#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <string>

void help() {
    std::cout << "Relative Lempel-Ziv decompression.\n\n";
    std::cout << "Usage: decompress -r <reference_file> -s <compressed_sequences_file> -o <output_path>\n";
    std::cout << "   <reference_file>   Path to reference file (not bwt!).\n";
    std::cout << "   <sequences_file>   Path to files containing compressed sequences.\n";
    std::cout << "   <output_path>      Path to directory where the decompressed files will be stored.\n";
    std::cout << "Both files are required.\n";
    std::cout << "The resulting decompressed files will carry the same name as input files, appended with 'decompressed.txt'.\n\n";
    std::cout << "Example: decompress -r reference.txt -s sequences.txt -o ../decompressed/" << std::endl;
    exit(0);
}

void decompress(const std::string& reference_file_path, const std::string& sequences_file_path, const std::string& output_path) {

    // Read the reference file
    std::ifstream ref_file(reference_file_path);
    std::string reference = "";
    std::string line;
    constexpr std::string_view suffix = "_phrases.bin";

    while (std::getline(ref_file, line)) {
        if (line[0] == '>' || line[0] == '\n')
            continue;
        reference += line;
    }

    // Save the file name of the compressed sequences
    std::ifstream seq_file(sequences_file_path);
    if (!seq_file) {
        throw std::runtime_error("Cannot open file: " + std::string(sequences_file_path));
    }

    std::vector<std::string> seq_files;
    std::string p;

    bool triple = false;

    while (seq_file >> p)
        seq_files.push_back(p);

    //#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < seq_files.size(); i++) {

        std::string filename =
            std::filesystem::path(seq_files[i]).filename().string();
        if (filename.ends_with(suffix))
            filename.resize(filename.size() - suffix.size());

        std::ofstream out_file(output_path + "decompressed_" + filename);

        std::fstream in_file;
        in_file.open(seq_files[i], std::ios::binary | std::ios::in);
        if (!in_file) {
            throw std::runtime_error("Cannot open file: " + std::string(seq_files[i]));
        }

        uint8_t des;
        in_file.read(reinterpret_cast<char*>(&des), 1);

        if (des > 0)
            triple = true;

        uint64_t length, position;
        uint64_t mask = 1ULL << 63;

        // Read the paris of length and position
        while (in_file.read(reinterpret_cast<char*>(&length),sizeof(uint64_t)) &&
            in_file.read(reinterpret_cast<char*>(&position),sizeof(uint64_t))) {

            if (triple)
                in_file.read(reinterpret_cast<char*>(&des),sizeof(uint8_t));

            if (mask & length) {
                length -= mask;
                for (uint64_t j = 0; j < length; j++) {
                    out_file << (char)position;
                }
                continue;
            }

            // For each pair, output the corresponding reference subsequence into the outputfile
            for (uint64_t j = 0; j < length; j++) {
                out_file << reference[(j+position)%reference.size()];
            }
            if (triple)
                out_file << (char)des;
        }

        //out_file << "\n";
        out_file.close();
    }
}

int main(int argc, char const* argv[]) {
    if (argc < 3) {
        std::cerr << "Input and pattern files are required\n" << std::endl;
        help();
    }
    std::string reference_file_path = "";
    std::string sequences_file_path = "";
    std::string output_path = "../";

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-r") == 0) {
            reference_file_path = argv[++i];
        } else if (strcmp(argv[i], "-s") == 0) {
            sequences_file_path = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0) {
            output_path = argv[++i];
        } else {
            help();
        }
    }


    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;
    auto start = high_resolution_clock::now();

    decompress(reference_file_path, sequences_file_path, output_path);

    auto end = high_resolution_clock::now();
    double time = duration_cast<milliseconds>(end - start).count();

    std::cerr << "Decompression time: " << time << "ms\n";
}