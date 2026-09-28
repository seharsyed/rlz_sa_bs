//#define VERB

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <filesystem>
#include <string>
#include <omp.h>
#include <bits/regex_constants.h>

#include "include/reader.hpp"
#include "include/types.hpp"

void help() {
    std::cout << "Relative Lempel-Ziv compression using BWT data structure with powered backward search.\n\n";
    std::cout << "Usage: rlz_parser [options] -r <reference_file> -s <sequences_file> -o <output_path>\n";
    std::cout << "   <reference_file>               Path to bwt file of the reference.\n";
    std::cout << "   <sequences_file>               Path to file containing the names of the sequence files. The names are separated with a newline character.\n";
    std::cout << "   <output_path>                  Path to directory where the phrase files will be stored.\n";
    std::cout << "   --rle                          The input file is run length encoded. This means the alphabet size is 4, not 256.\n";
    std::cout << "   --triples                      The output will contain triples (length, position, mismatching symbol), instead of tuples.\n";
    std::cout << "   --original                     Turn on when the bwt of the original alphabet is used, but it is not rle.\n";
    std::cout << "Bwt file and patterns file are required.\n\n";
    std::cout << "Example: rlz_parser -r bwt.bwt -s sequences.txt -o ../phrases/" << std::endl;
    exit(0);
}

template <class bwt_type, typename T>
std::pair<double, size_t> parse(const std::string& in_file_path, std::ifstream& sequences,
                        const std::string& output_path, uint8_t code_size, bool triple) {

    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::nanoseconds;
    using std::chrono::milliseconds;

    bwt_type bwt(in_file_path);

    std::vector<std::string> seq_files;
    std::string p;

    while (sequences >> p)
        seq_files.push_back(p);

    uint8_t description = triple == true;
    constexpr size_t phrase_size =
        std::tuple_size_v<T> == 3
            ? sizeof(uint64_t) * 2 + sizeof(uint8_t)
            : sizeof(uint64_t) * 2;

    auto total_start = high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < seq_files.size(); i++) {
        std::vector<T> parses;

        std::string filename =
            std::filesystem::path(seq_files[i]).filename().string();
        std::ofstream out_file(output_path + filename + "_phrases.bin", std::ios::binary);
        out_file.write(reinterpret_cast<const char*>(&description), sizeof(uint8_t));

        //read the sequence
        std::ifstream file(seq_files[i], std::ios::binary | std::ios::ate);
        if (!file) {
            throw std::runtime_error("Cannot open file: " + std::string(seq_files[i]));
        }
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);
        std::string sequence(size, '\0');
        file.read(sequence.data(), size);
        const auto sv = std::string_view(sequence.data(), size);

        //parses.reserve(size/1000);
        //parse the sequence
        if constexpr (std::tuple_size_v<T> == 3)
            bwt.parse_triples(sv, parses, code_size);
        else
            bwt.parse_tuples(sv, parses, code_size);

        // output the phrases
        std::vector<char> output;
        output.resize(parses.size() * phrase_size);

        char* dst = output.data();
        for (auto it = parses.rbegin(); it != parses.rend(); ++it) {
            const auto& phrase = *it;

            const uint64_t a = std::get<0>(phrase);
            const uint64_t b = std::get<1>(phrase);

            std::memcpy(dst, &a, sizeof(a));
            dst += sizeof(a);

            std::memcpy(dst, &b, sizeof(b));
            dst += sizeof(b);

            if constexpr (std::tuple_size_v<T> == 3) {
                const uint8_t c = std::get<2>(phrase);

                std::memcpy(dst, &c, sizeof(c));
                dst += sizeof(c);
            }
        }

        out_file.write(output.data(), output.size());
        out_file.close();
    }

    auto total_end = high_resolution_clock::now();
    double total_time = duration_cast<milliseconds>(total_end - total_start).count();

    return {total_time, seq_files.size()};
}

int main(int argc, char const* argv[]) {
    if (argc < 3) {
        std::cerr << "Input and pattern files are required\n" << std::endl;
        help();
    }
    std::string in_file_path = "";
    std::string sequences_file = "";
    std::string output_path = "";
    uint8_t code_size = 4;
    bool triple = false;
    bool rle = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--rle") == 0) {
            rle = true;
            code_size = 1;
        } else if (strcmp(argv[i], "--triples") == 0) {
            triple = true;
        }  else if (strcmp(argv[i], "--original") == 0) {
            code_size = 1;
        } else if (strcmp(argv[i], "-r") == 0) {
            in_file_path = argv[++i];
        } else if (strcmp(argv[i], "-s") == 0) {
            sequences_file = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0) {
            output_path = argv[++i];
        } else {
            help();
        }
    }

    std::ifstream s(sequences_file);
    if (!s) {
        throw std::runtime_error("Cannot open file: " + std::string(sequences_file));
    }

    std::pair<double, size_t> res;

    if (!rle) {
        if (triple) {
            res = parse<bbwt::non_rle<>, std::tuple<uint64_t, uint64_t, uint8_t>>(in_file_path, s, output_path, code_size, triple);
        } else {
            res = parse<bbwt::non_rle<>, std::tuple<uint64_t, uint64_t>>(in_file_path, s, output_path, code_size, triple);
        }
    } else {
        if (triple) {
            res = parse<bbwt::two_byte<>, std::tuple<uint64_t, uint64_t, uint8_t>>(in_file_path, s, output_path, code_size, triple);
        } else {
            res = parse<bbwt::two_byte<>, std::tuple<uint64_t, uint64_t>>(in_file_path, s, output_path, code_size, triple);
        }
    }
    std::cerr << "Mean compression time: " << res.first << " / " << res.second << " = " << (res.first / res.second) << "ms\n";

    s.close();
}