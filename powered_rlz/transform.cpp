#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <cstdlib>
#include <filesystem>

#include "include/types.hpp"
#include "include/reader.hpp"

void help() {
    std::cout << "Transform ordinary bwt to bwt with coded alphabet. \n\n";
    std::cout << "Usage: transform [options] -i <bwt_file> [-sa <sa_file>] -o <output_file>\n";
    std::cout << "   bwt_file       Path to bwt element.\n";
    std::cout << "   sa_file        Path to file containing sa (or gca) of the bwt file. Only required when -b flag is on.\n";
    std::cout << "   output_file    Path to file containing coded bwt text.\n";
    std::cout << "   -b             Build the bwt structure from the transformed text\n";
    std::cout << "   -e             The input BWT structure is not run length encoded.\n";
    std::cout << "Bwt file, output file and the size of alphabet are required.\n\n";
    std::cout << "Example: transform -i bwt.bwt -o bwt_coded.txt" << std::endl;
    exit(0);
}

int  code (const unsigned char ch) {
    return (ch >= 'C') + (ch >= 'G') + (ch >= 'T');
}

template <class bwt_type>
void transform(const std::string& in_file_path, const std::string& out_file_path) {

    // Read the FM-index
    bwt_type bwt(in_file_path);
    std::ofstream out_file(out_file_path);

    for (uint64_t i = 0; i < bwt.size(); i++) {
        uint8_t x = code(bwt.at(i));
        uint8_t mul = 2;
        uint64_t pos = bwt.LF(i);

        // Incode multiple characters into one meta-character
        for (int j = 1; j < 4; j++) {
            x += code(bwt.at(pos))<<mul;
            pos = bwt.LF(pos);
            mul += 2;
        }

        // Output the metacharacter into the output file
        out_file << x;
    }
}

void build_BWT(const std::string& input, const std::string& sa_file) {
    size_t loc = input.find_last_of('.');
    std::string output = input.substr(0, loc) + ".bwt";

    std::string cmd = "./make_bwt -i " + input + " -sa " + sa_file + " -o " + output;
    int ret = std::system(cmd.c_str());

    std::filesystem::remove(input);

    if (ret != 0) {
        throw std::runtime_error("BWT build failed");
    } else {
        std::cerr << "New transformed BWT written in file " << output << std::endl;
    }
}

int main(int argc, char const* argv[]) {
    if (argc < 5) {
        std::cerr << "Input and output files, and the number of coded letters are required\n" << std::endl;
        help();
    }

    std::string in_file_path = "";
    std::string sa_file_path = "";
    std::string out_file_path = "";
    bool build_bwt = false;
    bool rle = true;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0) {
            build_bwt = true;
        } else if (strcmp(argv[i], "-e") == 0) {
            rle = false;
        } else if (strcmp(argv[i], "-i") == 0) {
            in_file_path = argv[++i];
        } else if (strcmp(argv[i], "-sa") == 0) {
            sa_file_path = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0) {
            out_file_path = argv[++i];
        } else {
            help();
        }
    }

    if (!rle) {
        if (sa_file_path.size() == 0) {
            std::cerr << "SA file is required with -b flag on\n";
            exit(1);
        }

        transform<bbwt::non_rle<>>(in_file_path, out_file_path);
    } else {
        transform<bbwt::two_byte<>>(in_file_path, out_file_path);
    }

    if (build_bwt) {
        build_BWT(out_file_path, sa_file_path);
    } else {
        std::cerr << "New transformed BWT written in file " << out_file_path << std::endl;
    }
}