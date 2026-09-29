#include <cstring>
#include <fstream>
#include <iostream>

#include "include/debug.hpp"
#include "include/reader.hpp"
#include "include/types.hpp"
#include <filesystem>

void help() {
    std::cout << "Create BWT data structure.\n\n"
        << "Usage: make_bwt [options] -i <bwt_file_name> -sa <sa_file_name> -o <output_file_name>\n"
        << "   <bwt_file_name>      File containing bwt of the reference.\n"
        << "   <sa_file_name>       File containing suffix array (or GCA).\n"
        << "   <output_file_name>   File where the output will be written to.\n"
        << "   --rle                Use run length encoding.\n"
        << "   -n                   Strip new line characters from input.\n\n";
    std::cout 
        << "Output file is required.\n";
    std::cout << "Example: make_bwt -i bwt.txt -sa sa.txt -o bwt.bwt" << std::endl;
    std::cout << "Please run the command again if the number of dense blocks shows 0 and -r option is on." << std::endl;
    exit(0);
}

typedef bbwt::two_byte_build<> bwt_type_a;
typedef bbwt::non_rle_build<> bwt_type_e;

template <class bwt_t>
void build(char const* argv[], size_t in_file_loc, size_t sa_file_loc,
            size_t out_file_loc, bool strip_new_line) {

    const std::string input = argv[in_file_loc];
    std::ifstream in(input);

    typename bwt_t::builder b(argv[out_file_loc], argv[sa_file_loc] );

    bbwt::file_reader<typename bwt_t::alphabet_type> reader(&in);

    for (auto it : reader) {
        if (strip_new_line && (it.head == '\n')) {
            continue;
        }
        if (it.length > 0) {
            b.append(it.head, it.length);
        }
    }

    b.finalize();
}

int main(int argc, char const* argv[]) {
    if (argc < 2) {
        help();
    }

    size_t in_file_loc = 0;
    size_t out_file_loc = 0;
    size_t sa_file_loc = 0;

    bool strip_new_line = false;
    bool rle = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            in_file_loc = ++i;
        } else if (strcmp(argv[i], "-sa") == 0) {
            sa_file_loc = ++i;
        } else if (strcmp(argv[i], "-n") == 0) {
            strip_new_line = true;
        } else if (strcmp(argv[i], "--rle") == 0) {
             rle = true;
        } else if (strcmp(argv[i], "-o") == 0) {
            out_file_loc = ++i;
        } else {
            help();
        }
    }
    if (in_file_loc == 0) {
        throw std::runtime_error("Input file is required");
    }
    if (sa_file_loc == 0) {
        throw std::runtime_error("SA file is required");
    }
    if (out_file_loc == 0) {
        throw std::runtime_error("Output file is required");
    }

    std::string cmd;

    if (rle) {
        if (strip_new_line) {
            cmd = "make\n./make_alphabet_header -n -i " + (std::string)argv[in_file_loc] + " > include/custom_alphabet.hpp";
        } else {
            cmd = "make\n./make_alphabet_header -i " + (std::string)argv[in_file_loc] + " > include/custom_alphabet.hpp";
        }
        int ret = std::system(cmd.c_str());
        ret = std::system(cmd.c_str());
        if (ret != 0) {
            throw std::runtime_error("Alphabet build failed");
        }

        build<bwt_type_a>(argv, in_file_loc, sa_file_loc,  out_file_loc, strip_new_line);
    } else {
        build<bwt_type_e>(argv, in_file_loc, sa_file_loc, out_file_loc, strip_new_line);
    }

    return 0;
}