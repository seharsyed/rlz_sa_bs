#include "parser.hpp"
#include "pt16_build.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using Symbol = unsigned char;
using SAType = std::uint32_t;

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <reference> <suffix_array> <output_table>\n";
        return EXIT_FAILURE;
    }

    const std::string reference_path = argv[1];
    const std::string suffix_array_path = argv[2];
    const std::string output_path = argv[3];

    std::cerr << "Loading reference..." << std::endl;
    const auto reference = read_file<Symbol>(reference_path.c_str());

    std::cerr << "Loading suffix array..." << std::endl;
    const auto suffix_array = read_file<SAType>(suffix_array_path.c_str());

    std::cerr << "Building PT16 table..." << std::endl;

    const auto start = std::chrono::steady_clock::now();

    build_pt16_table(reference, suffix_array, output_path);

    const auto end = std::chrono::steady_clock::now();

    const double build_ms =
        std::chrono::duration<double, std::milli>(end - start).count();

    std::cout << "pt16_build_ms=" << build_ms << std::endl;
    std::cout << "pt16_build_seconds=" << build_ms / 1000.0 << std::endl;

    return EXIT_SUCCESS;
}