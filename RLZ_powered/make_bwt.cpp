#include <cstring>
#include <fstream>
#include <iostream>

#include "include/debug.hpp"
#include "include/reader.hpp"
#include "include/types.hpp"
#include <filesystem>

void help() {
    std::cout << "Create powered BWT data structure.\n\n"
        << "Usage: make_bwt [options] -i <bwt_file_name> -sa <sa_file_name> -o <output_file_name>\n"
        << "   <bwt_file_name>      File containing bwt of the reference.\n"
        << "   <sa_file_name>       File containing suffix array (or GCA).\n"
        << "   <output_file_name>   File where the output will be written to.\n"
        << "   --original           Keep the original alphabet.\n";
    std::cout 
        << "Output file is required.\n";
    std::cout << "Example: make_bwt -i bwt.txt -sa sa.txt -o bwt.bwt" << std::endl;
    std::cout << "Please run the command again if the number of dense blocks shows 0 and -r option is on." << std::endl;
    exit(0);
}

typedef bbwt::two_byte_build<> bwt_type_a;
typedef bbwt::non_rle_build<> bwt_type_e;

template <class bwt_t>
void build(const std::string& input_filename, const std::string& sa_filename, const std::string& output_filename) {

    std::ifstream in(input_filename);

    typename bwt_t::builder b(output_filename, sa_filename);

    bbwt::file_reader<typename bwt_t::alphabet_type> reader(&in);

    for (auto it : reader) {
        if (it.length > 0) {
            b.append(it.head, it.length);
        }
    }
    b.finalize();
}

template <class bwt_type>
void transform() {
    auto constexpr encode_base = [](unsigned char ch) noexcept ->
                                uint8_t { return static_cast<uint8_t>
                                    ( (ch > 'A') + (ch > 'C') + (ch > 'G')); };

    // Read the FM-index
    bwt_type bwt("tmp.bwt");
    std::ofstream out_file("tmp.txt");

    for (uint64_t i = 0; i < bwt.size(); i++) {
        uint8_t x = encode_base(bwt.at(i));
        uint64_t pos = bwt.LF(i);

        // Incode multiple characters into one meta-character
        for (int j = 1; j < 4; j++) {
            x += encode_base(bwt.at(pos))<<(j<<1);
            pos = bwt.LF(pos);
        }

        // Output the metacharacter into the output file
        out_file << x;
    }
    out_file.close();
}

int main(int argc, char const* argv[]) {
    if (argc < 2) {
        help();
    }

    std::string in_filename;
    std::string sa_filename;
    std::string out_filename;
    bool original = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            in_filename = argv[++i];
        } else if (strcmp(argv[i], "-sa") == 0) {
            sa_filename = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0) {
            out_filename = argv[++i];
        } else if (strcmp(argv[i], "--original") == 0) {
            original = true;
        } else {
            help();
        }
    }
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto start = high_resolution_clock::now();

    std::string cmd = "make\n./make_alphabet_header -i " + in_filename + " > include/custom_alphabet.hpp";

    int ret = std::system(cmd.c_str());
    ret = std::system(cmd.c_str());
    if (ret != 0) {
        throw std::runtime_error("Alphabet build failed");
    }

    if (original) {
        build<bwt_type_a>(in_filename, sa_filename, out_filename);

    } else {
        build<bwt_type_a>(in_filename, sa_filename, "tmp.bwt");

        transform<bbwt::two_byte<>>();

        build<bwt_type_e>("tmp.txt", sa_filename, out_filename);
        std::filesystem::remove("tmp.txt");
        std::filesystem::remove("tmp.bwt");
        std::filesystem::remove("tmp_data.bwt");
    }

    auto end = high_resolution_clock::now();
    double time = duration_cast<milliseconds>(end - start).count();

    std::cerr << "\nBuilding the FM-index in " << out_filename << " took " << time << "ms" << std::endl;
    return 0;
}