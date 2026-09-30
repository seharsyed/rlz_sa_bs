#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

bool compare_files(const std::string& file1, const std::string& file2)
{
    std::ifstream a(file1, std::ios::binary);
    std::ifstream b(file2, std::ios::binary);

    if (!a) {
        std::cerr << "Cannot open: " << file1 << '\n';
        return false;
    }

    if (!b) {
        std::cerr << "Cannot open: " << file2 << '\n';
        return false;
    }

    // Check sizes first
    a.seekg(0, std::ios::end);
    b.seekg(0, std::ios::end);

    const auto size_a = a.tellg();
    const auto size_b = b.tellg();

    if (size_a != size_b) {
        std::cout << "DIFFERENT: " << file1
                  << " (size " << size_a << ") vs "
                  << file2 << " (size " << size_b << ")\n";
        return false;
    }

    a.seekg(0);
    b.seekg(0);

    constexpr std::size_t BUFFER_SIZE = 1 << 20; // 1 MiB

    std::vector<char> buffer_a(BUFFER_SIZE);
    std::vector<char> buffer_b(BUFFER_SIZE);

    std::uint64_t offset = 0;

    while (a && b) {
        a.read(buffer_a.data(), BUFFER_SIZE);
        b.read(buffer_b.data(), BUFFER_SIZE);

        const std::streamsize read_a = a.gcount();
        const std::streamsize read_b = b.gcount();

        if (read_a != read_b) {
            std::cout << "DIFFERENT: " << file1 << '\n';
            return false;
        }

        for (std::streamsize i = 0; i < read_a; ++i) {
            if (buffer_a[i] != buffer_b[i]) {
                std::cout << "DIFFERENT: " << file1
                          << " vs " << file2
                          << " at byte " << (offset + i)
                          << " (0x"
                          << std::hex
                          << (static_cast<unsigned>(static_cast<unsigned char>(buffer_a[i])))
                          << " vs 0x"
                          << (static_cast<unsigned>(static_cast<unsigned char>(buffer_b[i])))
                          << std::dec
                          << ")\n";

                return false;
            }
        }

        offset += read_a;
    }

    return true;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_files.txt " << " decompressed_files.txt\n";
        return 1;
    }

    std::ifstream input(argv[1]);
    std::ifstream decompressed(argv[2]);

    if (!input) {
        std::cerr << "Cannot open input list: " << argv[1] << '\n';
        return 1;
    }

    std::string input_file, decompressed_file;
    std::size_t total = 0;
    std::size_t identical = 0;
    std::size_t different = 0;

    while (std::getline(input, input_file) && std::getline(decompressed, decompressed_file)) {

        ++total;

        if (compare_files(input_file, decompressed_file)) {
            ++identical;
        } else {
            std::cout << "DIFFERENT: " << input_file << '\n';
            ++different;
        }
    }

    std::cout << "\nSummary:\n";
    std::cout << "  Files checked : " << total << '\n';
    std::cout << "  Identical     : " << identical << '\n';
    std::cout << "  Different     : " << different << '\n';

    return different == 0 ? 0 : 2;
}
