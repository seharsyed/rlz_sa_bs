#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Phase 1: Read FASTA or plain input and create an A/C/G/T-only reference
// -----------------------------------------------------------------------------

struct CleaningStats {
    std::uint64_t sequence_length = 0;
    std::uint64_t replacement_count = 0;
    std::uint64_t header_lines_removed = 0;
    std::uint64_t whitespace_characters_removed = 0;
    std::uint64_t nul_bytes_removed = 0;
};

static std::vector<char> clean_reference(
    const std::string& input_path,
    CleaningStats& stats) {

    std::ifstream input(input_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open input reference: " + input_path);
    }

    std::vector<char> cleaned;
    input.seekg(0, std::ios::end);
    const std::streamoff input_size = input.tellg();
    input.seekg(0, std::ios::beg);

    if (input_size > 0) {
        cleaned.reserve(static_cast<std::size_t>(input_size));
    }

    bool at_line_start = true;
    bool inside_header = false;
    char raw = '\0';

    while (input.get(raw)) {
        const unsigned char symbol = static_cast<unsigned char>(raw);

        if (symbol == 0) {
            ++stats.nul_bytes_removed;
            continue;
        }

        if (inside_header) {
            if (raw == '\n' || raw == '\r') {
                inside_header = false;
                at_line_start = true;
            }
            continue;
        }

        if (at_line_start && raw == '>') {
            inside_header = true;
            ++stats.header_lines_removed;
            continue;
        }

        if (std::isspace(symbol) != 0) {
            ++stats.whitespace_characters_removed;
            if (raw == '\n' || raw == '\r') {
                at_line_start = true;
            }
            continue;
        }

        at_line_start = false;

        if (raw == 'A' || raw == 'C' || raw == 'G' || raw == 'T') {
            cleaned.push_back(raw);
        } else {
            cleaned.push_back('A');
            ++stats.replacement_count;
        }
    }

    stats.sequence_length = cleaned.size();

    if (cleaned.empty()) {
        throw std::runtime_error("The input contains no sequence characters.");
    }

    return cleaned;
}

// -----------------------------------------------------------------------------
// Phase 1 output: Write the cleaned reference and cleaning statistics
// -----------------------------------------------------------------------------

static void write_cleaned_reference(
    const std::string& output_path,
    const std::vector<char>& cleaned) {

    std::ofstream output(output_path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("Cannot create cleaned reference: " + output_path);
    }

    output.write(cleaned.data(), static_cast<std::streamsize>(cleaned.size()));
    if (!output) {
        throw std::runtime_error("Failed while writing cleaned reference: " + output_path);
    }
}

static void write_cleaning_csv(
    const std::string& output_path,
    const CleaningStats& stats) {

    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("Cannot create cleaning CSV: " + output_path);
    }

    output << "metric,value\n";
    output << "sequence_length," << stats.sequence_length << '\n';
    output << "replacement_count," << stats.replacement_count << '\n';
    output << "header_lines_removed," << stats.header_lines_removed << '\n';
    output << "whitespace_characters_removed,"
           << stats.whitespace_characters_removed << '\n';
    output << "nul_bytes_removed," << stats.nul_bytes_removed << '\n';
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr
            << "Usage: " << argv[0]
            << " <input_fasta_or_plain> <cleaned_reference> <cleaning_csv>\n";
        return EXIT_FAILURE;
    }

    try {
        CleaningStats stats;
        const std::vector<char> cleaned = clean_reference(argv[1], stats);

        write_cleaned_reference(argv[2], cleaned);
        write_cleaning_csv(argv[3], stats);

        std::cout << "cleaned_reference=" << argv[2] << '\n';
        std::cout << "sequence_length=" << stats.sequence_length << '\n';
        std::cout << "replacement_count=" << stats.replacement_count << '\n';
        std::cout << "header_lines_removed=" << stats.header_lines_removed << '\n';
        std::cout << "whitespace_characters_removed="
                  << stats.whitespace_characters_removed << '\n';
        std::cout << "nul_bytes_removed="
                  << stats.nul_bytes_removed << '\n';

        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
