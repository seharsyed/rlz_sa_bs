#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct CleaningStats {
    std::uint64_t sequence_length = 0;
    std::uint64_t replacement_count = 0;
    std::uint64_t header_lines_removed = 0;
    std::uint64_t whitespace_characters_removed = 0;
    std::uint64_t nul_bytes_removed = 0;
};

static std::vector<char> clean_sequence(
    const std::string& input_path,
    CleaningStats& stats) {

    std::ifstream input(input_path, std::ios::binary);

    if (!input) {
        throw std::runtime_error("Cannot open input file: " + input_path);
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
        throw std::runtime_error(
            "Input file contains no sequence characters: " + input_path);
    }

    return cleaned;
}

static void write_cleaned_sequence(
    const fs::path& output_path,
    const std::vector<char>& cleaned) {

    std::ofstream output(output_path, std::ios::binary);

    if (!output) {
        throw std::runtime_error(
            "Cannot create cleaned input file: " + output_path.string());
    }

    output.write(
        cleaned.data(),
        static_cast<std::streamsize>(cleaned.size()));

    if (!output) {
        throw std::runtime_error(
            "Failed while writing cleaned input file: " +
            output_path.string());
    }
}

static std::vector<std::string> read_input_list(
    const std::string& input_list_path) {

    std::ifstream input_list(input_list_path);

    if (!input_list) {
        throw std::runtime_error(
            "Cannot open input list: " + input_list_path);
    }

    std::vector<std::string> paths;
    std::string line;

    while (std::getline(input_list, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (!line.empty()) {
            paths.push_back(line);
        }
    }

    if (paths.empty()) {
        throw std::runtime_error(
            "Input list contains no file paths: " + input_list_path);
    }

    return paths;
}

static void write_csv_header(std::ofstream& output) {
    output
        << "input_file,"
        << "cleaned_file,"
        << "sequence_length,"
        << "replacement_count,"
        << "header_lines_removed,"
        << "whitespace_characters_removed,"
        << "nul_bytes_removed\n";
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr
            << "Usage: " << argv[0]
            << " <input_list.txt>"
            << " <cleaned_output_folder>"
            << " <cleaning_stats.csv>\n";

        return EXIT_FAILURE;
    }

    try {
        const std::string input_list_path = argv[1];
        const fs::path output_folder = argv[2];
        const fs::path statistics_path = argv[3];

        fs::create_directories(output_folder);

        const std::vector<std::string> input_paths =
            read_input_list(input_list_path);

        const fs::path cleaned_list_path =
            output_folder.parent_path() / "cleaned_input_list.txt";

        std::ofstream cleaned_list(cleaned_list_path);

        if (!cleaned_list) {
            throw std::runtime_error(
                "Cannot create cleaned input list: " +
                cleaned_list_path.string());
        }

        std::ofstream statistics(statistics_path);

        if (!statistics) {
            throw std::runtime_error(
                "Cannot create cleaning statistics CSV: " +
                statistics_path.string());
        }

        write_csv_header(statistics);

        std::uint64_t total_sequence_length = 0;
        std::uint64_t total_replacements = 0;

        for (const std::string& input_path : input_paths) {
            CleaningStats stats;

            const std::vector<char> cleaned =
                clean_sequence(input_path, stats);

            const fs::path input_filename =
                fs::path(input_path).filename();

            const fs::path output_path =
                output_folder / input_filename;

            write_cleaned_sequence(output_path, cleaned);

            cleaned_list << output_path.string() << '\n';

            statistics
                << '"' << input_path << '"' << ','
                << '"' << output_path.string() << '"' << ','
                << stats.sequence_length << ','
                << stats.replacement_count << ','
                << stats.header_lines_removed << ','
                << stats.whitespace_characters_removed << ','
                << stats.nul_bytes_removed << '\n';

            total_sequence_length += stats.sequence_length;
            total_replacements += stats.replacement_count;

            std::cout
                << "cleaned=" << output_path
                << " sequence_length=" << stats.sequence_length
                << " replacements=" << stats.replacement_count
                << '\n';
        }

        std::cout << "files_cleaned=" << input_paths.size() << '\n';
        std::cout
            << "total_sequence_length="
            << total_sequence_length << '\n';
        std::cout
            << "total_replacement_count="
            << total_replacements << '\n';
        std::cout
            << "cleaned_input_list="
            << cleaned_list_path << '\n';

        return EXIT_SUCCESS;

    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}