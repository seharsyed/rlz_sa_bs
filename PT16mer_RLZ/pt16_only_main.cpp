#include "parser.hpp"
#include "pt16_build.hpp"
#include "pt16_rlz.hpp"
#include "pt16_utils.hpp"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using Symbol = unsigned char;
using SAType = std::uint32_t;

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    try {
        std::cerr << "========================================" << std::endl;
        std::cerr << "PT16 RLZ" << std::endl;
        std::cerr << "========================================" << std::endl;

        std::cerr << "[1] Reading arguments..." << std::endl;
        const Args args = parse_args(argc, argv);

        std::cerr << "[2] Loading reference..." << std::endl;
        const auto reference = load_reference<Symbol>(args.reference);
        std::cerr << "    Reference loaded: " << reference.size() << " bytes" << std::endl;

        std::cerr << "[3] Loading suffix array..." << std::endl;
        const auto suffix_array = load_suffix_array<SAType>(args.suffix_array);
        std::cerr << "    Suffix array loaded: " << suffix_array.size() << " entries" << std::endl;

        std::cerr << "[4] Loading input list..." << std::endl;
        const auto files = load_input_list(args.filenames);
        std::cerr << "    Input files: " << files.size() << std::endl;

        const std::string pt16_path = get_pt16_path(args);

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[5] PT16 TABLE" << std::endl;
        std::cerr << "========================================" << std::endl;

        if (fs::exists(pt16_path)) {
            std::cerr << "Removing existing PT16 table..." << std::endl;
            fs::remove(pt16_path);
        }

        std::cerr << "Building H, L and interval starts..." << std::endl;

        const double pt16_build_ms = time_ms([&] {
            build_pt16_table(reference, suffix_array, pt16_path);
        });

        std::cerr << "PT16 table built." << std::endl;
        std::cerr << "PT16 build time: " << pt16_build_ms << " ms" << std::endl;
        std::cerr << "Table: " << pt16_path << std::endl;

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[6] LOAD PT16" << std::endl;
        std::cerr << "========================================" << std::endl;

        PT16RLZParser<Symbol, SAType> parser(reference, suffix_array, pt16_path);

        std::cerr << "PT16 loaded." << std::endl;
        std::cerr << "Entries: " << parser.stats().entries << std::endl;
        std::cerr << "PT16 memory: "
                  << static_cast<double>(parser.stats().approx_bytes) / (1024.0 * 1024.0)
                  << " MB" << std::endl;

        std::ofstream csv(args.results);

        if (!csv) {
            throw std::runtime_error("Could not open results file: " + args.results);
        }

        csv << "filename,input_bytes,pt16_ms,pt16_phrases\n";

        std::size_t total_input_bytes = 0;
        std::size_t total_pt16_phrases = 0;
        double total_pt16_ms = 0.0;

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[7] PT16 RLZ" << std::endl;
        std::cerr << "========================================" << std::endl;

        std::cerr << "PT16 files: 0/" << files.size() << std::flush;

        for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
            const auto input = load_input<Symbol>(files[file_index]);

            Triples pt16;

            const double pt16_ms = time_ms([&] {
                pt16 = parser.lzFactorize(input);
            });

            csv << fs::path(files[file_index]).filename().string() << ","
                << input.size() << ","
                << pt16_ms << ","
                << pt16.size() << "\n";

            total_input_bytes += input.size();
            total_pt16_phrases += pt16.size();
            total_pt16_ms += pt16_ms;

            std::cerr << "\rPT16 files: " << file_index + 1 << "/" << files.size() << std::flush;
        }

        std::cerr << std::endl;
        std::cerr << "PT16 complete." << std::endl;

        const std::size_t total_queries = parser.stats().hits + parser.stats().misses;

        const double total_hit_rate = total_queries == 0
            ? 0.0
            : static_cast<double>(parser.stats().hits) /
              static_cast<double>(total_queries);

        const std::size_t total_hit_types =
            parser.stats().singleton_hits + parser.stats().range_hits;

        const double singleton_hit_rate = total_hit_types == 0
            ? 0.0
            : static_cast<double>(parser.stats().singleton_hits) /
              static_cast<double>(total_hit_types);

        csv << "\n";
        csv << "summary\n";
        csv << "processed_files," << files.size() << "\n";
        csv << "total_input_bytes," << total_input_bytes << "\n";
        csv << "total_pt16_phrases," << total_pt16_phrases << "\n";
        csv << "total_pt16_ms," << total_pt16_ms << "\n";
        csv << "avg_pt16_ms,"
            << (files.empty() ? 0.0 : total_pt16_ms / static_cast<double>(files.size()))
            << "\n";
        csv << "pt16_build_ms," << pt16_build_ms << "\n";
        csv << "pt16_hits," << parser.stats().hits << "\n";
        csv << "pt16_misses," << parser.stats().misses << "\n";
        csv << "pt16_hit_rate," << total_hit_rate << "\n";
        csv << "pt16_entries," << parser.stats().entries << "\n";
        csv << "singleton_hits," << parser.stats().singleton_hits << "\n";
        csv << "range_hits," << parser.stats().range_hits << "\n";
        csv << "singleton_hit_rate," << singleton_hit_rate << "\n";
        csv << "pt16_MB,"
            << static_cast<double>(parser.stats().approx_bytes) / (1024.0 * 1024.0)
            << "\n";
        csv << "peak_RSS_MB," << peak_rss_mb() << "\n";

        std::cout << "results=" << args.results << std::endl;
        std::cout << "processed_files=" << files.size() << std::endl;
        std::cout << "total_input_bytes=" << total_input_bytes << std::endl;
        std::cout << "total_pt16_phrases=" << total_pt16_phrases << std::endl;
        std::cout << "total_pt16_ms=" << total_pt16_ms << std::endl;
        std::cout << "avg_pt16_ms="
                  << (files.empty() ? 0.0 : total_pt16_ms / static_cast<double>(files.size()))
                  << std::endl;
        std::cout << "pt16_hits=" << parser.stats().hits << std::endl;
        std::cout << "pt16_misses=" << parser.stats().misses << std::endl;
        std::cout << "pt16_hit_rate=" << total_hit_rate << std::endl;
        std::cout << "pt16_entries=" << parser.stats().entries << std::endl;
        std::cout << "singleton_hits=" << parser.stats().singleton_hits << std::endl;
        std::cout << "range_hits=" << parser.stats().range_hits << std::endl;
        std::cout << "singleton_hit_rate=" << singleton_hit_rate << std::endl;
        std::cout << "pt16_build_ms=" << pt16_build_ms << std::endl;
        std::cout << "pt16_MB="
                  << static_cast<double>(parser.stats().approx_bytes) / (1024.0 * 1024.0)
                  << std::endl;
        std::cout << "peak_RSS_MB=" << peak_rss_mb() << std::endl;

        return EXIT_SUCCESS;
    }
    catch (const std::exception& error) {
        std::cerr << std::endl;
        std::cerr << "ERROR: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
}