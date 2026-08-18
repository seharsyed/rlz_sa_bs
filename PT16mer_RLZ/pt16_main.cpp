#include "parser.hpp"
#include "pt16_build.hpp"
#include "pt16_rlz.hpp"
#include "pt16_utils.hpp"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using Symbol = unsigned char;
using SAType = std::uint32_t;

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    std::string temporary_directory;

    try {
        std::cerr << "========================================" << std::endl;
        std::cerr << "PT16 RLZ experiment" << std::endl;
        std::cerr << "========================================" << std::endl;

        // ---------- Load experiment inputs ----------

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

        // ---------- Baseline ----------

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[5] BASELINE RLZ" << std::endl;
        std::cerr << "========================================" << std::endl;

        temporary_directory = args.results + ".baseline_tmp";
        fs::remove_all(temporary_directory);
        fs::create_directories(temporary_directory);

        std::vector<BaselineResult> baseline_results;
        baseline_results.reserve(files.size());

        std::size_t total_input_bytes = 0;
        std::size_t total_baseline_phrases = 0;
        double total_baseline_ms = 0.0;

        std::cerr << "Baseline files: 0/" << files.size() << std::flush;

        for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
            const auto input = load_input<Symbol>(files[file_index]);

            Triples baseline;

            const double baseline_ms = time_ms([&] {
                baseline = lzFactorize<Symbol, SAType>(input, reference, suffix_array);
            });

            const std::string factor_file = temporary_directory + "/baseline_" + std::to_string(file_index) + ".bin";
            write_factor_file(factor_file, baseline);

            baseline_results.push_back({
                files[file_index],
                input.size(),
                baseline_ms,
                baseline.size(),
                factor_file
            });

            total_input_bytes += input.size();
            total_baseline_phrases += baseline.size();
            total_baseline_ms += baseline_ms;

            std::cerr << "\rBaseline files: " << file_index + 1 << "/" << files.size() << std::flush;
        }

        std::cerr << std::endl;
        std::cerr << "Baseline complete." << std::endl;
        std::cerr << "Total baseline time: " << total_baseline_ms << " ms" << std::endl;

        // ---------- PT16 preprocessing ----------

std::cerr << std::endl;
std::cerr << "========================================" << std::endl;
std::cerr << "[6] PT16 TABLE" << std::endl;
std::cerr << "========================================" << std::endl;

double pt16_build_ms = 0.0;

// Always rebuild the PT16 table for the current experiment.
if (fs::exists(pt16_path)) {
    std::cerr << "Removing existing PT16 table..." << std::endl;
    fs::remove(pt16_path);
}

std::cerr << "Building H, L and interval starts..." << std::endl;

pt16_build_ms = time_ms([&] {
    build_pt16_table(reference, suffix_array, pt16_path);
});

std::cerr << "PT16 table built." << std::endl;
std::cerr << "PT16 build time: " << pt16_build_ms << " ms" << std::endl;
std::cerr << "PT16 build time: " << pt16_build_ms / 1000.0 << " s" << std::endl;
std::cerr << "Table: " << pt16_path << std::endl;

        // ---------- Load PT16 ----------

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[7] LOAD PT16" << std::endl;
        std::cerr << "========================================" << std::endl;

        std::cerr << "Loading PT16 table..." << std::endl;

        PT16RLZParser<Symbol, SAType> parser(reference, suffix_array, pt16_path);

        std::cerr << "PT16 loaded." << std::endl;
        std::cerr << "Entries: " << parser.stats().entries << std::endl;
        std::cerr << "PT16 memory: "
                  << static_cast<double>(parser.stats().approx_bytes) / (1024.0 * 1024.0)
                  << " MB"
                  << std::endl;

        // ---------- CSV ----------

        std::cerr << std::endl;
        std::cerr << "[8] Opening results CSV..." << std::endl;

        CSVWriter csv(args.results);

        std::cerr << "    " << args.results << std::endl;

        std::size_t total_pt16_phrases = 0;
        double total_pt16_ms = 0.0;
        bool all_equal = true;

        auto previous_stats = parser.stats();

        // ---------- PT16 parsing ----------

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[9] PT16 RLZ" << std::endl;
        std::cerr << "========================================" << std::endl;

        std::cerr << "PT16 files: 0/" << baseline_results.size() << std::flush;

        for (std::size_t file_index = 0; file_index < baseline_results.size(); ++file_index) {
            const BaselineResult& baseline_result = baseline_results[file_index];

            const auto input = load_input<Symbol>(baseline_result.filename);

            Triples pt16;

            const double pt16_ms = time_ms([&] {
                pt16 = parser.lzFactorize(input);
            });

            const bool equal = factor_file_equals(baseline_result.factor_file, pt16);

            const auto current_stats = parser.stats();
            const PT16Delta file_stats = stats_difference(previous_stats, current_stats);

            csv.write_row(
                baseline_result,
                pt16_ms,
                pt16.size(),
                equal,
                file_stats,
                current_stats.entries,
                current_stats.approx_bytes
            );

            previous_stats = current_stats;

            all_equal = all_equal && equal;
            total_pt16_phrases += pt16.size();
            total_pt16_ms += pt16_ms;

            fs::remove(baseline_result.factor_file);

            std::cerr << "\rPT16 files: " << file_index + 1 << "/" << baseline_results.size() << std::flush;
        }

        std::cerr << std::endl;
        std::cerr << "PT16 complete." << std::endl;

        // ---------- Final summary ----------

        std::cerr << std::endl;
        std::cerr << "========================================" << std::endl;
        std::cerr << "[10] FINAL SUMMARY" << std::endl;
        std::cerr << "========================================" << std::endl;

        const BenchmarkSummary summary{
            baseline_results.size(),
            total_input_bytes,
            total_baseline_phrases,
            total_pt16_phrases,
            total_baseline_ms,
            total_pt16_ms,
            pt16_build_ms,
            all_equal
        };

        csv.write_summary(summary, parser.stats());

        fs::remove_all(temporary_directory);

        const double overall_speedup = total_pt16_ms == 0.0 ? 0.0 : total_baseline_ms / total_pt16_ms;

        const std::size_t total_queries = parser.stats().hits + parser.stats().misses;
        const double total_hit_rate = total_queries == 0
            ? 0.0
            : static_cast<double>(parser.stats().hits) / static_cast<double>(total_queries);




            const std::size_t total_hit_types =
            parser.stats().singleton_hits + parser.stats().range_hits;

            const double singleton_hit_rate =
            total_hit_types == 0
            ? 0.0
            : static_cast<double>(parser.stats().singleton_hits) /
            static_cast<double>(total_hit_types);


        std::cout << "results=" << args.results << std::endl;
        std::cout << "total_baseline_ms=" << total_baseline_ms << std::endl;
        std::cout << "total_pt16_ms=" << total_pt16_ms << std::endl;
        std::cout << "overall_speedup=" << overall_speedup << std::endl;
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
        std::cout << "all_outputs_equal=" << (all_equal ? "YES" : "NO") << std::endl;

        return all_equal ? EXIT_SUCCESS : 4;
    }
    catch (const std::exception& error) {
        if (!temporary_directory.empty()) {
            fs::remove_all(temporary_directory);
        }

        std::cerr << std::endl;
        std::cerr << "ERROR: " << error.what() << std::endl;

        return EXIT_FAILURE;
    }
}