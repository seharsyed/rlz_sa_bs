#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "parser.hpp"
#include "variants/pt16_build_sassy.hpp"
#include "pt16_build_v2.hpp"
#include "pt16_rlz_v2.hpp"
#include "variants/pt16_sassy.hpp"
#include "pt16_utils.hpp"

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
    std::cerr << "    Reference loaded: " << reference.size() << " bytes ("
              << static_cast<double>(reference.size()) / (1024.0 * 1024.0)
              << " MB)" << std::endl;

    std::cerr << "[3] Loading suffix array..." << std::endl;
    const auto suffix_array = load_suffix_array<SAType>(args.suffix_array);
    const std::size_t suffix_array_bytes =
        suffix_array.size() * sizeof(SAType);
    std::cerr << "    Suffix array loaded: " << suffix_array_bytes
              << " bytes ("
              << static_cast<double>(suffix_array_bytes) / (1024.0 * 1024.0)
              << " MB)" << std::endl;

    std::cerr << "[4] Loading input list..." << std::endl;
    const auto files = load_input_list(args.filenames);
    std::cerr << "    Input files: " << files.size() << std::endl;

    const std::string pt16_path = get_pt16_path(args);

    // A sibling file: the sassy table is a different, self-contained format
    // (see pt16_sassy.hpp), not interchangeable with the one at pt16_path.
    const std::string sassy_path = pt16_path + ".sassy";

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

      const std::string factor_file = temporary_directory + "/baseline_" +
                                      std::to_string(file_index) + ".bin";
      write_factor_file(factor_file, baseline);

      baseline_results.push_back({files[file_index], input.size(), baseline_ms,
                                  baseline.size(), factor_file});

      total_input_bytes += input.size();
      total_baseline_phrases += baseline.size();
      total_baseline_ms += baseline_ms;

      std::cerr << "\rBaseline files: " << file_index + 1 << "/" << files.size()
                << std::flush;
    }

    std::cerr << std::endl;
    std::cerr << "Baseline complete." << std::endl;
    std::cerr << "Total baseline time: " << total_baseline_ms << " ms"
              << std::endl;

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

    std::cerr << "Building H, H_sa and lower-level entries..." << std::endl;

    pt16_build_ms =
        time_ms([&] { build_pt16_table(reference, suffix_array, pt16_path); });

    std::cerr << "PT16 table built." << std::endl;
    std::cerr << "PT16 build time: " << pt16_build_ms << " ms" << std::endl;
    std::cerr << "PT16 build time: " << pt16_build_ms / 1000.0 << " s"
              << std::endl;
    std::cerr << "Table: " << pt16_path << std::endl;

    // ---------- Load PT16 ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[7] LOAD PT16" << std::endl;
    std::cerr << "========================================" << std::endl;

    std::cerr << "Loading PT16 table..." << std::endl;

    PT16RLZParser<Symbol, SAType> parser(reference, suffix_array, pt16_path);

    std::cerr << "PT16 loaded." << std::endl;
    std::cerr << "PT16 entries: " << parser.stats().entries << std::endl;
    std::cerr << "PT16 memory: "
              << static_cast<double>(parser.stats().approx_bytes) /
                     (1024.0 * 1024.0)
              << " MB" << std::endl;

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

    for (std::size_t file_index = 0; file_index < baseline_results.size();
         ++file_index) {
      const BaselineResult& baseline_result = baseline_results[file_index];

      const auto input = load_input<Symbol>(baseline_result.filename);

      Triples pt16;

      const double pt16_ms = time_ms([&] { pt16 = parser.lzFactorize(input); });

      // Write PT16 factors after timing so correctness testing does not affect
      // runtime.
      write_factor_file(
          args.results + ".pt16_" + std::to_string(file_index) + ".bin", pt16);

      const bool equal = factor_file_equals(baseline_result.factor_file, pt16,
                                            input, reference);

      const auto current_stats = parser.stats();
      const PT16Delta file_stats =
          stats_difference(previous_stats, current_stats);

      csv.write_row(baseline_result, pt16_ms, pt16.size(), equal, file_stats,
                    current_stats.entries, current_stats.approx_bytes);

      previous_stats = current_stats;

      all_equal = all_equal && equal;
      total_pt16_phrases += pt16.size();
      total_pt16_ms += pt16_ms;

      // The baseline factor file is also needed by the sassy correctness
      // check below, so it is not removed here; see the sassy loop.

      std::cerr << "\rPT16 files: " << file_index + 1 << "/"
                << baseline_results.size() << std::flush;
    }

    std::cerr << std::endl;
    std::cerr << "PT16 complete." << std::endl;
    std::cerr << "Total PT16 time: " << total_pt16_ms << " ms" << std::endl;

    std::cerr << "PT16 entries: " << parser.stats().entries << std::endl;
    std::cerr << "PT16 hits: " << parser.stats().hits << std::endl;
    std::cerr << "PT16 misses: " << parser.stats().misses << std::endl;
    std::cerr << "PT16 singleton hits: " << parser.stats().singleton_hits
              << std::endl;
    std::cerr << "PT16 range hits: " << parser.stats().range_hits
              << std::endl;
    std::cerr << "PT16 memory: "
              << static_cast<double>(parser.stats().approx_bytes) /
                     (1024.0 * 1024.0)
              << " MB" << std::endl;

    // ---------- Sassy preprocessing ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[10] SASSY TABLE" << std::endl;
    std::cerr << "========================================" << std::endl;

    double sassy_build_ms = 0.0;

    // Always rebuild the sassy table for the current experiment, exactly as
    // for the PT16 table above.
    if (fs::exists(sassy_path)) {
      std::cerr << "Removing existing sassy table..." << std::endl;
      fs::remove(sassy_path);
    }

    std::cerr << "Building H, H_sa, L, sampled_sa and short suffixes..."
              << std::endl;

    sassy_build_ms = time_ms(
        [&] { build_pt16_sassy_table(reference, suffix_array, sassy_path); });

    std::cerr << "Sassy table built." << std::endl;
    std::cerr << "Sassy build time: " << sassy_build_ms << " ms" << std::endl;
    std::cerr << "Sassy build time: " << sassy_build_ms / 1000.0 << " s"
              << std::endl;
    std::cerr << "Table: " << sassy_path << std::endl;

    // ---------- Load sassy ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[11] LOAD SASSY" << std::endl;
    std::cerr << "========================================" << std::endl;

    std::cerr << "Loading sassy table..." << std::endl;

    PT16SassyLookup sassy(sassy_path);

    std::cerr << "Sassy loaded." << std::endl;
    std::cerr << "Sassy entries: " << sassy.stats().entries << std::endl;
    std::cerr << "Sassy sampled entries: " << sassy.stats().sampled_entries
              << std::endl;
    std::cerr << "Sassy short suffixes: " << sassy.stats().short_suffixes
              << std::endl;
    std::cerr << "Sassy memory: "
              << static_cast<double>(sassy.stats().approx_bytes) /
                     (1024.0 * 1024.0)
              << " MB" << std::endl;

    // ---------- Sassy CSV ----------

    std::cerr << std::endl;
    std::cerr << "[12] Opening sassy results CSV..." << std::endl;

    const std::string sassy_results = args.results + ".sassy.csv";
    CSVWriter sassy_csv(sassy_results);

    std::cerr << "    " << sassy_results << std::endl;

    std::size_t total_sassy_phrases = 0;
    double total_sassy_ms = 0.0;
    bool sassy_all_equal = true;

    auto previous_sassy_stats = sassy.stats();

    // ---------- Sassy parsing ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[13] SASSY RLZ" << std::endl;
    std::cerr << "========================================" << std::endl;

    std::cerr << "Sassy files: 0/" << baseline_results.size() << std::flush;

    for (std::size_t file_index = 0; file_index < baseline_results.size();
         ++file_index) {
      const BaselineResult& baseline_result = baseline_results[file_index];

      const auto input = load_input<Symbol>(baseline_result.filename);

      Triples sassy_factors;

      const double sassy_ms = time_ms(
          [&] { sassy_factors = sassy.lzFactorize(input, reference); });

      // Write sassy factors after timing so correctness testing does not
      // affect runtime.
      write_factor_file(
          args.results + ".sassy_" + std::to_string(file_index) + ".bin",
          sassy_factors);

      const bool equal = factor_file_equals(baseline_result.factor_file,
                                            sassy_factors, input, reference);

      const auto current_sassy_stats = sassy.stats();
      const PT16Delta file_sassy_stats =
          stats_difference(previous_sassy_stats, current_sassy_stats);

      sassy_csv.write_row(baseline_result, sassy_ms, sassy_factors.size(),
                          equal, file_sassy_stats,
                          current_sassy_stats.entries,
                          current_sassy_stats.approx_bytes);

      previous_sassy_stats = current_sassy_stats;

      sassy_all_equal = sassy_all_equal && equal;
      total_sassy_phrases += sassy_factors.size();
      total_sassy_ms += sassy_ms;

      std::cerr << "\rSassy files: " << file_index + 1 << "/"
                << baseline_results.size() << std::flush;
    }

    std::cerr << std::endl;
    std::cerr << "Sassy complete." << std::endl;
    std::cerr << "Total sassy time: " << total_sassy_ms << " ms" << std::endl;

    std::cerr << "Sassy entries: " << sassy.stats().entries << std::endl;
    std::cerr << "Sassy sampled entries: " << sassy.stats().sampled_entries
              << std::endl;
    std::cerr << "Sassy short suffixes: " << sassy.stats().short_suffixes
              << std::endl;
    std::cerr << "Sassy hits: " << sassy.stats().hits << std::endl;
    std::cerr << "Sassy misses: " << sassy.stats().misses << std::endl;
    std::cerr << "Sassy singleton hits: " << sassy.stats().singleton_hits
              << std::endl;
    std::cerr << "Sassy range hits: " << sassy.stats().range_hits
              << std::endl;
    std::cerr << "Sassy memory: "
              << static_cast<double>(sassy.stats().approx_bytes) /
                     (1024.0 * 1024.0)
              << " MB" << std::endl;

    // The temporary baseline factor files are needed by both the PT16 and
    // the sassy correctness checks above, so they are only removed now.
    for (const BaselineResult& baseline_result : baseline_results) {
      fs::remove(baseline_result.factor_file);
    }

    // ---------- Final summary ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[14] FINAL SUMMARY" << std::endl;
    std::cerr << "========================================" << std::endl;

    const BenchmarkSummary summary{baseline_results.size(), total_input_bytes,
                                   total_baseline_phrases,  total_pt16_phrases,
                                   total_baseline_ms,       total_pt16_ms,
                                   pt16_build_ms,           all_equal};

    csv.write_summary(summary, parser.stats());

    const BenchmarkSummary sassy_summary{
        baseline_results.size(), total_input_bytes,     total_baseline_phrases,
        total_sassy_phrases,     total_baseline_ms,     total_sassy_ms,
        sassy_build_ms,          sassy_all_equal};

    sassy_csv.write_summary(sassy_summary, sassy.stats());

    fs::remove_all(temporary_directory);

    // Both methods are checked against the very same baseline run above, so
    // their timings, phrase counts and hit rates below are directly
    // comparable to each other, not just to the baseline.

    const double pt16_speedup =
        total_pt16_ms == 0.0 ? 0.0 : total_baseline_ms / total_pt16_ms;
    const double sassy_speedup =
        total_sassy_ms == 0.0 ? 0.0 : total_baseline_ms / total_sassy_ms;

    const std::size_t pt16_queries = parser.stats().hits + parser.stats().misses;
    const double pt16_hit_rate =
        pt16_queries == 0 ? 0.0
                         : static_cast<double>(parser.stats().hits) /
                               static_cast<double>(pt16_queries);

    const std::size_t pt16_hit_types =
        parser.stats().singleton_hits + parser.stats().range_hits;
    const double pt16_singleton_hit_rate =
        pt16_hit_types == 0
            ? 0.0
            : static_cast<double>(parser.stats().singleton_hits) /
                  static_cast<double>(pt16_hit_types);

    const std::size_t sassy_queries =
        sassy.stats().hits + sassy.stats().misses;
    const double sassy_hit_rate =
        sassy_queries == 0 ? 0.0
                          : static_cast<double>(sassy.stats().hits) /
                                static_cast<double>(sassy_queries);

    const std::size_t sassy_hit_types =
        sassy.stats().singleton_hits + sassy.stats().range_hits;
    const double sassy_singleton_hit_rate =
        sassy_hit_types == 0
            ? 0.0
            : static_cast<double>(sassy.stats().singleton_hits) /
                  static_cast<double>(sassy_hit_types);

    // Both correctness checks are already checked separately (each method's
    // own CSV summary above); this is what a script watching only the exit
    // code needs to catch either one diverging from the baseline.
    const bool overall_equal = all_equal && sassy_all_equal;

    std::cout << "results=" << args.results << std::endl;
    std::cout << "sassy_results=" << sassy_results << std::endl;
    std::cout << "total_baseline_ms=" << total_baseline_ms << std::endl;

    std::cout << "total_pt16_ms=" << total_pt16_ms << std::endl;
    std::cout << "pt16_speedup=" << pt16_speedup << std::endl;
    std::cout << "pt16_hits=" << parser.stats().hits << std::endl;
    std::cout << "pt16_misses=" << parser.stats().misses << std::endl;
    std::cout << "pt16_hit_rate=" << pt16_hit_rate << std::endl;
    std::cout << "pt16_entries=" << parser.stats().entries << std::endl;
    std::cout << "pt16_singleton_hits=" << parser.stats().singleton_hits
              << std::endl;
    std::cout << "pt16_range_hits=" << parser.stats().range_hits << std::endl;
    std::cout << "pt16_singleton_hit_rate=" << pt16_singleton_hit_rate
              << std::endl;
    std::cout << "pt16_build_ms=" << pt16_build_ms << std::endl;
    std::cout << "pt16_MB="
              << static_cast<double>(parser.stats().approx_bytes) /
                     (1024.0 * 1024.0)
              << std::endl;
    std::cout << "pt16_outputs_equal=" << (all_equal ? "YES" : "NO")
              << std::endl;

    std::cout << "total_sassy_ms=" << total_sassy_ms << std::endl;
    std::cout << "sassy_speedup=" << sassy_speedup << std::endl;
    std::cout << "sassy_hits=" << sassy.stats().hits << std::endl;
    std::cout << "sassy_misses=" << sassy.stats().misses << std::endl;
    std::cout << "sassy_hit_rate=" << sassy_hit_rate << std::endl;
    std::cout << "sassy_entries=" << sassy.stats().entries << std::endl;
    std::cout << "sassy_sampled_entries=" << sassy.stats().sampled_entries
              << std::endl;
    std::cout << "sassy_short_suffixes=" << sassy.stats().short_suffixes
              << std::endl;
    std::cout << "sassy_singleton_hits=" << sassy.stats().singleton_hits
              << std::endl;
    std::cout << "sassy_range_hits=" << sassy.stats().range_hits << std::endl;
    std::cout << "sassy_singleton_hit_rate=" << sassy_singleton_hit_rate
              << std::endl;
    std::cout << "sassy_build_ms=" << sassy_build_ms << std::endl;
    std::cout << "sassy_MB="
              << static_cast<double>(sassy.stats().approx_bytes) /
                     (1024.0 * 1024.0)
              << std::endl;
    std::cout << "sassy_outputs_equal=" << (sassy_all_equal ? "YES" : "NO")
              << std::endl;

    std::cout << "peak_RSS_MB=" << peak_rss_mb() << std::endl;
    std::cout << "all_outputs_equal=" << (overall_equal ? "YES" : "NO")
              << std::endl;

    return overall_equal ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    if (!temporary_directory.empty()) {
      fs::remove_all(temporary_directory);
    }

    std::cerr << std::endl;
    std::cerr << "ERROR: " << error.what() << std::endl;

    return EXIT_FAILURE;
  }
}