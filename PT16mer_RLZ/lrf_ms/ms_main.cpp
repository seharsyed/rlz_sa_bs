#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ms_utils.hpp"
#include "ms_variants.hpp"

using msbench::SAType;
using msbench::Symbol;

/**
 * Matching-statistics benchmark.
 *
 * Every implementation registered in ms_variants.hpp is built once, then
 * each input file is processed on its own: the baseline runs, its
 * lengths are kept, every variant runs and is compared against those
 * lengths, and the file is released before the next one is loaded. Times
 * accumulate per implementation, so the end of the run gives both the
 * per-file detail and the whole-collection total for each variant.
 *
 * Peak memory is therefore one input plus the baseline's lengths (four
 * bytes per position) plus one variant's MatchingStatistics, rather than
 * every implementation's output for every file at once.
 */

int main(int argc, char** argv) {
  try {
    std::cerr << "========================================" << std::endl;
    std::cerr << "Matching statistics benchmark" << std::endl;
    std::cerr << "========================================" << std::endl;

    // ---------- Load experiment inputs ----------

    std::cerr << "[1] Reading arguments..." << std::endl;
    const msbench::Args args = msbench::parse_args(argc, argv);

    std::cerr << "[2] Loading reference..." << std::endl;

    std::vector<Symbol> reference;
    const double reference_load_ms = msbench::time_ms(
        [&] { reference = msbench::load_reference<Symbol>(args.reference); });

    std::cerr << "    Reference loaded: " << reference.size() << " bytes in "
              << reference_load_ms << " ms" << std::endl;

    std::cerr << "[3] Loading suffix array..." << std::endl;

    std::vector<SAType> suffix_array;
    const double suffix_array_load_ms = msbench::time_ms([&] {
      suffix_array = msbench::load_suffix_array<SAType>(args.suffix_array);
    });

    std::cerr << "    Suffix array loaded: " << suffix_array.size()
              << " entries in " << suffix_array_load_ms << " ms" << std::endl;

    if (suffix_array.size() != reference.size()) {
      throw std::runtime_error(
          "suffix array has " + std::to_string(suffix_array.size()) +
          " entries but the reference has " + std::to_string(reference.size()) +
          " symbols; the suffix array must be built over the reference "
          "WITHOUT a terminating sentinel");
    }

    std::cerr << "[4] Loading input list..." << std::endl;
    const auto files = msbench::load_input_list(args.filenames);
    std::cerr << "    Input files: " << files.size() << std::endl;

    std::cerr << "[5] Scanning reference alphabet..." << std::endl;
    const msbench::SymbolTable<Symbol> alphabet(reference);
    std::cerr << "    Distinct symbols: " << alphabet.distinct() << std::endl;

    // ---------- Build ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[6] BUILD" << std::endl;
    std::cerr << "========================================" << std::endl;

    msbench::Implementations implementations =
        msbench::build_implementations(reference, suffix_array,
                                       msbench::pt16_table_path(args));

    if (implementations.empty()) {
      throw std::runtime_error("no implementations registered");
    }

    std::vector<msbench::ImplementationTotals> totals(implementations.size());

    for (std::size_t k = 0; k < implementations.size(); ++k) {
      totals[k].name = implementations[k]->name();
      totals[k].is_baseline = k == 0;
      totals[k].build_ms = implementations[k]->build_ms();

      std::cerr << "    " << implementations[k]->name()
                << (k == 0 ? "  (baseline)" : "") << ": build "
                << implementations[k]->build_ms() << " ms" << std::endl;
    }

    std::cerr << "Peak RSS after build: " << msbench::peak_rss_mb() << " MB"
              << std::endl;

    // ---------- Output files ----------

    std::cerr << std::endl;
    std::cerr << "[7] Opening results files..." << std::endl;

    msbench::CSVWriter csv(args.results);
    msbench::ChecksumWriter checksums(args.checksums);

    std::cerr << "    results:   " << args.results << std::endl;
    std::cerr << "    checksums: " << args.checksums << std::endl;

    if (!args.dump_directory.empty()) {
      std::cerr << "    dumps:     " << args.dump_directory << std::endl;
    }

    // ---------- Per-file runs ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[8] MATCHING STATISTICS" << std::endl;
    std::cerr << "========================================" << std::endl;

    std::size_t processed_files = 0;
    std::size_t total_input_bytes = 0;
    bool stopped_early = false;

    for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
      const std::string& filename = files[file_index];

      // Loading is deliberately outside every timed region.
      const auto input = msbench::load_input<Symbol>(filename);

      std::cerr << "\n[" << file_index + 1 << "/" << files.size() << "] "
                << filename << "  (" << input.size() << " bytes)" << std::endl;

      // The baseline's lengths for THIS file only, released with the
      // file. Four bytes per position rather than sixteen.
      std::vector<std::uint32_t> baseline_lengths;
      double baseline_min_ms = 0.0;
      bool file_diverged = false;

      for (std::size_t k = 0; k < implementations.size(); ++k) {
        msbench::MSImplementation& implementation = *implementations[k];
        const bool is_baseline = k == 0;

        msbench::FileRunResult result;
        result.filename = filename;
        result.input_bytes = input.size();
        result.implementation = implementation.name();
        result.is_baseline = is_baseline;
        result.build_ms = implementation.build_ms();

        {
          MatchingStatistics ms;

          result.timing = msbench::time_repeated(
              args.repeats, [&] { ms = implementation.compute(input); });

          // Everything below runs after timing, so checking never
          // affects the measurement.

          result.digest = msbench::digest_ms(ms);

          if (args.check_invariants) {
            result.invariants =
                msbench::validate_invariants(ms, input, reference, alphabet);
          }

          if (result.invariants.ok) {
            result.positions =
                args.verify_full
                    ? msbench::verify_all_positions(ms, input, reference,
                                                    suffix_array,
                                                    args.verify_maximality)
                    : msbench::verify_sampled(ms, input, reference,
                                              suffix_array, args.sample,
                                              args.seed,
                                              args.verify_maximality);
          }

          if (args.verify_brute && result.invariants.ok) {
            if (input.size() <= args.brute_limit &&
                reference.size() <= args.brute_limit) {
              const msbench::Validation brute =
                  msbench::verify_against_brute_force(ms, input, reference);

              if (!brute.ok) {
                result.invariants = brute;
              }
            } else if (file_index == 0 && k == 0) {
              std::cerr << "    (brute-force check skipped: over "
                           "--brute-limit)"
                        << std::endl;
            }
          }

          if (is_baseline) {
            baseline_lengths = msbench::extract_lengths(ms);
            baseline_min_ms = result.timing.min_ms;
          } else {
            result.lengths = msbench::compare_lengths(baseline_lengths, ms);
          }

          if (!args.dump_directory.empty()) {
            msbench::dump_ms(args.dump_directory, file_index, filename,
                             implementation.name(), ms);
          }

          // ms is released here, before the next implementation runs.
        }

        result.baseline_min_ms = baseline_min_ms;

        csv.write_row(result);
        checksums.write(result);
        totals[k].accumulate(result);

        // ---------- Per-implementation report for this file ----------

        std::cerr << "    " << std::left << std::setw(16)
                  << implementation.name() << std::right << std::fixed
                  << std::setprecision(3) << std::setw(10)
                  << result.timing.min_ms << " ms";

        if (!is_baseline) {
          std::cerr << "  speedup " << std::setprecision(2) << result.speedup()
                    << "x";
        }

        std::cerr << "  len_hash=" << msbench::to_hex(result.digest.len_hash);

        if (is_baseline) {
          std::cerr << "  total_len=" << result.digest.total_len
                    << " max_len=" << result.digest.max_len
                    << " zero_len=" << result.digest.zero_len_count;
        } else {
          std::cerr << "  lengths=" << (result.lengths.equal ? "EQUAL" : "DIFFER");
        }

        std::cerr << std::endl;

        if (!result.invariants.ok) {
          std::cerr << "        INVARIANT FAILURE: " << result.invariants
                    << std::endl;
        }

        if (!result.positions.ok) {
          std::cerr << "        POSITION FAILURE: " << result.positions
                    << std::endl;
        }

        if (result.lengths.checked && !result.lengths.equal) {
          std::cerr << "        LENGTH MISMATCH: " << result.lengths.describe()
                    << std::endl;
          file_diverged = true;
        }
      }

      ++processed_files;
      total_input_bytes += input.size();

      if (file_diverged && args.stop_on_mismatch) {
        std::cerr << "\nStopping: --stop-on-mismatch and this file diverged."
                  << std::endl;
        stopped_early = true;
        break;
      }
    }

    // ---------- Collection totals ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[9] COLLECTION TOTALS" << std::endl;
    std::cerr << "========================================" << std::endl;

    csv.write_summary(totals, processed_files, total_input_bytes,
                      reference.size(), alphabet.distinct());

    const double baseline_total = totals.front().total_min_ms;
    const double megabytes =
        static_cast<double>(total_input_bytes) / (1024.0 * 1024.0);

    std::cerr << std::left << std::setw(16) << "implementation" << std::right
              << std::setw(12) << "build ms" << std::setw(14) << "total ms"
              << std::setw(10) << "MB/s" << std::setw(10) << "speedup"
              << std::setw(16) << "lengths" << std::endl;

    for (const msbench::ImplementationTotals& implementation : totals) {
      const double seconds = implementation.total_min_ms / 1000.0;

      std::cerr << std::left << std::setw(16) << implementation.name
                << std::right << std::fixed << std::setprecision(2)
                << std::setw(12) << implementation.build_ms << std::setw(14)
                << implementation.total_min_ms << std::setw(10)
                << (seconds == 0.0 ? 0.0 : megabytes / seconds) << std::setw(10)
                << (implementation.total_min_ms == 0.0
                        ? 0.0
                        : baseline_total / implementation.total_min_ms)
                << std::setw(16);

      if (implementation.is_baseline) {
        std::cerr << "baseline";
      } else {
        std::cerr << (std::to_string(implementation.files_lengths_equal) + "/" +
                      std::to_string(implementation.files_compared) + " equal");
      }

      std::cerr << std::endl;
    }

    // ---------- Machine-readable summary ----------

    bool all_ok = !stopped_early;

    for (const msbench::ImplementationTotals& implementation : totals) {
      const double seconds = implementation.total_min_ms / 1000.0;
      const std::string prefix = implementation.name + ".";

      std::cout << prefix << "build_ms=" << implementation.build_ms
                << std::endl;
      std::cout << prefix << "total_min_ms=" << implementation.total_min_ms
                << std::endl;
      std::cout << prefix << "total_first_ms=" << implementation.total_first_ms
                << std::endl;
      std::cout << prefix << "total_entries=" << implementation.total_entries
                << std::endl;
      std::cout << prefix
                << "MB_per_s=" << (seconds == 0.0 ? 0.0 : megabytes / seconds)
                << std::endl;
      std::cout << prefix << "speedup_vs_baseline="
                << (implementation.total_min_ms == 0.0
                        ? 0.0
                        : baseline_total / implementation.total_min_ms)
                << std::endl;
      std::cout << prefix
                << "files_lengths_equal=" << implementation.files_lengths_equal
                << "/" << implementation.files_compared << std::endl;
      std::cout << prefix << "all_lengths_equal="
                << (implementation.is_baseline
                        ? "NA"
                        : (implementation.all_lengths_equal() ? "YES" : "NO"))
                << std::endl;
      std::cout << prefix
                << "invariant_failures=" << implementation.invariant_failures
                << std::endl;
      std::cout << prefix
                << "position_failures=" << implementation.position_failures
                << std::endl;

      all_ok = all_ok && implementation.ok();
    }

    std::cout << "results=" << args.results << std::endl;
    std::cout << "checksums=" << args.checksums << std::endl;
    std::cout << "implementations=" << implementations.size() << std::endl;
    std::cout << "processed_files=" << processed_files << std::endl;
    std::cout << "total_input_bytes=" << total_input_bytes << std::endl;
    std::cout << "reference_bytes=" << reference.size() << std::endl;
    std::cout << "peak_RSS_MB=" << msbench::peak_rss_mb() << std::endl;
    std::cout << "all_ok=" << (all_ok ? "YES" : "NO") << std::endl;

    return all_ok ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << std::endl;
    std::cerr << "ERROR: " << error.what() << std::endl;

    return EXIT_FAILURE;
  }
}