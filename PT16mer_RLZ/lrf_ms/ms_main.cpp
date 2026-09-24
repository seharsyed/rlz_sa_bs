#include <algorithm>
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
 * and diagnostics accumulate per implementation, so the end of the run
 * gives the whole-collection totals and diagnostics for each variant
 * (per-file rows still go to the results CSV); stderr only reports a
 * failure per file.
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

    // A structural property of the built table (e.g. singleton vs. range
    // entry counts), not of any query, so it does not vary between
    // implementations sharing the same table format -- print only the
    // first non-empty one, once, rather than once per implementation.
    for (const auto& implementation : implementations) {
      const std::string composition = implementation->indexComposition();

      if (!composition.empty()) {
        std::cerr << composition;
        break;
      }
    }

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

    // The raw-search classification (see search_composition below), summed
    // over every file for the summary.
    SearchComposition search_totals;
    bool stopped_early = false;

    for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
      const std::string& filename = files[file_index];

      // Loading is deliberately outside every timed region.
      const auto input = msbench::load_input<Symbol>(filename);

      std::cerr << "[" << file_index + 1 << "/" << files.size() << "] "
                << filename << "  (" << input.size() << " bytes)" << std::endl;

      // The baseline's lengths for THIS file only, released with the
      // file. Four bytes per position rather than sixteen.
      std::vector<std::uint32_t> baseline_lengths;
      double baseline_min_ms = 0.0;
      bool file_diverged = false;

      // The raw-search classification is identical across every scan-only
      // PT16 variant for this file (same table content, same queries), so
      // only the first one that offers it is kept, and added once to
      // search_totals after the loop below -- not once per implementation.
      // Only the classification itself, not a length: the scan-only
      // floor's own average length is not the true average phrase length
      // (every hit is capped at 16 there, whatever the true match reaches),
      // so it is not reported at all -- see the baseline's own
      // avg_phrase_length in the summary for that.
      bool search_composition_captured = false;
      SearchComposition search_composition;  // global namespace (pt16_utils.hpp)

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

          // Gated the same as validate_invariants above (not just on its
          // result): Validation defaults to .ok = true, so with
          // --no-invariants skipping the check above, result.invariants
          // was still reading as "ok" here and this ran anyway -- for
          // every position (or up to --sample of them), a binary search
          // over the whole suffix array. That made --no-invariants not
          // actually skip the expensive part.
          if (args.check_invariants && result.invariants.ok) {
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

          if (is_baseline) {
            baseline_lengths = msbench::extract_lengths(ms);
            baseline_min_ms = result.timing.min_ms;
          } else if (implementation.exactExpected()) {
            // Skipped entirely for a variant documented as NOT exact by
            // design (a scan-only floor, or one-step chainExtend): it is
            // never going to equal the baseline, so comparing it was
            // never a correctness question, only a per-file, per-run
            // O(n) scan (and a LENGTH MISMATCH report) that cost time
            // without telling us anything the doc comment doesn't
            // already say.
            result.lengths = msbench::compare_lengths(baseline_lengths, ms);
          }

          if (!args.dump_directory.empty()) {
            msbench::dump_ms(args.dump_directory, file_index, filename,
                             implementation.name(), ms);
          }

          // ms is released here, before the next implementation runs.
        }

        if (!search_composition_captured) {
          const SearchComposition composition = implementation.searchComposition();

          if (composition.available) {
            search_composition_captured = true;
            search_composition = composition;
          }
        }

        result.baseline_min_ms = baseline_min_ms;

        csv.write_row(result);
        checksums.write(result);
        totals[k].accumulate(result);

        // Summed over every file and printed once per implementation in
        // the summary at the end, rather than per file.
        totals[k].diagnostics.accumulate(implementation.diagnostics());

        if (!result.invariants.ok) {
          std::cerr << "    " << implementation.name()
                    << ": INVARIANT FAILURE: " << result.invariants
                    << std::endl;
        }

        if (!result.positions.ok) {
          std::cerr << "    " << implementation.name()
                    << ": POSITION FAILURE: " << result.positions
                    << std::endl;
        }

        if (result.lengths.checked && !result.lengths.equal) {
          std::cerr << "    " << implementation.name()
                    << ": LENGTH MISMATCH: " << result.lengths.describe()
                    << std::endl;
          file_diverged = true;
        }
      }

      if (search_composition_captured) {
        search_totals.available = true;
        search_totals.singleton_hits += search_composition.singleton_hits;
        search_totals.range_hits += search_composition.range_hits;
        search_totals.misses += search_composition.misses;
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

    // Implementation names vary a lot in length (e.g. "lrf-ms" vs.
    // "pt16-v2-bucket-chain-multi"); a fixed width overflows for the
    // longer ones and breaks every column after it. Width it to the
    // longest name actually registered instead, with a little breathing
    // room before the next column.
    std::size_t name_width = std::string("implementation").size();

    for (const msbench::ImplementationTotals& implementation : totals) {
      name_width = std::max(name_width, implementation.name.size());
    }

    name_width += 2;

    std::cerr << std::left << std::setw(static_cast<int>(name_width))
              << "implementation" << std::right << std::setw(12)
              << "build ms" << std::setw(14) << "total ms" << std::setw(10)
              << "MB/s" << std::setw(10) << "speedup" << std::setw(16)
              << "lengths" << std::endl;

    for (const msbench::ImplementationTotals& implementation : totals) {
      const double seconds = implementation.total_min_ms / 1000.0;

      std::cerr << std::left << std::setw(static_cast<int>(name_width))
                << implementation.name << std::right << std::fixed
                << std::setprecision(2)
                << std::setw(12) << implementation.build_ms << std::setw(14)
                << implementation.total_min_ms << std::setw(10)
                << (seconds == 0.0 ? 0.0 : megabytes / seconds) << std::setw(10)
                << (implementation.total_min_ms == 0.0
                        ? 0.0
                        : baseline_total / implementation.total_min_ms)
                << std::setw(16);

      if (implementation.is_baseline) {
        std::cerr << "baseline";
      } else if (implementation.files_compared == 0) {
        // Never checked: a variant documented as not exact by design
        // (see MSImplementation::exactExpected), not a file count of 0.
        std::cerr << "not compared";
      } else {
        std::cerr << (std::to_string(implementation.files_lengths_equal) + "/" +
                      std::to_string(implementation.files_compared) + " equal");
      }

      std::cerr << std::endl;
    }

    // ---------- Per-implementation diagnostics, summed over files ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[10] DIAGNOSTICS (summed over " << processed_files
              << " files)" << std::endl;
    std::cerr << "========================================" << std::endl;

    for (const msbench::ImplementationTotals& implementation : totals) {
      std::cerr << implementation.name << std::endl;

      if (implementation.is_baseline) {
        // The TRUE average phrase length (total length / input length):
        // the baseline computes complete, unextended matching statistics,
        // unlike the scan-only PT16 variants, whose hits are capped at 16.
        const double avg_phrase_length =
            total_input_bytes == 0
                ? 0.0
                : static_cast<double>(implementation.total_len) /
                      static_cast<double>(total_input_bytes);

        std::cerr << "        total_len=" << implementation.total_len
                  << " max_len=" << implementation.max_len
                  << " zero_len=" << implementation.zero_len_count
                  << " avg_phrase_length=" << std::fixed
                  << std::setprecision(3) << avg_phrase_length << std::endl;
      }

      std::cerr << implementation.diagnostics.format();
    }

    // How the raw search's lookups classified, over every file. Identical
    // across the scan-only PT16 variants (same table content, same
    // queries), so printed once. No length here -- a raw-search hit is
    // always reported at exactly 16 regardless of how far the true match
    // extends; see the baseline's avg_phrase_length above for that.
    if (search_totals.available) {
      const std::size_t total_queries = search_totals.singleton_hits +
                                        search_totals.range_hits +
                                        search_totals.misses;
      const auto percent = [&](const std::size_t count) {
        return total_queries == 0 ? 0.0
                                  : 100.0 * static_cast<double>(count) /
                                        static_cast<double>(total_queries);
      };

      std::cerr << "search (shared across pt16 variants): " << std::fixed
                << std::setprecision(1)
                << "singleton=" << percent(search_totals.singleton_hits)
                << "%  range=" << percent(search_totals.range_hits)
                << "%  short=" << percent(search_totals.misses) << "%"
                << std::endl;
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