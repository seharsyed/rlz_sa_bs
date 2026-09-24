#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ms_tools.hpp"
#include "ms_utils.hpp"
#include "ms_variants.hpp"
#include "probe_pipeline.hpp"

using msbench::SAType;
using msbench::Symbol;

/**
 * Matching-statistics benchmark.
 *
 * Two kinds of implementation are built once, then each input file is
 * processed on its own and released before the next one is loaded:
 *
 *   Full implementations (ms_variants.hpp, the baseline first) each
 *   compute complete matching statistics; the baseline's lengths are kept
 *   for the file and everything else is compared against them.
 *
 *   The PT16 table variants share one pipeline (probe_pipeline.hpp): the
 *   16-mer keys are rolled once, each probe order (bucket, sorted) is
 *   built once and every variant probes in it, and the chain runs once,
 *   on the first variant's lookup results. Only the probe is timed per
 *   variant; every other variant's lookups are checked against the
 *   chained ones, and the chain's output is checked like a full
 *   implementation's. A probe row's pipeline time (keys + its order + its
 *   probe + the chain) is what its speedup is computed from.
 *
 * Times and diagnostics accumulate over files, so the end of the run gives
 * the whole-collection totals and diagnostics; per-file rows go to the
 * results CSV, and stderr only reports a failure per file.
 *
 * Peak memory is one input plus the baseline's lengths (four bytes per
 * position) plus, for the pipeline, the keys and one order (four bytes
 * each), two lookup-result lists (the chained one and the one being
 * checked) and one MatchingStatistics.
 */

namespace {

// Percentage of `part` in `whole`, 0 if whole is 0.
double percent(const double part, const double whole) {
  return whole == 0.0 ? 0.0 : 100.0 * part / whole;
}

}  // namespace

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

    const std::string table_path = msbench::pt16_table_path(args);

    msbench::Implementations implementations =
        msbench::build_implementations(reference, suffix_array, table_path);

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

    msbench::ProberSet prober_set =
        msbench::build_probers(reference, suffix_array, table_path);
    auto& probers = prober_set.probers;

    for (const msbench::TableBuild& build : prober_set.table_builds) {
      std::cerr << "    " << build.name << ": write " << build.build_ms
                << " ms" << std::endl;
    }

    for (const auto& prober : probers) {
      std::cerr << "    " << prober->name() << ": load " << prober->build_ms()
                << " ms" << std::endl;
    }

    std::cerr << "Peak RSS after build: " << msbench::peak_rss_mb() << " MB"
              << std::endl;

    // A structural property of a built index, not of any query: printed
    // once, from the first implementation that offers one.
    for (const auto& implementation : implementations) {
      const std::string composition = implementation->indexComposition();

      if (!composition.empty()) {
        std::cerr << composition;
        break;
      }
    }

    // Pipeline totals: the shared stages, the chain's output (checked like
    // a full implementation), and one row per (variant, order).
    const bool run_pipeline = !probers.empty();

    std::vector<msbench::StageTotals> stage_totals;
    stage_totals.push_back({"keys"});
    for (const msbench::ProbeOrder order : msbench::probe_orders) {
      stage_totals.push_back({std::string(msbench::order_name(order)) +
                              "-order"});
    }
    stage_totals.push_back({"chain"});

    const std::size_t keys_stage = 0;
    const std::size_t chain_stage = stage_totals.size() - 1;

    msbench::ImplementationTotals chain_totals;
    chain_totals.name = "pt16-chain";

    // probe_totals[o * probers.size() + p]: prober p in order o.
    std::vector<msbench::ProbeTotals> probe_totals;
    for (const msbench::ProbeOrder order : msbench::probe_orders) {
      for (const auto& prober : probers) {
        msbench::ProbeTotals probe;
        probe.name = prober->name() + "-" + msbench::order_name(order);
        probe.build_ms = prober->build_ms();
        probe_totals.push_back(std::move(probe));
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
    bool stopped_early = false;

    // How the 16-mer lookups classified, summed over every file.
    SearchComposition search_totals;

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

      // Checks one complete MatchingStatistics (a full implementation's or
      // the chain's) and records it: digest, invariants, positions, lengths
      // against the baseline (or keeps them, for the baseline itself),
      // dump, CSV, checksums, totals, and a stderr line per failure. Runs
      // after timing, so checking never affects a measurement.
      const auto record_ms = [&](msbench::FileRunResult& result,
                                 const MatchingStatistics& ms,
                                 const bool compare_lengths,
                                 msbench::ImplementationTotals& row_totals) {
        result.digest = msbench::digest_ms(ms);

        if (args.check_invariants) {
          result.invariants =
              msbench::validate_invariants(ms, input, reference, alphabet);
        }

        // Gated on --no-invariants too: Validation defaults to ok, and
        // position verification is the expensive part.
        if (args.check_invariants && result.invariants.ok) {
          result.positions =
              args.verify_full
                  ? msbench::verify_all_positions(ms, input, reference,
                                                  suffix_array,
                                                  args.verify_maximality)
                  : msbench::verify_sampled(ms, input, reference,
                                            suffix_array, args.sample,
                                            args.seed, args.verify_maximality);
        }

        if (result.is_baseline) {
          baseline_lengths = msbench::extract_lengths(ms);
          baseline_min_ms = result.timing.min_ms;
        } else if (compare_lengths) {
          result.lengths = msbench::compare_lengths(baseline_lengths, ms);
        }

        if (!args.dump_directory.empty()) {
          msbench::dump_ms(args.dump_directory, file_index, filename,
                           result.implementation, ms);
        }

        result.baseline_min_ms = baseline_min_ms;

        csv.write_row(result);
        checksums.write(result);
        row_totals.accumulate(result);

        if (!result.invariants.ok) {
          std::cerr << "    " << result.implementation
                    << ": INVARIANT FAILURE: " << result.invariants
                    << std::endl;
        }

        if (!result.positions.ok) {
          std::cerr << "    " << result.implementation
                    << ": POSITION FAILURE: " << result.positions << std::endl;
        }

        if (result.lengths.checked && !result.lengths.equal) {
          std::cerr << "    " << result.implementation
                    << ": LENGTH MISMATCH: " << result.lengths.describe()
                    << std::endl;
          file_diverged = true;
        }
      };

      const auto new_result = [&](const std::string& name,
                                  const msbench::RowKind kind) {
        msbench::FileRunResult result;
        result.filename = filename;
        result.input_bytes = input.size();
        result.implementation = name;
        result.kind = kind;
        return result;
      };

      // ---------- Full implementations (the baseline first) ----------

      for (std::size_t k = 0; k < implementations.size(); ++k) {
        msbench::MSImplementation& implementation = *implementations[k];

        msbench::FileRunResult result =
            new_result(implementation.name(), msbench::RowKind::full);
        result.is_baseline = k == 0;
        result.build_ms = implementation.build_ms();

        MatchingStatistics ms;

        result.timing = msbench::time_repeated(
            args.repeats, [&] { ms = implementation.compute(input); });

        // A variant documented as not exact by design is never compared.
        record_ms(result, ms, implementation.exactExpected(), totals[k]);
        totals[k].diagnostics.accumulate(implementation.diagnostics());
      }

      // ---------- PT16 pipeline ----------

      if (run_pipeline) {
        const std::size_t n = input.size();
        const std::size_t kmer_positions = msbench::kmer_positions_of(n);

        const auto record_stage = [&](const std::size_t stage,
                                      const msbench::Timing& timing) {
          msbench::FileRunResult result =
              new_result(stage_totals[stage].name, msbench::RowKind::stage);
          result.timing = timing;
          result.baseline_min_ms = baseline_min_ms;
          csv.write_row(result);
          stage_totals[stage].accumulate(timing);
        };

        // Keys, once.
        std::vector<std::uint32_t> keys;
        const msbench::Timing keys_timing = msbench::time_repeated(
            args.repeats, [&] { msbench::roll_keys(input, keys); });
        record_stage(keys_stage, keys_timing);

        const std::uint32_t tail_key = msbench::first_tail_key(input, keys);

        // `chained` holds the first variant's first results: what the
        // chain runs on and every other probe is checked against (their
        // positions spans point into that variant's table, which lives
        // for the whole run). Every later probe writes into `results`.
        std::vector<KmerLookupResult> chained(n);
        std::vector<KmerLookupResult> results(n);
        bool have_chained = false;

        std::vector<std::uint32_t> order;
        std::vector<std::uint32_t> scratch;

        // Probe rows wait for the chain's time before they are written:
        // their pipeline time includes it.
        std::vector<std::pair<std::size_t, msbench::FileRunResult>> pending;

        for (std::size_t o = 0; o < msbench::probe_orders.size(); ++o) {
          const msbench::ProbeOrder probe_order = msbench::probe_orders[o];

          // The order, once.
          const msbench::Timing order_timing =
              msbench::time_repeated(args.repeats, [&] {
                if (probe_order == msbench::ProbeOrder::bucket) {
                  msbench::bucket_order(keys, order, scratch);
                } else {
                  msbench::sorted_order(keys, order, scratch);
                }
              });
          record_stage(1 + o, order_timing);

          // Every variant probes in it.
          for (std::size_t p = 0; p < probers.size(); ++p) {
            msbench::Prober& prober = *probers[p];
            const std::size_t row = o * probers.size() + p;

            std::vector<KmerLookupResult>& target =
                have_chained ? results : chained;

            msbench::FileRunResult result =
                new_result(probe_totals[row].name, msbench::RowKind::probe);
            result.build_ms = prober.build_ms();

            result.timing = msbench::time_repeated(args.repeats, [&] {
              prober.probe(keys, order, tail_key, target);
            });

            result.pipeline_ms = keys_timing.min_ms + order_timing.min_ms +
                                 result.timing.min_ms;

            // After timing: the chained results only need their positions
            // checked; every other probe is compared against them.
            result.lookups = msbench::compare_lookups(
                chained, target, input, reference);
            have_chained = true;

            probe_totals[row].diagnostics.accumulate(prober.diagnostics());
            pending.emplace_back(row, std::move(result));
          }

          // The order is dropped before the next one is built.
        }

        results = {};

        // The chain, once, on the chained variant's results.
        MatchingStatistics ms;
        const msbench::Timing chain_timing = msbench::time_repeated(
            args.repeats, [&] { ms = backwardChainExtend(chained); });
        record_stage(chain_stage, chain_timing);

        msbench::FileRunResult chain_result =
            new_result(chain_totals.name, msbench::RowKind::chain);
        chain_result.timing = chain_timing;
        record_ms(chain_result, ms, true, chain_totals);

        const SearchComposition composition =
            msbench::classify_lookups(chained, kmer_positions);
        search_totals.available = true;
        search_totals.singleton_hits += composition.singleton_hits;
        search_totals.range_hits += composition.range_hits;
        search_totals.misses += composition.misses;

        for (auto& [row, result] : pending) {
          result.pipeline_ms += chain_timing.min_ms;
          result.baseline_min_ms = baseline_min_ms;

          csv.write_row(result);
          probe_totals[row].accumulate(result);

          if (!result.lookups.equal) {
            std::cerr << "    " << result.implementation
                      << ": LOOKUP MISMATCH: " << result.lookups.describe()
                      << std::endl;
            file_diverged = true;
          }
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

    std::vector<msbench::ImplementationTotals> summary_totals = totals;
    if (run_pipeline) {
      summary_totals.push_back(chain_totals);
    }

    csv.write_summary(summary_totals, processed_files, total_input_bytes,
                      reference.size(), alphabet.distinct());

    const double baseline_total = totals.front().total_min_ms;
    const double megabytes =
        static_cast<double>(total_input_bytes) / (1024.0 * 1024.0);

    const auto mb_per_s = [&](const double total_ms) {
      return total_ms == 0.0 ? 0.0 : megabytes / (total_ms / 1000.0);
    };
    const auto speedup = [&](const double total_ms) {
      return total_ms == 0.0 ? 0.0 : baseline_total / total_ms;
    };
    const auto equal_count = [](std::size_t equal, std::size_t compared) {
      return std::to_string(equal) + "/" + std::to_string(compared) + " equal";
    };

    // Widths follow the longest name actually printed.
    std::size_t name_width = std::string("implementation").size();

    for (const auto& implementation : totals) {
      name_width = std::max(name_width, implementation.name.size());
    }
    for (const auto& probe : probe_totals) {
      name_width = std::max(name_width, probe.name.size());
    }
    for (const auto& stage : stage_totals) {
      name_width = std::max(name_width, stage.name.size());
    }

    name_width += 2;
    const int name_w = static_cast<int>(name_width);

    std::cerr << std::left << std::setw(name_w) << "implementation"
              << std::right << std::setw(12) << "build ms" << std::setw(14)
              << "total ms" << std::setw(10) << "MB/s" << std::setw(10)
              << "speedup" << std::setw(16) << "lengths" << std::endl;

    for (const auto& implementation : totals) {
      std::cerr << std::left << std::setw(name_w) << implementation.name
                << std::right << std::fixed << std::setprecision(2)
                << std::setw(12) << implementation.build_ms << std::setw(14)
                << implementation.total_min_ms << std::setw(10)
                << mb_per_s(implementation.total_min_ms) << std::setw(10)
                << speedup(implementation.total_min_ms) << std::setw(16);

      if (implementation.is_baseline) {
        std::cerr << "baseline";
      } else if (implementation.files_compared == 0) {
        std::cerr << "not compared";
      } else {
        std::cerr << equal_count(implementation.files_lengths_equal,
                                 implementation.files_compared);
      }

      std::cerr << std::endl;
    }

    if (run_pipeline) {
      double shared_total = 0.0;
      for (const auto& stage : stage_totals) {
        shared_total += stage.total_min_ms;
      }

      std::cerr << std::endl
                << "PT16 pipeline: shared stages, run once per file"
                << std::endl;
      std::cerr << std::left << std::setw(name_w) << "stage" << std::right
                << std::setw(14) << "total ms" << std::setw(10) << "share"
                << std::endl;

      for (std::size_t s = 0; s < stage_totals.size(); ++s) {
        const auto& stage = stage_totals[s];

        std::cerr << std::left << std::setw(name_w) << stage.name
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(14) << stage.total_min_ms << std::setw(9)
                  << std::setprecision(1)
                  << percent(stage.total_min_ms, shared_total) << "%";

        if (s == chain_stage) {
          std::cerr << "   lengths "
                    << equal_count(chain_totals.files_lengths_equal,
                                   chain_totals.files_compared);
        }

        std::cerr << std::endl;
      }

      std::cerr << std::endl
                << "PT16 pipeline: per variant and order (pipeline = keys + "
                   "order + probe + chain)"
                << std::endl;
      std::cerr << std::left << std::setw(name_w) << "variant" << std::right
                << std::setw(12) << "load ms" << std::setw(14) << "probe ms"
                << std::setw(15) << "pipeline ms" << std::setw(10) << "MB/s"
                << std::setw(10) << "speedup" << std::setw(16) << "lookups"
                << std::endl;

      for (const auto& probe : probe_totals) {
        std::cerr << std::left << std::setw(name_w) << probe.name
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(12) << probe.build_ms << std::setw(14)
                  << probe.total_min_ms << std::setw(15)
                  << probe.total_pipeline_ms << std::setw(10)
                  << mb_per_s(probe.total_pipeline_ms) << std::setw(10)
                  << speedup(probe.total_pipeline_ms) << std::setw(16)
                  << equal_count(probe.files_lookups_equal,
                                 probe.files_compared)
                  << std::endl;
      }

      std::cerr << "table files:";
      for (const auto& build : prober_set.table_builds) {
        std::cerr << "  " << build.name << " " << std::setprecision(2)
                  << build.build_ms << " ms";
      }
      std::cerr << std::endl;
    }

    // ---------- Diagnostics, summed over files ----------

    std::cerr << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "[10] DIAGNOSTICS (summed over " << processed_files
              << " files)" << std::endl;
    std::cerr << "========================================" << std::endl;

    for (const auto& implementation : totals) {
      std::cerr << implementation.name << std::endl;

      if (implementation.is_baseline) {
        // The true average phrase length (total length / input length):
        // the baseline computes complete matching statistics.
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

    for (const auto& probe : probe_totals) {
      if (!probe.diagnostics.empty()) {
        std::cerr << probe.name << std::endl << probe.diagnostics.format();
      }
    }

    // How the 16-mer lookups classified, over every file: the same for
    // every correct variant, so printed once.
    if (search_totals.available) {
      const double total_queries =
          static_cast<double>(search_totals.singleton_hits +
                              search_totals.range_hits + search_totals.misses);

      std::cerr << "16-mer lookups: " << std::fixed << std::setprecision(1)
                << "singleton="
                << percent(static_cast<double>(search_totals.singleton_hits),
                           total_queries)
                << "%  range="
                << percent(static_cast<double>(search_totals.range_hits),
                           total_queries)
                << "%  miss="
                << percent(static_cast<double>(search_totals.misses),
                           total_queries)
                << "%" << std::endl;
    }

    // ---------- Machine-readable summary ----------

    bool all_ok = !stopped_early;

    const auto print_implementation = [&](const auto& implementation) {
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
    };

    for (const auto& implementation : totals) {
      print_implementation(implementation);
      std::cout << implementation.name
                << ".MB_per_s=" << mb_per_s(implementation.total_min_ms)
                << std::endl;
      std::cout << implementation.name << ".speedup_vs_baseline="
                << speedup(implementation.total_min_ms) << std::endl;
    }

    if (run_pipeline) {
      // The chain row's time is the chain stage alone.
      print_implementation(chain_totals);

      for (const auto& stage : stage_totals) {
        std::cout << "stage." << stage.name
                  << ".total_min_ms=" << stage.total_min_ms << std::endl;
        csv.write_summary_value("stage." + stage.name, "total_min_ms",
                                stage.total_min_ms);
      }

      for (const auto& probe : probe_totals) {
        const std::string prefix = probe.name + ".";

        std::cout << prefix << "load_ms=" << probe.build_ms << std::endl;
        std::cout << prefix << "probe_min_ms=" << probe.total_min_ms
                  << std::endl;
        std::cout << prefix << "probe_first_ms=" << probe.total_first_ms
                  << std::endl;
        std::cout << prefix << "pipeline_ms=" << probe.total_pipeline_ms
                  << std::endl;
        std::cout << prefix
                  << "MB_per_s=" << mb_per_s(probe.total_pipeline_ms)
                  << std::endl;
        std::cout << prefix << "speedup_vs_baseline="
                  << speedup(probe.total_pipeline_ms) << std::endl;
        std::cout << prefix
                  << "files_lookups_equal=" << probe.files_lookups_equal << "/"
                  << probe.files_compared << std::endl;

        csv.write_summary_value(probe.name, "probe_min_ms",
                                probe.total_min_ms);
        csv.write_summary_value(probe.name, "pipeline_ms",
                                probe.total_pipeline_ms);
        csv.write_summary_value(probe.name, "speedup_vs_baseline",
                                speedup(probe.total_pipeline_ms));
        csv.write_summary_value(
            probe.name, "files_lookups_equal",
            std::to_string(probe.files_lookups_equal) + "/" +
                std::to_string(probe.files_compared));

        all_ok = all_ok && probe.ok();
      }

      for (const auto& build : prober_set.table_builds) {
        std::cout << "table." << build.name.substr(0, build.name.find(' '))
                  << ".write_ms=" << build.build_ms << std::endl;
      }
    }

    std::cout << "results=" << args.results << std::endl;
    std::cout << "checksums=" << args.checksums << std::endl;
    std::cout << "implementations=" << implementations.size() << std::endl;
    std::cout << "probers=" << probers.size() << std::endl;
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
