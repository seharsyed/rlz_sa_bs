// ms2 benchmark: the new matching-statistics methods, built up step by
// step. Separate from lrf_ms/ms_main.cpp on purpose.
//
//   g++ -std=c++20 -O3 ms2/ms2_main.cpp -o ms2_main
//   ./ms2_main --reference REF --suffix-array REF.sa --filenames LIST
//              [--results CSV] [--table PATH] [--max-files N]
//              [--max-bytes N] [--no-check] [--no-chain]
//              [--dump-lookups PATH] [--quiet] [--sample-step L]
//
// Built once: lrf-ms, and the PT16 tables (v2, sassy; each also searched
// with a finger on the sorted order). Then per input file:
//
//   lrf-ms         the baseline: full matching statistics (lrf_ms.hpp);
//   keys           one pass over the input (input_keys.hpp): the packed
//                  (key, position) values to query, with all-A / all-T
//                  16-mers left out and recorded as runs, the short
//                  queries (positions before a separator or the end) and
//                  the separator count;
//   bucket-order   the keys ordered by bucket (high 16 bits of the key);
//   sorted-order   the keys ordered fully by key;
//   results        the lookup-results arrays made ready for this file:
//                  sized to the input, separator and run positions
//                  emptied (every lookup writes the other positions);
//   <variant>-<order>
//                  the table lookups over the keys in bucket order and in
//                  sorted order (the finger variant: sorted only), into an
//                  array of lookup results by text position
//                  (probe_order.hpp); every row's results are checked
//                  against the first row's;
//   chain          the backward chain (chain.hpp, the same algorithm as
//                  backwardChainExtend) on the first row's results,
//                  checked against lrf-ms. Run positions
//                  have no lookup yet, so mismatches caused by them are
//                  counted apart from real ones. --no-chain skips it.
//
// After that, per file and table, the sampled method (sampled.hpp): one
// lookup at every L-th position in text order (default L = 16), a
// right-to-left pass linking neighbouring samples into stretches, and one
// midpoint lookup per hole -- each step timed, the stretches checked
// against the reference and lrf-ms.
//
// Last, per file and table, the sampled method end to end (sampled_ms.hpp):
// samples (sorted, finger lookups) linked into stretches, the remaining
// hole positions' keys rolled, sorted and looked up, then one right-to-left
// pass giving the full matching statistics -- checked against lrf-ms at
// every position.
//
// Each file's lines print as soon as they are measured, so progress on a
// large input is visible; totals follow at the end. Each order is checked
// after it is timed (unless --no-check).

#include <algorithm>
#include <cstdint>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "../lrf_ms/lrf_ms.hpp"
#include "../lrf_ms/ms_utils.hpp"  // loaders, SymbolTable
#include "input_keys.hpp"
#include "phrase_stats.hpp"
#include "sampled.hpp"
#include "sampled_ms.hpp"
#include "chain.hpp"
#include "probe_order.hpp"

using Symbol = unsigned char;
using SAType = std::uint32_t;

namespace {

struct Ms2Args {
  std::string reference;
  std::string suffix_array;
  std::string filenames;
  std::string results;  // optional CSV
  std::string table;    // PT16 table files; default beside the reference
  std::size_t max_files = 0;  // 0: all
  std::size_t max_bytes = 0;  // 0: whole file; otherwise a prefix
  bool check = true;
  bool chain = true;
  std::string dump_lookups;  // optional: the first file's lookup table
  bool quiet = false;        // no per-file lines, only totals
  std::size_t sample_step = 16;  // the sampled method's L
  bool keep_runs = false;    // all-A / all-T 16-mers as ordinary keys
};

void print_usage(const char* program) {
  std::cout << "Usage: " << program
            << " --reference PATH --suffix-array PATH --filenames PATH\n"
            << "  [--results PATH]    per-file CSV\n"
            << "  [--table PATH]      PT16 table files (default: "
               "<reference>.ms2_pt16, and .sassy)\n"
            << "  [--max-files N]     only the first N input files\n"
            << "  [--max-bytes N]     only the first N bytes of each input\n"
            << "  [--no-check]        skip checking the orders\n"
            << "  [--no-chain]        skip the chain (and its check)\n"
            << "  [--quiet]           no per-file lines, only the totals\n"
            << "  [--keep-runs]       look up all-A / all-T 16-mers like any "
               "other key\n"
            << "                      (default: set them aside as runs)\n"
            << "  [--sample-step L]   the sampled method's step (1..16, "
               "default 16)\n"
            << "  [--dump-lookups PATH]  write the first file's lookup table:\n"
            << "                      one line per input position, the size of\n"
            << "                      its 16-mer's SA interval (0 empty, 1\n"
            << "                      singleton, else the range size)\n";
}

Ms2Args parse_ms2_args(int argc, char** argv) {
  Ms2Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];

    if (option == "--reference") {
      args.reference = msbench::require_value(i, argc, argv);
    } else if (option == "--suffix-array") {
      args.suffix_array = msbench::require_value(i, argc, argv);
    } else if (option == "--filenames") {
      args.filenames = msbench::require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = msbench::require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = msbench::parse_number(
          option, msbench::require_value(i, argc, argv));
    } else if (option == "--max-bytes") {
      args.max_bytes = msbench::parse_number(
          option, msbench::require_value(i, argc, argv));
    } else if (option == "--table") {
      args.table = msbench::require_value(i, argc, argv);
    } else if (option == "--no-check") {
      args.check = false;
    } else if (option == "--no-chain") {
      args.chain = false;
    } else if (option == "--quiet") {
      args.quiet = true;
    } else if (option == "--keep-runs") {
      args.keep_runs = true;
    } else if (option == "--sample-step") {
      args.sample_step = msbench::parse_number(
          option, msbench::require_value(i, argc, argv));
      if (args.sample_step < 1 || args.sample_step > 16) {
        throw std::runtime_error("--sample-step must be 1..16");
      }
    } else if (option == "--dump-lookups") {
      args.dump_lookups = msbench::require_value(i, argc, argv);
    } else if (option == "--help" || option == "-h") {
      print_usage(argv[0]);
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }

  if (args.reference.empty() || args.suffix_array.empty() ||
      args.filenames.empty()) {
    throw std::runtime_error("reference, suffix-array and filenames are required");
  }

  if (args.table.empty()) {
    args.table = args.reference + ".ms2_pt16";
  }

  return args;
}

// ---------- Checks (untimed) ----------

// `order` holds exactly the values of `keys` (same count, same wrapping
// sum and xor), and every neighbouring pair is in `before` order.
template <typename Before>
bool is_ordering_of(const std::vector<std::uint64_t>& keys,
                    const std::vector<std::uint64_t>& order,
                    const Before& before) {
  if (keys.size() != order.size()) return false;

  std::uint64_t sum_keys = 0, sum_order = 0, xor_keys = 0, xor_order = 0;
  for (const std::uint64_t v : keys) {
    sum_keys += v;
    xor_keys ^= v;
  }
  for (const std::uint64_t v : order) {
    sum_order += v;
    xor_order ^= v;
  }
  if (sum_keys != sum_order || xor_keys != xor_order) return false;

  for (std::size_t i = 1; i < order.size(); ++i) {
    if (!before(order[i - 1], order[i])) return false;
  }

  return true;
}

// Bucket order: buckets non-decreasing, positions increasing within one.
bool bucket_before(const std::uint64_t a, const std::uint64_t b) {
  const std::uint64_t bucket_a = a >> 48;
  const std::uint64_t bucket_b = b >> 48;
  return bucket_a < bucket_b ||
         (bucket_a == bucket_b && ms2::position_of(a) < ms2::position_of(b));
}

// Sorted order: keys non-decreasing, positions increasing within one --
// which for packed values is just strictly increasing.
bool sorted_before(const std::uint64_t a, const std::uint64_t b) {
  return a < b;
}

// ---------- Lookup table dump ----------

// One line per input position: 0 for an empty entry (a miss, a short
// query, a separator or a run position -- found is false), 1 for a
// singleton, otherwise the range's occurrence count.
void dump_lookup_table(const std::string& path,
                       const std::vector<KmerLookupResult>& results) {
  std::ofstream output(path, std::ios::binary);
  if (!output) throw std::runtime_error("cannot create " + path);

  std::string buffer;
  buffer.reserve(1 << 20);
  char digits[16];

  for (const KmerLookupResult& result : results) {
    const std::uint32_t size = result.found ? result.count : 0;
    const auto [end, error] = std::to_chars(digits, digits + sizeof(digits), size);
    buffer.append(digits, end);
    buffer.push_back('\n');

    if (buffer.size() >= (1 << 20) - 32) {
      output.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      buffer.clear();
    }
  }

  output.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  if (!output) throw std::runtime_error("failed writing " + path);
}

// ---------- Lookup variants ----------

// One table and one way of searching it.
struct Variant {
  std::string name;
  bool sorted_only = false;  // e.g. a finger search: needs sorted keys
  std::function<void(const ms2::KeyedInput&, const std::vector<std::uint64_t>&,
                     std::vector<KmerLookupResult>&, Diagnostics*)>
      probe;
};

template <typename Policy>
Variant make_variant(std::string name,
                     std::shared_ptr<const typename Policy::Table> table,
                     const bool sorted_only = false) {
  Variant variant;
  variant.name = std::move(name);
  variant.sorted_only = sorted_only;
  variant.probe = [table](const ms2::KeyedInput& keyed,
                          const std::vector<std::uint64_t>& order,
                          std::vector<KmerLookupResult>& results,
                          Diagnostics* diagnostics) {
    ms2::probe_order<Policy>(*table, keyed, order, results, diagnostics);
  };
  return variant;
}

// One variant run on one order: a benchmark row, "<variant>-<order>".
struct Row {
  std::string name;
  std::size_t variant = 0;
  bool sorted = false;  // false: bucket order

  double total_ms = 0.0;
  Diagnostics diagnostics;  // summed over files
  std::size_t files_compared = 0;
  std::size_t files_equal = 0;
};

// ---------- The sampled method ----------

std::string percent(std::size_t part, std::size_t whole);  // below

// "stretches 92.96%, holes 1.35% (mean hole 9.1 positions), separators
// 5.69% of the input".
std::string coverage_line(const ms2::SampleStats& stats) {
  std::ostringstream text;
  const std::size_t holes = stats.holes;
  text << "positions: stretches "
       << percent(stats.stretch_characters, stats.positions) << " (strict "
       << percent(stats.strict_characters, stats.positions) << "), holes "
       << percent(stats.hole_positions(), stats.positions) << " (mean hole "
       << std::fixed << std::setprecision(1)
       << (holes == 0 ? 0.0
                      : static_cast<double>(stats.hole_positions()) /
                            static_cast<double>(holes))
       << " positions), separators "
       << percent(stats.separator_positions, stats.positions)
       << " of the input; overlapping stretches "
       << percent(stats.overlap_positions, stats.positions) << " of the input ("
       << percent(stats.overlap_positions, stats.stretch_characters)
       << " of the stretched positions)";
  return text.str();
}

// Its working arrays, reused across files and tables.
struct SampledBuffers {
  std::vector<KmerLookupResult> samples;
  std::vector<ms2::SampleKind> kinds;
  std::vector<std::uint32_t> chosen;
  std::vector<std::uint8_t> linked;
  std::vector<std::uint64_t> packed, sorted, scratch;  // sorted lookups
  Diagnostics lookup_phases;  // the last run's (sorted rows only)
};

// The sampled method on one table: steps 1-3 of sampled.hpp, each timed.
struct SampledRow {
  std::string name;
  std::function<void(const std::vector<Symbol>&, std::size_t, SampledBuffers&,
                     ms2::SampleStats&, double (&)[3])>
      run;

  double lookups_ms = 0.0, link_ms = 0.0, midpoints_ms = 0.0;
  double total_ms() const { return lookups_ms + link_ms + midpoints_ms; }

  // The sorted rows' lookup step split into collect / sort / lookups.
  Diagnostics lookup_phases;  // summed over files
};

// `LookupPolicy` does the sample lookups: in text order, or -- `sorted` --
// in sorted key order (a finger policy). The midpoint lookups always use
// `MidPolicy` (plain: they are few and in no useful order).
template <typename LookupPolicy, typename MidPolicy>
SampledRow make_sampled(std::string name,
                        std::shared_ptr<const typename MidPolicy::Table> table,
                        const bool sorted) {
  SampledRow row;
  row.name = std::move(name);
  row.run = [table, sorted](const std::vector<Symbol>& input,
                            const std::size_t step, SampledBuffers& b,
                            ms2::SampleStats& stats, double (&ms)[3]) {
    ms[0] = msbench::time_ms([&] {
      if (sorted) {
        ms2::sample_lookups_sorted<LookupPolicy>(*table, input, step,
                                                 b.samples, b.kinds, b.packed,
                                                 b.sorted, b.scratch,
                                                 &b.lookup_phases);
      } else {
        ms2::sample_lookups<LookupPolicy>(*table, input, step, b.samples,
                                          b.kinds);
        b.lookup_phases = Diagnostics{};
      }
    });
    ms[1] = msbench::time_ms([&] {
      ms2::link_samples(b.samples, b.kinds, step, b.chosen, b.linked, stats);
    });
    ms[2] = msbench::time_ms([&] {
      ms2::midpoint_lookups<MidPolicy>(*table, input, b.samples, b.kinds,
                                       b.chosen, b.linked, step, stats);
    });

    // Bookkeeping for the report (untimed).
    for (const ms2::SampleKind kind : b.kinds) {
      ++stats.samples;
      switch (kind) {
        case ms2::SampleKind::separator: ++stats.separators; break;
        case ms2::SampleKind::tail: ++stats.tails; break;
        case ms2::SampleKind::miss: ++stats.misses; break;
        case ms2::SampleKind::singleton: ++stats.singletons; break;
        case ms2::SampleKind::range: ++stats.ranges; break;
      }
    }
    if (!b.kinds.empty()) stats.pairs += b.kinds.size() - 1;
    for (std::size_t j = 0; j + 1 < b.linked.size(); ++j) {
      stats.linked += b.linked[j];
    }
    ms2::count_stretches(b.kinds, b.linked, step, input.size(), stats);
    ms2::strict_stretches(b.kinds, b.linked, step, input.size(), stats);
    ms2::overlap_stats(b.samples, b.kinds, b.chosen, b.linked, step,
                       input.size(), stats);
    stats.positions += input.size();
    for (const Symbol c : input) {
      stats.separator_positions += is_acgt(c) ? 0 : 1;
    }
  };
  return row;
}

// The sampled method end to end on one table.
struct SampledMsRow {
  std::string name;
  std::function<MatchingStatistics(const std::vector<Symbol>&, std::size_t,
                                   ms2::SampledMsBuffers&, ms2::SampledMsStats&,
                                   Diagnostics*)>
      run;

  double total_ms = 0.0;
  Diagnostics phases;  // summed over files
  std::size_t files_equal = 0;
  std::size_t files_compared = 0;
};

template <typename FingerPolicy, typename PlainPolicy>
SampledMsRow make_sampled_ms(
    std::string name, std::shared_ptr<const typename PlainPolicy::Table> table,
    const std::vector<Symbol>& reference) {
  SampledMsRow row;
  row.name = std::move(name);
  row.run = [table, &reference](const std::vector<Symbol>& input,
                                const std::size_t step,
                                ms2::SampledMsBuffers& buffers,
                                ms2::SampledMsStats& stats,
                                Diagnostics* phases) {
    return ms2::sampled_matching_statistics<FingerPolicy, PlainPolicy>(
        *table, input, reference, step, buffers, stats, phases);
  };
  return row;
}

// ---------- Totals ----------

struct Totals {
  std::size_t files = 0;
  std::size_t bytes = 0;
  std::size_t non_acgt = 0;
  std::uint64_t ms_length_sum = 0;

  double lrf_ms = 0.0;
  double keys_ms = 0.0;
  double bucket_ms = 0.0;
  double sorted_ms = 0.0;

  std::size_t keys = 0;
  std::size_t short_queries = 0;
  std::size_t separators = 0;
  std::size_t poly_a_runs = 0;
  std::size_t poly_a_positions = 0;
  std::size_t poly_t_runs = 0;
  std::size_t poly_t_positions = 0;
  std::size_t longest_run = 0;

  Diagnostics bucket_phases;
  Diagnostics sorted_phases;

  double results_ms = 0.0;
  double results_first_ms = 0.0;  // the first array only (the pipeline's one)
  Diagnostics results_phases;

  double chain_ms = 0.0;
  std::size_t chain_positions = 0;
  std::size_t chain_equal = 0;
  std::size_t chain_at_runs = 0;
  std::size_t chain_into_runs = 0;
  std::size_t chain_unexplained = 0;
  std::size_t chain_bad_positions = 0;

  bool checks_ok = true;
};

std::string percent(const std::size_t part, const std::size_t whole) {
  std::ostringstream text;
  text << std::fixed << std::setprecision(2)
       << (whole == 0 ? 0.0
                      : 100.0 * static_cast<double>(part) /
                            static_cast<double>(whole))
       << '%';
  return text.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::cerr << "========================================\n"
              << "ms2 benchmark\n"
              << "========================================\n";

    const Ms2Args args = parse_ms2_args(argc, argv);

    std::vector<Symbol> reference;
    const double reference_ms = msbench::time_ms(
        [&] { reference = msbench::load_reference<Symbol>(args.reference); });
    std::cerr << "[1] reference: " << reference.size() << " bytes in "
              << reference_ms << " ms\n";

    std::vector<SAType> suffix_array;
    const double suffix_array_ms = msbench::time_ms([&] {
      suffix_array = msbench::load_suffix_array<SAType>(args.suffix_array);
    });
    std::cerr << "[2] suffix array: " << suffix_array.size() << " entries in "
              << suffix_array_ms << " ms\n";

    if (suffix_array.size() != reference.size()) {
      throw std::runtime_error(
          "the suffix array must be built over the reference without a "
          "sentinel (sizes differ)");
    }

    std::vector<std::string> files = msbench::load_input_list(args.filenames);
    if (args.max_files != 0 && files.size() > args.max_files) {
      files.resize(args.max_files);
    }
    std::cerr << "[3] input files: " << files.size()
              << (args.max_bytes != 0
                      ? " (first " + std::to_string(args.max_bytes) +
                            " bytes of each)"
                      : std::string())
              << '\n';

    // Only A, C, G and T match: every other input byte is replaced by one
    // byte that does not occur in the reference, so lrf-ms never matches
    // it either (the same rule as the key pass, where it is a separator).
    const msbench::SymbolTable<Symbol> alphabet(reference);
    const Symbol separator = [&] {
      for (unsigned value = 1; value < 256; ++value) {
        const auto candidate = static_cast<Symbol>(value);
        if (!is_acgt(candidate) && !alphabet.contains(candidate)) {
          return candidate;
        }
      }
      throw std::runtime_error("no free byte left for a separator");
    }();

    // ---------- Build ----------

    std::cerr << "[4] build\n";

    std::unique_ptr<LRFMS<Symbol, SAType>> lrf;
    const double lrf_build_ms = msbench::time_ms([&] {
      lrf = std::make_unique<LRFMS<Symbol, SAType>>(reference, suffix_array);
    });

    Diagnostics::Line build_line{"build", true, {}};
    for (const auto& [name, ms] : lrf->buildPhases()) {
      build_line.values.emplace_back(name, ms);
    }
    const Diagnostics lrf_build{{build_line}};

    std::cerr << "    lrf-ms: " << lrf_build_ms << " ms\n"
              << lrf_build.format();

    // PT16 tables: always rebuilt, never reused from an earlier run.
    const std::string v2_path = args.table;
    const std::string sassy_path = args.table + ".sassy";

    // Plain v2 does not handle separators in the reference (sassy does);
    // on a reference with any non-ACGT byte it is left out.
    bool reference_is_acgt = true;
    for (unsigned value = 0; value < 256; ++value) {
      const auto symbol = static_cast<Symbol>(value);
      if (alphabet.contains(symbol) && !is_acgt(symbol)) {
        reference_is_acgt = false;
      }
    }

    std::vector<Variant> variants;

    std::vector<SampledRow> sampled_rows;
    std::shared_ptr<const PT16RLZParser<Symbol, SAType>> v2;

    if (reference_is_acgt) {
      const double write_ms = msbench::time_ms([&] {
        std::filesystem::remove(v2_path);
        build_pt16_table(reference, suffix_array, v2_path);
      });
      const double load_ms = msbench::time_ms([&] {
        v2 = std::make_shared<const PT16RLZParser<Symbol, SAType>>(
            reference, suffix_array, v2_path);
      });
      std::cerr << "    v2 table: write " << write_ms << " ms, load " << load_ms
                << " ms\n";
      variants.push_back(make_variant<msbench::V2Policy>("pt16-v2", v2));
      variants.push_back(make_variant<msbench::V2FingerPolicy>(
          "pt16-v2-finger", v2, /*sorted_only=*/true));
    } else {
      std::cerr << "    v2 table: skipped (the reference has non-ACGT bytes)\n";
    }

    RangeSizeStats range_stats;
    const double sassy_write_ms = msbench::time_ms([&] {
      std::filesystem::remove(sassy_path);
      range_stats = build_pt16_sassy_table(reference, suffix_array, sassy_path);
    });
    std::shared_ptr<const PT16SassyLookup> sassy;
    const double sassy_load_ms = msbench::time_ms(
        [&] { sassy = std::make_shared<const PT16SassyLookup>(sassy_path); });
    std::cerr << "    sassy table: write " << sassy_write_ms << " ms, load "
              << sassy_load_ms << " ms\n";

    variants.push_back(make_variant<msbench::SassyPolicy>("pt16-sassy", sassy));
    variants.push_back(make_variant<msbench::SassyFingerPolicy>(
        "pt16-sassy-finger", sassy, /*sorted_only=*/true));

    std::cerr << range_stats.format();

    if (v2) {
      sampled_rows.push_back(make_sampled<msbench::V2Policy, msbench::V2Policy>(
          "sampled-v2-textorder", v2, /*sorted=*/false));
      sampled_rows.push_back(
          make_sampled<msbench::V2FingerPolicy, msbench::V2Policy>(
              "sampled-v2-sorted", v2, /*sorted=*/true));
    }
    sampled_rows.push_back(
        make_sampled<msbench::SassyPolicy, msbench::SassyPolicy>(
            "sampled-sassy-textorder", sassy, /*sorted=*/false));
    sampled_rows.push_back(
        make_sampled<msbench::SassyFingerPolicy, msbench::SassyPolicy>(
            "sampled-sassy-sorted", sassy, /*sorted=*/true));

    std::vector<SampledMsRow> sampled_ms_rows;
    if (v2) {
      sampled_ms_rows.push_back(
          make_sampled_ms<msbench::V2FingerPolicy, msbench::V2Policy>(
              "sampled-ms-v2", v2, reference));
    }
    sampled_ms_rows.push_back(
        make_sampled_ms<msbench::SassyFingerPolicy, msbench::SassyPolicy>(
            "sampled-ms-sassy", sassy, reference));

    // Every variant on the bucket order (unless sorted-only), then every
    // variant on the sorted order. The first row's results are chained.
    std::vector<Row> rows;
    for (const bool sorted_rows : {false, true}) {
      for (std::size_t k = 0; k < variants.size(); ++k) {
        if (!sorted_rows && variants[k].sorted_only) continue;
        Row row;
        row.name = variants[k].name + (sorted_rows ? "-sorted" : "-bucket");
        row.variant = k;
        row.sorted = sorted_rows;
        rows.push_back(std::move(row));
      }
    }

    std::ofstream csv;
    if (!args.results.empty()) {
      csv.open(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "file,input_bytes,non_acgt,lrf_ms,keys_ms,bucket_ms,sorted_ms,"
             "results_ms";
      for (const Row& row : rows) {
        csv << ',' << row.name << "_ms";
      }
      csv << ",chain_ms,chain_equal,chain_at_runs,chain_into_runs,"
             "chain_unexplained,keys,short_queries,separators,poly_a_runs,"
             "poly_a_positions,poly_t_runs,poly_t_positions,longest_run,"
             "checks_ok\n";
    }

    // ---------- Per file ----------

    // --quiet: the per-file lines go nowhere (a stream without a buffer
    // discards what is written to it); only a file counter is shown.
    std::ostream discard(nullptr);
    std::ostream& out = args.quiet ? discard : std::cerr;

    out << "[5] files\n";

    Totals totals;
    SampledBuffers sampled_buffers;
    ms2::SampledMsBuffers sampled_ms_buffers;
    ms2::SampledMsStats total_sms_stats;
    ms2::SampleStats total_sample_stats;
    ms2::PhraseStats total_phrases;  // of lrf-ms's matching statistics
    ms2::KeyedInput keyed;
    std::vector<std::uint64_t> bucket, sorted, scratch;

    // The first row's lookup results (the ones chained), and the current
    // row's when it is not the first.
    std::vector<KmerLookupResult> first_results, other_results;

    for (std::size_t f = 0; f < files.size(); ++f) {
      if (args.quiet) {
        std::cerr << "\rfiles: " << f << "/" << files.size() << std::flush;
      }

      std::vector<Symbol> input;
      const double load_ms =
          msbench::time_ms([&] { input = msbench::load_input<Symbol>(files[f]); });

      if (args.max_bytes != 0 && input.size() > args.max_bytes) {
        input.resize(args.max_bytes);
      }

      std::size_t non_acgt = 0;
      for (Symbol& c : input) {
        if (!is_acgt(c)) {
          c = separator;
          ++non_acgt;
        }
      }

      out << "[" << f + 1 << "/" << files.size() << "] " << files[f]
                << "  (" << input.size() << " bytes, " << non_acgt
                << " non-ACGT, loaded in " << std::fixed
                << std::setprecision(1) << load_ms << " ms)\n";

      const auto step = [&out](const char* name) {
        out << "    " << std::left << std::setw(26) << name
                  << std::right << " ... " << std::flush;
      };

      // lrf-ms
      step("lrf-ms");
      MatchingStatistics ms;
      const double lrf_ms =
          msbench::time_ms([&] { ms = lrf->computeMatchingStatistics(input); });
      std::uint64_t length_sum = 0;
      for (const auto& [position, length] : ms) length_sum += length;
      const ms2::PhraseStats phrases = ms2::phrase_stats(ms);
      total_phrases.add(phrases);
      out << std::setw(10) << std::setprecision(2) << lrf_ms
                << " ms  " << phrases.format() << '\n';

      // keys
      step("keys");
      const double keys_ms =
          msbench::time_ms([&] { ms2::prepare_keys(input, keyed, args.keep_runs); });

      std::size_t longest_run = 0;
      for (const ms2::HomopolymerRun& run : keyed.runs) {
        longest_run = std::max<std::size_t>(longest_run, run.count + 15);
      }

      const std::size_t poly_a_runs = keyed.run_count('A');
      const std::size_t poly_t_runs = keyed.run_count('T');
      const std::size_t poly_a_positions = keyed.run_positions('A');
      const std::size_t poly_t_positions = keyed.run_positions('T');

      out << std::setw(10) << keys_ms << " ms  " << keyed.keys.size()
                << " keys, " << keyed.short_queries.size()
                << " short queries, " << keyed.separators
                << " separators (match 0, not searched)\n"
                << "        runs: poly-A " << poly_a_runs << " ("
                << poly_a_positions << " positions)  poly-T " << poly_t_runs
                << " (" << poly_t_positions << " positions)  = "
                << percent(poly_a_positions + poly_t_positions, input.size())
                << " of positions, longest " << longest_run
                << " characters\n";

      bool checks_ok = true;

      // bucket-order
      step("bucket-order");
      Diagnostics bucket_phases;
      const double bucket_ms = msbench::time_ms(
          [&] { ms2::bucket_order(keyed.keys, bucket, &bucket_phases); });
      out << std::setw(10) << bucket_ms << " ms";
      if (args.check) {
        const bool ok = is_ordering_of(keyed.keys, bucket, bucket_before);
        checks_ok = checks_ok && ok;
        out << (ok ? "  check OK" : "  CHECK FAILED");
      }
      out << '\n' << bucket_phases.format();

      // sorted-order
      step("sorted-order");
      Diagnostics sorted_phases;
      const double sorted_ms = msbench::time_ms([&] {
        ms2::sorted_order(keyed.keys, sorted, scratch, &sorted_phases);
      });
      out << std::setw(10) << sorted_ms << " ms";
      if (args.check) {
        const bool ok = is_ordering_of(keyed.keys, sorted, sorted_before);
        checks_ok = checks_ok && ok;
        out << (ok ? "  check OK" : "  CHECK FAILED");
      }
      out << '\n' << sorted_phases.format();

      // results arrays: prepared once per file, for every row
      step("results");
      Diagnostics results_phases;
      const double results_first_ms = msbench::time_ms([&] {
        ms2::prepare_results(input, keyed, first_results, &results_phases);
      });
      const double results_ms =
          results_first_ms + msbench::time_ms([&] {
            if (rows.size() > 1) {
              Diagnostics other;
              ms2::prepare_results(input, keyed, other_results, &other);
              results_phases.accumulate(other);
            }
          });
      out << std::setw(10) << results_ms << " ms  "
                << (rows.size() > 1 ? 2 : 1) << " arrays of "
                << keyed.input_size << " lookup results\n"
                << results_phases.format();
      totals.results_ms += results_ms;
      totals.results_first_ms += results_first_ms;
      totals.results_phases.accumulate(results_phases);

      // lookups, per row (variant x order)
      std::vector<double> probe_ms(rows.size(), 0.0);

      for (std::size_t k = 0; k < rows.size(); ++k) {
        Row& row = rows[k];
        const Variant& variant = variants[row.variant];
        const std::vector<std::uint64_t>& order = row.sorted ? sorted : bucket;
        std::vector<KmerLookupResult>& results =
            k == 0 ? first_results : other_results;

        step(row.name.c_str());
        Diagnostics diagnostics;
        probe_ms[k] = msbench::time_ms(
            [&] { variant.probe(keyed, order, results, &diagnostics); });
        out << std::setw(10) << probe_ms[k] << " ms";

        if (k > 0) {
          const msbench::LookupComparison comparison =
              msbench::compare_lookups(first_results, results, input,
                                       reference);
          ++row.files_compared;
          if (comparison.equal) {
            ++row.files_equal;
            out << "  lookups EQUAL";
          } else {
            checks_ok = false;
            out << "  lookups DIFFER: " << comparison.mismatches
                      << " positions, first at " << comparison.first_mismatch
                      << " (" << comparison.first_reason << ")";
          }
        }

        out << '\n' << diagnostics.format();
        row.total_ms += probe_ms[k];
        row.diagnostics.accumulate(diagnostics);
      }

      // the lookup table of the first file, once (untimed)
      if (f == 0 && !args.dump_lookups.empty() && !rows.empty()) {
        dump_lookup_table(args.dump_lookups, first_results);
        out << "    lookup table (" << rows.front().name << ") written to "
                  << args.dump_lookups << '\n';
      }

      // chain, on the first row's results
      double chain_ms = 0.0;
      ms2::ChainComparison chain;

      if (args.chain && !rows.empty()) {
        step("chain");

        MatchingStatistics chained;
        chain_ms = msbench::time_ms([&] {
          chained = ms2::chain_extend(first_results,
                                      [](std::size_t, std::size_t) {});
        });

        chain = ms2::compare_chain(ms, chained, keyed);
        const std::size_t bad_positions =
            ms2::bad_chain_positions(chained, input, reference, 10000);

        out << std::setw(10) << std::setprecision(2) << chain_ms
                  << " ms  lengths equal "
                  << chain.equal << "/" << chain.positions << "; differ: "
                  << chain.at_runs << " at run positions, " << chain.into_runs
                  << " reaching a run, " << chain.unexplained
                  << " unexplained";
        if (chain.unexplained != 0) {
          out << " (first at " << chain.first_unexplained << ": lrf-ms "
                    << chain.first_want << ", chain " << chain.first_got << ")";
        }
        out << "; sampled positions "
                  << (bad_positions == 0 ? "OK" : "BAD") << '\n';

        checks_ok = checks_ok && chain.unexplained == 0 && bad_positions == 0;

        totals.chain_ms += chain_ms;
        totals.chain_positions += chain.positions;
        totals.chain_equal += chain.equal;
        totals.chain_at_runs += chain.at_runs;
        totals.chain_into_runs += chain.into_runs;
        totals.chain_unexplained += chain.unexplained;
        totals.chain_bad_positions += bad_positions;
      }

      // The sampled method, per table; its statistics do not depend on
      // the table, so they are kept (and checked) for the first one.
      out << "    sampled method (step " << args.sample_step << "):\n";
      ms2::SampleStats file_sample_stats;
      for (std::size_t r = 0; r < sampled_rows.size(); ++r) {
        SampledRow& row = sampled_rows[r];
        ms2::SampleStats stats;
        double step_ms[3];
        row.run(input, args.sample_step, sampled_buffers, stats, step_ms);
        row.lookups_ms += step_ms[0];
        row.link_ms += step_ms[1];
        row.midpoints_ms += step_ms[2];
        row.lookup_phases.accumulate(sampled_buffers.lookup_phases);
        const Diagnostics file_lookup_phases = sampled_buffers.lookup_phases;

        step(row.name.c_str());
        out << std::setw(10) << step_ms[0] + step_ms[1] + step_ms[2] << " ms  (lookups "
            << step_ms[0] << ", link " << step_ms[1] << ", midpoints " << step_ms[2] << ")";

        // Every row's stretches are checked; and every row must find the
        // same links as the first (the lookups agree, so the links do).
        checks_ok = checks_ok && stats.overlap_into_singleton == 0;
        const std::size_t bad = ms2::check_stretches(
            sampled_buffers.kinds, sampled_buffers.chosen,
            sampled_buffers.linked, args.sample_step, input, reference, ms);
        bool same = true;
        if (r == 0) {
          file_sample_stats = stats;
        } else {
          same = stats.linked == file_sample_stats.linked &&
                 stats.stretches == file_sample_stats.stretches &&
                 stats.stretch_characters ==
                     file_sample_stats.stretch_characters &&
                 stats.holes == file_sample_stats.holes &&
                 stats.strict_characters == file_sample_stats.strict_characters;
        }
        if (bad == 0 && same) {
          out << "  stretches OK";
        } else {
          out << "  STRETCHES BAD (" << bad << " failing"
              << (same ? "" : ", differ from the first row") << ")";
        }
        checks_ok = checks_ok && bad == 0 && same;
        out << '\n' << file_lookup_phases.format();
      }
      out << "        " << coverage_line(file_sample_stats) << '\n'
          << file_sample_stats.diagnostics().format();
      total_sample_stats.add(file_sample_stats);

      // The sampled method end to end, per table (its statistics do not
      // depend on the table: kept for the first one).
      out << "    sampled method end to end (step " << args.sample_step
          << "):\n";
      ms2::SampledMsStats file_sms_stats;
      for (std::size_t r = 0; r < sampled_ms_rows.size(); ++r) {
        SampledMsRow& row = sampled_ms_rows[r];
        ms2::SampledMsStats stats;
        Diagnostics phases;
        MatchingStatistics got;
        const double t = msbench::time_ms([&] {
          got = row.run(input, args.sample_step, sampled_ms_buffers, stats,
                        &phases);
        });

        // Against lrf-ms at every position (untimed), and a sample of
        // reported positions really matching.
        std::size_t mismatches = 0, first = 0;
        for (std::size_t i = 0; i < ms.size(); ++i) {
          if (ms[i].second != got[i].second && mismatches++ == 0) first = i;
        }
        std::size_t bad = 0;
        const std::size_t stride = std::max<std::size_t>(1, got.size() / 10000);
        for (std::size_t i = 0; i < got.size(); i += stride) {
          const std::size_t p = got[i].first, len = got[i].second;
          if (len != 0 &&
              (p + len > reference.size() ||
               !std::equal(input.begin() + i, input.begin() + i + len,
                           reference.begin() +
                               static_cast<std::ptrdiff_t>(p)))) {
            ++bad;
          }
        }
        const bool ok = mismatches == 0 && bad == 0;
        checks_ok = checks_ok && ok;

        step(row.name.c_str());
        out << std::setw(10) << t << " ms  " << std::setprecision(3)
            << (t == 0.0 ? 0.0 : lrf_ms / t) << "x lrf-ms  "
            << (ok ? "EQUAL" : "DIFFER") << std::setprecision(2);
        if (mismatches != 0) {
          out << " (" << mismatches << " lengths, first at " << first
              << ": lrf-ms " << ms[first].second << ", sampled "
              << got[first].second << ")";
        }
        if (bad != 0) out << " (" << bad << " sampled positions do not match)";
        out << '\n' << phases.format();

        row.total_ms += t;
        row.phases.accumulate(phases);
        ++row.files_compared;
        row.files_equal += ok ? 1 : 0;
        if (r == 0) file_sms_stats = stats;
      }
      if (!sampled_ms_rows.empty()) {
        out << file_sms_stats.diagnostics().format();
        total_sms_stats.add(file_sms_stats);
      }

      // totals
      ++totals.files;
      totals.bytes += input.size();
      totals.non_acgt += non_acgt;
      totals.ms_length_sum += length_sum;
      totals.lrf_ms += lrf_ms;
      totals.keys_ms += keys_ms;
      totals.bucket_ms += bucket_ms;
      totals.sorted_ms += sorted_ms;
      totals.keys += keyed.keys.size();
      totals.short_queries += keyed.short_queries.size();
      totals.separators += keyed.separators;
      totals.poly_a_runs += poly_a_runs;
      totals.poly_a_positions += poly_a_positions;
      totals.poly_t_runs += poly_t_runs;
      totals.poly_t_positions += poly_t_positions;
      totals.longest_run = std::max(totals.longest_run, longest_run);
      totals.bucket_phases.accumulate(bucket_phases);
      totals.sorted_phases.accumulate(sorted_phases);
      totals.checks_ok = totals.checks_ok && checks_ok;

      if (csv.is_open()) {
        csv << files[f] << ',' << input.size() << ',' << non_acgt << ','
            << std::setprecision(3) << lrf_ms << ',' << keys_ms << ','
            << bucket_ms << ',' << sorted_ms << ',' << results_ms;
        for (const double ms_k : probe_ms) csv << ',' << ms_k;
        csv << ',' << chain_ms << ',' << chain.equal << ',' << chain.at_runs
            << ',' << chain.into_runs << ',' << chain.unexplained << ','
            << keyed.keys.size() << ','
            << keyed.short_queries.size() << ',' << keyed.separators << ','
            << poly_a_runs << ','
            << poly_a_positions << ',' << poly_t_runs << ','
            << poly_t_positions << ',' << longest_run << ','
            << (checks_ok || !args.check ? "YES" : "NO") << '\n';
      }
    }

    if (args.quiet) {
      std::cerr << "\rfiles: " << files.size() << "/" << files.size() << '\n';
    }

    // ---------- Totals ----------

    std::cerr << std::fixed;

    std::cerr << "[6] totals over " << totals.files << " files ("
              << totals.bytes << " bytes, " << totals.non_acgt
              << " non-ACGT)\n"
              << "    input: " << total_phrases.format() << '\n'
              << std::setprecision(2) << "    lrf-ms        " << std::setw(12)
              << totals.lrf_ms << " ms\n"
              << "    keys          " << std::setw(12) << totals.keys_ms
              << " ms\n"
              << "    bucket-order  " << std::setw(12) << totals.bucket_ms
              << " ms\n"
              << totals.bucket_phases.format()
              << "    sorted-order  " << std::setw(12) << totals.sorted_ms
              << " ms\n"
              << totals.sorted_phases.format()
              << "    results       " << std::setw(12) << totals.results_ms
              << " ms\n"
              << totals.results_phases.format();

    for (const Row& row : rows) {
      std::cerr << "    " << std::left << std::setw(26) << row.name
                << std::right << std::setw(12) << row.total_ms << " ms";
      if (row.files_compared != 0) {
        std::cerr << "  lookups equal " << row.files_equal << "/"
                  << row.files_compared << " files";
      }
      std::cerr << '\n' << row.diagnostics.format();
    }

    if (args.chain) {
      std::cerr << "    chain         " << std::setw(12) << totals.chain_ms
                << " ms  lengths equal " << totals.chain_equal << "/"
                << totals.chain_positions << "; differ: "
                << totals.chain_at_runs << " at run positions, "
                << totals.chain_into_runs << " reaching a run, "
                << totals.chain_unexplained << " unexplained\n";

      // Everything after lrf-ms, for each row: keys + its order + one
      // results array + its lookups + the chain.
      std::cerr << "    pipeline (keys + order + results + lookups + chain):\n";
      for (const Row& row : rows) {
        const double pipeline =
            totals.keys_ms + (row.sorted ? totals.sorted_ms : totals.bucket_ms) +
            totals.results_first_ms + row.total_ms + totals.chain_ms;
        std::cerr << "      " << std::left << std::setw(26) << row.name
                  << std::right << std::setw(12) << pipeline << " ms  lrf-ms / pipeline "
                  << std::setprecision(3)
                  << (pipeline == 0.0 ? 0.0 : totals.lrf_ms / pipeline) << "x\n"
                  << std::setprecision(2);
      }
    }

    std::cerr << "    keys " << totals.keys << " ("
              << percent(totals.keys, totals.bytes) << "), short queries "
              << totals.short_queries << ", separators " << totals.separators
              << ", poly-A runs " << totals.poly_a_runs
              << " (" << totals.poly_a_positions << " positions), poly-T runs "
              << totals.poly_t_runs << " (" << totals.poly_t_positions
              << " positions), longest run " << totals.longest_run << '\n'
              << "    peak RSS " << msbench::peak_rss_mb() << " MB\n";

    std::cerr << "    sampled method (step " << args.sample_step << "):\n";
    for (const SampledRow& row : sampled_rows) {
      std::cerr << "      " << std::left << std::setw(20) << row.name
                << std::right << std::setw(12) << row.total_ms()
                << " ms  (lookups " << row.lookups_ms << ", link "
                << row.link_ms << ", midpoints " << row.midpoints_ms
                << ")  " << std::setprecision(3)
                << (row.total_ms() == 0.0 ? 0.0 : row.total_ms() / totals.lrf_ms)
                << "x lrf-ms's time\n"
                << std::setprecision(2) << row.lookup_phases.format("          ");
    }
    std::cerr << "      " << coverage_line(total_sample_stats) << '\n'
              << total_sample_stats.diagnostics().format("        ");

    std::cerr << "    sampled method end to end (step " << args.sample_step
              << "):\n";
    for (const SampledMsRow& row : sampled_ms_rows) {
      std::cerr << "      " << std::left << std::setw(20) << row.name
                << std::right << std::setw(12) << row.total_ms << " ms  "
                << std::setprecision(3)
                << (row.total_ms == 0.0 ? 0.0 : totals.lrf_ms / row.total_ms)
                << "x lrf-ms  equal " << row.files_equal << "/"
                << row.files_compared << " files\n"
                << std::setprecision(2) << row.phases.format("          ");
    }
    std::cerr << total_sms_stats.diagnostics().format("        ");

    if (args.check) {
      std::cerr << "    checks: " << (totals.checks_ok ? "all OK" : "FAILED")
                << '\n';
    }

    // Machine-readable.
    std::cout << "files=" << totals.files << '\n'
              << "input_bytes=" << totals.bytes << '\n'
              << "lrf_build_ms=" << lrf_build_ms << '\n'
              << "lrf_ms=" << totals.lrf_ms << '\n'
              << "avg_match_length=" << total_phrases.average_match_length()
              << '\n'
              << "short_ms_percent=" << total_phrases.short_ms_percent() << '\n'
              << "avg_file_short_ms_percent="
              << total_phrases.average_file_short_ms_percent() << '\n'
              << "max_match_length=" << total_phrases.max_length << '\n'
              << "median_match_length=" << total_phrases.median_length() << '\n'
              << "unmatched_positions=" << total_phrases.unmatched << '\n'
              << "ms_length_sum=" << totals.ms_length_sum << '\n'
              << "keys_ms=" << totals.keys_ms << '\n'
              << "results_ms=" << totals.results_ms << '\n'
              << "bucket_ms=" << totals.bucket_ms << '\n'
              << "sorted_ms=" << totals.sorted_ms << '\n'
              << "keys=" << totals.keys << '\n'
              << "short_queries=" << totals.short_queries << '\n'
              << "separators=" << totals.separators << '\n'
              << "poly_a_runs=" << totals.poly_a_runs << '\n'
              << "poly_a_positions=" << totals.poly_a_positions << '\n'
              << "poly_t_runs=" << totals.poly_t_runs << '\n'
              << "poly_t_positions=" << totals.poly_t_positions << '\n'
              << "longest_run=" << totals.longest_run << '\n'
              << "chain_ms=" << totals.chain_ms << '\n'
              << "chain_equal=" << totals.chain_equal << '\n'
              << "chain_at_runs=" << totals.chain_at_runs << '\n'
              << "chain_into_runs=" << totals.chain_into_runs << '\n'
              << "chain_unexplained=" << totals.chain_unexplained << '\n';
    for (const Row& row : rows) {
      std::cout << row.name << "_ms=" << row.total_ms << '\n';
    }
    for (const SampledRow& row : sampled_rows) {
      std::cout << row.name << "_ms=" << row.total_ms() << '\n'
                << row.name << "_lookups_ms=" << row.lookups_ms << '\n'
                << row.name << "_link_ms=" << row.link_ms << '\n'
                << row.name << "_midpoints_ms=" << row.midpoints_ms << '\n';
    }
    for (const SampledMsRow& row : sampled_ms_rows) {
      std::cout << row.name << "_ms=" << row.total_ms << '\n'
                << row.name << "_equal_files=" << row.files_equal << '\n';
    }
    std::cout << "sampled_ms_lookups=" << total_sms_stats.lookups() << '\n'
              << "sampled_ms_stretch_positions="
              << total_sms_stats.stretch_positions << '\n';
    std::cout << "sample_step=" << args.sample_step << '\n'
              << "samples=" << total_sample_stats.samples << '\n'
              << "samples_linked=" << total_sample_stats.linked << '\n'
              << "stretches=" << total_sample_stats.stretches << '\n'
              << "stretch_characters=" << total_sample_stats.stretch_characters
              << '\n'
              << "strict_stretches=" << total_sample_stats.strict_stretches
              << '\n'
              << "strict_characters=" << total_sample_stats.strict_characters
              << '\n'
              << "ranges_dismissed=" << total_sample_stats.ranges_dismissed
              << '\n'
              << "ranges_kept=" << total_sample_stats.ranges_kept << '\n'
              << "overlapping_ends=" << total_sample_stats.overlapping_ends
              << '\n'
              << "overlap_positions=" << total_sample_stats.overlap_positions
              << '\n'
              << "overlap_into_singleton="
              << total_sample_stats.overlap_into_singleton << '\n'
              << "holes=" << total_sample_stats.holes << '\n'
              << "hole_positions=" << total_sample_stats.hole_positions() << '\n'
              << "hole_percent="
              << (total_sample_stats.positions == 0
                      ? 0.0
                      : 100.0 *
                            static_cast<double>(
                                total_sample_stats.hole_positions()) /
                            static_cast<double>(total_sample_stats.positions))
              << '\n';
    std::cout << "checks_ok=" << (totals.checks_ok ? "YES" : "NO") << '\n';

    return totals.checks_ok ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
