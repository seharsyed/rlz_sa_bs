// Experimental benchmark: PT16 lookups on the input's unique 16-mers only.
//
//   g++ -std=c++20 -O3 ms2/unique_main.cpp -o unique_main
//   ./unique_main --reference REF --suffix-array REF.sa --filenames LIST
//                 [--results CSV] [--table PATH] [--max-files N]
//                 [--max-bytes N] [--quiet]
//
// Built once: lrf-ms (for its time, as a yardstick) and the PT16 tables
// (v2 when the reference is all ACGT, sassy). Per input file:
//
//   lrf-ms         full matching statistics, for comparison;
//   keys           one pass over the input (input_keys.hpp);
//   sorted-order   the keys sorted by key (LSD radix);
//   unique         keep only the keys that occur once in the input
//                  (unique_keys.hpp): one pass over the sorted list;
//   results        the two lookup-result arrays (the first variant's and
//                  the current one's) sized for this file's unique keys,
//                  so no variant pays for allocating them;
//   <variant>      the lookups of the unique keys, in sorted order, into a
//                  compact array aligned with the unique list: pt16-v2,
//                  pt16-v2-finger, pt16-sassy, pt16-sassy-finger. Every
//                  variant's results are checked against the first one's;
//   text-order     the unique keys put back in input text order (radix
//                  sort of position << 32 | index; unique_chain.hpp);
//   chain          one right-to-left sweep over them (chain_unique):
//                  misses are exact, a found 16-mer that lines up with the
//                  entry to its right (at most 16 positions away) takes its
//                  length + the distance, anything else gets a lower bound
//                  of 16 and is marked approximate. Checked against lrf-ms:
//                  exact lengths must equal it, approximate ones must not
//                  exceed it.
// Each file's lines print as they are measured (unless
// --quiet); totals at the end.

#include <algorithm>
#include <cstdint>
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
#include "unique_chain.hpp"
#include "unique_keys.hpp"

using Symbol = unsigned char;
using SAType = std::uint32_t;

namespace {

struct UniqueArgs {
  std::string reference;
  std::string suffix_array;
  std::string filenames;
  std::string results;  // optional CSV
  std::string table;    // PT16 table files; default beside the reference
  std::size_t max_files = 0;
  std::size_t max_bytes = 0;
  bool quiet = false;
  bool keep_runs = false;  // all-A / all-T 16-mers as ordinary keys
};

void print_usage(const char* program) {
  std::cout << "Usage: " << program
            << " --reference PATH --suffix-array PATH --filenames PATH\n"
            << "  [--results PATH]    per-file CSV\n"
            << "  [--table PATH]      PT16 table files (default: "
               "<reference>.unique_pt16, and .sassy)\n"
            << "  [--max-files N]     only the first N input files\n"
            << "  [--max-bytes N]     only the first N bytes of each input\n"
            << "  [--quiet]           no per-file lines, only the totals\n"
            << "  [--keep-runs]       look up all-A / all-T 16-mers like any "
               "other key\n"
            << "                      (default: set them aside as runs)\n";
}

UniqueArgs parse_unique_args(int argc, char** argv) {
  UniqueArgs args;

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
    } else if (option == "--table") {
      args.table = msbench::require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = msbench::parse_number(
          option, msbench::require_value(i, argc, argv));
    } else if (option == "--max-bytes") {
      args.max_bytes = msbench::parse_number(
          option, msbench::require_value(i, argc, argv));
    } else if (option == "--quiet") {
      args.quiet = true;
    } else if (option == "--keep-runs") {
      args.keep_runs = true;
    } else if (option == "--help" || option == "-h") {
      print_usage(argv[0]);
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }

  if (args.reference.empty() || args.suffix_array.empty() ||
      args.filenames.empty()) {
    throw std::runtime_error(
        "reference, suffix-array and filenames are required");
  }

  if (args.table.empty()) {
    args.table = args.reference + ".unique_pt16";
  }

  return args;
}

// One table and one way of searching it, over the sorted unique keys.
struct Variant {
  std::string name;
  std::function<void(const std::vector<std::uint64_t>&,
                     std::vector<KmerLookupResult>&, Diagnostics*)>
      probe;

  double total_ms = 0.0;
  Diagnostics diagnostics;  // summed over files
  std::size_t files_compared = 0;
  std::size_t files_equal = 0;
};

template <typename Policy>
Variant make_variant(std::string name,
                     std::shared_ptr<const typename Policy::Table> table) {
  Variant variant;
  variant.name = std::move(name);
  variant.probe = [table](const std::vector<std::uint64_t>& order,
                          std::vector<KmerLookupResult>& results,
                          Diagnostics* diagnostics) {
    ms2::probe_compact<Policy>(*table, order, results, diagnostics);
  };
  return variant;
}

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
              << "unique 16-mer lookup benchmark\n"
              << "========================================\n";

    const UniqueArgs args = parse_unique_args(argc, argv);

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
    // byte that does not occur in the reference.
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
    std::cerr << "    lrf-ms: " << lrf_build_ms << " ms\n";

    const std::string v2_path = args.table;
    const std::string sassy_path = args.table + ".sassy";

    bool reference_is_acgt = true;
    for (unsigned value = 0; value < 256; ++value) {
      const auto symbol = static_cast<Symbol>(value);
      if (alphabet.contains(symbol) && !is_acgt(symbol)) {
        reference_is_acgt = false;
      }
    }

    std::vector<Variant> variants;

    if (reference_is_acgt) {
      const double write_ms = msbench::time_ms([&] {
        std::filesystem::remove(v2_path);
        build_pt16_table(reference, suffix_array, v2_path);
      });
      std::shared_ptr<const PT16RLZParser<Symbol, SAType>> v2;
      const double load_ms = msbench::time_ms([&] {
        v2 = std::make_shared<const PT16RLZParser<Symbol, SAType>>(
            reference, suffix_array, v2_path);
      });
      std::cerr << "    v2 table: write " << write_ms << " ms, load " << load_ms
                << " ms\n";
      variants.push_back(make_variant<msbench::V2Policy>("pt16-v2", v2));
      variants.push_back(
          make_variant<msbench::V2FingerPolicy>("pt16-v2-finger", v2));
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
              << sassy_load_ms << " ms\n"
              << range_stats.format();

    variants.push_back(make_variant<msbench::SassyPolicy>("pt16-sassy", sassy));
    variants.push_back(
        make_variant<msbench::SassyFingerPolicy>("pt16-sassy-finger", sassy));

    std::ofstream csv;
    if (!args.results.empty()) {
      csv.open(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "file,input_bytes,non_acgt,lrf_ms,keys_ms,sorted_ms,unique_ms,"
             "results_ms,keys,unique,repeated_keys,repeated_positions";
      for (const Variant& variant : variants) csv << ',' << variant.name << "_ms";
      csv << ",text_order_ms,chain_ms,chain_exact,chain_approximate,"
             "chain_approximate_right,checks_ok\n";
    }

    // ---------- Per file ----------

    // --quiet: the per-file lines go nowhere (a stream without a buffer
    // discards what is written to it); only a file counter is shown.
    std::ostream discard(nullptr);
    std::ostream& out = args.quiet ? discard : std::cerr;

    out << "[5] files\n";

    std::size_t total_bytes = 0, total_non_acgt = 0;
    double total_lrf_ms = 0.0, total_keys_ms = 0.0, total_sorted_ms = 0.0,
           total_unique_ms = 0.0, total_results_ms = 0.0,
           total_results_first_ms = 0.0;
    std::size_t total_keys = 0, total_unique = 0, total_repeated_keys = 0,
                total_repeated_positions = 0;
    Diagnostics total_sorted_phases;
    ms2::PhraseStats total_phrases;  // of lrf-ms's matching statistics
    Diagnostics total_unique_phases;
    bool all_ok = true;

    ms2::KeyedInput keyed;
    std::vector<std::uint64_t> sorted, scratch, unique, packed, by_position;
    std::vector<KmerLookupResult> first_results, other_results;
    std::vector<ms2::ChainEntry> chained;

    double total_text_ms = 0.0, total_chain_ms = 0.0;
    Diagnostics total_text_phases, total_composition;
    ms2::UniqueChainStats total_chain;
    std::size_t total_approx_right = 0;

    for (std::size_t f = 0; f < files.size(); ++f) {
      if (args.quiet) {
        std::cerr << "\rfiles: " << f << "/" << files.size() << std::flush;
      }

      std::vector<Symbol> input;
      const double load_ms = msbench::time_ms(
          [&] { input = msbench::load_input<Symbol>(files[f]); });

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

      out << "[" << f + 1 << "/" << files.size() << "] " << files[f] << "  ("
          << input.size() << " bytes, " << non_acgt << " non-ACGT, loaded in "
          << std::fixed << std::setprecision(1) << load_ms << " ms)\n";

      const auto step = [&out](const std::string& name) {
        out << "    " << std::left << std::setw(20) << name << std::right
            << " ... " << std::flush;
      };

      // lrf-ms, for comparison
      step("lrf-ms");
      MatchingStatistics want;
      const double lrf_ms = msbench::time_ms(
          [&] { want = lrf->computeMatchingStatistics(input); });

      const ms2::PhraseStats phrases = ms2::phrase_stats(want);
      total_phrases.add(phrases);
      out << std::setw(10) << std::setprecision(2) << lrf_ms << " ms  "
          << phrases.format() << '\n';

      // keys
      step("keys");
      const double keys_ms =
          msbench::time_ms([&] { ms2::prepare_keys(input, keyed, args.keep_runs); });
      out << std::setw(10) << keys_ms << " ms  " << keyed.keys.size()
          << " keys\n";

      // sorted-order
      step("sorted-order");
      Diagnostics sorted_phases;
      const double sorted_ms = msbench::time_ms([&] {
        ms2::sorted_order(keyed.keys, sorted, scratch, &sorted_phases);
      });
      out << std::setw(10) << sorted_ms << " ms\n" << sorted_phases.format();

      // unique
      step("unique");
      ms2::UniqueStats unique_stats;
      Diagnostics unique_phases;
      const double unique_ms = msbench::time_ms([&] {
        unique_stats = ms2::keep_unique(sorted, unique, &unique_phases);
      });
      out << std::setw(10) << unique_ms << " ms  " << unique_stats.unique
          << " unique keys (" << percent(unique_stats.unique, unique_stats.keys)
          << " of keys); dropped " << unique_stats.repeated_positions
          << " positions of " << unique_stats.repeated_keys
          << " repeated keys\n"
          << unique_phases.format();

      // results arrays, sized once per file for every variant
      step("results");
      const double results_first_ms =
          msbench::time_ms([&] { first_results.resize(unique.size()); });
      const double results_ms =
          results_first_ms +
          msbench::time_ms([&] { other_results.resize(unique.size()); });
      out << std::setw(10) << results_ms << " ms  2 arrays of " << unique.size()
          << " lookup results\n";

      // lookups, per variant
      bool checks_ok = true;
      std::vector<double> probe_ms(variants.size(), 0.0);

      for (std::size_t k = 0; k < variants.size(); ++k) {
        Variant& variant = variants[k];
        std::vector<KmerLookupResult>& results =
            k == 0 ? first_results : other_results;

        step(variant.name);
        Diagnostics diagnostics;
        probe_ms[k] = msbench::time_ms(
            [&] { variant.probe(unique, results, &diagnostics); });
        out << std::setw(10) << probe_ms[k] << " ms";

        // The first variant is checked for genuine positions against
        // itself; every other one against it.
        std::size_t at = 0;
        std::string reason;
        const std::size_t mismatches = ms2::compare_compact(
            first_results, results, unique, input, reference, at, reason);
        ++variant.files_compared;
        if (mismatches == 0) {
          ++variant.files_equal;
          out << (k == 0 ? "  positions OK" : "  lookups EQUAL");
        } else {
          checks_ok = false;
          out << "  " << (k == 0 ? "BAD POSITIONS" : "lookups DIFFER") << ": "
              << mismatches << ", first at " << at << " (" << reason << ")";
        }
        out << '\n';
        out << diagnostics.format();

        variant.total_ms += probe_ms[k];
        variant.diagnostics.accumulate(diagnostics);
      }

      // How the lookups resolved: the same for every variant, so once.
      const Diagnostics composition = ms2::lookup_composition(first_results);
      out << "    lookups (every variant):" << '\n' << composition.format();

      // text-order
      step("text-order");
      Diagnostics text_phases;
      const double text_ms = msbench::time_ms([&] {
        ms2::text_order(unique, by_position, packed, scratch, &text_phases);
      });
      out << std::setw(10) << text_ms << " ms\n" << text_phases.format();

      // chain, on the first variant's results
      step("chain");
      ms2::UniqueChainStats chain_stats;
      const double chain_ms = msbench::time_ms([&] {
        chain_stats = ms2::chain_unique(first_results, by_position, chained);
      });

      // Against lrf-ms (untimed): exact lengths equal, approximate ones not
      // above it; and how many approximate ones are right anyway.
      std::size_t exact_wrong = 0, approx_above = 0, approx_right = 0,
                  bad_positions = 0;
      for (std::size_t idx = 0; idx < chained.size(); ++idx) {
        const std::uint32_t i = ms2::key_of(by_position[idx]);
        const ms2::ChainEntry& entry = chained[idx];
        const std::uint32_t truth = want[i].second;

        if (ms2::is_exact(entry.status)) {
          exact_wrong += entry.length != truth;
        } else {
          approx_above += entry.length > truth;
          approx_right += entry.length == truth;
        }

        const std::size_t p = entry.reference_position;
        const std::size_t len = entry.length;
        if (p + len > reference.size() || i + len > input.size() ||
            !std::equal(input.begin() + i, input.begin() + i + len,
                        reference.begin() + static_cast<std::ptrdiff_t>(p))) {
          ++bad_positions;
        }
      }
      const bool chain_ok =
          exact_wrong == 0 && approx_above == 0 && bad_positions == 0;
      checks_ok = checks_ok && chain_ok;

      out << std::setw(10) << chain_ms << " ms  exact "
          << chain_stats.exact() << " ("
          << percent(chain_stats.exact(), chained.size()) << "), approximate "
          << chain_stats.approximate() << " (" << approx_right
          << " of them right anyway)  "
          << (chain_ok ? "check OK" : "CHECK FAILED");
      if (!chain_ok) {
        out << " (" << exact_wrong << " exact wrong, " << approx_above
            << " approximate above lrf-ms, " << bad_positions
            << " positions not matching)";
      }
      out << '\n' << chain_stats.diagnostics().format();

      total_text_ms += text_ms;
      total_chain_ms += chain_ms;
      total_text_phases.accumulate(text_phases);
      total_composition.accumulate(composition);
      total_chain.add(chain_stats);
      total_approx_right += approx_right;

      total_bytes += input.size();
      total_non_acgt += non_acgt;
      total_lrf_ms += lrf_ms;
      total_keys_ms += keys_ms;
      total_sorted_ms += sorted_ms;
      total_unique_ms += unique_ms;
      total_results_ms += results_ms;
      total_results_first_ms += results_first_ms;
      total_keys += unique_stats.keys;
      total_unique += unique_stats.unique;
      total_repeated_keys += unique_stats.repeated_keys;
      total_repeated_positions += unique_stats.repeated_positions;
      total_sorted_phases.accumulate(sorted_phases);
      total_unique_phases.accumulate(unique_phases);
      all_ok = all_ok && checks_ok;

      if (csv.is_open()) {
        csv << files[f] << ',' << input.size() << ',' << non_acgt << ','
            << std::setprecision(3) << lrf_ms << ',' << keys_ms << ','
            << sorted_ms << ',' << unique_ms << ',' << results_ms << ','
            << unique_stats.keys << ','
            << unique_stats.unique << ',' << unique_stats.repeated_keys << ','
            << unique_stats.repeated_positions;
        for (const double ms : probe_ms) csv << ',' << ms;
        csv << ',' << text_ms << ',' << chain_ms << ',' << chain_stats.exact()
            << ',' << chain_stats.approximate() << ',' << approx_right;
        csv << ',' << (checks_ok ? "YES" : "NO") << '\n';
      }
    }

    if (args.quiet) {
      std::cerr << "\rfiles: " << files.size() << "/" << files.size() << '\n';
    }

    // ---------- Totals ----------

    std::cerr << std::fixed << std::setprecision(2) << "[6] totals over "
              << files.size() << " files (" << total_bytes << " bytes, "
              << total_non_acgt << " non-ACGT)\n"
              << "    input: " << total_phrases.format() << '\n'
              << "    lrf-ms              " << std::setw(12) << total_lrf_ms
              << " ms  (build " << lrf_build_ms << " ms)\n"
              << "    keys                " << std::setw(12) << total_keys_ms
              << " ms  " << total_keys << " keys\n"
              << "    sorted-order        " << std::setw(12) << total_sorted_ms
              << " ms\n"
              << total_sorted_phases.format()
              << "    unique              " << std::setw(12) << total_unique_ms
              << " ms  " << total_unique << " unique keys ("
              << percent(total_unique, total_keys) << "); dropped "
              << total_repeated_positions << " positions of "
              << total_repeated_keys << " repeated keys\n"
              << total_unique_phases.format()
              << "    results             " << std::setw(12) << total_results_ms
              << " ms  (2 arrays)\n";

    for (const Variant& variant : variants) {
      std::cerr << "    " << std::left << std::setw(20) << variant.name
                << std::right << std::setw(12) << variant.total_ms << " ms  "
                << (variant.files_equal == variant.files_compared ? "equal "
                                                                  : "DIFFER ")
                << variant.files_equal << "/" << variant.files_compared
                << " files\n"
                << variant.diagnostics.format();
    }

    std::cerr << "    lookups (every variant):\n" << total_composition.format()
              << "    text-order          " << std::setw(12) << total_text_ms
              << " ms\n"
              << total_text_phases.format()
              << "    chain               " << std::setw(12) << total_chain_ms
              << " ms  exact " << total_chain.exact() << " ("
              << percent(total_chain.exact(), total_unique) << "), approximate "
              << total_chain.approximate() << " (" << total_approx_right
              << " of them right anyway)\n"
              << total_chain.diagnostics().format();

    std::cerr << "    pipeline (keys + sorted-order + unique + one results "
                 "array + lookups + text-order + chain):\n";
    for (const Variant& variant : variants) {
      const double pipeline = total_keys_ms + total_sorted_ms +
                              total_unique_ms + total_results_first_ms +
                              variant.total_ms + total_text_ms + total_chain_ms;
      std::cerr << "      " << std::left << std::setw(20) << variant.name
                << std::right << std::setw(12) << pipeline << " ms  ("
                << std::setprecision(3)
                << (total_lrf_ms == 0.0 ? 0.0 : pipeline / total_lrf_ms)
                << "x lrf-ms's time)\n"
                << std::setprecision(2);
    }

    std::cerr << "    peak RSS " << msbench::peak_rss_mb() << " MB\n";

    // Machine-readable.
    std::cout << "files=" << files.size() << '\n'
              << "input_bytes=" << total_bytes << '\n'
              << "lrf_ms=" << total_lrf_ms << '\n'
              << "avg_match_length=" << total_phrases.average_match_length()
              << '\n'
              << "short_ms_percent=" << total_phrases.short_ms_percent() << '\n'
              << "avg_file_short_ms_percent="
              << total_phrases.average_file_short_ms_percent() << '\n'
              << "max_match_length=" << total_phrases.max_length << '\n'
              << "median_match_length=" << total_phrases.median_length() << '\n'
              << "unmatched_positions=" << total_phrases.unmatched << '\n'
              << "keys_ms=" << total_keys_ms << '\n'
              << "sorted_ms=" << total_sorted_ms << '\n'
              << "unique_ms=" << total_unique_ms << '\n'
              << "results_ms=" << total_results_ms << '\n'
              << "keys=" << total_keys << '\n'
              << "unique=" << total_unique << '\n'
              << "repeated_keys=" << total_repeated_keys << '\n'
              << "repeated_positions=" << total_repeated_positions << '\n';
    for (const Variant& variant : variants) {
      std::cout << variant.name << "_ms=" << variant.total_ms << '\n';
    }
    std::cout << "text_order_ms=" << total_text_ms << '\n'
              << "chain_ms=" << total_chain_ms << '\n'
              << "chain_exact=" << total_chain.exact() << '\n'
              << "chain_exact_chained=" << total_chain.exact_chained << '\n'
              << "chain_approximate=" << total_chain.approximate() << '\n'
              << "chain_approximate_right=" << total_approx_right << '\n'
              << "all_ok=" << (all_ok ? "YES" : "NO") << '\n';

    return all_ok ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
