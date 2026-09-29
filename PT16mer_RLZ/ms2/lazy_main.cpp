// Lazy matching-statistics benchmark: backward left-extension with PT16
// lookups only at breaks (lazy_ms.hpp), against lrf-ms.
//
//   g++ -std=c++20 -O3 ms2/lazy_main.cpp -o lazy_main
//   ./lazy_main --reference REF --suffix-array REF.sa --filenames LIST
//               [--results CSV] [--table PATH] [--max-files N]
//               [--max-bytes N] [--sample N] [--quiet]
//
// Built once: lrf-ms and the PT16 tables (v2 when the reference is all
// ACGT, sassy). Per input file: lrf-ms, then each lazy variant; each
// variant's lengths are checked against lrf-ms at every position, and its
// reported positions at --sample evenly spaced positions (default 10000).
// Each file's lines print as soon as they are measured; totals at the end.

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
#include "../lrf_ms/ms_utils.hpp"        // loaders, SymbolTable
#include "../lrf_ms/probe_pipeline.hpp"  // table builders, lookup policies
#include "lazy_ms.hpp"
#include "phrase_stats.hpp"

using Symbol = unsigned char;
using SAType = std::uint32_t;

namespace {

struct LazyArgs {
  std::string reference;
  std::string suffix_array;
  std::string filenames;
  std::string results;  // optional CSV
  std::string table;    // PT16 table files; default beside the reference
  std::size_t max_files = 0;
  std::size_t max_bytes = 0;
  std::size_t sample = 10000;
  bool quiet = false;  // no per-file lines, only totals
};

void print_usage(const char* program) {
  std::cout << "Usage: " << program
            << " --reference PATH --suffix-array PATH --filenames PATH\n"
            << "  [--results PATH]    per-file CSV\n"
            << "  [--table PATH]      PT16 table files (default: "
               "<reference>.lazy_pt16, and .sassy)\n"
            << "  [--max-files N]     only the first N input files\n"
            << "  [--max-bytes N]     only the first N bytes of each input\n"
            << "  [--sample N]        positions whose reported match is "
               "checked (default 10000)\n"
            << "  [--quiet]           no per-file lines, only the totals\n";
}

LazyArgs parse_lazy_args(int argc, char** argv) {
  LazyArgs args;

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
    } else if (option == "--sample") {
      args.sample = msbench::parse_number(
          option, msbench::require_value(i, argc, argv));
    } else if (option == "--quiet") {
      args.quiet = true;
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
    args.table = args.reference + ".lazy_pt16";
  }

  return args;
}

// One lazy variant: a table, searched through its policy.
struct Variant {
  std::string name;
  std::function<MatchingStatistics(const std::vector<Symbol>&,
                                   const std::vector<Symbol>&,
                                   ms2::LazyStats*)>
      compute;

  double total_ms = 0.0;
  ms2::LazyStats stats;  // summed over files
  std::size_t files_equal = 0;
  std::size_t files_compared = 0;
  std::size_t length_mismatches = 0;
  std::size_t bad_positions = 0;
};

template <typename Policy>
Variant make_variant(std::string name,
                     std::shared_ptr<const typename Policy::Table> table) {
  Variant variant;
  variant.name = std::move(name);
  variant.compute = [table](const std::vector<Symbol>& input,
                            const std::vector<Symbol>& reference,
                            ms2::LazyStats* stats) {
    return ms2::lazy_matching_statistics<Policy>(*table, input, reference,
                                                 stats);
  };
  return variant;
}

// Positions whose length differs from the baseline, and the first one.
std::size_t length_mismatches(const MatchingStatistics& want,
                              const MatchingStatistics& got,
                              std::size_t& first) {
  std::size_t count = 0;
  for (std::size_t i = 0; i < want.size(); ++i) {
    if (want[i].second != got[i].second && count++ == 0) first = i;
  }
  return count;
}

// Of up to `samples` evenly spaced positions, how many report a position
// that does not really match the input for the reported length.
std::size_t bad_positions(const MatchingStatistics& got,
                          const std::vector<Symbol>& input,
                          const std::vector<Symbol>& reference,
                          const std::size_t samples) {
  if (got.empty() || samples == 0) return 0;

  const std::size_t stride = std::max<std::size_t>(1, got.size() / samples);
  std::size_t bad = 0;

  for (std::size_t i = 0; i < got.size(); i += stride) {
    const std::size_t p = got[i].first;
    const std::size_t len = got[i].second;
    if (len == 0) continue;

    if (p + len > reference.size() || i + len > input.size() ||
        !std::equal(input.begin() + static_cast<std::ptrdiff_t>(i),
                    input.begin() + static_cast<std::ptrdiff_t>(i + len),
                    reference.begin() + static_cast<std::ptrdiff_t>(p))) {
      ++bad;
    }
  }

  return bad;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::cerr << "========================================\n"
              << "lazy matching statistics benchmark\n"
              << "========================================\n";

    const LazyArgs args = parse_lazy_args(argc, argv);

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

    Diagnostics::Line build_line{"build", true, {}};
    for (const auto& [name, ms] : lrf->buildPhases()) {
      build_line.values.emplace_back(name, ms);
    }
    std::cerr << "    lrf-ms: " << lrf_build_ms << " ms\n"
              << Diagnostics{{build_line}}.format();

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
    double v2_build_ms = 0.0;

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
      v2_build_ms = write_ms + load_ms;
      std::cerr << "    v2 table: write " << write_ms << " ms, load " << load_ms
                << " ms\n";
      variants.push_back(make_variant<msbench::V2Policy>("lazy-v2", v2));
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

    variants.push_back(make_variant<msbench::SassyPolicy>("lazy-sassy", sassy));

    std::ofstream csv;
    if (!args.results.empty()) {
      csv.open(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "file,input_bytes,non_acgt,lrf_ms";
      for (const Variant& variant : variants) {
        csv << ',' << variant.name << "_ms," << variant.name << "_breaks,"
            << variant.name << "_equal";
      }
      csv << '\n';
    }

    // ---------- Per file ----------

    // --quiet: the per-file lines go nowhere (a stream without a buffer
    // discards what is written to it); only a file counter is shown.
    std::ostream discard(nullptr);
    std::ostream& out = args.quiet ? discard : std::cerr;

    out << "[5] files\n";

    std::size_t total_bytes = 0;
    std::size_t total_non_acgt = 0;
    double total_lrf_ms = 0.0;
    ms2::PhraseStats total_phrases;  // of lrf-ms's matching statistics
    bool all_ok = true;

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

      out << "[" << f + 1 << "/" << files.size() << "] " << files[f]
                << "  (" << input.size() << " bytes, " << non_acgt
                << " non-ACGT, loaded in " << std::fixed
                << std::setprecision(1) << load_ms << " ms)\n";

      const auto step = [&out](const std::string& name) {
        out << "    " << std::left << std::setw(12) << name
                  << std::right << " ... " << std::flush;
      };

      step("lrf-ms");
      MatchingStatistics want;
      const double lrf_ms = msbench::time_ms(
          [&] { want = lrf->computeMatchingStatistics(input); });
      const ms2::PhraseStats phrases = ms2::phrase_stats(want);
      total_phrases.add(phrases);
      out << std::setw(10) << std::setprecision(2) << lrf_ms
                << " ms  " << phrases.format() << '\n';

      if (csv.is_open()) {
        csv << files[f] << ',' << input.size() << ',' << non_acgt << ','
            << std::setprecision(3) << lrf_ms;
      }

      for (Variant& variant : variants) {
        step(variant.name);
        ms2::LazyStats stats;
        MatchingStatistics got;
        const double ms = msbench::time_ms(
            [&] { got = variant.compute(input, reference, &stats); });

        std::size_t first = 0;
        const std::size_t mismatches = length_mismatches(want, got, first);
        const std::size_t bad = bad_positions(got, input, reference, args.sample);
        const bool ok = mismatches == 0 && bad == 0;

        out << std::setw(10) << std::setprecision(2) << ms << " ms  "
                  << std::setprecision(3)
                  << (ms == 0.0 ? 0.0 : lrf_ms / ms) << "x lrf-ms  "
                  << (ok ? "EQUAL" : "DIFFER");
        if (mismatches != 0) {
          out << " (" << mismatches << " lengths, first at " << first
                    << ": lrf-ms " << want[first].second << ", lazy "
                    << got[first].second << ")";
        }
        if (bad != 0) {
          out << " (" << bad << " sampled positions do not match)";
        }
        out << std::setprecision(2) << '\n' << stats.diagnostics().format();

        variant.total_ms += ms;
        variant.stats.add(stats);
        ++variant.files_compared;
        variant.files_equal += ok ? 1 : 0;
        variant.length_mismatches += mismatches;
        variant.bad_positions += bad;
        all_ok = all_ok && ok;

        if (csv.is_open()) {
          csv << ',' << std::setprecision(3) << ms << ',' << stats.breaks << ','
              << (ok ? "YES" : "NO");
        }
      }

      if (csv.is_open()) csv << '\n';

      total_bytes += input.size();
      total_non_acgt += non_acgt;
      total_lrf_ms += lrf_ms;
    }

    if (args.quiet) {
      std::cerr << "\rfiles: " << files.size() << "/" << files.size() << '\n';
    }

    // ---------- Totals ----------

    std::cerr << std::fixed;

    std::cerr << "[6] totals over " << files.size() << " files ("
              << total_bytes << " bytes, " << total_non_acgt << " non-ACGT)\n"
              << "    input: " << total_phrases.format() << '\n'
              << std::setprecision(2) << "    lrf-ms      " << std::setw(12)
              << total_lrf_ms << " ms  (build " << lrf_build_ms << " ms)\n";

    for (const Variant& variant : variants) {
      std::cerr << "    " << std::left << std::setw(12) << variant.name
                << std::right << std::setw(12) << variant.total_ms << " ms  "
                << std::setprecision(3)
                << (variant.total_ms == 0.0 ? 0.0
                                            : total_lrf_ms / variant.total_ms)
                << "x lrf-ms  equal " << variant.files_equal << "/"
                << variant.files_compared << " files\n"
                << std::setprecision(2) << variant.stats.diagnostics().format();
    }

    std::cerr << "    peak RSS " << msbench::peak_rss_mb() << " MB\n";

    // Machine-readable.
    std::cout << "files=" << files.size() << '\n'
              << "input_bytes=" << total_bytes << '\n'
              << "lrf_build_ms=" << lrf_build_ms << '\n'
              << "lrf_ms=" << total_lrf_ms << '\n'
              << "avg_match_length=" << total_phrases.average_match_length()
              << '\n'
              << "short_ms_percent=" << total_phrases.short_ms_percent() << '\n'
              << "avg_file_short_ms_percent="
              << total_phrases.average_file_short_ms_percent() << '\n'
              << "max_match_length=" << total_phrases.max_length << '\n'
              << "median_match_length=" << total_phrases.median_length() << '\n'
              << "unmatched_positions=" << total_phrases.unmatched << '\n'
              << "v2_build_ms=" << v2_build_ms << '\n'
              << "sassy_build_ms=" << sassy_write_ms + sassy_load_ms << '\n';
    for (const Variant& variant : variants) {
      std::cout << variant.name << "_ms=" << variant.total_ms << '\n'
                << variant.name << "_breaks=" << variant.stats.breaks << '\n'
                << variant.name << "_extended=" << variant.stats.extended << '\n'
                << variant.name << "_length_mismatches="
                << variant.length_mismatches << '\n';
    }
    std::cout << "all_ok=" << (all_ok ? "YES" : "NO") << '\n';

    return all_ok ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
