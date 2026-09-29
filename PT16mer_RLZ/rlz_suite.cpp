// RLZ parsing benchmark suite: every parser on the same inputs, single-threaded,
// with comparable timings.
//
//   sa-binary-search  the baseline: rlz::lzFactorize (parser.hpp), binary
//             search over the suffix array
//   lrf-ms    a second baseline: the matching statistics by lrf-ms
//             (lrf_ms/lrf_ms.hpp; its ISA, LCP, LRF and RMQ are built once),
//             then one sweep turning them into the parse: the phrase at i is
//             MS[i] long (a literal when MS[i] <= 1), the next starts after it
//   pt16      the original PT16 table (pt16_build.hpp, pt16_rlz.hpp)
//   pt16-v2   the v2 table (pt16_build_v2.hpp, pt16_rlz_v2.hpp)
//   sassy     the self-contained sassy table (variants/pt16_*sassy*.hpp)
//   powered   powered backward search in an FM-index (../powered_rlz), only
//             when built with -DWITH_POWERED (x86-64, GCC) and given
//             --powered-index
//
// What is timed, for every parser alike: parsing one input that is already in
// memory into an in-memory list of phrases. Not timed: loading inputs, building
// or loading tables and indexes (reported separately), and the checks.
//
// Checks (untimed), for every parser alike: its phrases are decoded against the
// reference and must reproduce the input exactly. lrf-ms and the PT16 parsers
// must also give the baseline's phrase count (they compute the same greedy
// left-to-right parse). powered parses right to left (longest suffix of the
// remaining prefix), so its phrases may differ; its count is reported.
//
// Inputs and reference should be plain ACGT (powered maps any other byte to
// one of A/C/G/T); non-ACGT bytes are reported.
//
// Build (from PT16mer_RLZ/):
//   without powered (any platform):
//     g++ -std=c++20 -O3 rlz_suite.cpp -o rlz_suite
//   with powered (x86-64, GCC; block sizes as powered_rlz was built with):
//     g++ -std=c++2a -O3 -march=native -DNDEBUG -DWITH_POWERED \
//         -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 rlz_suite.cpp -o rlz_suite
//
// Run:
//   ./rlz_suite --reference REF --suffix-array REF.sa --filenames LIST
//               [--powered-index REF_four.bwt] [--table PATH] [--results CSV]
//               [--max-files N] [--quiet]

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "lrf_ms/lrf_ms.hpp"  // LRFMS (the lrf-ms baseline)
#include "parser.hpp"      // lzFactorize (baseline)
#include "pt16_utils.hpp"  // loaders, time_ms, shared constants
#include "rlz_common.hpp"

// The original PT16 format defines the same names as v2 (PT16RLZParser,
// build_pt16_table, ScanResult, ...), so it lives in its own namespace. Every
// header it includes is already included above, so nothing else ends up in it.
namespace pt16_v1 {
#include "pt16_build.hpp"
#include "pt16_rlz.hpp"
}  // namespace pt16_v1

#include "pt16_build_v2.hpp"
#include "pt16_rlz_v2.hpp"
#include "variants/pt16_build_sassy.hpp"
#include "variants/pt16_sassy.hpp"

#ifdef WITH_POWERED
#include "../powered_rlz/include/types.hpp"  // bbwt::non_rle
#endif

using Symbol = unsigned char;
using SAType = std::uint32_t;
namespace fs = std::filesystem;

namespace {

// ---------- Arguments ----------

struct SuiteArgs {
  std::string reference;
  std::string suffix_array;
  std::string filenames;
  std::string table;          // PT16 table files; default beside the reference
  std::string powered_index;  // powered_rlz index (.bwt from transform -b)
  std::string results;        // optional CSV
  std::size_t max_files = 0;
  bool quiet = false;
};

std::string require_value(int& i, int argc, char** argv) {
  if (i + 1 >= argc) {
    throw std::runtime_error(std::string("missing value for ") + argv[i]);
  }
  return argv[++i];
}

void print_usage(const char* program) {
  std::cout << "Usage: " << program
            << " --reference PATH --suffix-array PATH --filenames PATH\n"
               "  [--powered-index PATH]  powered_rlz index (REF_four.bwt)"
#ifndef WITH_POWERED
               " -- this build has no powered support"
#endif
               "\n"
               "  [--table PATH]          PT16 table files (default: "
               "<reference>.suite_pt16, .v2, .sassy)\n"
               "  [--results PATH]        per-file CSV\n"
               "  [--max-files N]         only the first N input files\n"
               "  [--quiet]               no per-file lines, only the totals\n";
}

SuiteArgs parse_suite_args(int argc, char** argv) {
  SuiteArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--reference") {
      args.reference = require_value(i, argc, argv);
    } else if (option == "--suffix-array") {
      args.suffix_array = require_value(i, argc, argv);
    } else if (option == "--filenames") {
      args.filenames = require_value(i, argc, argv);
    } else if (option == "--table") {
      args.table = require_value(i, argc, argv);
    } else if (option == "--powered-index") {
      args.powered_index = require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = std::stoull(require_value(i, argc, argv));
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
#ifndef WITH_POWERED
  if (!args.powered_index.empty()) {
    throw std::runtime_error(
        "--powered-index given, but this build has no powered support "
        "(rebuild with -DWITH_POWERED on x86-64)");
  }
#endif
  if (args.table.empty()) args.table = args.reference + ".suite_pt16";
  return args;
}

// ---------- Checks ----------

// Our parsers' phrases: (input position, reference position, length); a
// length-1 phrase is a literal whose "reference position" is the character.
// Returns the first position where the decoded text differs from the input,
// or input.size() when they are equal.
using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;

std::size_t decode_triples(const Triples& phrases,
                           const std::vector<Symbol>& reference,
                           const std::vector<Symbol>& input) {
  std::size_t at = 0;
  for (const auto& [i, position, length] : phrases) {
    if (i != at) return std::min(at, input.size());
    for (std::size_t j = 0; j < length; ++j, ++at) {
      const Symbol c = length == 1
                           ? static_cast<Symbol>(position)
                           : (position + j < reference.size()
                                  ? reference[position + j]
                                  : static_cast<Symbol>(0));
      if (at >= input.size() || input[at] != c) return std::min(at, input.size());
    }
  }
  return at == input.size() ? input.size() : std::min(at, input.size());
}

// ---------- Parsers ----------

// One parser: builds or loads what it needs once, then parses inputs.
struct Parser {
  virtual ~Parser() = default;
  virtual std::string name() const = 0;
  // Time to build its table (0 if none) and to load it.
  double build_ms = 0.0;
  double load_ms = 0.0;
  // Memory: its own structure (table or index), and whether parsing also
  // reads the suffix array and/or the reference.
  std::size_t own_bytes = 0;
  // Memory allocated beyond that, if the structure reserves more than it
  // fills (powered's super blocks); 0 otherwise.
  std::size_t reserved_bytes = 0;
  bool needs_suffix_array = false;
  bool needs_reference = false;
  std::size_t total_bytes(std::size_t n) const {
    return own_bytes + (needs_suffix_array ? n * sizeof(SAType) : 0) +
           (needs_reference ? n * sizeof(Symbol) : 0);
  }
  // Parses `input` (the only timed part, see parse_timed), keeping the result.
  virtual void parse(const std::vector<Symbol>& input) = 0;
  // Frees the previous result, so that every parser starts from an empty
  // phrase list (as powered_rlz's own rlz_parser does per file) and no
  // parser's timing includes freeing the previous file's result.
  virtual void release() = 0;
  virtual std::size_t phrase_count() const = 0;
  // First differing position of the decoded parse, or input.size().
  virtual std::size_t decode(const std::vector<Symbol>& reference,
                             const std::vector<Symbol>& input) const = 0;
  // True for parsers computing the same greedy left-to-right parse as the
  // baseline (their phrase count must match).
  virtual bool left_to_right() const { return true; }

  double parse_timed(const std::vector<Symbol>& input) {
    release();  // untimed
    return time_ms([&] { parse(input); });
  }
};

struct BaselineParser : Parser {
  const std::vector<Symbol>& reference;
  const std::vector<SAType>& suffix_array;
  Triples result;

  BaselineParser(const std::vector<Symbol>& ref, const std::vector<SAType>& sa)
      : reference(ref), suffix_array(sa) {
    needs_suffix_array = true;
    needs_reference = true;
  }

  std::string name() const override { return "sa-binary-search"; }
  void parse(const std::vector<Symbol>& input) override {
    result = lzFactorize<Symbol, SAType>(input, reference, suffix_array);
  }
  void release() override {
    std::remove_reference_t<decltype(result)>().swap(result);
  }
  std::size_t phrase_count() const override { return result.size(); }
  std::size_t decode(const std::vector<Symbol>& ref,
                     const std::vector<Symbol>& input) const override {
    return decode_triples(result, ref, input);
  }
};

// The lrf-ms baseline: matching statistics, then the greedy parse read off
// them in one sweep.
struct LrfMsParser : Parser {
  std::unique_ptr<LRFMS<Symbol, SAType>> lrf;
  Triples result;

  LrfMsParser(const std::vector<Symbol>& ref, const std::vector<SAType>& sa) {
    // Building its auxiliary structures (ISA, LCP, LRF, RMQ over LCP).
    build_ms = time_ms(
        [&] { lrf = std::make_unique<LRFMS<Symbol, SAType>>(ref, sa); });
    // ISA, LCP and LRF: 4 bytes per position each; the RMQ is a tree over
    // blocks of 2^7 LCP entries (2 * size ints, size the next power of two
    // with size * 128 >= n; see rmq_tree.h).
    std::size_t rmq_size = 1;
    while ((rmq_size << 7) < ref.size()) rmq_size <<= 1;
    own_bytes = (lrf->isa().size() + lrf->lcp().size() + lrf->lrf().size()) *
                    sizeof(std::int32_t) +
                2 * rmq_size * sizeof(int);
    needs_suffix_array = true;
    needs_reference = true;
  }
  std::string name() const override { return "lrf-ms"; }
  void release() override {
    std::remove_reference_t<decltype(result)>().swap(result);
  }
  void parse(const std::vector<Symbol>& input) override {
    // Both steps are the parse: the matching statistics, then the sweep.
    const MatchingStatistics ms = lrf->computeMatchingStatistics(input);
    std::size_t i = 0;
    while (i < input.size()) {
      const auto [position, length] = ms[i];
      if (length <= 1) {
        result.emplace_back(i, static_cast<std::size_t>(input[i]), 1);
        ++i;
      } else {
        result.emplace_back(i, position, length);
        i += length;
      }
    }
  }
  std::size_t phrase_count() const override { return result.size(); }
  std::size_t decode(const std::vector<Symbol>& ref,
                     const std::vector<Symbol>& input) const override {
    return decode_triples(result, ref, input);
  }
};

struct Pt16Parser : Parser {
  std::unique_ptr<pt16_v1::PT16RLZParser<Symbol, SAType>> parser;
  Triples result;

  Pt16Parser(const std::vector<Symbol>& ref, const std::vector<SAType>& sa,
             const std::string& path) {
    build_ms = time_ms([&] {
      fs::remove(path);
      pt16_v1::build_pt16_table(ref, sa, path);
    });
    load_ms = time_ms([&] {
      parser = std::make_unique<pt16_v1::PT16RLZParser<Symbol, SAType>>(ref, sa,
                                                                        path);
    });
    own_bytes = parser->stats().approx_bytes;
    needs_suffix_array = true;
    needs_reference = true;
  }
  std::string name() const override { return "pt16"; }
  void parse(const std::vector<Symbol>& input) override {
    result = parser->lzFactorize(input);
  }
  void release() override {
    std::remove_reference_t<decltype(result)>().swap(result);
  }
  std::size_t phrase_count() const override { return result.size(); }
  std::size_t decode(const std::vector<Symbol>& ref,
                     const std::vector<Symbol>& input) const override {
    return decode_triples(result, ref, input);
  }
};

struct V2Parser : Parser {
  std::unique_ptr<PT16RLZParser<Symbol, SAType>> parser;
  Triples result;

  V2Parser(const std::vector<Symbol>& ref, const std::vector<SAType>& sa,
           const std::string& path) {
    build_ms = time_ms([&] {
      fs::remove(path);
      build_pt16_table(ref, sa, path);
    });
    load_ms = time_ms([&] {
      parser = std::make_unique<PT16RLZParser<Symbol, SAType>>(ref, sa, path);
    });
    own_bytes = parser->stats().approx_bytes;
    needs_suffix_array = true;
    needs_reference = true;
  }
  std::string name() const override { return "pt16-v2"; }
  void parse(const std::vector<Symbol>& input) override {
    result = parser->lzFactorize(input);
  }
  void release() override {
    std::remove_reference_t<decltype(result)>().swap(result);
  }
  std::size_t phrase_count() const override { return result.size(); }
  std::size_t decode(const std::vector<Symbol>& ref,
                     const std::vector<Symbol>& input) const override {
    return decode_triples(result, ref, input);
  }
};

struct SassyParser : Parser {
  const std::vector<Symbol>& reference;
  std::unique_ptr<PT16SassyLookup> lookup;
  Triples result;

  SassyParser(const std::vector<Symbol>& ref, const std::vector<SAType>& sa,
              const std::string& path)
      : reference(ref) {
    build_ms = time_ms([&] {
      fs::remove(path);
      build_pt16_sassy_table(ref, sa, path);
    });
    load_ms = time_ms([&] { lookup = std::make_unique<PT16SassyLookup>(path); });
    own_bytes = lookup->stats().approx_bytes;
    needs_reference = true;  // to extend matches; no suffix array
  }
  std::string name() const override { return "sassy"; }
  void parse(const std::vector<Symbol>& input) override {
    result = lookup->lzFactorize(input, reference);
  }
  void release() override {
    std::remove_reference_t<decltype(result)>().swap(result);
  }
  std::size_t phrase_count() const override { return result.size(); }
  std::size_t decode(const std::vector<Symbol>& ref,
                     const std::vector<Symbol>& input) const override {
    return decode_triples(result, ref, input);
  }
};

#ifdef WITH_POWERED
// powered_rlz's parser (right to left, phrases (length, position) into the
// cyclic reference, kept last-to-first). Its index is built beforehand with
// powered_rlz's tools; here it is only loaded.
struct PoweredParser : Parser {
  std::unique_ptr<bbwt::non_rle<>> index;
  std::vector<std::tuple<std::uint64_t, std::uint64_t>> result;
  static constexpr std::uint8_t code_size = 4;  // the powered index

  explicit PoweredParser(const std::string& path) {
    // powered_rlz's loader only prints " -> Failed" and exits on a missing
    // file, so check both files first and name them.
    const std::size_t dot = path.find_last_of('.');
    const std::string data_file =
        dot == std::string::npos
            ? path + "_data"
            : path.substr(0, dot) + "_data" + path.substr(dot);
    for (const std::string& file : {path, data_file}) {
      if (!fs::exists(file)) {
        throw std::runtime_error("powered index file not found: " + file);
      }
    }

    load_ms = time_ms([&] { index = std::make_unique<bbwt::non_rle<>>(path); });

    // Its data: the header's data section and character counts (from the
    // index file), the GCA (held as 8-byte integers), and the super-block
    // data (the _data file, which the loader reads in full). bytes() counts
    // each super block's allocated capacity instead (read_super_block
    // mallocs `in_bytes` but fills only the stored words), so it is reported
    // apart as reserved memory.
    std::uint64_t data_bytes = 0;
    {
      std::ifstream header(path, std::ios::binary);
      header.read(reinterpret_cast<char*>(&data_bytes), sizeof(data_bytes));
    }
    own_bytes = sizeof(bbwt::non_rle<>) + data_bytes +
                257 * sizeof(std::uint64_t) +
                index->gca_.size() * sizeof(std::uint64_t) +
                fs::file_size(data_file);
    const std::size_t allocated =
        index->bytes() + index->gca_.size() * sizeof(std::uint64_t);
    reserved_bytes = allocated > own_bytes ? allocated - own_bytes : 0;
  }
  std::string name() const override { return "powered"; }
  bool left_to_right() const override { return false; }
  void parse(const std::vector<Symbol>& input) override {
    const std::string_view view(reinterpret_cast<const char*>(input.data()),
                                input.size());
    index->parse_tuples(view, result, code_size);
  }
  void release() override {
    std::remove_reference_t<decltype(result)>().swap(result);
  }
  std::size_t phrase_count() const override { return result.size(); }
  std::size_t decode(const std::vector<Symbol>& ref,
                     const std::vector<Symbol>& input) const override {
    const std::uint64_t literal = std::uint64_t{1} << 63;
    const std::size_t m = ref.size();
    std::size_t at = 0;
    for (auto it = result.rbegin(); it != result.rend(); ++it) {
      std::uint64_t length = std::get<0>(*it);
      const std::uint64_t position = std::get<1>(*it);
      const bool is_literal = (length & literal) != 0;
      if (is_literal) length -= literal;
      for (std::uint64_t j = 0; j < length; ++j, ++at) {
        const Symbol c = is_literal ? static_cast<Symbol>(position)
                                    : ref[(position + j) % m];
        if (at >= input.size() || input[at] != c) {
          return std::min(at, input.size());
        }
      }
    }
    return at == input.size() ? input.size() : std::min(at, input.size());
  }
};
#endif

double mb(const std::size_t bytes) {
  return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

// ---------- Totals ----------

struct Totals {
  double parse_ms = 0.0;
  std::size_t phrases = 0;
  std::size_t files_ok = 0;         // decoded exactly
  std::size_t count_mismatches = 0;  // phrase count differs from baseline
};

}  // namespace

int main(int argc, char** argv) {
  try {
    std::cerr << "========================================\n"
              << "RLZ parsing suite (single-threaded)\n"
              << "========================================\n";

    const SuiteArgs args = parse_suite_args(argc, argv);

    std::vector<Symbol> reference;
    const double reference_ms =
        time_ms([&] { reference = load_reference<Symbol>(args.reference); });
    std::cerr << "[1] reference: " << reference.size() << " bytes in "
              << reference_ms << " ms\n";

    std::size_t reference_non_acgt = 0;
    for (const Symbol c : reference) reference_non_acgt += is_acgt(c) ? 0 : 1;
    if (reference_non_acgt != 0) {
      std::cerr << "    WARNING: " << reference_non_acgt
                << " non-ACGT bytes in the reference\n";
    }

    std::vector<SAType> suffix_array;
    const double sa_ms = time_ms(
        [&] { suffix_array = load_suffix_array<SAType>(args.suffix_array); });
    std::cerr << "[2] suffix array: " << suffix_array.size() << " entries in "
              << sa_ms << " ms\n";
    if (suffix_array.size() != reference.size()) {
      throw std::runtime_error("suffix array and reference sizes differ");
    }

    std::vector<std::string> files = load_input_list(args.filenames);
    if (args.max_files != 0 && files.size() > args.max_files) {
      files.resize(args.max_files);
    }
    std::cerr << "[3] input files: " << files.size() << '\n';

    // ---------- Build / load every parser once ----------

    std::cerr << "[4] building and loading\n";
    std::vector<std::unique_ptr<Parser>> parsers;
    parsers.push_back(std::make_unique<BaselineParser>(reference, suffix_array));
    parsers.push_back(std::make_unique<LrfMsParser>(reference, suffix_array));
    parsers.push_back(
        std::make_unique<Pt16Parser>(reference, suffix_array, args.table));
    parsers.push_back(
        std::make_unique<V2Parser>(reference, suffix_array, args.table + ".v2"));
    parsers.push_back(std::make_unique<SassyParser>(reference, suffix_array,
                                                    args.table + ".sassy"));
#ifdef WITH_POWERED
    if (!args.powered_index.empty()) {
      parsers.push_back(std::make_unique<PoweredParser>(args.powered_index));
    }
#endif

    for (const auto& parser : parsers) {
      std::cerr << "    " << std::left << std::setw(17) << parser->name()
                << std::right << " build " << std::fixed << std::setprecision(1)
                << std::setw(9) << parser->build_ms << " ms, load "
                << std::setw(8) << parser->load_ms << " ms, own "
                << std::setw(8) << mb(parser->own_bytes) << " MB, needs "
                << std::setw(8) << mb(parser->total_bytes(reference.size()))
                << " MB"
                << (parser->reserved_bytes != 0
                        ? " [+" + std::to_string(static_cast<long long>(
                                      mb(parser->reserved_bytes))) +
                              " MB reserved, unused]"
                        : std::string())
                << (parser->needs_suffix_array || parser->needs_reference
                        ? std::string(" (own") +
                              (parser->needs_suffix_array ? " + SA" : "") +
                              (parser->needs_reference ? " + reference" : "") +
                              ")"
                        : std::string(" (own only)"))
                << '\n';
    }
    if (args.powered_index.empty()) {
      std::cerr << "    (powered: not run; give --powered-index"
#ifndef WITH_POWERED
                   " and build with -DWITH_POWERED"
#endif
                   ")\n";
    }

    std::ofstream csv;
    if (!args.results.empty()) {
      csv.open(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "file,input_bytes,non_acgt";
      for (const auto& parser : parsers) {
        csv << ',' << parser->name() << "_ms," << parser->name()
            << "_phrases," << parser->name() << "_ok";
      }
      csv << '\n';
    }

    // ---------- Per file ----------

    std::ostream discard(nullptr);
    std::ostream& out = args.quiet ? discard : std::cerr;
    out << "[5] files\n";

    std::vector<Totals> totals(parsers.size());
    std::size_t total_bytes = 0;

    for (std::size_t f = 0; f < files.size(); ++f) {
      if (args.quiet) {
        std::cerr << "\rfiles: " << f << "/" << files.size() << std::flush;
      }

      const std::vector<Symbol> input = load_input<Symbol>(files[f]);
      std::size_t non_acgt = 0;
      for (const Symbol c : input) non_acgt += is_acgt(c) ? 0 : 1;
      total_bytes += input.size();

      out << "[" << f + 1 << "/" << files.size() << "] " << files[f] << "  ("
          << input.size() << " bytes";
      if (non_acgt != 0) out << ", " << non_acgt << " non-ACGT";
      out << ")\n";

      if (csv.is_open()) {
        csv << files[f] << ',' << input.size() << ',' << non_acgt;
      }

      std::size_t baseline_phrases = 0;
      double baseline_ms = 0.0;

      for (std::size_t k = 0; k < parsers.size(); ++k) {
        Parser& parser = *parsers[k];

        // The only timed part: parsing the in-memory input.
        const double ms = parser.parse_timed(input);

        // Checks, untimed.
        const std::size_t difference = parser.decode(reference, input);
        const bool decoded = difference == input.size();
        const std::size_t phrases = parser.phrase_count();
        if (k == 0) {
          baseline_phrases = phrases;
          baseline_ms = ms;
        }
        const bool count_ok = !parser.left_to_right() || phrases == baseline_phrases;

        out << "    " << std::left << std::setw(17) << parser.name()
            << std::right << std::fixed << std::setprecision(2) << std::setw(10)
            << ms << " ms";
        if (k > 0) {
          out << "  " << std::setprecision(2) << std::setw(6)
              << (ms == 0.0 ? 0.0 : baseline_ms / ms) << "x";
        } else {
          out << "         ";
        }
        out << "  " << std::setw(9) << phrases << " phrases  "
            << (decoded ? "decode OK" : "DECODE DIFFERS");
        if (!decoded) out << " (first at " << difference << ")";
        if (!count_ok) out << "  PHRASE COUNT DIFFERS from baseline";
        out << '\n';

        Totals& t = totals[k];
        t.parse_ms += ms;
        t.phrases += phrases;
        t.files_ok += decoded ? 1 : 0;
        t.count_mismatches += count_ok ? 0 : 1;

        if (csv.is_open()) {
          csv << ',' << std::setprecision(3) << ms << ',' << phrases << ','
              << (decoded && count_ok ? "YES" : "NO");
        }
      }

      if (csv.is_open()) csv << '\n';
    }

    if (args.quiet) {
      std::cerr << "\rfiles: " << files.size() << "/" << files.size() << '\n';
    }

    // ---------- Totals ----------

    std::cerr << std::fixed << std::setprecision(2) << "[6] totals over "
              << files.size() << " files (" << total_bytes << " bytes)\n"
              << "    parser              build ms    load ms    parse ms  speedup"
                 "      phrases  decoded  count        own MB  needs MB\n";
    bool all_ok = true;
    for (std::size_t k = 0; k < parsers.size(); ++k) {
      const Parser& parser = *parsers[k];
      const Totals& t = totals[k];
      const bool ok = t.files_ok == files.size() && t.count_mismatches == 0;
      all_ok = all_ok && ok;
      std::cerr << "    " << std::left << std::setw(17) << parser.name()
                << std::right << std::setw(11) << parser.build_ms
                << std::setw(11) << parser.load_ms << std::setw(12)
                << t.parse_ms << std::setw(9)
                << (t.parse_ms == 0.0 ? 0.0 : totals[0].parse_ms / t.parse_ms)
                << std::setw(13) << t.phrases << std::setw(6) << t.files_ok
                << "/" << files.size() << "  "
                << std::left << std::setw(12)
                << (parser.left_to_right()
                        ? (t.count_mismatches == 0 ? "same" : "DIFFERS")
                        : "(own parse)")
                << std::right << std::setw(9) << mb(parser.own_bytes)
                << std::setw(10) << mb(parser.total_bytes(reference.size()))
                << '\n';
    }
    std::cerr << "    peak RSS " << peak_rss_mb() << " MB\n";

    std::cout << "files=" << files.size() << '\n'
              << "input_bytes=" << total_bytes << '\n';
    for (std::size_t k = 0; k < parsers.size(); ++k) {
      const std::string n = parsers[k]->name();
      std::cout << n << "_own_bytes=" << parsers[k]->own_bytes << '\n'
                << n << "_reserved_bytes=" << parsers[k]->reserved_bytes << '\n'
                << n << "_needed_bytes="
                << parsers[k]->total_bytes(reference.size()) << '\n'
                << n << "_build_ms=" << parsers[k]->build_ms << '\n'
                << n << "_load_ms=" << parsers[k]->load_ms << '\n'
                << n << "_parse_ms=" << totals[k].parse_ms << '\n'
                << n << "_phrases=" << totals[k].phrases << '\n'
                << n << "_files_ok=" << totals[k].files_ok << '\n';
    }
    std::cout << "all_ok=" << (all_ok ? "YES" : "NO") << '\n';

    return all_ok ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
