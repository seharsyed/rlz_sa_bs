#pragma once

#include <sys/resource.h>

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"  // time_ms, MatchingStatistics

/**
 * Utilities for the matching-statistics benchmark: arguments, loading,
 * timing, result checking and output files.
 *
 * Matching statistics are stored as (reference position, match length).
 * The position of a match is not unique, so the implementations are
 * compared on lengths. Positions are checked separately, by confirming
 * that the reported occurrence really matches the input.
 */

namespace msbench {

namespace fs = std::filesystem;

// ---------- Arguments ----------

struct Args {
  std::string reference;
  std::string suffix_array;
  std::string filenames;
  std::string results;
  std::string checksums;
  std::string pt16_table;
  std::string dump_directory;

  std::size_t repeats = 1;

  bool check_invariants = true;

  // Position verification: every position, or a random sample.
  bool verify_full = false;
  std::size_t sample = 10000;
  std::uint64_t seed = 1;

  // Also prove that no longer match exists (binary search in the SA).
  bool verify_maximality = false;

  bool stop_on_mismatch = false;
};

inline void print_usage(const char* program) {
  std::cout
      << "Usage: " << program << " --reference PATH --suffix-array PATH"
      << " --filenames PATH --results PATH\n"
      << "  [--checksums PATH]        default: <results>.checksums.csv\n"
      << "  [--pt16-table PATH]       default: next to the reference\n"
      << "  [--dump-dir DIR]          write every implementation's output\n"
      << "  [--repeats N]             timed runs per file, min is reported\n"
      << "  [--no-invariants]         skip the structural checks\n"
      << "  [--verify-full]           check every position's occurrence\n"
      << "  [--sample N --seed S]     check N random positions (default "
         "10000)\n"
      << "  [--verify-maximality]     prove the length cannot be extended\n"
      << "  [--stop-on-mismatch]      stop after the first divergent file\n";
}

inline std::string require_value(int& i, int argc, char** argv) {
  if (i + 1 >= argc) {
    throw std::runtime_error(std::string("missing value for ") + argv[i]);
  }

  return argv[++i];
}

inline std::uint64_t parse_number(const std::string& option,
                                  const std::string& text) {
  try {
    std::size_t used = 0;
    const std::uint64_t value = std::stoull(text, &used);

    if (used == text.size()) {
      return value;
    }
  } catch (const std::exception&) {
  }

  throw std::runtime_error("invalid number for " + option + ": " + text);
}

inline Args parse_args(int argc, char** argv) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];

    if (option == "--reference") {
      args.reference = require_value(i, argc, argv);
    } else if (option == "--suffix-array") {
      args.suffix_array = require_value(i, argc, argv);
    } else if (option == "--filenames") {
      args.filenames = require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = require_value(i, argc, argv);
    } else if (option == "--checksums") {
      args.checksums = require_value(i, argc, argv);
    } else if (option == "--pt16-table") {
      args.pt16_table = require_value(i, argc, argv);
    } else if (option == "--dump-dir") {
      args.dump_directory = require_value(i, argc, argv);
    } else if (option == "--repeats") {
      args.repeats = parse_number(option, require_value(i, argc, argv));
    } else if (option == "--no-invariants") {
      args.check_invariants = false;
    } else if (option == "--verify-full") {
      args.verify_full = true;
    } else if (option == "--sample") {
      args.sample = parse_number(option, require_value(i, argc, argv));
    } else if (option == "--seed") {
      args.seed = parse_number(option, require_value(i, argc, argv));
    } else if (option == "--verify-maximality") {
      args.verify_maximality = true;
    } else if (option == "--stop-on-mismatch") {
      args.stop_on_mismatch = true;
    } else if (option == "--help" || option == "-h") {
      print_usage(argv[0]);
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }

  if (args.reference.empty() || args.suffix_array.empty() ||
      args.filenames.empty() || args.results.empty()) {
    throw std::runtime_error(
        "reference, suffix-array, filenames and results are required");
  }

  if (args.repeats == 0) {
    throw std::runtime_error("--repeats must be at least 1");
  }

  if (args.checksums.empty()) {
    args.checksums = args.results + ".checksums.csv";
  }

  return args;
}

// PT16 table location: explicit, or "<dataset>_pt16_hl.bin" beside the
// reference, as in the RLZ benchmark.
inline std::string pt16_table_path(const Args& args) {
  if (!args.pt16_table.empty()) {
    return args.pt16_table;
  }

  const fs::path reference_path(args.reference);
  const std::string dataset = reference_path.parent_path().filename().string();

  return (reference_path.parent_path() / (dataset + "_pt16_hl.bin")).string();
}

// ---------- Loading ----------

template <typename T>
inline std::vector<T> read_binary(const std::string& path) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);

  if (!input) {
    throw std::runtime_error("cannot open file: " + path);
  }

  const std::streamoff bytes = input.tellg();
  input.seekg(0);

  std::vector<T> values(static_cast<std::size_t>(bytes) / sizeof(T));

  input.read(reinterpret_cast<char*>(values.data()),
             static_cast<std::streamsize>(values.size() * sizeof(T)));

  if (!input) {
    throw std::runtime_error("failed while reading file: " + path);
  }

  return values;
}

template <typename Symbol>
inline std::vector<Symbol> load_reference(const std::string& path) {
  return read_binary<Symbol>(path);
}

template <typename SAType>
inline std::vector<SAType> load_suffix_array(const std::string& path) {
  return read_binary<SAType>(path);
}

template <typename Symbol>
inline std::vector<Symbol> load_input(const std::string& path) {
  return read_binary<Symbol>(path);
}

inline std::vector<std::string> load_input_list(const std::string& path) {
  std::ifstream input(path);

  if (!input) {
    throw std::runtime_error("cannot open input list: " + path);
  }

  std::vector<std::string> files;
  std::string filename;

  while (std::getline(input, filename)) {
    if (!filename.empty() && filename.back() == '\r') {
      filename.pop_back();
    }

    if (!filename.empty()) {
      files.push_back(filename);
    }
  }

  return files;
}

// ---------- Alphabet ----------

// The set of symbols that occur in the reference.
template <typename Symbol>
class SymbolTable {
 public:
  explicit SymbolTable(const std::vector<Symbol>& reference) {
    for (const Symbol symbol : reference) {
      present_.set(static_cast<unsigned char>(symbol));
    }
  }

  bool contains(Symbol symbol) const {
    return present_.test(static_cast<unsigned char>(symbol));
  }

  std::size_t distinct() const { return present_.count(); }

 private:
  std::bitset<256> present_;
};

// ---------- Timing ----------

// time_ms is shared with the RLZ benchmark (pt16_utils.hpp).
using ::time_ms;

struct Timing {
  double first_ms = 0.0;
  double min_ms = 0.0;
  double total_ms = 0.0;
  std::size_t repeats = 0;
};

// Runs fn `repeats` times. The first run is reported separately because it
// includes cold caches; the minimum is the headline figure.
template <typename Fn>
inline Timing time_repeated(std::size_t repeats, Fn&& fn) {
  Timing timing;
  timing.repeats = std::max<std::size_t>(repeats, 1);

  for (std::size_t run = 0; run < timing.repeats; ++run) {
    const double elapsed = msbench::time_ms(fn);

    if (run == 0) {
      timing.first_ms = elapsed;
      timing.min_ms = elapsed;
    } else {
      timing.min_ms = std::min(timing.min_ms, elapsed);
    }

    timing.total_ms += elapsed;
  }

  return timing;
}

// peak_rss_mb is shared with the RLZ benchmark (pt16_utils.hpp).
using ::peak_rss_mb;

// ---------- Validation results ----------

struct Validation {
  bool ok = true;
  std::string message;

  static Validation fail(std::string message) {
    Validation validation;
    validation.ok = false;
    validation.message = std::move(message);
    return validation;
  }
};

inline std::ostream& operator<<(std::ostream& out,
                                const Validation& validation) {
  return out << (validation.ok ? std::string("ok") : validation.message);
}

// ---------- Digest ----------

inline std::string to_hex(std::uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << value;
  return out.str();
}

struct Digest {
  std::size_t entries = 0;
  std::uint64_t total_len = 0;
  std::uint32_t max_len = 0;
  std::size_t zero_len_count = 0;

  // Depends on the lengths only.
  std::uint64_t len_hash = 0;

  // Depends on the positions of non-empty matches. Positions are not
  // canonical, so two correct implementations may differ here.
  std::uint64_t pos_hash = 0;
};

inline void fnv_mix(std::uint64_t& hash, std::uint32_t value) {
  for (int byte = 0; byte < 4; ++byte) {
    hash ^= (value >> (8 * byte)) & 0xFFU;
    hash *= 1099511628211ULL;
  }
}

inline Digest digest_ms(const MatchingStatistics& ms) {
  Digest digest;
  digest.entries = ms.size();
  digest.len_hash = 1469598103934665603ULL;
  digest.pos_hash = 1469598103934665603ULL;

  for (const auto& entry : ms) {
    const std::uint32_t pos = entry.first;
    const std::uint32_t len = entry.second;

    digest.total_len += len;
    digest.max_len = std::max(digest.max_len, len);

    fnv_mix(digest.len_hash, len);

    if (len == 0) {
      ++digest.zero_len_count;
    } else {
      fnv_mix(digest.pos_hash, pos);
    }
  }

  return digest;
}

// ---------- Comparison with the baseline ----------

inline std::vector<std::uint32_t> extract_lengths(
    const MatchingStatistics& ms) {
  std::vector<std::uint32_t> lengths;
  lengths.reserve(ms.size());

  for (const auto& entry : ms) {
    lengths.push_back(entry.second);
  }

  return lengths;
}

struct LengthComparison {
  bool checked = false;
  bool equal = false;

  std::size_t baseline_count = 0;
  std::size_t variant_count = 0;

  std::size_t first_divergence = 0;
  std::size_t divergent_positions = 0;
  std::uint32_t baseline_value = 0;
  std::uint32_t variant_value = 0;

  std::string describe() const {
    std::ostringstream out;

    if (baseline_count != variant_count) {
      out << "entry counts differ (baseline " << baseline_count << ", variant "
          << variant_count << ")";

      if (divergent_positions > 0) {
        out << "; ";
      }
    }

    if (divergent_positions > 0) {
      out << divergent_positions << " divergent positions, first at "
          << first_divergence << " (baseline " << baseline_value << ", variant "
          << variant_value << ")";
    }

    return out.str();
  }
};

inline LengthComparison compare_lengths(
    const std::vector<std::uint32_t>& baseline, const MatchingStatistics& ms) {
  LengthComparison comparison;
  comparison.checked = true;
  comparison.baseline_count = baseline.size();
  comparison.variant_count = ms.size();

  const std::size_t common = std::min(baseline.size(), ms.size());

  for (std::size_t i = 0; i < common; ++i) {
    if (baseline[i] != ms[i].second) {
      if (comparison.divergent_positions == 0) {
        comparison.first_divergence = i;
        comparison.baseline_value = baseline[i];
        comparison.variant_value = ms[i].second;
      }

      ++comparison.divergent_positions;
    }
  }

  comparison.equal =
      comparison.divergent_positions == 0 && baseline.size() == ms.size();

  return comparison;
}

// ---------- Checks on a single implementation's output ----------

// Structural properties every correct matching-statistics array has.
template <typename Symbol>
inline Validation validate_invariants(const MatchingStatistics& ms,
                                      const std::vector<Symbol>& input,
                                      const std::vector<Symbol>& reference,
                                      const SymbolTable<Symbol>& alphabet) {
  if (ms.size() != input.size()) {
    return Validation::fail("entry count " + std::to_string(ms.size()) +
                            " does not match input size " +
                            std::to_string(input.size()));
  }

  std::size_t previous_len = 0;

  for (std::size_t i = 0; i < ms.size(); ++i) {
    const std::size_t pos = ms[i].first;
    const std::size_t len = ms[i].second;

    const std::string where = " at input position " + std::to_string(i);

    if (len > input.size() - i) {
      return Validation::fail("match runs past the end of the input" + where +
                              " (length " + std::to_string(len) + ")");
    }

    if (len > 0 && pos + len > reference.size()) {
      return Validation::fail("match runs past the end of the reference" +
                              where + " (position " + std::to_string(pos) +
                              ", length " + std::to_string(len) + ")");
    }

    // A symbol of the reference always matches at least itself; a symbol
    // that is absent from the reference cannot match at all.
    const bool present = alphabet.contains(input[i]);

    if (present && len == 0) {
      return Validation::fail("length 0 although the symbol occurs in the "
                              "reference" +
                              where);
    }

    if (!present && len > 0) {
      return Validation::fail("non-zero length although the symbol is absent "
                              "from the reference" +
                              where);
    }

    // Dropping the first character of a match leaves a match.
    if (i > 0 && len + 1 < previous_len) {
      return Validation::fail("length drops by more than one" + where +
                              " (previous " + std::to_string(previous_len) +
                              ", current " + std::to_string(len) + ")");
    }

    previous_len = len;
  }

  return {};
}

// True if input[start, start + length) occurs in the reference.
template <typename Symbol, typename SAType>
inline bool pattern_occurs(const std::vector<Symbol>& reference,
                           const std::vector<SAType>& suffix_array,
                           const std::vector<Symbol>& input, std::size_t start,
                           std::size_t length) {
  std::size_t low = 0;
  std::size_t high = suffix_array.size();

  while (low < high) {
    const std::size_t mid = low + (high - low) / 2;
    const std::size_t suffix = static_cast<std::size_t>(suffix_array[mid]);

    // Compare the suffix with the pattern; a suffix that ends first sorts
    // before the pattern.
    int order = 0;

    for (std::size_t k = 0; k < length && order == 0; ++k) {
      if (suffix + k >= reference.size()) {
        order = -1;
      } else if (reference[suffix + k] != input[start + k]) {
        order = reference[suffix + k] < input[start + k] ? -1 : 1;
      }
    }

    if (order == 0) {
      return true;
    }

    if (order < 0) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  return false;
}

// Checks one position. `implied` says the occurrence is the previous,
// already verified, one shifted by a character, so it needs no comparison.
template <typename Symbol, typename SAType>
inline Validation check_position(const MatchingStatistics& ms,
                                 const std::vector<Symbol>& input,
                                 const std::vector<Symbol>& reference,
                                 const std::vector<SAType>& suffix_array,
                                 std::size_t i, bool maximality,
                                 bool implied) {
  const std::size_t pos = ms[i].first;
  const std::size_t len = ms[i].second;

  if (len > input.size() - i) {
    return Validation::fail("match runs past the end of the input at input "
                            "position " +
                            std::to_string(i));
  }

  if (len > 0) {
    if (pos + len > reference.size()) {
      return Validation::fail("occurrence outside the reference at input "
                              "position " +
                              std::to_string(i) + " (position " +
                              std::to_string(pos) + ", length " +
                              std::to_string(len) + ")");
    }

    if (!implied) {
      for (std::size_t j = 0; j < len; ++j) {
        if (reference[pos + j] != input[i + j]) {
          return Validation::fail(
              "reported occurrence does not match the input at input "
              "position " +
              std::to_string(i) + " (position " + std::to_string(pos) +
              ", length " + std::to_string(len) + ", mismatch at offset " +
              std::to_string(j) + ")");
        }
      }
    }
  }

  if (maximality && i + len < input.size() &&
      pattern_occurs(reference, suffix_array, input, i, len + 1)) {
    return Validation::fail("match is not maximal at input position " +
                            std::to_string(i) + " (length " +
                            std::to_string(len) + " can be extended)");
  }

  return {};
}

// Every position.
template <typename Symbol, typename SAType>
inline Validation verify_all_positions(const MatchingStatistics& ms,
                                       const std::vector<Symbol>& input,
                                       const std::vector<Symbol>& reference,
                                       const std::vector<SAType>& suffix_array,
                                       bool maximality) {
  if (ms.size() != input.size()) {
    return Validation::fail("entry count does not match input size");
  }

  for (std::size_t i = 0; i < ms.size(); ++i) {
    const bool implied = i > 0 && ms[i - 1].second > 0 &&
                         ms[i].first == ms[i - 1].first + 1U &&
                         ms[i].second + 1U == ms[i - 1].second;

    const Validation result = check_position(ms, input, reference,
                                             suffix_array, i, maximality,
                                             implied);

    if (!result.ok) {
      return result;
    }
  }

  return {};
}

// A random sample of positions, plus the first and the last.
template <typename Symbol, typename SAType>
inline Validation verify_sampled(const MatchingStatistics& ms,
                                 const std::vector<Symbol>& input,
                                 const std::vector<Symbol>& reference,
                                 const std::vector<SAType>& suffix_array,
                                 std::size_t sample, std::uint64_t seed,
                                 bool maximality) {
  if (ms.size() != input.size()) {
    return Validation::fail("entry count does not match input size");
  }

  if (ms.empty()) {
    return {};
  }

  if (sample == 0 || sample >= ms.size()) {
    return verify_all_positions(ms, input, reference, suffix_array,
                                maximality);
  }

  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<std::size_t> pick(0, ms.size() - 1);

  std::vector<std::size_t> positions;
  positions.reserve(sample + 2);
  positions.push_back(0);
  positions.push_back(ms.size() - 1);

  for (std::size_t k = 0; k < sample; ++k) {
    positions.push_back(pick(rng));
  }

  std::sort(positions.begin(), positions.end());
  positions.erase(std::unique(positions.begin(), positions.end()),
                  positions.end());

  for (const std::size_t i : positions) {
    const Validation result = check_position(ms, input, reference,
                                             suffix_array, i, maximality,
                                             false);

    if (!result.ok) {
      return result;
    }
  }

  return {};
}

// Independent O(|input| * |reference|) computation of the lengths. Only for
// small inputs.
template <typename Symbol>
inline Validation verify_against_brute_force(
    const MatchingStatistics& ms, const std::vector<Symbol>& input,
    const std::vector<Symbol>& reference) {
  if (ms.size() != input.size()) {
    return Validation::fail("entry count does not match input size");
  }

  for (std::size_t i = 0; i < input.size(); ++i) {
    std::size_t best = 0;

    for (std::size_t start = 0; start < reference.size(); ++start) {
      std::size_t length = 0;

      while (i + length < input.size() && start + length < reference.size() &&
             input[i + length] == reference[start + length]) {
        ++length;
      }

      best = std::max(best, length);
    }

    if (ms[i].second != best) {
      return Validation::fail(
          "brute-force length differs at input position " +
          std::to_string(i) + " (brute force " + std::to_string(best) +
          ", reported " + std::to_string(ms[i].second) + ")");
    }
  }

  return {};
}

// ---------- Per-file results ----------

struct FileRunResult {
  std::string filename;
  std::size_t input_bytes = 0;
  std::string implementation;
  bool is_baseline = false;

  double build_ms = 0.0;
  Timing timing;
  double baseline_min_ms = 0.0;

  Digest digest;
  Validation invariants;
  Validation positions;
  LengthComparison lengths;

  double speedup() const {
    return timing.min_ms == 0.0 ? 0.0 : baseline_min_ms / timing.min_ms;
  }
};

// ---------- Per-implementation totals ----------

struct ImplementationTotals {
  std::string name;
  bool is_baseline = false;
  double build_ms = 0.0;

  double total_min_ms = 0.0;
  double total_first_ms = 0.0;
  std::size_t total_entries = 0;

  std::size_t files_compared = 0;
  std::size_t files_lengths_equal = 0;
  std::size_t invariant_failures = 0;
  std::size_t position_failures = 0;

  // Length digest summed over every file (max_len is the maximum).
  std::uint64_t total_len = 0;
  std::uint32_t max_len = 0;
  std::size_t zero_len_count = 0;

  // The implementation's per-call diagnostics (phase times, lookup
  // counters), summed over every file.
  Diagnostics diagnostics;

  void accumulate(const FileRunResult& result) {
    total_min_ms += result.timing.min_ms;
    total_first_ms += result.timing.first_ms;
    total_entries += result.digest.entries;
    total_len += result.digest.total_len;
    max_len = std::max(max_len, result.digest.max_len);
    zero_len_count += result.digest.zero_len_count;

    if (result.lengths.checked) {
      ++files_compared;

      if (result.lengths.equal) {
        ++files_lengths_equal;
      }
    }

    if (!result.invariants.ok) {
      ++invariant_failures;
    }

    if (!result.positions.ok) {
      ++position_failures;
    }
  }

  bool all_lengths_equal() const { return files_lengths_equal == files_compared; }

  bool ok() const {
    return invariant_failures == 0 && position_failures == 0 &&
           all_lengths_equal();
  }
};

// ---------- Output files ----------

class CSVWriter {
 public:
  explicit CSVWriter(const std::string& path) : output_(path) {
    if (!output_) {
      throw std::runtime_error("cannot create results file: " + path);
    }

    output_ << "file,input_bytes,implementation,is_baseline,build_ms,first_ms,"
               "min_ms,baseline_min_ms,speedup,MB_per_s,entries,total_len,"
               "max_len,zero_len,lengths_equal,divergent_positions,"
               "invariants_ok,positions_ok,peak_RSS_MB\n";
  }

  void write_row(const FileRunResult& result) {
    const double seconds = result.timing.min_ms / 1000.0;
    const double megabytes =
        static_cast<double>(result.input_bytes) / (1024.0 * 1024.0);

    output_ << result.filename << ',' << result.input_bytes << ','
            << result.implementation << ',' << (result.is_baseline ? 1 : 0)
            << ',' << std::fixed << std::setprecision(3) << result.build_ms
            << ',' << result.timing.first_ms << ',' << result.timing.min_ms
            << ',' << result.baseline_min_ms << ',' << std::setprecision(4)
            << result.speedup() << ',' << std::setprecision(3)
            << (seconds == 0.0 ? 0.0 : megabytes / seconds) << ','
            << result.digest.entries << ',' << result.digest.total_len << ','
            << result.digest.max_len << ',' << result.digest.zero_len_count
            << ','
            << (result.is_baseline ? "NA"
                                   : (result.lengths.equal ? "YES" : "NO"))
            << ',' << result.lengths.divergent_positions << ','
            << (result.invariants.ok ? "YES" : "NO") << ','
            << (result.positions.ok ? "YES" : "NO") << ','
            << std::setprecision(2) << peak_rss_mb() << '\n';

    output_.flush();
  }

  void write_summary(const std::vector<ImplementationTotals>& totals,
                     std::size_t processed_files,
                     std::size_t total_input_bytes,
                     std::size_t reference_bytes,
                     std::size_t distinct_symbols) {
    const double baseline_total =
        totals.empty() ? 0.0 : totals.front().total_min_ms;

    output_ << "summary,all,processed_files," << processed_files << '\n';
    output_ << "summary,all,total_input_bytes," << total_input_bytes << '\n';
    output_ << "summary,all,reference_bytes," << reference_bytes << '\n';
    output_ << "summary,all,distinct_symbols," << distinct_symbols << '\n';

    for (const ImplementationTotals& implementation : totals) {
      const std::string prefix = "summary," + implementation.name + ",";

      output_ << prefix << "build_ms," << std::fixed << std::setprecision(3)
              << implementation.build_ms << '\n';
      output_ << prefix << "total_min_ms," << implementation.total_min_ms
              << '\n';
      output_ << prefix << "total_first_ms," << implementation.total_first_ms
              << '\n';
      output_ << prefix << "speedup_vs_baseline," << std::setprecision(4)
              << (implementation.total_min_ms == 0.0
                      ? 0.0
                      : baseline_total / implementation.total_min_ms)
              << '\n';
      output_ << prefix << "total_entries," << implementation.total_entries
              << '\n';
      output_ << prefix << "files_lengths_equal,"
              << implementation.files_lengths_equal << '/'
              << implementation.files_compared << '\n';
      output_ << prefix << "invariant_failures,"
              << implementation.invariant_failures << '\n';
      output_ << prefix << "position_failures,"
              << implementation.position_failures << '\n';
      output_ << prefix << "ok," << (implementation.ok() ? "YES" : "NO")
              << '\n';
    }

    output_ << "summary,all,peak_RSS_MB," << std::setprecision(2)
            << peak_rss_mb() << '\n';
    output_.flush();
  }

 private:
  std::ofstream output_;
};

// One line per (file, implementation): what an implementation produced, in a
// form that can be diffed between runs and machines.
class ChecksumWriter {
 public:
  explicit ChecksumWriter(const std::string& path) : output_(path) {
    if (!output_) {
      throw std::runtime_error("cannot create checksum file: " + path);
    }

    output_ << "file,implementation,entries,total_len,max_len,zero_len,"
               "len_hash,pos_hash\n";
  }

  void write(const FileRunResult& result) {
    output_ << result.filename << ',' << result.implementation << ','
            << result.digest.entries << ',' << result.digest.total_len << ','
            << result.digest.max_len << ',' << result.digest.zero_len_count
            << ',' << to_hex(result.digest.len_hash) << ','
            << to_hex(result.digest.pos_hash) << '\n';

    output_.flush();
  }

 private:
  std::ofstream output_;
};

// Writes one implementation's output for one file: a 64-bit entry count, then
// (position, length) as pairs of 32-bit values.
inline void dump_ms(const std::string& directory, std::size_t file_index,
                    const std::string& filename,
                    const std::string& implementation,
                    const MatchingStatistics& ms) {
  fs::create_directories(directory);

  const std::string path = directory + "/" + std::to_string(file_index) + "_" +
                           fs::path(filename).filename().string() + "." +
                           implementation + ".ms.bin";

  std::ofstream output(path, std::ios::binary);

  if (!output) {
    throw std::runtime_error("cannot create dump file: " + path);
  }

  const std::uint64_t count = ms.size();
  output.write(reinterpret_cast<const char*>(&count), sizeof(count));

  for (const auto& entry : ms) {
    const std::uint32_t values[2] = {entry.first, entry.second};
    output.write(reinterpret_cast<const char*>(values), sizeof(values));
  }

  if (!output) {
    throw std::runtime_error("failed while writing dump file: " + path);
  }
}

}  // namespace msbench
