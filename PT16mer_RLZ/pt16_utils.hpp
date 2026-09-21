#pragma once

#include <sys/resource.h>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "parser.hpp"

namespace fs = std::filesystem;

using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;

using MatchingStatistics =
    std::vector<std::pair<std::uint32_t, std::uint32_t>>;

// ---------- 16-mer encoding ----------

constexpr std::uint32_t KMER_LENGTH = 16;

// Two-bit code of each DNA character; every other character maps to 0.
inline std::array<std::uint8_t, 256> build_alphatab() {
  std::array<std::uint8_t, 256> alphatab{};

  alphatab[static_cast<unsigned char>('A')] = 0;  // 00
  alphatab[static_cast<unsigned char>('C')] = 1;  // 01
  alphatab[static_cast<unsigned char>('G')] = 2;  // 10
  alphatab[static_cast<unsigned char>('T')] = 3;  // 11

  return alphatab;
}

inline const std::array<std::uint8_t, 256> alphatab = build_alphatab();

// Packs one 16-mer from the reference into a 32-bit key.

inline std::uint32_t encode_16mer(const std::vector<unsigned char>& reference,
                                  const std::uint32_t position) {
  std::uint32_t key = 0;

  for (std::uint32_t j = 0; j < KMER_LENGTH; ++j) {
    const std::uint8_t code =
        alphatab[static_cast<unsigned char>(reference[position + j])];
    key = (key << 2U) | code;
  }

  return key;
}

// Number of common leading bits between two 16-bit bucket prefixes.

inline std::uint32_t lcp_bits_16(const std::uint32_t first,
                                 const std::uint32_t second) {
  const std::uint32_t difference = (first ^ second) << 16U;
  return std::countl_zero(difference);
}

// ---------- PT16 buckets ----------

// A 16-mer key is split into a 16-bit bucket (its high bits) and a 16-bit
// low part.
constexpr std::uint32_t BUCKET_SIZE = 65536;
constexpr std::uint32_t LOW_BITS = 16;
constexpr std::uint32_t LOW_MASK = BUCKET_SIZE - 1;
constexpr std::uint32_t NUMBER_OF_BUCKETS = 65536;
constexpr std::uint32_t EMPTY_BUCKET_FLAG = 1U << 31;
constexpr std::uint16_t LARGE_OFFSET_FLAG = 65535;

/*
Builds the normal H prefix sums, then marks each empty bucket.
`count` is the number of PT16 entries in each bucket.
For an empty bucket b:
    MSB = 1
    lower 31 bits = non-empty bucket x with maximum LCP with b.
*/

inline std::vector<std::uint32_t> build_H(
    const std::array<std::uint32_t, NUMBER_OF_BUCKETS>& count) {
  std::vector<std::uint32_t> H(static_cast<std::size_t>(NUMBER_OF_BUCKETS) + 1,
                               0);
  std::vector<std::uint32_t> non_empty_buckets;

  H[0] = 0;

  // First construct the ordinary H directory.
  for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
    H[bucket + 1] = H[bucket] + count[bucket];

    if (count[bucket] != 0) {
      non_empty_buckets.push_back(bucket);
    }
  }

  // Encode every empty bucket.
  for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
    if (count[bucket] != 0) {
      continue;
    }

    // The best LCP candidate must be the nearest non-empty
    // bucket on the left or right in prefix order.
    const auto next = std::lower_bound(non_empty_buckets.begin(),
                                       non_empty_buckets.end(), bucket);

    std::uint32_t x;

    if (next == non_empty_buckets.begin()) {
      x = *next;
    } else if (next == non_empty_buckets.end()) {
      x = non_empty_buckets.back();
    } else {
      const std::uint32_t left = *(next - 1);
      const std::uint32_t right = *next;

      if (lcp_bits_16(bucket, left) >= lcp_bits_16(bucket, right)) {
        x = left;
      } else {
        x = right;
      }
    }

    H[bucket] = EMPTY_BUCKET_FLAG | x;
  }

  return H;
}

// ---------- Binary output ----------

// Writes the raw bytes of one value.
template <typename T>
inline void write_value(std::ofstream& output, const T& value) {
  output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

// Writes the raw bytes of a vector's elements, with no length prefix.
template <typename T>
inline void write_vector(std::ofstream& output, const std::vector<T>& values) {
  if (!values.empty()) {
    output.write(reinterpret_cast<const char*>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(T)));
  }
}

// ---------- Arguments ----------

struct Args {
  std::string reference;
  std::string suffix_array;
  std::string filenames;
  std::string pt16_table;
  std::string results;
};

inline std::string require_value(int& i, int argc, char** argv) {
  if (i + 1 >= argc) {
    throw std::runtime_error("missing command-line value");
  }

  return argv[++i];
}

inline void print_usage(const char* program) {
  std::cout << "Usage: " << program << " --reference PATH"
            << " --suffix-array PATH"
            << " --filenames PATH"
            << " --results PATH"
            << " [--pt16-table PATH]\n";
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
    } else if (option == "--pt16-table") {
      args.pt16_table = require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = require_value(i, argc, argv);
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

  return args;
}

// ---------- Loading ----------

template <typename Symbol>
inline std::vector<Symbol> load_reference(const std::string& path) {
  return read_file<Symbol>(path.c_str());
}

template <typename SAType>
inline std::vector<SAType> load_suffix_array(const std::string& path) {
  return read_file<SAType>(path.c_str());
}

template <typename Symbol>
inline std::vector<Symbol> load_input(const std::string& path) {
  return read_file<Symbol>(path.c_str());
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

// ---------- PT16 table path ----------

inline std::string get_pt16_path(const Args& args) {
  if (!args.pt16_table.empty()) {
    return args.pt16_table;
  }

  const fs::path reference_path(args.reference);
  const std::string dataset = reference_path.parent_path().filename().string();

  return (reference_path.parent_path() / (dataset + "_pt16_hl.bin")).string();
}

// ---------- Timing ----------

template <typename Fn>
inline double time_ms(Fn&& fn) {
  const auto start = std::chrono::steady_clock::now();
  fn();
  const auto end = std::chrono::steady_clock::now();

  return std::chrono::duration<double, std::milli>(end - start).count();
}

inline double peak_rss_mb() {
  struct rusage usage{};

  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return 0.0;
  }

  return static_cast<double>(usage.ru_maxrss) / 1024.0;
}

// ---------- Baseline results ----------

struct BaselineResult {
  std::string filename;
  std::size_t input_bytes;
  double baseline_ms;
  std::size_t baseline_phrases;
  std::string factor_file;
};

// ---------- Baseline factor storage ----------

inline void write_factor_file(const std::string& path, const Triples& factors) {
  std::ofstream output(path, std::ios::binary);

  if (!output) {
    throw std::runtime_error("cannot create baseline factor file");
  }

  const std::uint64_t count = factors.size();
  output.write(reinterpret_cast<const char*>(&count), sizeof(count));

  for (const auto& factor : factors) {
    const std::uint64_t input_position = std::get<0>(factor);
    const std::uint64_t reference_position = std::get<1>(factor);
    const std::uint64_t match_length = std::get<2>(factor);

    output.write(reinterpret_cast<const char*>(&input_position),
                 sizeof(input_position));
    output.write(reinterpret_cast<const char*>(&reference_position),
                 sizeof(reference_position));
    output.write(reinterpret_cast<const char*>(&match_length),
                 sizeof(match_length));
  }
}

template <typename Symbol>
inline bool factor_file_equals(const std::string& path, const Triples& factors,
                               const std::vector<Symbol>& input_sequence,
                               const std::vector<Symbol>& reference) {
  std::ifstream input(path, std::ios::binary);

  if (!input) {
    throw std::runtime_error("cannot open baseline factor file");
  }

  std::uint64_t count = 0;
  input.read(reinterpret_cast<char*>(&count), sizeof(count));

  // Check 1: same factor count.
  if (!input) {
    std::cerr << "CORRECTNESS FAILURE: could not read baseline factor count\n";
    return false;
  }

  if (count != factors.size()) {
    std::cerr << "CORRECTNESS FAILURE\n"
              << "check=factor_count\n"
              << "baseline=" << count << '\n'
              << "pt16=" << factors.size() << '\n';
    return false;
  }

  std::size_t expected_input_position = 0;

  for (std::size_t i = 0; i < factors.size(); ++i) {
    std::uint64_t baseline_input_position = 0;
    std::uint64_t baseline_reference_position = 0;
    std::uint64_t baseline_match_length = 0;

    input.read(reinterpret_cast<char*>(&baseline_input_position),
               sizeof(baseline_input_position));
    input.read(reinterpret_cast<char*>(&baseline_reference_position),
               sizeof(baseline_reference_position));
    input.read(reinterpret_cast<char*>(&baseline_match_length),
               sizeof(baseline_match_length));

    if (!input) {
      std::cerr << "CORRECTNESS FAILURE\n"
                << "check=baseline_factor_read\n"
                << "factor_index=" << i << '\n';
      return false;
    }

    const std::size_t pt16_input_position = std::get<0>(factors[i]);
    const std::size_t pt16_reference_position = std::get<1>(factors[i]);
    const std::size_t pt16_match_length = std::get<2>(factors[i]);

    // Check 2: same factor start.
    if (pt16_input_position != baseline_input_position) {
      std::cerr << "CORRECTNESS FAILURE\n"
                << "check=input_position\n"
                << "factor_index=" << i << '\n'
                << "baseline=" << baseline_input_position << '\n'
                << "pt16=" << pt16_input_position << '\n';
      return false;
    }

    // Check 3: same greedy factor length.
    if (pt16_match_length != baseline_match_length) {
      std::cerr << "CORRECTNESS FAILURE\n"
                << "check=match_length\n"
                << "factor_index=" << i << '\n'
                << "input_position=" << pt16_input_position << '\n'
                << "baseline=" << baseline_match_length << '\n'
                << "pt16=" << pt16_match_length << '\n';
      return false;
    }

    // Check 4: no gaps or overlaps.
    if (pt16_input_position != expected_input_position) {
      std::cerr << "CORRECTNESS FAILURE\n"
                << "check=tiling\n"
                << "factor_index=" << i << '\n'
                << "expected_input_position=" << expected_input_position << '\n'
                << "actual_input_position=" << pt16_input_position << '\n';
      return false;
    }

    if (pt16_match_length == 0) {
      std::cerr << "CORRECTNESS FAILURE\n"
                << "check=zero_length_factor\n"
                << "factor_index=" << i << '\n';
      return false;
    }

    expected_input_position += pt16_match_length;

    // Literal factor.
    if (pt16_match_length == 1) {
      if (pt16_reference_position !=
          static_cast<std::size_t>(input_sequence[pt16_input_position])) {
        std::cerr << "CORRECTNESS FAILURE\n"
                  << "check=literal_value\n"
                  << "factor_index=" << i << '\n'
                  << "input_position=" << pt16_input_position << '\n';
        return false;
      }

      continue;
    }

    // Check 5: copy stays inside reference.
    if (pt16_reference_position + pt16_match_length > reference.size()) {
      std::cerr << "CORRECTNESS FAILURE\n"
                << "check=reference_bounds\n"
                << "factor_index=" << i << '\n'
                << "ref_pos=" << pt16_reference_position << '\n'
                << "length=" << pt16_match_length << '\n'
                << "reference_size=" << reference.size() << '\n';
      return false;
    }

    // Check 6: PT16 reference occurrence really matches the input.
    for (std::size_t j = 0; j < pt16_match_length; ++j) {
      if (reference[pt16_reference_position + j] !=
          input_sequence[pt16_input_position + j]) {
        std::cerr << "CORRECTNESS FAILURE\n"
                  << "check=reference_match\n"
                  << "factor_index=" << i << '\n'
                  << "input_position=" << pt16_input_position << '\n'
                  << "ref_pos=" << pt16_reference_position << '\n'
                  << "length=" << pt16_match_length << '\n'
                  << "mismatch_offset=" << j << '\n';
        return false;
      }
    }
  }

  // Check 7: factors cover the complete input.
  if (expected_input_position != input_sequence.size()) {
    std::cerr << "CORRECTNESS FAILURE\n"
              << "check=final_coverage\n"
              << "covered=" << expected_input_position << '\n'
              << "input_size=" << input_sequence.size() << '\n';
    return false;
  }

  return true;
}
// ---------- PT16 per-file statistics ----------

struct PT16Delta {
  std::size_t hits;
  std::size_t misses;
  std::size_t singleton_hits;
  std::size_t range_hits;
};

template <typename Stats>
inline PT16Delta stats_difference(const Stats& previous, const Stats& current) {
  return {current.hits - previous.hits, current.misses - previous.misses,
          current.singleton_hits - previous.singleton_hits,
          current.range_hits - previous.range_hits};
}

// ---------- Final summary ----------

struct BenchmarkSummary {
  std::size_t processed_files;
  std::size_t total_input_bytes;
  std::size_t total_baseline_phrases;
  std::size_t total_pt16_phrases;
  double total_baseline_ms;
  double total_pt16_ms;
  double pt16_build_ms;
  bool all_equal;
};

// ---------- CSV ----------

class CSVWriter {
 public:
  explicit CSVWriter(const std::string& path) : output(path) {
    if (!output) {
      throw std::runtime_error("cannot create results file: " + path);
    }

    // output <<
    // "file,input_bytes,baseline_ms,pt16_ms,speedup,baseline_phrases,pt16_phrases,outputs_equal,pt16_hits,pt16_misses,pt16_hit_rate,pt16_entries,pt16_MB,peak_RSS_MB\n";

    output << "file,input_bytes,baseline_ms,pt16_ms,speedup,baseline_phrases,"
              "pt16_phrases,outputs_equal,pt16_hits,pt16_misses,pt16_hit_rate,"
              "singleton_hits,range_hits,singleton_hit_rate,pt16_entries,pt16_"
              "MB,peak_RSS_MB\n";
  }

  void write_row(const BaselineResult& baseline, double pt16_ms,
                 std::size_t pt16_phrases, bool equal, const PT16Delta& stats,
                 std::size_t entries, std::size_t approx_bytes) {
    const std::size_t queries = stats.hits + stats.misses;
    const double hit_rate = queries == 0 ? 0.0
                                         : static_cast<double>(stats.hits) /
                                               static_cast<double>(queries);
    const double speedup =
        pt16_ms == 0.0 ? 0.0 : baseline.baseline_ms / pt16_ms;
    const std::size_t hit_types = stats.singleton_hits + stats.range_hits;
    const double singleton_hit_rate =
        hit_types == 0 ? 0.0
                       : static_cast<double>(stats.singleton_hits) /
                             static_cast<double>(hit_types);

    output << baseline.filename << ',' << baseline.input_bytes << ','
           << std::fixed << std::setprecision(3) << baseline.baseline_ms << ','
           << pt16_ms << ',' << speedup << ',' << baseline.baseline_phrases
           << ',' << pt16_phrases << ',' << (equal ? "YES" : "NO") << ','
           << stats.hits << ',' << stats.misses << ',' << std::setprecision(6)
           << hit_rate << ',' << stats.singleton_hits << ',' << stats.range_hits
           << ',' << singleton_hit_rate << ',' << entries << ','
           << std::setprecision(3)
           << static_cast<double>(approx_bytes) / (1024.0 * 1024.0) << ','
           << std::setprecision(2) << peak_rss_mb() << '\n';
  }

  template <typename Stats>
  void write_summary(const BenchmarkSummary& summary, const Stats& stats) {
    const double speedup =
        summary.total_pt16_ms == 0.0
            ? 0.0
            : summary.total_baseline_ms / summary.total_pt16_ms;
    const std::size_t hit_types = stats.singleton_hits + stats.range_hits;
    const double singleton_hit_rate =
        hit_types == 0 ? 0.0
                       : static_cast<double>(stats.singleton_hits) /
                             static_cast<double>(hit_types);

    output << "summary,processed_files," << summary.processed_files << '\n';
    output << "summary,total_input_bytes," << summary.total_input_bytes << '\n';
    output << "summary,total_baseline_phrases,"
           << summary.total_baseline_phrases << '\n';
    output << "summary,total_pt16_phrases," << summary.total_pt16_phrases
           << '\n';
    output << "summary,total_baseline_ms," << std::fixed << std::setprecision(3)
           << summary.total_baseline_ms << '\n';
    output << "summary,total_pt16_ms," << summary.total_pt16_ms << '\n';
    output << "summary,overall_speedup," << speedup << '\n';
    output << "summary,all_outputs_equal," << (summary.all_equal ? "YES" : "NO")
           << '\n';
    output << "summary,total_pt16_hits," << stats.hits << '\n';
    output << "summary,total_pt16_misses," << stats.misses << '\n';
    output << "summary,pt16_entries," << stats.entries << '\n';
    output << "summary,pt16_MB," << std::setprecision(3)
           << static_cast<double>(stats.approx_bytes) / (1024.0 * 1024.0)
           << '\n';
    output << "summary,peak_RSS_MB," << std::setprecision(2) << peak_rss_mb()
           << '\n';
    output << "summary,total_singleton_hits," << stats.singleton_hits << '\n';
    output << "summary,total_range_hits," << stats.range_hits << '\n';
    output << "summary,singleton_hit_rate," << std::setprecision(6)
           << singleton_hit_rate << '\n';
    output << "summary,pt16_build_ms," << summary.pt16_build_ms << '\n';
  }

 private:
  std::ofstream output;
};