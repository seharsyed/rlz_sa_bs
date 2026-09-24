#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "../pt16_build_v2.hpp"
#include "../variants/pt16_rlz_v2_fastmiss.hpp"
#include "../pt16_utils.hpp"
#include "ms_tools.hpp"
#include "sorted_kmer_scan.hpp"

/**
 * Matching-statistics variants over PT16FastMissParser
 * (pt16_rlz_v2_fastmiss.hpp): the v2 table with a cheaper miss path. They
 * mirror the v2 sorted and bucketed rows exactly -- same scan order, same
 * backwardChainExtend -- so any difference between the rows is the lookup
 * itself.
 */

// Builds a fresh v2 table at `table_path` and loads it as a fast-miss
// parser. Each variant gets its own table file so it never depends on
// another variant's.
template <typename T1, typename T2>
std::unique_ptr<PT16FastMissParser<T1, T2>> make_fastmiss_parser(
    const std::vector<T1>& reference, const std::vector<T2>& suffix_array,
    const std::string& table_path) {
  fs::remove(table_path);
  build_pt16_table(reference, suffix_array, table_path);

  return std::make_unique<PT16FastMissParser<T1, T2>>(reference, suffix_array,
                                                      table_path);
}

/**
 * Scan only: every hit capped at 16, like pt16-v2-sorted-scan.
 *
 * Builds its own copy of the v2 table (at `table_path`, typically the v2
 * path plus a suffix) so it never depends on another variant's file.
 */
template <typename T1, typename T2>
class PT16FastMissSortedScanMS {
 public:
  PT16FastMissSortedScanMS(const std::vector<T1>& reference,
                           const std::vector<T2>& suffix_array,
                           const std::string& table_path)
      : parser_(make_fastmiss_parser(reference, suffix_array, table_path)) {}

  PT16FastMissSortedScanMS(const PT16FastMissSortedScanMS&) = delete;
  PT16FastMissSortedScanMS& operator=(const PT16FastMissSortedScanMS&) =
      delete;

  std::vector<KmerLookupResult> sortedScan(const std::vector<T1>& input) {
    const auto& parser = *parser_;
    const auto before = parser.stats();

    std::vector<KmerLookupResult> results = sortedKmerScan(
        input,
        [&](const std::vector<T1>&, std::size_t, std::uint32_t key) {
          return parser.lookupKmerByKey(key);
        },
        [&](const std::vector<T1>& in, std::size_t i, std::uint32_t key) {
          return parser.lookupTailByKey(
              key, static_cast<std::uint32_t>(in.size() - i));
        },
        &timings_);

    const auto after = parser.stats();

    diagnostics_ = sortedScanDiagnostics(timings_) +
                   miss_diagnostics(before, after);

    search_composition_ = {true, after.singleton_hits - before.singleton_hits,
                           after.range_hits - before.range_hits,
                           after.misses - before.misses};

    return results;
  }

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    const std::vector<KmerLookupResult> results = sortedScan(input);

    MatchingStatistics ms(results.size());
    for (std::size_t i = 0; i < results.size(); ++i) {
      ms[i] = {results[i].match_position, results[i].match_length};
    }

    return ms;
  }

  Diagnostics diagnostics() const { return diagnostics_; }

  SearchComposition searchComposition() const { return search_composition_; }

  // Scan-only: every hit is capped at exactly 16.
  bool exactExpected() const { return false; }

  const PT16FastMissParser<T1, T2>& parser() const { return *parser_; }

 private:
  std::unique_ptr<PT16FastMissParser<T1, T2>> parser_;
  SortedScanTimings timings_;
  SearchComposition search_composition_;
  Diagnostics diagnostics_;
};

/**
 * Same sorted scan, fed through multi-step chain extension -- the
 * counterpart of pt16-v2-sorted-chain-multi, expected to be exact.
 */
template <typename T1, typename T2>
class PT16FastMissSortedChainMS {
 public:
  PT16FastMissSortedChainMS(const std::vector<T1>& reference,
                            const std::vector<T2>& suffix_array,
                            const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16FastMissSortedChainMS(const PT16FastMissSortedChainMS&) = delete;
  PT16FastMissSortedChainMS& operator=(const PT16FastMissSortedChainMS&) =
      delete;

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    return timedBackwardChainExtend(impl_.sortedScan(input), chain_ms_);
  }

  // The scan's own phases, plus the chain extension run on its result.
  Diagnostics diagnostics() const {
    Diagnostics diagnostics = impl_.diagnostics();
    diagnostics.add_phase("chain", chain_ms_);
    return diagnostics;
  }

  SearchComposition searchComposition() const {
    return impl_.searchComposition();
  }

 private:
  PT16FastMissSortedScanMS<T1, T2> impl_;
  double chain_ms_ = 0.0;  // last backwardChainExtend call
};

/**
 * The fast-miss lookup with prebucketing instead of a full sort: the same
 * phases as pt16-v2-bucket-scan (PT16ScanMS::bucketedScan) -- roll keys
 * and count per table bucket, one counting-sort pass by bucket, probe in
 * bucket order, tail -- so its phase line compares directly with that row
 * and with pt16-v2-fastmiss-sorted-scan.
 */
template <typename T1, typename T2>
class PT16FastMissBucketScanMS {
 public:
  PT16FastMissBucketScanMS(const std::vector<T1>& reference,
                           const std::vector<T2>& suffix_array,
                           const std::string& table_path)
      : parser_(make_fastmiss_parser(reference, suffix_array, table_path)) {}

  PT16FastMissBucketScanMS(const PT16FastMissBucketScanMS&) = delete;
  PT16FastMissBucketScanMS& operator=(const PT16FastMissBucketScanMS&) =
      delete;

  std::vector<KmerLookupResult> bucketedScan(const std::vector<T1>& input) {
    using clock = std::chrono::steady_clock;
    const auto ms_since = [](const clock::time_point start) {
      return std::chrono::duration<double, std::milli>(clock::now() - start)
          .count();
    };

    const auto& parser = *parser_;
    const auto before = parser.stats();

    const std::size_t n = input.size();
    std::vector<KmerLookupResult> results(n);

    const std::size_t kmer_positions =
        n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;

    // ---------- Prebucket: roll every key, count per table bucket ----------

    auto phase_start = clock::now();

    std::vector<std::uint32_t> keys(kmer_positions);
    auto count =
        std::make_unique<std::array<std::uint32_t, NUMBER_OF_BUCKETS>>();
    count->fill(0);

    std::uint32_t key = kmer_positions > 0 ? encode_16mer(input, 0) : 0;

    for (std::size_t i = 0; i < kmer_positions; ++i) {
      if (i > 0) {
        key = (key << 2U) |
              alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
      }

      keys[i] = key;
      ++(*count)[key >> LOW_BITS];
    }

    phase_barrier(keys.data());
    phase_barrier(count->data());
    const double prebucket_ms = ms_since(phase_start);

    // ---------- Bucket: prefix sum + counting-sort placement ----------

    phase_start = clock::now();

    std::vector<std::uint32_t> offsets(
        static_cast<std::size_t>(NUMBER_OF_BUCKETS) + 1, 0);

    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
      offsets[bucket + 1] = offsets[bucket] + (*count)[bucket];
    }

    std::vector<std::uint32_t> cursor(offsets.begin(),
                                      offsets.begin() + NUMBER_OF_BUCKETS);
    std::vector<std::uint32_t> order(kmer_positions);

    for (std::size_t i = 0; i < kmer_positions; ++i) {
      const std::uint32_t bucket = keys[i] >> LOW_BITS;
      order[cursor[bucket]++] = static_cast<std::uint32_t>(i);
    }

    phase_barrier(order.data());
    const double bucket_ms = ms_since(phase_start);

    // ---------- Probe: the fast-miss lookups, in bucket order ----------

    phase_start = clock::now();

    for (std::size_t slot = 0; slot < kmer_positions; ++slot) {
      const std::uint32_t i = order[slot];
      results[i] = parser.lookupKmerByKey(keys[i]);
    }

    phase_barrier(results.data());
    const double probe_ms = ms_since(phase_start);

    // ---------- Tail: not reordered (at most 15 of them) ----------

    phase_start = clock::now();

    // The last 16-mer key rolled past the end: each `<< 2` drops the first
    // character and pads with a zero one, giving the padded tail key.
    std::uint32_t tail_key =
        kmer_positions > 0 ? keys[kmer_positions - 1] << 2U
                           : (n > 0 ? encode_tail(input, 0) : 0);

    for (std::size_t i = kmer_positions; i < n; ++i, tail_key <<= 2U) {
      const auto tail = parser.lookupTailByKey(
          tail_key, static_cast<std::uint32_t>(n - i));
      results[i].found = false;
      results[i].match_position = tail.match_position;
      results[i].match_length = tail.match_length;
    }

    phase_barrier(results.data());
    const double tail_ms = ms_since(phase_start);

    const auto after = parser.stats();

    diagnostics_ = phase_diagnostics({{"prebucket", prebucket_ms},
                                  {"bucket", bucket_ms},
                                  {"probe", probe_ms},
                                  {"tail", tail_ms}}) +
                   miss_diagnostics(before, after);

    search_composition_ = {true, after.singleton_hits - before.singleton_hits,
                           after.range_hits - before.range_hits,
                           after.misses - before.misses};

    return results;
  }

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    const std::vector<KmerLookupResult> results = bucketedScan(input);

    MatchingStatistics ms(results.size());
    for (std::size_t i = 0; i < results.size(); ++i) {
      ms[i] = {results[i].match_position, results[i].match_length};
    }

    return ms;
  }

  Diagnostics diagnostics() const { return diagnostics_; }

  SearchComposition searchComposition() const { return search_composition_; }

  // Scan-only: every hit is capped at exactly 16.
  bool exactExpected() const { return false; }

 private:
  std::unique_ptr<PT16FastMissParser<T1, T2>> parser_;
  SearchComposition search_composition_;
  Diagnostics diagnostics_;
};

/**
 * Same bucketed scan, fed through multi-step chain extension -- the
 * counterpart of pt16-v2-bucket-chain-multi, expected to be exact.
 */
template <typename T1, typename T2>
class PT16FastMissBucketChainMS {
 public:
  PT16FastMissBucketChainMS(const std::vector<T1>& reference,
                            const std::vector<T2>& suffix_array,
                            const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16FastMissBucketChainMS(const PT16FastMissBucketChainMS&) = delete;
  PT16FastMissBucketChainMS& operator=(const PT16FastMissBucketChainMS&) =
      delete;

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    return timedBackwardChainExtend(impl_.bucketedScan(input), chain_ms_);
  }

  // The scan's own phases, plus the chain extension run on its result.
  Diagnostics diagnostics() const {
    Diagnostics diagnostics = impl_.diagnostics();
    diagnostics.add_phase("chain", chain_ms_);
    return diagnostics;
  }

  SearchComposition searchComposition() const {
    return impl_.searchComposition();
  }

 private:
  PT16FastMissBucketScanMS<T1, T2> impl_;
  double chain_ms_ = 0.0;  // last backwardChainExtend call
};
