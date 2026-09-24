#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../pt16_build_v2.hpp"
#include "../pt16_rlz_v2.hpp"
#include "chain_extend.hpp"
#include "lrf_ms.hpp"
#include "ms_tools.hpp"
#include "ms_utils.hpp"
#include "pt16_fastmiss_ms.hpp"
#include "pt16_sassy_ms.hpp"
#include "sorted_kmer_scan.hpp"

/**
 * The registry of matching-statistics implementations under test.
 *
 * THIS IS THE ONLY FILE YOU EDIT WHEN A NEW VARIANT LANDS.
 *
 * An implementation qualifies if it has
 *
 *   Impl(const std::vector<Symbol>&, const std::vector<SAType>&, ...)
 *   MatchingStatistics computeMatchingStatistics(const std::vector<Symbol>&)
 *
 * Extra constructor arguments (a table path, a block size, a thread
 * count) are forwarded, so a variant that needs them does not need a
 * different adapter.
 *
 * The first registered implementation is the baseline: it is what every
 * other implementation's lengths are compared against, and what the
 * speedup column is relative to.
 */

namespace msbench {

using Symbol = unsigned char;
using SAType = std::uint32_t;

// ---------- Interface ----------

class MSImplementation {
 public:
  virtual ~MSImplementation() = default;

  virtual const std::string& name() const = 0;

  // Wall time of this implementation's construction, in milliseconds.
  virtual double build_ms() const = 0;

  virtual MatchingStatistics compute(const std::vector<Symbol>& input) = 0;

  // Diagnostics from the last compute() call (phase times, lookup
  // counters; empty if there are none). ms_main.cpp sums them over every
  // file and prints one summary per implementation at the end.
  virtual Diagnostics diagnostics() const { return {}; }

  // A structural summary of the built index (e.g. singleton vs. range
  // entry counts), independent of any query -- empty unless the wrapped
  // Index exposes indexComposition() as formatted text. This does not vary
  // between implementations sharing the same table format, so ms_main.cpp
  // prints only the first non-empty one it finds, once, right after
  // [6] BUILD, not once per implementation.
  virtual std::string indexComposition() const { return {}; }

  // The classification of THIS implementation's own last compute() call's
  // raw 16-mer lookups (singleton/range/miss), for the scan-only PT16
  // variants whose own reported lengths already ARE the raw search result.
  // .available is false unless the wrapped Index actually provides this
  // (PT16MS/pt16-brute does not: its own lengths are extended past the raw
  // search, so pairing them with these counts would be misleading).
  virtual SearchComposition searchComposition() const { return {}; }

  // True unless the wrapped Index is known, BY DESIGN, not to produce a
  // real matching-statistics answer (a scan-only floor, capped at 16, or
  // one-step chainExtend, which only ever extends by one character) --
  // see the scan-only/one-step wrapper classes below for the overrides.
  // ms_main.cpp uses this to skip compare_lengths (and the LENGTH
  // MISMATCH report) entirely for those: comparing them against the
  // baseline was never a correctness question -- they are documented
  // lower bounds, not candidate answers -- so doing it on every file,
  // every run, only cost time without ever telling us anything new.
  virtual bool exactExpected() const { return true; }
};

// ---------- Adapter ----------

/**
 * Wraps a concrete index type. The constructor is what gets timed, so
 * whatever preprocessing the variant does lands in build_ms().
 *
 * The index is held by unique_ptr because these classes are typically
 * neither copyable nor movable.
 */
template <typename Index>
class MSAdapter : public MSImplementation {
 public:
  template <typename... CtorArgs>
  MSAdapter(std::string name, CtorArgs&&... args) : name_(std::move(name)) {
    build_ms_ = msbench::time_ms([&] {
      index_ = std::make_unique<Index>(std::forward<CtorArgs>(args)...);
    });
  }

  const std::string& name() const override { return name_; }

  double build_ms() const override { return build_ms_; }

  MatchingStatistics compute(const std::vector<Symbol>& input) override {
    return index_->computeMatchingStatistics(input);
  }

  // Forwards to Index::diagnostics() when the wrapped type has one (the
  // PT16 variants do), otherwise the base class's empty default (LRFMS
  // and any future variant that has nothing to report).
  Diagnostics diagnostics() const override {
    if constexpr (requires(const Index& index) { index.diagnostics(); }) {
      return index_->diagnostics();
    } else {
      return {};
    }
  }

  std::string indexComposition() const override {
    if constexpr (requires(const Index& index) { index.indexComposition(); }) {
      return index_->indexComposition();
    } else {
      return {};
    }
  }

  SearchComposition searchComposition() const override {
    if constexpr (requires(const Index& index) { index.searchComposition(); }) {
      return index_->searchComposition();
    } else {
      return {};
    }
  }

  bool exactExpected() const override {
    if constexpr (requires(const Index& index) { index.exactExpected(); }) {
      return index_->exactExpected();
    } else {
      return true;
    }
  }

 private:
  std::string name_;
  double build_ms_ = 0.0;
  std::unique_ptr<Index> index_;
};

using Implementations = std::vector<std::unique_ptr<MSImplementation>>;

// ---------- Shared PT16 lookup diagnostics ----------

// One line: how many of this call's bucket dispatches searched their
// bucket linearly vs. with std::lower_bound (see PT16RLZParser::Stats /
// PT16SassyLookup::Stats: linear_bucket_searches/binary_bucket_searches).
// Templated on the Stats type so it works for both table formats' Stats
// struct without either needing to know about the other.
template <typename Stats>
inline Diagnostics bucketSearchDiagnostics(const Stats& before,
                                           const Stats& after) {
  return counter_diagnostics(
      "bucket search",
      {{"linear", after.linear_bucket_searches - before.linear_bucket_searches},
       {"binary",
        after.binary_bucket_searches - before.binary_bucket_searches}});
}

template <typename Index, typename... CtorArgs>
inline void add_implementation(Implementations& implementations,
                               std::string name, CtorArgs&&... args) {
  implementations.push_back(std::make_unique<MSAdapter<Index>>(
      std::move(name), std::forward<CtorArgs>(args)...));
}

// ---------- PT16 brute-force matching statistics ----------

/**
 * PT16 matching statistics: one PT16 longest-match query per input position
 * (PT16RLZParser::computeMS_brute).
 *
 * The constructor rebuilds the PT16 table at `table_path` and loads it, so
 * build_ms() covers both the table construction and the load.
 */
template <typename T1, typename T2>
class PT16MS {
 public:
  PT16MS(const std::vector<T1>& reference, const std::vector<T2>& suffix_array,
         const std::string& table_path) {
    // Always rebuild, so a stale table from another reference or another
    // version of the table format is never used.
    fs::remove(table_path);

    build_pt16_table(reference, suffix_array, table_path);

    parser_ = std::make_unique<PT16RLZParser<T1, T2>>(reference, suffix_array,
                                                      table_path);
  }

  PT16MS(const PT16MS&) = delete;
  PT16MS& operator=(const PT16MS&) = delete;

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    const auto before = parser_->stats();

    MatchingStatistics ms = parser_->computeMS_brute(input);

    diagnostics_ = bucketSearchDiagnostics(before, parser_->stats());

    return ms;
  }

  Diagnostics diagnostics() const { return diagnostics_; }

  // A structural property of the built table, not of any query -- printed
  // once for the whole run (see MSImplementation::indexComposition), from
  // whichever implementation is first to offer it. pt16-brute is
  // registered right after the baseline, so in practice this is the one
  // that supplies it.
  std::string indexComposition() const {
    const EntryComposition composition = parser_->entryComposition();
    const std::size_t total =
        composition.singleton_entries + composition.range_entries;

    std::ostringstream out;

    out << "PT16 index: " << total << " distinct 16-mers -- "
        << std::fixed << std::setprecision(1)
        << (total == 0 ? 0.0
                       : 100.0 * static_cast<double>(composition.singleton_entries) /
                             static_cast<double>(total))
        << "% singleton, "
        << (total == 0 ? 0.0
                       : 100.0 * static_cast<double>(composition.range_entries) /
                             static_cast<double>(total))
        << "% range\n";

    return out.str();
  }

 private:
  std::unique_ptr<PT16RLZParser<T1, T2>> parser_;
  Diagnostics diagnostics_;
};

// ---------- PT16 (non-sassy) 16-mer-only scan ----------

/**
 * Same idea as PT16SassyScanMS (pt16_sassy_ms.hpp): the mandatory first
 * scan alone, one 16-mer lookup per input position, no SA interval
 * narrowing/extension -- so the two formats' pure lookup cost is directly
 * comparable, apples to apples, with neither format's extension work
 * counted. Unlike PT16MS::computeMS_brute (which also narrows/extends
 * every match through the suffix array), this only calls
 * PT16RLZParser::lookupKmer/lookupKmerOrTail.
 *
 * Shares its table file with PT16MS (both build the same H/L format); see
 * that class's own comment about what that costs when both are
 * registered.
 */
template <typename T1, typename T2>
class PT16ScanMS {
 public:
  PT16ScanMS(const std::vector<T1>& reference,
            const std::vector<T2>& suffix_array,
            const std::string& table_path) {
    fs::remove(table_path);

    build_pt16_table(reference, suffix_array, table_path);

    parser_ = std::make_unique<PT16RLZParser<T1, T2>>(reference, suffix_array,
                                                      table_path);
  }

  PT16ScanMS(const PT16ScanMS&) = delete;
  PT16ScanMS& operator=(const PT16ScanMS&) = delete;

  // Exposed the same way PT16SassyMS::lookup() is: lets a caller run raw
  // lookups directly (e.g. to build a KmerLookupResult list for
  // chain_extend.hpp) without going through computeMatchingStatistics.
  const PT16RLZParser<T1, T2>& parser() const { return *parser_; }

  // The 16-mer key is a ROLLING window here, not a fresh pack_16mer per
  // position -- same reasoning as PT16SassyMS::scanPass1 (pt16_sassy_ms.hpp):
  // lookupKmer/pack_16mer exist for lzFactorize's skip-ahead calls, where
  // consecutive calls are NOT at consecutive positions; MS visits every
  // position in order, so 15 of the next position's 16 characters are
  // already in the current key. Only i=0 pays a full encode; every later
  // key is `(previous key << 2) | new trailing character`.
  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    const auto before = parser_->stats();
    const std::size_t n = input.size();

    MatchingStatistics ms;
    ms.reserve(n);

    std::uint32_t key = n >= KMER_LENGTH ? encode_16mer(input, 0)
                        : n > 0              ? encode_tail(input, 0)
                                             : 0;

    for (std::size_t i = 0; i < n; ++i) {
      if (n - i >= KMER_LENGTH) {
        if (i > 0) {
          key = (key << 2U) |
                alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
        }

        const auto result = parser_->lookupKmerByKey(input, i, key);
        ms.emplace_back(result.match_position, result.match_length);
      } else {
        // Tail: keep rolling past the end, shifting in zero (padding)
        // characters, so the key is the zero-padded tail.
        if (i > 0) {
          key <<= 2U;
        }

        const auto tail = parser_->lookupTailByKey(
            key, static_cast<std::uint32_t>(n - i));
        ms.emplace_back(tail.match_position, tail.match_length);
      }
    }

    const auto after = parser_->stats();
    diagnostics_ = bucketSearchDiagnostics(before, after);
    search_composition_ = {true, after.singleton_hits - before.singleton_hits,
                           after.range_hits - before.range_hits,
                           after.misses - before.misses};

    return ms;
  }

  /**
   * Same scan, but grouped and run in table-bucket order rather than
   * input order -- the counterpart of
   * PT16SassyMS::computeMatchingStatisticsScanOnlyBucketed, for this table
   * format. See that method's doc comment for the idea; this is the same
   * three-pass mechanics (group by bucket, counting-sort into `order`,
   * then look up bucket by bucket), just calling
   * PT16RLZParser::lookupKmerByKey for pass C instead of
   * PT16SassyLookup::lookup. A thin wrapper around bucketedScan below,
   * collapsing its full KmerLookupResult per position down to just
   * (position, length) -- exactly what bucketedScan's own results already
   * hold in .match_position/.match_length, so this changes nothing about
   * what gets reported (see lrf_ms/v2_bucket_equiv_test-style checks).
   */
  MatchingStatistics computeMatchingStatisticsBucketed(
      const std::vector<T1>& input) {
    const std::vector<KmerLookupResult> results = bucketedScan(input);

    MatchingStatistics ms(results.size());

    for (std::size_t i = 0; i < results.size(); ++i) {
      ms[i] = {results[i].match_position, results[i].match_length};
    }

    return ms;
  }

  /**
   * The bucketed scan itself, returning the FULL KmerLookupResult per
   * position (indexed by original text position, same order
   * computeMatchingStatisticsBucketed's own output is in) instead of
   * collapsing it to (position, length) -- so a caller can feed this
   * through chain_extend.hpp/ms_tools.hpp afterward and get this format's
   * fastest lookup floor (established earlier: bucketing beats both
   * sequential scanning and the sassy table on raw lookup speed) PLUS
   * chain extension on top, instead of just the scan-only answer.
   */
  std::vector<KmerLookupResult> bucketedScan(const std::vector<T1>& input) {
    using clock = std::chrono::steady_clock;

    const auto before = parser_->stats();

    const std::size_t n = input.size();
    std::vector<KmerLookupResult> results(n);

    const std::size_t kmer_positions =
        n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;

    // ---------- Prebucket: pack every key, count per table bucket ----------

    const auto prebucket_start = clock::now();

    std::vector<std::uint32_t> keys(kmer_positions);
    auto count = std::make_unique<std::array<std::uint32_t, NUMBER_OF_BUCKETS>>();
    count->fill(0);

    // Rolling key, same trick as the sequential scan above: only i=0 pays
    // a full encode.
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
    const double prebucket_ms = std::chrono::duration<double, std::milli>(
                                    clock::now() - prebucket_start)
                                    .count();

    // ---------- Bucket: prefix sum + counting-sort placement ----------

    const auto bucket_start = clock::now();

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
    const double bucket_ms = std::chrono::duration<double, std::milli>(
                                 clock::now() - bucket_start)
                                 .count();

    // ---------- Probe: the lookups themselves, in bucket order ----------

    const auto probe_start = clock::now();

    for (std::size_t slot = 0; slot < kmer_positions; ++slot) {
      const std::uint32_t i = order[slot];
      results[i] = parser_->lookupKmerByKey(input, i, keys[i]);
    }

    phase_barrier(results.data());
    const double probe_ms = std::chrono::duration<double, std::milli>(
                                clock::now() - probe_start)
                                .count();

    // ---------- Tail: not reordered (at most 15 of them) ----------

    const auto tail_start = clock::now();

    // The last 16-mer key rolled past the end: each `<< 2` drops the first
    // character and pads with a zero one, giving the padded tail key.
    std::uint32_t tail_key =
        kmer_positions > 0 ? keys[kmer_positions - 1] << 2U
                           : (n > 0 ? encode_tail(input, 0) : 0);

    for (std::size_t i = kmer_positions; i < n; ++i, tail_key <<= 2U) {
      const auto tail = parser_->lookupTailByKey(
          tail_key, static_cast<std::uint32_t>(n - i));
      results[i].found = false;
      results[i].match_position = tail.match_position;
      results[i].match_length = tail.match_length;
    }

    phase_barrier(results.data());
    const double tail_ms = std::chrono::duration<double, std::milli>(
                               clock::now() - tail_start)
                               .count();

    const auto after = parser_->stats();

    diagnostics_ = phase_diagnostics({{"prebucket", prebucket_ms},
                                      {"bucket", bucket_ms},
                                      {"probe", probe_ms},
                                      {"tail", tail_ms}}) +
                   bucketSearchDiagnostics(before, after);
    search_composition_ = {true, after.singleton_hits - before.singleton_hits,
                           after.range_hits - before.range_hits,
                           after.misses - before.misses};

    return results;
  }

  Diagnostics diagnostics() const { return diagnostics_; }

  // Both compute methods above populate this identically (the raw search
  // classification does not depend on lookup ORDER), so it does not matter
  // which one ran last.
  SearchComposition searchComposition() const { return search_composition_; }

 private:
  std::unique_ptr<PT16RLZParser<T1, T2>> parser_;
  Diagnostics diagnostics_;
  SearchComposition search_composition_;
};

/**
 * Benchmark-only adapter forwarding to PT16ScanMS::computeMatchingStatisticsBucketed,
 * so the bucketed non-sassy scan gets its own row, the same way
 * PT16SassyBucketScanMS does for the sassy one.
 */
template <typename T1, typename T2>
class PT16BucketScanMS {
 public:
  PT16BucketScanMS(const std::vector<T1>& reference,
                   const std::vector<T2>& suffix_array,
                   const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16BucketScanMS(const PT16BucketScanMS&) = delete;
  PT16BucketScanMS& operator=(const PT16BucketScanMS&) = delete;

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    return impl_.computeMatchingStatisticsBucketed(input);
  }

  Diagnostics diagnostics() const { return impl_.diagnostics(); }

  SearchComposition searchComposition() const {
    return impl_.searchComposition();
  }

  // Scan-only: every hit is capped at exactly 16, never extended, so this
  // never matches the baseline's real matching statistics by design.
  bool exactExpected() const { return false; }

 private:
  PT16ScanMS<T1, T2> impl_;
};

/**
 * The bucketed v2 scan (the fastest way established so far to build the
 * raw KmerLookupResult list -- it beat both sequential scanning and the
 * sassy table on raw lookup speed) fed through multi-step position-only
 * chain extension (backwardChainExtend, ms_tools.hpp): the v2-format
 * counterpart of PT16SassyBackwardChainMS below, expected to be both the
 * fastest and (like it) exact against brute force.
 */
template <typename T1, typename T2>
class PT16V2BucketChainMS {
 public:
  PT16V2BucketChainMS(const std::vector<T1>& reference,
                      const std::vector<T2>& suffix_array,
                      const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16V2BucketChainMS(const PT16V2BucketChainMS&) = delete;
  PT16V2BucketChainMS& operator=(const PT16V2BucketChainMS&) = delete;

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
  PT16ScanMS<T1, T2> impl_;
  double chain_ms_ = 0.0;  // last backwardChainExtend call
};

/**
 * sortedKmerScan (sorted_kmer_scan.hpp) over the v2 table: full sort by
 * the 32-bit key instead of bucket grouping alone, via a lookup callable
 * that adapts PT16RLZParser::lookupKmerByKey/lookupKmerOrTail to
 * sortedKmerScan's black-box interface -- no v2-specific scan logic lives
 * here, only that one-line adapter. Its own "-scan" row, the counterpart
 * of pt16-v2-bucket-scan, so the two orderings' raw lookup cost can be
 * read side by side.
 */
template <typename T1, typename T2>
class PT16V2SortedScanMS {
 public:
  PT16V2SortedScanMS(const std::vector<T1>& reference,
                     const std::vector<T2>& suffix_array,
                     const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16V2SortedScanMS(const PT16V2SortedScanMS&) = delete;
  PT16V2SortedScanMS& operator=(const PT16V2SortedScanMS&) = delete;

  std::vector<KmerLookupResult> sortedScan(const std::vector<T1>& input) {
    const auto& parser = impl_.parser();

    std::vector<KmerLookupResult> results = sortedKmerScan(
        input,
        [&](const std::vector<T1>& in, std::size_t i, std::uint32_t key) {
          return parser.lookupKmerByKey(in, i, key);
        },
        [&](const std::vector<T1>& in, std::size_t i, std::uint32_t key) {
          return parser.lookupTailByKey(
              key, static_cast<std::uint32_t>(in.size() - i));
        },
        &timings_);

    const std::size_t n = input.size();
    const std::size_t kmer_positions =
        n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;
    //search_composition_ = classifySearchComposition(results, kmer_positions);

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

  Diagnostics diagnostics() const {
    return sortedScanDiagnostics(timings_);
  }

  SearchComposition searchComposition() const { return search_composition_; }

  // Scan-only: every hit is capped at exactly 16, never extended, so this
  // never matches the baseline's real matching statistics by design.
  bool exactExpected() const { return false; }

 private:
  PT16ScanMS<T1, T2> impl_;
  SearchComposition search_composition_;
  SortedScanTimings timings_;
};

/**
 * Same sorted scan, fed through multi-step chain extension
 * (backwardChainExtend, ms_tools.hpp) -- the sorted-scan counterpart of
 * pt16-v2-bucket-chain-multi, exact against brute force for the same
 * reason (backwardChainExtend only ever reads `results`, never how it was
 * built).
 */
template <typename T1, typename T2>
class PT16V2SortedChainMS {
 public:
  PT16V2SortedChainMS(const std::vector<T1>& reference,
                      const std::vector<T2>& suffix_array,
                      const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16V2SortedChainMS(const PT16V2SortedChainMS&) = delete;
  PT16V2SortedChainMS& operator=(const PT16V2SortedChainMS&) = delete;

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
  PT16V2SortedScanMS<T1, T2> impl_;
  double chain_ms_ = 0.0;  // last backwardChainExtend call
};

// ---------- PT16 sassy + chain extension ----------

/**
 * A real matching-statistics answer, not capped at 16: the bucketed sassy
 * scan (PT16SassyMS::bucketedScan -- rolling-key, bucket-ordered, the
 * fastest way established so far to build the raw KmerLookupResult list
 * for this table format) fed through chainExtend (chain_extend.hpp) --
 * one-step position-only chain extension, no reference read anywhere.
 * Because it only looks one position ahead, a true match longer than 17 is
 * still only extended by one character here; see pt16-sassy-chain-multi
 * for the version that chases the whole chain.
 *
 * Previously built on a free-standing sequential scan (re-encoding every
 * 16-mer from scratch, not bucket-ordered): that was needlessly the
 * slowest way to build this list, not a property of chaining a sassy
 * result itself, so it has been dropped in favor of bucketedScan.
 */
class PT16SassyChainMS {
 public:
  PT16SassyChainMS(const std::vector<unsigned char>& reference,
                   const std::vector<std::uint32_t>& suffix_array,
                   const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16SassyChainMS(const PT16SassyChainMS&) = delete;
  PT16SassyChainMS& operator=(const PT16SassyChainMS&) = delete;

  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
    return chainExtend(impl_.bucketedScan(input));
  }

  Diagnostics diagnostics() const { return impl_.lastDiagnostics(); }

  SearchComposition searchComposition() const {
    return impl_.lastSearchComposition();
  }

  // One-step chainExtend only ever extends a hit by one character, so it
  // matches the baseline wherever the true match is <= 17 and falls short
  // for anything longer -- not expected to be exact by design.
  bool exactExpected() const { return false; }

 private:
  PT16SassyMS impl_;
};

/**
 * Same bucketed scan, but multi-step position-only chain extension
 * (backwardChainExtend, ms_tools.hpp) instead: proven exact against
 * brute force, not just a sound lower bound (lrf_ms/chain_extend_test.cpp),
 * so this is expected to match the baseline's lengths exactly -- the sassy
 * counterpart of pt16-v2-bucket-chain-multi, now built the same way (best
 * available scan first, chain on top), for a fair comparison between the
 * two table formats' chain-extension cost.
 */
class PT16SassyBackwardChainMS {
 public:
  PT16SassyBackwardChainMS(const std::vector<unsigned char>& reference,
                           const std::vector<std::uint32_t>& suffix_array,
                           const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16SassyBackwardChainMS(const PT16SassyBackwardChainMS&) = delete;
  PT16SassyBackwardChainMS& operator=(const PT16SassyBackwardChainMS&) = delete;

  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
    return timedBackwardChainExtend(impl_.bucketedScan(input), chain_ms_);
  }

  // The scan's own phases, plus the chain extension run on its result.
  Diagnostics diagnostics() const {
    Diagnostics diagnostics = impl_.lastDiagnostics();
    diagnostics.add_phase("chain", chain_ms_);
    return diagnostics;
  }

  SearchComposition searchComposition() const {
    return impl_.lastSearchComposition();
  }

 private:
  PT16SassyMS impl_;
  double chain_ms_ = 0.0;  // last backwardChainExtend call
};

/**
 * sortedKmerScan (sorted_kmer_scan.hpp) over the sassy table: same
 * function as PT16V2SortedScanMS above, only the adapter differs --
 * PT16SassyLookup::lookup ignores the `input`/`position` sortedKmerScan
 * passes it (its table is self-contained), and lookup_tail returns a
 * TailResult, not a KmerLookupResult, so the tail adapter wraps it. No
 * sassy-specific scan logic here either. Its own "-scan" row, the
 * counterpart of pt16-sassy-bucket-scan.
 */
class PT16SassySortedScanMS {
 public:
  PT16SassySortedScanMS(const std::vector<unsigned char>& reference,
                        const std::vector<std::uint32_t>& suffix_array,
                        const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16SassySortedScanMS(const PT16SassySortedScanMS&) = delete;
  PT16SassySortedScanMS& operator=(const PT16SassySortedScanMS&) = delete;

  std::vector<KmerLookupResult> sortedScan(
      const std::vector<unsigned char>& input) {
    const PT16SassyLookup& lookup = impl_.lookup();
    const PT16SassyLookup::Stats before = lookup.stats();

    std::vector<KmerLookupResult> results = sortedKmerScan(
        input,
        [&](const std::vector<unsigned char>&, std::size_t,
            std::uint32_t key) { return lookup.lookup(key); },
        [&](const std::vector<unsigned char>& in, std::size_t i,
            std::uint32_t key) -> KmerLookupResult {
          const auto tail = lookup.lookup_tail(
              key, static_cast<std::uint32_t>(in.size() - i));
          KmerLookupResult result;
          result.found = false;
          result.match_position = tail.match_position;
          result.match_length = tail.match_length;
          return result;
        },
        &timings_);

    diagnostics_ = sortedScanDiagnostics(timings_) +
                   miss_diagnostics(before, lookup.stats());

    const std::size_t n = input.size();
    const std::size_t kmer_positions =
        n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;
    // search_composition_ = classifySearchComposition(results, kmer_positions);

    return results;
  }

  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
    const std::vector<KmerLookupResult> results = sortedScan(input);

    MatchingStatistics ms(results.size());
    for (std::size_t i = 0; i < results.size(); ++i) {
      ms[i] = {results[i].match_position, results[i].match_length};
    }

    return ms;
  }

  Diagnostics diagnostics() const { return diagnostics_; }

  SearchComposition searchComposition() const { return search_composition_; }

  // Scan-only: every hit is capped at exactly 16, never extended, so this
  // never matches the baseline's real matching statistics by design.
  bool exactExpected() const { return false; }

 private:
  PT16SassyMS impl_;
  SearchComposition search_composition_;
  SortedScanTimings timings_;
  Diagnostics diagnostics_;
};

/**
 * Same sorted scan, fed through multi-step chain extension
 * (backwardChainExtend, ms_tools.hpp) -- the sorted-scan counterpart of
 * pt16-sassy-chain-multi.
 */
class PT16SassySortedChainMS {
 public:
  PT16SassySortedChainMS(const std::vector<unsigned char>& reference,
                         const std::vector<std::uint32_t>& suffix_array,
                         const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16SassySortedChainMS(const PT16SassySortedChainMS&) = delete;
  PT16SassySortedChainMS& operator=(const PT16SassySortedChainMS&) = delete;

  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
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
  PT16SassySortedScanMS impl_;
  double chain_ms_ = 0.0;  // last backwardChainExtend call
};

// ---------- Registry ----------

/**
 * Builds every implementation under test, in order. Index 0 is the
 * baseline.
 *
 * To add a variant, include its header above and add one line here:
 *
 *   add_implementation<MyVariant<Symbol, SAType>>(
 *       implementations, "my-variant", reference, suffix_array);
 *
 * Constructor arguments beyond the reference and suffix array are
 * forwarded as given, for example:
 *
 *   add_implementation<PT16MS<Symbol, SAType>>(
 *       implementations, "pt16-ms", reference, suffix_array, table_path);
 */
inline Implementations build_implementations(
    const std::vector<Symbol>& reference,
    const std::vector<SAType>& suffix_array,
    [[maybe_unused]] const std::string& pt16_table) {
  Implementations implementations;

  add_implementation<LRFMS<Symbol, SAType>>(implementations, "lrf-ms",
                                            reference, suffix_array);

  // --- Add full implementations below this line ---
  //
  // The PT16 table variants are not registered here: they share one
  // pipeline (keys, bucket/sorted order, probe, chain -- see
  // probe_pipeline.hpp), and ms_main runs the shared stages once per file
  // and only the probe per variant. Register a new table variant in
  // build_probers (probe_pipeline.hpp) instead. This list is for
  // implementations that compute complete matching statistics on their
  // own, like the baseline.
  //
  // The classes that used to be registered here (PT16V2BucketChainMS,
  // PT16SassyBackwardChainMS, PT16FastMissSortedChainMS, ... and their
  // scan-only counterparts) are still defined below and in
  // pt16_sassy_ms.hpp / pt16_fastmiss_ms.hpp, for the tests.

  return implementations;
}

}  // namespace msbench
