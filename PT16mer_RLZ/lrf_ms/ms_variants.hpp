#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../pt16_build_v2.hpp"
#include "../pt16_rlz_v2_interleaved.hpp"
#include "lrf_ms.hpp"
#include "ms_utils.hpp"
#include "pt16_sassy_ms.hpp"

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

  // Free-form diagnostic text from the last compute() call, already
  // terminated by its own trailing newline(s) (empty if there is none).
  // Printed by ms_main.cpp right after this implementation's own summary
  // row -- never before it, unlike printing straight to stderr from
  // inside compute() would (compute() always finishes before that row is
  // printed, so anything it prints directly ends up attached to whichever
  // row printed most recently, not its own).
  virtual std::string diagnostics() const { return {}; }
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
  std::string diagnostics() const override {
    if constexpr (requires(const Index& index) { index.diagnostics(); }) {
      return index_->diagnostics();
    } else {
      return {};
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
// struct without either needing to know about the other. Returned as a
// string (see MSImplementation::diagnostics) rather than printed directly,
// so the caller can attach it to the right implementation's own row.
template <typename Stats>
inline std::string formatBucketSearchCounts(const Stats& before,
                                            const Stats& after) {
  std::ostringstream out;

  out << "        bucket search: linear="
      << (after.linear_bucket_searches - before.linear_bucket_searches)
      << " binary="
      << (after.binary_bucket_searches - before.binary_bucket_searches)
      << "\n";

  return out.str();
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

    diagnostics_ = formatBucketSearchCounts(before, parser_->stats());

    return ms;
  }

  std::string diagnostics() const { return diagnostics_; }

 private:
  std::unique_ptr<PT16RLZParser<T1, T2>> parser_;
  std::string diagnostics_;
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

    std::uint32_t key = n >= KMER_LENGTH ? encode_16mer(input, 0) : 0;

    for (std::size_t i = 0; i < n; ++i) {
      if (n - i >= KMER_LENGTH) {
        if (i > 0) {
          key = (key << 2U) |
                alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
        }

        const auto result = parser_->lookupKmerByKey(input, i, key);
        ms.emplace_back(result.match_position, result.match_length);
      } else {
        const auto tail = parser_->lookupKmerOrTail(input, i);
        ms.emplace_back(tail.match_position, tail.match_length);
      }
    }

    diagnostics_ = formatBucketSearchCounts(before, parser_->stats());

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
   * PT16SassyLookup::lookup.
   */
  MatchingStatistics computeMatchingStatisticsBucketed(
      const std::vector<T1>& input) {
    using clock = std::chrono::steady_clock;

    const auto before = parser_->stats();

    const std::size_t n = input.size();
    MatchingStatistics ms(n);

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

    const double bucket_ms = std::chrono::duration<double, std::milli>(
                                 clock::now() - bucket_start)
                                 .count();

    // ---------- Probe: the lookups themselves, in bucket order ----------

    const auto probe_start = clock::now();

    for (std::size_t slot = 0; slot < kmer_positions; ++slot) {
      const std::uint32_t i = order[slot];
      const auto result = parser_->lookupKmerByKey(input, i, keys[i]);
      ms[i] = {result.match_position, result.match_length};
    }

    const double probe_ms = std::chrono::duration<double, std::milli>(
                                clock::now() - probe_start)
                                .count();

    // ---------- Tail: not reordered (at most 15 of them) ----------

    const auto tail_start = clock::now();

    for (std::size_t i = kmer_positions; i < n; ++i) {
      const auto tail = parser_->lookupKmerOrTail(input, i);
      ms[i] = {tail.match_position, tail.match_length};
    }

    const double tail_ms = std::chrono::duration<double, std::milli>(
                               clock::now() - tail_start)
                               .count();

    std::ostringstream out;

    out << "        phases: prebucket " << prebucket_ms << " ms"
        << "  bucket " << bucket_ms << " ms"
        << "  probe " << probe_ms << " ms"
        << "  tail " << tail_ms << " ms\n";

    out << formatBucketSearchCounts(before, parser_->stats());

    diagnostics_ = out.str();

    return ms;
  }

  std::string diagnostics() const { return diagnostics_; }

 private:
  std::unique_ptr<PT16RLZParser<T1, T2>> parser_;
  std::string diagnostics_;
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

  std::string diagnostics() const { return impl_.diagnostics(); }

 private:
  PT16ScanMS<T1, T2> impl_;
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
    const std::vector<SAType>& suffix_array, const std::string& pt16_table) {
  Implementations implementations;

  add_implementation<LRFMS<Symbol, SAType>>(implementations, "lrf-ms",
                                            reference, suffix_array);

  // --- Add variants below this line ---

  add_implementation<PT16MS<Symbol, SAType>>(implementations, "pt16-brute",
                                             reference, suffix_array,
                                             pt16_table);

  // The mandatory first scan alone (no SA narrowing/extension) over the
  // non-sassy H/L table -- the same idea as pt16-sassy-scan below, but for
  // the OTHER table format, so the two "just the lookup" costs can be
  // read side by side.
  add_implementation<PT16ScanMS<Symbol, SAType>>(
      implementations, "pt16-v2-scan", reference, suffix_array, pt16_table);

  // Same scan, grouped and run in table-bucket order -- the counterpart of
  // pt16-sassy-bucket-scan below, for the non-sassy table format.
  add_implementation<PT16BucketScanMS<Symbol, SAType>>(
      implementations, "pt16-v2-bucket-scan", reference, suffix_array,
      pt16_table);

  // The mandatory first scan alone (no chain extension) over the
  // self-contained sassy table -- a lower bound on any sassy-based MS
  // variant's cost. A different table file than pt16-brute's, since the
  // two formats are not interchangeable (see pt16_sassy.hpp).
  add_implementation<PT16SassyScanMS>(implementations, "pt16-sassy-scan",
                                      reference, suffix_array,
                                      pt16_table + ".sassy");

  // Same scan, but the lookups are grouped and run in table-bucket order
  // (by the input's own 8-mer prefix) instead of input order, to see
  // whether that ordering's cache locality pays for the grouping pass
  // itself. Shares the same table file as pt16-sassy-scan (both are
  // PT16SassyMS underneath).
  add_implementation<PT16SassyBucketScanMS>(
      implementations, "pt16-sassy-bucket-scan", reference, suffix_array,
      pt16_table + ".sassy");

  return implementations;
}

}  // namespace msbench
