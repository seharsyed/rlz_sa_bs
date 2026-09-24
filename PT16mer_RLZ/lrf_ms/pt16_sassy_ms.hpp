#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../variants/pt16_build_sassy.hpp"
#include "../variants/pt16_sassy.hpp"

/**
 * Table construction and raw lookup access over the self-contained sassy
 * PT16 table (pt16_sassy.hpp). Matching statistics themselves are never
 * computed inside this class: it only builds/loads the table and exposes
 * scan methods (scanPass1, bucketedScan) that produce a KmerLookupResult
 * per input position; turning that into matching statistics is left to
 * chain_extend.hpp/ms_tools.hpp (format-agnostic) so every chain-extension
 * variant runs on whichever scan is currently fastest, not a copy of it
 * hand-rolled in here that would drift as the scan gets tuned further.
 *
 * PT16SassyLookup does not keep its own copy of the reference (lookup and
 * find_longest_matching_factor both take it as an argument instead), so this
 * class keeps a pointer to it, exactly as PT16RLZParser does.
 */
class PT16SassyMS {
 public:
  PT16SassyMS(const std::vector<unsigned char>& reference,
             const std::vector<std::uint32_t>& suffix_array,
             const std::string& table_path)
      : reference_(&reference) {
    // Always rebuild, so a stale table from another reference or another
    // version of the table format is never used (matching PT16MS).
    std::filesystem::remove(table_path);

    build_pt16_sassy_table(reference, suffix_array, table_path);

    lookup_ = std::make_unique<PT16SassyLookup>(table_path);
  }

  PT16SassyMS(const PT16SassyMS&) = delete;
  PT16SassyMS& operator=(const PT16SassyMS&) = delete;

  /**
   * The mandatory first scan alone: one lookup()/lookup_tail() per input
   * position, no pass 2 (no chain extension), so every hit is reported at
   * exactly 16 -- whatever lookup itself found, untouched.
   *
   * This is not meant to be a good matching-statistics answer; it is a
   * lower bound on the cost ANY sassy-based MS variant must pay, since
   * that first scan is unavoidable. It exists so this floor can be
   * compared on its own, before spending effort on chaining/extension: if
   * this alone is already slower than the baseline, no amount of extension
   * work on top of it will change that.
   */
  MatchingStatistics computeMatchingStatisticsScanOnly(
      const std::vector<unsigned char>& input) {
    const PT16SassyLookup::Stats before = lookup_->stats();

    const std::vector<ScanEntry> entries = scanPass1(input);

    const PT16SassyLookup::Stats after = lookup_->stats();
    diagnostics_ = bucketSearchDiagnostics(before, after) +
                   miss_diagnostics(before, after);
    search_composition_ = {true, after.singleton_hits - before.singleton_hits,
                           after.range_hits - before.range_hits,
                           after.misses - before.misses};

    MatchingStatistics ms;
    ms.reserve(entries.size());

    for (const ScanEntry& entry : entries) {
      ms.emplace_back(entry.match_position, entry.match_length);
    }

    return ms;
  }

  /**
   * Same answer as computeMatchingStatisticsScanOnly (pass 1 alone, no
   * chain extension), but the lookups are done in a different ORDER: first
   * a linear scan over the input buckets every position by the top 16 bits
   * of its 16-mer (its leading 8 characters) -- exactly the bucket a
   * PT16SassyLookup dispatch on that key would use, since it is the same
   * `key >> LOW_BITS` -- storing the full 16-mer alongside so the input is
   * never read again. Then the lookups themselves run bucket by bucket,
   * so every lookup that lands in the same table bucket happens back to
   * back, instead of scattered across the input in whatever order those
   * 16-mers happened to occur. The results are still written back to
   * `entries` in ORIGINAL text-position order, so that write is the random
   * one here, not the lookup itself.
   *
   * The bet: since a table bucket's L_/sampled_sa_ slice is small (tens of
   * entries) but sits inside arrays far bigger than cache, the first
   * lookup into a given bucket is a cold miss either way, but grouping
   * means every LATER lookup into that same bucket -- from a different,
   * possibly far-away input position -- finds it still warm, instead of
   * paying a fresh cold miss each time as scanPass1's input-order lookups
   * do. This only pays off if the input actually revisits table buckets;
   * on top of whatever it saves, it now also pays the grouping pass
   * itself, and needs 8 extra bytes per 16-mer position (its key and its
   * slot) while it runs.
   */
  MatchingStatistics computeMatchingStatisticsScanOnlyBucketed(
      const std::vector<unsigned char>& input) {
    const std::vector<KmerLookupResult> results = bucketedScan(input);

    MatchingStatistics ms(results.size());

    for (std::size_t i = 0; i < results.size(); ++i) {
      ms[i] = {results[i].match_position, results[i].match_length};
    }

    return ms;
  }

  // Wall time spent in each phase of scanPass1Bucketed the last time it
  // ran: one clock::now() pair per phase (4 total), not per lookup, so
  // this stays negligible next to the work it measures.
  struct BucketScanTimings {
    std::int64_t prebucket_ns = 0;  // pack every key, count per bucket
    std::int64_t bucket_ns = 0;     // prefix sum + counting-sort placement
    std::int64_t probe_ns = 0;      // the lookups themselves, in bucket order
    std::int64_t tail_ns = 0;       // the last <16 characters, unreordered
  };

  const BucketScanTimings& lastBucketScanTimings() const {
    return bucket_scan_timings_;
  }

  // The last call's diagnostics (from computeMatchingStatisticsScanOnly/
  // ScanOnlyBucketed above), forwarded by PT16SassyScanMS/
  // PT16SassyBucketScanMS::diagnostics() for ms_main.cpp to sum per
  // implementation.
  const Diagnostics& lastDiagnostics() const { return diagnostics_; }

  // The last call's raw-lookup classification (singleton/range/miss),
  // identical whichever of the two compute methods above populated it
  // (lookup order does not change the classification). See
  // SearchComposition (pt16_utils.hpp) for what this is for.
  const SearchComposition& lastSearchComposition() const {
    return search_composition_;
  }

  // The reference held for a later pass (extending beyond what pure
  // position-chaining can resolve will need it, since PT16SassyLookup keeps
  // no copy of its own).
  const std::vector<unsigned char>& reference() const { return *reference_; }

  const PT16SassyLookup& lookup() const { return *lookup_; }

  /**
   * The bucketed scan itself, returning the FULL KmerLookupResult per
   * position (in original text-position order) instead of collapsing it to
   * (position, length) -- the sassy counterpart of PT16ScanMS::bucketedScan
   * (ms_variants.hpp), and now the ONLY scan a chain-extension variant
   * builds on for this table format: it is both bucket-ordered (established
   * fastest for raw lookup speed) and rolling-key (one encode_16mer total,
   * not one per position), unlike the free-standing sassyScan helper this
   * replaced, which re-encoded every 16-mer from scratch. See that method's
   * doc comment on computeMatchingStatisticsScanOnlyBucketed above for the
   * bucketing idea; this is the same three-pass mechanics (group by bucket,
   * counting-sort into `order`, then look up bucket by bucket), just
   * writing a KmerLookupResult per slot -- lookup_->lookup(key) already
   * returns one fully populated, hit or not, so no found/not-found
   * branching is needed here (contrast the old ScanEntry-based version this
   * replaced).
   */
  std::vector<KmerLookupResult> bucketedScan(
      const std::vector<unsigned char>& input) {
    using clock = std::chrono::steady_clock;

    const PT16SassyLookup::Stats before = lookup_->stats();

    const std::size_t n = input.size();
    std::vector<KmerLookupResult> results(n);

    bucket_scan_timings_ = BucketScanTimings{};

    const std::size_t kmer_positions =
        n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;

    // ---------- Prebucket: pack every key, count per table bucket ----------

    const auto prebucket_start = clock::now();

    std::vector<std::uint32_t> keys(kmer_positions);
    auto count = std::make_unique<std::array<std::uint32_t, NUMBER_OF_BUCKETS>>();
    count->fill(0);

    // Rolling key: only i=0 pays a full encode_16mer.
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
    bucket_scan_timings_.prebucket_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
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
    bucket_scan_timings_.bucket_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now() - bucket_start)
            .count();

    // ---------- Probe: the lookups themselves, in bucket order ----------

    const auto probe_start = clock::now();

    for (std::size_t slot = 0; slot < kmer_positions; ++slot) {
      const std::uint32_t i = order[slot];
      results[i] = lookup_->lookup(keys[i]);
    }

    phase_barrier(results.data());
    bucket_scan_timings_.probe_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now() - probe_start)
            .count();

    // ---------- Tail: not reordered (at most 15 of them) ----------

    const auto tail_start = clock::now();

    for (std::size_t i = kmer_positions; i < n; ++i) {
      const auto tail = lookup_->lookup_tail(input, i);
      results[i].found = false;
      results[i].match_position = tail.match_position;
      results[i].match_length = tail.match_length;
    }

    phase_barrier(results.data());
    bucket_scan_timings_.tail_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now() - tail_start)
            .count();

    const PT16SassyLookup::Stats after = lookup_->stats();
    diagnostics_ = bucketedScanDiagnostics(before, after);
    search_composition_ = {true, after.singleton_hits - before.singleton_hits,
                           after.range_hits - before.range_hits,
                           after.misses - before.misses};

    return results;
  }

 private:
  // One position's raw 16-mer/tail lookup, kept only long enough to check
  // whether the PREVIOUS position's match extends into this one (pass 2
  // above). `positions` is a view into PT16SassyLookup's own storage (valid
  // for as long as `lookup_` is), so nothing here is copied out of it.
  struct ScanEntry {
    // True once match_length/match_position are the final answer for this
    // position: a short match or a tail is resolved immediately in pass 1;
    // a hit becomes resolved once pass 2 has decided whether it extends.
    bool resolved = false;

    std::uint32_t match_length = 0;
    std::uint32_t match_position = 0;

    // This position's own occurrences, used only to check whether the
    // PREVIOUS position's match extends into it. 0/empty for a resolved
    // (short/tail) entry, since there is nothing to chain into.
    std::uint32_t count = 0;
    std::span<const std::uint32_t> positions;
  };

  // Pass 1: one raw lookup per position, no chaining yet. Shared by the
  // full computeMatchingStatistics (which runs pass 2 on top of this) and
  // computeMatchingStatisticsScanOnly (which returns exactly this). No
  // internal timing: whatever wraps the call (the benchmark's own
  // time_repeated) already measures this end to end, and adding a
  // clock::now() pair per position here only inflates that measurement
  // with its own overhead.
  //
  // The 16-mer key is a ROLLING window, not a fresh encode_16mer per
  // position: encode_16mer/pack_16mer exist for lzFactorize's skip-ahead
  // calls, where the next call is NOT at the next position, so nothing can
  // be reused. Here every position from 0 to n-KMER_LENGTH is visited in
  // order, so position i's key and position i+1's key share 15 of their 16
  // characters. Only i=0 pays the full 16-character encode; every later
  // key is `(previous key << 2) | new trailing character`: the shift
  // drops the oldest character off the top of the 32-bit key (nothing
  // masks it away -- the shift itself discards whatever no longer fits),
  // and the new character (input[i + KMER_LENGTH - 1], the last character
  // of the new window) fills the two bits the shift just freed at the
  // bottom. Same encoding encode_16mer produces, just built incrementally.
  std::vector<ScanEntry> scanPass1(const std::vector<unsigned char>& input) {
    const std::size_t n = input.size();
    std::vector<ScanEntry> entries(n);

    std::uint32_t key = n >= KMER_LENGTH ? encode_16mer(input, 0) : 0;

    for (std::size_t i = 0; i < n; ++i) {
      if (n - i >= KMER_LENGTH) {
        if (i > 0) {
          key = (key << 2U) |
                alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
        }

        const auto result = lookup_->lookup(key);

        if (!result.found) {
          // A short match: the 16-mer itself does not occur anywhere, so
          // match_length is already the longest possible match at i.
          entries[i] = {true, result.match_length, result.match_position,
                       0, {}};
        } else {
          entries[i].match_length = KMER_LENGTH;
          entries[i].match_position = result.match_position;
          entries[i].count = result.count;
          entries[i].positions = result.positions;
        }
      } else {
        // Fewer than 16 characters remain: the input ends here, so nothing
        // can extend this further, whatever lookup_tail finds.
        const auto tail = lookup_->lookup_tail(input, i);
        entries[i] = {true, tail.match_length, tail.match_position, 0, {}};
      }
    }

    return entries;
  }

  // One line: how many of this call's bucket dispatches searched their
  // bucket linearly vs. with std::lower_bound (see
  // PT16SassyLookup::Stats::linear_bucket_searches/binary_bucket_searches).
  static Diagnostics bucketSearchDiagnostics(
      const PT16SassyLookup::Stats& before,
      const PT16SassyLookup::Stats& after) {
    return counter_diagnostics(
        "bucket search",
        {{"linear",
          after.linear_bucket_searches - before.linear_bucket_searches},
         {"binary",
          after.binary_bucket_searches - before.binary_bucket_searches}});
  }

  // Same, plus the 4-phase timing breakdown from the bucketed scan.
  Diagnostics bucketedScanDiagnostics(
      const PT16SassyLookup::Stats& before,
      const PT16SassyLookup::Stats& after) const {
    const auto& t = bucket_scan_timings_;

    return phase_diagnostics(
               {{"prebucket", static_cast<double>(t.prebucket_ns) / 1e6},
                {"bucket", static_cast<double>(t.bucket_ns) / 1e6},
                {"probe", static_cast<double>(t.probe_ns) / 1e6},
                {"tail", static_cast<double>(t.tail_ns) / 1e6}}) +
           bucketSearchDiagnostics(before, after) +
           miss_diagnostics(before, after);
  }

  const std::vector<unsigned char>* reference_ = nullptr;
  std::unique_ptr<PT16SassyLookup> lookup_;
  BucketScanTimings bucket_scan_timings_;
  Diagnostics diagnostics_;
  SearchComposition search_composition_;
};

/**
 * Benchmark-only adapter around PT16SassyMS: same construction (builds and
 * loads the sassy table), but computeMatchingStatistics here calls the
 * scan-only pass, not the full one. Registering this as its own
 * implementation in ms_variants.hpp gives the scan floor its own row/timing
 * in the benchmark, separate from (and comparable to) whatever the full
 * chain-extending variant costs once that is registered too.
 */
class PT16SassyScanMS {
 public:
  PT16SassyScanMS(const std::vector<unsigned char>& reference,
                  const std::vector<std::uint32_t>& suffix_array,
                  const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16SassyScanMS(const PT16SassyScanMS&) = delete;
  PT16SassyScanMS& operator=(const PT16SassyScanMS&) = delete;

  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
    return impl_.computeMatchingStatisticsScanOnly(input);
  }

  Diagnostics diagnostics() const { return impl_.lastDiagnostics(); }

  SearchComposition searchComposition() const {
    return impl_.lastSearchComposition();
  }

  // Scan-only: every hit is capped at exactly 16, never extended, so this
  // never matches the baseline's real matching statistics by design.
  bool exactExpected() const { return false; }

 private:
  PT16SassyMS impl_;
};

/**
 * Same idea as PT16SassyScanMS, but computeMatchingStatistics here calls
 * the BUCKETED scan-only pass (computeMatchingStatisticsScanOnlyBucketed):
 * lookups grouped and run in table-bucket order rather than input order.
 * Registering it as its own implementation gives it its own row/timing,
 * directly comparable to "pt16-sassy-scan" on the same input.
 */
class PT16SassyBucketScanMS {
 public:
  PT16SassyBucketScanMS(const std::vector<unsigned char>& reference,
                        const std::vector<std::uint32_t>& suffix_array,
                        const std::string& table_path)
      : impl_(reference, suffix_array, table_path) {}

  PT16SassyBucketScanMS(const PT16SassyBucketScanMS&) = delete;
  PT16SassyBucketScanMS& operator=(const PT16SassyBucketScanMS&) = delete;

  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
    return impl_.computeMatchingStatisticsScanOnlyBucketed(input);
  }

  Diagnostics diagnostics() const { return impl_.lastDiagnostics(); }

  SearchComposition searchComposition() const {
    return impl_.lastSearchComposition();
  }

  // Scan-only: every hit is capped at exactly 16, never extended, so this
  // never matches the baseline's real matching statistics by design.
  bool exactExpected() const { return false; }

 private:
  PT16SassyMS impl_;
};
