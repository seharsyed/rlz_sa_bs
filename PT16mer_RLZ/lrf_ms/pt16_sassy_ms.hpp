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

#include "../pt16_build_sassy.hpp"
#include "../pt16_sassy.hpp"

/**
 * Matching statistics over the self-contained sassy PT16 table
 * (pt16_sassy.hpp). This is the initialization step only: the constructor
 * builds and loads the table, matching PT16MS's shape in ms_variants.hpp so
 * this can be registered the same way once computeMatchingStatistics is
 * written.
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
   * A first, deliberately incomplete sketch. Extends a hit by at most one
   * character, and only by following a chain of exact 16-mer occurrences:
   * it compares the reference POSITIONS lookup() returns, never reference
   * CHARACTERS -- the reference is not read at all here.
   *
   * Why comparing positions is enough: if the 16-mer at i occurs at
   * reference position p, and the 16-mer at i+1 occurs at p+1, then
   * input[i..i+16) == ref[p..p+16) (the first fact) and in particular
   * input[i+16] == ref[p+16] (the last character of the second fact), so
   * input[i..i+17) == ref[p..p+17): the match at i is at least 17 long.
   * This holds regardless of how many OTHER occurrences either 16-mer has,
   * so the successor check also looks inside a range, not just a singleton.
   *
   * What is NOT yet done (left for a later pass; every length reported here
   * is always a valid lower bound on the true matching statistic, so none
   * of this can make an already-reported answer wrong):
   *   - extending a chain by more than one step;
   *   - resolving a RANGE at the current position (left at length 16).
   */
  MatchingStatistics computeMatchingStatistics(
      const std::vector<unsigned char>& input) {
    const std::size_t n = input.size();

    std::vector<ScanEntry> entries = scanPass1(input);

    // ---------- Pass 2: extend a singleton hit by one character, if the
    // ---------- next position's own occurrences include position + 1.

    for (std::size_t i = 0; i < n; ++i) {
      ScanEntry& entry = entries[i];

      if (!entry.resolved && entry.count == 1 && i + 1 < n) {
        const std::uint32_t wanted = entry.match_position + 1;
        const ScanEntry& next = entries[i + 1];
        bool extends = false;

        if (next.count == 1) {
          extends = next.match_position == wanted;
        } else {
          for (const std::uint32_t candidate : next.positions) {
            if (candidate == wanted) {
              extends = true;
              break;
            }
          }
        }

        if (extends) {
          ++entry.match_length;
        }
      }

      entry.resolved = true;
    }

    // ---------- Collect ----------

    MatchingStatistics ms;
    ms.reserve(n);

    for (const ScanEntry& entry : entries) {
      ms.emplace_back(entry.match_position, entry.match_length);
    }

    return ms;
  }

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

    diagnostics_ = formatBucketSearchCounts(before, lookup_->stats());

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
    const PT16SassyLookup::Stats before = lookup_->stats();

    const std::vector<ScanEntry> entries = scanPass1Bucketed(input);

    diagnostics_ = formatBucketedScanTimings(before, lookup_->stats());

    MatchingStatistics ms;
    ms.reserve(entries.size());

    for (const ScanEntry& entry : entries) {
      ms.emplace_back(entry.match_position, entry.match_length);
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

  // The last call's diagnostic text (built, not printed, by
  // computeMatchingStatisticsScanOnly/ScanOnlyBucketed above), ready for
  // whoever prints the benchmark's per-implementation row to print this
  // right after it -- see PT16SassyScanMS/PT16SassyBucketScanMS::
  // diagnostics() and ms_main.cpp. Building a string instead of printing
  // directly means it can be attached to the right implementation's row
  // instead of always leaking out before it (this runs inside
  // computeMatchingStatistics, which finishes before that row is ever
  // printed).
  const std::string& lastDiagnostics() const { return diagnostics_; }

  // The reference held for a later pass (extending beyond what pure
  // position-chaining can resolve will need it, since PT16SassyLookup keeps
  // no copy of its own).
  const std::vector<unsigned char>& reference() const { return *reference_; }

  const PT16SassyLookup& lookup() const { return *lookup_; }

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
          entries[i].match_position =
              result.count == 1 ? result.position : result.positions.front();
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

  // Grouping pass + bucket-ordered lookup pass behind
  // computeMatchingStatisticsScanOnlyBucketed. See that method's doc
  // comment for the idea; this is the mechanics. No internal timing, for
  // the same reason as scanPass1.
  //
  // Pass A: one linear scan over the 16-mer positions, packing each key
  // (encode_16mer, same as scanPass1) and counting how many land in each
  // of the NUMBER_OF_BUCKETS table buckets (key >> LOW_BITS -- the same
  // split PT16SassyLookup::lookup uses internally).
  //
  // Pass B: prefix-sum those counts into slot offsets, then a second
  // linear scan places each position into its bucket's slice of `order`
  // (an ordinary counting sort, the same scheme build_H uses over the
  // reference at build time, here run over the input instead).
  //
  // Pass C: walk `order` bucket by bucket and look up the key already
  // stored for each slot -- the input itself is never read again -- and
  // write each result back into `entries` at its ORIGINAL position.
  //
  // Timed in 4 phases (one clock::now() pair each, not per lookup):
  // prebucket (pack + count), bucket (prefix sum + placement), probe (the
  // lookups), tail. See BucketScanTimings.
  std::vector<ScanEntry> scanPass1Bucketed(
      const std::vector<unsigned char>& input) {
    using clock = std::chrono::steady_clock;

    const std::size_t n = input.size();
    std::vector<ScanEntry> entries(n);

    bucket_scan_timings_ = BucketScanTimings{};

    const std::size_t kmer_positions =
        n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;

    // ---------- Prebucket: pack every key, count per table bucket ----------

    const auto prebucket_start = clock::now();

    std::vector<std::uint32_t> keys(kmer_positions);
    auto count = std::make_unique<std::array<std::uint32_t, NUMBER_OF_BUCKETS>>();
    count->fill(0);

    // Rolling key, same trick as scanPass1: only i=0 pays a full encode.
    std::uint32_t key = kmer_positions > 0 ? encode_16mer(input, 0) : 0;

    for (std::size_t i = 0; i < kmer_positions; ++i) {
      if (i > 0) {
        key = (key << 2U) |
              alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
      }

      keys[i] = key;
      ++(*count)[key >> LOW_BITS];
    }

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

    bucket_scan_timings_.bucket_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now() - bucket_start)
            .count();

    // ---------- Probe: the lookups themselves, in bucket order ----------

    const auto probe_start = clock::now();

    for (std::size_t slot = 0; slot < kmer_positions; ++slot) {
      const std::uint32_t i = order[slot];
      const auto result = lookup_->lookup(keys[i]);

      if (!result.found) {
        entries[i] = {true, result.match_length, result.match_position, 0,
                     {}};
      } else {
        entries[i].match_length = KMER_LENGTH;
        entries[i].match_position =
            result.count == 1 ? result.position : result.positions.front();
        entries[i].count = result.count;
        entries[i].positions = result.positions;
      }
    }

    bucket_scan_timings_.probe_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now() - probe_start)
            .count();

    // ---------- Tail: same as scanPass1, not reordered (at most 15 of them) ----------

    const auto tail_start = clock::now();

    for (std::size_t i = kmer_positions; i < n; ++i) {
      const auto tail = lookup_->lookup_tail(input, i);
      entries[i] = {true, tail.match_length, tail.match_position, 0, {}};
    }

    bucket_scan_timings_.tail_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock::now() - tail_start)
            .count();

    return entries;
  }

  // One line: how many of this call's bucket dispatches searched their
  // bucket linearly vs. with std::lower_bound (see
  // PT16SassyLookup::Stats::linear_bucket_searches/binary_bucket_searches).
  static std::string formatBucketSearchCounts(
      const PT16SassyLookup::Stats& before,
      const PT16SassyLookup::Stats& after) {
    std::ostringstream out;

    out << "        bucket search: linear="
        << (after.linear_bucket_searches - before.linear_bucket_searches)
        << " binary="
        << (after.binary_bucket_searches - before.binary_bucket_searches)
        << "\n";

    return out.str();
  }

  // Same, plus the 4-phase timing breakdown from the bucketed scan.
  std::string formatBucketedScanTimings(
      const PT16SassyLookup::Stats& before,
      const PT16SassyLookup::Stats& after) const {
    const auto& t = bucket_scan_timings_;
    std::ostringstream out;

    out << "        phases: prebucket "
        << static_cast<double>(t.prebucket_ns) / 1e6 << " ms"
        << "  bucket " << static_cast<double>(t.bucket_ns) / 1e6 << " ms"
        << "  probe " << static_cast<double>(t.probe_ns) / 1e6 << " ms"
        << "  tail " << static_cast<double>(t.tail_ns) / 1e6 << " ms\n";

    out << formatBucketSearchCounts(before, after);

    return out.str();
  }

  const std::vector<unsigned char>* reference_ = nullptr;
  std::unique_ptr<PT16SassyLookup> lookup_;
  BucketScanTimings bucket_scan_timings_;
  std::string diagnostics_;
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

  std::string diagnostics() const { return impl_.lastDiagnostics(); }

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

  std::string diagnostics() const { return impl_.lastDiagnostics(); }

 private:
  PT16SassyMS impl_;
};
