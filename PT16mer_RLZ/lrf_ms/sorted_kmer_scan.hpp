#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "../pt16_utils.hpp"  // KMER_LENGTH, LOW_BITS, NUMBER_OF_BUCKETS,
                              // KmerLookupResult, SearchComposition,
                              // alphatab, encode_16mer

/**
 * Wall time spent in each phase of one sortedKmerScan call, in ms: one
 * clock::now() pair per phase, not per lookup, so measuring it costs
 * nothing noticeable.
 */
struct SortedScanTimings {
  double keys_ms = 0.0;        // rolling every 16-mer key
  double radix_low_ms = 0.0;   // pass 1: sort by the low 16 bits
  double radix_high_ms = 0.0;  // pass 2: stable sort by the high 16 bits
  double probe_ms = 0.0;       // the lookups, in sorted order
  double tail_ms = 0.0;        // the last < KMER_LENGTH positions
};

inline Diagnostics sortedScanDiagnostics(const SortedScanTimings& t) {
  return phase_diagnostics({{"keys", t.keys_ms},
                        {"radix-low", t.radix_low_ms},
                        {"radix-high", t.radix_high_ms},
                        {"probe", t.probe_ms},
                        {"tail", t.tail_ms}});
}

/**
 * sortedKmerScan: builds the raw KmerLookupResult list (in original text
 * position order) the same way every "-scan"/"-bucket-scan" variant in
 * ms_variants.hpp does, but oblivious to which PT16 table format is behind
 * it. Earlier scan builders (PT16ScanMS::bucketedScan,
 * PT16SassyMS::bucketedScan) each hard-code their own format's lookup call,
 * so improving the scan strategy meant editing both, by hand, in lockstep --
 * exactly the kind of copy that drifts. This function takes the lookup
 * itself as a black box instead: `lookup` and `lookup_tail` are the only
 * places a table format's own call appears, supplied once by the caller as
 * ordinary callables, so a future scan-strategy change is written here
 * ONCE and every caller -- whichever table format -- gets it automatically.
 *
 *   lookup(input, position, key) -> KmerLookupResult
 *     A real 16-mer lookup (>= KMER_LENGTH characters remain at
 *     `position`), `key` already packed by this function (see the rolling
 *     comment below). PT16SassyLookup::lookup ignores `input`/`position`
 *     entirely (its table is self-contained), and so does PT16RLZParser::
 *     lookupKmerByKey (its short-suffix fallback works from `key` too).
 *
 *   lookup_tail(input, position, key) -> KmerLookupResult
 *     Fewer than KMER_LENGTH characters remain; `key` is the tail packed
 *     from the top bit and padded with zero bits, rolled on from the last
 *     16-mer key by this function (PT16RLZParser::lookupTailByKey,
 *     PT16SassyLookup::lookup_tail(key, length)). found=false, with
 *     match_position/match_length already the final answer.
 *
 * Ordering: unlike the earlier bucketed scans, which only grouped queries
 * by their table BUCKET (the top LOW_BITS of the key -- so a bucket's
 * queries land together, but in no particular order within it), this
 * function fully SORTS every query by its complete 32-bit key before
 * probing. A table bucket's own entries are themselves kept sorted by
 * exactly those low bits (that is what makes std::lower_bound valid
 * inside a bucket once it is large enough -- see BINARY_SEARCH_THRESHOLD),
 * so sorting the queries the same way means the probe order now tracks
 * the table's own layout throughout a bucket, not just up to which bucket
 * it lands in -- and identical keys (a query repeated at another text
 * position) become exactly adjacent, guaranteeing the second one is warm.
 *
 * Done as two stable counting-sort passes (classic LSD radix sort: least
 * significant digit first) instead of a comparison sort, for the same
 * reason the bucket pass always has been one -- linear time, and NUMBER_OF_
 * BUCKETS already happens to be exactly the right size (2^LOW_BITS) for
 * both halves of a 32-bit key: pass 1 sorts by the LOW 16 bits, pass 2
 * stably sorts that result by the HIGH 16 bits, and the combination is a
 * full ascending sort by the whole key.
 *
 * If `timings` is given, it receives the time spent in each phase.
 */
template <typename LookupFn, typename TailFn>
std::vector<KmerLookupResult> sortedKmerScan(
    const std::vector<unsigned char>& input, LookupFn&& lookup,
    TailFn&& lookup_tail, SortedScanTimings* timings = nullptr) {
  using clock = std::chrono::steady_clock;
  const auto ms_since = [](const clock::time_point start) {
    return std::chrono::duration<double, std::milli>(clock::now() - start)
        .count();
  };

  SortedScanTimings local_timings;
  SortedScanTimings& t = timings != nullptr ? *timings : local_timings;

  const std::size_t n = input.size();
  std::vector<KmerLookupResult> results(n);

  const std::size_t kmer_positions =
      n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;

  // ---------- Roll every key: only position 0 pays a full encode_16mer;
  // ---------- every later key is `(previous key << 2) | new character`
  // ---------- (same trick as every other scan in this codebase). ----------

  auto phase_start = clock::now();

  std::vector<std::uint32_t> keys(kmer_positions);
  std::uint32_t key = kmer_positions > 0 ? encode_16mer(input, 0) : 0;

  for (std::size_t i = 0; i < kmer_positions; ++i) {
    if (i > 0) {
      key = (key << 2U) |
            alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
    }
    keys[i] = key;
  }

  phase_barrier(keys.data());
  t.keys_ms = ms_since(phase_start);

  // ---------- Two-pass LSD radix sort by the full 32-bit key ----------

  auto radix_pass = [&](const std::vector<std::uint32_t>& order_in,
                        auto digit_of) {
    auto count =
        std::make_unique<std::array<std::uint32_t, NUMBER_OF_BUCKETS>>();
    count->fill(0);

    for (const std::uint32_t i : order_in) {
      ++(*count)[digit_of(keys[i])];
    }

    std::vector<std::uint32_t> offsets(
        static_cast<std::size_t>(NUMBER_OF_BUCKETS) + 1, 0);
    for (std::uint32_t d = 0; d < NUMBER_OF_BUCKETS; ++d) {
      offsets[d + 1] = offsets[d] + (*count)[d];
    }

    std::vector<std::uint32_t> cursor(offsets.begin(),
                                      offsets.begin() + NUMBER_OF_BUCKETS);
    std::vector<std::uint32_t> order_out(order_in.size());

    for (const std::uint32_t i : order_in) {
      const std::uint32_t d = digit_of(keys[i]);
      order_out[cursor[d]++] = i;
    }

    return order_out;
  };

  phase_start = clock::now();

  std::vector<std::uint32_t> identity(kmer_positions);
  for (std::size_t i = 0; i < kmer_positions; ++i) {
    identity[i] = static_cast<std::uint32_t>(i);
  }

  // Pass 1 (least significant digit): the LOW 16 bits -- the same split a
  // bucket's own internal ordering is sorted by.
  const std::vector<std::uint32_t> by_low = radix_pass(
      identity, [](std::uint32_t k) { return k & (NUMBER_OF_BUCKETS - 1); });

  phase_barrier(by_low.data());
  t.radix_low_ms = ms_since(phase_start);
  phase_start = clock::now();

  // Pass 2 (most significant digit), stable on top of pass 1's result: the
  // HIGH 16 bits -- the same split the table bucket itself is chosen by.
  // Stability is what makes the combination a full sort, not just a
  // by-bucket grouping.
  const std::vector<std::uint32_t> order =
      radix_pass(by_low, [](std::uint32_t k) { return k >> LOW_BITS; });

  phase_barrier(order.data());
  t.radix_high_ms = ms_since(phase_start);

  // ---------- Probe: the lookups themselves, in fully sorted order ----------

  phase_start = clock::now();

  for (const std::uint32_t i : order) {
    results[i] = lookup(input, static_cast<std::size_t>(i), keys[i]);
  }

  phase_barrier(results.data());
  t.probe_ms = ms_since(phase_start);

  // ---------- Tail: not reordered (at most KMER_LENGTH - 1 of them) ----------

  phase_start = clock::now();

  // The last 16-mer key rolled past the end: each `<< 2` drops the first
  // character and pads with a zero one, giving the padded tail key.
  std::uint32_t tail_key = kmer_positions > 0
                               ? keys[kmer_positions - 1] << 2U
                               : (n > 0 ? encode_tail(input, 0) : 0);

  for (std::size_t i = kmer_positions; i < n; ++i, tail_key <<= 2U) {
    results[i] = lookup_tail(input, i, tail_key);

    // A tail position has no 16-mer, so it has no .count/.positions to
    // chain through -- .found must be false regardless of what
    // `lookup_tail` itself set it to: every consumer of KmerLookupResult
    // (chain_extend.hpp, ms_tools.hpp) requires found to mean "a full
    // 16-mer hit, so .count/.positions are populated". Enforced once,
    // here, so a future lookup_tail callable never has to get this right
    // itself.
    results[i].found = false;
  }

  phase_barrier(results.data());
  t.tail_ms = ms_since(phase_start);

  return results;
}

/**
 * The raw-lookup classification (singleton/range/miss) a search over
 * `results` produced, computed directly from the results themselves rather
 * than from a table-format-specific Stats counter -- another consequence
 * of `lookup`/`lookup_tail` being an opaque black box here: this function
 * has no handle on whatever internal stats object the caller's table kept,
 * so it reads the same fact back out of what it already returned instead.
 * Only the first `kmer_positions` entries of `results` are real 16-mer
 * probes (see sortedKmerScan); the tail entries after that are not counted,
 * matching what every earlier scan's own Stats-based composition counted.
 */
inline SearchComposition classifySearchComposition(
    const std::vector<KmerLookupResult>& results,
    const std::size_t kmer_positions) {
  SearchComposition composition;
  composition.available = true;

  for (std::size_t i = 0; i < kmer_positions; ++i) {
    const KmerLookupResult& entry = results[i];

    if (!entry.found) {
      ++composition.misses;
    } else if (entry.count == 1) {
      ++composition.singleton_hits;
    } else {
      ++composition.range_hits;
    }
  }

  return composition;
}
