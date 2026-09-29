#pragma once

// Experimental: matching statistics from sampled lookups -- first steps.
//
// The idea is to avoid redundant lookups: on an input close to the
// reference, a 16-mer at every position is mostly wasted work, because
// consecutive 16-mers sit at consecutive reference positions. So only
// every l-th position is looked up, and lookups are added only where the
// samples do not line up.
//
//   1. sample_lookups: in text order, one lookup at every l-th position
//      s_j = j * l, into an array of about n / l results: a 16-mer lookup
//      where the window is all ACGT; a tail lookup (a short phrase,
//      exact) where fewer than 16 ACGT characters are left before a
//      separator or the end; an empty entry at a separator.
//
//   2. link_samples: right to left over the samples. A hit at s_j links to
//      the sample at s_{j+1} = s_j + l when one of its occurrences p has
//      p + l == q_{j+1} (the chosen reference position of sample j + 1):
//      then input[s_j .. s_j + l + 16) matches ref[p ..] contiguously
//      (for l <= 16 the two 16-mers overlap or touch). A range keeps only
//      the occurrence that links; a range that does not link is left
//      undecided, and the sample to its left may still link into it (then
//      that decides it). Misses and tails (short phrases) never link.
//      Maximal runs of linked samples are stretches: one reference
//      alignment, no lookups needed inside.
//
//   3. midpoint_lookups: for every pair of neighbouring samples that did
//      not link (a hole), one 16-mer lookup at s_j + l / 2, recording
//      whether it lines up with the left sample, the right one, both or
//      neither. Nothing is inferred from it yet.

#include <cstdint>
#include <limits>
#include <vector>

#include "../pt16_utils.hpp"  // alphatab, is_acgt, KMER_LENGTH, Diagnostics
#include "input_keys.hpp"     // sorted_order, key_of

namespace ms2 {

inline constexpr std::uint32_t UNDECIDED =
    std::numeric_limits<std::uint32_t>::max();

enum class SampleKind : std::uint8_t {
  separator,  // a non-ACGT byte: match length 0
  tail,       // fewer than 16 ACGT characters left: a short phrase (exact)
  miss,       // 16-mer not in the reference: a short phrase (exact)
  singleton,
  range,
};

struct SampleStats {
  std::size_t samples = 0;
  std::size_t separators = 0, tails = 0, misses = 0, singletons = 0,
              ranges = 0;

  std::size_t pairs = 0;         // neighbouring sample pairs
  std::size_t linked = 0;        // of them, lined up
  std::size_t range_decided_by_left = 0;  // an undecided range fixed later
  std::size_t occurrences_scanned = 0;

  std::size_t stretches = 0;           // maximal runs of linked samples
  std::size_t stretch_characters = 0;  // input characters they cover

  // Strict stretches (strict_stretches): a chain of linked samples cut to
  // run from its first singleton to its last one (two different samples).
  // Ranges strictly inside are dismissed (the singleton path fixes their
  // occurrence); ranges outside -- at a chain's ends, or in a chain with
  // fewer than two singletons -- are kept for later evaluation.
  std::size_t strict_stretches = 0;
  std::size_t strict_characters = 0;
  std::size_t ranges_dismissed = 0;
  std::size_t ranges_kept = 0;
  std::size_t chains_without_singleton = 0;
  std::size_t chains_one_singleton = 0;

  // Overlaps (overlap_stats): a chain continued past either end along its
  // own alignment, through samples whose occurrences allow it (only
  // ranges can: they may fit two alignments).
  std::size_t chain_ends = 0;          // ends tried (two per chain)
  std::size_t overlapping_ends = 0;    // ends that reach into another chain
  std::size_t overlap_samples = 0;     // distinct samples of another chain
  std::size_t overlap_positions = 0;   // the input positions they cover
  std::size_t hole_extension_samples = 0;  // distinct samples outside chains
  std::size_t overlap_into_singleton = 0;  // should be 0

  // Every input position is in a stretch, a separator (non-ACGT byte), or
  // a hole (everything else: what still has to be resolved).
  std::size_t positions = 0;
  std::size_t separator_positions = 0;
  std::size_t hole_positions() const {
    return positions - stretch_characters - separator_positions;
  }

  // Midpoint lookups, one per hole (holes between two separator samples,
  // inside a run of non-ACGT bytes, are only counted).
  std::size_t separator_holes = 0;
  std::size_t holes = 0;
  std::size_t mid_no_window = 0;  // no full ACGT window at the midpoint
  std::size_t mid_miss = 0;
  std::size_t mid_left = 0, mid_right = 0, mid_both = 0, mid_neither = 0;

  void add(const SampleStats& o) {
    samples += o.samples;
    separators += o.separators;
    tails += o.tails;
    misses += o.misses;
    singletons += o.singletons;
    ranges += o.ranges;
    pairs += o.pairs;
    linked += o.linked;
    range_decided_by_left += o.range_decided_by_left;
    occurrences_scanned += o.occurrences_scanned;
    stretches += o.stretches;
    stretch_characters += o.stretch_characters;
    strict_stretches += o.strict_stretches;
    strict_characters += o.strict_characters;
    ranges_dismissed += o.ranges_dismissed;
    ranges_kept += o.ranges_kept;
    chains_without_singleton += o.chains_without_singleton;
    chains_one_singleton += o.chains_one_singleton;
    chain_ends += o.chain_ends;
    overlapping_ends += o.overlapping_ends;
    overlap_samples += o.overlap_samples;
    overlap_positions += o.overlap_positions;
    hole_extension_samples += o.hole_extension_samples;
    overlap_into_singleton += o.overlap_into_singleton;
    positions += o.positions;
    separator_positions += o.separator_positions;
    separator_holes += o.separator_holes;
    holes += o.holes;
    mid_no_window += o.mid_no_window;
    mid_miss += o.mid_miss;
    mid_left += o.mid_left;
    mid_right += o.mid_right;
    mid_both += o.mid_both;
    mid_neither += o.mid_neither;
  }

  Diagnostics diagnostics() const {
    return counter_diagnostics("samples", {{"total", samples},
                                           {"singletons", singletons},
                                           {"ranges", ranges},
                                           {"misses", misses},
                                           {"tails", tails},
                                           {"separators", separators}}) +
           counter_diagnostics("links", {{"pairs", pairs},
                                         {"linked", linked},
                                         {"range-decided-by-left",
                                          range_decided_by_left},
                                         {"occurrences-scanned",
                                          occurrences_scanned}}) +
           counter_diagnostics("stretches",
                               {{"count", stretches},
                                {"characters", stretch_characters}}) +
           counter_diagnostics("strict",
                               {{"stretches", strict_stretches},
                                {"characters", strict_characters},
                                {"ranges-dismissed", ranges_dismissed},
                                {"ranges-kept", ranges_kept},
                                {"chains-without-singleton",
                                 chains_without_singleton},
                                {"chains-one-singleton",
                                 chains_one_singleton}}) +
           counter_diagnostics("overlaps",
                               {{"chain-ends", chain_ends},
                                {"overlapping-ends", overlapping_ends},
                                {"samples", overlap_samples},
                                {"positions", overlap_positions},
                                {"into-singleton", overlap_into_singleton},
                                {"hole-extension-samples",
                                 hole_extension_samples}}) +
           counter_diagnostics("positions",
                               {{"in-stretches", stretch_characters},
                                {"in-holes", hole_positions()},
                                {"separators", separator_positions}}) +
           counter_diagnostics("midpoints", {{"holes", holes},
                                             {"separator-holes",
                                              separator_holes},
                                             {"left", mid_left},
                                             {"right", mid_right},
                                             {"both", mid_both},
                                             {"neither", mid_neither},
                                             {"miss", mid_miss},
                                             {"no-window", mid_no_window}});
  }
};

// The ACGT run starting at `s` (capped at 16), and its characters packed
// like a key (unused bits 0).
template <typename Symbol>
std::uint32_t window_at(const std::vector<Symbol>& input, const std::size_t s,
                        std::uint32_t& key) {
  key = 0;
  std::uint32_t run = 0;
  while (run < KMER_LENGTH && s + run < input.size() &&
         is_acgt(static_cast<unsigned char>(input[s + run]))) {
    key |= static_cast<std::uint32_t>(
               alphatab[static_cast<unsigned char>(input[s + run])])
           << (30U - 2U * run);
    ++run;
  }
  return run;
}

inline SampleKind kind_of(const KmerLookupResult& result,
                          const std::uint32_t run) {
  if (run == 0) return SampleKind::separator;
  if (run < KMER_LENGTH) return SampleKind::tail;
  if (!result.found) return SampleKind::miss;
  return result.count == 1 ? SampleKind::singleton : SampleKind::range;
}

// Step 1. `samples[j]` is the lookup at j * step; `kinds[j]` what it was.
template <typename Policy, typename Symbol>
void sample_lookups(const typename Policy::Table& table,
                    const std::vector<Symbol>& input, const std::size_t step,
                    std::vector<KmerLookupResult>& samples,
                    std::vector<SampleKind>& kinds) {
  const std::size_t m = (input.size() + step - 1) / step;
  samples.resize(m);
  kinds.resize(m);

  typename Policy::State state{};

  for (std::size_t j = 0; j < m; ++j) {
    std::uint32_t key;
    const std::uint32_t run = window_at(input, j * step, key);

    KmerLookupResult result;
    if (run == KMER_LENGTH) {
      result = Policy::lookup(table, state, key);
    } else if (run > 0) {
      result = Policy::tail(table, key, run);
      result.found = false;
    }  // separator: empty

    samples[j] = result;
    kinds[j] = kind_of(result, run);
  }
}

inline bool contains(const KmerLookupResult& result, const std::uint32_t wanted,
                     std::size_t& scanned) {
  if (result.count == 1) return result.match_position == wanted;
  scanned += result.positions.size();
  for (const std::uint32_t position : result.positions) {
    if (position == wanted) return true;
  }
  return false;
}

// Step 1, sorted: the same samples, but the 16-mer lookups are done in
// sorted key order (so a finger search walks through the table): the keys
// of the full windows are packed as key << 32 | j, radix sorted, looked up
// in that order, and each result written to samples[j]. Tails and
// separators are answered directly. Phases: collect (keys, tails), sort,
// lookups.
template <typename Policy, typename Symbol>
void sample_lookups_sorted(const typename Policy::Table& table,
                           const std::vector<Symbol>& input,
                           const std::size_t step,
                           std::vector<KmerLookupResult>& samples,
                           std::vector<SampleKind>& kinds,
                           std::vector<std::uint64_t>& packed,
                           std::vector<std::uint64_t>& sorted,
                           std::vector<std::uint64_t>& scratch,
                           Diagnostics* phases = nullptr) {
  const std::size_t m = (input.size() + step - 1) / step;

  const double collect_ms = time_ms([&] {
    samples.resize(m);
    kinds.resize(m);
    packed.clear();

    for (std::size_t j = 0; j < m; ++j) {
      std::uint32_t key;
      const std::uint32_t run = window_at(input, j * step, key);

      if (run == KMER_LENGTH) {
        packed.push_back(static_cast<std::uint64_t>(key) << 32 | j);
        continue;  // looked up below
      }

      KmerLookupResult result;
      if (run > 0) {
        result = Policy::tail(table, key, run);
        result.found = false;
      }
      samples[j] = result;
      kinds[j] = kind_of(result, run);
    }
  });

  const double sort_ms =
      time_ms([&] { sorted_order(packed, sorted, scratch); });

  const double lookups_ms = time_ms([&] {
    typename Policy::State state{};
    for (const std::uint64_t value : sorted) {
      const std::uint32_t j = static_cast<std::uint32_t>(value);
      samples[j] = Policy::lookup(table, state, key_of(value));
      kinds[j] = kind_of(samples[j], KMER_LENGTH);
    }
  });

  if (phases) {
    *phases = phase_diagnostics(
        {{"collect", collect_ms}, {"sort", sort_ms}, {"lookups", lookups_ms}});
  }
}

// Step 2. `chosen[j]`: the reference position of hit j's alignment
// (UNDECIDED for a range nothing linked to, and for non-hits);
// `linked[j]`: sample j lines up with sample j + 1.
inline void link_samples(const std::vector<KmerLookupResult>& samples,
                         const std::vector<SampleKind>& kinds,
                         const std::size_t step,
                         std::vector<std::uint32_t>& chosen,
                         std::vector<std::uint8_t>& linked,
                         SampleStats& stats) {
  const std::size_t m = samples.size();
  chosen.assign(m, UNDECIDED);
  linked.assign(m, 0);

  const auto is_hit = [&](std::size_t j) {
    return kinds[j] == SampleKind::singleton || kinds[j] == SampleKind::range;
  };

  for (std::size_t back = 0; back < m; ++back) {
    const std::size_t j = m - 1 - back;
    const KmerLookupResult& result = samples[j];

    if (!is_hit(j)) continue;

    if (kinds[j] == SampleKind::singleton) chosen[j] = result.match_position;

    if (j + 1 >= m || !is_hit(j + 1)) continue;

    const std::uint32_t d = static_cast<std::uint32_t>(step);

    if (chosen[j + 1] != UNDECIDED) {
      // The right sample's alignment is known: does ours continue it?
      const std::uint32_t q = chosen[j + 1];
      if (q >= d && contains(result, q - d, stats.occurrences_scanned)) {
        chosen[j] = q - d;
        linked[j] = 1;
      }
    } else if (kinds[j] == SampleKind::singleton) {
      // The right sample is an undecided range: does it contain ours + d?
      if (contains(samples[j + 1], chosen[j] + d, stats.occurrences_scanned)) {
        chosen[j + 1] = chosen[j] + d;
        linked[j] = 1;
        ++stats.range_decided_by_left;
      }
    }
    // Two undecided ranges next to each other: not linked (rare).
  }
}

// Stretches: maximal runs of linked samples, and the characters they cover
// (from the first sample to the end of the last one's 16-mer).
inline void count_stretches(const std::vector<SampleKind>& kinds,
                            const std::vector<std::uint8_t>& linked,
                            const std::size_t step, const std::size_t n,
                            SampleStats& stats) {
  const std::size_t m = kinds.size();
  std::size_t j = 0;
  while (j < m) {
    if (j + 1 < m && linked[j]) {
      std::size_t end = j;
      while (end + 1 < m && linked[end]) ++end;
      ++stats.stretches;
      stats.stretch_characters +=
          std::min(end * step + KMER_LENGTH, n) - j * step;
      j = end + 1;
    } else {
      ++j;
    }
  }
}

// Strict stretches, from the same links: each chain (maximal run of
// linked samples, at least two) is cut to run from its first singleton to
// its last. See SampleStats.
inline void strict_stretches(const std::vector<SampleKind>& kinds,
                             const std::vector<std::uint8_t>& linked,
                             const std::size_t step, const std::size_t n,
                             SampleStats& stats) {
  const std::size_t m = kinds.size();
  std::size_t j = 0;

  while (j < m) {
    if (!(j + 1 < m && linked[j])) {
      ++j;
      continue;
    }

    std::size_t end = j;
    while (end + 1 < m && linked[end]) ++end;

    // First and last singleton of the chain j..end.
    std::size_t first = m, last = m, ranges = 0;
    for (std::size_t k = j; k <= end; ++k) {
      if (kinds[k] == SampleKind::singleton) {
        if (first == m) first = k;
        last = k;
      } else if (kinds[k] == SampleKind::range) {
        ++ranges;
      }
    }

    if (first == m) {
      ++stats.chains_without_singleton;
      stats.ranges_kept += ranges;
    } else if (first == last) {
      ++stats.chains_one_singleton;
      stats.ranges_kept += ranges;
    } else {
      ++stats.strict_stretches;
      stats.strict_characters +=
          std::min(last * step + KMER_LENGTH, n) - first * step;
      for (std::size_t k = j; k <= end; ++k) {
        if (kinds[k] != SampleKind::range) continue;
        if (k > first && k < last) ++stats.ranges_dismissed;
        else ++stats.ranges_kept;
      }
    }

    j = end + 1;
  }
}

// Overlapping stretches: every chain is continued past each end along its
// own alignment -- to the right while the next sample contains the
// alignment's position + step, to the left while the previous one contains
// position - step. Samples of another chain reached this way overlap (only
// a range can be reached: a singleton reached would mean the two chains
// were one alignment and should have linked); samples outside any chain
// reached this way extend the chain into a hole. See SampleStats.
inline void overlap_stats(const std::vector<KmerLookupResult>& samples,
                          const std::vector<SampleKind>& kinds,
                          const std::vector<std::uint32_t>& chosen,
                          const std::vector<std::uint8_t>& linked,
                          const std::size_t step, const std::size_t n,
                          SampleStats& stats) {
  const std::size_t m = kinds.size();
  const auto is_hit = [&](std::size_t k) {
    return kinds[k] == SampleKind::singleton || kinds[k] == SampleKind::range;
  };
  const auto in_chain = [&](std::size_t k) {
    return (k > 0 && linked[k - 1]) || (k + 1 < m && linked[k]);
  };

  // 1: reached as an overlap, 2: reached as a hole extension.
  std::vector<std::uint8_t> reached(m, 0);
  const std::uint32_t d = static_cast<std::uint32_t>(step);
  std::size_t scanned = 0;

  const auto reach = [&](std::size_t k, bool& overlapped) {
    if (in_chain(k)) {
      if (kinds[k] == SampleKind::singleton) ++stats.overlap_into_singleton;
      reached[k] |= 1;
      overlapped = true;
    } else {
      reached[k] |= 2;
    }
  };

  for (std::size_t k = 0; k < m; ++k) {
    if (!in_chain(k) || chosen[k] == UNDECIDED) continue;

    const bool is_end = !(k + 1 < m && linked[k]);
    const bool is_start = !(k > 0 && linked[k - 1]);

    if (is_end) {  // continue to the right
      ++stats.chain_ends;
      bool overlapped = false;
      std::uint32_t position = chosen[k];
      for (std::size_t next = k + 1;
           next < m && is_hit(next) &&
           contains(samples[next], position + d, scanned);
           ++next) {
        position += d;
        reach(next, overlapped);
      }
      stats.overlapping_ends += overlapped ? 1 : 0;
    }

    if (is_start) {  // continue to the left
      ++stats.chain_ends;
      bool overlapped = false;
      std::uint32_t position = chosen[k];
      for (std::size_t previous = k;
           previous > 0 && is_hit(previous - 1) && position >= d &&
           contains(samples[previous - 1], position - d, scanned);
           --previous) {
        position -= d;
        reach(previous - 1, overlapped);
      }
      stats.overlapping_ends += overlapped ? 1 : 0;
    }
  }

  for (std::size_t k = 0; k < m; ++k) {
    if (reached[k] & 1) {
      ++stats.overlap_samples;
      stats.overlap_positions += std::min(k * step + step, n) - k * step;
    }
    if (reached[k] & 2) ++stats.hole_extension_samples;
  }
  stats.occurrences_scanned += scanned;
}

// Step 3: one lookup at the midpoint of every hole.
template <typename Policy, typename Symbol>
void midpoint_lookups(const typename Policy::Table& table,
                      const std::vector<Symbol>& input,
                      const std::vector<KmerLookupResult>& samples,
                      const std::vector<SampleKind>& kinds,
                      const std::vector<std::uint32_t>& chosen,
                      const std::vector<std::uint8_t>& linked,
                      const std::size_t step, SampleStats& stats) {
  const std::size_t m = samples.size();
  const std::uint32_t half = static_cast<std::uint32_t>(step / 2);
  typename Policy::State state{};

  for (std::size_t j = 0; j + 1 < m; ++j) {
    if (linked[j]) continue;

    if (kinds[j] == SampleKind::separator &&
        kinds[j + 1] == SampleKind::separator) {
      ++stats.separator_holes;
      continue;
    }
    ++stats.holes;

    std::uint32_t key;
    const std::size_t mid = j * step + half;
    if (window_at(input, mid, key) < KMER_LENGTH) {
      ++stats.mid_no_window;
      continue;
    }

    const KmerLookupResult result = Policy::lookup(table, state, key);
    if (!result.found) {
      ++stats.mid_miss;
      continue;
    }

    const bool left = chosen[j] != UNDECIDED &&
                      contains(result, chosen[j] + half,
                               stats.occurrences_scanned);
    const bool right = chosen[j + 1] != UNDECIDED && chosen[j + 1] >= half &&
                       contains(result, chosen[j + 1] - half,
                                stats.occurrences_scanned);

    if (left && right) ++stats.mid_both;
    else if (left) ++stats.mid_left;
    else if (right) ++stats.mid_right;
    else ++stats.mid_neither;
  }
}

// Untimed checks: every stretch's alignment really matches the reference
// over its span, and lrf-ms's match at the stretch start is at least that
// long. Returns the number of stretches failing either.
template <typename Symbol>
std::size_t check_stretches(const std::vector<SampleKind>& kinds,
                            const std::vector<std::uint32_t>& chosen,
                            const std::vector<std::uint8_t>& linked,
                            const std::size_t step,
                            const std::vector<Symbol>& input,
                            const std::vector<Symbol>& reference,
                            const MatchingStatistics& want) {
  const std::size_t m = kinds.size();
  const std::size_t n = input.size();
  std::size_t bad = 0;
  std::size_t j = 0;

  while (j < m) {
    if (!(j + 1 < m && linked[j])) {
      ++j;
      continue;
    }

    std::size_t end = j;
    while (end + 1 < m && linked[end]) ++end;

    const std::size_t start = j * step;
    const std::size_t stop = std::min(end * step + KMER_LENGTH, n);
    const std::size_t p = chosen[j];
    bool ok = p != UNDECIDED && p + (stop - start) <= reference.size();

    for (std::size_t k = j; ok && k <= end; ++k) {
      ok = chosen[k] == p + (k - j) * step;
    }
    for (std::size_t i = start; ok && i < stop; ++i) {
      ok = input[i] == reference[p + (i - start)];
    }
    ok = ok && want[start].second >= stop - start;

    bad += ok ? 0 : 1;
    j = end + 1;
  }

  return bad;
}

}  // namespace ms2
