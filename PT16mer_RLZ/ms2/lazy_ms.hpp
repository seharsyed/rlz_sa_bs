#pragma once

// Matching statistics by backward left-extension, with lazy PT16 lookups
// ("option B"): the reference point against lrf-ms.
//
// MS[i] is the longest prefix of input[i..] that occurs in the reference,
// with one reference position where it does. Walking the input backward,
// with (q, L) the answer at i + 1:
//
//   extend   MS[i] <= MS[i+1] + 1 always (a match at i without its first
//            character is a match at i + 1). So if ref[q - 1] == input[i],
//            then input[i .. i + L] occurs at q - 1 and reaches that bound:
//            MS[i] = (q - 1, L + 1), exactly, from one comparison.
//
//   break    otherwise, MS[i] is computed from position i alone, with one
//            table lookup (the only place a lookup happens):
//              - fewer than 16 ACGT characters before a separator or the
//                end: the tail lookup is the answer;
//              - a miss (the 16-mer does not occur): the lookup's closest
//                match is the answer;
//              - a singleton at p: every match of 16 or more starts at p,
//                so extend from p;
//              - a range: narrow among its occurrences (in suffix order)
//                by binary search, one character at a time, as the RLZ
//                parser does, then extend once one occurrence is left.
//
//   separator  a non-ACGT byte never matches: MS[i] = (0, 0).
//
// On an input close to the reference, almost every position extends, and
// lookups happen only at breaks (variants, repeat boundaries).
//
// The table is used through a lookup policy (../lrf_ms/probe_pipeline.hpp:
// V2Policy, SassyPolicy), so any PT16 table variant works.

#include <cstdint>
#include <vector>

#include "../pt16_utils.hpp"  // encode_16mer, alphatab, is_acgt, KMER_LENGTH
#include "../rlz_common.hpp"  // rlz::binarySearchLB / RB

namespace ms2 {

struct LazyStats {
  std::size_t extended = 0;    // answered by one comparison
  std::size_t separators = 0;  // non-ACGT: length 0
  std::size_t breaks = 0;      // needed a lookup (the rest)

  // How the breaks were resolved.
  std::size_t tails = 0;       // fewer than 16 ACGT characters left
  std::size_t misses = 0;      // 16-mer not in the reference
  std::size_t singletons = 0;  // one occurrence, then extended
  std::size_t ranges = 0;      // several occurrences, narrowed

  std::size_t singleton_extension = 0;  // characters compared past 16
  std::size_t range_occurrences = 0;    // occurrences of the ranges met
  std::size_t range_narrowing = 0;      // characters narrowed past 16

  void add(const LazyStats& other) {
    extended += other.extended;
    separators += other.separators;
    breaks += other.breaks;
    tails += other.tails;
    misses += other.misses;
    singletons += other.singletons;
    ranges += other.ranges;
    singleton_extension += other.singleton_extension;
    range_occurrences += other.range_occurrences;
    range_narrowing += other.range_narrowing;
  }

  Diagnostics diagnostics() const {
    return counter_diagnostics("positions", {{"extended", extended},
                                             {"breaks", breaks},
                                             {"separators", separators}}) +
           counter_diagnostics("breaks", {{"tails", tails},
                                          {"misses", misses},
                                          {"singletons", singletons},
                                          {"ranges", ranges}}) +
           counter_diagnostics(
               "work", {{"singleton-extension", singleton_extension},
                        {"range-occurrences", range_occurrences},
                        {"range-narrowing", range_narrowing}});
  }
};

template <typename Policy, typename Symbol>
MatchingStatistics lazy_matching_statistics(
    const typename Policy::Table& table, const std::vector<Symbol>& input,
    const std::vector<Symbol>& reference, LazyStats* stats_out = nullptr) {
  const std::size_t n = input.size();
  MatchingStatistics ms(n);
  LazyStats stats;

  typename Policy::State state{};

  // ACGT characters from i up to the next separator or the end.
  std::size_t run = 0;

  for (std::size_t step = 0; step < n; ++step) {
    const std::size_t i = n - 1 - step;
    const unsigned char c = static_cast<unsigned char>(input[i]);

    if (!is_acgt(c)) {
      ms[i] = {0, 0};
      run = 0;
      ++stats.separators;
      continue;
    }

    ++run;

    // ---------- Extend the match at i + 1 by one character ----------

    if (i + 1 < n) {
      const auto [q, L] = ms[i + 1];

      if (L > 0 && q > 0 && reference[q - 1] == input[i]) {
        ms[i] = {q - 1, L + 1};
        ++stats.extended;
        continue;
      }
    }

    // ---------- Break: compute MS[i] from position i ----------

    ++stats.breaks;

    if (run < KMER_LENGTH) {
      // Pack the run's characters like a key; the rest stays 0.
      std::uint32_t key = 0;
      for (std::size_t j = 0; j < run; ++j) {
        key |= static_cast<std::uint32_t>(
                   alphatab[static_cast<unsigned char>(input[i + j])])
               << (30U - 2U * j);
      }

      const KmerLookupResult tail =
          Policy::tail(table, key, static_cast<std::uint32_t>(run));
      ms[i] = {tail.match_position, tail.match_length};
      ++stats.tails;
      continue;
    }

    const KmerLookupResult hit = Policy::lookup(
        table, state, encode_16mer(input, static_cast<std::uint32_t>(i)));

    if (!hit.found) {
      ms[i] = {hit.match_position, hit.match_length};
      ++stats.misses;
      continue;
    }

    std::size_t offset = KMER_LENGTH;
    std::size_t j = i + KMER_LENGTH;
    std::size_t match = hit.match_position;

    if (hit.count > 1) {
      ++stats.ranges;
      stats.range_occurrences += hit.count;

      std::int64_t nlb = 0;
      std::int64_t nrb = static_cast<std::int64_t>(hit.positions.size()) - 1;

      while (nlb < nrb && j < n) {
        const auto lb = rlz::binarySearchLB(reference, hit.positions, nlb, nrb,
                                            static_cast<std::int64_t>(offset),
                                            input[j]);
        if (!lb) break;

        const auto rb =
            rlz::binarySearchRB(reference, hit.positions, lb.value(), nrb,
                                static_cast<std::int64_t>(offset), input[j]);
        if (!rb) break;

        nlb = lb.value();
        nrb = rb.value();
        ++j;
        ++offset;
        ++stats.range_narrowing;
      }

      match = hit.positions[static_cast<std::size_t>(nlb)];

      // Still several occurrences: the input (or its run) ended first, or
      // the next character matches none of them; either way done.
      if (nlb != nrb) {
        ms[i] = {static_cast<std::uint32_t>(match),
                 static_cast<std::uint32_t>(offset)};
        continue;
      }
    } else {
      ++stats.singletons;
    }

    const std::size_t before = offset;

    while (j < n && match + offset < reference.size() &&
           reference[match + offset] == input[j]) {
      ++j;
      ++offset;
    }

    if (hit.count == 1) {
      stats.singleton_extension += offset - before;
    }

    ms[i] = {static_cast<std::uint32_t>(match),
             static_cast<std::uint32_t>(offset)};
  }

  if (stats_out) *stats_out = stats;
  return ms;
}

}  // namespace ms2
