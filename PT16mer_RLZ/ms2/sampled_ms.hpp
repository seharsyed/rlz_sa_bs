#pragma once

// Experimental, end to end: matching statistics from sampled lookups.
//
//   Phase 1  samples: one 16-mer lookup at every step-th position, keys
//            sorted, finger lookups (sampled.hpp); neighbouring samples
//            linked into chains. Each chain gives a (strict) stretch: from
//            its first singleton sample to its last one. The input matches
//            the reference along one path between them (a range in between
//            takes the path's occurrence -- the linking chooses it); ranges
//            at a chain's ends are left outside, like the holes.
//
//   Phase 2  holes: every position outside the stretches (not a sample,
//            which is already looked up, and not a separator) gets a key,
//            rolled through each hole one character per position; a
//            stretch is skipped and the rolling restarts after it. Keys
//            with a full 16-character window are packed with their index,
//            radix sorted and looked up with a finger; positions with
//            fewer than 16 ACGT characters left get a tail lookup. So
//            every position outside a stretch has a lookup result.
//
//   Phase 3  matching statistics, one pass right to left, with (q, L) the
//            answer at i + 1, from the lookup results only:
//              separator        0;
//              lookup at i      (outside a stretch, and at a stretch's last
//                               singleton): chained when one of its
//                               occurrences is q - 1 -- (q - 1, L + 1),
//                               exact since MS[i] <= MS[i+1] + 1; otherwise
//                               from the result alone: a short phrase
//                               (miss, tail) as is, a singleton extended, a
//                               range narrowed -- exact;
//              stretch interior after its last singleton (whose answer lies
//                               on the path: its only occurrence), every
//                               position to the left up to the first
//                               singleton is the previous answer + 1 at
//                               path - 1 -- exact: the path proves the
//                               character matches, so it reaches the bound.
//                               Filled in one loop, no lookups, no
//                               character read.

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "../rlz_common.hpp"  // rlz::binarySearchLB / RB
#include "sampled.hpp"

namespace ms2 {

struct SampledMsStats {
  std::size_t positions = 0;  // input positions

  // Phase 1
  std::size_t samples = 0;
  std::size_t stretches = 0;          // strict stretches
  std::size_t stretch_positions = 0;  // their positions, both ends included

  // Phase 2 (positions outside the stretches)
  std::size_t hole_keys = 0;        // full windows: sorted lookups
  std::size_t hole_tails = 0;       // fewer than 16 ACGT characters
  std::size_t hole_samples = 0;     // already looked up in phase 1
  std::size_t hole_separators = 0;

  // Phase 3: how every position was answered
  std::size_t chained = 0;       // a lookup's occurrence continuing i + 1
  std::size_t by_short = 0;      // a lookup's short phrase (miss or tail)
  std::size_t by_singleton = 0;  // a lookup's singleton, extended
  std::size_t by_range = 0;      // a lookup's range, narrowed
  std::size_t by_stretch = 0;    // a stretch interior, filled (+1 per step)
  std::size_t separators = 0;
  std::size_t stretch_ends = 0;  // stretches entered (also counted above)
  std::size_t hit_extension = 0;  // characters compared/narrowed at hits

  void add(const SampledMsStats& o) {
    positions += o.positions;
    samples += o.samples;
    stretches += o.stretches;
    stretch_positions += o.stretch_positions;
    hole_keys += o.hole_keys;
    hole_tails += o.hole_tails;
    hole_samples += o.hole_samples;
    hole_separators += o.hole_separators;
    chained += o.chained;
    by_short += o.by_short;
    by_singleton += o.by_singleton;
    by_range += o.by_range;
    by_stretch += o.by_stretch;
    separators += o.separators;
    stretch_ends += o.stretch_ends;
    hit_extension += o.hit_extension;
  }

  std::size_t lookups() const { return samples + hole_keys + hole_tails; }

  Diagnostics diagnostics() const {
    return counter_diagnostics("phase-1", {{"samples", samples},
                                           {"stretches", stretches},
                                           {"stretch-positions",
                                            stretch_positions}}) +
           counter_diagnostics("phase-2", {{"hole-keys", hole_keys},
                                           {"hole-tails", hole_tails},
                                           {"hole-samples", hole_samples},
                                           {"hole-separators",
                                            hole_separators}}) +
           counter_diagnostics("phase-3", {{"chained", chained},
                                           {"by-short", by_short},
                                           {"by-singleton", by_singleton},
                                           {"by-range", by_range},
                                           {"by-stretch", by_stretch},
                                           {"separators", separators},
                                           {"stretch-ends", stretch_ends}}) +
           counter_diagnostics("work", {{"positions", positions},
                                        {"lookups", lookups()},
                                        {"lookups-skipped",
                                         positions - separators - lookups()},
                                        {"hit-extension", hit_extension}});
  }
};

struct SampledMsBuffers {
  std::vector<KmerLookupResult> samples;
  std::vector<SampleKind> kinds;
  std::vector<std::uint32_t> chosen;
  std::vector<std::uint8_t> linked;
  std::vector<std::uint64_t> packed, sorted, scratch;

  // A strict stretch: positions [start, last] (both singleton samples),
  // matching the reference at path + (i - start).
  struct Stretch {
    std::uint32_t start, last;
    std::uint32_t path;
  };
  std::vector<Stretch> stretches;

  std::vector<std::uint32_t> hole_positions;  // text order
  std::vector<KmerLookupResult> hole_results;  // aligned with them
};

// The exact MS at i from a lookup made at i (lazy_ms.hpp's break step).
template <typename Symbol>
std::pair<std::uint32_t, std::uint32_t> resolve_lookup(
    const KmerLookupResult& hit, const std::vector<Symbol>& input,
    const std::vector<Symbol>& reference, const std::size_t i,
    std::size_t& work) {
  if (!hit.found) return {hit.match_position, hit.match_length};

  const std::size_t n = input.size();
  std::size_t offset = KMER_LENGTH;
  std::size_t j = i + KMER_LENGTH;
  std::size_t match = hit.match_position;

  if (hit.count > 1) {
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
      ++work;
    }

    match = hit.positions[static_cast<std::size_t>(nlb)];
    if (nlb != nrb) {
      return {static_cast<std::uint32_t>(match),
              static_cast<std::uint32_t>(offset)};
    }
  }

  while (j < n && match + offset < reference.size() &&
         reference[match + offset] == input[j]) {
    ++j;
    ++offset;
    ++work;
  }

  return {static_cast<std::uint32_t>(match), static_cast<std::uint32_t>(offset)};
}

// The whole method on one table. `FingerPolicy` does the sorted lookups
// (phases 1 and 2), `PlainPolicy` the tails. The phase times go to
// `phases`: 1-collect, 1-sort, 1-lookups, 1-link, 2-collect, 2-sort,
// 2-lookups, 3-ms.
template <typename FingerPolicy, typename PlainPolicy, typename Symbol>
MatchingStatistics sampled_matching_statistics(
    const typename PlainPolicy::Table& table, const std::vector<Symbol>& input,
    const std::vector<Symbol>& reference, const std::size_t step,
    SampledMsBuffers& b, SampledMsStats& stats, Diagnostics* phases) {
  const std::size_t n = input.size();
  MatchingStatistics ms(n);
  stats.positions += n;

  // ---------- Phase 1: samples, sorted lookups, links, stretches ----------

  Diagnostics sample_phases;
  sample_lookups_sorted<FingerPolicy>(table, input, step, b.samples, b.kinds,
                                      b.packed, b.sorted, b.scratch,
                                      &sample_phases);
  stats.samples += b.samples.size();

  const double link_ms = time_ms([&] {
    SampleStats unused;
    link_samples(b.samples, b.kinds, step, b.chosen, b.linked, unused);

    // Strict stretches: each chain from its first singleton to its last.
    b.stretches.clear();
    const std::size_t m = b.samples.size();
    std::size_t j = 0;
    while (j < m) {
      if (!(j + 1 < m && b.linked[j])) {
        ++j;
        continue;
      }
      std::size_t end = j;
      while (end + 1 < m && b.linked[end]) ++end;

      std::size_t first = m, last = m;
      for (std::size_t k = j; k <= end; ++k) {
        if (b.kinds[k] == SampleKind::singleton) {
          if (first == m) first = k;
          last = k;
        }
      }
      if (first != m && first < last) {
        const std::uint32_t start = static_cast<std::uint32_t>(first * step);
        const std::uint32_t stop = static_cast<std::uint32_t>(last * step);
        b.stretches.push_back({start, stop, b.chosen[first]});
        stats.stretch_positions += stop - start + 1;
      }
      j = end + 1;
    }
    stats.stretches += b.stretches.size();
  });

  // ---------- Phase 2: hole keys (rolled), sorted lookups ----------

  const double collect_ms = time_ms([&] {
    b.hole_positions.clear();
    b.hole_results.clear();
    b.packed.clear();

    // A hole position gets a slot: a full window's key is queued for the
    // sorted lookups, a shorter one gets its tail lookup now.
    const auto add_position = [&](std::size_t i, bool full, std::uint32_t key) {
      if (!is_acgt(static_cast<unsigned char>(input[i]))) {
        ++stats.hole_separators;
        return;
      }
      if (i % step == 0) {
        ++stats.hole_samples;
        return;
      }
      const std::size_t k = b.hole_positions.size();
      b.hole_positions.push_back(static_cast<std::uint32_t>(i));
      b.hole_results.emplace_back();

      if (full) {
        b.packed.push_back(static_cast<std::uint64_t>(key) << 32 | k);
        ++stats.hole_keys;
      } else {
        std::uint32_t tail_key;
        const std::uint32_t run = window_at(input, i, tail_key);
        KmerLookupResult result = PlainPolicy::tail(table, tail_key, run);
        result.found = false;
        b.hole_results[k] = result;
        ++stats.hole_tails;
      }
    };

    // Roll through one hole [a, b_end): the key of position i is complete
    // once character i + 15 has been added.
    const auto roll_hole = [&](std::size_t a, std::size_t b_end) {
      std::uint32_t key = 0;
      std::uint32_t run = 0;
      std::size_t next = a;  // next hole position to emit

      for (std::size_t j = a; j < std::min(b_end + KMER_LENGTH - 1, n); ++j) {
        const unsigned char c = static_cast<unsigned char>(input[j]);
        if (is_acgt(c)) {
          key = (key << 2) | alphatab[c];
          if (run < KMER_LENGTH) ++run;
        } else {
          key = 0;
          run = 0;
        }
        if (j + 1 >= a + KMER_LENGTH) {
          const std::size_t i = j + 1 - KMER_LENGTH;
          if (i >= b_end) break;
          add_position(i, run == KMER_LENGTH, key);
          next = i + 1;
        }
      }
      // Positions whose window runs past the input's end.
      for (std::size_t i = next; i < b_end; ++i) add_position(i, false, 0);
    };

    std::size_t at = 0;
    for (const auto& stretch : b.stretches) {
      if (at < stretch.start) roll_hole(at, stretch.start);
      at = stretch.last + 1;
    }
    if (at < n) roll_hole(at, n);
  });

  const double sort_ms =
      time_ms([&] { sorted_order(b.packed, b.sorted, b.scratch); });

  const double lookups_ms = time_ms([&] {
    typename FingerPolicy::State state{};
    for (const std::uint64_t value : b.sorted) {
      b.hole_results[static_cast<std::uint32_t>(value)] =
          FingerPolicy::lookup(table, state, key_of(value));
    }
  });

  // ---------- Phase 3: matching statistics, right to left ----------

  const double ms_ms = time_ms([&] {
    std::size_t hole = b.hole_positions.size();  // cursor, one past
    std::size_t stretch = b.stretches.size();    // cursor, one past

    std::size_t i = n;
    while (i > 0) {
      --i;

      if (!is_acgt(static_cast<unsigned char>(input[i]))) {
        ms[i] = {0, 0};
        ++stats.separators;
        continue;
      }

      std::uint32_t q = 0, L = 0;
      if (i + 1 < n) {
        q = ms[i + 1].first;
        L = ms[i + 1].second;
      }

      // From the lookup result at i: chained to i + 1 if one of its
      // occurrences is q - 1, else from the result alone.
      const auto from_lookup = [&](const KmerLookupResult& result) {
        if (result.found && L > 0 && q > 0) {
          std::size_t scanned = 0;
          if (contains(result, q - 1, scanned)) {
            ms[i] = {q - 1, L + 1};
            ++stats.chained;
            return;
          }
        }
        ms[i] = resolve_lookup(result, input, reference, i, stats.hit_extension);
        if (!result.found) ++stats.by_short;
        else if (result.count == 1) ++stats.by_singleton;
        else ++stats.by_range;
      };

      // A stretch's last singleton: its lookup, then the whole interior.
      if (stretch > 0 && b.stretches[stretch - 1].last == i) {
        const auto& s = b.stretches[--stretch];
        from_lookup(b.samples[i / step]);
        ++stats.stretch_ends;

        // The singleton's answer lies on the path (its only occurrence);
        // each position to the left is one longer, one earlier.
        auto [position, length] = ms[i];
        for (std::size_t k = i; k > s.start;) {
          --k;
          ms[k] = {--position, ++length};
        }
        stats.by_stretch += i - s.start;
        i = s.start;
        continue;
      }

      // Outside the stretches: every position has a lookup result.
      if (i % step == 0) {
        from_lookup(b.samples[i / step]);
        continue;
      }
      while (hole > 0 && b.hole_positions[hole - 1] > i) --hole;
      if (hole > 0 && b.hole_positions[hole - 1] == i) {
        from_lookup(b.hole_results[hole - 1]);
        continue;
      }

      // Unreachable: phase 2 gives every such position a result.
      throw std::logic_error("sampled_ms: position " + std::to_string(i) +
                             " has no lookup result");
    }
  });

  if (phases) {
    Diagnostics::Line line{"phases", true, {}};
    for (const auto& [name, t] : sample_phases.lines.front().values) {
      line.values.emplace_back(std::string("1-") + name, t);
    }
    line.values.emplace_back("1-link", link_ms);
    line.values.emplace_back("2-collect", collect_ms);
    line.values.emplace_back("2-sort", sort_ms);
    line.values.emplace_back("2-lookups", lookups_ms);
    line.values.emplace_back("3-ms", ms_ms);
    *phases = Diagnostics{{line}};
  }

  return ms;
}

}  // namespace ms2
