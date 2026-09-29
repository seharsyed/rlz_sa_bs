#pragma once

// Experimental: chaining the lookups of the input's unique 16-mers
// (unique_keys.hpp) in one right-to-left sweep.
//
// text_order puts the unique keys back in input text order, without an
// input-sized array: each unique key's (position, index k into the compact
// results) is packed as position << 32 | k and radix sorted.
//
// chain_unique then walks those positions right to left, keeping in a
// scratchbook the last entry seen (its input position P, the reference
// position q of its match, its length L, and whether L is exact). For the
// entry at input position i, with d = P - i:
//
//   miss        the lookup's answer (a short factor, under 16) is exact;
//   found, d <= 16 and q - d is one of the entry's occurrences
//               the entry's 16-mer and the previous match line up (they
//               overlap or touch), so its match is (q - d, L + d). If the
//               previous entry is exact, so is this one -- MS[i] <=
//               MS[i + d] + d, and this reaches it; otherwise it is an
//               approximate lower bound;
//   found, d > 16
//               a gap: unresolved, lower bound 16, approximate;
//   found, d <= 16 but no occurrence lines up
//               unresolved, lower bound 16, approximate.
//
// Only the previous entry's chosen position q is checked (not all of a
// range's occurrences), so each entry costs at most one scan of its own
// occurrences.

#include <cstdint>
#include <vector>

#include "input_keys.hpp"  // sorted_order, key_of / position_of

namespace ms2 {

// Puts the unique keys in input text order: `by_position` holds
// position << 32 | k, sorted by position, where k indexes `unique` (and the
// results aligned with it). Phases: pack, then sorted_order's own.
inline void text_order(const std::vector<std::uint64_t>& unique,
                       std::vector<std::uint64_t>& by_position,
                       std::vector<std::uint64_t>& packed,
                       std::vector<std::uint64_t>& scratch,
                       Diagnostics* phases = nullptr) {
  const double pack_ms = time_ms([&] {
    packed.resize(unique.size());
    for (std::size_t k = 0; k < unique.size(); ++k) {
      packed[k] = static_cast<std::uint64_t>(position_of(unique[k])) << 32 | k;
    }
  });

  Diagnostics sort_phases;
  sorted_order(packed, by_position, scratch, &sort_phases);

  if (phases) {
    *phases = phase_diagnostics({{"pack", pack_ms}});
    for (const auto& [name, ms] : sort_phases.lines.front().values) {
      phases->lines.front().values.emplace_back(name, ms);
    }
  }
}

enum class ChainStatus : std::uint8_t {
  exact_miss,       // a short factor: the lookup's own answer
  exact_chained,    // lined up with an exact entry to its right
  approx_chained,   // lined up with an approximate entry: a lower bound
  approx_gap,       // the next entry is more than 16 positions away: 16
  approx_no_link,   // within 16, but no occurrence lines up: 16
};

inline bool is_exact(const ChainStatus status) {
  return status == ChainStatus::exact_miss ||
         status == ChainStatus::exact_chained;
}

// One unique position's chained answer, aligned with `by_position`.
struct ChainEntry {
  std::uint32_t reference_position;
  std::uint32_t length;  // exact, or a lower bound when approximate
  ChainStatus status;
};

struct UniqueChainStats {
  std::size_t exact_miss = 0;
  std::size_t exact_chained = 0;
  std::size_t approx_chained = 0;
  std::size_t approx_gap = 0;
  std::size_t approx_no_link = 0;
  std::size_t occurrences_scanned = 0;  // range occurrences checked

  std::size_t exact() const { return exact_miss + exact_chained; }
  std::size_t approximate() const {
    return approx_chained + approx_gap + approx_no_link;
  }

  void add(const UniqueChainStats& other) {
    exact_miss += other.exact_miss;
    exact_chained += other.exact_chained;
    approx_chained += other.approx_chained;
    approx_gap += other.approx_gap;
    approx_no_link += other.approx_no_link;
    occurrences_scanned += other.occurrences_scanned;
  }

  Diagnostics diagnostics() const {
    return counter_diagnostics("exact", {{"misses", exact_miss},
                                         {"chained", exact_chained}}) +
           counter_diagnostics("approximate",
                               {{"chained", approx_chained},
                                {"gap", approx_gap},
                                {"no-link", approx_no_link}}) +
           counter_diagnostics("work",
                               {{"range-occurrences-scanned",
                                 occurrences_scanned}});
  }
};

// Is `wanted` one of `result`'s occurrences? A range is scanned (its
// occurrences are in suffix order, not by position).
inline bool has_occurrence(const KmerLookupResult& result,
                           const std::uint32_t wanted,
                           std::size_t& scanned) {
  if (result.count == 1) return result.match_position == wanted;

  scanned += result.positions.size();
  for (const std::uint32_t position : result.positions) {
    if (position == wanted) return true;
  }
  return false;
}

// The sweep. `results[k]` is the lookup of unique key k; `by_position`
// comes from text_order. Fills `chained` (aligned with by_position).
inline UniqueChainStats chain_unique(
    const std::vector<KmerLookupResult>& results,
    const std::vector<std::uint64_t>& by_position,
    std::vector<ChainEntry>& chained) {
  UniqueChainStats stats;
  const std::size_t u = by_position.size();
  chained.resize(u);

  // The scratchbook: the last entry seen (to the right).
  bool have_previous = false;
  std::uint32_t previous_position = 0;  // P, in the input
  std::uint32_t previous_match = 0;     // q, in the reference
  std::uint32_t previous_length = 0;    // L
  bool previous_exact = false;

  for (std::size_t step = 0; step < u; ++step) {
    const std::size_t idx = u - 1 - step;
    const std::uint32_t i = key_of(by_position[idx]);  // the position half
    const KmerLookupResult& result =
        results[static_cast<std::uint32_t>(by_position[idx])];

    ChainEntry entry;

    if (!result.found) {
      entry = {result.match_position, result.match_length,
               ChainStatus::exact_miss};
      ++stats.exact_miss;
    } else if (have_previous && previous_position - i <= KMER_LENGTH) {
      const std::uint32_t d = previous_position - i;

      if (previous_match >= d &&
          has_occurrence(result, previous_match - d,
                         stats.occurrences_scanned)) {
        const ChainStatus status = previous_exact ? ChainStatus::exact_chained
                                                  : ChainStatus::approx_chained;
        entry = {previous_match - d, previous_length + d, status};
        ++(previous_exact ? stats.exact_chained : stats.approx_chained);
      } else {
        entry = {result.match_position, KMER_LENGTH,
                 ChainStatus::approx_no_link};
        ++stats.approx_no_link;
      }
    } else {
      entry = {result.match_position, KMER_LENGTH, ChainStatus::approx_gap};
      ++stats.approx_gap;
    }

    chained[idx] = entry;

    have_previous = true;
    previous_position = i;
    previous_match = entry.reference_position;
    previous_length = entry.length;
    previous_exact = is_exact(entry.status);
  }

  return stats;
}

}  // namespace ms2
