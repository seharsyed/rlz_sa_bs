#pragma once

#include <chrono>
#include <cstdint>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"  // KMER_LENGTH, KmerLookupResult, MatchingStatistics

/**
 * backwardChainExtend: multi-step position-only chain extension.
 *
 * chain_extend.hpp's chainExtend() only ever looks one position ahead: a
 * singleton at i extends by exactly one character if i+1's own occurrence
 * set contains position + 1, and stops there even if i+1 ITSELF turned
 * out to extend further right (chainExtend never finds out -- it never
 * looks at i+2). This function chases that all the way: it walks the
 * input BACKWARD, so by the time it reaches position i, it already knows
 * -- for every one of i+1's occurrence positions that itself extended
 * beyond the base 16 -- exactly how deep each one goes. Extending i then
 * means checking i's own candidates against that already-resolved
 * information, not just against i+1's raw occurrence set.
 *
 * The induction (still no reference read anywhere): if input[i..i+16)
 * occurs at reference position p (p is one of i's own occurrences), and
 * input[i+1..i+1+L) == ref[p+1..p+1+L) for some already-established
 * L >= 16 (i.e. position p+1 at i+1 is KNOWN, from having processed i+1
 * already, to extend to length L), then input[i..i+1+L) ==
 * ref[p..p+1+L): p's match at i is at least L+1. The one-step version is
 * the L=16 base case of this; here L can be anything already resolved
 * for p+1, at any earlier (rightward) step.
 *
 * The scratch below carries exactly that: (position, length) pairs for
 * positions at the most-recently-processed index (i+1 relative to the
 * position now being handled) that extended past the base 16 -- i.e.
 * their own depth is already resolved, not just yes/no. A position that
 * stayed at the base 16 is NOT kept in the scratch; it costs nothing to
 * re-check "is p+1 one of i+1's own raw occurrences" directly against
 * results[i+1] when needed, since that's already sitting right there in
 * the array -- no need to carry it forward separately.
 *
 * Why a position can need MORE than one length at once: a range's several
 * occurrence positions are independent -- some of them may inherit a long
 * chain from the scratch, others may not extend at all, and there is no
 * reason they all agree. So every one of a range's own candidates is
 * checked on its own, and the position's reported length is the maximum
 * over all of them. This is deliberately a brute-force, small-scratch
 * approach (linear search within the scratch, no hashing): the bet is
 * that singleton hits dominate in practice (a singleton has exactly one
 * candidate, so no branching at all), keeping both the scratch and the
 * per-position candidate loop small. If ranges turn out to be common and
 * this becomes a real cost, that is the first thing to revisit.
 *
 * Same input contract as chainExtend (chain_extend.hpp): results[i] must
 * be position i's own lookup result, in text order, with a tail position
 * (fewer than 16 characters left) already resolved (found=false,
 * match_length/match_position set to its own answer) by the caller.
 */
inline MatchingStatistics backwardChainExtend(
    const std::vector<KmerLookupResult>& results) {
  const std::size_t n = results.size();
  MatchingStatistics ms(n);

  // (position, length) pairs surviving from the position just processed
  // (one step to the right of whichever position is being handled right
  // now), restricted to those that extended past the base 16. Swapped
  // with next_scratch at the end of every iteration instead of being
  // rebuilt from scratch.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> scratch;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> next_scratch;

  for (std::size_t step = 0; step < n; ++step) {
    const std::size_t i = n - 1 - step;
    const KmerLookupResult& entry = results[i];

    if (!entry.found) {
      // A short factor: already fully resolved, and nothing can chain
      // through a position whose 16-mer does not even occur.
      ms[i] = {entry.match_position, entry.match_length};
      scratch.clear();
      continue;
    }

    const bool next_exists = i + 1 < n;
    const KmerLookupResult* next = next_exists ? &results[i + 1] : nullptr;

    // Is `wanted` one of i+1's own raw occurrences (the plain one-step
    // case, giving length 17)? Only meaningful when i+1 is itself a hit --
    // a miss's match_position is an LCP neighbour, not a genuine
    // occurrence, so it can never be chained into.
    auto extends_into_next_raw = [&](const std::uint32_t wanted) {
      if (!next_exists || !next->found) {
        return false;
      }

      if (next->count == 1) {
        return next->match_position == wanted;
      }

      for (const std::uint32_t candidate : next->positions) {
        if (candidate == wanted) {
          return true;
        }
      }

      return false;
    };

    // Is `wanted` in the scratch (already known to extend past 16, at a
    // specific depth)? Returns that depth, or 0 (never a real length --
    // the shortest possible entry is 17) if not found.
    auto deep_length = [&](const std::uint32_t wanted) -> std::uint32_t {
      for (const auto& [position, length] : scratch) {
        if (position == wanted) {
          return length;
        }
      }

      return 0;
    };

    next_scratch.clear();

    std::uint32_t best_length = KMER_LENGTH;
    std::uint32_t best_position =
        entry.count == 1 ? entry.match_position : entry.positions.front();

    const auto consider = [&](const std::uint32_t p) {
      const std::uint32_t wanted = p + 1;
      std::uint32_t length = 0;

      const std::uint32_t deep = deep_length(wanted);

      if (deep != 0) {
        length = deep + 1;
      } else if (extends_into_next_raw(wanted)) {
        length = KMER_LENGTH + 1;
      }

      if (length == 0) {
        return;  // this candidate does not extend past the base 16
      }

      next_scratch.emplace_back(p, length);

      if (length > best_length) {
        best_length = length;
        best_position = p;
      }
    };

    if (entry.count == 1) {
      consider(entry.match_position);
    } else {
      for (const std::uint32_t p : entry.positions) {
        consider(p);
      }
    }

    ms[i] = {best_position, best_length};
    scratch.swap(next_scratch);
  }

  return ms;
}

// backwardChainExtend, with its wall time written to `chain_ms`, for the
// chain variants' "chain" phase (see Diagnostics::add_phase).
inline MatchingStatistics timedBackwardChainExtend(
    const std::vector<KmerLookupResult>& results, double& chain_ms) {
  const auto start = std::chrono::steady_clock::now();
  MatchingStatistics ms = backwardChainExtend(results);
  phase_barrier(ms.data());
  chain_ms = std::chrono::duration<double, std::milli>(
                 std::chrono::steady_clock::now() - start)
                 .count();
  return ms;
}
