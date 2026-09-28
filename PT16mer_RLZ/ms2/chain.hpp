#pragma once

// The chain for the ms2 benchmark: the same algorithm as
// backwardChainExtend (../lrf_ms/ms_tools.hpp, see the explanation
// there), copied so it can report progress without touching the original.
//
// `progress(done, total)` is called at every 1/marks of the positions
// (the last call with done == total), so a caller can print a mark with the elapsed time;
// the chain walks the input backward, `done` counts positions handled.
// A long gap between marks points at a slow region (long runs of range
// hits with many occurrences).

#include <cstdint>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"  // KMER_LENGTH, KmerLookupResult, MatchingStatistics

namespace ms2 {

template <typename Progress>
MatchingStatistics chain_extend(const std::vector<KmerLookupResult>& results,
                                Progress&& progress,
                                const std::size_t marks = 20) {
  const std::size_t n = results.size();
  MatchingStatistics ms(n);

  // The next mark is at mark * n / marks positions.
  std::size_t mark = 1;
  std::size_t next_mark = n / marks;

  // (position, length) pairs from the position just handled (one to the
  // right of the current one) that extended past the base 16.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> scratch;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> next_scratch;

  for (std::size_t step = 0; step < n; ++step) {
    while (mark < marks && step == next_mark) {
      progress(step, n);
      ++mark;
      next_mark = mark * n / marks;
    }

    const std::size_t i = n - 1 - step;
    const KmerLookupResult& entry = results[i];

    if (!entry.found) {
      ms[i] = {entry.match_position, entry.match_length};
      scratch.clear();
      continue;
    }

    const bool next_exists = i + 1 < n;
    const KmerLookupResult* next = next_exists ? &results[i + 1] : nullptr;

    // Is `wanted` one of i+1's own occurrences (a hit only)?
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

    // The depth `wanted` is known to extend to, or 0.
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
        return;
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

  progress(n, n);
  return ms;
}

}  // namespace ms2
