#pragma once

// Experimental: look up only the input's 16-mers that occur once in the
// input.
//
// After the keys are sorted (input_keys.hpp, sorted_order), equal keys are
// next to each other, so one sequential pass keeps the keys whose
// neighbours differ (keep_unique). The lookups then run over that sorted
// list, into a compact array of results aligned with it: results[k]
// belongs to unique[k] (its position is position_of(unique[k])).

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "../lrf_ms/probe_pipeline.hpp"  // policies, occurrences helpers
#include "input_keys.hpp"

namespace ms2 {

struct UniqueStats {
  std::size_t keys = 0;             // keys before filtering
  std::size_t unique = 0;           // keys that occur once
  std::size_t repeated_keys = 0;    // distinct keys that occur more than once
  std::size_t repeated_positions = 0;  // positions dropped (their copies)
};

// Keeps the keys of the sorted packed list `sorted` (key << 32 | position)
// that occur exactly once, in order, into `unique` (resized to fit).
// Phases: alloc (growing `unique`: only when it is smaller than needed,
// e.g. on the first file), scan.
inline UniqueStats keep_unique(const std::vector<std::uint64_t>& sorted,
                               std::vector<std::uint64_t>& unique,
                               Diagnostics* phases = nullptr) {
  UniqueStats stats;
  const std::size_t n = sorted.size();
  stats.keys = n;

  const double alloc_ms = time_ms([&] { unique.resize(n); });
  std::size_t kept = 0;

  const double scan_ms = time_ms([&] {
  std::size_t i = 0;

  while (i < n) {
    const std::uint32_t key = key_of(sorted[i]);
    std::size_t j = i + 1;
    while (j < n && key_of(sorted[j]) == key) ++j;

    if (j - i == 1) {
      unique[kept++] = sorted[i];
    } else {
      ++stats.repeated_keys;
      stats.repeated_positions += j - i;
    }

    i = j;
  }
  });

  unique.resize(kept);
  stats.unique = kept;

  if (phases) {
    *phases = phase_diagnostics({{"alloc", alloc_ms}, {"scan", scan_ms}});
  }

  return stats;
}

// Looks up every key of `order` with one table, into `results` (resized
// to order.size(); results[k] for order[k]). Phases: alloc (growing
// `results`: only when it is smaller than needed, e.g. its first use),
// lookups (one per key, in `order`).
template <typename Policy>
void probe_compact(const typename Policy::Table& table,
                   const std::vector<std::uint64_t>& order,
                   std::vector<KmerLookupResult>& results,
                   Diagnostics* diagnostics = nullptr) {
  const auto before = table.stats();
  const double alloc_ms = time_ms([&] { results.resize(order.size()); });

  const double keys_ms = time_ms([&] {
    typename Policy::State state{};

    for (std::size_t k = 0; k < order.size(); ++k) {
      results[k] = Policy::lookup(table, state, key_of(order[k]));
    }
  });

  phase_barrier(results.data());

  if (diagnostics) {
    *diagnostics = phase_diagnostics({{"alloc", alloc_ms},
                                      {"lookups", keys_ms}}) +
                   Policy::counters(before, table.stats());
  }
}

// How the lookups resolved.
inline Diagnostics lookup_composition(
    const std::vector<KmerLookupResult>& results) {
  std::size_t singletons = 0, ranges = 0, misses = 0;
  for (const KmerLookupResult& result : results) {
    if (!result.found) ++misses;
    else if (result.count == 1) ++singletons;
    else ++ranges;
  }
  return counter_diagnostics(
      "lookups", {{"singletons", singletons}, {"ranges", ranges},
                  {"misses", misses}});
}

// `results` against `expected` (both aligned with `order`): same found,
// length, count and occurrence set; and every reported position really
// matches the input for the reported length. Returns the number of
// differing entries; `first` / `reason` describe the first one.
template <typename Symbol>
std::size_t compare_compact(const std::vector<KmerLookupResult>& expected,
                            const std::vector<KmerLookupResult>& results,
                            const std::vector<std::uint64_t>& order,
                            const std::vector<Symbol>& input,
                            const std::vector<Symbol>& reference,
                            std::size_t& first, std::string& reason) {
  std::size_t mismatches = 0;
  const auto mismatch = [&](std::size_t k, const char* why) {
    if (mismatches++ == 0) {
      first = k;
      reason = why;
    }
  };

  if (expected.size() != results.size()) {
    mismatch(0, "different number of results");
    return mismatches;
  }

  for (std::size_t k = 0; k < results.size(); ++k) {
    const KmerLookupResult& want = expected[k];
    const KmerLookupResult& got = results[k];

    if (got.found != want.found || got.match_length != want.match_length) {
      mismatch(k, "found/match_length");
      continue;
    }

    if (got.found &&
        (got.count != want.count ||
         !msbench::same_occurrences(msbench::occurrences(got),
                                    msbench::occurrences(want)))) {
      mismatch(k, "occurrences");
      continue;
    }

    const std::size_t i = position_of(order[k]);
    const std::size_t p = got.match_position;
    const std::size_t len = got.match_length;

    if (p + len > reference.size() || i + len > input.size() ||
        !std::equal(input.begin() + static_cast<std::ptrdiff_t>(i),
                    input.begin() + static_cast<std::ptrdiff_t>(i + len),
                    reference.begin() + static_cast<std::ptrdiff_t>(p))) {
      mismatch(k, "match_position is not a real match");
    }
  }

  return mismatches;
}

}  // namespace ms2
