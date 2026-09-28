#pragma once

// 16-mer table lookups over an input's keys in a probe order (bucket or
// sorted, input_keys.hpp), for the ms2 benchmark. The table side is reused as is: the lookup policies
// (V2Policy, SassyPolicy, SassyFingerPolicy) and compare_lookups from
// ../lrf_ms/probe_pipeline.hpp.
//
// The result of every lookup is stored at its own text position, in one
// array of KmerLookupResult per input:
//
//   key positions     the 16-mer lookup, in the given order (in sorted
//                     order a finger search can walk forward through the
//                     table);
//   short queries     the tail lookup of their 1 to 15 characters, with
//                     found = false (the chain treats them as resolved);
//   separators        left empty (found = false, length 0): the correct
//                     answer, no search;
//   run positions     left empty for now: all-A / all-T 16-mers are not
//                     looked up yet (see input_keys.hpp).
//
// Every lookup of a file writes the same positions (the keys and short
// queries), whatever the order and the table. So a results array is
// prepared once per file (prepare_results: sized to the input, separator
// and run positions emptied), and then each probe only writes its
// lookups.

#include <algorithm>
#include <cstdint>
#include <vector>

#include "../lrf_ms/probe_pipeline.hpp"  // policies, compare_lookups
#include "input_keys.hpp"

namespace ms2 {

// Sizes `results` to the input and empties the positions no lookup writes:
// the separators (non-ACGT bytes of `input`) and the run positions. What
// the array held before (an earlier file's results) does not matter:
// every other position is overwritten by the file's lookups. Phases:
// resize (only when the size changes; new entries start empty), clear.
template <typename Symbol>
void prepare_results(const std::vector<Symbol>& input, const KeyedInput& keyed,
                     std::vector<KmerLookupResult>& results,
                     Diagnostics* phases = nullptr) {
  const double resize_ms =
      time_ms([&] { results.resize(keyed.input_size); });

  const double clear_ms = time_ms([&] {
    for (std::size_t i = 0; i < input.size(); ++i) {
      if (!is_acgt(static_cast<unsigned char>(input[i]))) {
        results[i] = KmerLookupResult{};
      }
    }

    for (const HomopolymerRun& run : keyed.runs) {
      for (std::uint32_t k = 0; k < run.count; ++k) {
        results[run.first + k] = KmerLookupResult{};
      }
    }
  });

  if (phases) {
    *phases = phase_diagnostics({{"resize", resize_ms}, {"clear", clear_ms}});
  }
}

// Writes the lookups of one table into `results`, which prepare_results
// has made ready for this file, looking the keys up in the order of
// `order` (packed keys: bucket_order or sorted_order). Phases: keys (the
// 16-mer lookups), short (the tail lookups).
template <typename Policy>
void probe_order(const typename Policy::Table& table, const KeyedInput& keyed,
                 const std::vector<std::uint64_t>& order,
                 std::vector<KmerLookupResult>& results,
                 Diagnostics* diagnostics = nullptr) {
  const auto before = table.stats();

  const double keys_ms = time_ms([&] {
    typename Policy::State state{};

    for (const std::uint64_t packed : order) {
      results[position_of(packed)] =
          Policy::lookup(table, state, key_of(packed));
    }
  });

  const double short_ms = time_ms([&] {
    for (const ShortQuery& query : keyed.short_queries) {
      KmerLookupResult& result = results[query.position];
      result = Policy::tail(table, query.key, query.length);
      result.found = false;
    }
  });

  phase_barrier(results.data());

  if (diagnostics) {
    *diagnostics = phase_diagnostics({{"keys", keys_ms},
                                      {"short", short_ms}}) +
                   Policy::counters(before, table.stats());
  }
}

// The chained lengths against the baseline's. A mismatch is "explained by
// runs" when the baseline's match at i covers a run position: those
// positions have no lookup yet, so the chain cannot carry a match through
// them. Anything else is a real error.
struct ChainComparison {
  std::size_t positions = 0;
  std::size_t equal = 0;
  std::size_t at_runs = 0;       // mismatches at a run position itself
  std::size_t into_runs = 0;     // mismatches whose true match reaches a run
  std::size_t unexplained = 0;   // everything else: a bug
  std::size_t first_unexplained = 0;
  std::uint32_t first_want = 0;
  std::uint32_t first_got = 0;
};

inline ChainComparison compare_chain(const MatchingStatistics& want,
                                     const MatchingStatistics& got,
                                     const KeyedInput& keyed) {
  ChainComparison comparison;
  const std::size_t n = want.size();
  comparison.positions = n;

  // Walk backward, keeping the nearest run position at or after i.
  std::size_t run = keyed.runs.size();  // runs[run] is the next run after i
  std::size_t next_run_position = n;    // n: none

  for (std::size_t step = 0; step < n; ++step) {
    const std::size_t i = n - 1 - step;

    while (run > 0 && keyed.runs[run - 1].first + keyed.runs[run - 1].count > i) {
      --run;
    }
    // runs[run - 1] (if any) ends at or before i; runs[run] starts after
    // or contains i.
    bool inside = false;
    if (run < keyed.runs.size()) {
      const HomopolymerRun& r = keyed.runs[run];
      inside = r.first <= i;
      next_run_position = inside ? i : r.first;
    } else {
      next_run_position = n;
    }

    const std::uint32_t want_length = want[i].second;
    const std::uint32_t got_length = got[i].second;

    if (want_length == got_length) {
      ++comparison.equal;
    } else if (inside) {
      ++comparison.at_runs;
    } else if (next_run_position < n && i + want_length > next_run_position) {
      ++comparison.into_runs;
    } else {
      if (comparison.unexplained++ == 0 || i < comparison.first_unexplained) {
        comparison.first_unexplained = i;
        comparison.first_want = want_length;
        comparison.first_got = got_length;
      }
    }
  }

  return comparison;
}

// Checks up to `samples` evenly spaced positions: the chain's reported
// position must really match the input for the reported length.
inline std::size_t bad_chain_positions(const MatchingStatistics& got,
                                       const std::vector<unsigned char>& input,
                                       const std::vector<unsigned char>& reference,
                                       const std::size_t samples) {
  if (got.empty()) return 0;

  const std::size_t stride = std::max<std::size_t>(1, got.size() / samples);
  std::size_t bad = 0;

  for (std::size_t i = 0; i < got.size(); i += stride) {
    const std::size_t p = got[i].first;
    const std::size_t len = got[i].second;

    if (len == 0) continue;

    if (p + len > reference.size() || i + len > input.size() ||
        !std::equal(input.begin() + static_cast<std::ptrdiff_t>(i),
                    input.begin() + static_cast<std::ptrdiff_t>(i + len),
                    reference.begin() + static_cast<std::ptrdiff_t>(p))) {
      ++bad;
    }
  }

  return bad;
}

}  // namespace ms2
