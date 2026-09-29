#pragma once

// Preparing an input for 16-mer table lookups, for the ms2 benchmark
// (ms2_main.cpp). Separate from lrf_ms/probe_pipeline.hpp on purpose: the
// new matching-statistics methods are built up here.
//
// One pass over the input (prepare_keys) splits its positions three ways:
//
//   keys          positions whose 16-character window is all ACGT, as
//                 packed (key << 32 | position) values in text order --
//                 except the two below;
//
//   runs          positions whose 16-mer is all A (key 0) or all T (key
//                 0xFFFFFFFF), recorded as runs of consecutive positions
//                 instead of as keys. They are not queried (yet): a
//                 homopolymer 16-mer has a huge range in the table, and it
//                 is the chain's worst case, so they are dealt with apart.
//                 With keep_runs they are ordinary keys instead;
//
//   short queries positions without a full ACGT window: the last 15 before
//                 a separator (any non-ACGT byte) or the input's end, with
//                 their 1 to 15 ACGT characters packed like a key;
//
//   separators    positions holding a non-ACGT byte. Nothing is stored for
//                 them, only their count: their match length is 0, so they
//                 need no search.
//
// The keys are then ordered for probing, by bucket (the high 16 bits of
// the key, bucket_order) or fully by key (sorted_order). Both keep equal
// keys (or equal buckets) in text order.

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"  // KMER_LENGTH, alphatab, is_acgt, Diagnostics

namespace ms2 {

// A run of consecutive positions whose 16-mer is all A or all T: positions
// first .. first + count - 1. The characters span count + 15 of them.
struct HomopolymerRun {
  std::uint32_t first;
  std::uint32_t count;
  char base;  // 'A' or 'T'
};

// A position without a full ACGT window: its `length` (1 to 15) ACGT
// characters, packed from the top bit down like a key (unused bits 0).
struct ShortQuery {
  std::uint32_t position;
  std::uint32_t key;
  std::uint32_t length;
};

inline constexpr std::uint32_t POLY_A_KEY = 0;
inline constexpr std::uint32_t POLY_T_KEY = 0xFFFFFFFFU;

inline std::uint32_t key_of(const std::uint64_t packed) {
  return static_cast<std::uint32_t>(packed >> 32);
}

inline std::uint32_t position_of(const std::uint64_t packed) {
  return static_cast<std::uint32_t>(packed);
}

struct KeyedInput {
  std::size_t input_size = 0;

  // key << 32 | position, in text order.
  std::vector<std::uint64_t> keys;

  // In text order.
  std::vector<HomopolymerRun> runs;
  std::vector<ShortQuery> short_queries;

  // Positions holding a non-ACGT byte: match length 0, nothing stored.
  std::size_t separators = 0;

  std::size_t run_positions(const char base) const {
    std::size_t total = 0;
    for (const HomopolymerRun& run : runs) {
      if (run.base == base) total += run.count;
    }
    return total;
  }

  std::size_t run_count(const char base) const {
    std::size_t total = 0;
    for (const HomopolymerRun& run : runs) {
      if (run.base == base) ++total;
    }
    return total;
  }
};

/**
 * The single pass. Rolls the key of the 16 characters ending at each
 * input position, and the length of the ACGT run ending there; once the
 * run is 16 long, the window starting 15 positions back is complete, and
 * that position is a key or part of a homopolymer run. A position whose
 * window never completes (a separator comes, or the input ends, first)
 * is a short query, answered from the characters it does have.
 *
 * `keep_runs`: all-A / all-T 16-mers are kept as ordinary keys (looked up
 * like any other), and no runs are recorded.
 */
template <typename Symbol>
void prepare_keys(const std::vector<Symbol>& input, KeyedInput& out,
                  const bool keep_runs = false) {
  if (input.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("input too long for 32-bit positions");
  }

  const std::size_t n = input.size();

  out.input_size = n;
  out.keys.clear();
  out.runs.clear();
  out.short_queries.clear();
  out.separators = 0;
  out.keys.reserve(n);

  std::uint32_t key = 0;
  std::uint32_t run = 0;  // ACGT characters ending here, capped at 16

  // The positions before the current ACGT run's window became complete
  // are the pending short-query candidates: when a separator comes they
  // are the last (up to 15) positions of the run, answered as short
  // queries. See flush_short below.
  std::size_t run_start = 0;  // first position of the current ACGT run

  const auto flush_short = [&](const std::size_t end) {
    // Positions max(run_start, end - 15) .. end - 1 have fewer than 16
    // ACGT characters before `end` (a separator or the input's end).
    const std::size_t first =
        end - run_start >= KMER_LENGTH ? end - (KMER_LENGTH - 1) : run_start;

    for (std::size_t p = first; p < end; ++p) {
      const std::uint32_t length = static_cast<std::uint32_t>(end - p);
      std::uint32_t packed = 0;

      for (std::uint32_t j = 0; j < length; ++j) {
        packed |= static_cast<std::uint32_t>(
                      alphatab[static_cast<unsigned char>(input[p + j])])
                  << (30U - 2U * j);
      }

      out.short_queries.push_back(
          {static_cast<std::uint32_t>(p), packed, length});
    }
  };

  for (std::size_t j = 0; j < n; ++j) {
    const unsigned char c = static_cast<unsigned char>(input[j]);

    if (!is_acgt(c)) {
      flush_short(j);
      ++out.separators;
      run = 0;
      key = 0;
      run_start = j + 1;
      continue;
    }

    key = (key << 2) | alphatab[c];

    if (run < KMER_LENGTH) {
      ++run;
    }

    if (run < KMER_LENGTH) {
      continue;
    }

    const std::uint32_t position =
        static_cast<std::uint32_t>(j + 1 - KMER_LENGTH);

    if (!keep_runs && (key == POLY_A_KEY || key == POLY_T_KEY)) {
      const char base = key == POLY_A_KEY ? 'A' : 'T';

      if (!out.runs.empty() && out.runs.back().base == base &&
          out.runs.back().first + out.runs.back().count == position) {
        ++out.runs.back().count;
      } else {
        out.runs.push_back({position, 1, base});
      }

      continue;
    }

    out.keys.push_back(static_cast<std::uint64_t>(key) << 32 | position);
  }

  flush_short(n);
}

// ---------- Bucket order: by the high 16 bits of the key ----------

// A stable counting sort of `keys` by bucket into `out` (resized to fit).
// Phases: alloc (growing `out`: only when it is smaller than before, e.g.
// on the first file), count, prefix, scatter.
inline void bucket_order(const std::vector<std::uint64_t>& keys,
                         std::vector<std::uint64_t>& out,
                         Diagnostics* phases = nullptr) {
  constexpr std::size_t buckets = std::size_t{1} << 16;
  std::vector<std::uint32_t> start(buckets + 1, 0);

  const double alloc_ms = time_ms([&] { out.resize(keys.size()); });

  const double count_ms = time_ms([&] {
    for (const std::uint64_t packed : keys) {
      ++start[(packed >> 48) + 1];
    }
  });

  const double prefix_ms = time_ms([&] {
    for (std::size_t b = 0; b < buckets; ++b) {
      start[b + 1] += start[b];
    }
  });

  const double scatter_ms = time_ms([&] {
    for (const std::uint64_t packed : keys) {
      out[start[packed >> 48]++] = packed;
    }
  });

  if (phases) {
    *phases = phase_diagnostics({{"alloc", alloc_ms},
                                 {"count", count_ms},
                                 {"prefix", prefix_ms},
                                 {"scatter", scatter_ms}});
  }
}

// ---------- Sorted order: fully by key ----------

// An LSD radix sort of `keys` by their 32-bit key, in three stable passes
// of 11, 11 and 10 bits, into `out` (resized to fit); `scratch` is the
// other buffer. The histograms of all three passes are counted in one
// read. Stable, so equal keys stay in text order. Phases: alloc (growing
// `out` and `scratch`, as in bucket_order), count, pass-1, pass-2, pass-3.
inline void sorted_order(const std::vector<std::uint64_t>& keys,
                         std::vector<std::uint64_t>& out,
                         std::vector<std::uint64_t>& scratch,
                         Diagnostics* phases = nullptr) {
  constexpr std::array<unsigned, 3> shift{32, 43, 54};
  constexpr std::array<unsigned, 3> width{11, 11, 10};

  std::array<std::vector<std::uint32_t>, 3> start;
  for (std::size_t p = 0; p < 3; ++p) {
    start[p].assign((std::size_t{1} << width[p]) + 1, 0);
  }

  const double alloc_ms = time_ms([&] {
    out.resize(keys.size());
    scratch.resize(keys.size());
  });

  const double count_ms = time_ms([&] {
    for (const std::uint64_t packed : keys) {
      for (std::size_t p = 0; p < 3; ++p) {
        const std::uint64_t digit =
            (packed >> shift[p]) & ((std::uint64_t{1} << width[p]) - 1);
        ++start[p][digit + 1];
      }
    }

    for (std::size_t p = 0; p < 3; ++p) {
      for (std::size_t d = 0; d + 1 < start[p].size(); ++d) {
        start[p][d + 1] += start[p][d];
      }
    }
  });

  // keys -> scratch -> out -> scratch would end in the wrong buffer, so:
  // keys -> out -> scratch -> out.
  const std::array<const std::vector<std::uint64_t>*, 3> from{&keys, &out,
                                                              &scratch};
  const std::array<std::vector<std::uint64_t>*, 3> to{&out, &scratch, &out};

  std::array<double, 3> pass_ms{};

  for (std::size_t p = 0; p < 3; ++p) {
    pass_ms[p] = time_ms([&] {
      const std::uint64_t mask = (std::uint64_t{1} << width[p]) - 1;
      std::vector<std::uint32_t>& next = start[p];
      std::vector<std::uint64_t>& target = *to[p];

      for (const std::uint64_t packed : *from[p]) {
        target[next[(packed >> shift[p]) & mask]++] = packed;
      }
    });
  }

  if (phases) {
    *phases = phase_diagnostics({{"alloc", alloc_ms},
                                 {"count", count_ms},
                                 {"pass-1", pass_ms[0]},
                                 {"pass-2", pass_ms[1]},
                                 {"pass-3", pass_ms[2]}});
  }
}

}  // namespace ms2
