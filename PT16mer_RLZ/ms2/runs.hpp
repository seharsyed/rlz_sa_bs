#pragma once

// Shared by lazy_ms.hpp and sampled_ms.hpp: the reference's longest run of
// each character, and the exact MS at i from a lookup made at i.
//
// Long single-character runs (e.g. N runs replaced by A when cleaning) are
// the hard case for both methods: inside an input run longer than the
// reference's longest run of that character, every position needs its own
// answer, and each lookup of c^16 meets a range with one occurrence per
// position of the reference's runs. Two things keep it cheap:
//
//   run shortcut  if input[i .. i + k) is a run of c with k > R_c, the
//                 reference's longest run of c, then MS[i] = (that run's
//                 start, R_c) exactly: c^(R_c + 1) does not occur and c^R_c
//                 does. O(1), no lookup.
//
//   mlr search    a range is resolved with rlz::mlrLongestMatch (Manber and
//                 Myers' binary search keeping the LCP with both bounds):
//                 O(log c + L) character comparisons for the longest match
//                 of length L among c occurrences, instead of two binary
//                 searches per character (O(L log c) random reads).

#include <array>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"  // KmerLookupResult, KMER_LENGTH
#include "../rlz_common.hpp"  // rlz::mlrLongestMatch

namespace ms2 {

struct ReferenceRuns {
  std::array<std::uint32_t, 256> longest{};  // per byte value
  std::array<std::uint32_t, 256> start{};    // where that run starts

  ReferenceRuns() = default;

  template <typename Symbol>
  explicit ReferenceRuns(const std::vector<Symbol>& reference) {
    std::size_t i = 0;
    while (i < reference.size()) {
      std::size_t j = i + 1;
      while (j < reference.size() && reference[j] == reference[i]) ++j;
      const auto c = static_cast<unsigned char>(reference[i]);
      if (j - i > longest[c]) {
        longest[c] = static_cast<std::uint32_t>(j - i);
        start[c] = static_cast<std::uint32_t>(i);
      }
      i = j;
    }
  }

  // MS at a position where a run of `run_length` characters `c` starts, if
  // the run shortcut applies (run_length > the reference's longest run).
  bool shortcut(const unsigned char c, const std::size_t run_length,
                std::pair<std::uint32_t, std::uint32_t>& answer) const {
    if (run_length <= longest[c]) return false;
    answer = {start[c], longest[c]};
    return true;
  }
};

// The exact MS at i from a lookup made at i: a miss or tail as is, a
// singleton extended, a range resolved by the mlr search. `work` counts the
// characters matched past 16.
template <typename Symbol>
std::pair<std::uint32_t, std::uint32_t> resolve_hit(
    const KmerLookupResult& hit, const std::vector<Symbol>& input,
    const std::vector<Symbol>& reference, const std::size_t i,
    std::size_t& work) {
  if (!hit.found) return {hit.match_position, hit.match_length};

  if (hit.count > 1) {
    const auto [position, length] = rlz::mlrLongestMatch(
        reference, hit.positions, 0, hit.positions.size() - 1, input,
        i + KMER_LENGTH, KMER_LENGTH);
    work += length - KMER_LENGTH;
    return {static_cast<std::uint32_t>(position),
            static_cast<std::uint32_t>(length)};
  }

  const std::size_t match = hit.match_position;
  std::size_t j = i + KMER_LENGTH;
  std::size_t offset = KMER_LENGTH;
  while (j < input.size() && match + offset < reference.size() &&
         reference[match + offset] == input[j]) {
    ++j;
    ++offset;
  }
  work += offset - KMER_LENGTH;
  return {static_cast<std::uint32_t>(match), static_cast<std::uint32_t>(offset)};
}

// Whether `wanted` is one of a hit's occurrences. All occurrences of a
// found 16-mer share it, so for a range this compares the 16 characters at
// `wanted` with those at its first occurrence instead of scanning them all.
template <typename Symbol>
bool occurs_at(const KmerLookupResult& result, const std::uint32_t wanted,
              const std::vector<Symbol>& reference) {
  if (result.count == 1) return result.match_position == wanted;
  if (result.count == 0 || result.positions.empty()) return false;
  const std::size_t first = result.positions[0];
  if (static_cast<std::size_t>(wanted) + KMER_LENGTH > reference.size()) return false;
  for (std::size_t k = 0; k < KMER_LENGTH; ++k) {
    if (reference[wanted + k] != reference[first + k]) return false;
  }
  return true;
}

}  // namespace ms2
