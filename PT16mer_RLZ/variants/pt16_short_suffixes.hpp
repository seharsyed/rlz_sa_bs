#pragma once

// Separators and short suffixes, shared by the PT16 table variants
// (pt16_rlz_v2_fastmiss.hpp, pt16_sassy.hpp) and their builders.
//
// SEPARATORS. Only A, C, G and T are sequence characters. Any other byte
// (N, IUPAC codes, lowercase, ...) is a separator: it never matches
// anything, in the reference or in the input. So a 16-mer window that
// contains one has no table entry, and a match never runs across one.
//
// SHORT SUFFIXES. A reference position whose run of ACGT characters is
// shorter than 16 -- it is within 15 characters of the reference's end or
// of a separator -- has no 16-mer, so it has no table entry either; but a
// query can still match it, for up to that run's length. Every such
// position is kept as a short suffix record (PackedShortSuffix), and a
// lookup that misses the table matches the records too, which keeps the
// longest match it reports exact. The reference's end is one such case;
// each separator adds up to 15 more positions just before it.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"      // KMER_LENGTH, LOW_BITS, alphatab, acgt_runs
#include "pt16_sassy_format.hpp"  // PackedShortSuffix

// Every short suffix of `reference`: each position whose ACGT run is 1 to
// 15 characters long, packed like a 16-mer key (unused low bits 0). Records
// with the same characters are kept once (the first position), since any
// one occurrence is a correct answer. Sorted by length, then characters.
inline std::vector<PackedShortSuffix> collect_short_suffixes(
    const std::vector<unsigned char>& reference) {
  if (reference.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("reference too long for 32-bit positions");
  }

  const std::vector<std::uint8_t> run = acgt_runs(reference);
  std::vector<PackedShortSuffix> suffixes;

  for (std::size_t p = 0; p < reference.size(); ++p) {
    const std::uint32_t length = run[p];

    if (length == 0 || length >= KMER_LENGTH) {
      continue;
    }

    std::uint32_t packed = 0;

    for (std::uint32_t j = 0; j < length; ++j) {
      const std::uint32_t code =
          alphatab[static_cast<unsigned char>(reference[p + j])];
      packed |= code << (30U - 2U * j);
    }

    suffixes.push_back({packed, static_cast<std::uint32_t>(p), length});
  }

  std::sort(suffixes.begin(), suffixes.end(),
            [](const PackedShortSuffix& a, const PackedShortSuffix& b) {
              return std::tie(a.length, a.packed, a.ref_pos) <
                     std::tie(b.length, b.packed, b.ref_pos);
            });

  suffixes.erase(
      std::unique(suffixes.begin(), suffixes.end(),
                  [](const PackedShortSuffix& a, const PackedShortSuffix& b) {
                    return a.length == b.length && a.packed == b.packed;
                  }),
      suffixes.end());

  return suffixes;
}

/**
 * The short suffixes, arranged for the two questions a lookup asks:
 *
 *   best_in_bucket(bucket)  for an empty bucket, precomputed at load time:
 *                           the longest match any short suffix has with
 *                           the bucket's 8 characters (0 to 8), and where.
 *                           One table per prefix length 1..8 (4^l entries)
 *                           answers it without scanning the records.
 *
 *   raise_long(key, ...)    at query time: a short suffix longer than 8
 *                           characters in the key's own bucket may match
 *                           more than 8; only those can beat a miss that
 *                           already shares the bucket's 8 characters.
 *                           has_long(bucket) says whether there are any.
 */
class ShortSuffixIndex {
 public:
  ShortSuffixIndex() = default;

  explicit ShortSuffixIndex(const std::vector<PackedShortSuffix>& suffixes)
      : count_(suffixes.size()) {
    for (std::uint32_t l = 1; l <= 8; ++l) {
      prefix_position_[l].assign(std::size_t{1} << (2 * l), none);
    }

    for (const PackedShortSuffix& suffix : suffixes) {
      for (std::uint32_t l = 1; l <= std::min<std::uint32_t>(suffix.length, 8);
           ++l) {
        std::uint32_t& slot = prefix_position_[l][suffix.packed >> (32 - 2 * l)];

        if (slot == none) {
          slot = suffix.ref_pos;
        }
      }

      if (suffix.length > 8) {
        const std::uint32_t bucket = suffix.packed >> LOW_BITS;
        long_.push_back(suffix);
        long_buckets_[bucket / 64] |= std::uint64_t{1} << (bucket % 64);
      }
    }

    // By bucket, and longest first within a bucket.
    std::sort(long_.begin(), long_.end(),
              [](const PackedShortSuffix& a, const PackedShortSuffix& b) {
                const std::uint32_t a_bucket = a.packed >> LOW_BITS;
                const std::uint32_t b_bucket = b.packed >> LOW_BITS;
                return a_bucket != b_bucket ? a_bucket < b_bucket
                                            : a.length > b.length;
              });
  }

  // The longest match (0 to 8 characters) of any short suffix with the
  // bucket's 8 characters, and one reference position with it.
  std::pair<std::uint32_t, std::uint32_t> best_in_bucket(
      const std::uint32_t bucket) const {
    for (std::uint32_t l = 8; l >= 1; --l) {
      const std::uint32_t position =
          prefix_position_[l][bucket >> (2 * (8 - l))];

      if (position != none) {
        return {l, position};
      }
    }

    return {0, 0};
  }

  bool has_long(const std::uint32_t bucket) const {
    return (long_buckets_[bucket / 64] >> (bucket % 64)) & 1U;
  }

  // Raises (length, position) if a short suffix longer than 8 in the key's
  // bucket shares a longer prefix with `key`. Longest first, stopping at
  // the first suffix no longer than the match already found.
  template <typename Length, typename Position>
  void raise_long(const std::uint32_t key, Length& length,
                  Position& position) const {
    const std::uint32_t bucket = key >> LOW_BITS;

    auto it = std::lower_bound(
        long_.begin(), long_.end(), bucket,
        [](const PackedShortSuffix& suffix, const std::uint32_t value) {
          return (suffix.packed >> LOW_BITS) < value;
        });

    for (; it != long_.end() && (it->packed >> LOW_BITS) == bucket &&
           it->length > length;
         ++it) {
      // Capped at the suffix's own length: its unused low bits are 0 and
      // must not count as matches.
      const std::uint32_t shared = std::min<std::uint32_t>(
          it->length,
          static_cast<std::uint32_t>(std::countl_zero(key ^ it->packed) / 2));

      if (shared > length) {
        length = static_cast<Length>(shared);
        position = static_cast<Position>(it->ref_pos);
      }
    }
  }

  std::size_t size() const { return count_; }

  std::size_t memory_bytes() const {
    std::size_t bytes = sizeof(long_buckets_) +
                        long_.size() * sizeof(PackedShortSuffix);

    for (const auto& level : prefix_position_) {
      bytes += level.size() * sizeof(std::uint32_t);
    }

    return bytes;
  }

 private:
  static constexpr std::uint32_t none =
      std::numeric_limits<std::uint32_t>::max();

  std::size_t count_ = 0;

  // prefix_position_[l][p]: a reference position of a short suffix at
  // least l long whose first l characters are p (packed), or `none`.
  std::array<std::vector<std::uint32_t>, 9> prefix_position_;

  // The suffixes longer than 8, by bucket then longest first, and one bit
  // per bucket that has any.
  std::vector<PackedShortSuffix> long_;
  std::array<std::uint64_t, NUMBER_OF_BUCKETS / 64> long_buckets_{};
};
