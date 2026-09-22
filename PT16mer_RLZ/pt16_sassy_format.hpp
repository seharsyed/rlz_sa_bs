#pragma once

// The parts of the sassy PT16 table's binary format that both the builder
// (pt16_build_sassy.hpp) and the reader (pt16_sassy.hpp) need to agree on:
// the short-suffix records, and the L entry bit layout. Neither is used by
// the older PT16 formats (pt16_build.hpp/pt16_rlz.hpp,
// pt16_build_v2.hpp/pt16_rlz_v2*.hpp), which is why they live apart from the
// rest of pt16_utils.hpp. Keeping the encode and decode sides of each here,
// used by nothing else, means the two sides cannot drift apart the way the
// pre-PT16SA03 format's duplicated bit math once did.

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "pt16_utils.hpp"  // KMER_LENGTH, LOW_BITS, alphatab

// ---------- Short suffixes ----------

// A suffix of the reference shorter than KMER_LENGTH has no 16-mer, so it gets
// no table entry. The sassy table keeps each of them (down to
// SHORT_SUFFIX_MIN_LENGTH characters) as a record, to be matched when a
// lookup misses the table. Keeping every length makes the longest match found
// by a lookup exact: a short suffix can be the best match for a 16-mer whose
// bucket is empty in the table, and there is no reason to leave the very
// shortest ones out, since matching one is a single XOR.
constexpr std::uint32_t SHORT_SUFFIX_MIN_LENGTH = 1;

struct PackedShortSuffix {
  // The characters, 2 bits each from the top bit down, exactly as in a 16-mer
  // key; the unused low bits are 0. So packed >> LOW_BITS is the bucket, and
  // the number of leading characters shared with a 16-mer key q is
  // min(length, countl_zero(q ^ packed) / 2).
  std::uint32_t packed;

  // Where the suffix starts in the reference.
  std::uint32_t ref_pos;

  std::uint32_t length;
};

static_assert(sizeof(PackedShortSuffix) == 3 * sizeof(std::uint32_t),
              "PackedShortSuffix is written to the table as raw bytes");

// One record for each length from SHORT_SUFFIX_MIN_LENGTH up to
// KMER_LENGTH - 1 that the reference is long enough to have, shortest first.
inline std::vector<PackedShortSuffix> build_packed_short_suffixes(
    const std::vector<unsigned char>& reference) {
  if (reference.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("reference too long for 32-bit positions");
  }

  std::vector<PackedShortSuffix> suffixes;

  for (std::uint32_t length = SHORT_SUFFIX_MIN_LENGTH; length < KMER_LENGTH;
       ++length) {
    // A reference shorter than `length` has no suffix of that length, and
    // neither does it have a longer one.
    if (length > reference.size()) {
      break;
    }

    const std::uint32_t ref_pos =
        static_cast<std::uint32_t>(reference.size() - length);

    std::uint32_t packed = 0;

    for (std::uint32_t j = 0; j < length; ++j) {
      const std::uint32_t code =
          alphatab[static_cast<unsigned char>(reference[ref_pos + j])];
      packed |= code << (30U - 2U * j);
    }

    suffixes.push_back({packed, ref_pos, length});
  }

  return suffixes;
}

// ---------- Sassy L entry ----------
//
// Each 64-bit L entry (pt16_build_sassy.hpp / pt16_sassy.hpp) has the low 16
// bits of the 16-mer in its top 16 bits; the bucket (the high 16 bits of the
// key) selects which bucket the entry belongs to. The remaining 48 bits hold
// either a singleton's text position or a range's slice into sampled_sa,
// chosen by the lowest bit:
//
//   bit 0        0: one occurrence (singleton)   1: several (a range)
//
//   singleton:   bits 1..47   the text position (47 bits)
//
//   range:       bits 1..16   the occurrence count, or SASSY_COUNT_ESCAPE
//                             (0xFFFF) if the real count does not fit
//                bits 17..47  the offset of this entry's slice into
//                             sampled_sa, relative to H_sa[bucket] (31 bits)
//
// A range whose true count is SASSY_COUNT_ESCAPE or more stores that count as
// a plain uint32_t at the start of its slice, then its positions right after
// it: the escape trades one extra sampled_sa slot for room to record any
// count up to 2^32 - 1, which covers every count a real genome can produce (a
// single 16-mer cannot occur more than reference.size() times, and
// reference.size() already fits in uint32_t elsewhere in this codebase).
//
// 31 bits for a range's offset and 47 bits for a singleton's position are
// both far beyond what a real genome needs: a bucket would need close to
// 2^31 occurrences accumulated in it, which means essentially the whole
// reference sharing one 8-character prefix, and that does not happen in
// practice. Only a single 16-mer's own occurrence count can plausibly need
// more than 16 bits (a highly conserved repeat family), which is why only it
// has an escape.

constexpr int SASSY_COUNT_BITS = 16;
constexpr int SASSY_OFFSET_BITS = 31;
constexpr int SASSY_POSITION_BITS = 47;
constexpr std::uint16_t SASSY_COUNT_ESCAPE = 0xFFFFU;

constexpr std::uint64_t sassy_position_mask =
    (std::uint64_t{1} << SASSY_POSITION_BITS) - 1;
constexpr std::uint64_t sassy_count_mask =
    (std::uint64_t{1} << SASSY_COUNT_BITS) - 1;
constexpr std::uint64_t sassy_offset_mask =
    (std::uint64_t{1} << SASSY_OFFSET_BITS) - 1;

inline std::uint64_t sassy_encode_low(const std::uint16_t low) {
  return static_cast<std::uint64_t>(low) << 48;
}

inline std::uint16_t sassy_decode_low(const std::uint64_t entry) {
  return static_cast<std::uint16_t>(entry >> 48);
}

inline bool sassy_is_range(const std::uint64_t entry) {
  return (entry & 1U) != 0;
}

// A singleton entry: one occurrence, at `position`.
inline std::uint64_t sassy_encode_singleton(const std::uint16_t low,
                                            const std::uint64_t position) {
  if (position > sassy_position_mask) {
    throw std::runtime_error(
        "PT16 sassy: text position does not fit in an entry");
  }

  return sassy_encode_low(low) | (position << 1);
}

inline std::uint64_t sassy_decode_position(const std::uint64_t entry) {
  return (entry >> 1) & sassy_position_mask;
}

// A range entry, whose slice into sampled_sa starts at `offset` relative to
// the bucket's H_sa. `count_field` is what to store in the 16-bit count
// field: the real count if it fits, or SASSY_COUNT_ESCAPE if the caller has
// arranged for the real count to be readable from sampled_sa instead (see
// sassy_decode_count_field and the layout above).
inline std::uint64_t sassy_encode_range(const std::uint16_t low,
                                        const std::uint32_t offset,
                                        const std::uint16_t count_field) {
  if (offset > sassy_offset_mask) {
    throw std::runtime_error(
        "PT16 sassy: range offset does not fit in an entry");
  }

  return sassy_encode_low(low) | (static_cast<std::uint64_t>(offset) << 17) |
         (static_cast<std::uint64_t>(count_field) << 1) | 1U;
}

inline std::uint32_t sassy_decode_offset(const std::uint64_t entry) {
  return static_cast<std::uint32_t>((entry >> 17) & sassy_offset_mask);
}

inline std::uint16_t sassy_decode_count_field(const std::uint64_t entry) {
  return static_cast<std::uint16_t>((entry >> 1) & sassy_count_mask);
}
