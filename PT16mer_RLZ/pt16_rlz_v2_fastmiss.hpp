#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pt16_utils.hpp"  // KmerLookupResult, EntryComposition, constants
#include "rlz_common.hpp"

/**
 * PT16FastMissParser: the v2 H/L table (same file, same lookup results as
 * PT16RLZParser::lookupKmerByKey), rearranged at load time so that a MISS
 * -- often the majority of lookups -- is as cheap as possible:
 *
 *   1. Every L entry also stores its first reference position (ref_pos),
 *      the value v2 recovers with a random suffix-array read
 *      (sa_[H_sa_[bucket] + sa_offset]) on every miss and every hit. Here
 *      that read is gone: the entry the bucket search just landed on
 *      already holds it. Costs 4 more bytes per entry (8 instead of 4).
 *
 *   2. Empty buckets: the whole answer (table LCP neighbour AND the best
 *      short suffix) depends only on the bucket, so it is computed once
 *      per empty bucket at load time. A query into an empty bucket is one
 *      array read.
 *
 *   3. Short-suffix check skipped on almost every miss. A non-empty-bucket
 *      miss already shares >= 8 characters with its neighbour, so only a
 *      short suffix of length >= 9 whose first 8 characters ARE this
 *      bucket could beat it -- at most 7 buckets in the whole table. Those
 *      buckets are flagged in a 65536-bit bitmap (8 KB); every other miss
 *      skips the check with one bit test. The same flag covers the empty
 *      buckets whose precomputed answer (capped at 8 characters) could
 *      still be beaten by such a suffix.
 *
 * Hits are unchanged from v2 apart from (1): count/positions still come
 * from the suffix array itself.
 */
template <typename T1, typename T2>
class PT16FastMissParser {
 public:
  using input_type = std::vector<T1>;
  using reference_type = std::vector<T1>;
  using suffix_array_type = std::vector<T2>;
  using KmerLookupResult = ::KmerLookupResult;

  static constexpr std::uint32_t kmer_length = KMER_LENGTH;
  static constexpr std::uint32_t low_bits = LOW_BITS;
  static constexpr std::uint32_t low_mask = LOW_MASK;
  static constexpr std::uint32_t number_of_buckets = NUMBER_OF_BUCKETS;
  static constexpr std::uint32_t empty_bucket_flag = EMPTY_BUCKET_FLAG;
  static constexpr std::uint32_t empty_bucket_mask = empty_bucket_flag - 1;
  static constexpr std::uint32_t binary_search_threshold =
      BINARY_SEARCH_THRESHOLD;
  static constexpr std::uint16_t large_offset_flag = LARGE_OFFSET_FLAG;

  struct Stats {
    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t singleton_hits = 0;
    std::size_t range_hits = 0;
    std::size_t entries = 0;
    std::size_t approx_bytes = 0;

    std::size_t linear_bucket_searches = 0;
    std::size_t binary_bucket_searches = 0;

    // Misses that fell into an empty bucket (answered from the
    // precomputed table), and misses that still had to run the full
    // short-suffix check (flagged buckets only).
    std::size_t empty_bucket_misses = 0;
    std::size_t short_suffix_checks = 0;
  };

  PT16FastMissParser(const reference_type& ref, const suffix_array_type& sa,
                     const std::string& pt16_path)
      : ref_(&ref), sa_(&sa) {
    load_hl(pt16_path);
    build_short_suffixes();
    build_empty_answers();

    stats_.entries = L_.size();
    stats_.approx_bytes =
        H_.size() * sizeof(std::uint32_t) +
        H_sa_.size() * sizeof(std::uint32_t) + L_.size() * sizeof(Entry) +
        large_offsets_.size() *
            sizeof(std::pair<const std::uint32_t, std::uint32_t>) +
        empty_ref_pos_.size() * sizeof(std::uint32_t) +
        empty_length_.size() * sizeof(std::uint8_t) +
        sizeof(long_short_buckets_);
  }

  PT16FastMissParser(const PT16FastMissParser&) = delete;
  PT16FastMissParser& operator=(const PT16FastMissParser&) = delete;

  /**
   * Same result as PT16RLZParser::lookupKmerByKey for the same key (found,
   * count, positions, match_length; match_position too, as long as ties
   * between equally long short suffixes break the same way). Needs only
   * the key: nothing here reads the input.
   */
  KmerLookupResult lookupKmerByKey(const std::uint32_t key) const {
    const std::uint32_t bucket = key >> low_bits;
    const std::uint32_t h = H_[bucket];

    KmerLookupResult result;

    // ---------- Empty bucket: precomputed answer ----------

    if (h & empty_bucket_flag) {
      std::size_t ref_pos = empty_ref_pos_[bucket];
      std::size_t match_length = empty_length_[bucket];

      if (has_long_short_suffix(bucket)) {
        ++stats_.short_suffix_checks;
        check_short_suffixes(key, ref_pos, match_length);
      }

      ++stats_.misses;
      ++stats_.empty_bucket_misses;
      result.match_position = static_cast<std::uint32_t>(ref_pos);
      result.match_length = static_cast<std::uint32_t>(match_length);
      return result;
    }

    // ---------- Non-empty bucket: search L ----------

    const std::uint16_t low = static_cast<std::uint16_t>(key & low_mask);
    const std::uint32_t begin = h;
    const std::uint32_t end = H_[next_nonempty_bucket(bucket)];

    std::uint32_t position = begin;

    if (end - begin < binary_search_threshold) {
      ++stats_.linear_bucket_searches;

      while (position < end && L_[position].low < low) {
        ++position;
      }
    } else {
      ++stats_.binary_bucket_searches;

      const auto it = std::lower_bound(
          L_.begin() + begin, L_.begin() + end, low,
          [](const Entry& entry, const std::uint16_t value) {
            return entry.low < value;
          });

      position = static_cast<std::uint32_t>(it - L_.begin());
    }

    // ---------- Hit ----------

    if (position < end && L_[position].low == low) {
      ++stats_.hits;

      const std::uint32_t sa_start = sa_start_at(bucket, position);
      const std::uint32_t sa_end = interval_end(bucket, position, end);

      result.found = true;
      result.match_length = kmer_length;
      result.count = sa_end - sa_start + 1;
      result.match_position = L_[position].ref_pos;
      result.positions = std::span<const std::uint32_t>(
          sa_->data() + sa_start, static_cast<std::size_t>(result.count));

      if (sa_start == sa_end) {
        ++stats_.singleton_hits;
      } else {
        ++stats_.range_hits;
      }

      return result;
    }

    // ---------- Miss: neighbour with the longest common prefix ----------

    ++stats_.misses;

    std::uint32_t best;
    std::size_t lcp_chars;

    if (position == begin) {
      best = begin;
      lcp_chars = lcp_with(key, bucket, best);
    } else if (position == end) {
      best = end - 1;
      lcp_chars = lcp_with(key, bucket, best);
    } else {
      const std::size_t predecessor_lcp = lcp_with(key, bucket, position - 1);
      const std::size_t successor_lcp = lcp_with(key, bucket, position);

      if (predecessor_lcp >= successor_lcp) {
        best = position - 1;
        lcp_chars = predecessor_lcp;
      } else {
        best = position;
        lcp_chars = successor_lcp;
      }
    }

    std::size_t ref_pos = L_[best].ref_pos;

    // lcp_chars >= 8 here, so only a flagged bucket can do better.
    if (has_long_short_suffix(bucket)) {
      ++stats_.short_suffix_checks;
      check_short_suffixes(key, ref_pos, lcp_chars);
    }

    result.match_position = static_cast<std::uint32_t>(ref_pos);
    result.match_length = static_cast<std::uint32_t>(lcp_chars);
    return result;
  }

  // Same as PT16RLZParser::lookupTailByKey: `key` is the 1 .. 15
  // character tail packed from the top bit and padded with zero bits; the
  // padded key is looked up as usual and the match capped at `length`.
  KmerLookupResult lookupTailByKey(const std::uint32_t key,
                                   const std::uint32_t length) const {
    const KmerLookupResult padded = lookupKmerByKey(key);

    KmerLookupResult result;
    result.match_position = padded.match_position;
    result.match_length = std::min(padded.match_length, length);
    return result;
  }

  // Fewer than 16 characters left: pack the padded tail and use
  // lookupTailByKey.
  KmerLookupResult lookupKmerOrTail(const input_type& input,
                                    const std::size_t input_pos) const {
    if (input.size() - input_pos < kmer_length) {
      return lookupTailByKey(
          encode_tail(input, input_pos),
          static_cast<std::uint32_t>(input.size() - input_pos));
    }

    return lookupKmerByKey(encode_16mer(input, static_cast<std::uint32_t>(
                                                   input_pos)));
  }

  const Stats& stats() const { return stats_; }

 private:
  // v2's on-disk L entry.
  struct LowerEntry {
    std::uint16_t low;
    std::uint16_t sa_offset;
  };

  // In memory: the same, plus the entry's first reference position.
  struct Entry {
    std::uint16_t low;
    std::uint16_t sa_offset;
    std::uint32_t ref_pos;
  };

  static_assert(sizeof(LowerEntry) == 4);
  static_assert(sizeof(Entry) == 8);

  const reference_type* ref_ = nullptr;
  const suffix_array_type* sa_ = nullptr;

  std::vector<std::uint32_t> H_;
  std::vector<std::uint32_t> H_sa_;
  std::vector<Entry> L_;
  std::unordered_map<std::uint32_t, std::uint32_t> large_offsets_;

  // short_suffix_keys_[L], L = 1 .. short_suffix_count_: the reference's
  // last L characters packed like a 16-mer key, unused low bits 0.
  std::array<std::uint32_t, KMER_LENGTH> short_suffix_keys_{};
  std::size_t short_suffix_count_ = 0;

  // Bit b set: some short suffix of length >= 9 starts with bucket b's 8
  // characters, so a query in bucket b may match it for more than 8.
  std::array<std::uint64_t, NUMBER_OF_BUCKETS / 64> long_short_buckets_{};

  // Precomputed answer per empty bucket (unused for non-empty ones).
  std::vector<std::uint32_t> empty_ref_pos_;
  std::vector<std::uint8_t> empty_length_;

  mutable Stats stats_;

  // ---------- Loading (same file format as PT16RLZParser::load_hl) ----------

  template <typename T>
  static void read_value(std::ifstream& input, T& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(T));

    if (!input) {
      throw std::runtime_error("Failed while reading PT16 H/L table");
    }
  }

  void load_hl(const std::string& path) {
    std::ifstream input(path, std::ios::binary);

    if (!input) {
      throw std::runtime_error("Cannot open PT16 H/L table: " + path);
    }

    char magic[8]{};
    input.read(magic, sizeof(magic));

    const char expected_magic[8] = {'P', 'T', '1', '6', 'H', 'L', '0', '1'};

    if (!input || std::memcmp(magic, expected_magic, sizeof(magic)) != 0) {
      throw std::runtime_error("Invalid PT16 H/L table");
    }

    std::uint64_t entry_count;
    read_value(input, entry_count);

    H_.resize(number_of_buckets + 1);
    H_sa_.resize(number_of_buckets);
    std::vector<LowerEntry> lower(entry_count);

    input.read(reinterpret_cast<char*>(H_.data()),
               static_cast<std::streamsize>(H_.size() * sizeof(std::uint32_t)));
    input.read(
        reinterpret_cast<char*>(H_sa_.data()),
        static_cast<std::streamsize>(H_sa_.size() * sizeof(std::uint32_t)));
    input.read(reinterpret_cast<char*>(lower.data()),
               static_cast<std::streamsize>(lower.size() * sizeof(LowerEntry)));

    for (std::uint32_t position = 0; position < lower.size(); ++position) {
      if (lower[position].sa_offset == large_offset_flag) {
        std::uint32_t offset;
        read_value(input, offset);
        large_offsets_.emplace(position, offset);
      }
    }

    if (!input) {
      throw std::runtime_error("Failed while loading PT16 H/L table");
    }

    if (H_.back() != lower.size()) {
      throw std::runtime_error("PT16 H directory does not end at m");
    }

    // Widen every entry with its first reference position. Entries are in
    // SA order, so this walks the suffix array front to back.
    L_.resize(lower.size());

    for (std::uint32_t position = 0; position < lower.size(); ++position) {
      L_[position].low = lower[position].low;
      L_[position].sa_offset = lower[position].sa_offset;
    }

    for (std::uint32_t bucket = 0; bucket < number_of_buckets; ++bucket) {
      if (H_[bucket] & empty_bucket_flag) {
        continue;
      }

      const std::uint32_t end = H_[next_nonempty_bucket(bucket)];

      for (std::uint32_t position = H_[bucket]; position < end; ++position) {
        L_[position].ref_pos =
            static_cast<std::uint32_t>((*sa_)[sa_start_at(bucket, position)]);
      }
    }
  }

  // ---------- Short suffixes ----------

  void build_short_suffixes() {
    short_suffix_count_ =
        std::min(static_cast<std::size_t>(kmer_length - 1), ref_->size());

    for (std::size_t length = 1; length <= short_suffix_count_; ++length) {
      const std::size_t start = ref_->size() - length;
      std::uint32_t key = 0;

      for (std::size_t j = 0; j < length; ++j) {
        key = (key << 2U) |
              alphatab[static_cast<unsigned char>((*ref_)[start + j])];
      }

      short_suffix_keys_[length] =
          key << (32U - 2U * static_cast<std::uint32_t>(length));

      if (length > 8) {
        const std::uint32_t bucket = short_suffix_keys_[length] >> low_bits;
        long_short_buckets_[bucket / 64] |= std::uint64_t{1} << (bucket % 64);
      }
    }
  }

  bool has_long_short_suffix(const std::uint32_t bucket) const {
    return (long_short_buckets_[bucket / 64] >> (bucket % 64)) & 1U;
  }

  // Same as PT16RLZParser::check_short_suffixes_empty_bucket: longest first, stopping
  // once no remaining suffix could be longer than match_length.
  void check_short_suffixes(const std::uint32_t key, std::size_t& ref_pos,
                            std::size_t& match_length) const {
    for (std::size_t length = short_suffix_count_; length > match_length;
         --length) {
      const std::size_t shared = std::min<std::size_t>(
          length, static_cast<std::size_t>(
                      std::countl_zero(key ^ short_suffix_keys_[length])) /
                      2);

      if (shared > match_length) {
        match_length = shared;
        ref_pos = ref_->size() - length;
      }
    }
  }

  // ---------- Empty-bucket answers ----------

  // For every empty bucket: the v2 answer (LCP with the nearest non-empty
  // bucket, first position of that bucket), raised by any short suffix.
  // A short suffix's match is capped at 8 here because only the bucket's
  // 8 characters are known; if a suffix could go past 8, the bucket is
  // flagged in long_short_buckets_ and the query finishes the check.
  void build_empty_answers() {
    empty_ref_pos_.assign(number_of_buckets, 0);
    empty_length_.assign(number_of_buckets, 0);

    for (std::uint32_t bucket = 0; bucket < number_of_buckets; ++bucket) {
      if (!(H_[bucket] & empty_bucket_flag)) {
        continue;
      }

      const std::uint32_t matching_bucket = H_[bucket] & empty_bucket_mask;
      const std::uint32_t difference = (bucket ^ matching_bucket) << 16U;

      std::size_t match_length = std::countl_zero(difference) / 2;
      std::size_t ref_pos = L_[H_[matching_bucket]].ref_pos;

      const std::uint32_t bucket_key = bucket << low_bits;

      for (std::size_t length = short_suffix_count_; length > match_length;
           --length) {
        const std::size_t shared = std::min<std::size_t>(
            {length, 8,
             static_cast<std::size_t>(
                 std::countl_zero(bucket_key ^ short_suffix_keys_[length])) /
                 2});

        if (shared > match_length) {
          match_length = shared;
          ref_pos = ref_->size() - length;
        }
      }

      empty_ref_pos_[bucket] = static_cast<std::uint32_t>(ref_pos);
      empty_length_[bucket] = static_cast<std::uint8_t>(match_length);
    }
  }

  // ---------- Table helpers (same as PT16RLZParser) ----------

  std::uint32_t next_nonempty_bucket(std::uint32_t bucket) const {
    ++bucket;

    while (bucket < number_of_buckets && (H_[bucket] & empty_bucket_flag)) {
      ++bucket;
    }

    return bucket;
  }

  std::size_t lcp_with(const std::uint32_t key, const std::uint32_t bucket,
                       const std::uint32_t position) const {
    const std::uint32_t candidate =
        (bucket << low_bits) | static_cast<std::uint32_t>(L_[position].low);
    return static_cast<std::size_t>(std::countl_zero(key ^ candidate)) / 2;
  }

  std::uint32_t sa_offset(const std::uint32_t position) const {
    const std::uint16_t offset = L_[position].sa_offset;

    if (offset != large_offset_flag) {
      return offset;
    }

    return large_offsets_.at(position);
  }

  std::uint32_t sa_start_at(const std::uint32_t bucket,
                            const std::uint32_t position) const {
    return H_sa_[bucket] + sa_offset(position);
  }

  std::uint32_t interval_end(const std::uint32_t bucket,
                             const std::uint32_t position,
                             const std::uint32_t end) const {
    std::size_t next_start;

    if (position + 1 < end) {
      next_start = sa_start_at(bucket, position + 1);
    } else {
      const std::uint32_t next_bucket = next_nonempty_bucket(bucket);

      if (next_bucket < number_of_buckets) {
        next_start = H_sa_[next_bucket];
      } else {
        next_start = sa_->size();
      }
    }

    std::size_t sa_end = next_start - 1;

    // Suffixes shorter than 16 have no entry but can sit between two
    // intervals in the SA; drop them from this interval's end.
    while (static_cast<std::size_t>((*sa_)[sa_end]) + kmer_length >
           ref_->size()) {
      --sa_end;
    }

    return static_cast<std::uint32_t>(sa_end);
  }
};
