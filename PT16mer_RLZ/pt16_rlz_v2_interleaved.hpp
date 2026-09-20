#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "parser.hpp"

template <typename T1, typename T2>
class PT16RLZParser {
 public:
  using input_type = std::vector<T1>;
  using reference_type = std::vector<T1>;
  using suffix_array_type = std::vector<T2>;

  using factor_type = std::tuple<std::size_t, std::size_t>;
  using phrase_type = std::tuple<std::size_t, std::size_t, std::size_t>;
  using phrase_vector_type = std::vector<phrase_type>;

  static constexpr std::uint32_t kmer_length = 16;
  static constexpr std::uint32_t bucket_size = 65536;
  static constexpr std::uint32_t low_bits = 16;
  static constexpr std::uint32_t low_mask = bucket_size - 1;
  static constexpr std::uint32_t number_of_buckets = 65536;
  static constexpr std::uint32_t empty_bucket_flag = 1U << 31;
  static constexpr std::uint32_t empty_bucket_mask = empty_bucket_flag - 1;
  static constexpr std::uint32_t binary_search_threshold = 64;
  static constexpr std::uint16_t large_offset_flag = 65535;

  struct Stats {
    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t singleton_hits = 0;
    std::size_t range_hits = 0;
    std::size_t entries = 0;
    std::size_t approx_bytes = 0;
  };

 private:
  struct LookupResult {
    bool found = false;
    std::uint32_t sa_start = 0;
    std::uint32_t sa_end = 0;
    std::size_t ref_pos = 0;
    std::size_t match_length = 0;
  };

  struct LowerEntry {
    std::uint16_t low;        // current lower 16 bits of the 16-mer
    std::uint16_t sa_offset;  // relative SA interval start
  };
  struct ShortSuffix {
    std::uint32_t bucket;
    std::size_t ref_pos;
    std::size_t length;
  };

  std::vector<ShortSuffix> short_suffixes_;

  const reference_type* ref_ = nullptr;
  const suffix_array_type* sa_ = nullptr;

  std::array<std::uint8_t, 256> alphatab_{};

  // ---------- PT16 H/L representation ----------

  // One H array interleaving starting position of each bucket in the lower
  // level and starting SA position corresponding to each bucket.
  std::vector<std::uint32_t> H_interleaved_;

  // Lower-level entries containing the 16-mer value and relative SA offset.
  std::vector<LowerEntry> L_;

  // Full relative SA offsets for entries marked with the 65535 escape value.
  std::unordered_map<std::uint32_t, std::uint32_t> large_offsets_;

  mutable Stats stats_;

 public:
  PT16RLZParser(const reference_type& ref, const suffix_array_type& sa,
                const std::string& pt16_path)
      : ref_(&ref), sa_(&sa) {
    initialise_alphatab();
    load_hl(pt16_path);
    build_short_suffixes();
  }

  factor_type computeLZFactorAt(const input_type& input,
                                const std::size_t input_pos) {
    // Fewer than 16 characters remain: use ordinary RLZ.
    if (input.size() - input_pos < kmer_length) {
      return ::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);
    }

    // Pack the next 16 characters and select the H bucket.
    const std::uint32_t key = pack_16mer(input, input_pos);
    const std::uint32_t bucket = key >> low_bits;

    // Empty bucket: compute the short factor directly.
    if (H_interleaved_[2 * bucket] & empty_bucket_flag) {
      // Lower 31 bits store the non-empty bucket with maximum LCP.
      const std::uint32_t matching_bucket =
          H_interleaved_[2 * bucket] & empty_bucket_mask;

      // Compute the common prefix length in bits, then DNA characters.
      const std::uint32_t difference = (bucket ^ matching_bucket) << 16U;
      const std::size_t lcp_bits = std::countl_zero(difference);
      const std::size_t lcp_chars = lcp_bits / 2;

      // Take a reference position from the matching bucket.
      const std::size_t interval_position = H_interleaved_[2 * matching_bucket];
      // Recover the first SA interval start of the matching bucket.
      const std::size_t sa_position =
          sa_start_at(matching_bucket, interval_position);

      std::size_t ref_pos = static_cast<std::size_t>((*sa_)[sa_position]);
      std::size_t match_length = lcp_chars;

      check_short_suffixes(input, input_pos, bucket, ref_pos, match_length);
      ++stats_.misses;

      return {ref_pos, match_length};
    }

    // Non-empty bucket: use the existing PT16 lookup.
    const LookupResult result = lookup(input, input_pos, key);

    // Non-empty bucket miss: lookup() has already computed the short factor.
    if (!result.found) {
      return {result.ref_pos, result.match_length};
    }

    if (result.sa_start == result.sa_end) {
      ++stats_.singleton_hits;
    } else {
      ++stats_.range_hits;
    }

    std::size_t offset = kmer_length;
    std::size_t j = input_pos + kmer_length;
    std::size_t nlb = result.sa_start;
    std::size_t nrb = result.sa_end;

    // Range case: narrow the SA interval from character 17 onward.
    while (nlb < nrb && j < input.size()) {
      const auto lb =
          ::binarySearchLB<T1, T2>(*ref_, *sa_, nlb, nrb, offset, input[j]);

      if (!lb) {
        break;
      }

      const auto rb = ::binarySearchRB<T1, T2>(
          *ref_, *sa_, static_cast<std::size_t>(lb.value()), nrb, offset,
          input[j]);

      if (!rb) {
        break;
      }

      nlb = static_cast<std::size_t>(lb.value());
      nrb = static_cast<std::size_t>(rb.value());

      ++j;
      ++offset;
    }

    std::size_t match = static_cast<std::size_t>((*sa_)[nlb]);

    // Singleton case: extend directly from character 17 onward.
    if (nlb == nrb) {
      while (j < input.size() && match + offset < ref_->size() &&
             (*ref_)[match + offset] == input[j]) {
        ++j;
        ++offset;
      }
    }

    return {match, offset};
  }

  phrase_vector_type lzFactorize(const input_type& input) {
    phrase_vector_type spl_vec;
    std::size_t i = 0;

    while (i < input.size()) {
      auto [pos, len] = computeLZFactorAt(input, i);

      if (len <= 1) {
        pos = static_cast<std::size_t>(input[i]);
        len = 1;
      }

      spl_vec.push_back({i, pos, len});
      i += len;
    }

    return spl_vec;
  }

  const Stats& stats() const { return stats_; }

 private:
  // ---------- Build alphatab once ----------

  void initialise_alphatab() {
    alphatab_[static_cast<unsigned char>('A')] = 0;  // 00
    alphatab_[static_cast<unsigned char>('C')] = 1;  // 01
    alphatab_[static_cast<unsigned char>('G')] = 2;  // 10
    alphatab_[static_cast<unsigned char>('T')] = 3;  // 11
  }

  void build_short_suffixes() {
    for (std::size_t length = 8; length < kmer_length; ++length) {
      const std::size_t ref_pos = ref_->size() - length;

      std::uint32_t bucket = 0;

      for (std::size_t j = 0; j < 8; ++j) {
        const std::uint8_t code =
            alphatab_[static_cast<unsigned char>((*ref_)[ref_pos + j])];

        bucket = (bucket << 2U) | code;
      }

      short_suffixes_.push_back({bucket, ref_pos, length});
    }
  }

  void check_short_suffixes(const input_type& input,
                            const std::size_t input_pos,
                            const std::uint32_t bucket, std::size_t& ref_pos,
                            std::size_t& match_length) const {
    for (const ShortSuffix& suffix : short_suffixes_) {
      if (suffix.bucket != bucket) {
        continue;
      }

      std::size_t length = 8;

      while (length < suffix.length &&
             input[input_pos + length] == (*ref_)[suffix.ref_pos + length]) {
        ++length;
      }

      if (length > match_length) {
        match_length = length;
        ref_pos = suffix.ref_pos;
      }
    }
  }

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

    // local arrays just for reading the input
    std::vector<std::uint32_t> H_;
    std::vector<std::uint32_t> H_sa_;

    H_.resize(number_of_buckets + 1);
    H_sa_.resize(number_of_buckets);
    H_interleaved_.resize(number_of_buckets * 2 + 1);
    L_.resize(entry_count);

    input.read(reinterpret_cast<char*>(H_.data()),
               static_cast<std::streamsize>(H_.size() * sizeof(std::uint32_t)));

    input.read(
        reinterpret_cast<char*>(H_sa_.data()),
        static_cast<std::streamsize>(H_sa_.size() * sizeof(std::uint32_t)));

    // interleave the H arrays
    for (size_t i = 0; i < number_of_buckets; ++i) {
      H_interleaved_.push_back(H_[i]);
      H_interleaved_.push_back(H_sa_[i]);
    }
    H_interleaved_.push_back(H_[number_of_buckets]);

    input.read(reinterpret_cast<char*>(L_.data()),
               static_cast<std::streamsize>(L_.size() * sizeof(LowerEntry)));

    // Read the full offsets for lower-level entries marked with 65535.
    for (std::uint32_t position = 0; position < L_.size(); ++position) {
      if (L_[position].sa_offset == large_offset_flag) {
        std::uint32_t offset;
        read_value(input, offset);
        large_offsets_.emplace(position, offset);
      }
    }

    if (!input) {
      throw std::runtime_error("Failed while loading PT16 H/L table");
    }

    if (H_.back() != L_.size()) {
      throw std::runtime_error("PT16 H directory does not end at m");
    }

    stats_.entries = L_.size();

    stats_.approx_bytes =
        H_.size() * sizeof(std::uint32_t) +
        H_sa_.size() * sizeof(std::uint32_t) + L_.size() * sizeof(LowerEntry) +
        large_offsets_.size() *
            sizeof(std::pair<const std::uint32_t, std::uint32_t>);
  }

  std::uint32_t pack_16mer(const input_type& input,
                           const std::size_t position) const {
    std::uint32_t key = 0;

    for (std::uint32_t j = 0; j < kmer_length; ++j) {
      const std::uint8_t code =
          alphatab_[static_cast<unsigned char>(input[position + j])];
      key = (key << 2U) | code;
    }

    return key;
  }

  // Recover the inclusive SA interval end from the next interval start.
  // For the final entry in a bucket, use the next non-empty bucket's SA start.
  std::uint32_t interval_end(const std::uint32_t bucket,
                             const std::uint32_t position,
                             const std::uint32_t end) const {
    std::size_t next_start;

    if (position + 1 < end) {
      next_start = sa_start_at(bucket, position + 1);
    } else {
      std::uint32_t next_bucket = bucket + 1;

      while (next_bucket < number_of_buckets &&
             (H_interleaved_[2 * next_bucket] & empty_bucket_flag)) {
        ++next_bucket;
      }

      if (next_bucket < number_of_buckets) {
        next_start = H_interleaved_[2 * next_bucket + 1];
      } else {
        next_start = sa_->size();
      }
    }

    std::size_t sa_end = next_start - 1;

    /*
    build_entries() skips suffixes shorter than 16. Such suffixes can
    occur between two valid PT16 intervals, so remove them from the
    calculated end of the current interval.
    */
    while (static_cast<std::size_t>((*sa_)[sa_end]) + kmer_length >
           ref_->size()) {
      --sa_end;
    }

    return static_cast<std::uint32_t>(sa_end);
  }

  // Return the relative SA offset for a lower-level entry.
  // The value 65535 indicates that the full offset is stored in large_offsets_.

  std::uint32_t sa_offset(const std::uint32_t position) const {
    const std::uint16_t offset = L_[position].sa_offset;

    if (offset != large_offset_flag) {
      return offset;
    }

    return large_offsets_.at(position);
  }

  // Recover the absolute SA interval start from the bucket's SA start
  // and the relative offset stored for the lower-level entry.
  std::uint32_t sa_start_at(const std::uint32_t bucket,
                            const std::uint32_t position) const {
    return H_interleaved_[2 * bucket + 1] + sa_offset(position);
  }

  LookupResult lookup(const input_type& input, const std::size_t input_pos,
                      const std::uint32_t key) const {
    const std::uint32_t bucket = key >> low_bits;
    const std::uint16_t low = static_cast<std::uint16_t>(key & low_mask);

    const std::uint32_t begin = H_interleaved_[2 * bucket];

    // Find the next non-empty H entry to obtain the end of this L bucket.
    std::uint32_t next_bucket = bucket + 1;

    while (next_bucket < number_of_buckets &&
           (H_interleaved_[2 * next_bucket] & empty_bucket_flag)) {
      ++next_bucket;
    }

    const std::uint32_t end = H_interleaved_[2 * next_bucket];
    const std::uint32_t bucket_entries = end - begin;

    std::uint32_t insertion_position = begin;

    // Small bucket: linear scan.

    if (bucket_entries < binary_search_threshold) {
      while (insertion_position < end && L_[insertion_position].low < low) {
        ++insertion_position;
      }

    }

    // Large bucket: binary search.
    else {
      const auto it = std::lower_bound(
          L_.begin() + begin, L_.begin() + end, low,
          [](const LowerEntry& entry, const std::uint16_t value) {
            return entry.low < value;
          });

      insertion_position = static_cast<std::uint32_t>(it - L_.begin());
    }

    // Exact 16-mer hit.
    if (insertion_position < end && L_[insertion_position].low == low) {
      ++stats_.hits;

      // Recover the SA interval from the bucket start and relative offsets.
      const std::uint32_t sa_start = sa_start_at(bucket, insertion_position);

      const std::uint32_t sa_end =
          interval_end(bucket, insertion_position, end);

      return {true, sa_start, sa_end, 0, 0};
    }

    ++stats_.misses;

    // Non-empty bucket miss:
    // choose the neighbouring key with the longest common prefix.
    std::uint32_t best_position;
    std::size_t lcp_chars;

    // Key is smaller than everything in the bucket.
    if (insertion_position == begin) {
      best_position = begin;

      const std::uint32_t candidate_key =
          (bucket << low_bits) |
          static_cast<std::uint32_t>(L_[best_position].low);

      lcp_chars = std::countl_zero(key ^ candidate_key) / 2;
    }

    // Key is larger than everything in the bucket.
    else if (insertion_position == end) {
      best_position = end - 1;

      const std::uint32_t candidate_key =
          (bucket << low_bits) |
          static_cast<std::uint32_t>(L_[best_position].low);

      lcp_chars = std::countl_zero(key ^ candidate_key) / 2;
    }

    // Key lies between two entries.
    else {
      const std::uint32_t predecessor = insertion_position - 1;
      const std::uint32_t successor = insertion_position;

      const std::uint32_t predecessor_key =
          (bucket << low_bits) |
          static_cast<std::uint32_t>(L_[predecessor].low);

      const std::uint32_t successor_key =
          (bucket << low_bits) | static_cast<std::uint32_t>(L_[successor].low);

      const std::size_t predecessor_lcp =
          std::countl_zero(key ^ predecessor_key) / 2;

      const std::size_t successor_lcp =
          std::countl_zero(key ^ successor_key) / 2;

      if (predecessor_lcp >= successor_lcp) {
        best_position = predecessor;
        lcp_chars = predecessor_lcp;
      } else {
        best_position = successor;
        lcp_chars = successor_lcp;
      }
    }

    // Recover the SA start of the selected neighbouring 16-mer.
    const std::size_t sa_position = sa_start_at(bucket, best_position);

    std::size_t ref_pos = static_cast<std::size_t>((*sa_)[sa_position]);

    check_short_suffixes(input, input_pos, bucket, ref_pos, lcp_chars);

    return {false, 0, 0, ref_pos, lcp_chars};
  }
};
