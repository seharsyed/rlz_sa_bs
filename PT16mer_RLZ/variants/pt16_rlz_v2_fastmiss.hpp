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

#include "../pt16_utils.hpp"  // KmerLookupResult, EntryComposition, constants
#include "pt16_short_suffixes.hpp"  // collect_short_suffixes, ShortSuffixIndex
#include "../rlz_common.hpp"

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

    // Finger lookups (lookupKmerByKey(key, finger)): how many started a
    // fresh bucket search (a new bucket, or a key below the previous one),
    // how many continued from the previous insertion point in the same
    // bucket, and how many L entries those continuations stepped over in
    // total.
    std::size_t finger_restarts = 0;
    std::size_t finger_continues = 0;
    std::size_t finger_steps = 0;
  };

  /**
   * A finger into the table, for a run of lookups with non-decreasing keys
   * (e.g. an input's 16-mers sorted by key): see lookupKmerByKey(key,
   * finger). Start each run with a fresh Finger; it is only a position in
   * this table, so it is cheap to copy and holds nothing to release.
   */
  class Finger {
   public:
    Finger() = default;

   private:
    friend class PT16FastMissParser;

    // The previous lookup's bucket (number_of_buckets: none yet), whether
    // it is empty, its entries [begin, end) in L when it is not, and the
    // previous key's low part and insertion point there.
    std::uint32_t bucket = NUMBER_OF_BUCKETS;
    bool empty = false;
    std::uint32_t begin = 0;
    std::uint32_t end = 0;
    std::uint32_t at = 0;
    std::uint16_t low = 0;
  };

  PT16FastMissParser(const reference_type& ref, const suffix_array_type& sa,
                     const std::string& pt16_path)
      : ref_(&ref), sa_(&sa) {
    load_hl(pt16_path);
    short_index_ = ShortSuffixIndex(collect_short_suffixes(*ref_));
    build_empty_answers();
    build_trimmed_ends();

    stats_.entries = L_.size();
    stats_.approx_bytes =
        H_.size() * sizeof(std::uint32_t) +
        H_sa_.size() * sizeof(std::uint32_t) + L_.size() * sizeof(Entry) +
        large_offsets_.size() *
            sizeof(std::pair<const std::uint32_t, std::uint32_t>) +
        empty_ref_pos_.size() * sizeof(std::uint32_t) +
        empty_length_.size() * sizeof(std::uint8_t) +
        short_index_.memory_bytes() +
        trimmed_.size() * sizeof(std::uint64_t) +
        trimmed_end_.size() *
            sizeof(std::pair<const std::uint32_t, std::uint32_t>);
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

    if (h & empty_bucket_flag) {
      return empty_bucket_result(key, bucket);
    }

    const std::uint16_t low = static_cast<std::uint16_t>(key & low_mask);
    const std::uint32_t begin = h;
    const std::uint32_t end = H_[next_nonempty_bucket(bucket)];

    return bucket_result(key, bucket, begin, end,
                         lower_bound_low(begin, end, low));
  }

  /**
   * Same result as lookupKmerByKey(key), for a run of lookups whose keys
   * never decrease (the input's 16-mers in sorted order): `finger`
   * remembers where the previous lookup landed, so a lookup does not
   * search its bucket from scratch.
   *
   * In the previous lookup's bucket, the search walks forward in L from
   * the previous insertion point -- the key is not smaller, so its
   * insertion point is not earlier -- and ends at this key's insertion
   * point, which becomes the next start. Only a key in a new bucket
   * restarts: its bucket's range is looked up in H and searched as by
   * lookupKmerByKey(key) (linearly or by binary search, by size). An empty
   * bucket is answered as by lookupKmerByKey(key).
   *
   * A key below the previous one also restarts, so the result is correct
   * in any order; it is only fast in sorted order. (Same scheme as
   * PT16SassyLookup::lookup(key, finger).)
   */
  KmerLookupResult lookupKmerByKey(const std::uint32_t key,
                                   Finger& finger) const {
    const std::uint32_t bucket = key >> low_bits;
    const std::uint16_t low = static_cast<std::uint16_t>(key & low_mask);

    if (bucket == finger.bucket && low >= finger.low) {
      finger.low = low;

      if (finger.empty) {
        return empty_bucket_result(key, bucket);
      }

      ++stats_.finger_continues;

      std::uint32_t at = finger.at;

      while (at < finger.end && L_[at].low < low) {
        ++at;
      }

      stats_.finger_steps += at - finger.at;
      finger.at = at;

      return bucket_result(key, bucket, finger.begin, finger.end, at);
    }

    // ---------- Restart: a new bucket (or a smaller key) ----------

    ++stats_.finger_restarts;
    finger.bucket = bucket;
    finger.low = low;
    finger.empty = (H_[bucket] & empty_bucket_flag) != 0;

    if (finger.empty) {
      return empty_bucket_result(key, bucket);
    }

    finger.begin = H_[bucket];
    finger.end = H_[next_nonempty_bucket(bucket)];
    finger.at = lower_bound_low(finger.begin, finger.end, low);

    return bucket_result(key, bucket, finger.begin, finger.end, finger.at);
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

  /**
   * The longest prefix of input[input_pos..] that occurs in the reference,
   * as (reference position, length): lookupKmerOrTail, then extended past
   * the 16-mer on a hit, as PT16SassyLookup::find_longest_matching_factor
   * does:
   *
   *   - a singleton hit extends character by character from the entry's
   *     own ref_pos (no suffix array access);
   *   - a range hit narrows within its SA slice (hit.positions) with
   *     rlz::binarySearchLB/RB, then extends once one occurrence is left;
   *   - a miss or a tail is already the longest match.
   */
  std::pair<std::size_t, std::size_t> findLongestMatchingFactor(
      const input_type& input, const std::size_t input_pos) const {
    const KmerLookupResult hit = lookupKmerOrTail(input, input_pos);

    if (!hit.found) {
      return {hit.match_position, hit.match_length};
    }

    std::size_t offset = kmer_length;
    std::size_t j = input_pos + kmer_length;
    std::size_t match = hit.match_position;

    if (hit.count > 1) {
      std::int64_t nlb = 0;
      std::int64_t nrb = static_cast<std::int64_t>(hit.positions.size()) - 1;

      while (nlb < nrb && j < input.size()) {
        const auto lb = rlz::binarySearchLB(*ref_, hit.positions, nlb, nrb,
                                            static_cast<std::int64_t>(offset),
                                            input[j]);

        if (!lb) {
          break;
        }

        const auto rb =
            rlz::binarySearchRB(*ref_, hit.positions, lb.value(), nrb,
                                static_cast<std::int64_t>(offset), input[j]);

        if (!rb) {
          break;
        }

        nlb = lb.value();
        nrb = rb.value();
        ++j;
        ++offset;
      }

      match = hit.positions[static_cast<std::size_t>(nlb)];

      // Still several occurrences when the input ran out: any is correct.
      if (nlb != nrb) {
        return {match, offset};
      }
    }

    while (j < input.size() && match + offset < ref_->size() &&
           (*ref_)[match + offset] == input[j]) {
      ++j;
      ++offset;
    }

    return {match, offset};
  }

  // The same greedy parse as PT16RLZParser::lzFactorize and
  // PT16SassyLookup::lzFactorize: a match of 0 or 1 characters becomes a
  // literal factor.
  Triples lzFactorize(const input_type& input) const {
    Triples factors;
    std::size_t i = 0;

    while (i < input.size()) {
      auto [pos, len] = findLongestMatchingFactor(input, i);

      if (len <= 1) {
        pos = static_cast<std::size_t>(input[i]);
        len = 1;
      }

      factors.push_back({i, pos, len});
      i += len;
    }

    return factors;
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

  // Every short suffix of the reference (see pt16_short_suffixes.hpp): the
  // positions within 15 characters of its end or of a separator, arranged
  // for the empty-bucket precomputation and for the query-time check of
  // suffixes longer than 8 in the query's bucket.
  ShortSuffixIndex short_index_;

  // An entry's SA interval ends where the next entry's begins, except when
  // suffixes without an entry (too short, or containing a separator) lie
  // in between: bit `position` of trimmed_ is set for such an entry, and
  // trimmed_end_ holds its real (inclusive) end. Worked out once at load
  // time, so a lookup never walks back over them -- a separator run in the
  // reference can put millions there.
  std::vector<std::uint64_t> trimmed_;
  std::unordered_map<std::uint32_t, std::uint32_t> trimmed_end_;

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

  // ---------- Interval ends ----------

  // The last SA index before the next entry's interval: the entry's end,
  // unless suffixes without an entry lie in between (see trimmed_).
  std::uint32_t raw_interval_end(const std::uint32_t bucket,
                                 const std::uint32_t position,
                                 const std::uint32_t end) const {
    if (position + 1 < end) {
      return sa_start_at(bucket, position + 1) - 1;
    }

    const std::uint32_t next_bucket = next_nonempty_bucket(bucket);

    return static_cast<std::uint32_t>(
        (next_bucket < number_of_buckets ? H_sa_[next_bucket] : sa_->size()) -
        1);
  }

  void build_trimmed_ends() {
    trimmed_.assign((L_.size() + 63) / 64, 0);

    for (std::uint32_t bucket = 0; bucket < number_of_buckets; ++bucket) {
      if (H_[bucket] & empty_bucket_flag) {
        continue;
      }

      const std::uint32_t end = H_[next_nonempty_bucket(bucket)];

      for (std::uint32_t position = H_[bucket]; position < end; ++position) {
        std::uint32_t sa_end = raw_interval_end(bucket, position, end);

        if (window_is_acgt(*ref_, (*sa_)[sa_end])) {
          continue;
        }

        // Every suffix with a 16-mer between this entry and the next
        // belongs to this entry, so the ones without are all at the end.
        while (!window_is_acgt(*ref_, (*sa_)[sa_end])) {
          --sa_end;
        }

        trimmed_[position / 64] |= std::uint64_t{1} << (position % 64);
        trimmed_end_.emplace(position, sa_end);
      }
    }
  }

  // ---------- Empty-bucket answers ----------

  // For every empty bucket: the v2 answer (LCP with the nearest non-empty
  // bucket, first position of that bucket), raised by any short suffix.
  // A short suffix's match is capped at 8 here because only the bucket's
  // 8 characters are known; if a suffix could go past 8, the bucket is
  // flagged (ShortSuffixIndex::has_long) and the query finishes the check.
  void build_empty_answers() {
    empty_ref_pos_.assign(number_of_buckets, 0);
    empty_length_.assign(number_of_buckets, 0);

    for (std::uint32_t bucket = 0; bucket < number_of_buckets; ++bucket) {
      if (!(H_[bucket] & empty_bucket_flag)) {
        continue;
      }

      const std::uint32_t matching_bucket = H_[bucket] & empty_bucket_mask;

      // No non-empty bucket at all (a table without entries): only the
      // short suffixes can match.
      std::size_t match_length = 0;
      std::size_t ref_pos = 0;

      if (matching_bucket < number_of_buckets) {
        const std::uint32_t difference = (bucket ^ matching_bucket) << 16U;
        match_length = std::countl_zero(difference) / 2;
        ref_pos = L_[H_[matching_bucket]].ref_pos;
      }

      const auto [suffix_length, suffix_position] =
          short_index_.best_in_bucket(bucket);

      if (suffix_length > match_length) {
        match_length = suffix_length;
        ref_pos = suffix_position;
      }

      empty_ref_pos_[bucket] = static_cast<std::uint32_t>(ref_pos);
      empty_length_[bucket] = static_cast<std::uint8_t>(match_length);
    }
  }

  // ---------- Table helpers (same as PT16RLZParser) ----------

  // ---------- Building a lookup's result ----------

  // An empty bucket: the precomputed answer, finished by the full
  // short-suffix check only in a flagged bucket.
  KmerLookupResult empty_bucket_result(const std::uint32_t key,
                                       const std::uint32_t bucket) const {
    std::size_t ref_pos = empty_ref_pos_[bucket];
    std::size_t match_length = empty_length_[bucket];

    if (short_index_.has_long(bucket)) {
      ++stats_.short_suffix_checks;
      short_index_.raise_long(key, match_length, ref_pos);
    }

    ++stats_.misses;
    ++stats_.empty_bucket_misses;

    KmerLookupResult result;
    result.match_position = static_cast<std::uint32_t>(ref_pos);
    result.match_length = static_cast<std::uint32_t>(match_length);
    return result;
  }

  // Index of the first entry in [begin, end) whose low part is not below
  // `low`, or `end`: linear below binary_search_threshold entries,
  // std::lower_bound at or above it.
  std::uint32_t lower_bound_low(const std::uint32_t begin,
                                const std::uint32_t end,
                                const std::uint16_t low) const {
    if (end - begin < binary_search_threshold) {
      ++stats_.linear_bucket_searches;

      std::uint32_t position = begin;

      while (position < end && L_[position].low < low) {
        ++position;
      }

      return position;
    }

    ++stats_.binary_bucket_searches;

    const auto it = std::lower_bound(
        L_.begin() + begin, L_.begin() + end, low,
        [](const Entry& entry, const std::uint16_t value) {
          return entry.low < value;
        });

    return static_cast<std::uint32_t>(it - L_.begin());
  }

  // A non-empty bucket, with entries [begin, end) in L, and `position` the
  // key's insertion point there: a hit if that entry is the key, otherwise
  // a miss answered by the neighbouring entry with the longest common
  // prefix.
  KmerLookupResult bucket_result(const std::uint32_t key,
                                 const std::uint32_t bucket,
                                 const std::uint32_t begin,
                                 const std::uint32_t end,
                                 const std::uint32_t position) const {
    const std::uint16_t low = static_cast<std::uint16_t>(key & low_mask);

    KmerLookupResult result;

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
    if (short_index_.has_long(bucket)) {
      ++stats_.short_suffix_checks;
      short_index_.raise_long(key, lcp_chars, ref_pos);
    }

    result.match_position = static_cast<std::uint32_t>(ref_pos);
    result.match_length = static_cast<std::uint32_t>(lcp_chars);
    return result;
  }

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
    if ((trimmed_[position / 64] >> (position % 64)) & 1U) {
      return trimmed_end_.at(position);
    }

    return raw_interval_end(bucket, position, end);
  }
};
