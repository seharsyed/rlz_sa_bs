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
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../pt16_utils.hpp"  // EntryComposition, KmerLookupResult
#include "../rlz_common.hpp"

template <typename T1, typename T2>
class PT16RLZParser {
 public:
  using input_type = std::vector<T1>;
  using reference_type = std::vector<T1>;
  using suffix_array_type = std::vector<T2>;

  using factor_type = std::tuple<std::size_t, std::size_t>;
  using phrase_type = std::tuple<std::size_t, std::size_t, std::size_t>;
  using phrase_vector_type = std::vector<phrase_type>;
  using factor_vector_type = std::vector<factor_type>;
  using ms_vector_type = std::vector<std::pair<std::uint32_t, std::uint32_t>>;

  static constexpr std::uint32_t kmer_length = 16;
  static constexpr std::uint32_t bucket_size = 65536;
  static constexpr std::uint32_t low_bits = 16;
  static constexpr std::uint32_t low_mask = bucket_size - 1;
  static constexpr std::uint32_t number_of_buckets = 65536;
  static constexpr std::uint32_t empty_bucket_flag = 1U << 31;
  static constexpr std::uint32_t empty_bucket_mask = empty_bucket_flag - 1;
  static constexpr std::uint32_t binary_search_threshold =
      BINARY_SEARCH_THRESHOLD;
  static constexpr std::uint16_t large_offset_flag = 65535;

  struct Stats {
    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t singleton_hits = 0;
    std::size_t range_hits = 0;
    std::size_t entries = 0;
    std::size_t approx_bytes = 0;

    // How every non-empty-bucket dispatch searched its bucket: a linear
    // scan below binary_search_threshold entries, std::lower_bound at or
    // above it. Updated by lookup(), so this covers both hits and
    // non-empty-bucket misses (an empty-bucket miss never searches L_ at
    // all, so it touches neither counter).
    std::size_t linear_bucket_searches = 0;
    std::size_t binary_bucket_searches = 0;
  };

  // The shared, format-agnostic type (pt16_utils.hpp) -- PT16SassyLookup's
  // lookup returns the same one, so code consuming either format's
  // results (e.g. a chain-extension pass) never needs to know which table
  // produced them.
  using KmerLookupResult = ::KmerLookupResult;

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
  // short_suffix_keys_[L], for L = 1 .. short_suffix_count_: the
  // reference's last L characters, packed like a 16-mer key (2 bits each
  // from the top bit down, unused low bits 0). These positions are too
  // close to the end to have an H/L entry, so check_short_suffixes compares
  // the query's key against them directly.
  std::array<std::uint32_t, kmer_length> short_suffix_keys_{};
  std::size_t short_suffix_count_ = 0;

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
      return rlz::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);
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

      check_short_suffixes_empty_bucket(key, ref_pos, match_length);
      ++stats_.misses;

      return {ref_pos, match_length};
    }

    // Non-empty bucket: use the existing PT16 lookup.
    const LookupResult result = lookup(key);

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
          rlz::binarySearchLB(*ref_, *sa_, nlb, nrb, offset, input[j]);

      if (!lb) {
        break;
      }

      const auto rb = rlz::binarySearchRB(
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

  // Brute-force matching statistics: one PT16 longest-match query per input
  // position, returned as (reference position, match length) pairs.
  ms_vector_type computeMS_brute(const input_type& input) {
    ms_vector_type ms;
    ms.reserve(input.size());

    for (std::size_t i = 0; i < input.size(); ++i) {
      const auto [pos, len] = computeLZFactorAt(input, i);

      ms.emplace_back(static_cast<std::uint32_t>(pos),
                      static_cast<std::uint32_t>(len));
    }

    return ms;
  }

  /**
   * A 16-mer-only lookup: the same bucket dispatch, bucket search and
   * short-suffix fallback as computeLZFactorAt/lookup, but -- like
   * PT16SassyLookup::lookup, not find_longest_matching_factor -- no SA
   * interval narrowing beyond the raw H/L entry on a hit.
   *
   * This exists to compare the two table formats' pure lookup cost head to
   * head: this format still needs one suffix-array access to turn a hit's
   * SA interval into an actual reference position (sa_start_at, then
   * (*sa_)[...] below), which the self-contained sassy table does not --
   * it decodes a position (or its first occurrence, for a range) directly
   * out of the L entry itself.
   *
   * Requires input.size() - input_pos >= kmer_length; use
   * lookupKmerOrTail for a shorter tail.
   */
  KmerLookupResult lookupKmer(const input_type& input,
                              const std::size_t input_pos) const {
    return lookupKmerByKey(input, input_pos, pack_16mer(input, input_pos));
  }

  /**
   * Same as lookupKmer, but with the key already packed elsewhere -- for a
   * caller that bucketed several positions by table bucket first and
   * cached each one's key, so the input never needs to be re-read to
   * repack it (see PT16ScanMS::computeMatchingStatisticsBucketed in
   * ms_variants.hpp, the counterpart of
   * PT16SassyMS::computeMatchingStatisticsScanOnlyBucketed).
   * `input`/`input_pos` are unused (the short-suffix check works from `key`
   * too); they are kept so this has the same shape as
   * PT16SassyLookup::lookup.
   */
  KmerLookupResult lookupKmerByKey([[maybe_unused]] const input_type& input,
                                   [[maybe_unused]] const std::size_t input_pos,
                                   const std::uint32_t key) const {
    const std::uint32_t bucket = key >> low_bits;

    KmerLookupResult result;

    // ---------- Empty bucket: same short factor as computeLZFactorAt ----------

    if (H_interleaved_[2 * bucket] & empty_bucket_flag) {
      const std::uint32_t matching_bucket =
          H_interleaved_[2 * bucket] & empty_bucket_mask;

      const std::uint32_t difference = (bucket ^ matching_bucket) << 16U;
      const std::size_t lcp_chars = std::countl_zero(difference) / 2;

      const std::size_t interval_position =
          H_interleaved_[2 * matching_bucket];
      const std::size_t sa_position = sa_start_at(
          matching_bucket, static_cast<std::uint32_t>(interval_position));

      std::size_t ref_pos = static_cast<std::size_t>((*sa_)[sa_position]);
      std::size_t match_length = lcp_chars;

      check_short_suffixes_empty_bucket(key, ref_pos, match_length);

      result.match_position = static_cast<std::uint32_t>(ref_pos);
      result.match_length = static_cast<std::uint32_t>(match_length);
      ++stats_.misses;
      return result;
    }

    // ---------- Non-empty bucket: the existing PT16 lookup ----------

    const LookupResult inner = lookup(key);

    if (!inner.found) {
      result.match_position = static_cast<std::uint32_t>(inner.ref_pos);
      result.match_length = static_cast<std::uint32_t>(inner.match_length);
      return result;
    }

    result.found = true;
    result.match_length = kmer_length;
    result.count = inner.sa_end - inner.sa_start + 1;
    result.match_position = static_cast<std::uint32_t>((*sa_)[inner.sa_start]);

    // The occurrences of a range are exactly sa_[sa_start..sa_end]
    // (inclusive) -- already contiguous in the suffix array itself,
    // thanks to lexicographic sorting grouping identical-prefix suffixes
    // together, so exposing them needs no extra storage, unlike sassy's
    // sampled_sa_ (this format has no equivalent scratch array; the
    // suffix array already IS one, for this purpose).
    result.positions = std::span<const std::uint32_t>(
        sa_->data() + inner.sa_start, static_cast<std::size_t>(result.count));

    // lookup() above only counts stats_.hits (generic); the
    // singleton/range split is computeLZFactorAt's job normally, which
    // this bypasses, so it is repeated here -- same condition, same
    // stats_ fields, so PT16SassyLookup::lookup's own singleton/range
    // counters and these stay comparable.
    if (inner.sa_start == inner.sa_end) {
      ++stats_.singleton_hits;
    } else {
      ++stats_.range_hits;
    }

    return result;
  }

  /**
   * lookupKmer for a query with all 16 characters, or the same suffix-array
   * fallback computeLZFactorAt uses for the last few characters of an
   * input: this format has no self-contained short-tail scheme like
   * PT16SassyLookup::lookup_tail's zero-padding, so a tail here bypasses
   * the PT16 table entirely (as it always has) rather than staying inside
   * it. Tails are rare (at most kmer_length - 1 per input, always at the
   * very end), so this does not affect what the bulk lookup comparison is
   * measuring.
   */
  KmerLookupResult lookupKmerOrTail(const input_type& input,
                                    const std::size_t input_pos) const {
    if (input.size() - input_pos < kmer_length) {
      const auto [pos, len] =
          rlz::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);

      KmerLookupResult result;
      result.found = len > 0;
      result.match_position = static_cast<std::uint32_t>(pos);
      result.match_length = static_cast<std::uint32_t>(len);
      return result;
    }

    return lookupKmer(input, input_pos);
  }

  const Stats& stats() const { return stats_; }

  /**
   * Classifies every entry in the table by whether its SA interval has
   * exactly one occurrence (singleton) or more than one (range). A
   * structural property of the table itself, independent of any query --
   * computed once (walks every bucket and every entry, same cost as a
   * single full scan of the table), not per lookup.
   */
  EntryComposition entryComposition() const {
    EntryComposition result;

    for (std::uint32_t bucket = 0; bucket < number_of_buckets; ++bucket) {
      if (H_interleaved_[2 * bucket] & empty_bucket_flag) {
        continue;
      }

      const std::uint32_t begin = H_interleaved_[2 * bucket];

      std::uint32_t next_bucket = bucket + 1;

      while (next_bucket < number_of_buckets &&
             (H_interleaved_[2 * next_bucket] & empty_bucket_flag)) {
        ++next_bucket;
      }

      const std::uint32_t end = H_interleaved_[2 * next_bucket];

      for (std::uint32_t position = begin; position < end; ++position) {
        const std::uint32_t sa_start = sa_start_at(bucket, position);
        const std::uint32_t sa_end = interval_end(bucket, position, end);

        if (sa_start == sa_end) {
          ++result.singleton_entries;
        } else {
          ++result.range_entries;
        }
      }
    }

    return result;
  }

 private:
  // ---------- Build alphatab once ----------

  void initialise_alphatab() {
    alphatab_[static_cast<unsigned char>('A')] = 0;  // 00
    alphatab_[static_cast<unsigned char>('C')] = 1;  // 01
    alphatab_[static_cast<unsigned char>('G')] = 2;  // 10
    alphatab_[static_cast<unsigned char>('T')] = 3;  // 11
  }

  // Packs the reference's last 1 .. kmer_length-1 characters once (see
  // short_suffix_keys_), so a miss never reads the reference or the input
  // to check them.
  void build_short_suffixes() {
    short_suffix_count_ =
        std::min(static_cast<std::size_t>(kmer_length - 1), ref_->size());

    for (std::size_t length = 1; length <= short_suffix_count_; ++length) {
      const std::size_t start = ref_->size() - length;
      std::uint32_t key = 0;

      for (std::size_t j = 0; j < length; ++j) {
        key = (key << 2U) |
              alphatab_[static_cast<unsigned char>((*ref_)[start + j])];
      }

      short_suffix_keys_[length] =
          key << (32U - 2U * static_cast<std::uint32_t>(length));
    }
  }

  // Empty-bucket miss: raises the match if a short suffix of the reference
  // shares a longer prefix with the query `key` than the table did. The
  // table's LCP is < 8 here, so a suffix of any length can beat it and every
  // one is checked (not only those in the query's bucket), but longest
  // first: a suffix of length L shares at most L characters, so the loop
  // stops as soon as L <= match_length.
  void check_short_suffixes_empty_bucket(const std::uint32_t key,
                                         std::size_t& ref_pos,
                                         std::size_t& match_length) const {
    for (std::size_t length = short_suffix_count_; length > match_length;
         --length) {
      // Capped at the suffix's own length: its unused low bits are 0 and
      // must not count as matches.
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

  // Non-empty-bucket miss: the neighbouring table entry already shares the
  // query's whole 8-character bucket prefix, so match_length >= 8. Only a
  // suffix of length >= 8 that lies in the same bucket can beat that --
  // any other shares < 8 characters -- so the rest are skipped.
  void check_short_suffixes_in_bucket(const std::uint32_t key,
                                      std::size_t& ref_pos,
                                      std::size_t& match_length) const {
    const std::uint32_t bucket = key >> low_bits;
    const std::size_t min_length = std::max<std::size_t>(match_length, 7);

    for (std::size_t length = short_suffix_count_; length > min_length;
         --length) {
      if ((short_suffix_keys_[length] >> low_bits) != bucket) {
        continue;
      }

      // Capped at the suffix's own length: its unused low bits are 0 and
      // must not count as matches.
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
    H_interleaved_.reserve(number_of_buckets * 2 + 1);
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

  LookupResult lookup(const std::uint32_t key) const {
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
      ++stats_.linear_bucket_searches;

      while (insertion_position < end && L_[insertion_position].low < low) {
        ++insertion_position;
      }

    }

    // Large bucket: binary search.
    else {
      ++stats_.binary_bucket_searches;

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

    check_short_suffixes_in_bucket(key, ref_pos, lcp_chars);

    return {false, 0, 0, ref_pos, lcp_chars};
  }
};
