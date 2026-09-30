#pragma once

// A PT16 table for powered_rlz's backward search, in the sassy layout
// (variants/pt16_sassy.hpp, variants/pt16_build_sassy.hpp), mirrored for a
// parse that runs right to left.
//
// powered's phrase starting at its right end i is the longest SUFFIX of
// seq[0, i) that occurs in the (cyclic) reference, found by backward search.
// So the table is keyed by the 16-mer read from its right end (the reversed
// key: its last character in the top two bits): then, exactly as a forward
// PT16 table gives the longest common PREFIX on a miss, this one gives the
// longest common SUFFIX -- which is powered's phrase when it is shorter
// than 16.
//
// A lookup of the 16-mer P = seq[i-16, i):
//   singleton hit: P occurs once; returns its reference position (from the
//         entry itself). The match is unique, so the parse extends it by
//         comparing characters with the reference -- no rank steps.
//   range hit: P occurs several times; returns its row interval [a, b) in the
//         powered index -- the interval powered's backward search holds after
//         matching P -- and one reference position.
//   miss: P does not occur, so the phrase is shorter than 16; returns its
//         length (the longest common suffix with a neighbouring key) and a
//         reference position of it (from that entry, or precomputed for an
//         empty bucket). No backward search at all.
// A singleton hit and a miss on a singleton neighbour read one L entry and
// nothing else: no row array, no gca_.
//
// Built from the reference and the powered index itself (REF_four.bwt), not
// from a suffix array: the index's rows are the rotations of the cyclic
// reference in sorted order and gca_[row] is the start of the rotation at
// that row. Every row has a full cyclic 16-mer, so every occurrence of any
// string ends a table 16-mer -- there are no short suffixes (sassy's short
// suffix records have no counterpart here).
//
// Layout, as the sassy table (entries in reversed-key order):
//   H           NUMBER_OF_BUCKETS + 1: index in L of each bucket's first
//               entry; an empty bucket is EMPTY_BUCKET_FLAG | the non-empty
//               bucket with the longest common prefix (build_H)
//   H_ranges    NUMBER_OF_BUCKETS: where each bucket's slices of
//               range_data begin (sassy's H_sa)
//   L           per entry 64 bits, sassy's L entry (pt16_sassy_format.hpp):
//               top 16 bits the reversed key's low 16 bits; bit 0 singleton
//               or range; a singleton's reference position, or a range's
//               count and offset into range_data (relative to H_ranges)
//   range_data  per range: [its count, if the entry's count field escaped],
//               its first row, the reference position of that row (sassy's
//               sampled_sa holds all positions; powered needs the rows)
// Derived at load, as in sassy: per empty bucket its miss answer.
//
// File: "PT16PS01", uint64 n (rows), uint64 entry_count, uint64
// range_data_count, H, H_ranges, L, range_data.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../pt16_utils.hpp"                // alphatab, buckets, build_H
#include "../variants/pt16_sassy_format.hpp"  // sassy_encode_*, sassy_decode_*

struct PoweredStart {
  bool found = false;
  bool singleton = false;
  // range hit: the row interval [a, b)
  std::uint64_t a = 0;
  std::uint64_t b = 0;
  // miss: the phrase length (< 16); hit: 16
  std::uint32_t match_length = 0;
  // hit: a reference position of the 16-mer; miss: of the phrase
  std::uint64_t position = 0;
};

// The shape of a built table.
struct PoweredTableStats {
  std::uint64_t rows = 0;
  std::uint64_t entries = 0;
  std::uint64_t singletons = 0;
  std::uint64_t ranges = 0;
  std::uint64_t largest_range = 0;
  std::uint64_t wrapping_rows = 0;  // rows whose 16-mer crosses the end
  std::uint64_t empty_buckets = 0;
};

// The 2-bit characters of a 16-mer key in reverse order: the key of the
// 16-mer read from its right end.
inline std::uint32_t reverse_16mer_key(std::uint32_t x) {
  x = ((x >> 2U) & 0x33333333U) | ((x & 0x33333333U) << 2U);
  x = ((x >> 4U) & 0x0F0F0F0FU) | ((x & 0x0F0F0F0FU) << 4U);
  x = ((x >> 8U) & 0x00FF00FFU) | ((x & 0x00FF00FFU) << 8U);
  return (x >> 16U) | (x << 16U);
}

// The forward key of the 16-mer of the cyclic reference at `position`.
inline std::uint32_t encode_cyclic_16mer(
    const std::vector<unsigned char>& reference, const std::uint64_t position) {
  const std::uint64_t n = reference.size();
  if (position + KMER_LENGTH <= n) {
    return encode_16mer(reference, static_cast<std::uint32_t>(position));
  }
  std::uint32_t key = 0;
  for (std::uint32_t j = 0; j < KMER_LENGTH; ++j) {
    key = (key << 2U) | alphatab[reference[(position + j) % n]];
  }
  return key;
}

class PT16PoweredTable {
 public:
  static constexpr std::uint32_t empty_bucket_mask = EMPTY_BUCKET_FLAG - 1;

  // One entry, decoded (for checks).
  struct EntryInfo {
    std::uint32_t reversed_key = 0;
    bool singleton = false;
    std::uint64_t a = 0, b = 0;  // range only
    std::uint64_t position = 0;
  };

  PT16PoweredTable() = default;

  // Builds the table from the reference and the index's rotation order
  // (`gca[row]` = start of that row's rotation). Throws if the reference is
  // not plain ACGT, does not match the index size, or the rows are not in
  // sorted 16-mer order.
  template <typename GCA>
  static PT16PoweredTable build(const std::vector<unsigned char>& reference,
                                const GCA& gca, PoweredTableStats* stats);

  static PT16PoweredTable load(const std::string& path);
  void write(const std::string& path) const;

  // Looks up the 16-mer text[0, 16).
  PoweredStart lookup(const unsigned char* text) const {
    return lookup_key(encode_text_reversed(text));
  }

  // The reversed key of text[0, 16): text[15] in the top two bits.
  static std::uint32_t encode_text_reversed(const unsigned char* text) {
    std::uint32_t key = 0;
    for (std::uint32_t j = KMER_LENGTH; j-- > 0;) {
      key = (key << 2U) | alphatab[text[j]];
    }
    return key;
  }

  PoweredStart lookup_key(const std::uint32_t key) const;

  std::uint64_t rows() const { return rows_; }
  std::uint64_t entries() const { return L_.size(); }
  // Calls f(EntryInfo) for every entry, in table order (for checks).
  template <typename F>
  void for_each_entry(F&& f) const {
    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
      if (H_[bucket] & EMPTY_BUCKET_FLAG) continue;
      for (std::uint32_t index = H_[bucket]; index < bucket_end(bucket); ++index) {
        const std::uint64_t e = L_[index];
        EntryInfo info;
        info.reversed_key = key_of(bucket, e);
        info.singleton = !sassy_is_range(e);
        if (info.singleton) {
          info.position = sassy_decode_position(e);
        } else {
          const RangeSlice slice = range_slice(bucket, e);
          info.a = slice.first_row;
          info.b = slice.first_row + slice.count;
          info.position = slice.position;
        }
        f(info);
      }
    }
  }
  std::uint64_t bytes() const {
    return (H_.size() + H_ranges_.size() + range_data_.size() +
            empty_position_.size()) *
               sizeof(std::uint32_t) +
           L_.size() * sizeof(std::uint64_t) + empty_length_.size();
  }

 private:
  // A range entry's slice of range_data: [first row, its position]; its row
  // count comes from the entry, or from range_data if it escaped.
  struct RangeSlice {
    std::uint64_t first_row;
    std::uint64_t count;
    std::uint64_t position;
  };

  RangeSlice range_slice(const std::uint32_t bucket,
                         const std::uint64_t entry) const {
    std::size_t start =
        static_cast<std::size_t>(H_ranges_[bucket]) + sassy_decode_offset(entry);
    std::uint64_t count = sassy_decode_count_field(entry);
    if (count == SASSY_COUNT_ESCAPE) {
      count = range_data_[start];
      ++start;
    }
    return {range_data_[start], count, range_data_[start + 1]};
  }

  // The full reversed key of an entry of `bucket`.
  static std::uint32_t key_of(const std::uint32_t bucket,
                              const std::uint64_t entry) {
    return (bucket << LOW_BITS) | sassy_decode_low(entry);
  }

  // One reference position of the 16-mer of the entry at `index` in L.
  std::uint64_t first_position(const std::uint32_t bucket,
                               const std::uint32_t index) const {
    const std::uint64_t entry = L_[index];
    return sassy_is_range(entry) ? range_slice(bucket, entry).position
                                 : sassy_decode_position(entry);
  }

  // A miss answered by entry `index` sharing `shared` characters (the common
  // suffix): the phrase is the last `shared` characters of its 16-mer.
  PoweredStart miss(const std::uint32_t bucket, const std::uint32_t index,
                    const std::uint32_t shared) const {
    PoweredStart result;
    result.match_length = shared;
    result.position =
        (first_position(bucket, index) + KMER_LENGTH - shared) % rows_;
    return result;
  }

  // Index in L one past the last entry of a non-empty bucket.
  std::uint32_t bucket_end(const std::uint32_t bucket) const {
    std::uint32_t next = bucket + 1;
    while (next < NUMBER_OF_BUCKETS && (H_[next] & EMPTY_BUCKET_FLAG)) ++next;
    return H_[next];
  }

  // Per empty bucket: the answer for any query in it (sassy's
  // build_empty_answers): the longest common prefix with the nearest
  // non-empty bucket, and a position from that bucket's first entry.
  void build_empty_answers();

  std::uint64_t rows_ = 0;
  std::vector<std::uint32_t> H_;
  std::vector<std::uint32_t> H_ranges_;
  std::vector<std::uint64_t> L_;
  std::vector<std::uint32_t> range_data_;

  std::vector<std::uint32_t> empty_position_;
  std::vector<std::uint8_t> empty_length_;
};

inline PoweredStart PT16PoweredTable::lookup_key(const std::uint32_t key) const {
  const std::uint32_t bucket = key >> LOW_BITS;

  // Empty bucket: no 16-mer ends with these 8 characters (precomputed).
  if (H_[bucket] & EMPTY_BUCKET_FLAG) {
    PoweredStart result;
    result.match_length = empty_length_[bucket];
    result.position = empty_position_[bucket];
    return result;
  }

  const std::uint16_t low = static_cast<std::uint16_t>(key & LOW_MASK);
  const std::uint32_t begin = H_[bucket];
  const std::uint32_t end = bucket_end(bucket);

  // Small bucket: linear scan. Large bucket: binary search.
  std::uint32_t at = begin;
  if (end - begin < BINARY_SEARCH_THRESHOLD) {
    while (at < end && sassy_decode_low(L_[at]) < low) ++at;
  } else {
    at = static_cast<std::uint32_t>(
        std::lower_bound(L_.begin() + begin, L_.begin() + end, low,
                         [](const std::uint64_t entry, const std::uint16_t value) {
                           return sassy_decode_low(entry) < value;
                         }) -
        L_.begin());
  }

  // Exact hit.
  if (at < end && sassy_decode_low(L_[at]) == low) {
    const std::uint64_t entry = L_[at];
    PoweredStart result;
    result.found = true;
    result.match_length = KMER_LENGTH;
    if (sassy_is_range(entry)) {
      const RangeSlice slice = range_slice(bucket, entry);
      result.a = slice.first_row;
      result.b = slice.first_row + slice.count;
      result.position = slice.position;
    } else {
      result.singleton = true;
      result.position = sassy_decode_position(entry);
    }
    return result;
  }

  // Miss: the neighbouring key with the longest common prefix (reversed: the
  // longest common suffix, >= 8 characters).
  std::uint32_t best;
  if (at == begin) {
    best = begin;
  } else if (at == end) {
    best = end - 1;
  } else {
    const int predecessor_lcp = std::countl_zero(key ^ key_of(bucket, L_[at - 1]));
    const int successor_lcp = std::countl_zero(key ^ key_of(bucket, L_[at]));
    best = predecessor_lcp >= successor_lcp ? at - 1 : at;
  }
  const auto shared = static_cast<std::uint32_t>(
      std::countl_zero(key ^ key_of(bucket, L_[best])) / 2);
  return miss(bucket, best, shared);
}

inline void PT16PoweredTable::build_empty_answers() {
  empty_position_.assign(NUMBER_OF_BUCKETS, 0);
  empty_length_.assign(NUMBER_OF_BUCKETS, 0);
  for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
    if ((H_[bucket] & EMPTY_BUCKET_FLAG) == 0) continue;
    const std::uint32_t matching_bucket = H_[bucket] & empty_bucket_mask;
    if (matching_bucket >= NUMBER_OF_BUCKETS) continue;  // no entries at all
    const std::uint32_t difference = (bucket ^ matching_bucket) << 16U;
    const auto shared =
        static_cast<std::uint32_t>(std::countl_zero(difference)) / 2;
    const PoweredStart answer = miss(matching_bucket, H_[matching_bucket], shared);
    empty_length_[bucket] = static_cast<std::uint8_t>(answer.match_length);
    empty_position_[bucket] = static_cast<std::uint32_t>(answer.position);
  }
}

template <typename GCA>
PT16PoweredTable PT16PoweredTable::build(
    const std::vector<unsigned char>& reference, const GCA& gca,
    PoweredTableStats* stats) {
  const std::uint64_t n = reference.size();

  if (gca.size() != n) {
    throw std::runtime_error(
        "index has " + std::to_string(gca.size()) + " rows, reference " +
        std::to_string(n) + " characters: not the index of this reference");
  }
  if (n < KMER_LENGTH || n >= (std::uint64_t{1} << 32)) {
    throw std::runtime_error("reference length must be in [16, 2^32)");
  }
  for (std::uint64_t i = 0; i < n; ++i) {
    if (!is_acgt(reference[i])) {
      throw std::runtime_error("reference is not plain ACGT (position " +
                               std::to_string(i) + ")");
    }
  }

  PoweredTableStats local;
  local.rows = n;

  // ---------- Pass 1: the row intervals, in row (forward key) order ----------
  //
  // gca is read in row order here, and each entry keeps the position of its
  // first row: pass 2 visits the entries in reversed-key order, which is a
  // random order of rows, and reading gca there would miss the cache once
  // per entry. The reference windows are prefetched as in the other
  // builders (prefetch_window).
  struct Entry {
    std::uint32_t reversed_key;
    std::uint32_t start;
    std::uint32_t end;
    std::uint32_t position;  // gca[start]
  };
  std::vector<Entry> entries;
  entries.reserve(n / 2);
  std::array<std::uint32_t, NUMBER_OF_BUCKETS> count{};

  std::uint32_t previous = 0;
  std::uint64_t range_start = 0;
  std::uint64_t range_position = 0;

  auto close_entry = [&](const std::uint64_t end) {
    const std::uint32_t reversed_key = reverse_16mer_key(previous);
    entries.push_back({reversed_key, static_cast<std::uint32_t>(range_start),
                       static_cast<std::uint32_t>(end),
                       static_cast<std::uint32_t>(range_position)});
    ++count[reversed_key >> LOW_BITS];
    const std::uint64_t size = end - range_start;
    if (size == 1) {
      ++local.singletons;
    } else {
      ++local.ranges;
      local.largest_range = std::max(local.largest_range, size);
    }
  };

  for (std::uint64_t row = 0; row < n; ++row) {
    if (BUILD_PREFETCH_DISTANCE != 0 && row + BUILD_PREFETCH_DISTANCE < n) {
      prefetch_window(reference, gca[row + BUILD_PREFETCH_DISTANCE]);
    }

    const std::uint64_t position = gca[row];
    if (position >= n) {
      throw std::runtime_error("gca value out of range at row " +
                               std::to_string(row));
    }
    if (position + KMER_LENGTH > n) ++local.wrapping_rows;

    const std::uint32_t key = encode_cyclic_16mer(reference, position);

    if (row != 0 && key == previous) continue;
    if (row != 0 && key < previous) {
      throw std::runtime_error(
          "rows are not in sorted 16-mer order at row " + std::to_string(row) +
          ": the index is not a sorted rotation order of this reference");
    }

    if (row != 0) close_entry(row);
    previous = key;
    range_start = row;
    range_position = position;
  }
  close_entry(n);

  // ---------- Pass 2: reversed-key order, sassy entries ----------
  //
  // Sorted by a counting sort on the bucket (the reversed key's high 16 bits,
  // counted in pass 1), then each bucket by its low 16 bits: the same order
  // as sorting by the whole reversed key (the keys are distinct).
  {
    std::vector<std::uint32_t> next(NUMBER_OF_BUCKETS);
    std::uint32_t total = 0;
    for (std::uint32_t b = 0; b < NUMBER_OF_BUCKETS; ++b) {
      next[b] = total;
      total += count[b];
    }
    std::vector<Entry> sorted(entries.size());
    for (const Entry& e : entries) sorted[next[e.reversed_key >> LOW_BITS]++] = e;
    std::vector<Entry>().swap(entries);
    std::uint32_t begin = 0;
    for (std::uint32_t b = 0; b < NUMBER_OF_BUCKETS; ++b) {
      std::sort(sorted.begin() + begin, sorted.begin() + next[b],
                [](const Entry& x, const Entry& y) {
                  return x.reversed_key < y.reversed_key;
                });
      begin = next[b];
    }
    entries.swap(sorted);
  }

  PT16PoweredTable table;
  table.rows_ = n;
  table.H_ranges_.assign(NUMBER_OF_BUCKETS, 0);
  table.L_.reserve(entries.size());

  count.fill(0);
  for (const Entry& e : entries) {
    const std::uint32_t bucket = e.reversed_key >> LOW_BITS;
    const auto low = static_cast<std::uint16_t>(e.reversed_key & LOW_MASK);

    // The first entry in the bucket gives the bucket's range_data start.
    if (count[bucket] == 0) {
      table.H_ranges_[bucket] =
          static_cast<std::uint32_t>(table.range_data_.size());
    }

    const std::uint32_t occurrences = e.end - e.start;
    if (occurrences == 1) {
      table.L_.push_back(sassy_encode_singleton(low, e.position));
    } else {
      const auto offset = static_cast<std::uint32_t>(
          table.range_data_.size() - table.H_ranges_[bucket]);
      std::uint16_t count_field;
      if (occurrences < SASSY_COUNT_ESCAPE) {
        count_field = static_cast<std::uint16_t>(occurrences);
      } else {
        count_field = SASSY_COUNT_ESCAPE;
        table.range_data_.push_back(occurrences);
      }
      table.range_data_.push_back(e.start);
      table.range_data_.push_back(e.position);
      table.L_.push_back(sassy_encode_range(low, offset, count_field));
    }
    ++count[bucket];
  }
  table.H_ = build_H(count);
  table.build_empty_answers();

  for (const std::uint32_t c : count) local.empty_buckets += c == 0 ? 1 : 0;
  local.entries = table.L_.size();
  if (stats != nullptr) *stats = local;
  return table;
}

inline void PT16PoweredTable::write(const std::string& path) const {
  std::ofstream output(path, std::ios::binary);
  if (!output) throw std::runtime_error("cannot create table: " + path);

  const char magic[8] = {'P', 'T', '1', '6', 'P', 'S', '0', '1'};
  const std::uint64_t entry_count = L_.size();
  const std::uint64_t range_data_count = range_data_.size();
  output.write(magic, sizeof(magic));
  write_value(output, rows_);
  write_value(output, entry_count);
  write_value(output, range_data_count);
  write_vector(output, H_);
  write_vector(output, H_ranges_);
  write_vector(output, L_);
  write_vector(output, range_data_);

  if (!output) throw std::runtime_error("failed while writing table: " + path);
}

inline PT16PoweredTable PT16PoweredTable::load(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open table: " + path);

  char magic[8];
  input.read(magic, sizeof(magic));
  if (!input || std::string(magic, 8) != "PT16PS01") {
    throw std::runtime_error(
        "not a (sassy) powered PT16 table: " + path +
        " -- rebuild it with pt16_powered_build");
  }

  PT16PoweredTable table;
  std::uint64_t entry_count = 0, range_data_count = 0;
  input.read(reinterpret_cast<char*>(&table.rows_), sizeof(table.rows_));
  input.read(reinterpret_cast<char*>(&entry_count), sizeof(entry_count));
  input.read(reinterpret_cast<char*>(&range_data_count),
             sizeof(range_data_count));

  table.H_.resize(static_cast<std::size_t>(NUMBER_OF_BUCKETS) + 1);
  table.H_ranges_.resize(NUMBER_OF_BUCKETS);
  table.L_.resize(entry_count);
  table.range_data_.resize(range_data_count);
  input.read(reinterpret_cast<char*>(table.H_.data()),
             table.H_.size() * sizeof(std::uint32_t));
  input.read(reinterpret_cast<char*>(table.H_ranges_.data()),
             table.H_ranges_.size() * sizeof(std::uint32_t));
  input.read(reinterpret_cast<char*>(table.L_.data()),
             table.L_.size() * sizeof(std::uint64_t));
  input.read(reinterpret_cast<char*>(table.range_data_.data()),
             table.range_data_.size() * sizeof(std::uint32_t));

  if (!input) throw std::runtime_error("truncated table: " + path);
  if (table.H_.back() != table.L_.size()) {
    throw std::runtime_error("PT16 H directory does not end at the entry count");
  }
  table.build_empty_answers();
  return table;
}
