#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "pt16_sassy_format.hpp"  // PackedShortSuffix, sassy_decode_*
#include "pt16_utils.hpp"         // KMER_LENGTH, buckets
#include "rlz_common.hpp"         // binarySearchLB, binarySearchRB

/**
 * Lookup over the PT16 table written by pt16_build_sassy.hpp
 * (write_hl_table). The table is self-contained: neither the reference nor
 * the suffix array is needed to read it or to use it.
 *
 * !!! CAVEAT: A LOOKUP NEEDS 16 CHARACTERS, A TAIL NEEDS lookup_tail !!!
 *
 * lookup takes the 16 characters at the start of a query. A query with fewer
 * than 16 characters left (the tail of an input) has no 16-mer key. The RLZ
 * parsers handled such a tail with plain RLZ, which binary-searches the suffix
 * array; this table has no suffix array to fall back on, so a tail must go
 * through lookup_tail instead, which pads it with zero characters ('A') and
 * does an ordinary 16-mer lookup. lookup(text, position) refuses a query with
 * fewer than 16 characters left instead of guessing, and lookup(key) has no
 * way to know, so its caller must choose the right function.
 *
 * The table knows only the first 16 characters of each suffix, so no lookup
 * reports a match longer than 16. Extending a match beyond that needs the
 * reference (or the suffix array), which is not part of this index.
 *
 * File layout (all integers little-endian, as written by the builder):
 *
 *   magic          8 bytes    "PT16SA03"
 *   entry_count    uint64     number of L entries
 *   sampled_count  uint64     number of sampled_sa entries
 *   short_count    uint64     number of short suffix records
 *   H              (NUMBER_OF_BUCKETS + 1) x uint32
 *   H_sa           NUMBER_OF_BUCKETS x uint32
 *   L              entry_count x uint64
 *   sampled_sa     sampled_count x uint32
 *   short suffixes short_count x PackedShortSuffix (3 x uint32)
 *
 * H[b] is the index in L of the first entry of bucket b. An empty bucket has
 * EMPTY_BUCKET_FLAG set and the non-empty bucket with the longest common
 * prefix in the low 31 bits. H[NUMBER_OF_BUCKETS] equals entry_count.
 *
 * H_sa[b] is the index in sampled_sa where the slices of bucket b begin.
 *
 * Each L entry is 64 bits, and the entries of a bucket are sorted by `low`.
 * See the "Sassy L entry" section of pt16_utils.hpp for the exact bit layout
 * (sassy_decode_low, sassy_is_range, sassy_decode_position,
 * sassy_decode_offset, sassy_decode_count_field): in short, bit 0 says
 * whether the entry is a singleton or a range, a singleton's text position
 * uses the other 47 bits, and a range uses 31 bits for its offset into
 * sampled_sa and 16 bits for its count, escaping to a count stored in
 * sampled_sa itself when 16 bits are not enough.
 *
 * Suffixes of the reference shorter than KMER_LENGTH have no table entry,
 * because an entry needs a full 16-mer. The builder stores every one of them
 * as a short suffix record (see PackedShortSuffix in pt16_utils.hpp), and a
 * lookup that misses the table matches them too.
 */
class PT16SassyLookup {
 public:
  struct Stats {
    std::size_t entries = 0;
    std::size_t sampled_entries = 0;
    std::size_t short_suffixes = 0;
    std::size_t approx_bytes = 0;

    // Per-query counts, updated by every call to lookup(const std::uint32_t)
    // (matching the older PT16RLZParser's Stats shape, so this class plugs
    // into the same PT16Delta/CSVWriter reporting unchanged). Unlike that
    // class, lookup_tail's own padded lookup counts too, since it genuinely
    // queries the table via a padded key rather than bypassing it the way
    // the older parsers' separate <16-character fallback did.
    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t singleton_hits = 0;
    std::size_t range_hits = 0;
  };

  /**
   * The answer to looking up one 16-mer.
   *
   * found = true:  the 16-mer occurs in the reference.
   *   count == 1:  `position` is its text position.
   *   count  > 1:  `positions` holds all `count` text positions, in suffix
   *                array order, as a view into the table.
   *
   * found = false: it does not occur.
   *
   * In both cases match_length is the length of the longest prefix of the
   * 16-mer that occurs in the reference (16 when found) and match_position is
   * one place where that prefix occurs (the first occurrence when found).
   * The length is exact; the position is one of possibly many.
   */
  struct LookupResult {
    bool found = false;

    std::uint32_t count = 0;
    std::uint32_t position = 0;
    std::span<const std::uint32_t> positions;

    std::uint32_t match_position = 0;
    std::uint32_t match_length = 0;
  };

  /**
   * The answer to looking up a tail: a query of 1 to 15 characters.
   *
   * match_length is the length of the longest prefix of the tail that occurs
   * in the reference, and match_position is one place where that prefix
   * occurs. found = true means the whole tail occurs (match_length equals the
   * tail's length).
   *
   * There is no count or list of occurrences: which 16-mer was looked up is an
   * implementation detail (see lookup_tail), so its occurrences are only some
   * of the tail's.
   */
  struct TailResult {
    bool found = false;
    std::uint32_t match_position = 0;
    std::uint32_t match_length = 0;
  };

  explicit PT16SassyLookup(const std::string& path) { load(path); }

  const Stats& stats() const { return stats_; }

  /**
   * Looks up a 16-mer packed as by encode_16mer: 2 bits per character, the
   * first character in the top bits. See the caveat above: this is only for
   * a query that has all 16 characters.
   */
  LookupResult lookup(const std::uint32_t key) const {
    const std::uint32_t bucket = key >> LOW_BITS;
    const std::uint16_t low = static_cast<std::uint16_t>(key & LOW_MASK);

    LookupResult result;

    // ---------- Empty bucket: no 16-mer starts with these 8 characters ----------

    if (H_[bucket] & EMPTY_BUCKET_FLAG) {
      // The non-empty bucket sharing the most leading bits; every entry of
      // it shares that many characters with the query.
      const std::uint32_t matching = H_[bucket] & empty_bucket_mask;

      result.match_length = lcp_bits_16(bucket, matching) / 2;
      result.match_position = first_position(matching, H_[matching]);

      match_short_suffixes(key, result);
      ++stats_.misses;
      return result;
    }

    // ---------- Non-empty bucket: search it for the low 16 bits ----------

    const std::uint32_t begin = H_[bucket];
    const std::uint32_t end = bucket_end(bucket);
    const std::uint32_t at = lower_bound_low(begin, end, low);

    // Exact hit.
    if (at < end && sassy_decode_low(L_[at]) == low) {
      const std::uint64_t entry = L_[at];

      result.found = true;
      result.match_length = KMER_LENGTH;

      if (sassy_is_range(entry)) {
        result.positions = entry_slice(bucket, entry);
        result.count = static_cast<std::uint32_t>(result.positions.size());
        result.match_position = result.positions.front();
        ++stats_.range_hits;
      } else {
        result.count = 1;
        result.position = static_cast<std::uint32_t>(sassy_decode_position(entry));
        result.match_position = result.position;
        ++stats_.singleton_hits;
      }

      ++stats_.hits;
      return result;
    }

    // ---------- Miss: the neighbouring key with the longest common prefix ----------

    std::uint32_t best;

    if (at == begin) {
      best = begin;  // smaller than everything in the bucket
    } else if (at == end) {
      best = end - 1;  // larger than everything in the bucket
    } else {
      const std::uint32_t predecessor = at - 1;
      const std::uint32_t successor = at;

      const int predecessor_lcp =
          std::countl_zero(key ^ key_of(bucket, L_[predecessor]));
      const int successor_lcp =
          std::countl_zero(key ^ key_of(bucket, L_[successor]));

      best = predecessor_lcp >= successor_lcp ? predecessor : successor;
    }

    result.match_length =
        static_cast<std::uint32_t>(std::countl_zero(key ^ key_of(bucket, L_[best])) / 2);
    result.match_position = first_position(bucket, best);

    match_short_suffixes(key, result);
    ++stats_.misses;
    return result;
  }

  /**
   * Looks up the 16 characters starting at text[position].
   *
   * Throws std::out_of_range if fewer than 16 characters remain: use
   * lookup_tail for such a query (see the caveat above). It is better to
   * refuse it than to answer it wrongly.
   */
  LookupResult lookup(const std::vector<unsigned char>& text,
                      const std::size_t position) const {
    if (position > text.size() || text.size() - position < KMER_LENGTH) {
      throw std::out_of_range(
          "PT16 sassy lookup needs 16 characters; fewer remain at position " +
          std::to_string(position) + " (use lookup_tail)");
    }

    return lookup(encode_16mer(text, static_cast<std::uint32_t>(position)));
  }

  /**
   * Looks up the tail text[position, text.size()), which must be 1 to 15
   * characters long; throws std::out_of_range otherwise.
   *
   * The tail is padded with zero characters ('A') to a 16-mer P, and P is
   * looked up as usual. The longest prefix of P that occurs is at least as
   * long as the tail exactly when the whole tail occurs, and otherwise it is
   * the longest prefix of the tail that occurs, because the tail is the first
   * characters of P. So the result is min(that length, tail length), which is
   * exact. If P itself is in the table, the tail occurs where P does. If not,
   * every 16-mer that starts with the tail sorts at or just after P (the
   * padding is the smallest possible), so the closest 16-mer the lookup finds
   * covers the tail whenever any 16-mer does. A tail that occurs only at the
   * very end of the reference is found among the short suffixes.
   */
  TailResult lookup_tail(const std::vector<unsigned char>& text,
                         const std::size_t position) const {
    if (position >= text.size() || text.size() - position >= KMER_LENGTH) {
      throw std::out_of_range(
          "PT16 sassy lookup_tail needs 1 to 15 characters; found " +
          std::to_string(position > text.size() ? 0
                                                : text.size() - position) +
          " at position " + std::to_string(position));
    }

    const std::uint32_t length =
        static_cast<std::uint32_t>(text.size() - position);

    // Pack the tail from the top bit, leaving the padding as zero bits.
    std::uint32_t key = 0;

    for (std::uint32_t j = 0; j < length; ++j) {
      const std::uint32_t code =
          alphatab[static_cast<unsigned char>(text[position + j])];
      key |= code << (30U - 2U * j);
    }

    const LookupResult padded = lookup(key);

    // A match of the padded key may run into the padding; the tail ends
    // after `length` characters.
    TailResult result;
    result.match_length = std::min(padded.match_length, length);
    result.match_position = padded.match_position;
    result.found = result.match_length == length;

    return result;
  }

  /**
   * Finds the longest match of text[position, text.size()) against
   * `reference`: an LZ factor, as (reference position, length).
   *
   * Requires text.size() - position >= KMER_LENGTH (throws std::out_of_range
   * otherwise, from the lookup(text, position) below; use lookup_tail
   * directly for a shorter tail, see the caveat above).
   *
   * Looks up the 16-mer at text[position, position + 16), then extends the
   * match beyond those 16 characters by reading `reference` directly:
   *
   *   - a singleton hit is extended character by character, comparing text
   *     against reference from character 17 onward;
   *   - a range hit is narrowed with the same left/right binary search the
   *     plain RLZ parser runs over the whole suffix array (rlz::
   *     binarySearchLB/RB, rlz_common.hpp), but scoped to just this 16-mer's
   *     own occurrences (hit.positions) instead: since that slice is already
   *     exactly the suffix-array range for this key, searching it directly
   *     needs no suffix array at all. Once it narrows to one occurrence, the
   *     rest is extended character by character, as for a singleton;
   *   - a miss returns lookup's own closest match unchanged: match_length is
   *     already the longest prefix that occurs anywhere in the reference (of
   *     the 16-mer itself, or of a short suffix), so there is nothing to
   *     extend.
   *
   * The caller is responsible for passing the same reference the table was
   * built from; this class does not keep one, so it cannot check that.
   */
  std::pair<std::uint32_t, std::uint32_t> find_longest_matching_factor(
      const std::vector<unsigned char>& text, const std::size_t position,
      const std::vector<unsigned char>& reference) const {
    const LookupResult hit = lookup(text, position);

    if (!hit.found) {
      return {hit.match_position, hit.match_length};
    }

    std::size_t offset = KMER_LENGTH;
    std::size_t j = position + KMER_LENGTH;

    if (hit.count == 1) {
      const std::uint32_t match = hit.position;

      while (j < text.size() && match + offset < reference.size() &&
             reference[match + offset] == text[j]) {
        ++j;
        ++offset;
      }

      return {match, static_cast<std::uint32_t>(offset)};
    }

    // Range: narrow among the 16-mer's own occurrences, exactly as the plain
    // RLZ parser narrows [nlb, nrb] within the suffix array, but with
    // hit.positions standing in for it.
    std::int64_t nlb = 0;
    std::int64_t nrb = static_cast<std::int64_t>(hit.positions.size()) - 1;

    while (nlb < nrb && j < text.size()) {
      const auto lb = rlz::binarySearchLB(reference, hit.positions, nlb, nrb,
                                          static_cast<std::int64_t>(offset),
                                          text[j]);

      if (!lb) {
        break;
      }

      const auto rb =
          rlz::binarySearchRB(reference, hit.positions, lb.value(), nrb,
                              static_cast<std::int64_t>(offset), text[j]);

      if (!rb) {
        break;
      }

      nlb = lb.value();
      nrb = rb.value();
      ++j;
      ++offset;
    }

    const std::uint32_t match = hit.positions[static_cast<std::size_t>(nlb)];

    // Narrowed to one occurrence: extend directly, as for a singleton. If
    // several occurrences are still tied when the text runs out, `match` is
    // arbitrarily the first of them; any of them is an equally valid answer.
    if (nlb == nrb) {
      while (j < text.size() && match + offset < reference.size() &&
             reference[match + offset] == text[j]) {
        ++j;
        ++offset;
      }
    }

    return {match, static_cast<std::uint32_t>(offset)};
  }

  /**
   * Runs the full greedy LZ factorization of `input` against `reference`:
   * the same greedy longest-match parse as rlz::lzFactorize (parser.hpp,
   * the baseline) and every PT16RLZParser::lzFactorize, using
   * find_longest_matching_factor for each position with 16 or more
   * characters left, and lookup_tail for the last few (see the caveat at
   * the top of this file). A match of 0 or 1 characters becomes a literal
   * factor, exactly as the other parsers do.
   *
   * The caller is responsible for passing the same reference the table was
   * built from, as for find_longest_matching_factor.
   */
  Triples lzFactorize(const std::vector<unsigned char>& input,
                      const std::vector<unsigned char>& reference) const {
    Triples factors;
    std::size_t i = 0;

    while (i < input.size()) {
      std::size_t pos;
      std::size_t len;

      if (input.size() - i < KMER_LENGTH) {
        const TailResult tail = lookup_tail(input, i);
        pos = tail.match_position;
        len = tail.match_length;
      } else {
        const auto [factor_pos, factor_len] =
            find_longest_matching_factor(input, i, reference);
        pos = factor_pos;
        len = factor_len;
      }

      if (len <= 1) {
        pos = static_cast<std::size_t>(input[i]);
        len = 1;
      }

      factors.push_back({i, pos, len});
      i += len;
    }

    return factors;
  }

  // The loaded arrays, exactly as the builder wrote them.
  const std::vector<std::uint32_t>& H() const { return H_; }
  const std::vector<std::uint32_t>& H_sa() const { return H_sa_; }
  const std::vector<std::uint64_t>& L() const { return L_; }
  const std::vector<std::uint32_t>& sampled_sa() const { return sampled_sa_; }
  const std::vector<PackedShortSuffix>& short_suffixes() const {
    return short_suffixes_;
  }

 private:
  // Must match the magic written by write_hl_table in pt16_build_sassy.hpp.
  static constexpr char magic_[8] = {'P', 'T', '1', '6', 'S', 'A', '0', '3'};

  static constexpr std::uint32_t empty_bucket_mask = EMPTY_BUCKET_FLAG - 1;

  // Buckets with fewer entries than this are searched linearly.
  static constexpr std::uint32_t binary_search_threshold = 64;

  // There is at most one short suffix per length in
  // [SHORT_SUFFIX_MIN_LENGTH, KMER_LENGTH).
  static constexpr std::uint64_t max_short_suffixes =
      KMER_LENGTH - SHORT_SUFFIX_MIN_LENGTH;

  // Starting index in L of each bucket, with the empty-bucket encoding.
  std::vector<std::uint32_t> H_;

  // Starting index in sampled_sa of each bucket's slices.
  std::vector<std::uint32_t> H_sa_;

  // One 64-bit entry per distinct 16-mer.
  std::vector<std::uint64_t> L_;

  // Text positions of every 16-mer that occurs more than once.
  std::vector<std::uint32_t> sampled_sa_;

  // Suffixes of the reference of length SHORT_SUFFIX_MIN_LENGTH..KMER_LENGTH-1.
  std::vector<PackedShortSuffix> short_suffixes_;

  // hits/misses/singleton_hits/range_hits are updated inside lookup(), which
  // is logically read-only (querying the table does not change it).
  mutable Stats stats_;

  // ---------- Decoding an L entry ----------
  //
  // The low-level bit decoding (sassy_decode_low, sassy_is_range,
  // sassy_decode_position, sassy_decode_offset, sassy_decode_count_field) is
  // shared with the builder, in pt16_utils.hpp. What is left here needs the
  // table's own arrays: turning a range's count field into its actual
  // position count and slice start, which means following the escape into
  // sampled_sa when the count field says to.

  // The text positions of a range entry, as a view into sampled_sa_. The
  // slice was checked at load time (see validate), so this never reads out
  // of bounds.
  std::span<const std::uint32_t> entry_slice(const std::uint32_t bucket,
                                             const std::uint64_t entry) const {
    std::size_t start =
        static_cast<std::size_t>(H_sa_[bucket]) + sassy_decode_offset(entry);
    std::uint32_t count = sassy_decode_count_field(entry);

    if (count == SASSY_COUNT_ESCAPE) {
      // The real count does not fit in the entry; it was pushed onto
      // sampled_sa right before the positions themselves.
      count = sampled_sa_[start];
      ++start;
    }

    return {sampled_sa_.data() + start, count};
  }

  // The full 16-mer key of an entry of `bucket`.
  static std::uint32_t key_of(const std::uint32_t bucket,
                              const std::uint64_t entry) {
    return (bucket << LOW_BITS) | sassy_decode_low(entry);
  }

  // One text position of the entry at index `index` of L, which is in
  // `bucket`.
  std::uint32_t first_position(const std::uint32_t bucket,
                               const std::uint32_t index) const {
    const std::uint64_t entry = L_[index];

    return sassy_is_range(entry)
              ? entry_slice(bucket, entry).front()
              : static_cast<std::uint32_t>(sassy_decode_position(entry));
  }

  // ---------- Searching a bucket ----------

  // Index in L one past the last entry of a non-empty bucket: the first entry
  // of the next non-empty bucket, or entry_count after the last one.
  std::uint32_t bucket_end(const std::uint32_t bucket) const {
    std::uint32_t next = bucket + 1;

    while (next < NUMBER_OF_BUCKETS && (H_[next] & EMPTY_BUCKET_FLAG)) {
      ++next;
    }

    return H_[next];
  }

  // Index of the first entry in [begin, end) whose low part is not below
  // `low`, or `end`.
  std::uint32_t lower_bound_low(const std::uint32_t begin,
                                const std::uint32_t end,
                                const std::uint16_t low) const {
    if (end - begin < binary_search_threshold) {
      std::uint32_t at = begin;

      while (at < end && sassy_decode_low(L_[at]) < low) {
        ++at;
      }

      return at;
    }

    const auto it = std::lower_bound(
        L_.begin() + begin, L_.begin() + end, low,
        [](const std::uint64_t entry, const std::uint16_t value) {
          return sassy_decode_low(entry) < value;
        });

    return static_cast<std::uint32_t>(it - L_.begin());
  }

  // ---------- Short suffixes ----------

  // Raises the match in `result` if a short suffix shares a longer prefix with
  // the query than the table does. A short suffix is at most 15 characters,
  // so this cannot turn a miss into a hit.
  void match_short_suffixes(const std::uint32_t key,
                            LookupResult& result) const {
    for (const PackedShortSuffix& suffix : short_suffixes_) {
      // Leading characters shared with the query, capped at the suffix's own
      // length (its unused low bits are 0 and must not count).
      const std::uint32_t shared = std::min<std::uint32_t>(
          suffix.length,
          static_cast<std::uint32_t>(std::countl_zero(key ^ suffix.packed) / 2));

      if (shared > result.match_length) {
        result.match_length = shared;
        result.match_position = suffix.ref_pos;
      }
    }
  }

  // ---------- Loading ----------

  template <typename T>
  static void read_value(std::ifstream& input, T& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(T));

    if (!input) {
      throw std::runtime_error("Failed while reading PT16 sassy table");
    }
  }

  template <typename T>
  static void read_vector(std::ifstream& input, std::vector<T>& values,
                          const std::uint64_t count) {
    values.resize(static_cast<std::size_t>(count));

    input.read(reinterpret_cast<char*>(values.data()),
               static_cast<std::streamsize>(values.size() * sizeof(T)));

    if (!input) {
      throw std::runtime_error("Failed while reading PT16 sassy table");
    }
  }

  void load(const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);

    if (!input) {
      throw std::runtime_error("Cannot open PT16 sassy table: " + path);
    }

    const std::uint64_t file_size = static_cast<std::uint64_t>(input.tellg());
    input.seekg(0);

    // ---------- Header ----------

    char magic[sizeof(magic_)]{};
    input.read(magic, sizeof(magic));

    if (!input || std::memcmp(magic, magic_, sizeof(magic)) != 0) {
      throw std::runtime_error("Not a PT16 sassy table (bad magic): " + path);
    }

    std::uint64_t entry_count = 0;
    std::uint64_t sampled_count = 0;
    std::uint64_t short_count = 0;

    read_value(input, entry_count);
    read_value(input, sampled_count);
    read_value(input, short_count);

    // Check the sizes before allocating anything, so that a corrupt count
    // cannot cause a huge allocation or an overflow below.
    if (entry_count > file_size / sizeof(std::uint64_t) ||
        sampled_count > file_size / sizeof(std::uint32_t) ||
        short_count > max_short_suffixes) {
      throw std::runtime_error(
          "PT16 sassy table counts are out of range: " + path);
    }

    const std::uint64_t expected_size =
        sizeof(magic_) + 3 * sizeof(std::uint64_t) +
        (static_cast<std::uint64_t>(NUMBER_OF_BUCKETS) + 1) *
            sizeof(std::uint32_t) +
        static_cast<std::uint64_t>(NUMBER_OF_BUCKETS) * sizeof(std::uint32_t) +
        entry_count * sizeof(std::uint64_t) +
        sampled_count * sizeof(std::uint32_t) +
        short_count * sizeof(PackedShortSuffix);

    if (file_size != expected_size) {
      throw std::runtime_error(
          "PT16 sassy table has the wrong size: expected " +
          std::to_string(expected_size) + " bytes, found " +
          std::to_string(file_size) + " (" + path + ")");
    }

    // ---------- Arrays, in the order the builder wrote them ----------

    read_vector(input, H_, static_cast<std::uint64_t>(NUMBER_OF_BUCKETS) + 1);
    read_vector(input, H_sa_, NUMBER_OF_BUCKETS);
    read_vector(input, L_, entry_count);
    read_vector(input, sampled_sa_, sampled_count);
    read_vector(input, short_suffixes_, short_count);

    validate(entry_count, sampled_count);

    stats_.entries = L_.size();
    stats_.sampled_entries = sampled_sa_.size();
    stats_.short_suffixes = short_suffixes_.size();
    stats_.approx_bytes = H_.size() * sizeof(std::uint32_t) +
                          H_sa_.size() * sizeof(std::uint32_t) +
                          L_.size() * sizeof(std::uint64_t) +
                          sampled_sa_.size() * sizeof(std::uint32_t) +
                          short_suffixes_.size() * sizeof(PackedShortSuffix);
  }

  // Checks everything a lookup relies on, so that a corrupt table is rejected
  // at load time and cannot make a lookup read out of bounds.
  void validate(const std::uint64_t entry_count,
                const std::uint64_t sampled_count) const {
    if (H_.back() != entry_count) {
      throw std::runtime_error("PT16 sassy H directory does not end at "
                               "entry_count");
    }

    if (entry_count == 0) {
      throw std::runtime_error("PT16 sassy table has no entries");
    }

    for (const std::uint32_t start : H_sa_) {
      if (start > sampled_count) {
        throw std::runtime_error(
            "PT16 sassy H_sa points past the end of sampled_sa");
      }
    }

    for (const PackedShortSuffix& suffix : short_suffixes_) {
      if (suffix.length < SHORT_SUFFIX_MIN_LENGTH ||
          suffix.length >= KMER_LENGTH) {
        throw std::runtime_error("PT16 sassy short suffix has a bad length");
      }

      // The characters occupy the top 2 * length bits; the rest must be 0.
      const std::uint32_t unused_bits = (1U << (32U - 2U * suffix.length)) - 1U;

      if ((suffix.packed & unused_bits) != 0) {
        throw std::runtime_error(
            "PT16 sassy short suffix has bits beyond its length");
      }
    }

    // An empty bucket must point at a non-empty one.
    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
      if (!(H_[bucket] & EMPTY_BUCKET_FLAG)) {
        continue;
      }

      const std::uint32_t matching = H_[bucket] & empty_bucket_mask;

      if (matching >= NUMBER_OF_BUCKETS || (H_[matching] & EMPTY_BUCKET_FLAG)) {
        throw std::runtime_error(
            "PT16 sassy empty bucket does not point at a non-empty bucket");
      }
    }

    // The non-empty buckets must tile L, and each bucket's entries must be
    // sorted by `low` with their slices inside sampled_sa.
    std::uint32_t expected_begin = 0;

    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
      if (H_[bucket] & EMPTY_BUCKET_FLAG) {
        continue;
      }

      const std::uint32_t begin = H_[bucket];
      const std::uint32_t end = bucket_end(bucket);

      if (begin != expected_begin || end <= begin || end > entry_count) {
        throw std::runtime_error("PT16 sassy H directory is inconsistent");
      }

      expected_begin = end;

      for (std::uint32_t i = begin; i < end; ++i) {
        const std::uint64_t entry = L_[i];

        if (i > begin &&
            sassy_decode_low(entry) <= sassy_decode_low(L_[i - 1])) {
          throw std::runtime_error(
              "PT16 sassy bucket entries are not sorted by low part");
        }

        if (sassy_is_range(entry)) {
          std::uint64_t start = static_cast<std::uint64_t>(H_sa_[bucket]) +
                                sassy_decode_offset(entry);
          const std::uint16_t count_field = sassy_decode_count_field(entry);

          std::uint64_t count = count_field;

          if (count_field == SASSY_COUNT_ESCAPE) {
            // The real count lives at sampled_sa[start]; it must be there to
            // read, and it must be why the entry escaped in the first place.
            if (start >= sampled_count) {
              throw std::runtime_error(
                  "PT16 sassy escaped range entry has no count marker");
            }

            count = sampled_sa_[start];
            ++start;

            if (count < SASSY_COUNT_ESCAPE) {
              throw std::runtime_error(
                  "PT16 sassy range entry escapes unnecessarily");
            }
          }

          if (count == 0 || start + count > sampled_count) {
            throw std::runtime_error(
                "PT16 sassy range entry lies outside sampled_sa");
          }
        }
      }
    }

    if (expected_begin != entry_count) {
      throw std::runtime_error("PT16 sassy H directory does not cover L");
    }
  }
};
