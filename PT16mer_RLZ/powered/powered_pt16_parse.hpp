#pragma once

// powered_rlz's parse_tuples (code size 4), with a PT16 lookup
// (pt16_powered.hpp) at every phrase start.
//
// powered parses right to left: a phrase starts at its right end i, with the
// interval of the metacharacter seq[i-4, i) taken from the C array, and then
// extends four characters per rank step. Here, while i >= 16, the 16-mer
// seq[i-16, i) is looked up first (sassy-layout table, pt16_powered.hpp):
//   singleton hit: the match is unique; it is extended leftwards by comparing
//         characters with the reference (see escape below);
//   range hit: the backward search continues from its row interval with 16
//         characters matched -- the state powered itself reaches after its
//         first three rank steps;
//   miss: the phrase is shorter than 16 and the table already gives it (its
//         length is the longest common suffix, see pt16_powered.hpp): it is
//         emitted with no backward search, and the next phrase starts.
// Below 16 characters the phrase starts as in powered.
//
// Escape on singleton (whenever the reference is given): as soon as the
// interval is a single row (b - a == 1) the match is unique, so instead of
// further rank steps the phrase is extended leftwards by comparing characters
// with the (cyclic) reference, as the SA and PT16 parsers do once their
// interval is a singleton.
//
// Variants: no table and no reference is powered itself; the reference alone
// is powered with the escape; the table needs the reference (a singleton
// entry has no row to go on ranking from).
//
// The phrase LENGTHS are identical to powered's (both are the greedy longest
// suffix); a miss's position may be another occurrence than the one powered
// reports. powered also emits a zero-length phrase when a phrase ends
// exactly at the start of the input; this parse does not.
//
// Otherwise it follows parse_tuples. It uses only public members of the
// index (rank, largest_extension, gca_, count); powered's C array
// (char_counts_, private) is rebuilt once from count().

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <vector>

#include "pt16_powered.hpp"

template <class Index>
class PoweredPT16Parser {
 public:
  using Phrase = std::tuple<std::uint64_t, std::uint64_t>;

  struct Stats {
    std::uint64_t lookups = 0;  // phrase starts with i >= 16
    std::uint64_t hits = 0;     // ... whose 16-mer occurs
    std::uint64_t misses = 0;   // ... resolved by the table alone
    std::uint64_t escapes = 0;  // phrases finished by character comparison
  };

  PoweredPT16Parser(const Index& index, const PT16PoweredTable* table,
                    const std::vector<unsigned char>* reference)
      : index_(index), table_(table), reference_(reference) {
    // C[c] = rows whose rotation starts with a metacharacter smaller than c;
    // count(m, 4, 4) of a 4-character m is the size of m's interval.
    const char dna[4] = {'A', 'C', 'G', 'T'};
    C_[0] = 0;
    for (unsigned c = 0; c < 256; ++c) {
      const char meta[4] = {dna[(c >> 6) & 3U], dna[(c >> 4) & 3U],
                            dna[(c >> 2) & 3U], dna[c & 3U]};
      C_[c + 1] = C_[c] + index.count(std::string_view(meta, 4), 4, 4);
    }
    if (C_[256] != index.size()) {
      throw std::runtime_error("rebuilt C array does not sum to the index size");
    }
    if (reference != nullptr && reference->size() != index.size()) {
      throw std::runtime_error("reference and powered index differ in size");
    }
    if (table != nullptr && reference == nullptr) {
      throw std::runtime_error("the PT16 table needs the reference");
    }
    if (table != nullptr && table->rows() != index.size()) {
      throw std::runtime_error("PT16 table and powered index differ in size");
    }
  }

  void parse(const std::string_view& seq, std::vector<Phrase>& parses,
             Stats& stats) const {
    constexpr std::uint8_t code_size = 4;
    auto constexpr encode_base = [](unsigned char ch) noexcept -> std::uint8_t {
      return static_cast<std::uint8_t>((ch > 'A') + (ch > 'C') + (ch > 'G'));
    };
    const auto* text = reinterpret_cast<const unsigned char*>(seq.data());
    const auto& gca = index_.gca_;
    const std::uint64_t size = index_.size();

    std::uint8_t c = 0;
    std::pair<std::uint64_t, std::uint64_t> tuple;

    // case if sequence is too short
    if (seq.size() <= code_size) {
      for (std::uint8_t j = 0; j < seq.size(); j++) {
        c <<= 2;
        c += encode_base(seq[j]);
      }
      c <<= ((code_size - seq.size()) * 2);

      while (C_[c + 1] == C_[c]) c++;
      parses.emplace_back(seq.size(), gca[C_[c]]);
      return;
    }

    std::uint64_t length = 0;
    std::uint64_t a = 0;
    std::uint64_t b = 0;

    // Emits the phrase whose match seq[i, i + matched) occurs only at
    // reference position r, extended leftwards character by character, and
    // moves i to its left end.
    const unsigned char* ref =
        reference_ != nullptr ? reference_->data() : nullptr;
    const auto extend_unique = [&](std::size_t& i, std::uint64_t r,
                                   const std::uint64_t matched) {
      std::size_t j = i;
      while (j > 0) {
        const std::uint64_t previous = (r == 0 ? size : r) - 1;
        if (ref[previous] != text[j - 1]) break;
        r = previous;
        --j;
      }
      ++stats.escapes;
      parses.emplace_back(matched + (i - j), r);
      i = j;
    };

    // Starts the phrase whose right end is i. Returns false if fewer than
    // code_size characters remain (the tail below parses them); otherwise
    // [a, b) is the interval of the matched seq[i, i + length) after i was
    // moved to its left end.
    const auto start_phrase = [&](std::size_t& i) -> bool {
      while (table_ != nullptr && i >= KMER_LENGTH) {
        ++stats.lookups;
        const PoweredStart s = table_->lookup(text + i - KMER_LENGTH);
        if (s.found && s.singleton) {
          ++stats.hits;
          i -= KMER_LENGTH;
          extend_unique(i, s.position, KMER_LENGTH);
          continue;
        }
        if (s.found) {
          ++stats.hits;
          a = s.a;
          b = s.b;
          length = KMER_LENGTH;
          i -= KMER_LENGTH;
          return true;
        }
        if (s.match_length == 0) break;  // (only if a base is absent)
        ++stats.misses;
        parses.emplace_back(s.match_length, s.position);
        i -= s.match_length;
      }
      if (i < code_size) return false;

      // as powered: the interval of the metacharacter seq[i-4, i)
      c = encode_base(seq[i - 4]) << 6 | encode_base(seq[i - 3]) << 4 |
          encode_base(seq[i - 2]) << 2 | encode_base(seq[i - 1]);
      length = code_size;
      a = C_[c];
      b = C_[c + 1];
      i -= code_size;
      return true;
    };

    std::size_t i = seq.size();
    if (!start_phrase(i)) length = 0;

    while (length > 0 && i >= code_size) {
      // escape on singleton: the match seq[i, i + length) occurs only at
      // p = gca[a]; extend it leftwards character by character
      if (ref != nullptr && b - a == 1) {
        extend_unique(i, gca[a], length);

        if (!start_phrase(i)) {
          length = 0;
          break;
        }
        continue;
      }

      // load the next metacharacter
      c = encode_base(seq[i - 4]) << 6 | encode_base(seq[i - 3]) << 4 |
          encode_base(seq[i - 2]) << 2 | encode_base(seq[i - 1]);

      // compute the rank of this metacharacter
      const std::uint64_t a_tmp = index_.rank(a, c);
      const std::uint64_t b_tmp = index_.rank(b, c);

      // check if the suffix can be extended by a whole metacharacter, else
      // extend by a partial metacharacter and start the next phrase
      if (b_tmp <= a_tmp) {
        tuple = index_.largest_extension(a, b - 1, c);
        parses.emplace_back(length + tuple.first, tuple.second);
        i -= tuple.first;  // < 4 <= i

        if (!start_phrase(i)) {
          length = 0;
          break;
        }
        continue;
      }

      // continue with extension
      a = a_tmp + C_[c];
      b = a - a_tmp + b_tmp;

      // update length
      length += code_size;
      i -= code_size;
    }

    // check if the loop covered the whole sequence or there is a partial
    // metacharacter left to cover
    if (i > 0) {
      c = 0;
      for (std::uint8_t j = 0; j < i; j++) {
        c += encode_base(seq[i - j - 1]) << (j << 1);
      }

      // extend the suffix
      if (length > 0) {
        tuple = index_.largest_extension(a, b - 1, c);
        if (tuple.first >= i) {
          parses.emplace_back(length + i,
                              (tuple.second + (tuple.first - i)) % size);
          return;
        }
        parses.emplace_back(length, gca[a]);
      }
      c <<= ((code_size - i) << 1);

      while (C_[c + 1] == C_[c]) c++;
      parses.emplace_back(i, gca[C_[c]]);
    } else if (length > 0) {
      parses.emplace_back(length, gca[a]);
    }
  }

 private:
  const Index& index_;
  const PT16PoweredTable* table_;
  const std::vector<unsigned char>* reference_;
  std::array<std::uint64_t, 257> C_{};
};
