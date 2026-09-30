#pragma once

// RLZ_powered's parse_tuples (code size 4), with escape on small intervals,
// matches only inside the linear reference, and an optional PT16 lookup
// (pt16_powered.hpp) at every phrase start.
//
// powered parses right to left: a phrase starts at its right end i, with the
// interval of the metacharacter seq[i-4, i) taken from the C array, and then
// extends four characters per rank step.
//
// Linear matches only. powered's index is cyclic (its rows are rotations of
// the reference), so a match may run over the end of the reference into its
// start -- a match the linear parsers (sa-binary-search, the PT16 tables,
// Varki) cannot make. This parse accepts only occurrences inside the
// reference, so its phrases are those of a linear index:
//   - the occurrences of a match of length L that run over the end start in
//     the last L-1 positions, so an interval of at least L rows always holds
//     an occurrence inside the reference: powered's rank steps go on;
//   - an interval of fewer than L rows (including the singleton the escape
//     used to handle) is resolved by comparing characters: for each row whose
//     occurrence lies inside the reference, the match is extended leftwards
//     with the reference (stopping at its start), and the longest wins -- the
//     exact longest linear match (the SA and PT16 parsers do the same once
//     their interval is small);
//   - if no row of such an interval lies inside the reference, the match of
//     the previous step (which does) is resolved that way instead; a partial
//     extension (largest_extension) or a table answer that points over the
//     end is resolved the same way, or by starting the phrase again without
//     the table.
//
// The table: while i >= 16, the 16-mer seq[i-16, i) is looked up first:
//   singleton hit: its one occurrence; extended by comparing characters;
//   range hit: the backward search continues from its row interval with 16
//         characters matched -- the state powered itself reaches after its
//         first three rank steps;
//   miss: the phrase is shorter than 16 and the table already gives it (its
//         length is the longest common suffix, see pt16_powered.hpp): it is
//         emitted with no backward search, and the next phrase starts.
// Below 16 characters, or without a table, the phrase starts as in powered.
//
// Forward mode (forward = true): the index (and the table, and the reference
// given here) are of the REVERSED reference R^r, and the input S is read
// reversed -- virtually, through at(k) = S[n-1-k], no copy. powered's
// right-to-left parse of S^r against R^r is then the greedy LEFT-TO-RIGHT
// parse of S against R (a longest suffix of a prefix of S^r is the reverse of
// a longest prefix of a suffix of S): with linear matches only, exactly the
// phrases of sa-binary-search and Varki. Positions are converted back to R,
// n - p - len, and the phrases come out in input order (not last to first).
//
// The phrase lengths are those of the greedy longest linear match (right to
// left; left to right in forward mode); a position may be another occurrence
// than the one powered reports. Unlike powered, the tail is greedy too
// (powered drops a partial extension before its last short phrase), and no
// zero-length phrase is emitted.
//
// It uses only public members of the index (rank, largest_extension, sa_,
// count); powered's C array (char_counts_, private) is rebuilt once from
// count(). The reference is required (the character comparisons).

#include <algorithm>
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
    std::uint64_t lookups = 0;    // phrase starts with i >= 16
    std::uint64_t hits = 0;       // ... whose 16-mer occurs
    std::uint64_t misses = 0;     // ... resolved by the table alone
    std::uint64_t escapes = 0;    // phrases finished by character comparison
    std::uint64_t fallbacks = 0;  // matches that only ran over the end
  };

  PoweredPT16Parser(const Index& index, const PT16PoweredTable* table,
                    const std::vector<unsigned char>* reference,
                    const bool forward = false)
      : index_(index), table_(table), reference_(reference), forward_(forward) {
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
    if (reference == nullptr) {
      throw std::runtime_error("the powered parse needs the reference");
    }
    if (reference->size() != index.size()) {
      throw std::runtime_error("reference and powered index differ in size");
    }
    if (table != nullptr && table->rows() != index.size()) {
      throw std::runtime_error("PT16 table and powered index differ in size");
    }
  }

  // Phrases (length, position): last to first normally, in input order in
  // forward mode (see above).
  void parse(const std::string_view& seq, std::vector<Phrase>& parses,
             Stats& stats) const {
    if (forward_) {
      parse_direction<true>(seq, parses, stats);
    } else {
      parse_direction<false>(seq, parses, stats);
    }
  }

  bool forward() const { return forward_; }

 private:
  template <bool Forward>
  void parse_direction(const std::string_view& seq, std::vector<Phrase>& parses,
                       Stats& stats) const {
    constexpr std::uint8_t code_size = 4;
    auto constexpr encode_base = [](unsigned char ch) noexcept -> std::uint8_t {
      return static_cast<std::uint8_t>((ch > 'A') + (ch > 'C') + (ch > 'G'));
    };
    const auto* text = reinterpret_cast<const unsigned char*>(seq.data());
    const std::size_t n = seq.size();
    const auto& gca = index_.sa_;
    const std::uint64_t size = index_.size();
    const unsigned char* ref = reference_->data();

    // The input as the parse sees it: S, or S reversed (forward mode).
    const auto at = [&](const std::size_t k) -> unsigned char {
      if constexpr (Forward) {
        return text[n - 1 - k];
      } else {
        return text[k];
      }
    };
    // Emits a phrase found at position p of the indexed text (R, or R^r in
    // forward mode, where it is converted back to R). p + len <= size.
    const auto emit = [&](const std::uint64_t len, const std::uint64_t p) {
      if constexpr (Forward) {
        parses.emplace_back(len, size - p - len);
      } else {
        parses.emplace_back(len, p);
      }
    };
    // The table lookup of the 16-mer seen at [i-16, i).
    const auto lookup = [&](const std::size_t i) {
      if constexpr (Forward) {
        // Its reversed key is the forward key of S[n-i, n-i+16).
        std::uint32_t key = 0;
        for (std::size_t j = n - i; j < n - i + KMER_LENGTH; ++j) {
          key = (key << 2U) | alphatab[text[j]];
        }
        return table_->lookup_key(key);
      } else {
        return table_->lookup(text + i - KMER_LENGTH);
      }
    };
    const auto metacharacter = [&](const std::size_t i) -> std::uint8_t {
      return static_cast<std::uint8_t>(
          encode_base(at(i - 4)) << 6 | encode_base(at(i - 3)) << 4 |
          encode_base(at(i - 2)) << 2 | encode_base(at(i - 1)));
    };

    // ---------- Resolving a match by comparing characters ----------

    // The longest leftward extension of the match seq[i, i + len) over the
    // occurrences at rows [lo, hi) that lie inside the reference.
    struct Best {
      bool found = false;
      std::uint64_t ext = 0;
      std::uint64_t p = 0;  // where the match (unextended) occurs
    };
    const auto extend_at = [&](const std::uint64_t p, const std::size_t i) {
      std::uint64_t k = 0;
      while (k < i && k < p && ref[p - 1 - k] == at(i - 1 - k)) ++k;
      return k;
    };
    const auto best_extension = [&](const std::uint64_t lo, const std::uint64_t hi,
                                    const std::uint64_t len, const std::size_t i) {
      Best best;
      for (std::uint64_t row = lo; row < hi; ++row) {
        const std::uint64_t p = gca[row];
        if (p + len > size) continue;  // runs over the reference's end
        const std::uint64_t k = extend_at(p, i);
        if (!best.found || k > best.ext) best = {true, k, p};
      }
      return best;
    };
    // Emits the match seq[i, i + len) extended as `best` says; i moves to
    // the phrase's left end.
    const auto emit_best = [&](std::size_t& i, const std::uint64_t len,
                               const Best& best) {
      ++stats.escapes;
      emit(len + best.ext, best.p - best.ext);
      i -= best.ext;
    };
    // The longest suffix of seq[0, i) of at most min(i, 3) characters that
    // occurs inside the reference (the rows of a 1-3 character prefix are
    // the C range of the metacharacters starting with it).
    const auto short_phrase = [&](std::size_t& i) {
      for (std::uint64_t len = std::min<std::size_t>(i, 3); len >= 1; --len) {
        unsigned prefix = 0;
        for (std::uint64_t j = 0; j < len; ++j) {
          prefix = (prefix << 2U) | encode_base(at(i - len + j));
        }
        const unsigned shift = 2U * static_cast<unsigned>(code_size - len);
        for (std::uint64_t row = C_[prefix << shift];
             row < C_[(prefix + 1) << shift]; ++row) {
          if (gca[row] + len <= size) {
            emit(len, gca[row]);
            i -= len;
            return;
          }
        }
      }
      // A character that does not occur in the reference: a literal, as
      // powered_rlz's decompress reads it (the length's top bit set).
      parses.emplace_back((std::uint64_t{1} << 63) | 1, at(i - 1));
      i -= 1;
    };

    // ---------- The search state ----------

    std::uint64_t length = 0;  // characters matched: seq[i, i + length)
    std::uint64_t a = 0, b = 0;  // their row interval
    // The previous step's state (its match has an occurrence inside the
    // reference), for a step whose rows all run over the end.
    bool have_previous = false;
    std::uint64_t previous_a = 0, previous_b = 0, previous_length = 0;
    std::size_t previous_i = 0;
    std::size_t phrase_end = 0;  // the right end of the current phrase
    bool from_table = false;     // the current match came from a range hit

    // Starts phrases at the right end i until one needs the backward search
    // (returns true, with its state set) or the input is parsed (false).
    // Phrases the table or the comparisons resolve are emitted here.
    const auto start_phrase = [&](std::size_t& i, bool use_table) -> bool {
      while (i > 0) {
        have_previous = false;
        phrase_end = i;

        if (use_table && table_ != nullptr && i >= KMER_LENGTH) {
          ++stats.lookups;
          const PoweredStart s = lookup(i);
          if (s.found && s.singleton) {
            ++stats.hits;
            if (s.position + KMER_LENGTH <= size) {
              i -= KMER_LENGTH;
              emit_best(i, KMER_LENGTH, Best{true, extend_at(s.position, i),
                                             s.position});
              continue;
            }
          } else if (s.found) {
            ++stats.hits;
            a = s.a;
            b = s.b;
            length = KMER_LENGTH;
            i -= KMER_LENGTH;
            from_table = true;
            return true;
          } else if (s.match_length != 0 &&
                     s.position + s.match_length <= size) {
            ++stats.misses;
            emit(s.match_length, s.position);
            i -= s.match_length;
            continue;
          }
          // The answer runs over the reference's end: this phrase again,
          // without the table.
          ++stats.fallbacks;
        }
        use_table = true;  // (for the next phrase)

        if (i < code_size) {
          short_phrase(i);
          continue;
        }
        // as powered: the interval of the metacharacter seq[i-4, i)
        const std::uint8_t c = metacharacter(i);
        if (C_[c + 1] == C_[c]) {  // (a metacharacter absent from the reference)
          short_phrase(i);
          continue;
        }
        a = C_[c];
        b = C_[c + 1];
        length = code_size;
        i -= code_size;
        from_table = false;
        return true;
      }
      return false;
    };

    // Resolves a small interval (fewer rows than characters matched): emits
    // the phrase and starts the next one. Returns start_phrase's result.
    const auto resolve = [&](std::size_t& i) -> bool {
      Best best = best_extension(a, b, length, i);
      if (!best.found) {
        // Every occurrence runs over the end.
        ++stats.fallbacks;
        if (have_previous) {
          i = previous_i;
          length = previous_length;
          best = best_extension(previous_a, previous_b, length, i);
        } else {
          i = phrase_end;
          if (from_table) return start_phrase(i, false);
          short_phrase(i);
          return start_phrase(i, true);
        }
      }
      emit_best(i, length, best);
      return start_phrase(i, true);
    };

    // ---------- The parse ----------

    std::size_t i = n;
    bool active = start_phrase(i, true);

    while (active) {
      if (b - a < length) {
        active = resolve(i);
        continue;
      }
      // Here the interval holds an occurrence inside the reference.
      if (i < code_size) break;

      const std::uint8_t c = metacharacter(i);
      const std::uint64_t a_tmp = index_.rank(a, c);
      const std::uint64_t b_tmp = index_.rank(b, c);

      if (b_tmp <= a_tmp) {
        // No whole metacharacter extends the match: a partial one, then the
        // next phrase. largest_extension's occurrence must lie inside the
        // reference; otherwise compare characters over the interval.
        const std::pair<std::uint64_t, std::uint64_t> tuple =
            index_.largest_extension(a, b - 1, c);
        if (tuple.second + length + tuple.first <= size) {
          emit(length + tuple.first, tuple.second);
          i -= tuple.first;  // < 4 <= i
        } else {
          ++stats.fallbacks;
          emit_best(i, length, best_extension(a, b, length, i));
        }
        active = start_phrase(i, true);
        continue;
      }

      have_previous = true;
      previous_a = a;
      previous_b = b;
      previous_length = length;
      previous_i = i;

      a = a_tmp + C_[c];
      b = a - a_tmp + b_tmp;
      length += code_size;
      i -= code_size;
    }

    // The tail: fewer than four characters left of a match whose interval
    // holds an occurrence inside the reference. Extend it greedily by
    // comparing characters, then parse what is left in short phrases.
    if (active) {
      emit_best(i, length, best_extension(a, b, length, i));
      while (i > 0) short_phrase(i);
    }
  }

  const Index& index_;
  const PT16PoweredTable* table_;
  const std::vector<unsigned char>* reference_;
  bool forward_;
  std::array<std::uint64_t, 257> C_{};
};
