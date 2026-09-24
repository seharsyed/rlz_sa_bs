#pragma once

#include <cstdint>
#include <vector>

#include "../pt16_utils.hpp"  // KMER_LENGTH, KmerLookupResult, MatchingStatistics

/**
 * Turns one 16-mer lookup per input position (in text order) into
 * matching statistics, by extending a singleton hit's match by one
 * character whenever the immediately following position's own lookup
 * result -- singleton or range -- contains this position's occurrence + 1.
 *
 * Why comparing positions is enough (no reference read needed): if the
 * 16-mer at i occurs at reference position p, and the 16-mer at i+1
 * occurs at p+1, then input[i..i+16) == ref[p..p+16) (the first fact) and
 * in particular input[i+16] == ref[p+16] (the last character of the
 * second fact), so input[i..i+17) == ref[p..p+17): the match at i is at
 * least 17 long. This holds regardless of how many OTHER occurrences
 * either 16-mer has, so the check also looks inside a range, not just a
 * singleton -- but only the CURRENT position needs to be a singleton for
 * this to fire (a range's own several candidate positions at i are not
 * resolved to one, so nothing here can be said about which of them, if
 * any, extends).
 *
 * One-step, single-position lookahead only -- a chain longer than 17 is
 * still extended by just one character here; see ms_tools.hpp's
 * backwardChainExtend for the multi-step version that chases the whole
 * chain and is proven exact against brute force. Generalized to work on
 * either table format's lookup results: both PT16SassyLookup::
 * lookup and PT16RLZParser::lookupKmerByKey return the same
 * KmerLookupResult type (pt16_utils.hpp), so this function never touches
 * a lookup table, the reference, or a suffix array itself -- it only
 * reads `results`, whichever format built it.
 *
 * `results[i]` must be position i's own lookup result, for every i in
 * [0, results.size()): a real 16-mer lookup where 16 or more characters
 * remain, or (found=false, match_length/match_position already set to
 * the short/tail answer) where fewer than 16 remain. A tail position's
 * answer must be computed separately (PT16SassyLookup::lookup_tail /
 * PT16RLZParser::lookupKmerOrTail's own tail branch) and placed into
 * `results` the same way -- this function cannot compute one itself,
 * since doing so needs the input text, which it never sees; it works
 * purely from what a caller already looked up.
 */
inline MatchingStatistics chainExtend(
    const std::vector<KmerLookupResult>& results) {
  const std::size_t n = results.size();
  MatchingStatistics ms(n);

  for (std::size_t i = 0; i < n; ++i) {
    const KmerLookupResult& entry = results[i];

    // A hit's match_length is exactly 16 before extension (that is what
    // "found" means here); a miss/tail's own match_length is already the
    // final answer, whatever produced `results` computed it.
    std::uint32_t match_length = entry.found ? KMER_LENGTH : entry.match_length;

    if (entry.found && entry.count == 1 && i + 1 < n) {
      const std::uint32_t wanted = entry.match_position + 1;
      const KmerLookupResult& next = results[i + 1];
      bool extends = false;

      if (next.found) {
        if (next.count == 1) {
          extends = next.match_position == wanted;
        } else {
          for (const std::uint32_t candidate : next.positions) {
            if (candidate == wanted) {
              extends = true;
              break;
            }
          }
        }
      }

      if (extends) {
        ++match_length;
      }
    }

    ms[i] = {entry.match_position, match_length};
  }

  return ms;
}
