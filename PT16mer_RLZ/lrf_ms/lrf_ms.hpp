#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "rmq_tree.h"

using MatchingStatistics = std::vector<std::pair<std::uint32_t, std::uint32_t>>;

/**
 * LRF-MS matching statistics.
 *
 * The reference sequence and suffix array are supplied by the caller.
 * No files are read and no suffix array is constructed here.
 *
 * The constructor builds:
 *   - ISA
 *   - LCP
 *   - LRF
 *   - RMQ over LCP
 *
 * The input sequence is a single sequence. There are no sequence
 * separators and no sentinel is added to either sequence.
 *
 * REQUIREMENTS ON THE CALLER
 *
 *   - suffix_array must be the suffix array of `reference` as given,
 *     i.e. built WITHOUT a terminating sentinel, using the usual
 *     convention that a suffix which is a proper prefix of another
 *     suffix sorts before it. If you build the SA with libsais over a
 *     sentinel-terminated string, drop the sentinel entry and pass the
 *     reference without it, or sizes will not match and the
 *     constructor will throw.
 *
 *   - The reference and the suffix array must outlive this object.
 *     Both are held by reference; binding a temporary is rejected at
 *     compile time by the deleted rvalue overload below.
 *
 * No assumption is made about the alphabet: suffixes shorter than the
 * current comparison depth are treated as sorting before every real
 * symbol, which is exactly what a sentinel-free suffix array encodes.
 * Symbol value 0 is therefore an ordinary symbol here.
 *
 * NOTE ON THE "NO MATCH" POSITION
 *
 * A factor of length 0 is reported with reference position
 * reference.size() - 1. In the original implementation that position
 * was the sentinel; here it is a real reference position. Downstream
 * code must use the length, not the position, to detect an empty
 * match.
 */

template <typename T1, typename T2>
class LRFMS {
 public:
  using Symbol = T1;
  using Index = std::int32_t;
  using input_type = std::vector<T1>;
  using reference_type = std::vector<T1>;
  using suffix_array_type = std::vector<T2>;

  LRFMS(const reference_type& reference, const suffix_array_type& suffix_array)
      : reference_(reference),
        suffix_array_(suffix_array),
        isa_(reference.size()),
        lcp_(reference.size()),
        lrf_(reference.size()) {
    if (reference_.empty()) {
      throw std::invalid_argument("LRFMS: reference is empty");
    }

    if (suffix_array_.size() != reference_.size()) {
      throw std::invalid_argument(
          "LRFMS: suffix array size does not match reference size");
    }

    if (reference_.size() >
        static_cast<std::size_t>((std::numeric_limits<Index>::max)())) {
      throw std::invalid_argument("LRFMS: reference is too large for Index");
    }

    constructISA();
    constructLCP();
    constructLRF();

    rmq_ = std::make_unique<rmq_tree<Index>>(lcp_.data(),
                                             static_cast<int>(lcp_.size()), 7);
  }

  // The reference and the suffix array are held by reference; temporaries
  // would dangle.
  LRFMS(reference_type&&, const suffix_array_type&) = delete;
  LRFMS(const reference_type&, suffix_array_type&&) = delete;
  LRFMS(reference_type&&, suffix_array_type&&) = delete;

  LRFMS(const LRFMS&) = delete;
  LRFMS& operator=(const LRFMS&) = delete;

  LRFMS(LRFMS&&) = delete;
  LRFMS& operator=(LRFMS&&) = delete;

  /**
   * Compute matching statistics for one input sequence.
   *
   * Each result is:
   *   (reference_position, match_length)
   *
   * The result at position i corresponds to input_sequence[i].
   */
  MatchingStatistics computeMatchingStatistics(
      const input_type& input_sequence) const {
    MatchingStatistics matching_statistics;
    matching_statistics.reserve(input_sequence.size());

    if (input_sequence.empty()) {
      return matching_statistics;
    }

    const Index n = static_cast<Index>(reference_.size());

    Index leftB = 0;
    Index rightB = n - 1;
    Index pos = n - 1;
    Index len = 0;

    std::size_t i = 0;

    while (i < input_sequence.size()) {
      computeMatchingFactor(input_sequence, i, &pos, &len, leftB, rightB);

      matching_statistics.emplace_back(static_cast<std::size_t>(pos),
                                       static_cast<std::size_t>(len));

      /*
       * Same left-extension step as the original implementation: after
       * emitting the factor beginning at i, the next iteration starts at
       * i + 1 while the current match is contracted by one.
       */
      --len;

      /*
       * len reaches 0 when the input symbol does not occur in the
       * reference at all, and the contraction above then makes it -1.
       * The original single-threaded routine fed that -1 straight back
       * into the next call, where it was read as an unsigned offset;
       * the CMS variant clamped it. We clamp.
       */
      if (len < 0) {
        len = 0;
      }

      if (leftB == rightB) {
        /*
         * LRF skip. The invariant is that on entry to each iteration
         * `len` is the matching statistic of input position i + 1, so
         * len >= 1 implies i + 1 is a valid input position and the loop
         * cannot run past the end of the input. The explicit bounds
         * checks below are defensive only and never change the result.
         */
        while (pos + 1 < n && len > lrf_[static_cast<std::size_t>(pos) + 1] &&
               i + 1 < input_sequence.size()) {
          matching_statistics.emplace_back(static_cast<std::size_t>(pos + 1),
                                           static_cast<std::size_t>(len));

          ++i;
          --len;
          ++pos;
        }

        if (i < input_sequence.size() && pos + 1 < n) {
          const auto interval =
              adjustInterval(isa_[static_cast<std::size_t>(pos) + 1],
                             isa_[static_cast<std::size_t>(pos) + 1], len);

          leftB = interval.first;
          rightB = interval.second;
        } else {
          // We reached the end of the reference suffix: back to the root.
          leftB = 0;
          rightB = n - 1;
        }
      } else {
        const auto interval = contractLeft(leftB, rightB, len);
        leftB = interval.first;
        rightB = interval.second;
      }

      ++i;
    }

    return matching_statistics;
  }

  const std::vector<Index>& isa() const noexcept { return isa_; }

  const std::vector<Index>& lcp() const noexcept { return lcp_; }

  const std::vector<Index>& lrf() const noexcept { return lrf_; }

 private:
  const reference_type& reference_;
  const suffix_array_type& suffix_array_;

  std::vector<Index> isa_;
  std::vector<Index> lcp_;
  std::vector<Index> lrf_;

  std::unique_ptr<rmq_tree<Index>> rmq_;

  void constructISA() {
    for (std::size_t i = 0; i < suffix_array_.size(); ++i) {
      const auto suffix = suffix_array_[i];

      if (suffix < 0 ||
          static_cast<std::size_t>(suffix) >= reference_.size()) {
        throw std::invalid_argument(
            "LRFMS: suffix array contains an invalid position");
      }

      isa_[static_cast<std::size_t>(suffix)] = static_cast<Index>(i);
    }
  }

  /**
   * Kasai LCP construction.
   *
   * lcp_[rank] is the LCP of SA[rank] and SA[rank - 1]; lcp_[0] is zero.
   * This matches the convention produced by libsais_plcp + libsais_lcp.
   */
  void constructLCP() {
    std::size_t h = 0;

    lcp_[0] = 0;

    for (std::size_t i = 0; i < reference_.size(); ++i) {
      const std::size_t rank = static_cast<std::size_t>(isa_[i]);

      if (rank == 0) {
        h = 0;
        continue;
      }

      const std::size_t j = static_cast<std::size_t>(suffix_array_[rank - 1]);

      while (i + h < reference_.size() && j + h < reference_.size() &&
             reference_[i + h] == reference_[j + h]) {
        ++h;
      }

      lcp_[rank] = static_cast<Index>(h);

      if (h > 0) {
        --h;
      }
    }
  }

  /**
   * LRF[i] = max(LCP[ISA[i]], LCP[ISA[i] + 1])
   *
   * The original code read one past the end of the LCP array at the
   * largest rank because it relied on the sentinel suffix being present.
   * Here the out-of-range neighbour contributes 0, which is the correct
   * value.
   */
  void constructLRF() {
    for (std::size_t i = 0; i < reference_.size(); ++i) {
      const std::size_t rank = static_cast<std::size_t>(isa_[i]);

      Index left_lcp = 0;
      Index right_lcp = 0;

      if (rank > 0) {
        left_lcp = lcp_[rank];
      }

      if (rank + 1 < lcp_.size()) {
        right_lcp = lcp_[rank + 1];
      }

      lrf_[i] = std::max(left_lcp, right_lcp);
    }
  }

  std::pair<Index, Index> adjustInterval(Index lo, Index hi,
                                         Index offset) const {
    Index psv = rmq_->psv(lo, offset);

    if (psv == -1) {
      psv = 0;
    }

    Index nsv;

    if (static_cast<std::size_t>(hi) + 1 >= lcp_.size()) {
      nsv = static_cast<Index>(reference_.size()) - 1;
    } else {
      nsv = rmq_->nsv(hi + 1, offset);

      if (nsv == -1) {
        nsv = static_cast<Index>(reference_.size()) - 1;
      } else {
        --nsv;
      }
    }

    return {psv, nsv};
  }

  std::pair<Index, Index> contractLeft(Index lo, Index hi, Index offset) const {
    const std::size_t suflo = static_cast<std::size_t>(suffix_array_[lo]);
    const std::size_t sufhi = static_cast<std::size_t>(suffix_array_[hi]);

    /*
     * Without a reference sentinel, reaching the last reference position
     * means the suffix cannot be contracted further, so we fall back to
     * the root. The original hit this case only for the sentinel suffix;
     * here it is reachable for the genuine last position as well.
     */
    if (suflo + 1 >= reference_.size() || sufhi + 1 >= reference_.size()) {
      return {0, static_cast<Index>(reference_.size()) - 1};
    }

    const Index tmplo = isa_[suflo + 1];
    const Index tmphi = isa_[sufhi + 1];

    return adjustInterval(tmplo, tmphi, offset);
  }

  /**
   * Symbol of the suffix at rank `rank`, `offset` characters in.
   *
   * Returns false if that suffix is shorter than `offset`, in which case
   * it is a proper prefix of every suffix it is being compared against
   * and therefore sorts before all of them. This is what replaces the
   * sentinel of the original implementation, and it holds for any
   * alphabet.
   */
  inline bool symbolAt(Index rank, std::size_t offset, Symbol* out) const {
    const std::size_t suffix = static_cast<std::size_t>(suffix_array_[rank]);

    if (suffix + offset >= reference_.size()) {
      return false;
    }

    *out = reference_[suffix + offset];
    return true;
  }

  /**
   * Leftmost occurrence of c at depth `offset` within [lo, hi], or
   * -(x + 1) where x is the insertion point, matching the original.
   */
  inline Index binarySearchLB(Index lo, Index hi, std::size_t offset,
                              Symbol c) const {
    Index low = lo;
    Index high = hi;

    while (low <= high) {
      const Index mid = low + ((high - low) >> 1);

      Symbol midVal;
      const bool midInRange = symbolAt(mid, offset, &midVal);

      if (!midInRange || midVal < c) {
        low = mid + 1;
      } else if (midVal > c) {
        high = mid - 1;
      } else {
        if (mid == lo) {
          return mid;
        }

        Symbol midValLeft;
        const bool leftInRange = symbolAt(mid - 1, offset, &midValLeft);

        if (leftInRange && midValLeft == midVal) {
          high = mid - 1;
        } else {
          return mid;
        }
      }
    }

    return -(low + 1);
  }

  inline Index binarySearchRB(Index lo, Index hi, std::size_t offset,
                              Symbol c) const {
    Index low = lo;
    Index high = hi;

    while (low <= high) {
      const Index mid = low + ((high - low) >> 1);

      Symbol midVal;
      const bool midInRange = symbolAt(mid, offset, &midVal);

      if (!midInRange || midVal < c) {
        low = mid + 1;
      } else if (midVal > c) {
        high = mid - 1;
      } else {
        if (mid == hi) {
          return mid;
        }

        Symbol midValRight;
        const bool rightInRange = symbolAt(mid + 1, offset, &midValRight);

        if (rightInRange && midValRight == midVal) {
          low = mid + 1;
        } else {
          return mid;
        }
      }
    }

    return -(low + 1);
  }

  /**
   * Scans the input from i + *len onwards until the maximal prefix
   * shared with the reference is found, updating the SA interval.
   *
   * On return *pos is the reference position of the match and *len its
   * length. This mirrors the original routine one-for-one, including the
   * depth `offset` advancing with every matched character.
   */
  void computeMatchingFactor(const input_type& input_sequence, std::size_t i,
                             Index* pos, Index* len, Index& leftB,
                             Index& rightB) const {
    std::size_t offset = static_cast<std::size_t>(*len);
    std::size_t j = i + offset;

    Index nlb = leftB;
    Index nrb = rightB;

    Index match = suffix_array_[nlb];

    while (j < input_sequence.size()) {
      if (nlb == nrb) {
        Symbol midVal;

        // A suffix shorter than the current depth cannot match.
        if (!symbolAt(nlb, offset, &midVal) || midVal != input_sequence[j]) {
          break;
        }

        leftB = nlb;
        rightB = nrb;
      } else {
        // Refine the bucket holding the match, from the left and then
        // from the right.
        nlb = binarySearchLB(nlb, nrb, offset, input_sequence[j]);

        if (nlb < 0) {
          Index maxMatch = -nlb - 1;

          if (maxMatch == nrb + 1) {
            --maxMatch;
          }

          match = suffix_array_[maxMatch];
          break;
        }

        nrb = binarySearchRB(nlb, nrb, offset, input_sequence[j]);

        leftB = nlb;
        rightB = nrb;
      }

      match = suffix_array_[nlb];

      ++j;
      ++offset;
    }

    *pos = match;
    *len = static_cast<Index>(offset);
  }
};
