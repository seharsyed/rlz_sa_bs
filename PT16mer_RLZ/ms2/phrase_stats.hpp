#pragma once

// Statistics of an input's matching statistics (e.g. lrf-ms's): over the
// positions with a match (MS >= 1; a separator has MS = 0 and is counted
// apart), the average MS, the share of positions whose MS is below 16, the
// maximum and the median. Summed over files, with the per-file average of
// the below-16 share.

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../pt16_utils.hpp"  // MatchingStatistics, KMER_LENGTH

namespace ms2 {

struct PhraseStats {
  std::size_t files = 0;
  std::size_t positions = 0;          // all input positions
  std::size_t matched_positions = 0;  // MS >= 1
  std::size_t unmatched = 0;          // MS = 0 (separators)
  std::uint64_t match_length_sum = 0;
  std::size_t short_ms_positions = 0;  // 1 <= MS < 16
  std::uint32_t max_length = 0;

  // histogram[L]: matched positions with MS = L, for the median.
  std::vector<std::uint64_t> histogram;

  // Sum over files of each file's below-16 share (per-file average).
  double short_ms_share_sum = 0.0;

  double average_match_length() const {
    return matched_positions == 0
               ? 0.0
               : static_cast<double>(match_length_sum) /
                     static_cast<double>(matched_positions);
  }

  double short_ms_percent() const {
    return matched_positions == 0
               ? 0.0
               : 100.0 * static_cast<double>(short_ms_positions) /
                     static_cast<double>(matched_positions);
  }

  double average_file_short_ms_percent() const {
    return files == 0 ? 0.0
                      : 100.0 * short_ms_share_sum / static_cast<double>(files);
  }

  std::uint32_t median_length() const {
    if (matched_positions == 0) return 0;
    const std::uint64_t half = (matched_positions + 1) / 2;
    std::uint64_t seen = 0;
    for (std::size_t length = 1; length < histogram.size(); ++length) {
      seen += histogram[length];
      if (seen >= half) return static_cast<std::uint32_t>(length);
    }
    return max_length;
  }

  void add(const PhraseStats& other) {
    files += other.files;
    positions += other.positions;
    matched_positions += other.matched_positions;
    unmatched += other.unmatched;
    match_length_sum += other.match_length_sum;
    short_ms_positions += other.short_ms_positions;
    max_length = std::max(max_length, other.max_length);
    if (histogram.size() < other.histogram.size()) {
      histogram.resize(other.histogram.size(), 0);
    }
    for (std::size_t length = 0; length < other.histogram.size(); ++length) {
      histogram[length] += other.histogram[length];
    }
    short_ms_share_sum += other.short_ms_share_sum;
  }

  // One line, e.g. "average match length 99.07, MS < 16 at 11.74% of
  // positions, max 1205, median 88".
  std::string format() const {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << "average match length "
        << average_match_length() << ", MS < 16 at " << short_ms_percent()
        << "% of positions";
    if (files > 1) {
      out << " (average per file " << average_file_short_ms_percent() << "%)";
    }
    out << ", max " << max_length << ", median " << median_length();
    if (unmatched != 0) {
      out << " (" << unmatched << " separator positions excluded)";
    }
    return out.str();
  }
};

inline PhraseStats phrase_stats(const MatchingStatistics& ms) {
  PhraseStats stats;
  stats.files = 1;
  stats.positions = ms.size();

  for (const auto& [position, length] : ms) {
    if (length == 0) {
      ++stats.unmatched;
      continue;
    }
    ++stats.matched_positions;
    stats.match_length_sum += length;
    if (length < KMER_LENGTH) ++stats.short_ms_positions;
    stats.max_length = std::max(stats.max_length, length);
    if (stats.histogram.size() <= length) {
      stats.histogram.resize(std::max<std::size_t>(length + 1,
                                                   2 * stats.histogram.size()),
                             0);
    }
    ++stats.histogram[length];
  }

  stats.short_ms_share_sum =
      stats.matched_positions == 0
          ? 0.0
          : static_cast<double>(stats.short_ms_positions) /
                static_cast<double>(stats.matched_positions);
  return stats;
}

}  // namespace ms2
