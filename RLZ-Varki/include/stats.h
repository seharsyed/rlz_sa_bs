/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef STATS_H
#define STATS_H

#include <CLI11.hpp>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "spdlog/spdlog.h"

template<typename int_t>
void calculate_stats(const std::string& parse_file, const uintmax_t ref_size)
{

    // Reading in parse file
    std::ifstream parse(parse_file, std::ios::binary);
    if (!parse) {
        spdlog::error("Error opening {}", parse_file);
        std::exit(EXIT_FAILURE);
    }

    uint64_t num_pairs = 0;
    parse.read(reinterpret_cast<char*>(&num_pairs), sizeof(num_pairs));

    int_t pos; 
    int_t len;

    // Welford's Algorithm state variables
    uint64_t count = 0;
    double mean = 0.0;
    double M2 = 0.0;

    // Sequence stats
    uint64_t total_chars = 0;
    
    // Interval tracking for unique positions covered
    // We store [start, end) ranges.
    std::vector<std::pair<int_t, int_t>> intervals;
    intervals.reserve(num_pairs);

    for (uint64_t i = 0; i < num_pairs; ++i) {
        parse.read(reinterpret_cast<char*>(&pos), sizeof(pos));
        parse.read(reinterpret_cast<char*>(&len), sizeof(len));

        if (!parse) {
            spdlog::error("Early EOF or read error at pair index {}", i);
            break;
        }

        // Keep track of the total amount of chars in the sequence file
        total_chars += len;

        // Welford's algorithm for length mean and variance
        count++;
        double delta = len - mean;
        mean += delta / count;
        double delta2 = len - mean;
        M2 += delta * delta2;

        // Store interval for unique position coverage
        // pos is 0-indexed, so the phrase covers [pos, pos + len)
        intervals.push_back({pos, pos + len});
    }

    // Calculate Final Standard Deviation (Sample StdDev)
    double variance = (count > 1) ? M2 / (count - 1) : 0.0;
    double std_dev = std::sqrt(variance);

    // Calculate Unique Positions Covered via Interval Merging
    uint64_t unique_covered = 0;
    if (!intervals.empty()) {
        // Sort intervals by start position
        std::sort(intervals.begin(), intervals.end());
        
        int_t current_start = intervals[0].first;
        int_t current_end = intervals[0].second;

        for (size_t i = 1; i < intervals.size(); ++i) {
            if (intervals[i].first <= current_end) {
                // Intervals overlap or are contiguous; extend the current end
                current_end = std::max(current_end, intervals[i].second);
            } else {
                // Disjoint interval found; commit the length of the previous merged interval
                unique_covered += (current_end - current_start);
                current_start = intervals[i].first;
                current_end = intervals[i].second;
            }
        }
        // Commit the final interval
        unique_covered += (current_end - current_start);
    }

    // Output the calculated statistics
    spdlog::info("--- RLZ Parse Statistics ---");
    spdlog::info("Total Phrases (pairs): {}", count);
    spdlog::info("Original Sequence Chars: {}", total_chars);
    spdlog::info("Phrase Length Mean: {:.2f}", mean);
    spdlog::info("Phrase Length Std. Dev: {:.2f}", std_dev);
    double coverage_percent = (static_cast<double>(unique_covered) / ref_size) * 100.0;
    spdlog::info("Percentage of Ref Positions Covered: {:.2f}%", coverage_percent);
}

#endif // STATS_H