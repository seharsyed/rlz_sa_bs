/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef BENCHMARK_LOGGER_H
#define BENCHMARK_LOGGER_H

#include <string>

void write_sort_benchmark_jsonl(const std::string& input_path, 
                                const std::string& config_name,
                                size_t text_size,
                                double sort_time,
                                size_t suffix_comps, 
                                size_t unit_comps,
                                double avg_unit_per_comp,
                                // RLZ-specific metric defaults
                                bool is_rlz = false,
                                size_t total_factors = 0,
                                size_t indicative = 0,
                                size_t not_indicative = 0,
                                size_t interval_hits = 0,
                                double interval_percentage = 0.0,
                                size_t backbone_hits = 0,
                                double backbone_percentage = 0.0,
                                double preprocess_time = 0.0,
                                double resync_time = 0.0,
                                size_t resync = 0,
                                size_t resync_indicative = 0,
                                size_t resync_not_indicative = 0);

#endif