/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef SORT_ALGO_TEXT_H
#define SORT_ALGO_TEXT_H

#include <vector>
#include <string>
#include <string_view>

class TEXT_SORT 
{
    public:
        std::string seq_content;
        std::vector<size_t> suffix_array;

        TEXT_SORT(const std::string seq_file);
        ~TEXT_SORT();

        bool comparator(std::string_view a, std::string_view b, size_t& char_count);
        void build_sa();
        void write_json(const std::string seq_file);
        void write_sa(const std::string seq_file);
    
    private:
        // Metrics to record throughout
        size_t metric_text_size = 0;
        size_t metric_char_hits = 0;
        size_t metric_suffix_comps = 0;
        double metric_sort_time = 0;
        double metric_avg_char_per_comp = 0;
};

#endif  // SORT_ALGO_TEXT_H