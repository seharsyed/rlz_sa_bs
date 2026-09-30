/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef FM_WRAPPER_H
#define FM_WRAPPER_H

#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <fstream>
#include <vector>
#include <tuple>
#include <map>
#include "utility.h"

class FM_Wrapper 
{
    public:
        FM_Wrapper();
        ~FM_Wrapper();
        std::tuple<size_t, size_t> backward_match(const rlz_fm_index_t& fm_index,
                                                const std::vector<size_t>& occs, 
                                                const std::tuple<size_t, size_t>& prev_backward_range,
                                                const char next_char);
        
        size_t get_suffix_array_value(const rlz_fm_index_t& fm_index, const size_t location);
};


#endif // FM_WRAPPER_H