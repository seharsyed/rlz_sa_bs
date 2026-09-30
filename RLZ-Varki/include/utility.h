/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef UTILITY_H
#define UTILITY_H

#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/suffix_arrays.hpp> 
#include <sdsl/lcp.hpp>           
#include <sdsl/int_vector.hpp> 
#include <sdsl/rmq_support.hpp>

// SA_SAMPLE_RATE, ISA_SAMPLE_RATE
using rlz_fm_index_t = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector>, 1, 4096>; // Only want the SA, not the ISA

using sort_csa_index_t = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector>, 1, 1>; // Want the SA only for resynchronization, ISA always needed
using sort_lcp_index_t = sdsl::lcp_bitcompressed<>; 
using sort_rmq_index_t = sdsl::rmq_succinct_sada<>;


#endif // UTILITY_H