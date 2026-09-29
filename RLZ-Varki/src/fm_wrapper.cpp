/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include "fm_wrapper.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/int_vector.hpp>
#include <tuple>
#include <vector>
#include <map>
#include <stdexcept>


/**
* @brief FM Wrapper Constructor
*
* Does nothing 
*
*/

FM_Wrapper::FM_Wrapper(){}

/**
* @brief FM Wrapper Destructor
*
* Does nothing 
*
*/

FM_Wrapper::~FM_Wrapper(){}

/**
* @brief Extend backward match with FM-index (LF Mapping)
*
* A wrapper around rlz_fm_index fm_index.
* Calculates the new backwards search range after trying to match the next char (backwards)
* Can continue from previous character match so do not have to keep redoing previously done backwards matches.
*
* @param[in] fm_index [rlz_fm_index] the fm_index queried
* @param[in] occs [std::vector<size_t>] essentially the compressed version of F column of the BWT.
* @param[in] prev_backward_range [sdsl::range] the previous backwards search range.
* @param[in] next_char [char] the next character to match 
* 
* @return the backwards search range of next_char 
*/

std::tuple<size_t, size_t> FM_Wrapper::backward_match(const rlz_fm_index_t& fm_index,
                                                    const std::vector<size_t>& occs,
                                                    const std::tuple<size_t, size_t>& prev_backward_range,
                                                    const char next_char)
{
    std::tuple<size_t, size_t> next_char_backward_range;
    std::size_t next_left = fm_index.bwt.rank(std::get<0>(prev_backward_range), next_char);
    std::size_t next_right = fm_index.bwt.rank(std::get<1>(prev_backward_range), next_char);

    // There is a special (smaller) character appended so have to add 1 to the offset for both
    unsigned char u_next_char = static_cast<unsigned char>(next_char);
    try{
        if (u_next_char < 255 && occs[u_next_char] == occs[u_next_char + 1]) {
            throw std::runtime_error("Missing character");
        }
        next_left += occs[u_next_char] + 1;
    } catch (const std::runtime_error& e) {
        std::cerr << "Character code: " << static_cast<int>(u_next_char) << " ('" << next_char << "') not found in reference text!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    try{ 
        if (u_next_char < 255 && occs[u_next_char] == occs[u_next_char + 1]) {
            throw std::runtime_error("Missing character");
        }
        next_right += occs[u_next_char] + 1;
    } catch (const std::runtime_error& e) {
        std::cerr << "Character code: " << static_cast<int>(u_next_char) << " ('" << next_char << "') not found in reference text!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    return std::make_tuple(next_left,next_right);
}

/**
* @brief Get corresponding location in reference via suffix array
*
* A wrapper around rlz_fm_index fm_index.
* When the backwards search range becomes empty it means there is no perfect match with the current pattern that is being processed.
* With the prior non-empty range we find one location of the previous perfect pattern match. This should always be next_left of backward_match
*
* @param[in] fm_index [rlz_fm_index] the fm_index queried
* @param[in] location [size_t] I think should always be next_left of backward_match.
* 
* @return the suffix array index
*/

size_t FM_Wrapper::get_suffix_array_value(const rlz_fm_index_t& fm_index, const size_t location)
{
    return fm_index[location];
}