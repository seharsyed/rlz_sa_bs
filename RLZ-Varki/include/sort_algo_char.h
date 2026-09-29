/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef SORT_ALGO_CHAR_H
#define SORT_ALGO_CHAR_H

#include <string>
#include "utility.h"
#include "benchmark_logger.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <fstream>
#include <system_error>
#include <filesystem>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>

template<typename int_t>
class RLZ_CHAR_SORT
{
    public:
        struct RLZ_Factor { int_t p; int_t l; };

        // Represents a suffix starting at a specific character in the parsed text
        // factor_idx: The 0-based index of the factor within the sequence of factors where the suffix begins
        // offset: How far within the factor where the suffix starts
        struct SuffixID { size_t factor_idx; int_t offset; };

        // Represents a resynchronized state for the sort
        // id: The underlying SuffixID coordinate
        // first_factor: The factor after resynchronization or the original factor if option not used
        // borrowed: Tracks characters consumed from subsequent factors during resync
        // sa_range: Caches the [start, end] interval
        struct SortableSuffix { SuffixID id; RLZ_Factor first_factor; int_t borrowed; std::pair<int_t, int_t> sa_range; };

        std::string ref_content;
        sort_csa_index_t csa_ref; 
        sort_lcp_index_t lcp_ref; 
        sort_rmq_index_t lcp_rmq; 
        std::vector<RLZ_Factor> rlz_factors;
        std::vector<SortableSuffix> sorted_suffixes; // Stores the suffix array 

        RLZ_CHAR_SORT(const std::string ref_file, const std::string parse_file);
        ~RLZ_CHAR_SORT();

        void sort_naive(bool apply_resync = false);
        void sort_lcp_interval(bool apply_resync = false);
        void sort_induced(bool apply_resync = false);
        void sort_factors_only(bool apply_resync = false);
        void write_json(const std::string out_file, const std::string configuration);
        void stream_sa_to_file(const std::string out_file);
        void log_memory_estimate(const std::string& sort_method, uintmax_t ref_bytes) const;
    
    private:
        int_t get_lce(int_t i, int_t j); 
        bool compare_suffixes(const SortableSuffix& a, const SortableSuffix& b, const std::vector<size_t>* factor_ranks = nullptr);
        bool is_indicative(const RLZ_Factor& f);
        std::pair<int_t, int_t> get_sa_range(const RLZ_Factor& f);
        RLZ_Factor apply_resynchronization(const RLZ_Factor& f_i, const RLZ_Factor& f_next);
        std::vector<SortableSuffix> generate_all_suffixes(bool apply_resync);
        std::vector<SortableSuffix> generate_factor_boundaries(bool apply_resync);
        // Metrics to record throughout
        size_t metric_text_size = 0;
        size_t metric_factor_size = 0;  
        size_t metric_boundary_hits = 0;
        size_t metric_suffix_comps = 0;
        double metric_avg_boundary_per_comp = 0;
        size_t metric_interval_hits = 0;
        double interval_percentage = 0;
        size_t metric_backbone_hits = 0;
        double backbone_percentage = 0;
        size_t metric_indicative = 0;
        size_t metric_not_indicative = 0;
        size_t metric_resync = 0;
        size_t metric_resync_indicative = 0;
        size_t metric_resync_not_indicative = 0;
        double metric_sort_time = 0;
        double metric_preprocess_time = 0;
        double metric_resync_time = 0;
};

/**
 * @brief Initializes the RLZ sequence sorting framework.
 * Loads the reference text and the binary sequence of RLZ factors. 
 * Automatically constructs the necessary succinct data structures 
 * (Compressed Suffix Array, LCP Array, and RMQ support) over the 
 * reference string to enable O(1) LCE queries.
 * @param [in] ref_file   [std::string] Path to the raw reference text file.
 * @param [in] parse_file [std::string] Path to the binary RLZ parse file.
 */

template <typename int_t>
RLZ_CHAR_SORT<int_t>::RLZ_CHAR_SORT(const std::string ref_file, const std::string parse_file)
{
    spdlog::debug("Reading in reference file");
    spdlog::stopwatch sw_convert;

    std::error_code ec;
    uintmax_t ref_size = std::filesystem::file_size(ref_file, ec);
    if (ec) {
        spdlog::error("Error getting file size for {}: {}", ref_file, ec.message());
        std::exit(EXIT_FAILURE);
    }

    std::ifstream ref(ref_file, std::ios::binary);
    if (!ref) {
        spdlog::error("Error opening {}", ref_file);
        std::exit(EXIT_FAILURE);
    }

    ref_content.resize(ref_size);

    if (!ref.read(&ref_content[0], ref_size)) {
        spdlog::error("Error reading data from {}", ref_file);
        std::exit(EXIT_FAILURE);
    }
    ref.close();

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished reading reference file in {:.3f} seconds", sw_convert_elapsed.count());
    
    spdlog::debug("Constructing ISA, LCP, and RMQ from reference");
    spdlog::stopwatch sw_build;

    sdsl::construct_im(csa_ref, ref_content, 1);
    sdsl::construct_im(lcp_ref, ref_content, 1);
    lcp_rmq = sdsl::rmq_succinct_sada<>(&lcp_ref);

    auto sw_build_elapsed = sw_build.elapsed();
    spdlog::debug("Finished building data structures in {:.3f} seconds", sw_build_elapsed.count());
    
    spdlog::debug("Reading in RLZ parse");
    spdlog::stopwatch sw_parse;

    uint64_t num_pairs; 

    std::ifstream parse(parse_file, std::ios::binary);
    if (!parse) {
        spdlog::error("Error opening {}", parse_file);
        std::exit(EXIT_FAILURE);
    }

    parse.read(reinterpret_cast<char*>(&num_pairs), sizeof(num_pairs));

    metric_factor_size = num_pairs;

    spdlog::debug("The RLZ parse contains {} factors", num_pairs);

    uintmax_t parse_size = std::filesystem::file_size(parse_file);
    assert(parse_size - sizeof(uint64_t) == num_pairs * sizeof(RLZ_Factor));

    rlz_factors.resize(num_pairs);
    parse.read(reinterpret_cast<char*>(rlz_factors.data()), num_pairs * sizeof(RLZ_Factor));

    for (auto factor : rlz_factors){ metric_text_size += factor.l; }

    spdlog::debug("The RLZ parse represents {} characters", metric_text_size);

    auto sw_parse_elapsed = sw_parse.elapsed();
    spdlog::debug("Finished reading RLZ parse in {:.3f} seconds", sw_parse_elapsed.count());
}

/**
 * @brief Default destructor for the RLZ_CHAR_SORT class.
 */

template<typename int_t>
RLZ_CHAR_SORT<int_t>::~RLZ_CHAR_SORT(){}

/**
 * @brief Computes the Longest Common Extension (LCE) between two reference positions.
 * Implements Equation 2.2 from the paper. Utilizes the Inverse Suffix Array (ISA) 
 * and a Range Minimum Query (RMQ) over the LCP array to determine the length of 
 * the longest common prefix between two suffixes of R in O(1) time (assuming no compression).
 * @param [in] i [int_t] First position in the reference text.
 * @param [in] j [int_t] Second position in the reference text.
 * @return [int_t] The LCE length.
 */

template<typename int_t>
int_t RLZ_CHAR_SORT<int_t>::get_lce(int_t i, int_t j) {
    if (i == j) return ref_content.size() - i;
    
    int_t ri = csa_ref.isa[i];
    int_t rj = csa_ref.isa[j];
    
    if (ri > rj) std::swap(ri, rj);
    
    size_t min_idx = lcp_rmq(ri + 1, rj);
    return lcp_ref[min_idx];
}

/**
 * @brief Lexicographically compares two factor-aligned suffixes of the text.
 * Implements Algorithm 3.1. Evaluates the lexicographical order of two suffixes 
 * by performing LCE queries on their corresponding reference segments. This allows 
 * the comparison to "jump" over matching segments of the reference text in O(1) time 
 * per factor transition, avoiding decompression and character-by-character scans.
 * @param [in] a [SortableSuffix] The first suffix object.
 * @param [in] b [SortableSuffix] The second suffix object.
 * @param [in] factor_ranks [std::vector<size_t>*] Optional: The ISA of the sorted complete factors 
 * @return [bool] True if Suffix A is strictly less than Suffix B.
 */

template<typename int_t>
bool RLZ_CHAR_SORT<int_t>::compare_suffixes(const SortableSuffix& a, const SortableSuffix& b, const std::vector<size_t>* factor_ranks) {
    int_t pa = a.first_factor.p;
    int_t pb = b.first_factor.p;
    int_t la = a.first_factor.l;
    int_t lb = b.first_factor.l;
    
    size_t idx_a = a.id.factor_idx + 1; // The next factor following a 
    size_t idx_b = b.id.factor_idx + 1; // The next factor following b

    int_t skip_a = a.borrowed;
    int_t skip_b = b.borrowed;

    // Helper lambda to fetch the next valid block of characters, safely skipping 'borrowed' segments
    auto fetch_next = [&](size_t& idx, int_t& p, int_t& l, int_t& skip) 
    {
        while (idx < this->rlz_factors.size()) 
        {
            int_t next_l = this->rlz_factors[idx].l;
            if (skip >= next_l) {
                skip -= next_l; // Entire factor was swallowed
                idx++;
            } 
            else {
                p = this->rlz_factors[idx].p + skip;
                l = next_l - skip;
                skip = 0;
                idx++;
                return true;
            }
        }
        l = 0;
        return false;
    };
    
    while (true) {

        spdlog::trace("Comparing Factor A: ({},{}) with Factor B: ({},{})", pa, la, pb, lb);

        int_t k = get_lce(pa, pb);
        metric_boundary_hits++;

        spdlog::trace("LCE value: {} ---- min len: {}", k, std::min(la,lb));

        if (k < std::min(la, lb)) {
            spdlog::trace("Mismatch: A: ['{}'] ---- B: ['{}']", ref_content[pa +k], ref_content[pb + k]);
            return ref_content[pa + k] < ref_content[pb + k];
        } 
        
        if (la < lb) {
            spdlog::trace("Factor A fully consumed --- Factor B partially consumed");
            // Adjust Factor B by la
            pb += la; 
            lb -= la;
            // Move Suffix A to next factor if possible
            if (!fetch_next(idx_a, pa, la, skip_a)){
                spdlog::trace("Suffix A has no more factors"); 
                return true; 
            }
        } else if (lb < la) {
            spdlog::trace("Factor A partially consumed --- Factor B fully consumed");
            // Adjust Factor A by lb
            pa += lb; 
            la -= lb;
            // Move Suffix B to next factor if possible
            if (!fetch_next(idx_b, pb, lb, skip_b)){
                spdlog::trace("Suffix B has no more factors");  
                return false;
            }
        } else {
            spdlog::trace("Factor A fully consumed --- Factor B fully consumed");

            // THE O(1) BACKBONE SHORTCUT
            // If we have access to the sorted ranks, and neither factor has a 
            // resynchronization debt (skip == 0), we are sitting exactly on the 
            // boundaries of two complete factors. We can terminate immediately!
            
            if (factor_ranks != nullptr && skip_a == 0 && skip_b == 0) {
                spdlog::trace("Can use complete factor shortcut to finish comparison");
                metric_backbone_hits++;

                // Handle end-of-text boundaries
                if (idx_a >= rlz_factors.size() && idx_b >= rlz_factors.size()) {
                    spdlog::error("Suffixes A and B have no more factors! Should not be possible.");
                    return false;
                }
                if (idx_a >= rlz_factors.size()) {
                    spdlog::trace("Suffix A has no more factors"); 
                    return true;
                } 
                if (idx_b >= rlz_factors.size()) {
                    spdlog::trace("Suffix B has no more factors"); 
                    return false;
                }
                
                // O(1) instant resolution using the sorted backbone
                return (*factor_ranks)[idx_a] < (*factor_ranks)[idx_b];
            }

            // Fallback for suffixes carrying a resync debt
            // Move Suffixes A and B to next factor if possible
            bool a_has_more = fetch_next(idx_a, pa, la, skip_a);
            bool b_has_more = fetch_next(idx_b, pb, lb, skip_b);

            if (!a_has_more && !b_has_more) {
                spdlog::error("Suffixes A and B have no more factors! Should not be possible.");
                return false;
            }
            if (!a_has_more){
                spdlog::trace("Suffix A has no more factors"); 
                return true;
            }
            if (!b_has_more){
                spdlog::trace("Suffix B has no more factors"); 
                return false;
            }
        }
    }
}

/**
 * @brief Determines if an RLZ factor uniquely identifies a match (Algorithm 3.2).
 * An indicative factor has a length strictly greater than the longest common prefix 
 * it shares with its immediate lexicographical neighbors in the reference. 
 * Therefore, it corresponds to a unique match (a single SA range) in R.
 * @param [in] f [RLZ_Factor] The RLZ factor to evaluate.
 * @return [bool] True if the factor is indicative, false otherwise.
 */

template<typename int_t>
bool RLZ_CHAR_SORT<int_t>::is_indicative(const RLZ_Factor& f) {
    size_t i = csa_ref.isa[f.p];
    
    int_t v_prev = (i > 0) ? lcp_ref[i] : 0;
    int_t v_next = (i < csa_ref.size() - 1) ? lcp_ref[i + 1] : 0;
    
    return f.l > std::max(v_prev, v_next);
}

/**
 * @brief Finds the Suffix Array (SA) range for a factor using O(log N) RMQ Binary Search.
 * Utilizes the succinct Range Minimum Query (RMQ) structure to binary search the exact boundaries 
 * where the minimum LCP drops below the factor length.
 * @param [in] f [RLZ_Factor] The RLZ factor.
 * @return [std::pair<int_t, int_t>] The [start, end] indices in the Suffix Array.
 */
template<typename int_t>
std::pair<int_t, int_t> RLZ_CHAR_SORT<int_t>::get_sa_range(const RLZ_Factor& f) {
    size_t i = csa_ref.isa[f.p];
    
    // Binary Search for the Start Pointer (sp)
    size_t low_sp = 0;
    size_t high_sp = i;
    
    while (low_sp < high_sp) {
        size_t mid = low_sp + (high_sp - low_sp) / 2;
        
        // O(1) query: Find the minimum LCP in the range [mid + 1, i]
        size_t min_idx = lcp_rmq(mid + 1, i);
        
        if (lcp_ref[min_idx] >= f.l) {
            // The whole range from mid to i has LCP >= f.l
            // So the drop-off boundary must be further left (or exactly at mid)
            high_sp = mid; 
        } else {
            // The minimum is less than f.l, so mid is too far left
            low_sp = mid + 1; 
        }
    }
    size_t sp = low_sp;

    // Binary Search for the End Pointer (ep)
    size_t low_ep = i;
    size_t high_ep = csa_ref.size() - 1;
    
    while (low_ep < high_ep) {
        // Ceiling division (+1) prevents infinite loops when low == high - 1
        size_t mid = low_ep + (high_ep - low_ep + 1) / 2; 
        
        // O(1) query: Find the minimum LCP in the range [i + 1, mid]
        size_t min_idx = lcp_rmq(i + 1, mid);
        
        if (lcp_ref[min_idx] >= f.l) {
            // The whole range from i to mid has LCP >= f.l
            // So the drop-off boundary must be further right (or exactly at mid)
            low_ep = mid; 
        } else {
            // The minimum is less than f.l, so mid is too far right
            high_ep = mid - 1; 
        }
    }
    size_t ep = high_ep;

    return {static_cast<int_t>(sp), static_cast<int_t>(ep)};
}

/**
 * @brief Restores right-maximality to incomplete factors (Algorithm 3.4).
 * Optimized to use O(log N) binary search over the Suffix Array. 
 * Because all suffixes in range_i share the prefix f_i, their 
 * continuations are perfectly sorted, allowing us to find the maximum overlap 
 * instantly without scanning millions of candidate locations.
 * @param [in] f_i    [RLZ_Factor] The incomplete factor.
 * @param [in] f_next [RLZ_Factor] The factor immediately following f_i.
 * @return [RLZ_Factor] The updated right-maximal factor (or original if no overlap).
 */
template<typename int_t>
typename RLZ_CHAR_SORT<int_t>::RLZ_Factor RLZ_CHAR_SORT<int_t>::apply_resynchronization(const RLZ_Factor& f_i, const RLZ_Factor& f_next) {

    std::pair<int_t, int_t> range_i = get_sa_range(f_i);
    
    // Safety check: if factor doesn't exist, return original
    if (range_i.first > range_i.second) return f_i;

    int_t best_p = f_i.p; 
    int_t max_k = 0;
    bool perfect_match = false;
    
    // Binary search for the lexicographical insertion point of f_next
    size_t low = range_i.first;
    size_t high = static_cast<size_t>(range_i.second) + 1;
    
    while (low < high) {
        size_t mid = low + (high - low) / 2;
        int_t candidate_p = csa_ref[mid];
        
        // If the candidate suffix is too close to the end of the text, treat as smaller
        if (static_cast<size_t>(candidate_p + f_i.l) >= ref_content.size()) {
            low = mid + 1;
            continue;
        }

        int_t k = get_lce(candidate_p + f_i.l, f_next.p);
        
        if (k >= f_next.l) {
            // We found a perfect overlap that consumes all of f_next. Stop searching.
            max_k = f_next.l;
            best_p = candidate_p;
            perfect_match = true;
            break;
        }
        
        // Look at the differing character to guide the binary search
        char c_candidate = (static_cast<size_t>(candidate_p + f_i.l + k) < ref_content.size()) ? ref_content[candidate_p + f_i.l + k] : '\0';
                           
        char c_target = (static_cast<size_t>(f_next.p + k) < ref_content.size()) ? ref_content[f_next.p + k] : '\0';
        
        if (c_candidate < c_target) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    
    // If a perfect match wasn't found, the longest match is guaranteed to be 
    // at either the insertion point (low) or the element right before it (low - 1).
    if (!perfect_match) {
        // Evaluate 'low'
        if (low >= static_cast<size_t>(range_i.first) && low <= static_cast<size_t>(range_i.second)) {
            int_t candidate_p = csa_ref[low];
            if (static_cast<size_t>(candidate_p + f_i.l) < ref_content.size()) {
                int_t current_k = get_lce(candidate_p + f_i.l, f_next.p);
                if (current_k > max_k) {
                    max_k = std::min(current_k, f_next.l);
                    best_p = candidate_p;
                }
            }
        }
        
        // Evaluate 'low - 1'
        if (low > static_cast<size_t>(range_i.first) && (low - 1) <= static_cast<size_t>(range_i.second)) {
            int_t candidate_p = csa_ref[low - 1];
            if (static_cast<size_t>(candidate_p + f_i.l) < ref_content.size()) {
                int_t current_k = get_lce(candidate_p + f_i.l, f_next.p);
                if (current_k > max_k) {
                    max_k = std::min(current_k, f_next.l);
                    best_p = candidate_p;
                }
            }
        }
    }
    
    if (max_k > 0) {
        spdlog::trace("Resynchronization occured");
        return {best_p, static_cast<int_t>(f_i.l + max_k)};
    }
    
    spdlog::trace("Resynchronization did not occur");
    return f_i;
}

/**
 * @brief Generates the full set of character-level suffixes for the original text.
 * Unrolls the 1D array of RLZ factors into individual Suffix objects, assigning 
 * an internal offset to every character position. If resynchronization is enabled, 
 * it dynamically applies Algorithm 3.4 to incomplete suffixes (offset > 0) to 
 * bypass artificial boundaries before the sort begins.
 * @param [in] apply_resync [bool] Toggle to enable or bypass boundary repair.
 * @return [std::vector<SortableSuffix>] The pre-processed array of all suffixes.
 */

template<typename int_t>
std::vector<typename RLZ_CHAR_SORT<int_t>::SortableSuffix> RLZ_CHAR_SORT<int_t>::generate_all_suffixes(bool apply_resync) {
    spdlog::info("Generating all character-level suffixes (Resync: {})", apply_resync);
    spdlog::stopwatch sw_preprocess;
    
    // Get number of suffixes
    size_t total_suffixes = 0;
    std::vector<size_t> factor_starts(rlz_factors.size());
    for (size_t i = 0; i < rlz_factors.size(); ++i) {
        factor_starts[i] = total_suffixes;
        total_suffixes += rlz_factors[i].l;
    }

    // Pre-allocate the exact amount of required memory once
    std::vector<SortableSuffix> all_suffixes(total_suffixes);
    
    // Thread-safe metric variables
    size_t loc_ind = 0, loc_not_ind = 0, loc_resync_ind = 0;
    size_t loc_resync_not_ind = 0, loc_resync = 0;
    double loc_resync_time = 0.0;

    // Figure out if this loop is actually parallelized
    #ifdef _OPENMP
    #pragma omp parallel
    {
        // #pragma omp single ensures only ONE thread prints this message, 
        #pragma omp single 
        spdlog::info("OpenMP is ACTIVE! Running {} parallel threads.", omp_get_num_threads());
    }
    #else
        spdlog::warn("OpenMP is NOT active. Running sequentially on 1 thread.");
    #endif

    // Parallel Loop 
    #pragma omp parallel for schedule(dynamic, 1024) \
        reduction(+:loc_ind, loc_not_ind, loc_resync_ind, loc_resync_not_ind, loc_resync, loc_resync_time)
    for (size_t i = 0; i < rlz_factors.size(); ++i) {
        size_t base_idx = factor_starts[i];
        
        for (size_t offset = 0; offset < rlz_factors[i].l; ++offset) {
            SortableSuffix suf;
            suf.id = {i, static_cast<int_t>(offset)};
            
            int_t orig_l = static_cast<int_t>(rlz_factors[i].l - offset);
            RLZ_Factor effective_first = {
                static_cast<int_t>(rlz_factors[i].p + offset), 
                orig_l
            };

            bool local_is_ind = is_indicative(effective_first);

            if (local_is_ind) { loc_ind++; } else { loc_not_ind++; }
            
            // Thread-safe resynchronization timing
            if (!local_is_ind && offset > 0 && apply_resync && i + 1 < rlz_factors.size()) {
                #ifdef _OPENMP
                double start_time = omp_get_wtime(); // High-res OpenMP timer
                #endif

                effective_first = apply_resynchronization(effective_first, rlz_factors[i+1]);

                #ifdef _OPENMP
                loc_resync_time += (omp_get_wtime() - start_time);
                #endif

                local_is_ind = is_indicative(effective_first); 
            } 

            if (local_is_ind) { loc_resync_ind++; } else { loc_resync_not_ind++; }
            
            suf.first_factor = effective_first;
            suf.borrowed = static_cast<int_t>(effective_first.l - orig_l);

            if (suf.borrowed > 0) { loc_resync++; } 

            // Calculate the sa_range
            if (local_is_ind) {
                int_t rank = csa_ref.isa[effective_first.p];
                suf.sa_range = {rank, rank}; 
            } else {
                suf.sa_range = get_sa_range(effective_first);
            }
            
            // Lock-free direct memory write
            all_suffixes[base_idx + offset] = suf;
        }
    }
    
    // Aggregate parallel results back into the global class metrics
    metric_indicative += loc_ind;
    metric_not_indicative += loc_not_ind;
    metric_resync_indicative += loc_resync_ind;
    metric_resync_not_indicative += loc_resync_not_ind;
    metric_resync += loc_resync;
    metric_resync_time += loc_resync_time;
    
    spdlog::debug("Prior to resynchronization there were {} indicative factors", metric_indicative);
    spdlog::debug("Prior to resynchronization there were {} non-indicative factors", metric_not_indicative);
    spdlog::debug("After resynchronization there were {} indicative factors", metric_resync_indicative);
    spdlog::debug("After resynchronization there were {} non-indicative factors", metric_resync_not_indicative);
    spdlog::debug("{} factors were resynchronized", metric_resync);
    spdlog::debug("Resynchronization of factors took {:.3f} seconds", metric_resync_time);

    metric_preprocess_time = sw_preprocess.elapsed().count();
    spdlog::info("Finished generating all character-level suffixes in {:.3f} seconds", metric_preprocess_time);

    return all_suffixes;
}

/**
 * @brief Generates complete factor boundaries (offset == 0) for partial SA construction.
 * Supports optional resynchronization, which is highly valuable for non-greedy 
 * parsing strategies where boundaries are structurally fixed but not strictly right-maximal.
 * @param [in] apply_resync [bool] Toggle to enable boundary right-extension.
 * @return [std::vector<SortableSuffix>] Array of complete factors with cached intervals.
 */
template<typename int_t>
std::vector<typename RLZ_CHAR_SORT<int_t>::SortableSuffix> RLZ_CHAR_SORT<int_t>::generate_factor_boundaries(bool apply_resync) {
    spdlog::info("Generating ONLY factor boundaries (Resync enabled: {})", apply_resync);
    spdlog::stopwatch sw_preprocess;

    std::vector<SortableSuffix> boundaries(rlz_factors.size());
    
    // Thread-safe metric variables
    size_t loc_ind = 0, loc_not_ind = 0, loc_resync_ind = 0;
    size_t loc_resync_not_ind = 0, loc_resync = 0;
    double loc_resync_time = 0.0;

    // Figure out if this loop is actually parallelized
    #ifdef _OPENMP
    #pragma omp parallel
    {
        // #pragma omp single ensures only ONE thread prints this message, 
        #pragma omp single 
        spdlog::info("OpenMP is ACTIVE! Running {} parallel threads.", omp_get_num_threads());
    }
    #else
        spdlog::warn("OpenMP is NOT active. Running sequentially on 1 thread.");
    #endif

    // Parallel Loop 
    #pragma omp parallel for schedule(dynamic, 1024) \
        reduction(+:loc_ind, loc_not_ind, loc_resync_ind, loc_resync_not_ind, loc_resync, loc_resync_time)
    for (size_t i = 0; i < rlz_factors.size(); ++i) {
        SortableSuffix suf;
        suf.id = {i, 0}; 
        
        RLZ_Factor effective_first = rlz_factors[i];
        int_t orig_l = effective_first.l;
        
        bool local_is_ind = is_indicative(effective_first);

        if (local_is_ind) { loc_ind++; } else { loc_not_ind++; }
        
        // Thread-safe resynchronization timing
        if (!local_is_ind && apply_resync && i + 1 < rlz_factors.size()) {
            #ifdef _OPENMP
            double start_time = omp_get_wtime();
            #endif

            effective_first = apply_resynchronization(effective_first, rlz_factors[i+1]);

            #ifdef _OPENMP
            loc_resync_time += (omp_get_wtime() - start_time);
            #endif
            
            local_is_ind = is_indicative(effective_first); 
        } 

        if (local_is_ind) { loc_resync_ind++; } else { loc_resync_not_ind++; }
        
        suf.first_factor = effective_first;
        suf.borrowed = static_cast<int_t>(effective_first.l - orig_l);

        if (suf.borrowed > 0) { loc_resync++; } 
        
        // Calculate the sa_range
        if (local_is_ind) {
            int_t rank = csa_ref.isa[effective_first.p];
            suf.sa_range = {rank, rank}; 
        } else {
            suf.sa_range = get_sa_range(effective_first);
        }
        
        // Lock-free direct memory write
        boundaries[i] = suf;
    }

    // Aggregate parallel results back into the global class metrics
    metric_indicative += loc_ind;
    metric_not_indicative += loc_not_ind;
    metric_resync_indicative += loc_resync_ind;
    metric_resync_not_indicative += loc_resync_not_ind;
    metric_resync += loc_resync;
    metric_resync_time += loc_resync_time;

    spdlog::debug("Prior to resynchronization there were {} indicative factors", metric_indicative);
    spdlog::debug("Prior to resynchronization there were {} non-indicative factors", metric_not_indicative);
    spdlog::debug("After resynchronization there were {} indicative factors", metric_resync_indicative);
    spdlog::debug("After resynchronization there were {} non-indicative factors", metric_resync_not_indicative);
    spdlog::debug("{} factors were resynchronized", metric_resync);
    spdlog::debug("Resynchronization of factors took {:.3f} seconds", metric_resync_time);
    
    metric_preprocess_time = sw_preprocess.elapsed().count();
    spdlog::info("Finished generating all factor-level suffixes in {:.3f} seconds", metric_preprocess_time);

    return boundaries;
}

/**
 * @brief Baseline sorting strategy (O(nlogn) factor comparisons).
 * Evaluates all suffixes indiscriminately using introsort combined with Algorithm 3.1. 
 * Does not utilize LCP interval pruning or induced sorting mechanics. Acts as the 
 * performance control for the benchmark.
 * @param [in] apply_resync [bool] Toggle to enable resynchronization preprocessing.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_naive(bool apply_resync) {

    // Preprocessing 
    std::vector<SortableSuffix> all_suffixes = generate_all_suffixes(apply_resync);

    // Actual Sorting Logic
    spdlog::info("Executing Naive Sort (Resync: {})", apply_resync);
    spdlog::stopwatch sw_sort;

    std::sort(all_suffixes.begin(), all_suffixes.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        metric_suffix_comps++;
        return compare_suffixes(a, b);
    });

    sorted_suffixes = std::move(all_suffixes); // O(1) pointer swap

    all_suffixes.clear();
    all_suffixes.shrink_to_fit();

    metric_sort_time = sw_sort.elapsed().count();
    spdlog::debug("Number of suffix comparison performed: {}", metric_suffix_comps);
    interval_percentage = (static_cast<double>(metric_interval_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with LCP-intervals: {} ({:.2f}%)", metric_interval_hits, interval_percentage);
    backbone_percentage = (static_cast<double>(metric_backbone_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with complete factor backbone: {} ({:.2f}%)", metric_backbone_hits, backbone_percentage);
    spdlog::debug("Number of LCE queries performed during sorting: {}", metric_boundary_hits);
    metric_avg_boundary_per_comp = static_cast<double>(metric_boundary_hits) / static_cast<double>(metric_suffix_comps);
    spdlog::debug("Average number of LCE queries per suffix comparison: {:.3f}", metric_avg_boundary_per_comp);
    spdlog::info("Naive Sort completed in {:.3f} seconds", metric_sort_time);
}

/**
 * @brief Sorting strategy utilizing LCP Interval Pruning.
 * Sorts suffixes by first checking their precomputed indicativity flags. If two 
 * indicative suffixes are compared, their lexicographical order is resolved 
 * immediately in O(1) time via their ISA ranks, bypassing LCE jumps entirely.
 * @param [in] apply_resync [bool] Toggle to enable resynchronization preprocessing.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_lcp_interval(bool apply_resync) {

    // Preprocessing
    std::vector<SortableSuffix> all_suffixes = generate_all_suffixes(apply_resync);

    // Actual Sorting Logic
    spdlog::info("Executing LCP Interval Sort (Resync: {})", apply_resync);
    spdlog::stopwatch sw_sort;

    std::sort(all_suffixes.begin(), all_suffixes.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        metric_suffix_comps++;
        // The Disjoint Interval Shortcut (O(1) resolution)
        if (a.sa_range.second < b.sa_range.first) {
            metric_interval_hits++;
            return true;
        }
        if (b.sa_range.second < a.sa_range.first) {
            metric_interval_hits++;
            return false;
        }

        // The Overlap Fallback (LCE Evaluation)
        return compare_suffixes(a, b);
    });

    sorted_suffixes = std::move(all_suffixes); // O(1) pointer swap

    all_suffixes.clear();
    all_suffixes.shrink_to_fit();

    metric_sort_time = sw_sort.elapsed().count();
    spdlog::debug("Number of suffix comparison performed: {}", metric_suffix_comps);
    interval_percentage = (static_cast<double>(metric_interval_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with LCP-intervals: {} ({:.2f}%)", metric_interval_hits, interval_percentage);
    backbone_percentage = (static_cast<double>(metric_backbone_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with complete factor backbone: {} ({:.2f}%)", metric_backbone_hits, backbone_percentage);
    spdlog::debug("Number of LCE queries performed during sorting: {}", metric_boundary_hits);
    metric_avg_boundary_per_comp = static_cast<double>(metric_boundary_hits) / static_cast<double>(metric_suffix_comps);
    spdlog::debug("Average number of LCE queries per suffix comparison: {:.3f}", metric_avg_boundary_per_comp);
    spdlog::info("LCP Interval Sort completed in {:.3f} seconds", metric_sort_time);
}

/**
 * @brief Highly optimized sorting strategy leveraging Lemma 3.2.
 * Bins suffixes into two categories: Complete (offset == 0) and Incomplete (offset > 0).
 * Sorts the complete, highly indicative backbone first. It then leverages this stable 
 * structural backbone to dramatically reduce the search space and LCE jump overhead 
 * when resolving the more difficult incomplete suffixes.
 * @param [in] apply_resync [bool] Toggle to enable resynchronization preprocessing.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_induced(bool apply_resync) {

    // Preprocessing
    std::vector<SortableSuffix> all_suffixes = generate_all_suffixes(apply_resync);
    
    std::vector<SortableSuffix> complete_bin;
    std::vector<SortableSuffix> incomplete_bin;
    
    for (auto& suf : all_suffixes) {
        if (suf.id.offset == 0) complete_bin.push_back(std::move(suf));
        else incomplete_bin.push_back(std::move(suf)); // To avoid duplicating the data
    }

    all_suffixes.clear();
    all_suffixes.shrink_to_fit();

    // Actual Sorting Logic
    spdlog::info("Executing Induced Sort (Resync: {})", apply_resync);
    spdlog::stopwatch sw_sort;

    // Sort the complete factor suffixes first since they are the most likely to be indicative
    std::sort(complete_bin.begin(), complete_bin.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        metric_suffix_comps++;        
        // The Disjoint Interval Shortcut (O(1) resolution)
        if (a.sa_range.second < b.sa_range.first) {
            metric_interval_hits++;
            return true;
        }
        if (b.sa_range.second < a.sa_range.first) {
            metric_interval_hits++;
            return false;
        }

        // The Overlap Fallback (LCE Evaluation)
        return compare_suffixes(a, b);
    });

    // Build the ISA of the sorted complete factors
    std::vector<size_t> factor_ranks(rlz_factors.size(), 0);
    for (size_t r = 0; r < complete_bin.size(); ++r) {
        factor_ranks[complete_bin[r].id.factor_idx] = r;
    }

    // Sort the incomplete factors suffixes using the sorted complete factor suffixes knowledge
    // Note this speedup only occurs if two factors are fully consumed at the same time
    std::sort(incomplete_bin.begin(), incomplete_bin.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        metric_suffix_comps++;
        // The Disjoint Interval Shortcut (O(1) resolution)
        if (a.sa_range.second < b.sa_range.first) {
            metric_interval_hits++;
            return true;
        }
        if (b.sa_range.second < a.sa_range.first) {
            metric_interval_hits++;
            return false;
        }

        // The Overlap Fallback (LCE Evaluation)
        return compare_suffixes(a, b, &factor_ranks);
    });

    sorted_suffixes.reserve(complete_bin.size() + incomplete_bin.size());
    
    // Merge the sorted complete and incomplete factors
    size_t i = 0, j = 0;
    while (i < complete_bin.size() && j < incomplete_bin.size()) {
        SortableSuffix& a = complete_bin[i];
        SortableSuffix& b = incomplete_bin[j];
        
        bool a_is_smaller;

        // The Disjoint Interval Shortcut (O(1) resolution)
        if (a.sa_range.second < b.sa_range.first) {
            a_is_smaller = true;
        } 
        else if (b.sa_range.second < a.sa_range.first) {
            a_is_smaller = false;
        } 
        // The Overlap Fallback (LCE Evaluation)
        else {
            a_is_smaller = compare_suffixes(a, b, &factor_ranks);
        }

        // Push the winner to the final array
        if (a_is_smaller) {
            sorted_suffixes.push_back(std::move(complete_bin[i++]));
        } else {
            sorted_suffixes.push_back(std::move(incomplete_bin[j++]));
        }
    }
    // Flush the remaining elements
    if (i < complete_bin.size()) {
        this->sorted_suffixes.insert(
            this->sorted_suffixes.end(),
            std::make_move_iterator(complete_bin.begin() + i),
            std::make_move_iterator(complete_bin.end())
        );
    }
    if (j < incomplete_bin.size()) {
        this->sorted_suffixes.insert(
            this->sorted_suffixes.end(),
            std::make_move_iterator(incomplete_bin.begin() + j),
            std::make_move_iterator(incomplete_bin.end())
        );
    }

    complete_bin.clear();
    complete_bin.shrink_to_fit();
    incomplete_bin.clear();
    incomplete_bin.shrink_to_fit();

    metric_sort_time = sw_sort.elapsed().count();
    spdlog::debug("Number of suffix comparison performed: {}", metric_suffix_comps);
    interval_percentage = (static_cast<double>(metric_interval_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with LCP-intervals: {} ({:.2f}%)", metric_interval_hits, interval_percentage);
    backbone_percentage = (static_cast<double>(metric_backbone_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with complete factor backbone: {} ({:.2f}%)", metric_backbone_hits, backbone_percentage);
    spdlog::debug("Number of LCE queries performed during sorting: {}", metric_boundary_hits);
    metric_avg_boundary_per_comp = static_cast<double>(metric_boundary_hits) / static_cast<double>(metric_suffix_comps);
    spdlog::debug("Average number of LCE queries per suffix comparison: {:.3f}", metric_avg_boundary_per_comp);
    spdlog::info("Induced Sort completed in {:.3f} seconds", metric_sort_time);
}

/**
 * @brief Sorts only the complete RLZ factors to output a partial Suffix Array.
 * @param [in] apply_resync [bool] Toggle to enable or bypass boundary repair.
 * @note apply_resync would only be needed if user specified match-length during initial RLZ parsing.
 * @return void
 */
template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_factors_only(bool apply_resync) {
    
    // Preprocessing
    std::vector<SortableSuffix> complete_bin = generate_factor_boundaries(apply_resync);

    // Actual Sorting Logic
    spdlog::info("Executing Factor-Only Sort (Resync flag: {})", apply_resync);
    spdlog::stopwatch sw_sort;

    std::sort(complete_bin.begin(), complete_bin.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        metric_suffix_comps++;
        // The Disjoint Interval Shortcut (O(1) resolution)
        if (a.sa_range.second < b.sa_range.first) {
            metric_interval_hits++;
            return true;
        }
        if (b.sa_range.second < a.sa_range.first) {
            metric_interval_hits++;
            return false;
        }
        
        // The Overlap Fallback (LCE Evaluation)
        return compare_suffixes(a, b);
    });

    sorted_suffixes = std::move(complete_bin); // O(1) pointer swap

    complete_bin.clear();
    complete_bin.shrink_to_fit();

    metric_sort_time = sw_sort.elapsed().count();
    spdlog::debug("Number of suffix comparison performed: {}", metric_suffix_comps);
    interval_percentage = (static_cast<double>(metric_interval_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with LCP-intervals: {} ({:.2f}%)", metric_interval_hits, interval_percentage);
    backbone_percentage = (static_cast<double>(metric_backbone_hits) / static_cast<double>(metric_suffix_comps)) * 100.0;
    spdlog::debug("Number of comparisons resolved with complete factor backbone: {} ({:.2f}%)", metric_backbone_hits, backbone_percentage);
    spdlog::debug("Number of LCE queries performed during sorting: {}", metric_boundary_hits);
    metric_avg_boundary_per_comp = static_cast<double>(metric_boundary_hits) / static_cast<double>(metric_suffix_comps);
    spdlog::debug("Average number of LCE queries per suffix comparison: {:.3f}", metric_avg_boundary_per_comp);
    spdlog::info("Factor-Only Sort completed in {:.3f} seconds", metric_sort_time);
}

/**
 * @brief Dumps internal sorting metrics to the JSON Lines benchmark log.
 * Acts as a class-level adapter for the global benchmark logger. It automatically 
 * aggregates all tracked private metrics (O(1) shortcut hits, parsing stats, 
 * and execution timings) and pushes them to disk using the RLZ-specific JSON schema.
 * Must be called at the very end of a sorting routine after final percentages are calculated.
 *
 * @param parse_file    Base path/identifier for the input data, used to name the output log file.
 * @param configuration Identifier for the specific sorting strategy executed (e.g., "Induced_Sort").
 */


template<typename int_t>
void RLZ_CHAR_SORT<int_t>::write_json(const std::string parse_file, const std::string configuration)
{
    write_sort_benchmark_jsonl(
        parse_file,                      // input_path 
        configuration,                   // config_name 
        metric_text_size,                // text_size
        metric_sort_time,                // sort_time
        metric_suffix_comps,             // suffix_comps
        metric_boundary_hits,            // unit_comps (LCE boundary checks)
        metric_avg_boundary_per_comp,    // avg_unit_per_comp
        true,                            // is_rlz
        metric_factor_size,              // total_factors
        metric_indicative,               // indicative
        metric_not_indicative,           // not_indicative
        metric_interval_hits,            // interval_hits
        interval_percentage,             // interval_percentage
        metric_backbone_hits,            // backbone_hits
        backbone_percentage,             // backbone_percentage
        metric_preprocess_time,          // preprocess_time
        metric_resync_time,              // resync_time
        metric_resync,                   // resync
        metric_resync_indicative,        // resync_indicative
        metric_resync_not_indicative     // resync_not_indicative
    );
}


/**
 * @brief Streams the 2D Suffix Array to disk as a 1D Suffix Array in text format.
 * Bypasses full 1D vector allocation to save memory (preventing OOM on large datasets).
 * Calculates the absolute 1D character index on the fly using a prefix sum of factor lengths.
 * Writes each absolute index as a standard text string followed by a newline.
 * @param [in] parse_file [std::string] The input parse file path used to generate output path.
 */
template<typename int_t>
void RLZ_CHAR_SORT<int_t>::stream_sa_to_file(const std::string parse_file) {

    std::string out_file = parse_file + ".sa";

    spdlog::info("Streaming 1D suffix array directly to text file {}", out_file);
    spdlog::stopwatch sw_stream;

    std::ofstream out(out_file);
    if (!out) {
        spdlog::error("Error opening output file {}", out_file);
        std::exit(EXIT_FAILURE);
    }

    // Precompute absolute starting positions of each factor in T
    std::vector<size_t> factor_starts;
    factor_starts.reserve(rlz_factors.size());
    
    size_t current_text_pos = 0;
    for (const auto& factor : rlz_factors) {
        factor_starts.push_back(current_text_pos);
        current_text_pos += factor.l;
    }

    // Stream elements directly to the file, one per line.
    for (const auto& suf : this->sorted_suffixes) {
        size_t absolute_index = factor_starts[suf.id.factor_idx] + suf.id.offset;
        out << absolute_index << '\n';
    }

    out.close();
    spdlog::info("Finished streaming suffix array to file in {:.3f} seconds", sw_stream.elapsed().count());
}

/**
 * @brief Calculates and logs a highly accurate peak memory estimate and safe SLURM allocation request.
 * This function dynamically estimates the maximum RAM footprint required during the sorting pipeline. 
 * 
 * It accounts for two primary memory consumers:
 * 1. The target sequence: Scaled by the exact byte-size of the SortableSuffix struct. If performing a 
 * factor-only sort, this scales with number of factors. Otherwise, it scales with total characters.
 * 2. The reference sequence: Assumes a ~3.5x multiplier on the raw file size to account for the construction 
 * of the SDSL Compressed Suffix Array, bit-compressed LCP, ISA samples, and RMQ succinct structures.
 * 
 * Finally, a 15% safety buffer is applied to account for Linux OS file buffers and OpenMP thread overhead 
 * to generate a guaranteed-safe #SBATCH --mem directive.
 * 
 * @param [in] sort_method [std::string] The exact sorting algorithm being executed (e.g., "Induced_Sort").
 * @param [in] ref_bytes [uintmax_t] The uncompressed size of the reference sequence file in bytes.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::log_memory_estimate(const std::string& sort_method, uintmax_t ref_bytes) const {
    
    // Calculate total target characters and total factors
    size_t total_chars = 0;
    for (const auto& factor : rlz_factors) {
        total_chars += factor.l;
    }
    size_t total_factors = rlz_factors.size();
    
    // Struct size
    double bytes_per_element = sizeof(SortableSuffix);
    
    // Calculate target sequence RAM footprint
    double target_memory_bytes = 0.0;
    if (sort_method == "Factor_Only_Sort") {
        target_memory_bytes = total_factors * bytes_per_element;
    } else {
        target_memory_bytes = total_chars * bytes_per_element;
    }
    
    // Calculate SDSL Reference RAM footprint
    // CSA, LCP, ISA, and RMQ safely peak at roughly 3.5 bytes per reference character
    double sdsl_memory_bytes = static_cast<double>(ref_bytes) * 3.5; 

    // Combine and convert to Gigabytes
    double peak_memory_bytes = target_memory_bytes + sdsl_memory_bytes;
    double mem_gb = peak_memory_bytes / (1024.0 * 1024.0 * 1024.0);
    
    // Calculate a safe SLURM request (Peak + 15% for the Linux OS and OpenMP buffers)
    double slurm_safe_gb = std::ceil(mem_gb * 1.15);
    
    // Print the precise dashboard
    spdlog::info("================ MEMORY ESTIMATE ================");
    
    spdlog::info("Sortable Suffix is {} bytes", bytes_per_element);

    // Format Reference Size
    if (ref_bytes < 1024 * 1024 * 1024) {
        spdlog::info("Reference Size  : {:.2f} MB (SDSL RAM: ~{:.2f} MB)", 
                     static_cast<double>(ref_bytes) / (1024.0 * 1024.0), 
                     sdsl_memory_bytes / (1024.0 * 1024.0));
    } else {
        spdlog::info("Reference Size  : {:.2f} GB (SDSL RAM: ~{:.2f} GB)", 
                     static_cast<double>(ref_bytes) / (1024.0 * 1024.0 * 1024.0), 
                     sdsl_memory_bytes / (1024.0 * 1024.0 * 1024.0));
    }

    // Format Target Sequence
    if (total_chars < 1e9) {
        spdlog::info("Target Sequence : ~{:.2f} Million Characters", static_cast<double>(total_chars) / 1e6);
    } else {
        spdlog::info("Target Sequence : ~{:.2f} Billion Characters", static_cast<double>(total_chars) / 1e9);
    }

    spdlog::info("Sorting Method  : {}", sort_method);

    // Format Peak RAM
    if (mem_gb < 1.0) {
        spdlog::info("Peak RAM Needed : {:.2f} MB", peak_memory_bytes / (1024.0 * 1024.0));
    } else {
        spdlog::info("Peak RAM Needed : {:.2f} GB", mem_gb);
    }
    
    if (slurm_safe_gb > 64.0) {
        spdlog::warn("HIGH MEMORY WARNING: This job requires a large node.");
    }
    
    // Ensure we always request at least 1G from SLURM
    double final_slurm_request = std::max(1.0, slurm_safe_gb);
    spdlog::info("SLURM Request   : #SBATCH --mem={:.0f}G", final_slurm_request);
    spdlog::info("=================================================");
}

#endif  // SORT_ALGO_CHAR_H