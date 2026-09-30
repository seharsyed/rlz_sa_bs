/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include "sort_algo_text.h"
#include "benchmark_logger.h"
#include <iostream>
#include <vector>
#include <string>
#include <string_view> //GCC 7
#include <algorithm>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <fstream>
#include <system_error>
#include <filesystem>


/**
 * @brief Constructor of the TEXT_SORT
 * @param [in] seq_file [string] Path to the sequence file from which to construct the Suffix Array
 */

TEXT_SORT::TEXT_SORT(const std::string seq_file)
{
    spdlog::debug("Reading in sequence file");

    spdlog::stopwatch sw_convert;

    // Getting size of sequence file
    std::error_code ec;
    uintmax_t seq_size = std::filesystem::file_size(seq_file, ec);
    if (ec) {
        spdlog::error("Error getting file size for {}: {}", seq_file, ec.message());
        std::exit(EXIT_FAILURE);
    }

    // Opening sequence file
    std::ifstream seq(seq_file, std::ios::binary);
    if (!seq) {
        spdlog::error("Error opening {}", seq_file);
        std::exit(EXIT_FAILURE);
    }

    // Preloading size of sequence buffer
    seq_content.resize(seq_size);

    // Directly loading sequence into buffer
    if (!seq.read(&seq_content[0], seq_size)) {
        spdlog::error("Error reading data from {}", seq_file);
        std::exit(EXIT_FAILURE);
    }
    seq.close();

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished reading sequence file in {:.3f} seconds", sw_convert_elapsed.count());
}

/**
* @brief Destructor of TEXT_SORT class.
*
* Currently does nothing.
*
*/
TEXT_SORT::~TEXT_SORT(){}


/**
 * @brief Comparison function to use for sorting
 * 
 * @param [in] a [string_view] suffix a
 * @param [in] b [string_view] suffix b
 * @param [in] char_count [size_t] the number of character comparsions at this timepoint
 * 
 * @return a less than b
 */

bool TEXT_SORT::comparator(std::string_view a, std::string_view b, size_t& char_count)
{ 
    size_t min_len = std::min(a.length(), b.length());

    for (size_t i = 0; i < min_len; ++i) {
        char_count++; 
        if (a[i] != b[i]) {
            return a[i] < b[i];
        }
    }

    char_count++; 
    return a.length() < b.length();
} 


/**
 * @brief Creates the suffix array of a text with naive method using custom comparison operator
 * 
 * Sorts the suffixes using a custom comparison operator. Worst case is O(N^2log(N)) since 
 * there are O(Nlog(N)) comparisons and each comparison takes O(N) time worst case. There 
 * are more efficient manners to create the suffix array, but the purpose is to compare 
 * sorting time with RLZ version, therefore an inefficient but compatible method with RLZ was chosen.
 * 
 */

void TEXT_SORT::build_sa() 
{
    spdlog::stopwatch sw_sort;

    metric_text_size = seq_content.size();
    suffix_array.resize(metric_text_size);
    
    // Fill the SA with the indices to start
    for (size_t i = 0; i < metric_text_size; i++) {
        suffix_array[i] = i;
    }

    // Read only view into original text
    std::string_view view(seq_content);

    std::sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b) {
        metric_suffix_comps++; 
        return comparator(view.substr(a), view.substr(b), metric_char_hits); 
    });


    metric_sort_time = sw_sort.elapsed().count();
    metric_avg_char_per_comp = static_cast<double>(metric_char_hits) / static_cast<double>(metric_suffix_comps);

    spdlog::debug("Number of suffix comparison performed: {}", metric_suffix_comps);
    spdlog::debug("Number of character comparisons performed during sorting: {}", metric_char_hits);
    spdlog::debug("Average number of characters compared per suffix comparison: {:.3f}", metric_avg_char_per_comp);

    spdlog::info("Finished sorting text suffixes in {:.3f} seconds", metric_sort_time);
}

/**
 * @brief Dumps internal sorting metrics to the JSON Lines benchmark log.
 * Acts as a class-level adapter for the global benchmark logger. It automatically 
 * aggregates all tracked private metrics and pushes them to disk using the Text-specific JSON schema.
 * Must be called at the very end of the sorting routine.
 *
 * @param seq_file    Base path/identifier for the input data, used to name the output log file.
 */

void TEXT_SORT::write_json(const std::string seq_file)
{
    write_sort_benchmark_jsonl(
        seq_file, // input path
        "Text_Sort", // config name
        metric_text_size, // text size
        metric_sort_time, // sort time
        metric_suffix_comps, // suffix comps
        metric_char_hits, // unit comps
        metric_avg_char_per_comp // avg unit per comp
    );   
}


/**
 * @brief Writes the suffix array of the text to file
 * 
 * @param [in] seq_file [string] The sequence file name is used to create output filename
 */

void TEXT_SORT::write_sa(const std::string seq_file)
{
    spdlog::stopwatch sw_write;

    std::string out_file = seq_file + ".sa";

    std::ofstream out(out_file);
    if (!out) {
        spdlog::error("Error opening {}", out_file);
        std::exit(EXIT_FAILURE);
    }

    for (size_t offset : suffix_array){
        out << offset << "\n";
    }

    out.close();

    auto sw_write_elapsed = sw_write.elapsed();
    spdlog::info("Finished writing suffix array to file in {:.3f} seconds", sw_write_elapsed.count());
}