/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef RLZ_ALGO_CHAR_H
#define RLZ_ALGO_CHAR_H

#include "fm_wrapper.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/int_vector.hpp>
#include <fstream>
#include <system_error>
#include <filesystem>
#include <vector>
#include <tuple>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <stack>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

std::chrono::duration<double> backward_match_time_char{0.0};
std::chrono::duration<double> sa_time_char{0.0};
std::chrono::duration<double> serialize_time_char{0.0};

template <typename int_t>
class RLZ_CHAR {
    public:
        std::string ref_file;
        std::string ref_content;

        RLZ_CHAR(const std::string ref_file);
        ~RLZ_CHAR();

        void compress(const std::string& seq_file, int threads, size_t max_len);

        void parse(const rlz_fm_index_t& fm_index,
            FM_Wrapper& fm_support,
            const std::vector<size_t>& occs,
            const std::string& seq_file,
            std::vector<std::vector<std::tuple<int_t, int_t>>>& seq_parse_vec_vec,
            size_t num_bits_to_process,
            size_t loop_iter,
            size_t num_threads,
            size_t max_len);

        void decompress(const std::string& parse_file);
        void load_reference(const std::string& ref_file, std::string& ref_content);
        void load_reverse_reference(const std::string& ref_file, std::string& ref_content);
        void calculate_occs(std::string& ref_content, std::vector<size_t>& occs);
        void serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file);
        std::vector<std::tuple<int_t, int_t>> deserialize(const std::string& parse_file);

        void print_serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file);
};

/**
* @brief Constructor of RLZ_CHAR class
* @param[in] ref_file [string] Path to reference file
*/
template<typename int_t>
RLZ_CHAR<int_t>::RLZ_CHAR(const std::string ref_file): ref_file(ref_file){}

/**
* @brief Destructor of RLZ_CHAR class.
*
* Currently does nothing.
*
*/
template<typename int_t>
RLZ_CHAR<int_t>::~RLZ_CHAR(){}

/**
* @brief Reads the reference and stores its content
*
* Loads the reference directly into a string. Obtains the file size
* using filesystem metadata and performs an efficient bulk binary read.   
*
* @param[in] ref_file Path to the reference file
* @param[in] ref_content Place where the reference content is stored.
* @return void
*/
template<typename int_t>
void RLZ_CHAR<int_t>::load_reference(const std::string& ref_file, std::string& ref_content)
{
    spdlog::debug("Reading reference content");
    spdlog::stopwatch sw_convert;

    // Getting size of reference
    std::error_code ec;
    uintmax_t ref_size = std::filesystem::file_size(ref_file, ec);
    if (ec) {
        spdlog::error("Error getting file size for {}: {}", ref_file, ec.message());
        std::exit(EXIT_FAILURE);
    }

    // Opening reference
    std::ifstream ref(ref_file, std::ios::binary);
    if (!ref) {
        spdlog::error("Error opening {}", ref_file);
        std::exit(EXIT_FAILURE);
    }

    // Preloading size of reference buffer
    ref_content.resize(ref_size);

    // Directly loading the reference into buffer
    if (!ref.read(&ref_content[0], ref_size)) {
        spdlog::error("Error reading data from {}", ref_file);
        std::exit(EXIT_FAILURE);
    }
    ref.close();

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished reading file in {:.3} seconds", sw_convert_elapsed.count());
}


/**
* @brief Reads the reference and stores the reversed content.
*
* Loads the reversed reference directly into a string. Obtains the file size
* using filesystem metadata, performs a bulk binary read, and reverses the 
* sequence efficiently in memory.
*
* @param[in] ref_file Path to the reference file
* @param[in] ref_content Place where reversed content is stored.
* @return void
*/
template<typename int_t>
void RLZ_CHAR<int_t>::load_reverse_reference(const std::string& ref_file, std::string& ref_content)
{
    spdlog::debug("Reading reversed reference content");
    spdlog::stopwatch sw_convert;

    // Getting size of reference
    std::error_code ec;
    uintmax_t ref_size = std::filesystem::file_size(ref_file, ec);
    if (ec) {
        spdlog::error("Error getting file size for {}: {}", ref_file, ec.message());
        std::exit(EXIT_FAILURE);
    }

    // Opening reference
    std::ifstream ref(ref_file, std::ios::binary);
    if (!ref) {
        spdlog::error("Error opening {}", ref_file);
        std::exit(EXIT_FAILURE);
    }

    // Preloading size of reference buffer
    ref_content.resize(ref_size);

    // Directly loading the reference into buffer
    if (!ref.read(&ref_content[0], ref_size)) {
        spdlog::error("Error reading data from {}", ref_file);
        std::exit(EXIT_FAILURE);
    }
    ref.close();

    // Reversing reference
    std::reverse(ref_content.begin(), ref_content.end());

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished reading file in {:.3} seconds", sw_convert_elapsed.count());
}

/**
* @brief Parses the sequence file in relation to the reference file
*
* This function does the RLZ_CHAR parsing of the sequence file. It parses relative to the original file content of the reference.
* Fails if there is a character in the sequence that is not present in the reference.
*
* RLZ algorithm tries to greedily find the longest sequence substring match within the reference.
* The RLZ parse in the end contains (pos, len) pairs
* in relation to the reference such that the sequence can be reconstructed from only the RLZ parse and the reference
* file. It is a O(n) algorithm. The size is the reference file + the RLZ parse.
*
* The algorithm implemented here is as follows.
* 1. Starting from the last char of the reversed sequence file or sequence file chunk, check if char matches the reversed reference
* (via backwards match with FM-index) 
* 2a. If match, check if next char also matches (ex. aab. I know that b matches then check if ab matches etc...)
* 2b. If match and end of sequence file or sequence file chunk, push current (pos,len) pair to parse vector
* 2c. If mismatch, push (prev pos, len - 1) to parse stack. Reset search from char that caused mismatch.
*
*
* @param [in] fm_index [rlz_fm_index_t] the fm-index of the reference
* @param [in] fm_support [FM_Wrapper] Utility object that allows us to do search and locate queries with fm-index.
* @param [in] occs [const std::vector<size_t>&] The compressed F column of the fm-index
* @param [in] seq_file [const std::string&] the sequence file
* @param [in] seq_parse_vec_vec [std::vector<std::vector<std::tuple<int_t, int_t>>>] empty RLZ_CHAR parse vectors equal to number of threads
* @param [in] num_char_to_process [size_t] the number of chars that should be processed. Useful for the OpenMP parallelization.
* @param [in] loop_iter [size_t] the loop iteration. Useful for OpenMP and making sure we are thread-safe.
* @param [in] num_threads [size_t] the total number of threads allocated.
* @param [in] max_len [size_t] the maximum length of the match
*
* @return void
*/
template<typename int_t>
void RLZ_CHAR<int_t>::parse(const rlz_fm_index_t& fm_index,
        FM_Wrapper& fm_support,
        const std::vector<size_t>& occs, 
        const std::string& seq_file,
        std::vector<std::vector<std::tuple<int_t, int_t>>>& seq_parse_vec_vec,
        size_t num_char_to_process,
        size_t loop_iter,
        size_t num_threads,
        size_t max_len)
{
    size_t pattern_len = 0;
    size_t prev_left = 0;
    size_t prev_right = fm_index.bwt.size();
    size_t next_left = 0;
    size_t next_right = fm_index.bwt.size();

    std::ifstream sfile(seq_file, std::ios::binary | std::ios::ate);
    if (!sfile) {
        spdlog::error("Error opening {}", seq_file);
        std::exit(EXIT_FAILURE);
    }
    
    size_t seq_size = static_cast<long long int>(sfile.tellg());
    size_t start_loc = loop_iter * num_char_to_process;

    // Move file this many characters
    sfile.seekg(start_loc, std::ios::beg);
    
    // Process the file in reverse for backwards matching with FM-index.
    bool retry = false;
    char next_char;
    size_t count = 0; // Keep track of how many characters processed
    while (sfile)
    {
        if (!retry) {  // Read a new character only if we're not retrying a char
            sfile.get(next_char);
            count++;
            if (count % 10000 == 0){
                spdlog::debug("******** Thread {}: Processed {} unique chars in sequence file. ********", loop_iter + 1, count);
            }
            if (sfile.eof()) break; // Exit if end of file
            if (count == num_char_to_process + 1) break; // Have processed all the characters this thread should process (The +1 because we increment before processing)
        }

        pattern_len++;

        std::tuple<size_t,size_t> previous_ranges = std::make_tuple(prev_left, prev_right);
        auto back_start = std::chrono::high_resolution_clock::now();
        std::tuple<size_t,size_t> next_ranges = fm_support.backward_match(fm_index, occs, previous_ranges, next_char);
        auto back_end = std::chrono::high_resolution_clock::now();
        backward_match_time_char += back_end - back_start;
        next_left = std::get<0>(next_ranges);
        next_right = std::get<1>(next_ranges);

        // If same then that means no perfect match so we reset.
        if (next_left == next_right){
            pattern_len--; // -1 due to not matching the last character successfully
            auto sa_start = std::chrono::high_resolution_clock::now();
            size_t sa_pos = fm_support.get_suffix_array_value(fm_index, prev_left);
            // spdlog::trace("Reversed SA pos: {}", sa_pos);
            size_t mirrored_sa_pos = fm_index.bwt.size() - 1 - sa_pos; // 0 based involution formula of sa position to correct for the reverse string matching (will give pos in ref where pattern ends)
            // spdlog::trace("Forward SA pos: {}", mirrored_sa_pos);
            // Note: Normally would subtract pattern_len - 1 from mirrored position since mirrored position is match of len 1.
            // However, $ is added to the end of the reversed text which means it appears at the front of the original text. 
            // This $ is purely virtual, so it creates a 1 offset for all matches which we account for by subtracting pattern_len
            size_t adjusted_sa_pos = mirrored_sa_pos - pattern_len; 
            // spdlog::trace("Adjusted SA pos: {}", adjusted_sa_pos);
            auto sa_end = std::chrono::high_resolution_clock::now();
            sa_time_char += sa_end - sa_start;
            spdlog::trace("Mismatch: Putting ({},{}) at end of vector", adjusted_sa_pos, pattern_len);
            seq_parse_vec_vec[loop_iter].emplace_back(std::make_tuple(adjusted_sa_pos, pattern_len));
            prev_left = 0;
            prev_right = fm_index.bwt.size();
            next_left = 0;
            next_right = fm_index.bwt.size();
            pattern_len = 0;
            retry = true;
        }
        // If at the end we are still in a perfect match, we save what we have. 
        else if (sfile.peek() == EOF || count == num_char_to_process)
        {
            auto sa_start = std::chrono::high_resolution_clock::now();
            size_t sa_pos = fm_support.get_suffix_array_value(fm_index, next_left);
            // spdlog::trace("Reversed SA pos: {}", sa_pos);
            size_t mirrored_sa_pos = fm_index.bwt.size() - 1 - sa_pos;
            // spdlog::trace("Forward SA pos: {}", mirrored_sa_pos);
            // Note: Normally would subtract pattern_len - 1 from mirrored position since mirrored position is match of len 1.
            // However, $ is added to the end of the reversed text which means it appears at the front of the original text. 
            // This $ is purely virtual, so it creates a 1 offset for all matches which we account for by subtracting pattern_len
            size_t adjusted_sa_pos = mirrored_sa_pos - pattern_len;
            // spdlog::trace("Adjusted SA pos: {}", adjusted_sa_pos);
            auto sa_end = std::chrono::high_resolution_clock::now();
            sa_time_char += sa_end - sa_start;
            spdlog::trace("END: Putting ({},{}) at end of vector", adjusted_sa_pos, pattern_len);
            seq_parse_vec_vec[loop_iter].emplace_back(std::make_tuple(adjusted_sa_pos, pattern_len));
            retry = false;
        }
        // If a match has reached the user-defined max len
        else if (max_len > 0 && pattern_len == max_len)
        {
            auto sa_start = std::chrono::high_resolution_clock::now();
            size_t sa_pos = fm_support.get_suffix_array_value(fm_index, next_left);
            // spdlog::trace("Reversed SA pos: {}", sa_pos);
            size_t mirrored_sa_pos = fm_index.bwt.size() - 1 - sa_pos;
            // spdlog::trace("Forward SA pos: {}", mirrored_sa_pos);
            // Note: Normally would subtract pattern_len - 1 from mirrored position since mirrored position is match of len 1.
            // However, $ is added to the end of the reversed text which means it appears at the front of the original text. 
            // This $ is purely virtual, so it creates a 1 offset for all matches which we account for by subtracting pattern_len
            size_t adjusted_sa_pos = mirrored_sa_pos - pattern_len;
            // spdlog::trace("Adjusted SA pos: {}", adjusted_sa_pos);
            auto sa_end = std::chrono::high_resolution_clock::now();
            sa_time_char += sa_end - sa_start;
            spdlog::trace("Max Len: Putting ({},{}) at end of vector", adjusted_sa_pos, pattern_len);
            seq_parse_vec_vec[loop_iter].emplace_back(std::make_tuple(adjusted_sa_pos, pattern_len));
            prev_left = 0;
            prev_right = fm_index.bwt.size();
            next_left = 0;
            next_right = fm_index.bwt.size();
            pattern_len = 0;
            retry = false;
        }
        // Currently in a perfect match
        else{
            prev_left = next_left;
            prev_right = next_right;
            retry = false;
        }
    }
    sfile.close();
}


/** 
* @brief Builds the compressed F column of BWT matrix of reference
* @param [in] ref_content [string] The reference file 
* @return void
*/
template<typename int_t>
void RLZ_CHAR<int_t>::calculate_occs(std::string& ref_content, std::vector<size_t>& occs)
{
    spdlog::debug("Constructing compressed F column of reference");
    spdlog::stopwatch sw_occs;

    for (char c : ref_content) {
        occs[static_cast<unsigned char>(c)]++;
    }
    size_t running_total = 0;
    for (int i = 0; i < 256; ++i) {
        size_t current_frequency = occs[i];
        occs[i] = running_total;      
        running_total += current_frequency; 
    }

    auto sw_occs_elapsed = sw_occs.elapsed();
    spdlog::debug("Finished building compressed F column in {:.3} seconds", sw_occs_elapsed.count());
}


/**
* @brief Compresses the sequence file in relation to the reference file.
*
* Creates a FM-index from the reversed reference string which we query using the reversed sequence string in order to simulate forward matching.
* We first create the FM-index from the reveresed string representation of the reference.
* We query the index one char at time from the reversed sequence string. When the sequence char does not have a match, 
* we add the last matching ref position of the sequence and the length of the match to the parse. Then we 
* restart the match at the last mismatch position. The parse is ultimately
* stored in a vector in the correct order. The parse at the end is serialized to a file.
*
* @param [in] seq_file [const std::string&] The sequence file to compress
* @param [in] threads [int] The number of threads provided by the user.
* @param [in] max_len [size_t] The maximum length of a match
*
* @return void
*
* @warning Providing multiple threads changes the output of the RLZ_CHAR parse slightly. 
* Might create two phrases at chunk boundaries if phrase spans chunk boundary. For proper RLZ_CHAR parse should run with 1 thread.
* * @warning Will fail if the sequence file contains a char not present in the reference file
*
*/
template<typename int_t>
void RLZ_CHAR<int_t>::compress(const std::string& seq_file, int threads, size_t max_len)
{

    rlz_fm_index_t fm_index;
    
    // Creates the FM-index
    spdlog::debug("Building FM-index of reversed reference");
    spdlog::stopwatch sw_fm_index;
    construct_im(fm_index, ref_content, 1);
    assert(ref_content.size() + 1 == fm_index.bwt.size()); // SDSL should add the sentinel to the end
    auto sw_fm_index_elapsed = sw_fm_index.elapsed();
    spdlog::debug("Finished building FM-index in {:.3} seconds", sw_fm_index_elapsed.count());

    // Get the number of occurances of each char in lexicographical order
    std::vector<size_t> occs(256, 0);
    calculate_occs(ref_content, occs);

    FM_Wrapper fm_support;

    std::ifstream s(seq_file);
    if (!s) {
        throw std::runtime_error("Cannot open file: " + std::string(seq_file));
    }
    std::vector<std::string> seq_files;
    std::string p;

    while (s >> p)
        seq_files.push_back(p);

    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto start = high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < seq_files.size(); i++) {

        std::vector<std::vector<std::tuple<int_t, int_t>>> seq_parse_vec_vec(threads);

        std::ifstream sfile(seq_files[i], std::ios::binary | std::ios::ate); // Opens the file in binary mode and moves indicator to end
        if (!sfile) {
            spdlog::error("Error opening {}", seq_files[i]);
            std::exit(EXIT_FAILURE);
        }
        size_t seq_size = static_cast<size_t>(sfile.tellg());
        sfile.close();

        size_t num_char_to_process;
        if (threads == 1){
            num_char_to_process = seq_size;
        }
        else{
            num_char_to_process = seq_size / (threads-1);  // Integer division
        }
        //spdlog::debug("Each thread will process {} characters", num_char_to_process);
        //spdlog::stopwatch sw_compress;

        #pragma omp parallel for schedule(static)
        for (size_t j = 0; j < threads; j++)
        {
            parse(fm_index, fm_support, occs, seq_files[i], seq_parse_vec_vec, num_char_to_process, j, threads, max_len);
        }

        //spdlog::debug("Total time spent processing occurrences (s): {:.6f}", std::chrono::duration<double>(backward_match_time_char).count());
        //spdlog::debug("Total time spent processing locations (s): {:.6f}", std::chrono::duration<double>(sa_time_char).count());

        // Store tuples of (pos,len) in correct order in vector
        size_t chars_stored = 0;
        std::vector<std::tuple<int_t, int_t>> seq_parse;
        // Can process the parse vectors sequentially since the first vector contains the parse of the start of the non-reversed sequence.
        for (int j = 0; j < threads; j++)
        {
            for (size_t k = 0; k < seq_parse_vec_vec[j].size(); k++)
            {
                chars_stored += std::get<1>(seq_parse_vec_vec[j][k]);
                seq_parse.emplace_back(seq_parse_vec_vec[j][k]);
            }
        }

        //spdlog::debug("The sequence file contained {} chars", seq_size);
        //spdlog::debug("The rlz parse encodes for {} chars", chars_stored);

        // Serialize the RLZ parse
        serialize(seq_parse, seq_files[i]);
        //auto sw_compress_elapsed = sw_compress.elapsed();
        //spdlog::info("Compression finished in {:.3} seconds", sw_compress_elapsed.count());
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = duration_cast<milliseconds>(end - start).count();
    std::cout << "Total Time : " << time << "ms" << std::endl;
    

}

/**
* @brief Serializes the parse of the sequence file
*
* The sequence parse contains tuples (ref pos, size) that can reconstruct the sequence file given the reference.
* We serialize the parse vector into binary file called seq_file_name.rlz
*
* File content of the .rlz file
* (uint64_t byte: size num of pair) (int_t byte: size pos) (int_t byte: size len) (int_t byte: size pos) (int_t byte: size len) ...
* @param[in] seq_parse [std::vector<std::tuple<int_t, int_t>>] The parse of the seq <(ref pos,len),(ref pos,len),(ref pos,len)... >
* @param[in] seq_file [std::string] the sequence file name
*
* @return void
*/
template<typename int_t>
void RLZ_CHAR<int_t>::serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file)
{
    spdlog::debug("Serializing RLZ parse");
    spdlog::stopwatch sw_serialize;

    std::ofstream ofs(seq_file + ".rlz", std::ios::binary);
    if (!ofs) {
        spdlog::error("Error opening {}", seq_file + ".rlz");
        std::exit(EXIT_FAILURE);
    }
    uint64_t size = seq_parse.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(uint64_t));
    for (size_t i = 0; i < size; i++)
    {
        ofs.write(reinterpret_cast<const char*>(&std::get<0>(seq_parse[i])), sizeof(int_t));
        ofs.write(reinterpret_cast<const char*>(&std::get<1>(seq_parse[i])), sizeof(int_t));
    }
    ofs.close();

    auto sw_serialize_elapsed = sw_serialize.elapsed();
    spdlog::debug("Serializing finished in {:.6f} seconds", sw_serialize_elapsed.count());
}


/**
* @brief Deserializes the parse of the sequence file
*
* Decompress seq_file_name.rlz into tuple vector <(ref pos,len),(ref pos,len),(ref pos,len)... > . 
* Return the vector.
* @param[in] parse_file [const std::string&] The parse filename (with .rlz extension)
*
* @return Return the vector.
*/
template<typename int_t>
std::vector<std::tuple<int_t, int_t>> RLZ_CHAR<int_t>::deserialize(const std::string& parse_file)
{
    spdlog::debug("Deserializing RLZ parse");
    spdlog::stopwatch sw_deserialize;

    std::ifstream ifs(parse_file, std::ios::binary);
    if (!ifs) {
        spdlog::error("Error opening {}", parse_file);
        std::exit(EXIT_FAILURE);
    }
    uint64_t size;
    std::vector<std::tuple<int_t, int_t>> seq_parse;

    ifs.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    seq_parse.reserve(size);
    std::tuple<int_t, int_t> elem;
    int_t val;

    for (size_t i = 0; i < 2 * size; i++)
    {
        ifs.read(reinterpret_cast<char*>(&val), sizeof(int_t));
        if (i % 2 == 0){
            std::get<0>(elem) = val;
        }
        else{
            std::get<1>(elem) = val;
            seq_parse.emplace_back(elem);
        }
    }
    ifs.close();

    auto sw_deserialize_elapsed = sw_deserialize.elapsed();
    spdlog::debug("Deserializing finished in {:.6f} seconds", sw_deserialize_elapsed.count());

    return seq_parse;
}

/**
* @brief Decompresses the sequence parse back to the original sequence file
*
* The sequence parse contains tuples (ref pos, size) that can reconstruct the sequence file given the reference.
* We read the ref position and length from each tuple in the sequence parse and get the corresponding chars from the reference to reconstruct the sequence
*
* @param[in] parse_file [const std::string&] The parse filename (with .rlz extension)
*
* @warning might have to change int to long long int depending on size
*/
template<typename int_t>
void RLZ_CHAR<int_t>::decompress(const std::string& parse_file)
{
    spdlog::stopwatch sw_decompress;

    std::vector<std::tuple<int_t, int_t>> seq_parse = deserialize(parse_file);
    
    size_t char_size = 0;
    for (const auto& [pos, len] : seq_parse){
        char_size += len;
    }

    spdlog::debug("The compressed sequence file encodes for {} chars", char_size);

    // Resize the array to be equal to the number of bits to be stored
    std::string seq_content;
    seq_content.resize(char_size);

    // Get the sequence chars from the parse + reference chars
    size_t prev_pos = 0;
    for (const auto& [pos, len] : seq_parse){
        size_t curr_pos = prev_pos + len;
        size_t len_count = 0;
        for (size_t i = prev_pos; i < curr_pos; ++i){
            seq_content[i] = ref_content[pos + len_count];
            len_count++;
        }
        prev_pos = curr_pos;
    }

    // Remove extension from parse file
    size_t last_dot = parse_file.find_last_of(".");
    std::string seq_file = parse_file.substr(0, last_dot);
    
    std::ofstream output_file(seq_file + ".out");
    if (!output_file) {
        spdlog::error("Error opening {}", seq_file + ".out");
        std::exit(EXIT_FAILURE);
    }

    output_file << seq_content;
    output_file.close();

    auto sw_decompress_elapsed = sw_decompress.elapsed();
    spdlog::info("Decompression finished in {:.3} seconds", sw_decompress_elapsed.count());
}


/**
* @brief Write the non-binary serialization of the sequence parse to a file.
*
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We write the non-binary serialization to a file. Testing purposes only.
* * @param[in] seq_parse [std::vector<std::tuple<int_t, int_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
* @param[in] seq_file [const std::string&] The sequence filename
*
* @return void
*/
template<typename int_t>
void RLZ_CHAR<int_t>::print_serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file)
{
    std::ofstream ofs(seq_file + ".readable.rlz");
    if (!ofs) {
        spdlog::error("Error opening {}", seq_file + ".readable.rlz");
        std::exit(EXIT_FAILURE);
    }

    for (const auto& [pos, len] : seq_parse){
        ofs << "Position: " << pos << "  Length: " << len << std::endl;
    }

    ofs.close();
}

#endif  // RLZ_ALGO_CHAR_H
