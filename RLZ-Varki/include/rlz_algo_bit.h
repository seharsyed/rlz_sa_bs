/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef RLZ_ALGO_BIT_H
#define RLZ_ALGO_BIT_H

#include "fm_wrapper.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/int_vector.hpp>
#include <fstream>
#include <filesystem>
#include <system_error>
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <stack>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

std::chrono::duration<double> backward_match_time{0.0};
std::chrono::duration<double> sa_time{0.0};
std::chrono::duration<double> serialize_time{0.0};

template <typename int_t>
class RLZ_BIT {
    public:
        std::string ref_file;
        sdsl::bit_vector ref_bit_array;

        RLZ_BIT(const std::string ref_file);
        ~RLZ_BIT();

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
        void load_reference_bit(const std::string& ref_file, sdsl::bit_vector& ref_bit_array);
        void load_reverse_reference_bit(const std::string& ref_file, sdsl::bit_vector& ref_bit_array);
        void calculate_occs(std::string& content, std::vector<size_t>& occs);
        void serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file);
        std::vector<std::tuple<int_t, int_t>> deserialize(const std::string& parse_file);

        void bits_to_str(sdsl::bit_vector bit_array, std::string prefix);
        void print_serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file);
};


/**
* @brief Constuctor of RLZ_BIT class.
* @param[in] ref_file [string] Path to reference file 
*/
template<typename int_t>
RLZ_BIT<int_t>::RLZ_BIT(const std::string ref_file): ref_file(ref_file){}

/**
* @brief Destructor of RLZ_BIT class.
*
* Currently does nothing.
*
*/
template<typename int_t>
RLZ_BIT<int_t>::~RLZ_BIT(){}

/**
* @brief Reads the "bits" of the reference and stores its content.
*
* Loads the reference content directly into an SDSL bit vector. Obtains the file size
* using filesystem metadata, performs a bulk binary read into an in-memory buffer, 
* and populates the bit vector efficiently from RAM.   
*
* @param[in] ref_file Path to the reference file
* @param[in] ref_bit_array Place where the reference content is stored.
* @return void
*/
template<typename int_t>
void RLZ_BIT<int_t>::load_reference_bit(const std::string& ref_file, sdsl::bit_vector& ref_bit_array)
{
    spdlog::debug("Reading reference content bits");
    spdlog::stopwatch sw_convert;

    // Getting size of the reference
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

    // Preloading size and loading reference
    std::vector<char> ref_buffer(ref_size);
    if (!ref.read(ref_buffer.data(), ref_size)) {
        spdlog::error("Error reading data from {}", ref_file);
        std::exit(EXIT_FAILURE);
    }
    ref.close();

    // Resize the SDSL bit array to hold the total number of bits required
    ref_bit_array.resize(ref_size * 8);

    // Populate the bit vector entirely from memory
    std::size_t bit_index = 0;
    for (char byte : ref_buffer) {
        for (int i = 7; i >= 0; --i) {
            ref_bit_array[bit_index++] = (byte >> i) & 1;
        }
    }

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished creating bit array in {:.3} seconds", sw_convert_elapsed.count());
}

/**
* @brief Reads the "bits" of the reversed reference and stores its content.
*
* Loads the reversed file content directly into an SDSL bit vector. Obtains the file size
* using filesystem metadata, performs a bulk binary read into an in-memory buffer, 
* and populates the bit vector in reverse order efficiently from RAM.   
*
* @param[in] ref_file Path to the reference file 
* @param[in] ref_bit_array Place where the reference content is stored.
* @return void
*/
template<typename int_t>
void RLZ_BIT<int_t>::load_reverse_reference_bit(const std::string& ref_file, sdsl::bit_vector& ref_bit_array)
{
    spdlog::debug("Reading reversed reference content bits");
    spdlog::stopwatch sw_convert;

    // Getting size of the reference
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
    std::vector<char> ref_buffer(ref_size);
    if (!ref.read(ref_buffer.data(), ref_size)) {
        spdlog::error("Error reading data from {}", ref_file);
        std::exit(EXIT_FAILURE);
    }
    ref.close();

    // Resize the SDSL bit array to hold the total number of bits required
    ref_bit_array.resize(ref_size * 8);

    // Directly loading the reference into buffer
    std::size_t bit_index = (ref_size * 8) - 1;
    for (char byte : ref_buffer) {
        for (int i = 7; i >= 0; --i) {
            ref_bit_array[bit_index--] = (byte >> i) & 1;
        }
    }

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished creating bit array in {:.3} seconds", sw_convert_elapsed.count());
}


/**
* @brief Parses the sequence file in relation to the reference file
*
* This function does the RLZ parsing of the sequence file. It currently works at the
* "psuedo" bit level. To clarify, it currently processes the string representation of the bits of the 
* sequence file. Working at bit level allows us to compress all types of files.
*
* RLZ algorithm tries to greedily find the longest sequence substring match within the reference.
* The RLZ parse in the end contains (pos, len) pairs
* in relation to the reference such that the sequence can be reconstructed from only the RLZ parse and the reference
* file. It is a O(n) algorithm. The size is the reference file + the RLZ parse.
*
* The algorithm implemented here is as follows.
* 1. Starting from the last bit of the reversed sequence file or sequence file chunk, check if bit matches the reversed reference
* (via backwards match with FM-index) 
* 2a. If match, check if next bit also matches (ex. 001. I know that 1 matches then check if 01 matches etc...)
* 2b. If match and end of sequence file or sequence file chunk, push current (pos,len) pair to parse vector
* 2c. If mismatch, push (prev pos, len - 1) to parse stack. Reset search from bit that caused mismatch.
*
* Push to parse stack since we process the string in reverse. Popping from stack gives correct order.
*
* @param [in] fm_index [rlz_fm_index_t] the fm-index of the reference
* @param [in] fm_support [FM_Wrapper] Utility object that allows us to do search and locate queries with fm-index.
* @param [in] occs [const std::vector<size_t>&] The compressed F column of the fm-index
* @param [in] seq_file [const std::string&] the sequence file
* @param [in] seq_parse_vec_vec [std::vector<std::vector<std::tuple<int_t, int_t>>>] empty RLZ parse vectors equal to number of threads
* @param [in] num_bits_to_process [size_t] the number of bits that should be processed. Useful for the OpenMP parallelization.
* @param [in] loop_iter [size_t] the loop iteration. Useful for OpenMP and making sure we are thread-safe.
* @param [in] num_threads [size_t] the total number of threads allocated.
* @param [in] max_len [size_t] the maximum length of the match
*
* @return void
*/
template<typename int_t>
void RLZ_BIT<int_t>::parse(const rlz_fm_index_t& fm_index,
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
    char byte;
    sdsl::bit_vector seq_bit_array;
    seq_bit_array.resize(8);
    size_t count = 0; // Keep track of how many chars processed
    while (sfile)
    {
        if (!retry) {  // Read a new character only if we're not retrying a char
            sfile.get(byte);
            count++;
            if (count % 10000 == 0){
                spdlog::debug("******** Thread {}: Processed {} bits in sequence file. ********", loop_iter + 1, count * 8);
            }
            if (sfile.eof()) break; // Exit if end of file
            if (count == num_char_to_process + 1) break; // Have processed all the characters this thread should process (The +1 because we increment before processing)
            for (int i = 7; i >= 0; --i) {
                seq_bit_array[7-i] = (byte >> i) & 1; // Here it gets stored in reverse order (0 pos = most sig bit, 1 pos = 2nd most sig bit)
            }
        }
        for (int i = 0; i < 8; i++)
        {
            char next_char = seq_bit_array[i] ? '1' : '0';

            pattern_len++;

            std::tuple<size_t,size_t> previous_ranges = std::make_tuple(prev_left, prev_right);
            auto back_start = std::chrono::high_resolution_clock::now();
            std::tuple<size_t,size_t> next_ranges = fm_support.backward_match(fm_index, occs, previous_ranges, next_char);
            auto back_end = std::chrono::high_resolution_clock::now();
            backward_match_time += back_end - back_start;
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
                sa_time += sa_end - sa_start;
                spdlog::trace("Mismatch: Putting ({},{}) at end of vector", adjusted_sa_pos, pattern_len);
                seq_parse_vec_vec[loop_iter].emplace_back(std::make_tuple(adjusted_sa_pos, pattern_len));
                prev_left = 0;
                prev_right = fm_index.bwt.size();
                next_left = 0;
                next_right = fm_index.bwt.size();
                pattern_len = 0;
                --i; // Do not increment if no perfect match
            }
            // If at the end we are still in a perfect match, we save what we have. 
            else if (i == 7 && count == num_char_to_process || i == 7 && sfile.peek() == EOF)
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
                sa_time += sa_end - sa_start;
                spdlog::trace("END: Putting ({},{}) at end of vector", adjusted_sa_pos, pattern_len);
                seq_parse_vec_vec[loop_iter].emplace_back(std::make_tuple(adjusted_sa_pos, pattern_len));
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
                sa_time += sa_end - sa_start;
                spdlog::trace("Max Len: Putting ({},{}) at end of vector", adjusted_sa_pos, pattern_len);
                seq_parse_vec_vec[loop_iter].emplace_back(std::make_tuple(adjusted_sa_pos, pattern_len));
                prev_left = 0;
                prev_right = fm_index.bwt.size();
                next_left = 0;
                next_right = fm_index.bwt.size();
                pattern_len = 0;
            }
            // Currently in a perfect match
            else{
                prev_left = next_left;
                prev_right = next_right;
            }
        }
    }
    sfile.close();
}


/**  
* @brief Builds the compressed F column of BWT matrix of reference "bits"
* @param [in] content [string] [string] The reference file "bits"
* @return void
*/
template<typename int_t>
void RLZ_BIT<int_t>::calculate_occs(std::string& ref_content, std::vector<size_t>& occs)
{
    spdlog::debug("Constructing compressed F column of reference 'bits' ");
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
* Creates a FM-index from the reversed reference bit array which we query using the reversed sequence bit array in order to simulate forward matching.
* We first convert the reversed reference bit array into its string representation so that we can create the FM-index.
* We query the index one bit at time from the reversed sequence bit array. When the sequence bit does not have a match, 
* we add the last matching ref position of the sequence and the length of the match to the parse. Then we 
* restart the match at the last mismatch position. The parse is ultimately
* stored in a vector in the correct order. The parse at the end is serialized to a file.
*
* @param [in] seq_file [const std::string&] The sequence file to compress
* @param [in] threads [int] The number of threads provided by the user.
* @param [in] max_len [size_t] The maximum length of the match 
*
* @return void
*
* @warning Providing multiple threads changes the output of the RLZ parse slightly. 
* Might create two phrases at chunk boundaries if phrase spans chunk boundary. For proper RLZ parse should run with 1 thread.
* * @note Supposedly cannot create a FM-index directly from bit array.
* Have to first convert into the reference bits into their string representation and then create the FM-index.
* Likely a bottleneck in the code as have to store a bit as a byte. [check if there is a way to build bit level FM-index]
*
*/
template<typename int_t>
void RLZ_BIT<int_t>::compress(const std::string& seq_file, int threads, size_t max_len)
{

    spdlog::stopwatch sw_compress;

    rlz_fm_index_t fm_index;
    std::string binary_reference_text;

    // Convert the reference bit array into its string representation
    for (size_t i = 0; i < ref_bit_array.size(); ++i) {
        binary_reference_text += (ref_bit_array[i] ? '1' : '0');
    }

    // Creates the FM-index
    spdlog::debug("Building FM-index of reversed reference 'bits'");
    spdlog::stopwatch sw_fm_index;
    construct_im(fm_index, binary_reference_text, 1);
    assert(binary_reference_text.size() + 1 == fm_index.bwt.size()); // SDSL should add the sentinel to the end
    auto sw_fm_index_elapsed = sw_fm_index.elapsed();
    spdlog::debug("Finished building FM-index in {:.3} seconds", sw_fm_index_elapsed.count());

    // Get the number of occurances of each bit in lexicographical order
    std::vector<size_t> occs(256, 0);
    calculate_occs(binary_reference_text, occs);

    FM_Wrapper fm_support;
    
    std::vector<std::vector<std::tuple<int_t, int_t>>> seq_parse_vec_vec(threads);

    std::ifstream sfile(seq_file, std::ios::binary | std::ios::ate); // Opens the file in binary mode and moves indicator to end
    if (!sfile) {
        spdlog::error("Error opening {}", seq_file);
        std::exit(EXIT_FAILURE);
    }
    size_t seq_size = static_cast<size_t>(sfile.tellg());
    sfile.close();

    size_t num_char_to_process;
    if (threads == 1){
        num_char_to_process = seq_size;
    }
    else{
        num_char_to_process = seq_size / (threads-1);  // Integer division (threads-1 to ensure last thread does not process more chars than num_char_to_process)
    }
    spdlog::debug("Each thread will process {} bits", num_char_to_process * 8);
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < threads; i++)
    {
        parse(fm_index, fm_support, occs, seq_file, seq_parse_vec_vec, num_char_to_process, i, threads, max_len);
    }

    spdlog::debug("Total time spent processing occurrences (s): {:.6f}", std::chrono::duration<double>(backward_match_time).count());
    spdlog::debug("Total time spent processing locations (s): {:.6f}", std::chrono::duration<double>(sa_time).count());

    // Store tuples of (pos,len) in correct order in vector
    size_t bits_stored = 0;
    std::vector<std::tuple<int_t, int_t>> seq_parse;
    // Can process the parse vectors sequentially since the first vector contains the parse of the start of the non-reversed sequence.
    for (int i = 0; i < threads; i++)
    {
        for (size_t j = 0; j < seq_parse_vec_vec[i].size(); j++)
        {
            bits_stored += std::get<1>(seq_parse_vec_vec[i][j]);
            seq_parse.emplace_back(seq_parse_vec_vec[i][j]);
            // spdlog::debug("Ref Pos: {}, Len: {}", std::get<0>(seq_parse.back()), std::get<1>(seq_parse.back()));
        }
    }

    spdlog::debug("The sequence file contained {} bits", seq_size * 8);
    spdlog::debug("The rlz parse encodes for {} bits", bits_stored);

    // Serialize the RLZ parse
    serialize(seq_parse, seq_file);

    auto sw_compress_elapsed = sw_compress.elapsed();
    spdlog::info("Compression finished in {:.3} seconds", sw_compress_elapsed.count());
}

/**
* @brief Serializes the parse of the sequence file
*
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We serialize the parse vector into binary file called seq_file_name.rlz
*
* File content of the .rlz file
* (uint64_t byte size num of pair) (int_t byte size pos) (int_t byte size len) (int_t byte size pos) (int_t byte size len) ...
* @param[in] seq_parse [std::vector<std::tuple<int_t, int_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
* @param[in] seq_file [std::string] the sequence file name
*
* @return void
*/
template<typename int_t>
void RLZ_BIT<int_t>::serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file)
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
* Decompress seq_file_name.rlz into tuple vector <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... > . 
* Return the vector.
*
* @param[in] parse_file [const std::string&] The parse filename (with .rlz extension)
*
* @return Return the vector.
*/
template<typename int_t>
std::vector<std::tuple<int_t, int_t>> RLZ_BIT<int_t>::deserialize(const std::string& parse_file)
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
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We read the ref position and length from each tuple in the sequence parse and get the corresponding bits from the reference bit array
* We then convert the bits into bytes and write the ASCII equivalent to file which sbould be the original sequence file.
*
* @param[in] parse_file [const std::string&] The parse filename (with .rlz extension)
*
* @warning might have to change int to long long int depending on size
*/
template<typename int_t>
void RLZ_BIT<int_t>::decompress(const std::string& parse_file)
{
    spdlog::stopwatch sw_decompress;

    std::vector<std::tuple<int_t, int_t>> seq_parse = deserialize(parse_file);
    
    size_t bit_size = 0;
    for (const auto& [pos, len] : seq_parse){
        bit_size += len;
    }

    spdlog::debug("The compressed sequence file encodes for {} bits", bit_size);

    // Resize the array to be equal to the number of bits to be stored
    sdsl::bit_vector seq_bit_array;
    seq_bit_array.resize(bit_size);

    // Get the sequence bits from the parse + reference bits
    size_t prev_pos = 0;
    for (const auto& [pos, len] : seq_parse){
        size_t curr_pos = prev_pos + len;
        size_t len_count = 0;
        for (size_t i = prev_pos; i < curr_pos; ++i){
            seq_bit_array[i] = ref_bit_array[pos + len_count];
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

    // Convert the sequence bits back into string
    std::string uncompressed_seq;

    uint8_t byte = 0;
    size_t bit_count = 0;

    // Write the sequence bits back into characters
    for (size_t i = 0; i < seq_bit_array.size(); ++i) {
        byte <<= 1; // Shift byte left by 1 bit
        byte |= seq_bit_array[i]; // Add the current bit to the byte
        ++bit_count;

        if (bit_count == 8) { // If 8 bits have been processed
            uncompressed_seq += static_cast<char>(byte); // Convert byte to char and append to result
            byte = 0; // Reset byte
            bit_count = 0; // Reset bit count
        }
    }

    output_file << uncompressed_seq;
    output_file.close();

    auto sw_decompress_elapsed = sw_decompress.elapsed();
    spdlog::info("Decompression finished in {:.3} seconds", sw_decompress_elapsed.count());
}


/**
* @brief Converts the bits to its corresponding character representation
*
* Prints out debug info about how many 0 and 1s are in the bit array.
* Writes the string representation of the bits to a file.
* Testing purposes only.
* * @param [in] bit_array [sdsl::bit_vector] the bit array of interest
* @param [in] ext [std::string] the extension of the file
*
*
*/
template<typename int_t>
void RLZ_BIT<int_t>::bits_to_str(sdsl::bit_vector bit_array, std::string ext)
{
    std::string bitstr;
    for (size_t i = 0; i < bit_array.size(); ++i) {
        bitstr += (bit_array[i] ? '1' : '0');
    }

    sdsl::rank_support_v<1> b_rank1(&bit_array);
    sdsl::rank_support_v<0> b_rank0(&bit_array);
    spdlog::debug("^^^^^^^^^^^^^^^^^^^^^^^^^");
    spdlog::debug("Number of 1s: {}", b_rank1(bit_array.size()));
    spdlog::debug("Number of 0s: {}", b_rank0(bit_array.size()));
    spdlog::debug("^^^^^^^^^^^^^^^^^^^^^^^^^");

    std::ofstream ofs("bits_to_str_debug" + ext);
    if (!ofs) {
        spdlog::error("Error opening {}", "bits_to_str_debug" + ext);
        std::exit(EXIT_FAILURE);
    }
    ofs << bitstr;
    ofs.close();
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
void RLZ_BIT<int_t>::print_serialize(const std::vector<std::tuple<int_t, int_t>>& seq_parse, const std::string& seq_file)
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

#endif  // RLZ_ALGO_BIT_H