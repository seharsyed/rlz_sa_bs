/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <CLI11.hpp>
#include "rlz_algo_bit.h"
#include "rlz_algo_char.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <cstdint>
#include <filesystem> // Note that this requires at least gcc 9
#include <limits>
#include <omp.h>
#include <iostream>
#include <fstream>

template <typename int_t>
void run_bit_decompression(const std::string& ref_file, const std::string& parse_file)
{
    spdlog::info("The reference file provided: {}", ref_file);
    spdlog::info("The parse file provided: {}", parse_file);

    RLZ_BIT<int_t> main_parser(ref_file);
    main_parser.load_reference_bit(ref_file, main_parser.ref_bit_array);
    main_parser.decompress(parse_file);
}

template <typename int_t>
void run_bit_compression(const std::string& ref_file, const std::string& seq_file, int threads, size_t max_len)
{
    spdlog::debug("The reference file provided: {}", ref_file);
    spdlog::debug("The sequence file provided: {}", seq_file);

    RLZ_BIT<int_t> main_parser(ref_file);
    main_parser.load_reverse_reference_bit(ref_file, main_parser.ref_bit_array);
    main_parser.compress(seq_file, threads, max_len);
    
    spdlog::info("#############################################################");
    spdlog::info("File Size Statistics:");
    uintmax_t ref_size = std::filesystem::file_size(ref_file); //bytes
    uintmax_t seq_size = std::filesystem::file_size(seq_file); //bytes
    uintmax_t parse_size = std::filesystem::file_size(seq_file + ".rlz"); //bytes
    double comp_ratio = static_cast<double>(ref_size + parse_size) / 
                static_cast<double>(ref_size + seq_size) * 100;
    spdlog::info("The reference (ref) file provided: {} is {} bytes", ref_file, ref_size);
    spdlog::info("The sequence (seq) file provided: {} is {} bytes", seq_file, seq_size);
    spdlog::info("The parse (parse) file created: {} is {} bytes", seq_file + ".rlz", parse_size);
    spdlog::info("Compression ratio [((ref + parse)/(ref + seq)) * 100]: {:.3}%", comp_ratio);
}

template <typename int_t>
void run_char_decompression(const std::string& ref_file, const std::string& parse_file)
{
    spdlog::info("The reference file provided: {}", ref_file);
    spdlog::info("The parse file provided: {}", parse_file);
    
    RLZ_CHAR<int_t> main_parser(ref_file);
    main_parser.load_reference(ref_file, main_parser.ref_content);
    main_parser.decompress(parse_file);
}

template <typename int_t>
void run_char_compression(const std::string& ref_file, const std::string& seq_file, int threads, size_t max_len)
{
    spdlog::info("The reference file provided: {}", ref_file);
    spdlog::info("The sequence file provided: {}", seq_file);
    
    RLZ_CHAR<int_t> main_parser(ref_file);
    main_parser.load_reverse_reference(ref_file, main_parser.ref_content);

    spdlog::info("#############################################################");

    main_parser.compress(seq_file, threads, max_len);

/*
    spdlog::info("#############################################################");
    spdlog::info("File Size Statistics:");
    uintmax_t ref_size = std::filesystem::file_size(ref_file); //bytes
    uintmax_t seq_size = std::filesystem::file_size(seq_files[i]); //bytes
    uintmax_t parse_size = std::filesystem::file_size(seq_files[i] + ".rlz"); //bytes
    double comp_ratio = static_cast<double>(ref_size + parse_size) / 
                static_cast<double>(ref_size + seq_size) * 100;
    spdlog::info("The reference (ref) file provided: {} is {} bytes", ref_file, ref_size);
    spdlog::info("The sequence (seq) file provided: {} is {} bytes", seq_files[i], seq_size);
    spdlog::info("The parse (parse) file created: {} is {} bytes", seq_files[i] + ".rlz", parse_size);
    spdlog::info("Compression ratio [((ref + parse)/(ref + seq)) * 100]: {:.3}%", comp_ratio);*/
}

int main(int argc, char **argv) 
{
    CLI::App app("rlz - An implementation of RLZ that compresses a sequence file using a reference file.\n\nImplemented by Rahul Varki");

    std::string ref_file;
    std::string seq_file;
    std::string parse_file;
    bool decompress = false;
    int verbosity = 0;
    bool clean = false;
    bool bit = false;
    int threads = 1;
    size_t max_len = 0; // 0 means not set
    bool rlz_repair = false;
    std::string version = "Version: 1.1.0";
    
    // Compress Subcommand
    auto* compress_cmd = app.add_subcommand("compress", "Compress a sequence file using RLZ");
    compress_cmd->add_option("-r,--ref", ref_file, "Reference file")->required();
    compress_cmd->add_option("-s,--seq", seq_file, "Sequence file to compress")->required();
    compress_cmd->add_option("-t,--threads", threads, "Number of threads to use")->default_val(1);
    compress_cmd->add_option("-l, --len", max_len, "Maximum length a match can span")->check(CLI::Range(static_cast<size_t>(1), std::numeric_limits<size_t>::max()));
    compress_cmd->add_flag("--bit", bit, "Experimental: Set if ref lacks unique sequence chars");
    compress_cmd->add_flag("--repair", rlz_repair, "Set if running RLZ-RePair");
    compress_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);   

    // Decompress Subcommand
    auto* decompress_cmd = app.add_subcommand("decompress", "Decompress an RLZ parse file");
    decompress_cmd->add_option("-r,--ref", ref_file, "Reference file")->required();
    decompress_cmd->add_option("-p,--parse", parse_file, "RLZ parse file to decompress")->required();
    decompress_cmd->add_option("-l, --len", max_len, "Maximum length a match can span (must be specified if used for compression)")->check(CLI::Range(1UL, UINT64_MAX));
    decompress_cmd->add_flag("--bit", bit, "Experimental: Set if ref lacks unique sequence chars (must be specified if used for compression)");
    decompress_cmd->add_flag("--repair", rlz_repair, "Set if running RLZ-RePair");
    decompress_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);   

    // Choose between compression or decompression subcommand
    app.require_subcommand(1, 1); 

    // Set version flag
    app.set_version_flag("--version", version);

    // Footer Updates
    app.footer("Example usage:\n"
               "  ./rlz compress -r reference.fasta -s sequence.fasta [--bit] [--len [int] ]\n"
               "  ./rlz decompress -r reference.fasta -p sequence.fasta.rlz [--bit] [--len [int] ]\n");

    CLI11_PARSE(app, argc, argv);

    if (verbosity == 2) {
        spdlog::set_level(spdlog::level::trace);
    }
    else if (verbosity == 1){
        spdlog::set_level(spdlog::level::debug);
    }
    else if (verbosity == 0){ 
        spdlog::set_level(spdlog::level::info);
    }
    
    if (compress_cmd->parsed())
    {
        // Encode with "bit" level compression
        if (bit)
        {
            spdlog::info("Bit alphabet compression enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference 
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
            uintmax_t ref_size_bits = ref_size * 8;
            
            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair){
                if (ref_size_bits < std::numeric_limits<int>::max()){
                    spdlog::info("Encoding entries with int");
                    run_bit_compression<int>(ref_file, seq_file, threads, max_len * 8);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    exit(1);
                }
                return 0;
            }

            if (ref_size_bits <= UINT8_MAX) { spdlog::info("Encoding entries with uint8_t"); run_bit_compression<uint8_t>(ref_file, seq_file, threads, max_len * 8); }
            else if (ref_size_bits <= UINT16_MAX) { spdlog::info("Encoding entries with uint16_t"); run_bit_compression<uint16_t>(ref_file, seq_file, threads, max_len * 8); }
            else if (ref_size_bits <= UINT32_MAX) { spdlog::info("Encoding entries with uint32_t"); run_bit_compression<uint32_t>(ref_file, seq_file, threads, max_len * 8); }
            else if (ref_size_bits <= UINT64_MAX) { spdlog::info("Encoding entries with uint64_t"); run_bit_compression<uint64_t>(ref_file, seq_file, threads, max_len * 8); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                exit(1);
            }
        }
        // Encode with regular alphabet compression
        else
        {
            spdlog::info("Original alphabet compression enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes

            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair){
                if (ref_size < std::numeric_limits<int>::max()){
                    spdlog::info("Encoding entries with int");
                    run_char_compression<int>(ref_file, seq_file, threads, max_len);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    exit(1);
                }
                return 0;
            }
            
            if (ref_size <= UINT8_MAX) { spdlog::info("Encoding entries with uint8_t"); run_char_compression<uint8_t>(ref_file, seq_file, threads, max_len); }
            else if (ref_size <= UINT16_MAX) { spdlog::info("Encoding entries with uint16_t"); run_char_compression<uint16_t>(ref_file, seq_file, threads, max_len); }
            else if (ref_size <= UINT32_MAX) { spdlog::info("Encoding entries with uint32_t"); run_char_compression<uint32_t>(ref_file, seq_file, threads, max_len); }
            else if (ref_size <= UINT64_MAX) { spdlog::info("Encoding entries with uint64_t"); run_char_compression<uint64_t>(ref_file, seq_file, threads, max_len); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                exit(1);
            }
        }
    }
    else if (decompress_cmd->parsed())
    {
        // Encoded with "bit" level compression 
        if (bit)
        {
            spdlog::info("Bit alphabet decompression enabled");
            
            // Cannot use max len to determine size because position of match can be anywhere on the reference
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
            uintmax_t ref_size_bits = ref_size * 8;
            
            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair){
                if (ref_size_bits < std::numeric_limits<int>::max()){
                    spdlog::info("Encoding entries with int");
                    run_bit_decompression<int>(ref_file, parse_file);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    exit(1);
                }
                return 0;
            }

            // Entries are decoded dynamically by upper bound specified
            if (ref_size_bits <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_bit_decompression<uint8_t>(ref_file, parse_file); }
            else if (ref_size_bits <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_bit_decompression<uint16_t>(ref_file, parse_file); }
            else if (ref_size_bits <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_bit_decompression<uint32_t>(ref_file, parse_file); }
            else if (ref_size_bits <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_bit_decompression<uint64_t>(ref_file, parse_file); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                exit(1);
            }
        }
        // Encoded with regular alphabet compression
        else
        {
            spdlog::info("Original alphabet decompression enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
            
            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair){
                if (ref_size < std::numeric_limits<int>::max()){
                    spdlog::info("Encoding entries with int");
                    run_char_decompression<int>(ref_file, parse_file);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    exit(1);
                }
                return 0;
            }

            // Entries are decoded dynamically by upper bound specified
            if (ref_size <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_char_decompression<uint8_t>(ref_file, parse_file); }
            else if (ref_size <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_char_decompression<uint16_t>(ref_file, parse_file); }
            else if (ref_size <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_char_decompression<uint32_t>(ref_file, parse_file); }
            else if (ref_size <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_char_decompression<uint64_t>(ref_file, parse_file); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                exit(1);
            }
        }
    }
    else{ spdlog::error("Neither compression or decompression. This condition should be impossible!"); }
    
    return 0;
}
