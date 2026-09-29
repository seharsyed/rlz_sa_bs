/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <CLI11.hpp>
#include "sort_algo_bit.h"
#include "sort_algo_char.h"
#include "sort_algo_text.h"
#include "benchmark_logger.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <cstdint>
#include <filesystem> // Note that this requires at least gcc 9
#include <limits>
#include <omp.h>

template <typename int_t>
void run_rlz_bit_sort(const std::string& ref_file, const std::string& parse_file, const std::string& out_prefix, uintmax_t ref_size, bool naive, bool interval, bool induced, bool only_factor, bool resync, bool json)
{
    spdlog::info("Sorting is currently not supported for bit-level RLZ parses");
}

template <typename int_t>
void run_rlz_char_sort(const std::string& ref_file, const std::string& parse_file, const std::string& out_prefix, uintmax_t ref_size, bool naive, bool interval, bool induced, bool only_factor, bool resync, bool json)
{
    RLZ_CHAR_SORT<int_t> main_parser(ref_file, parse_file);

    if (naive){
        main_parser.log_memory_estimate("Naive_Sort", ref_size); 
        main_parser.sort_naive(resync);
        if (json) { main_parser.write_json(out_prefix, "Naive_Sort"); } 
    }
    else if (interval) {
        main_parser.log_memory_estimate("Interval_Sort", ref_size); 
        main_parser.sort_lcp_interval(resync); 
        if (json) { main_parser.write_json(out_prefix, "Interval_Sort"); }
    }
    else if (induced) {
        main_parser.log_memory_estimate("Induced_Sort", ref_size); 
        main_parser.sort_induced(resync);
        if (json) { main_parser.write_json(out_prefix, "Induced_Sort"); }
    }
    else if (only_factor) { 
        main_parser.log_memory_estimate("Factor_Only_Sort", ref_size);
        main_parser.sort_factors_only(resync); 
        if (json) { main_parser.write_json(out_prefix, "Factor_Only_Sort"); }
    }
    else { 
        spdlog::error("Compression option not chosen!");
        std::exit(1); 
    }

    main_parser.stream_sa_to_file(out_prefix);
}


void run_text_sort(const std::string& seq_file, const std::string& out_prefix, bool json)
{
    TEXT_SORT main_parser(seq_file);
    main_parser.build_sa();
    if (json) { main_parser.write_json(out_prefix); }
    main_parser.write_sa(out_prefix);
}



int main(int argc, char **argv) 
{
    CLI::App app("sort - Sorting RLZ factors.\n\nImplemented by Rahul Varki");

    std::string ref_file;
    std::string seq_file;
    std::string parse_file;
    std::string out_prefix;
    bool bit = false;
    bool rlz_repair = false;
    int verbosity = 0;
    bool json = false; // Write out stats in JSON lines format
    bool resync = false; // Apply resycnronization to the RLZ parse
    int threads = 1; // Number of threads to use for preprocessing OPENMP loop

    // Compression groups for RLZ
    bool naive = false;
    bool interval = false;
    bool induced = false;
    bool only_factor = false;

    std::string version = "Version: 1.1.0";

    // RLZ sorting
    auto* rlz_cmd = app.add_subcommand("rlz", "Sorting suffixes directly from RLZ factors");
    rlz_cmd->add_option("-r,--ref", ref_file, "Reference file")->required()->check(CLI::ExistingFile);
    rlz_cmd->add_option("-p,--parse", parse_file, "RLZ parse file to sort")->required()->check(CLI::ExistingFile);
    rlz_cmd->add_option("-o,--output", out_prefix, "Output prefix name (default is parse filename)");
    rlz_cmd->add_flag("--bit", bit, "Set if used during compression");
    rlz_cmd->add_flag("--repair", rlz_repair, "Set if used during compression");

    auto compression_group = rlz_cmd->add_option_group("Compression Method", "Choose exactly one compression method");
    compression_group->add_flag("--naive", naive, "Sort all suffixes with naive method");
    compression_group->add_flag("--interval", interval, "Sort all suffixes with LCP-interval awareness");
    compression_group->add_flag("--induced", induced, "Sort all suffixes using sorted complete factor backbone");
    compression_group->add_flag("--factors-only", only_factor, "Sort only the suffixes corresponding to complete factors");
    compression_group->require_option(1);

    rlz_cmd->add_flag("--resync", resync, "Apply resynchronization prior to sorting");
    rlz_cmd->add_flag("--json", json, "Output JSON Lines file containing sorting statistics");
    rlz_cmd->add_option("-t,--threads", threads, "Number of OpenMP threads to use during pre-processing step (default: 1)");
    rlz_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);
    
    // Text sorting
    auto* text_cmd = app.add_subcommand("text", "Sorting suffixes directly from text");
    text_cmd->add_option("-s,--seq", seq_file, "Sequence file to sort")->required()->check(CLI::ExistingFile);
    text_cmd->add_option("-o,--output", out_prefix, "Output prefix name (default is sequence filename)");
    text_cmd->add_flag("--json", json, "Output JSON Lines file containing sorting statistics");
    text_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);

    // Choose between rlz or text sorting
    app.require_subcommand(1, 1); 

    // Set version flag
    app.set_version_flag("--version", version);

    // Footer Updates
    app.footer("Example usage:\n"
               "  ./sort rlz -r reference.fasta -p sequence.fasta.rlz [--bit] [--repair] [--naive] [--interval] [--induced] [--factors-only] [--resync] [--json]\n"
               "  ./sort text -s sequence.fasta [--json]\n");

    CLI11_PARSE(app, argc, argv);

    omp_set_num_threads(threads); // Can only modify from rlz sort where it is used

    if (verbosity == 2) {
        spdlog::set_level(spdlog::level::trace);
    }
    else if (verbosity == 1){
        spdlog::set_level(spdlog::level::debug);
    }
    else if (verbosity == 0){ 
        spdlog::set_level(spdlog::level::info);
    }

    if (rlz_cmd->parsed())
    {
        if (out_prefix.empty()) { out_prefix = parse_file; }

        if (bit)
        {
            spdlog::info("Bit alphabet sorting enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference 
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
            uintmax_t ref_size_bits = ref_size * 8;

            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair)
            {
                if (ref_size_bits < std::numeric_limits<int>::max()){
                    spdlog::info("Assuming entries encoded with int");
                    run_rlz_bit_sort<int>(ref_file, parse_file, out_prefix, ref_size_bits, naive, interval, induced, only_factor, resync, json);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    std::exit(1);
                }
                return 0;
            }
            // Entry size is determined by the size of the reference
            if (ref_size_bits <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_rlz_bit_sort<uint8_t>(ref_file, parse_file, out_prefix, ref_size_bits, naive, interval, induced, only_factor, resync, json); }
            else if (ref_size_bits <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_rlz_bit_sort<uint16_t>(ref_file, parse_file, out_prefix, ref_size_bits, naive, interval, induced, only_factor, resync, json); }
            else if (ref_size_bits <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_rlz_bit_sort<uint32_t>(ref_file, parse_file, out_prefix, ref_size_bits, naive, interval, induced, only_factor, resync, json); }
            else if (ref_size_bits <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_rlz_bit_sort<uint64_t>(ref_file, parse_file, out_prefix, ref_size_bits, naive, interval, induced, only_factor, resync, json); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                std::exit(1);
            }

            return 0;
        }
        else
        {
            spdlog::info("Original alphabet sorting enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes

            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair)
            {
                if (ref_size < std::numeric_limits<int>::max()){
                    spdlog::info("Assuming entries encoded with int");
                    run_rlz_char_sort<int>(ref_file, parse_file, out_prefix, ref_size, naive, interval, induced, only_factor, resync, json);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    std::exit(1);
                }
                return 0;
            }
            // Entries is determined by the size of the reference
            if (ref_size <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_rlz_char_sort<uint8_t>(ref_file, parse_file, out_prefix, ref_size, naive, interval, induced, only_factor, resync, json); }
            else if (ref_size <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_rlz_char_sort<uint16_t>(ref_file, parse_file, out_prefix, ref_size, naive, interval, induced, only_factor, resync, json); }
            else if (ref_size <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_rlz_char_sort<uint32_t>(ref_file, parse_file, out_prefix, ref_size, naive, interval, induced, only_factor, resync, json); }
            else if (ref_size <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_rlz_char_sort<uint64_t>(ref_file, parse_file, out_prefix, ref_size, naive, interval, induced, only_factor, resync, json); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                std::exit(1);
            }

            return 0;
        }
    }
    else if (text_cmd->parsed())
    {
        if (out_prefix.empty()) { out_prefix = seq_file; }
        spdlog::info("Text sorting enabled");
        run_text_sort(seq_file, out_prefix, json);
        return 0;
    }
    else{ 
        spdlog::error("Neither rlz or text sorting. This condition should be impossible!"); 
        std::exit(1); 
    }

    return 0;
}