/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <CLI11.hpp>
#include <cstdint>
#include <filesystem> // Note that this requires at least gcc 9
#include <limits>
#include "spdlog/spdlog.h"
#include "stats.h"

int main(int argc, char **argv) 
{
    CLI::App app("stats - Print stats related to the RLZ parsing.\n\nImplemented by Rahul Varki");
    std::string ref_file;
    std::string parse_file;
    bool bit = false;
    bool rlz_repair = false;
    int verbosity = 0;

    std::string version = "Version: 1.1.0";

    app.add_option("-r,--ref", ref_file, "Reference file")->required()->check(CLI::ExistingFile);
    app.add_option("-p,--parse", parse_file, "RLZ parse file")->required()->check(CLI::ExistingFile);
    app.add_flag("--repair", rlz_repair, "Set if used during compression");
    app.add_flag("--bit", bit, "Set if used during compression");
    app.add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);

    // Set version flag
    app.set_version_flag("--version", version);

    // Footer Updates
    app.footer("Example usage:\n"
               "  ./stats -r reference.fasta -p sequence.fasta.rlz [--bit] [--repair]\n");

    CLI11_PARSE(app, argc, argv);

    if (bit)
    {
        spdlog::info("Bit alphabet parsing was enabled");

        spdlog::info("Using the reference size to determine entry size");
        uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
        uintmax_t ref_size_bits = ref_size * 8;

        if (rlz_repair)
        {
            if (ref_size_bits < std::numeric_limits<int>::max()){
                spdlog::info("Assuming entries encoded with int");
                calculate_stats<int>(parse_file, ref_size_bits);
            }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                std::exit(1);
            }
            return 0;
        }

        // Entry size is determined by the size of the reference
        if (ref_size_bits <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); calculate_stats<uint8_t>(parse_file, ref_size_bits); }
        else if (ref_size_bits <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); calculate_stats<uint16_t>(parse_file, ref_size_bits); }
        else if (ref_size_bits <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); calculate_stats<uint32_t>(parse_file, ref_size_bits); }
        else if (ref_size_bits <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); calculate_stats<uint64_t>(parse_file, ref_size_bits); }
        else{
            spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
            std::exit(1);
        }

        return 0;
    }
    else
    {
        spdlog::info("Original alphabet sorting enabled");

        spdlog::info("Using the reference size to determine entry size");
        uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes

        if (rlz_repair)
        {
            if (ref_size < std::numeric_limits<int>::max()){
                spdlog::info("Assuming entries encoded with int");
                calculate_stats<int>(parse_file, ref_size);
            }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                std::exit(1);
            }
            return 0;
        }

        // Entries is determined by the size of the reference
        if (ref_size <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); calculate_stats<uint8_t>(parse_file, ref_size); }
        else if (ref_size <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); calculate_stats<uint16_t>(parse_file, ref_size); }
        else if (ref_size <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); calculate_stats<uint32_t>(parse_file, ref_size); }
        else if (ref_size <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); calculate_stats<uint64_t>(parse_file, ref_size); }
        else{
            spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
            std::exit(1);
        }

        return 0;
    }
}
