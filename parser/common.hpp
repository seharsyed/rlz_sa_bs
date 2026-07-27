#pragma once

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "parser.hpp"

using Symbol = unsigned char;
using SAType = unsigned int;
using Triples =
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;

struct Args {
    std::string reference;
    std::string suffix_array;
    std::string filenames;
    std::size_t p = 64;
};

struct LoadedData {
    std::vector<std::string> files;
    std::vector<Symbol> ref;
    std::vector<SAType> sa;
};

struct FileResult {
    std::string filename;
    std::size_t input_bytes = 0;
    double rlz_ms = 0.0;
    std::size_t phrases = 0;
    std::size_t cache_hits = 0;
    std::size_t cache_misses = 0;
    double cache_hit_rate = 0.0;
    std::size_t cache_entries = 0;
    std::size_t table_slots = 0;
    std::size_t max_probe_cluster = 0;
    double load_factor = 0.0;
    double cache_mb = 0.0;
    double peak_rss_mb = 0.0;
};

struct RunResults {
    std::size_t processed_files = 0;
    std::size_t skipped_files = 0;
    std::size_t total_input_bytes = 0;
    std::size_t total_phrases = 0;
    double total_rlz_ms = 0.0;

    void add(const FileResult& result) {
        ++processed_files;
        total_input_bytes += result.input_bytes;
        total_phrases += result.phrases;
        total_rlz_ms += result.rlz_ms;
    }
};

inline void print_usage(const char* program) {
    std::cerr
        << "Usage:\n"
        << "  " << program
        << " --reference PATH --filenames PATH "
        << "[--suffix-array PATH] [--p N]\n\n"
        << "If --suffix-array is omitted, <reference>.sa is used.\n"
        << "Default --p is 64.\n";
}

inline Args parse_args(int argc, char** argv) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];

        auto require_value = [&](const std::string& name) {
            if (i + 1 >= argc) {
                throw std::runtime_error(
                    "missing value for " + name
                );
            }

            return std::string(argv[++i]);
        };

        if (option == "--reference") {
            args.reference = require_value(option);
        } else if (option == "--suffix-array") {
            args.suffix_array = require_value(option);
        } else if (option == "--filenames") {
            args.filenames = require_value(option);
        } else if (option == "--p") {
            args.p = static_cast<std::size_t>(
                std::stoull(require_value(option))
            );
        } else if (option == "--help" || option == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error(
                "unknown argument: " + option
            );
        }
    }

    if (args.reference.empty()) {
        throw std::runtime_error("missing --reference");
    }

    if (args.filenames.empty()) {
        throw std::runtime_error("missing --filenames");
    }

    if (args.suffix_array.empty()) {
        args.suffix_array = args.reference + ".sa";
    }

    if (args.p == 0) {
        throw std::runtime_error(
            "--p must be greater than zero"
        );
    }

    return args;
}

inline std::string trim(std::string value) {
    while (
        !value.empty() &&
        (
            value.back() == '\n' ||
            value.back() == '\r' ||
            value.back() == ' ' ||
            value.back() == '\t'
        )
    ) {
        value.pop_back();
    }

    std::size_t first = 0;

    while (
        first < value.size() &&
        (
            value[first] == ' ' ||
            value[first] == '\t'
        )
    ) {
        ++first;
    }

    return value.substr(first);
}

inline std::vector<std::string> read_filename_list(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "cannot open filenames file: " + path
        );
    }

    std::vector<std::string> files;
    std::string line;

    while (std::getline(input, line)) {
        line = trim(line);

        if (!line.empty() && line.front() != '#') {
            files.push_back(line);
        }
    }

    if (files.empty()) {
        throw std::runtime_error(
            "filenames file is empty: " + path
        );
    }

    return files;
}

inline LoadedData load_data(const Args& args) {
    LoadedData data;

    std::cerr
        << "Loading input file list..."
        << std::endl;

    data.files = read_filename_list(args.filenames);

    std::cerr
        << "Input files loaded: "
        << data.files.size()
        << std::endl;

    std::cerr
        << "Loading reference..."
        << std::endl;

    data.ref = read_file<Symbol>(
        args.reference.c_str()
    );

    std::cerr
        << "Reference loaded: "
        << data.ref.size()
        << " bytes"
        << std::endl;

    std::cerr
        << "Loading suffix array..."
        << std::endl;

    data.sa = read_file<SAType>(
        args.suffix_array.c_str()
    );

    std::cerr
        << "Suffix array loaded: "
        << data.sa.size()
        << " entries"
        << std::endl;

    if (data.ref.empty() || data.sa.empty()) {
        throw std::runtime_error(
            "reference and suffix array must be non-empty"
        );
    }

    if (data.sa.size() != data.ref.size()) {
        std::cerr
            << "Warning: SA entries ("
            << data.sa.size()
            << ") != reference bytes ("
            << data.ref.size()
            << ")"
            << std::endl;
    }

    return data;
}

template <typename Function>
double measure_ms(Function&& function) {
    const auto start =
        std::chrono::steady_clock::now();

    function();

    const auto end =
        std::chrono::steady_clock::now();

    return std::chrono::duration<double, std::milli>(
        end - start
    ).count();
}

inline double peak_rss_mb() {
#if defined(__unix__) || defined(__APPLE__)
    struct rusage usage{};

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<double>(
            usage.ru_maxrss
        ) / (1024.0 * 1024.0);
#else
        return static_cast<double>(
            usage.ru_maxrss
        ) / 1024.0;
#endif
    }
#endif

    return 0.0;
}

template <typename CacheInfo>
void print_parser_ready(
    const CacheInfo& cache_info
) {
    std::cerr
        << "LP parser created"
        << std::endl;

    std::cerr
        << "Initial table slots: "
        << cache_info.table_slots
        << std::endl;
}

inline void print_run_details(
    const Args& args,
    const LoadedData& data
) {
    std::cout
        << "reference,"
        << args.reference
        << '\n';

    std::cout
        << "suffix_array,"
        << args.suffix_array
        << '\n';

    std::cout
        << "files_requested,"
        << data.files.size()
        << '\n';

    std::cout
        << "p,"
        << args.p
        << '\n';

    std::cout
        << "cache_mode,warm_across_files\n";

    std::cout
        << "cache_layout,linear_probing_dynamic_table\n";

    std::cout
        << "hash,polynomial_lb_rb_offset_symbol\n";

    std::cout
        << "resize_load_factor,0.70\n\n";

    std::cout.flush();
}

inline void print_result_header() {
    std::cout
        << "file,input_bytes,rlz_ms,phrases,"
        << "file_cache_hits,file_cache_misses,"
        << "file_cache_hit_rate,total_cache_entries,"
        << "table_slots,max_probe_cluster,load_factor,"
        << "approx_cache_MB,peak_RSS_MB\n";

    std::cout.flush();
}

inline void print_file_result(
    const FileResult& result
) {
    std::cout
        << result.filename << ','
        << result.input_bytes << ','
        << std::fixed
        << std::setprecision(3)
        << result.rlz_ms << ','
        << result.phrases << ','
        << result.cache_hits << ','
        << result.cache_misses << ','
        << std::fixed
        << std::setprecision(6)
        << result.cache_hit_rate << ','
        << result.cache_entries << ','
        << result.table_slots << ','
        << result.max_probe_cluster << ','
        << std::fixed
        << std::setprecision(6)
        << result.load_factor << ','
        << std::fixed
        << std::setprecision(3)
        << result.cache_mb << ','
        << std::fixed
        << std::setprecision(2)
        << result.peak_rss_mb
        << '\n';

    std::cout.flush();
}

template <typename Parser>
RunResults run_rlz(
    const LoadedData& data,
    Parser& parser
) {
    RunResults results;

    std::size_t previous_hits = 0;
    std::size_t previous_misses = 0;

    for (
        std::size_t i = 0;
        i < data.files.size();
        ++i
    ) {
        const std::string& filename =
            data.files[i];

        std::cerr
            << "Processing "
            << i + 1
            << '/'
            << data.files.size()
            << ": "
            << filename
            << std::endl;

        const auto input =
            read_file<Symbol>(filename.c_str());

        if (input.empty()) {
            ++results.skipped_files;

            std::cerr
                << "Skipped empty input: "
                << filename
                << std::endl;

            continue;
        }

        Triples phrases;

        FileResult result;
        result.filename = filename;
        result.input_bytes = input.size();

        result.rlz_ms = measure_ms([&] {
            phrases = parser.lzFactorize(input);
        });

        result.phrases = phrases.size();

        const auto after =
            parser.cache_info();

        result.cache_hits =
            after.hits - previous_hits;

        result.cache_misses =
            after.misses - previous_misses;

        previous_hits = after.hits;
        previous_misses = after.misses;

        const std::size_t queries =
            result.cache_hits +
            result.cache_misses;

        result.cache_hit_rate =
            queries == 0
                ? 0.0
                : static_cast<double>(
                      result.cache_hits
                  ) /
                  static_cast<double>(
                      queries
                  );

        result.cache_entries =
            after.current_size;

        result.table_slots =
            after.table_slots;

        result.max_probe_cluster =
            after.max_probe_cluster;

        result.load_factor =
            after.load_factor;

        result.cache_mb =
            static_cast<double>(
                after.approx_bytes
            ) /
            (1024.0 * 1024.0);

        result.peak_rss_mb =
            peak_rss_mb();

        results.add(result);
        print_file_result(result);

        std::cerr
            << "Completed "
            << i + 1
            << '/'
            << data.files.size()
            << ": "
            << std::fixed
            << std::setprecision(3)
            << result.rlz_ms
            << " ms"
            << ", phrases="
            << result.phrases
            << ", entries="
            << result.cache_entries
            << ", slots="
            << result.table_slots
            << ", max_cluster="
            << result.max_probe_cluster
            << std::endl;
    }

    return results;
}

template <typename CacheInfo>
void print_lp_summary(
    const RunResults& results,
    const CacheInfo& cache_info
) {
    const double average_ms =
        results.processed_files == 0
            ? 0.0
            : results.total_rlz_ms /
              static_cast<double>(
                  results.processed_files
              );

    std::cout
        << "\nsummary,processed_files,"
        << results.processed_files
        << '\n';

    std::cout
        << "summary,skipped_files,"
        << results.skipped_files
        << '\n';

    std::cout
        << "summary,total_input_bytes,"
        << results.total_input_bytes
        << '\n';

    std::cout
        << "summary,total_phrases,"
        << results.total_phrases
        << '\n';

    std::cout
        << "summary,total_rlz_ms,"
        << std::fixed
        << std::setprecision(3)
        << results.total_rlz_ms
        << '\n';

    std::cout
        << "summary,avg_rlz_ms,"
        << std::fixed
        << std::setprecision(3)
        << average_ms
        << '\n';

    std::cout
        << "summary,total_cache_hits,"
        << cache_info.hits
        << '\n';

    std::cout
        << "summary,total_cache_misses,"
        << cache_info.misses
        << '\n';

    std::cout
        << "summary,final_cache_hit_rate,"
        << std::fixed
        << std::setprecision(6)
        << cache_info.hit_rate
        << '\n';

    std::cout
        << "summary,final_cache_entries,"
        << cache_info.current_size
        << '\n';

    std::cout
        << "summary,final_table_slots,"
        << cache_info.table_slots
        << '\n';

    std::cout
        << "summary,final_max_probe_cluster,"
        << cache_info.max_probe_cluster
        << '\n';

    std::cout
        << "summary,final_load_factor,"
        << std::fixed
        << std::setprecision(6)
        << cache_info.load_factor
        << '\n';

    std::cout
        << "summary,final_cache_MB,"
        << std::fixed
        << std::setprecision(3)
        << static_cast<double>(
               cache_info.approx_bytes
           ) /
           (1024.0 * 1024.0)
        << '\n';

    std::cout
        << "summary,peak_RSS_MB,"
        << std::fixed
        << std::setprecision(2)
        << peak_rss_mb()
        << '\n';

    std::cout.flush();
}