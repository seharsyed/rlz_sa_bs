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
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "parser.hpp"

using Symbol = unsigned char;
using SAType = unsigned int;
using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;

struct Args {
    std::string reference;
    std::string suffix_array;
    std::string filenames;
    std::size_t div_p = 64;
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
    double cache_mb = 0.0;
    double peak_rss_mb = 0.0;
};

struct RunResults {
    std::vector<FileResult> files;
    std::size_t skipped_files = 0;
    std::size_t total_input_bytes = 0;
    std::size_t total_phrases = 0;
    double total_rlz_ms = 0.0;

    void add(FileResult result) {
        total_input_bytes += result.input_bytes;
        total_phrases += result.phrases;
        total_rlz_ms += result.rlz_ms;
        files.push_back(std::move(result));
    }
};

inline void print_usage(const char* program) {
    std::cerr << "Usage:\n"
              << "  " << program
              << " --reference PATH --filenames PATH [--suffix-array PATH] [--div_p N]\n\n"
              << "If --suffix-array is omitted, <reference>.sa is used.\n"
              << "Default --div_p is 64.\n";
}

inline Args parse_args(int argc, char** argv) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];

        auto require_value = [&](const std::string& name) {
            if (i + 1 >= argc) {
                throw std::runtime_error("missing value for " + name);
            }
            return std::string(argv[++i]);
        };

        if (option == "--reference") {
            args.reference = require_value(option);
        } else if (option == "--suffix-array") {
            args.suffix_array = require_value(option);
        } else if (option == "--filenames") {
            args.filenames = require_value(option);
        } else if (option == "--div_p") {
            args.div_p = static_cast<std::size_t>(std::stoull(require_value(option)));
        } else if (option == "--help" || option == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + option);
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

    if (args.div_p == 0) {
        throw std::runtime_error("--div_p must be greater than zero");
    }

    return args;
}

inline std::string trim(std::string value) {
    while (!value.empty() &&
           (value.back() == '\n' || value.back() == '\r' ||
            value.back() == ' ' || value.back() == '\t')) {
        value.pop_back();
    }

    std::size_t first = 0;

    while (first < value.size() &&
           (value[first] == ' ' || value[first] == '\t')) {
        ++first;
    }

    return value.substr(first);
}

inline std::vector<std::string> read_filename_list(const std::string& path) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error("cannot open filenames file: " + path);
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
        throw std::runtime_error("filenames file is empty: " + path);
    }

    return files;
}

inline LoadedData load_data(const Args& args) {
    LoadedData data;

    data.files = read_filename_list(args.filenames);
    data.ref = read_file<Symbol>(args.reference.c_str());
    data.sa = read_file<SAType>(args.suffix_array.c_str());

    if (data.ref.empty() || data.sa.empty()) {
        throw std::runtime_error("reference and suffix array must be non-empty");
    }

    if (data.sa.size() != data.ref.size()) {
        std::cerr << "Warning: SA entries (" << data.sa.size()
                  << ") != reference bytes (" << data.ref.size() << ")\n";
    }

    return data;
}

template <typename Function>
double measure_ms(Function&& function) {
    const auto start = std::chrono::steady_clock::now();
    function();
    const auto end = std::chrono::steady_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

inline double peak_rss_mb() {
#if defined(__unix__) || defined(__APPLE__)
    struct rusage usage{};

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
        return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
    }
#endif

    return 0.0;
}

template <typename Parser>
RunResults run_rlz(const LoadedData& data, Parser& parser) {
    RunResults results;

    for (const auto& filename : data.files) {
        const auto input = read_file<Symbol>(filename.c_str());

        if (input.empty()) {
            ++results.skipped_files;
            continue;
        }

        const auto before = parser.cache_info();
        Triples phrases;

        FileResult result;
        result.filename = filename;
        result.input_bytes = input.size();
        result.rlz_ms = measure_ms([&] { phrases = parser.lzFactorize(input); });
        result.phrases = phrases.size();

        const auto after = parser.cache_info();

        result.cache_hits = after.hits - before.hits;
        result.cache_misses = after.misses - before.misses;

        const std::size_t queries = result.cache_hits + result.cache_misses;

        result.cache_hit_rate = queries == 0
            ? 0.0
            : static_cast<double>(result.cache_hits) /
              static_cast<double>(queries);

        result.cache_entries = after.current_size;
        result.cache_mb = static_cast<double>(after.approx_bytes) / (1024.0 * 1024.0);
        result.peak_rss_mb = peak_rss_mb();

        results.add(std::move(result));
    }

    return results;
}

inline void print_results(
    const Args& args,
    const LoadedData& data,
    const RunResults& results
) {
    std::cout << "reference," << args.reference << '\n';
    std::cout << "suffix_array," << args.suffix_array << '\n';
    std::cout << "files_requested," << data.files.size() << '\n';
    std::cout << "div_p," << args.div_p << '\n';
    std::cout << "cache_mode,warm_across_files\n";
    std::cout << "cache_layout,linear_probing_fixed_table\n\n";

    std::cout << "file,input_bytes,rlz_ms,phrases,file_cache_hits,"
              << "file_cache_misses,file_cache_hit_rate,total_cache_entries,"
              << "approx_cache_MB,peak_RSS_MB\n";

    for (const auto& result : results.files) {
        std::cout << result.filename << ','
                  << result.input_bytes << ','
                  << std::fixed << std::setprecision(3) << result.rlz_ms << ','
                  << result.phrases << ','
                  << result.cache_hits << ','
                  << result.cache_misses << ','
                  << std::fixed << std::setprecision(6) << result.cache_hit_rate << ','
                  << result.cache_entries << ','
                  << std::fixed << std::setprecision(3) << result.cache_mb << ','
                  << std::fixed << std::setprecision(2) << result.peak_rss_mb
                  << '\n';
    }
}

template <typename CacheInfo>
void print_lp_summary(const RunResults& results, const CacheInfo& cache_info) {
    const double average_ms = results.files.empty()
        ? 0.0
        : results.total_rlz_ms / static_cast<double>(results.files.size());

    std::cout << "\nsummary,processed_files," << results.files.size() << '\n';
    std::cout << "summary,skipped_files," << results.skipped_files << '\n';
    std::cout << "summary,total_input_bytes," << results.total_input_bytes << '\n';
    std::cout << "summary,total_phrases," << results.total_phrases << '\n';

    std::cout << "summary,total_rlz_ms,"
              << std::fixed << std::setprecision(3)
              << results.total_rlz_ms << '\n';

    std::cout << "summary,avg_rlz_ms,"
              << std::fixed << std::setprecision(3)
              << average_ms << '\n';

    std::cout << "summary,total_cache_hits," << cache_info.hits << '\n';
    std::cout << "summary,total_cache_misses," << cache_info.misses << '\n';

    std::cout << "summary,final_cache_hit_rate,"
              << std::fixed << std::setprecision(6)
              << cache_info.hit_rate << '\n';

    std::cout << "summary,final_cache_entries,"
              << cache_info.current_size << '\n';

    std::cout << "summary,final_table_slots,"
              << cache_info.table_slots << '\n';

    std::cout << "summary,final_load_factor,"
              << std::fixed << std::setprecision(6)
              << cache_info.load_factor << '\n';

    std::cout << "summary,final_cache_MB,"
              << std::fixed << std::setprecision(3)
              << static_cast<double>(cache_info.approx_bytes) /
                 (1024.0 * 1024.0)
              << '\n';

    std::cout << "summary,table_full,"
              << (cache_info.table_full ? "YES" : "NO")
              << '\n';

    std::cout << "summary,peak_RSS_MB,"
              << std::fixed << std::setprecision(2)
              << peak_rss_mb()
              << '\n';
}