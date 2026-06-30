#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "parser.hpp"
#include "rlz_lp.hpp"

using Symbol = unsigned char;
using SAType = unsigned int;
using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;

struct Args {
    std::string reference;
    std::string suffix_array;
    std::string filenames;
    std::size_t div_p = 64;
    std::string mode = "both";
};

template <typename Fn>
double time_ms(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

std::size_t peak_rss_kb() {
#if defined(__unix__) || defined(__APPLE__)
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<std::size_t>(usage.ru_maxrss / 1024);
#else
        return static_cast<std::size_t>(usage.ru_maxrss);
#endif
    }
#endif
    return 0;
}

void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " --reference PATH --filenames PATH [--suffix-array PATH] "
        << "[--div_p N] [--mode both|baseline|cached]\n\n"
        << "If --suffix-array is omitted, the program uses <reference>.sa\n"
        << "Default --div_p is 64.\n";
}

std::string trim(std::string s) {
    while (!s.empty() &&
           (s.back() == '\n' || s.back() == '\r' ||
            s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }

    std::size_t first = 0;
    while (first < s.size() && (s[first] == ' ' || s[first] == '\t')) {
        ++first;
    }

    return s.substr(first);
}

std::vector<std::string> read_filename_list(const std::string& path) {
    std::ifstream in(path);

    if (!in) {
        throw std::runtime_error("cannot open filenames file: " + path);
    }

    std::vector<std::string> files;
    std::string line;

    while (std::getline(in, line)) {
        line = trim(line);

        if (line.empty() || line[0] == '#') {
            continue;
        }

        files.push_back(line);
    }

    if (files.empty()) {
        throw std::runtime_error("filenames file is empty: " + path);
    }

    return files;
}

Args parse_args(int argc, char** argv) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];

        auto require_value = [&](const std::string& option) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("missing value for " + option);
            }
            return argv[++i];
        };

        if (key == "--reference") {
            args.reference = require_value(key);
        } else if (key == "--suffix-array") {
            args.suffix_array = require_value(key);
        } else if (key == "--filenames") {
            args.filenames = require_value(key);
        } else if (key == "--div_p") {
            args.div_p =
                static_cast<std::size_t>(std::stoull(require_value(key)));
        } else if (key == "--mode") {
            args.mode = require_value(key);
        } else if (key == "--help" || key == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + key);
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
        throw std::runtime_error("--div_p must be > 0");
    }

    if (args.mode != "both" && args.mode != "baseline" && args.mode != "cached") {
        throw std::runtime_error("--mode must be one of: both, baseline, cached");
    }

    return args;
}

int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);
        const auto files = read_filename_list(args.filenames);

        std::cerr << "Loading reference and suffix array...\n";

        auto ref = read_file<Symbol>(args.reference.c_str());
        std::cerr << "Reference loaded: " << ref.size() << " bytes\n";

        auto sa = read_file<SAType>(args.suffix_array.c_str());
        std::cerr << "Suffix array loaded: " << sa.size() << " entries\n";

        if (ref.empty() || sa.empty()) {
            throw std::runtime_error("reference and suffix array must be non-empty");
        }

        if (sa.size() != ref.size()) {
            std::cerr << "Warning: SA entries (" << sa.size()
                      << ") != reference bytes (" << ref.size() << ")\n";
        }

        std::unique_ptr<RLZLPParser<Symbol, SAType>> parser;

        if (args.mode != "baseline") {
            parser = std::make_unique<RLZLPParser<Symbol, SAType>>(
                ref, sa, args.div_p
            );
            std::cerr << "LP parser created\n";
        }

        std::cout << "reference," << args.reference << "\n";
        std::cout << "suffix_array," << args.suffix_array << "\n";
        std::cout << "files," << files.size() << "\n";
        std::cout << "div_p," << args.div_p << "\n";
        std::cout << "run_mode," << args.mode << "\n";
        std::cout << "cache_mode,"
                  << (args.mode == "baseline" ? "none" : "warm_across_files")
                  << "\n";
        std::cout << "cache_layout,linear_probing_dynamic_table\n";
        std::cout << "\n";

        std::cout
            << "file,input_bytes,baseline_ms,cached_ms,baseline_phrases,cached_phrases,outputs_equal,"
            << "file_cache_hits,file_cache_misses,file_cache_hit_rate,total_cache_entries,approx_cache_MB,peak_RSS_MB\n";

        std::cout.flush();

        bool all_equal = true;
        std::size_t processed = 0;
        std::size_t total_input_bytes = 0;
        std::size_t total_baseline_phrases = 0;
        std::size_t total_cached_phrases = 0;
        double total_baseline_ms = 0.0;
        double total_cached_ms = 0.0;

        for (const auto& file : files) {
            std::cerr << "Processing " << (processed + 1) << "/" << files.size()
                      << ": " << file << "\n";

            auto input = read_file<Symbol>(file.c_str());

            if (input.empty()) {
                std::cerr << "Warning: empty input skipped: " << file << "\n";
                continue;
            }

            Triples baseline;
            Triples cached;
            double baseline_ms = 0.0;
            double cached_ms = 0.0;
            std::size_t baseline_phrases = 0;
            std::size_t cached_phrases = 0;
            std::string equal_text = "NA";
            std::size_t file_hits = 0;
            std::size_t file_misses = 0;
            double file_hit_rate = 0.0;
            std::size_t cache_entries = 0;
            double cache_mb = 0.0;

            if (args.mode != "cached") {
                baseline_ms = time_ms([&] {
                    baseline = lzFactorize<Symbol, SAType>(input, ref, sa);
                });

                baseline_phrases = baseline.size();
            }

            if (args.mode != "baseline") {
                const auto before = parser->cache_info();

                cached_ms = time_ms([&] {
                    cached = parser->lzFactorize(input);
                });

                const auto after = parser->cache_info();

                file_hits = after.hits - before.hits;
                file_misses = after.misses - before.misses;

                const std::size_t file_queries = file_hits + file_misses;
                file_hit_rate = file_queries == 0
                    ? 0.0
                    : static_cast<double>(file_hits) / static_cast<double>(file_queries);

                cached_phrases = cached.size();
                cache_entries = after.current_size;
                cache_mb = after.approx_bytes / (1024.0 * 1024.0);
            }

            if (args.mode == "both") {
                const bool equal = (baseline == cached);
                all_equal = all_equal && equal;
                equal_text = equal ? "YES" : "NO";
            }

            ++processed;
            total_input_bytes += input.size();
            total_baseline_phrases += baseline_phrases;
            total_cached_phrases += cached_phrases;
            total_baseline_ms += baseline_ms;
            total_cached_ms += cached_ms;

            std::cout
                << file << ','
                << input.size() << ','
                << std::fixed << std::setprecision(3) << baseline_ms << ','
                << std::fixed << std::setprecision(3) << cached_ms << ','
                << baseline_phrases << ','
                << cached_phrases << ','
                << equal_text << ','
                << file_hits << ','
                << file_misses << ','
                << std::fixed << std::setprecision(6) << file_hit_rate << ','
                << cache_entries << ','
                << std::fixed << std::setprecision(3) << cache_mb << ','
                << std::fixed << std::setprecision(2)
                << peak_rss_kb() / 1024.0
                << '\n';

            std::cout.flush();
        }

        RLZLPParser<Symbol, SAType>::CacheInfo final_info{};

        if (parser) {
            final_info = parser->cache_info();
        }

        const double avg_baseline_ms = processed == 0
            ? 0.0
            : total_baseline_ms / static_cast<double>(processed);

        const double avg_cached_ms = processed == 0
            ? 0.0
            : total_cached_ms / static_cast<double>(processed);

        const double speedup = total_cached_ms == 0.0
            ? 0.0
            : total_baseline_ms / total_cached_ms;

        std::cout << "\nsummary,processed_files," << processed << "\n";
        std::cout << "summary,total_input_bytes," << total_input_bytes << "\n";
        std::cout << "summary,total_baseline_phrases," << total_baseline_phrases << "\n";
        std::cout << "summary,total_cached_phrases," << total_cached_phrases << "\n";
        std::cout << "summary,total_baseline_ms,"
                  << std::fixed << std::setprecision(3) << total_baseline_ms << "\n";
        std::cout << "summary,total_cached_ms,"
                  << std::fixed << std::setprecision(3) << total_cached_ms << "\n";
        std::cout << "summary,avg_baseline_ms,"
                  << std::fixed << std::setprecision(3) << avg_baseline_ms << "\n";
        std::cout << "summary,avg_cached_ms,"
                  << std::fixed << std::setprecision(3) << avg_cached_ms << "\n";
        std::cout << "summary,overall_speedup,"
                  << std::fixed << std::setprecision(3) << speedup << "\n";
        std::cout << "summary,all_outputs_equal,"
                  << (all_equal ? "YES" : "NO") << "\n";
        std::cout << "summary,total_cache_hits," << final_info.hits << "\n";
        std::cout << "summary,total_cache_misses," << final_info.misses << "\n";
        std::cout << "summary,final_cache_hit_rate,"
                  << std::fixed << std::setprecision(6)
                  << final_info.hit_rate << "\n";
        std::cout << "summary,final_cache_entries,"
                  << final_info.current_size << "\n";
        std::cout << "summary,final_table_slots,"
                  << final_info.table_slots << "\n";
        std::cout << "summary,final_load_factor,"
                  << std::fixed << std::setprecision(6)
                  << final_info.load_factor << "\n";
        std::cout << "summary,final_cache_MB,"
                  << std::fixed << std::setprecision(3)
                  << final_info.approx_bytes / (1024.0 * 1024.0) << "\n";
        std::cout << "summary,peak_RSS_MB,"
                  << std::fixed << std::setprecision(2)
                  << peak_rss_kb() / 1024.0 << "\n";

        std::cout.flush();

        return all_equal ? 0 : 4;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        print_usage(argv[0]);
        return 1;
    }
}