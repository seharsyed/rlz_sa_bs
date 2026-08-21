#pragma once

#include "parser.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;


// ---------- Arguments ----------

struct Args {
    std::string reference;
    std::string suffix_array;
    std::string filenames;
    std::string pt16_table;
    std::string results;
};

inline std::string require_value(int& i, int argc, char** argv) {
    if (i + 1 >= argc) {
        throw std::runtime_error("missing command-line value");
    }

    return argv[++i];
}

inline void print_usage(const char* program) {
    std::cout << "Usage: " << program
              << " --reference PATH"
              << " --suffix-array PATH"
              << " --filenames PATH"
              << " --results PATH"
              << " [--pt16-table PATH]\n";
}

inline Args parse_args(int argc, char** argv) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];

        if (option == "--reference") {
            args.reference = require_value(i, argc, argv);
        } else if (option == "--suffix-array") {
            args.suffix_array = require_value(i, argc, argv);
        } else if (option == "--filenames") {
            args.filenames = require_value(i, argc, argv);
        } else if (option == "--pt16-table") {
            args.pt16_table = require_value(i, argc, argv);
        } else if (option == "--results") {
            args.results = require_value(i, argc, argv);
        } else if (option == "--help" || option == "-h") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            throw std::runtime_error("unknown argument: " + option);
        }
    }

    if (args.reference.empty() || args.suffix_array.empty() || args.filenames.empty() || args.results.empty()) {
        throw std::runtime_error("reference, suffix-array, filenames and results are required");
    }

    return args;
}


// ---------- Loading ----------

template <typename Symbol>
inline std::vector<Symbol> load_reference(const std::string& path) {
    return read_file<Symbol>(path.c_str());
}

template <typename SAType>
inline std::vector<SAType> load_suffix_array(const std::string& path) {
    return read_file<SAType>(path.c_str());
}

template <typename Symbol>
inline std::vector<Symbol> load_input(const std::string& path) {
    return read_file<Symbol>(path.c_str());
}

inline std::vector<std::string> load_input_list(const std::string& path) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error("cannot open input list: " + path);
    }

    std::vector<std::string> files;
    std::string filename;

    while (std::getline(input, filename)) {
        if (!filename.empty() && filename.back() == '\r') {
            filename.pop_back();
        }

        if (!filename.empty()) {
            files.push_back(filename);
        }
    }

    return files;
}


// ---------- PT16 table path ----------

inline std::string get_pt16_path(const Args& args) {
    if (!args.pt16_table.empty()) {
        return args.pt16_table;
    }

    const fs::path reference_path(args.reference);
    const std::string dataset = reference_path.parent_path().filename().string();

    return (reference_path.parent_path() / (dataset + "_pt16_hl.bin")).string();
}


// ---------- Timing ----------

template <typename Fn>
inline double time_ms(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

inline double peak_rss_mb() {
    struct rusage usage{};

    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0.0;
    }

    return static_cast<double>(usage.ru_maxrss) / 1024.0;
}


// ---------- Baseline results ----------

struct BaselineResult {
    std::string filename;
    std::size_t input_bytes;
    double baseline_ms;
    std::size_t baseline_phrases;
    std::string factor_file;
};


// ---------- Baseline factor storage ----------

inline void write_factor_file(const std::string& path, const Triples& factors) {
    std::ofstream output(path, std::ios::binary);

    if (!output) {
        throw std::runtime_error("cannot create baseline factor file");
    }

    const std::uint64_t count = factors.size();
    output.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (const auto& factor : factors) {
        const std::uint64_t input_position = std::get<0>(factor);
        const std::uint64_t reference_position = std::get<1>(factor);
        const std::uint64_t match_length = std::get<2>(factor);

        output.write(reinterpret_cast<const char*>(&input_position), sizeof(input_position));
        output.write(reinterpret_cast<const char*>(&reference_position), sizeof(reference_position));
        output.write(reinterpret_cast<const char*>(&match_length), sizeof(match_length));
    }
}

static bool factor_file_equals(
    const std::string& path,
    const Triples& factors
) {
    std::ifstream input(path, std::ios::binary);

    if (!input) {
        throw std::runtime_error("cannot open temporary factor file: " + path);
    }

    std::uint64_t count = 0;
    input.read(reinterpret_cast<char*>(&count), sizeof(count));

    if (!input) {
        throw std::runtime_error("failed while reading factor count: " + path);
    }

    if (count != static_cast<std::uint64_t>(factors.size())) {
        std::cerr
            << "Factor count mismatch: baseline="
            << count
            << ", PT16="
            << factors.size()
            << '\n';

        return false;
    }

    return true;
}
// ---------- PT16 per-file statistics ----------

struct PT16Delta {
    std::size_t hits;
    std::size_t misses;
    std::size_t singleton_hits;
    std::size_t range_hits;
};


template <typename Stats>
inline PT16Delta stats_difference(const Stats& previous, const Stats& current) {
    return {
        current.hits - previous.hits,
        current.misses - previous.misses,
        current.singleton_hits - previous.singleton_hits,
        current.range_hits - previous.range_hits
    };
}


// ---------- Final summary ----------

struct BenchmarkSummary {
    std::size_t processed_files;
    std::size_t total_input_bytes;
    std::size_t total_baseline_phrases;
    std::size_t total_pt16_phrases;
    double total_baseline_ms;
    double total_pt16_ms;
    double pt16_build_ms;
    bool all_equal;
};


// ---------- CSV ----------

class CSVWriter {
public:
    explicit CSVWriter(const std::string& path) : output(path) {
        if (!output) {
            throw std::runtime_error("cannot create results file: " + path);
        }

        //output << "file,input_bytes,baseline_ms,pt16_ms,speedup,baseline_phrases,pt16_phrases,outputs_equal,pt16_hits,pt16_misses,pt16_hit_rate,pt16_entries,pt16_MB,peak_RSS_MB\n";
    
        output << "file,input_bytes,baseline_ms,pt16_ms,speedup,baseline_phrases,pt16_phrases,outputs_equal,pt16_hits,pt16_misses,pt16_hit_rate,singleton_hits,range_hits,singleton_hit_rate,pt16_entries,pt16_MB,peak_RSS_MB\n";
    
    }

    void write_row(const BaselineResult& baseline, double pt16_ms, std::size_t pt16_phrases, bool equal, const PT16Delta& stats, std::size_t entries, std::size_t approx_bytes) {
        const std::size_t queries = stats.hits + stats.misses;
        const double hit_rate = queries == 0 ? 0.0 : static_cast<double>(stats.hits) / static_cast<double>(queries);
        const double speedup = pt16_ms == 0.0 ? 0.0 : baseline.baseline_ms / pt16_ms;
        const std::size_t hit_types = stats.singleton_hits + stats.range_hits;
        const double singleton_hit_rate = hit_types == 0 ? 0.0 : static_cast<double>(stats.singleton_hits) / static_cast<double>(hit_types);


        output << baseline.filename << ','
               << baseline.input_bytes << ','
               << std::fixed << std::setprecision(3)
               << baseline.baseline_ms << ','
               << pt16_ms << ','
               << speedup << ','
               << baseline.baseline_phrases << ','
               << pt16_phrases << ','
               << (equal ? "YES" : "NO") << ','
               << stats.hits << ','
               << stats.misses << ','
               << std::setprecision(6)
               << hit_rate << ','
               << stats.singleton_hits << ','
                << stats.range_hits << ','
                << singleton_hit_rate << ','
               << entries << ','
               << std::setprecision(3)
               << static_cast<double>(approx_bytes) / (1024.0 * 1024.0) << ','
               << std::setprecision(2)
               << peak_rss_mb()
               << '\n';
    }

    template <typename Stats>
    void write_summary(const BenchmarkSummary& summary, const Stats& stats) {
        const double speedup = summary.total_pt16_ms == 0.0 ? 0.0 : summary.total_baseline_ms / summary.total_pt16_ms;
        const std::size_t hit_types = stats.singleton_hits + stats.range_hits;
        const double singleton_hit_rate = hit_types == 0 ? 0.0 : static_cast<double>(stats.singleton_hits) / static_cast<double>(hit_types);
        
        output << "summary,processed_files," << summary.processed_files << '\n';
        output << "summary,total_input_bytes," << summary.total_input_bytes << '\n';
        output << "summary,total_baseline_phrases," << summary.total_baseline_phrases << '\n';
        output << "summary,total_pt16_phrases," << summary.total_pt16_phrases << '\n';
        output << "summary,total_baseline_ms," << std::fixed << std::setprecision(3) << summary.total_baseline_ms << '\n';
        output << "summary,total_pt16_ms," << summary.total_pt16_ms << '\n';
        output << "summary,overall_speedup," << speedup << '\n';
        output << "summary,all_outputs_equal," << (summary.all_equal ? "YES" : "NO") << '\n';
        output << "summary,total_pt16_hits," << stats.hits << '\n';
        output << "summary,total_pt16_misses," << stats.misses << '\n';
        output << "summary,pt16_entries," << stats.entries << '\n';
        output << "summary,pt16_MB," << std::setprecision(3) << static_cast<double>(stats.approx_bytes) / (1024.0 * 1024.0) << '\n';
        output << "summary,peak_RSS_MB," << std::setprecision(2) << peak_rss_mb() << '\n';
        output << "summary,total_singleton_hits," << stats.singleton_hits << '\n';
        output << "summary,total_range_hits," << stats.range_hits << '\n';
        output << "summary,singleton_hit_rate," << std::setprecision(6) << singleton_hit_rate << '\n';
        output << "summary,pt16_build_ms," << summary.pt16_build_ms << '\n';
    }

private:
    std::ofstream output;
};