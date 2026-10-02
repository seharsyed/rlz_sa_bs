// Statistics of a reference and an input collection, from lrf-ms
// (lrf_ms/lrf_ms.hpp): the reference's average LCP, and every input's
// average matching-statistics length against the reference.
//
//   MS[i] = the length of the longest prefix of input[i, ...) that occurs in
//           the reference (one value per input position; not the LZ phrases)
//   LCP[r] = the longest common prefix of the suffixes at suffix-array ranks
//           r - 1 and r of the reference (r = 1 .. n - 1)
//
// Timing does not matter here: the inputs are processed by several threads,
// each file on its own.
//
// Build (from PT16mer_RLZ/):
//   g++ -std=c++20 -O3 -march=native -pthread ms_stats.cpp -o ms_stats
//
// Run:
//   ./ms_stats --reference REF --suffix-array REF.sa --filenames LIST
//              [--threads T] [--per-file CSV]
//
// Prints key=value lines on stdout: reference_bytes, lcp_mean (over the n - 1
// adjacent pairs), lcp_max, files, input_bytes, ms_mean (over all positions
// of all inputs), ms_mean_of_files (the mean of the per-file means). With
// --per-file: one CSV row per input (file, input_bytes, ms_mean, ms_median,
// ms_max).

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "lrf_ms/lrf_ms.hpp"
#include "pt16_utils.hpp"  // loaders, time_ms

namespace {

using Symbol = unsigned char;
using SAType = std::uint32_t;

struct MsArgs {
  std::string reference, suffix_array, filenames, per_file;
  unsigned threads = 0;
};

MsArgs parse(int argc, char** argv) {
  MsArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    const auto value = [&]() -> std::string {
      if (i + 1 >= argc) throw std::runtime_error("missing value for " + option);
      return argv[++i];
    };
    if (option == "--reference") args.reference = value();
    else if (option == "--suffix-array") args.suffix_array = value();
    else if (option == "--filenames") args.filenames = value();
    else if (option == "--per-file") args.per_file = value();
    else if (option == "--threads") args.threads = static_cast<unsigned>(std::stoul(value()));
    else if (option == "--quiet") {}
    else throw std::runtime_error("unknown argument: " + option);
  }
  if (args.reference.empty() || args.suffix_array.empty() || args.filenames.empty()) {
    throw std::runtime_error("--reference, --suffix-array and --filenames are required");
  }
  if (args.threads == 0) args.threads = std::max(1u, std::thread::hardware_concurrency());
  return args;
}

struct FileStats {
  std::size_t bytes = 0;
  double sum = 0;  // of the MS lengths
  double mean = 0;
  double median = 0;
  std::uint32_t max = 0;
};

}  // namespace

int main(int argc, char** argv) {
  try {
    const MsArgs args = parse(argc, argv);

    const std::vector<Symbol> reference = load_reference<Symbol>(args.reference);
    const std::vector<SAType> suffix_array = load_suffix_array<SAType>(args.suffix_array);
    if (suffix_array.size() != reference.size()) {
      throw std::runtime_error("suffix array and reference sizes differ");
    }
    const std::vector<std::string> files = load_input_list(args.filenames);
    std::cerr << "reference " << reference.size() << " bytes, " << files.size()
              << " inputs, " << args.threads << " threads\n";

    std::unique_ptr<LRFMS<Symbol, SAType>> lrf;
    const double build_ms =
        time_ms([&] { lrf = std::make_unique<LRFMS<Symbol, SAType>>(reference, suffix_array); });
    std::cerr << "lrf-ms built in " << build_ms / 1000.0 << " s\n";

    // The reference's average LCP over the n - 1 adjacent suffix pairs.
    const auto& lcp = lrf->lcp();
    double lcp_sum = 0;
    std::int64_t lcp_max = 0;
    for (std::size_t r = 1; r < lcp.size(); ++r) {
      lcp_sum += lcp[r];
      lcp_max = std::max<std::int64_t>(lcp_max, lcp[r]);
    }
    const double lcp_mean = lcp.size() > 1 ? lcp_sum / static_cast<double>(lcp.size() - 1) : 0;

    // The inputs' matching statistics, several files at a time.
    std::vector<FileStats> stats(files.size());
    std::atomic<std::size_t> next{0}, done{0};
    const auto worker = [&] {
      for (std::size_t f = next.fetch_add(1); f < files.size(); f = next.fetch_add(1)) {
        const std::vector<Symbol> input = load_input<Symbol>(files[f]);
        const MatchingStatistics ms = lrf->computeMatchingStatistics(input);
        FileStats& s = stats[f];
        s.bytes = input.size();
        std::vector<std::uint32_t> lengths(ms.size());
        for (std::size_t i = 0; i < ms.size(); ++i) {
          lengths[i] = ms[i].second;
          s.sum += ms[i].second;
          s.max = std::max(s.max, ms[i].second);
        }
        s.mean = ms.empty() ? 0 : s.sum / static_cast<double>(ms.size());
        if (!lengths.empty()) {
          auto middle = lengths.begin() + static_cast<std::ptrdiff_t>(lengths.size() / 2);
          std::nth_element(lengths.begin(), middle, lengths.end());
          s.median = *middle;
        }
        const std::size_t k = done.fetch_add(1) + 1;
        if (k % 10 == 0 || k == files.size()) {
          std::cerr << "\r" << k << "/" << files.size() << " inputs" << std::flush;
        }
      }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 1; t < args.threads; ++t) pool.emplace_back(worker);
    worker();
    for (std::thread& thread : pool) thread.join();
    std::cerr << '\n';

    double total_sum = 0, total_positions = 0, mean_of_files = 0;
    for (const FileStats& s : stats) {
      total_sum += s.sum;
      total_positions += static_cast<double>(s.bytes);
      mean_of_files += s.mean;
    }
    mean_of_files /= std::max<std::size_t>(1, stats.size());

    if (!args.per_file.empty()) {
      std::ofstream csv(args.per_file);
      if (!csv) throw std::runtime_error("cannot create " + args.per_file);
      csv << "file,input_bytes,ms_mean,ms_median,ms_max\n" << std::fixed << std::setprecision(4);
      for (std::size_t f = 0; f < files.size(); ++f) {
        csv << std::filesystem::path(files[f]).filename().string() << ',' << stats[f].bytes
            << ',' << stats[f].mean << ',' << stats[f].median << ',' << stats[f].max << '\n';
      }
    }

    std::cout << std::fixed << std::setprecision(4)
              << "reference_bytes=" << reference.size() << '\n'
              << "lcp_mean=" << lcp_mean << '\n'
              << "lcp_max=" << lcp_max << '\n'
              << "files=" << files.size() << '\n'
              << "input_bytes=" << static_cast<std::uint64_t>(total_positions) << '\n'
              << "ms_mean=" << (total_positions > 0 ? total_sum / total_positions : 0) << '\n'
              << "ms_mean_of_files=" << mean_of_files << '\n';
    std::cerr << "average LCP " << lcp_mean << ", average MS "
              << (total_positions > 0 ? total_sum / total_positions : 0) << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
