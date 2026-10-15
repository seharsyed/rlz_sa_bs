// The shape of a reference's PT16 table (the sassy builder,
// variants/pt16_build_sassy.hpp): how many of its 16-mers are singletons (one
// occurrence) and how many are ranges, and how large the ranges are. Only the
// table's entries are built; nothing is written.
//
// Build (from PT16mer_RLZ/):
//   g++ -std=c++20 -O3 -march=native pt16_stats.cpp -o pt16_stats
//
// Run:
//   ./pt16_stats --reference REF --suffix-array REF.sa
//
// Prints key=value lines on stdout:
//   reference_bytes
//   entries             distinct 16-mers in the reference (entries of L)
//   singletons, ranges  entries with one occurrence / several
//   positions           suffixes covered by the entries (those of length >= 16)
//   range_positions     of them, those in a range
//   mean_range_size     range_positions / ranges
//   expected_range_size the range size of a random covered position (a
//                       singleton's counts as 0): sum of size^2 / positions
//   max_range_size
//   ranges_at_least_T   ranges of size >= T (T = 1000, 10000, 100000, 1000000)
//   bin_LOW_HIGH        ranges of size LOW..HIGH (powers of two, non-empty bins)

#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "variants/pt16_build_sassy.hpp"

namespace {

struct StatsArgs {
  std::string reference, suffix_array;
};

StatsArgs parse(int argc, char** argv) {
  StatsArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    const auto value = [&]() -> std::string {
      if (i + 1 >= argc) throw std::runtime_error("missing value for " + option);
      return argv[++i];
    };
    if (option == "--reference") args.reference = value();
    else if (option == "--suffix-array") args.suffix_array = value();
    else if (option == "--quiet") {}
    else throw std::runtime_error("unknown argument: " + option);
  }
  if (args.reference.empty() || args.suffix_array.empty()) {
    throw std::runtime_error("--reference and --suffix-array are required");
  }
  return args;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const StatsArgs args = parse(argc, argv);
    const std::vector<unsigned char> reference = load_reference<unsigned char>(args.reference);
    const std::vector<std::uint32_t> suffix_array =
        load_suffix_array<std::uint32_t>(args.suffix_array);
    if (suffix_array.size() != reference.size()) {
      throw std::runtime_error("suffix array and reference sizes differ");
    }

    const RangeSizeStats stats = build_sassy_entries(reference, suffix_array).range_stats;
    std::cerr << stats.format();

    const auto ratio = [](const double part, const std::size_t whole) {
      return whole == 0 ? 0.0 : part / static_cast<double>(whole);
    };
    std::cout << std::fixed << std::setprecision(4)
              << "reference_bytes=" << reference.size() << '\n'
              << "entries=" << stats.entries << '\n'
              << "singletons=" << stats.singletons << '\n'
              << "ranges=" << stats.ranges << '\n'
              << "positions=" << stats.positions << '\n'
              << "range_positions=" << stats.range_positions << '\n'
              << "mean_range_size="
              << ratio(static_cast<double>(stats.range_positions), stats.ranges) << '\n'
              << "expected_range_size=" << ratio(stats.range_size_squares, stats.positions)
              << '\n'
              << "max_range_size=" << (stats.top.empty() ? 1 : stats.top.front().first) << '\n';
    for (std::size_t t = 0; t < RangeSizeStats::large_thresholds.size(); ++t) {
      std::cout << "ranges_at_least_" << RangeSizeStats::large_thresholds[t] << '='
                << stats.large[t] << '\n';
    }
    for (std::size_t b = 0; b < RangeSizeStats::bin_count; ++b) {
      if (stats.size_bins[b] == 0) continue;
      std::cout << "bin_" << (std::size_t{1} << b) + 1 << '_' << (std::size_t{1} << (b + 1))
                << '=' << stats.size_bins[b] << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
