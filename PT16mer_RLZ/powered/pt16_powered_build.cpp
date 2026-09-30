// Builds the powered PT16 table (pt16_powered.hpp: keyed by the 16-mer read
// from its right end, for powered's backward search) from a reference and its
// powered_rlz index -- the same REF_four.bwt (+ REF_four_data.bwt) the
// parsing programs take. No suffix array: the rows and their rotation
// starts come from the index itself.
//
// Build (x86-64, GCC; block sizes as for powered_rlz's tools), from
// PT16mer_RLZ/:
//
//   g++ -std=c++2a -O3 -march=native -DNDEBUG \
//       -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 \
//       powered/pt16_powered_build.cpp -o pt16_powered_build
//
// Run:
//
//   ./pt16_powered_build --reference REF --index REF_four.bwt [--output PATH]
//
// The default output is the index path with .bwt replaced by .pt16
// (REF_four.pt16).

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "../../powered_rlz/include/types.hpp"  // bbwt::non_rle
#include "pt16_powered.hpp"

namespace {

struct BuildArgs {
  std::string reference, index, output;
};

BuildArgs parse_build_args(int argc, char** argv) {
  BuildArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--reference") {
      args.reference = require_value(i, argc, argv);
    } else if (option == "--index") {
      args.index = require_value(i, argc, argv);
    } else if (option == "--output") {
      args.output = require_value(i, argc, argv);
    } else if (option == "--help" || option == "-h") {
      std::cout << "Usage: " << argv[0]
                << " --reference REF --index REF_four.bwt [--output PATH]\n";
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }
  if (args.reference.empty() || args.index.empty()) {
    throw std::runtime_error("--reference and --index are required");
  }
  if (args.output.empty()) {
    args.output = fs::path(args.index).replace_extension(".pt16").string();
  }
  return args;
}

// Checks the index file and its _data file exist (the index loader exits
// without saying which file it could not open).
void require_index_files(const std::string& index) {
  const fs::path path(index);
  const fs::path data = path.parent_path() /
                        (path.stem().string() + "_data" + path.extension().string());
  for (const fs::path& p : {path, data}) {
    if (!fs::exists(p)) {
      throw std::runtime_error("powered index file not found: " + p.string());
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const BuildArgs args = parse_build_args(argc, argv);
    require_index_files(args.index);

    std::vector<unsigned char> reference;
    const double reference_ms = time_ms(
        [&] { reference = load_reference<unsigned char>(args.reference); });
    std::cerr << "[1] reference: " << reference.size() << " characters ("
              << std::fixed << std::setprecision(1) << reference_ms
              << " ms)\n";

    std::unique_ptr<bbwt::non_rle<>> index;
    const double index_ms = time_ms(
        [&] { index = std::make_unique<bbwt::non_rle<>>(args.index); });
    std::cerr << "[2] powered index: " << index->size() << " rows ("
              << index_ms << " ms)\n";

    PT16PoweredTable table;
    PoweredTableStats stats;
    const double build_ms = time_ms([&] {
      table = PT16PoweredTable::build(reference, index->gca_, &stats);
    });

    const double write_ms = time_ms([&] { table.write(args.output); });

    const auto percent = [&](std::uint64_t x) {
      return stats.entries == 0 ? 0.0 : 100.0 * x / stats.entries;
    };
    std::cerr << "[3] table built in " << build_ms << " ms, written in "
              << write_ms << " ms: " << args.output << '\n'
              << "    entries " << stats.entries << " (singletons "
              << stats.singletons << ", " << percent(stats.singletons)
              << "%; ranges " << stats.ranges << ", " << percent(stats.ranges)
              << "%), largest range " << stats.largest_range << " rows\n"
              << "    empty buckets " << stats.empty_buckets << " of "
              << NUMBER_OF_BUCKETS << '\n'
              << "    rows whose 16-mer wraps past the end: "
              << stats.wrapping_rows << '\n'
              << "    size " << table.bytes() / (1024.0 * 1024.0) << " MB\n";

    std::cout << "reference_chars=" << reference.size() << '\n'
              << "pt16_powered_entries=" << stats.entries << '\n'
              << "pt16_powered_singletons=" << stats.singletons << '\n'
              << "pt16_powered_ranges=" << stats.ranges << '\n'
              << "pt16_powered_build_ms=" << build_ms << '\n'
              << "pt16_powered_bytes=" << table.bytes() << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
