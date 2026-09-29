// RLZ parsing, parallel over input files: one parser at a time, T threads each
// taking the next input file and parsing it against the shared, read-only
// table or index.
//
// Parsers (as in rlz_suite.cpp): sa-binary-search, lrf-ms, pt16, pt16-v2,
// sassy, and powered (built with -DWITH_POWERED, given --powered-index).
//
// Per parser: build or load its structure once (timed apart), then parse all
// files with T threads, then free the structure before the next parser, so
// only one parser's structures are in memory at a time.
//
// Timed, per file: the parse call only (input in memory -> phrase list in
// memory); the phrase list is freed after the timer stops. Reported per
// parser: the wall time of the whole parallel run, and the sum of the
// per-file parse times (the work; sum / wall is the effective parallelism).
// Inputs are loaded before the parsers run (--preload, the default), so the
// wall time is parsing only; with --no-preload each thread loads its own
// files, which then counts in the wall time (not in the per-file times).
//
// No correctness checks here (rlz_suite does those); phrase totals are
// reported.
//
// Build (from PT16mer_RLZ/):
//   g++ -std=c++20 -O3 -pthread rlz_parallel.cpp -o rlz_parallel
//   with powered (x86-64, GCC):
//   g++ -std=c++2a -O3 -march=native -pthread -DNDEBUG -DWITH_POWERED \
//       -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 rlz_parallel.cpp -o rlz_parallel
//
// Run:
//   ./rlz_parallel --reference REF --suffix-array REF.sa --filenames LIST
//                  [--threads T] [--parsers a,b,...] [--powered-index REF_four.bwt]
//                  [--table PATH] [--max-files N] [--no-preload] [--results CSV]

// The PT16 parsers' lookup counters would be written by every thread: turn
// them off (see pt16_counter in pt16_utils.hpp).
#define PT16_NO_STATS

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "lrf_ms/lrf_ms.hpp"
#include "parser.hpp"
#include "pt16_utils.hpp"
#include "rlz_common.hpp"

namespace pt16_v1 {
#include "pt16_build.hpp"
#include "pt16_rlz.hpp"
}  // namespace pt16_v1

#include "pt16_build_v2.hpp"
#include "pt16_rlz_v2.hpp"
#include "variants/pt16_build_sassy.hpp"
#include "variants/pt16_sassy.hpp"

#ifdef WITH_POWERED
#include "../powered_rlz/include/types.hpp"
#endif

using Symbol = unsigned char;
using SAType = std::uint32_t;
using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;
namespace fs = std::filesystem;

namespace {

// ---------- Arguments ----------

const std::vector<std::string> kAllParsers = {
    "sa-binary-search", "lrf-ms", "pt16", "pt16-v2", "sassy", "powered"};

struct ParallelArgs {
  std::string reference, suffix_array, filenames, table, powered_index, results;
  std::size_t max_files = 0;
  unsigned threads = 0;  // 0: all hardware threads
  bool preload = true;
  std::set<std::string> parsers;  // empty: all
};

std::string require_value(int& i, int argc, char** argv) {
  if (i + 1 >= argc) {
    throw std::runtime_error(std::string("missing value for ") + argv[i]);
  }
  return argv[++i];
}

void print_usage(const char* program) {
  std::cout << "Usage: " << program
            << " --reference PATH --suffix-array PATH --filenames PATH\n"
               "  [--threads T]           threads (default: all hardware threads)\n"
               "  [--parsers a,b,...]     subset of: sa-binary-search, lrf-ms, "
               "pt16, pt16-v2, sassy, powered\n"
               "  [--powered-index PATH]  powered_rlz index (REF_four.bwt)\n"
               "  [--table PATH]          PT16 table files (default: "
               "<reference>.parallel_pt16, .v2, .sassy)\n"
               "  [--max-files N]         only the first N input files\n"
               "  [--no-preload]          threads load their own files (the "
               "loading then counts in the wall time)\n"
               "  [--results PATH]        CSV, one row per parser\n";
}

ParallelArgs parse_parallel_args(int argc, char** argv) {
  ParallelArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--reference") {
      args.reference = require_value(i, argc, argv);
    } else if (option == "--suffix-array") {
      args.suffix_array = require_value(i, argc, argv);
    } else if (option == "--filenames") {
      args.filenames = require_value(i, argc, argv);
    } else if (option == "--table") {
      args.table = require_value(i, argc, argv);
    } else if (option == "--powered-index") {
      args.powered_index = require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = std::stoull(require_value(i, argc, argv));
    } else if (option == "--threads") {
      args.threads = static_cast<unsigned>(std::stoul(require_value(i, argc, argv)));
    } else if (option == "--no-preload") {
      args.preload = false;
    } else if (option == "--parsers") {
      std::stringstream list(require_value(i, argc, argv));
      std::string name;
      while (std::getline(list, name, ',')) {
        if (std::find(kAllParsers.begin(), kAllParsers.end(), name) ==
            kAllParsers.end()) {
          throw std::runtime_error("unknown parser: " + name);
        }
        args.parsers.insert(name);
      }
    } else if (option == "--help" || option == "-h") {
      print_usage(argv[0]);
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }
  if (args.reference.empty() || args.suffix_array.empty() ||
      args.filenames.empty()) {
    throw std::runtime_error("reference, suffix-array and filenames are required");
  }
  if (args.table.empty()) args.table = args.reference + ".parallel_pt16";
  if (args.threads == 0) {
    args.threads = std::max(1u, std::thread::hardware_concurrency());
  }
#ifndef WITH_POWERED
  if (!args.powered_index.empty()) {
    throw std::runtime_error(
        "--powered-index given, but this build has no powered support");
  }
#endif
  return args;
}

double mb(const std::size_t bytes) {
  return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

// ---------- The inputs ----------

struct Inputs {
  std::vector<std::string> files;
  std::vector<std::vector<Symbol>> loaded;  // empty when not preloaded
  std::size_t total_bytes = 0;

  // The input of file f: the preloaded copy, or loaded now into `buffer`.
  const std::vector<Symbol>& get(std::size_t f,
                                 std::vector<Symbol>& buffer) const {
    if (!loaded.empty()) return loaded[f];
    buffer = load_input<Symbol>(files[f]);
    return buffer;
  }
};

// ---------- One parser, run in parallel ----------

struct Result {
  std::string name;
  double build_ms = 0.0, load_ms = 0.0;
  std::size_t own_bytes = 0, needed_bytes = 0;
  double wall_ms = 0.0;      // the whole parallel run
  double work_ms = 0.0;      // sum of the per-file parse times
  double longest_ms = 0.0;   // the slowest single file
  std::size_t phrases = 0;
};

// Runs `parse` (input -> phrase list, returned by value) on every file with
// `threads` threads; each thread takes the next unparsed file.
template <typename Parse>
void run_parallel(const Inputs& inputs, const unsigned threads, Parse parse,
                  Result& result) {
  const std::size_t n = inputs.files.size();
  std::vector<double> times(n, 0.0);
  std::vector<std::size_t> phrases(n, 0);
  std::atomic<std::size_t> next{0};

  const auto worker = [&] {
    std::vector<Symbol> buffer;
    for (std::size_t f = next.fetch_add(1); f < n; f = next.fetch_add(1)) {
      const std::vector<Symbol>& input = inputs.get(f, buffer);
      const auto start = std::chrono::steady_clock::now();
      auto list = parse(input);  // timed: the parse only
      times[f] = std::chrono::duration<double, std::milli>(
                     std::chrono::steady_clock::now() - start)
                     .count();
      phrases[f] = list.size();
    }  // the phrase list is freed here, untimed
  };

  result.wall_ms = time_ms([&] {
    std::vector<std::thread> pool;
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(worker);
    worker();  // the main thread works too
    for (std::thread& thread : pool) thread.join();
  });

  for (std::size_t f = 0; f < n; ++f) {
    result.work_ms += times[f];
    result.longest_ms = std::max(result.longest_ms, times[f]);
    result.phrases += phrases[f];
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const ParallelArgs args = parse_parallel_args(argc, argv);
    const auto wanted = [&](const std::string& name) {
      return args.parsers.empty() || args.parsers.count(name) != 0;
    };

    // What the machine offers (std::thread::hardware_concurrency: logical
    // cores, 0 if unknown), next to what this run uses.
    const unsigned available = std::thread::hardware_concurrency();

    std::cerr << "========================================\n"
              << "RLZ parsing, parallel over files (" << args.threads
              << " threads used, " << available << " available)\n"
              << "========================================\n";

    std::vector<Symbol> reference;
    const double reference_ms =
        time_ms([&] { reference = load_reference<Symbol>(args.reference); });
    std::cerr << "[1] reference: " << reference.size() << " bytes in "
              << reference_ms << " ms\n";

    std::vector<SAType> suffix_array;
    const double sa_ms = time_ms(
        [&] { suffix_array = load_suffix_array<SAType>(args.suffix_array); });
    std::cerr << "[2] suffix array: " << suffix_array.size() << " entries in "
              << sa_ms << " ms\n";
    if (suffix_array.size() != reference.size()) {
      throw std::runtime_error("suffix array and reference sizes differ");
    }
    const std::size_t n = reference.size();

    Inputs inputs;
    inputs.files = load_input_list(args.filenames);
    if (args.max_files != 0 && inputs.files.size() > args.max_files) {
      inputs.files.resize(args.max_files);
    }
    if (args.preload) {
      const double load_ms = time_ms([&] {
        for (const std::string& file : inputs.files) {
          inputs.loaded.push_back(load_input<Symbol>(file));
          inputs.total_bytes += inputs.loaded.back().size();
        }
      });
      std::cerr << "[3] input files: " << inputs.files.size() << ", "
                << mb(inputs.total_bytes) << " MB preloaded in " << load_ms
                << " ms\n";
    } else {
      for (const std::string& file : inputs.files) {
        inputs.total_bytes += fs::file_size(file);
      }
      std::cerr << "[3] input files: " << inputs.files.size()
                << " (loaded by the threads)\n";
    }

    std::vector<Result> results;
    const auto report = [&](const Result& r) {
      std::cerr << "    " << std::left << std::setw(17) << r.name << std::right
                << std::fixed << std::setprecision(1) << " build "
                << std::setw(8) << r.build_ms << " ms, load " << std::setw(7)
                << r.load_ms << " ms | wall " << std::setw(9) << r.wall_ms
                << " ms, work " << std::setw(10) << r.work_ms
                << " ms (x" << std::setprecision(2)
                << (r.wall_ms == 0.0 ? 0.0 : r.work_ms / r.wall_ms)
                << " parallel), slowest file " << std::setprecision(1)
                << r.longest_ms << " ms, " << r.phrases << " phrases\n";
      results.push_back(r);
    };

    std::cerr << "[4] parsers, one at a time\n";

    // ---------- sa-binary-search ----------
    if (wanted("sa-binary-search")) {
      Result r;
      r.name = "sa-binary-search";
      r.needed_bytes = n * sizeof(SAType) + n;
      run_parallel(inputs, args.threads, [&](const std::vector<Symbol>& input) {
        return lzFactorize<Symbol, SAType>(input, reference, suffix_array);
      }, r);
      report(r);
    }

    // ---------- lrf-ms ----------
    if (wanted("lrf-ms")) {
      Result r;
      r.name = "lrf-ms";
      std::unique_ptr<LRFMS<Symbol, SAType>> lrf;
      r.build_ms = time_ms(
          [&] { lrf = std::make_unique<LRFMS<Symbol, SAType>>(reference, suffix_array); });
      std::size_t rmq_size = 1;
      while ((rmq_size << 7) < n) rmq_size <<= 1;
      r.own_bytes = 3 * n * sizeof(std::int32_t) + 2 * rmq_size * sizeof(int);
      r.needed_bytes = r.own_bytes + n * sizeof(SAType) + n;
      run_parallel(inputs, args.threads, [&](const std::vector<Symbol>& input) {
        const MatchingStatistics ms = lrf->computeMatchingStatistics(input);
        Triples phrases;
        std::size_t i = 0;
        while (i < input.size()) {
          const auto [position, length] = ms[i];
          if (length <= 1) {
            phrases.emplace_back(i, static_cast<std::size_t>(input[i]), 1);
            ++i;
          } else {
            phrases.emplace_back(i, position, length);
            i += length;
          }
        }
        return phrases;
      }, r);
      report(r);
    }

    // ---------- pt16 (original) ----------
    if (wanted("pt16")) {
      Result r;
      r.name = "pt16";
      const std::string path = args.table;
      std::unique_ptr<pt16_v1::PT16RLZParser<Symbol, SAType>> parser;
      r.build_ms = time_ms([&] {
        fs::remove(path);
        pt16_v1::build_pt16_table(reference, suffix_array, path);
      });
      r.load_ms = time_ms([&] {
        parser = std::make_unique<pt16_v1::PT16RLZParser<Symbol, SAType>>(
            reference, suffix_array, path);
      });
      r.own_bytes = parser->stats().approx_bytes;
      r.needed_bytes = r.own_bytes + n * sizeof(SAType) + n;
      run_parallel(inputs, args.threads, [&](const std::vector<Symbol>& input) {
        return parser->lzFactorize(input);
      }, r);
      report(r);
    }

    // ---------- pt16-v2 ----------
    if (wanted("pt16-v2")) {
      Result r;
      r.name = "pt16-v2";
      const std::string path = args.table + ".v2";
      std::unique_ptr<PT16RLZParser<Symbol, SAType>> parser;
      r.build_ms = time_ms([&] {
        fs::remove(path);
        build_pt16_table(reference, suffix_array, path);
      });
      r.load_ms = time_ms([&] {
        parser = std::make_unique<PT16RLZParser<Symbol, SAType>>(reference,
                                                                 suffix_array, path);
      });
      r.own_bytes = parser->stats().approx_bytes;
      r.needed_bytes = r.own_bytes + n * sizeof(SAType) + n;
      run_parallel(inputs, args.threads, [&](const std::vector<Symbol>& input) {
        return parser->lzFactorize(input);
      }, r);
      report(r);
    }

    // ---------- sassy ----------
    if (wanted("sassy")) {
      Result r;
      r.name = "sassy";
      const std::string path = args.table + ".sassy";
      std::unique_ptr<PT16SassyLookup> lookup;
      r.build_ms = time_ms([&] {
        fs::remove(path);
        build_pt16_sassy_table(reference, suffix_array, path);
      });
      r.load_ms = time_ms([&] { lookup = std::make_unique<PT16SassyLookup>(path); });
      r.own_bytes = lookup->stats().approx_bytes;
      r.needed_bytes = r.own_bytes + n;  // no suffix array
      run_parallel(inputs, args.threads, [&](const std::vector<Symbol>& input) {
        return lookup->lzFactorize(input, reference);
      }, r);
      report(r);
    }

    // The suffix array is not needed by powered; free it before loading it.
    std::vector<SAType>().swap(suffix_array);

    // ---------- powered ----------
#ifdef WITH_POWERED
    if (wanted("powered") && !args.powered_index.empty()) {
      Result r;
      r.name = "powered";
      const std::string path = args.powered_index;
      const std::size_t dot = path.find_last_of('.');
      const std::string data_file =
          dot == std::string::npos ? path + "_data"
                                   : path.substr(0, dot) + "_data" + path.substr(dot);
      for (const std::string& file : {path, data_file}) {
        if (!fs::exists(file)) {
          throw std::runtime_error("powered index file not found: " + file);
        }
      }
      std::unique_ptr<bbwt::non_rle<>> index;
      r.load_ms = time_ms([&] { index = std::make_unique<bbwt::non_rle<>>(path); });
      std::uint64_t data_bytes = 0;
      {
        std::ifstream header(path, std::ios::binary);
        header.read(reinterpret_cast<char*>(&data_bytes), sizeof(data_bytes));
      }
      r.own_bytes = sizeof(bbwt::non_rle<>) + data_bytes +
                    257 * sizeof(std::uint64_t) +
                    index->gca_.size() * sizeof(std::uint64_t) +
                    fs::file_size(data_file);
      r.needed_bytes = r.own_bytes;
      run_parallel(inputs, args.threads, [&](const std::vector<Symbol>& input) {
        std::vector<std::tuple<std::uint64_t, std::uint64_t>> phrases;
        const std::string_view view(reinterpret_cast<const char*>(input.data()),
                                    input.size());
        index->parse_tuples(view, phrases, 4);
        return phrases;
      }, r);
      report(r);
    }
#endif
    if (wanted("powered") && args.powered_index.empty()) {
      std::cerr << "    (powered: not run; give --powered-index"
#ifndef WITH_POWERED
                   " and build with -DWITH_POWERED"
#endif
                   ")\n";
    }

    // ---------- Totals ----------

    const double base_wall = results.empty() ? 0.0 : results.front().wall_ms;
    std::cerr << "[5] summary (" << args.threads << " of " << available
              << " available threads, "
              << inputs.files.size() << " files, " << mb(inputs.total_bytes)
              << " MB)\n"
              << "    parser              wall ms   vs first      MB/s   parallel"
                 "      phrases    own MB  needs MB\n";
    for (const Result& r : results) {
      std::cerr << "    " << std::left << std::setw(17) << r.name << std::right
                << std::fixed << std::setprecision(1) << std::setw(10)
                << r.wall_ms << std::setprecision(2) << std::setw(10)
                << (r.wall_ms == 0.0 ? 0.0 : base_wall / r.wall_ms)
                << std::setprecision(1) << std::setw(10)
                << (r.wall_ms == 0.0 ? 0.0
                                     : mb(inputs.total_bytes) / (r.wall_ms / 1000.0))
                << std::setprecision(2) << std::setw(10)
                << (r.wall_ms == 0.0 ? 0.0 : r.work_ms / r.wall_ms)
                << std::setw(13) << r.phrases << std::setprecision(1)
                << std::setw(10) << mb(r.own_bytes) << std::setw(10)
                << mb(r.needed_bytes) << '\n';
    }
    std::cerr << "    peak RSS " << peak_rss_mb() << " MB\n";

    if (!args.results.empty()) {
      std::ofstream csv(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "parser,threads,threads_available,files,input_bytes,build_ms,"
             "load_ms,wall_ms,work_ms,"
             "slowest_file_ms,phrases,own_bytes,needed_bytes\n";
      for (const Result& r : results) {
        csv << r.name << ',' << args.threads << ',' << available << ','
            << inputs.files.size() << ','
            << inputs.total_bytes << ',' << std::setprecision(3) << r.build_ms
            << ',' << r.load_ms << ',' << r.wall_ms << ',' << r.work_ms << ','
            << r.longest_ms << ',' << r.phrases << ',' << r.own_bytes << ','
            << r.needed_bytes << '\n';
      }
    }

    std::cout << "threads=" << args.threads << '\n'
              << "threads_available=" << available << '\n'
              << "files=" << inputs.files.size() << '\n'
              << "input_bytes=" << inputs.total_bytes << '\n';
    for (const Result& r : results) {
      std::cout << r.name << "_wall_ms=" << r.wall_ms << '\n'
                << r.name << "_work_ms=" << r.work_ms << '\n'
                << r.name << "_build_ms=" << r.build_ms << '\n'
                << r.name << "_load_ms=" << r.load_ms << '\n'
                << r.name << "_phrases=" << r.phrases << '\n'
                << r.name << "_own_bytes=" << r.own_bytes << '\n'
                << r.name << "_needed_bytes=" << r.needed_bytes << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
