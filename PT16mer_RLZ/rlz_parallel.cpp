// RLZ parsing, parallel over input files: one parser at a time, T threads each
// taking the next input file and parsing it against the shared, read-only
// table or index.
//
// Parsers (as in rlz_suite.cpp): sa-binary-search, lrf-ms, pt16, pt16-v2,
// sassy, varki (built with -DWITH_VARKI; one chunk per file, see
// varki/varki_rlz.hpp), and powered-escape and powered-pt16-escape (built with
// -DWITH_POWERED, given --powered-index; see rlz_suite.cpp).
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
//                  [--threads T | --threads-list 1,2,4,...] [--per-file CSV]
//                  [--parsers a,b,...] [--powered-index REF_four.bwt]
//                  [--table PATH] [--max-files N] [--no-preload] [--results CSV]
//                  [--quiet]
//
// There is no per-file output; --quiet is accepted (and changes nothing) so
// that every benchmark takes the same flag.

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

#ifdef WITH_VARKI
#include "varki/varki_rlz.hpp"  // VarkiRLZ (sdsl-lite)
#endif
#ifdef WITH_POWERED
#include "../RLZ_powered/include/types.hpp"
#include "powered/powered_pt16_parse.hpp"  // PoweredPT16Parser
#include "powered/pt16_powered.hpp"        // PT16PoweredTable
#endif

using Symbol = unsigned char;
using SAType = std::uint32_t;
using Triples = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;
namespace fs = std::filesystem;

namespace {

// ---------- Arguments ----------

const std::vector<std::string> kAllParsers = {
    "sa-binary-search", "lrf-ms", "pt16", "pt16-v2", "sassy", "varki",
    "powered-escape", "powered-pt16-escape", "powered-fwd-escape",
    "powered-pt16-fwd-escape"};

struct ParallelArgs {
  std::string reference, suffix_array, filenames, table, powered_index, results;
  std::string powered_fwd_index;  // index of the reversed reference
  std::size_t max_files = 0;
  unsigned threads = 0;  // 0: all hardware threads
  // --threads-list: one parallel pass per thread count, on structures built
  // once (--threads T is the list {T}).
  std::vector<unsigned> threads_list;
  std::string per_file;  // CSV of every file's parse time (per parser, threads)
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
               "  [--threads-list 1,2,4]  one parallel pass per thread count, on\n"
               "                          structures built once (instead of --threads)\n"
               "  [--per-file PATH]       CSV: every file's parse time, per parser and\n"
               "                          thread count\n"
               "  [--parsers a,b,...]     subset of: sa-binary-search, lrf-ms, "
               "pt16, pt16-v2, sassy, varki, powered-escape, powered-pt16-escape,\n"
               "                          powered-fwd-escape, powered-pt16-fwd-escape\n"
               "  [--powered-index PATH]  RLZ_powered index (REF_four.bwt)\n"
               "  [--powered-fwd-index PATH]  the same, of the reversed reference "
               "(REF.rev_four.bwt)\n"
               "  [--table PATH]          PT16 table files (default: "
               "<reference>.parallel_pt16, .v2, .sassy)\n"
               "  [--max-files N]         only the first N input files\n"
               "  [--no-preload]          threads load their own files (the "
               "loading then counts in the wall time)\n"
               "  [--results PATH]        CSV, one row per parser\n"
               "  [--quiet]               accepted for uniformity (no per-file "
               "output anyway)\n";
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
    } else if (option == "--powered-fwd-index") {
      args.powered_fwd_index = require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = std::stoull(require_value(i, argc, argv));
    } else if (option == "--threads") {
      args.threads = static_cast<unsigned>(std::stoul(require_value(i, argc, argv)));
    } else if (option == "--threads-list") {
      std::stringstream list(require_value(i, argc, argv));
      std::string count;
      while (std::getline(list, count, ',')) {
        const unsigned t = static_cast<unsigned>(std::stoul(count));
        if (t == 0) throw std::runtime_error("--threads-list: 0 threads");
        args.threads_list.push_back(t);
      }
    } else if (option == "--per-file") {
      args.per_file = require_value(i, argc, argv);
    } else if (option == "--no-preload") {
      args.preload = false;
    } else if (option == "--quiet") {
      // no per-file output to suppress
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
  if (args.threads_list.empty()) {
    args.threads_list.push_back(args.threads);
  } else {
    args.threads = *std::max_element(args.threads_list.begin(),
                                     args.threads_list.end());
  }
#ifndef WITH_POWERED
  if (!args.powered_index.empty() || !args.powered_fwd_index.empty()) {
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
  unsigned threads = 0;
  std::vector<double> file_ms;            // per file (input order)
  std::vector<std::size_t> file_phrases;  // per file
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

  result.threads = threads;
  for (std::size_t f = 0; f < n; ++f) {
    result.work_ms += times[f];
    result.longest_ms = std::max(result.longest_ms, times[f]);
    result.phrases += phrases[f];
  }
  result.file_ms = std::move(times);
  result.file_phrases = std::move(phrases);
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
              << "RLZ parsing, parallel over files (up to " << args.threads
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
      std::cerr << "    " << std::left << std::setw(24) << r.name << std::right
                << std::setw(3) << r.threads << "t" << std::fixed
                << std::setprecision(1) << " build "
                << std::setw(8) << r.build_ms << " ms, load " << std::setw(7)
                << r.load_ms << " ms | wall " << std::setw(9) << r.wall_ms
                << " ms, work " << std::setw(10) << r.work_ms
                << " ms (x" << std::setprecision(2)
                << (r.wall_ms == 0.0 ? 0.0 : r.work_ms / r.wall_ms)
                << " parallel), slowest file " << std::setprecision(1)
                << r.longest_ms << " ms, " << r.phrases << " phrases\n";
      results.push_back(r);
    };

    // One parallel pass per thread count, on the structures already built
    // (`base` carries the parser's name, build/load times and sizes).
    const auto run_threads = [&](const Result& base, auto&& parse) {
      for (const unsigned t : args.threads_list) {
        Result r = base;
        run_parallel(inputs, t, parse, r);
        report(r);
      }
    };

    std::cerr << "[4] parsers, one at a time\n";

    // ---------- sa-binary-search ----------
    if (wanted("sa-binary-search")) {
      Result r;
      r.name = "sa-binary-search";
      r.needed_bytes = n * sizeof(SAType) + n;
      run_threads(r, [&](const std::vector<Symbol>& input) {
        return lzFactorize<Symbol, SAType>(input, reference, suffix_array);
      });
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
      run_threads(r, [&](const std::vector<Symbol>& input) {
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
      });
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
      run_threads(r, [&](const std::vector<Symbol>& input) {
        return parser->lzFactorize(input);
      });
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
      run_threads(r, [&](const std::vector<Symbol>& input) {
        return parser->lzFactorize(input);
      });
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
      run_threads(r, [&](const std::vector<Symbol>& input) {
        return lookup->lzFactorize(input, reference);
      });
    }

    // The suffix array is not needed by powered; free it before loading it.
    std::vector<SAType>().swap(suffix_array);

    // ---------- varki ----------
#ifdef WITH_VARKI
    if (wanted("varki")) {
      Result r;
      r.name = "varki";
      std::unique_ptr<VarkiRLZ> varki;
      r.build_ms = time_ms([&] { varki = std::make_unique<VarkiRLZ>(reference); });
      r.own_bytes = varki->bytes();
      r.needed_bytes = r.own_bytes;  // parsing reads the index only
      run_threads(r, [&](const std::vector<Symbol>& input) {
        std::vector<VarkiRLZ::Phrase> phrases;
        varki->parse(input.data(), input.size(), phrases);
        return phrases;
      });
    }
#else
    if (wanted("varki") && !args.parsers.empty()) {
      std::cerr << "    (varki: not run; build with -DWITH_VARKI, see "
                   "varki/build_sdsl.sh)\n";
    }
#endif

    // ---------- powered: escape, pt16-escape, and their fwd variants ----------
    // One index per direction for both of its parsers; each parse has its own
    // (unused) counters, so the parsers are safe to share between threads.
    // The fwd variants use the index of the reversed reference
    // (--powered-fwd-index) and the reference reversed in memory, and give
    // the greedy left-to-right parse (see powered/powered_pt16_parse.hpp).
    struct PoweredDirection {
      std::string index_path;
      bool forward;
      std::string escape_name, pt16_name;
    };
    const std::vector<PoweredDirection> directions = {
        {args.powered_index, false, "powered-escape", "powered-pt16-escape"},
        {args.powered_fwd_index, true, "powered-fwd-escape",
         "powered-pt16-fwd-escape"}};
    for (const PoweredDirection& direction : directions) {
      const bool any_wanted =
          wanted(direction.escape_name) || wanted(direction.pt16_name);
      if (!any_wanted) continue;
      if (direction.index_path.empty()) {
        std::cerr << "    (" << direction.escape_name << ", "
                  << direction.pt16_name << ": not run; give "
                  << (direction.forward ? "--powered-fwd-index"
                                        : "--powered-index")
#ifndef WITH_POWERED
                  << " and build with -DWITH_POWERED"
#endif
                  << ")\n";
        continue;
      }
#ifdef WITH_POWERED
      using Parse = PoweredPT16Parser<bbwt::non_rle<>>;
      const std::string& path = direction.index_path;
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
      const double index_ms =
          time_ms([&] { index = std::make_unique<bbwt::non_rle<>>(path); });
      std::uint64_t data_bytes = 0;
      {
        std::ifstream header(path, std::ios::binary);
        header.read(reinterpret_cast<char*>(&data_bytes), sizeof(data_bytes));
      }
      const std::size_t index_bytes =
          sizeof(bbwt::non_rle<>) + data_bytes + 257 * sizeof(std::uint64_t) +
          index->sa_.size() * sizeof(std::uint64_t) + fs::file_size(data_file);

      // The reference the table and the escape read: the reversed one for
      // the fwd variants (built at load, reported with the load).
      std::vector<Symbol> reversed;
      const double reverse_ms = time_ms([&] {
        if (direction.forward) reversed.assign(reference.rbegin(), reference.rend());
      });
      const std::vector<Symbol>& ref = direction.forward ? reversed : reference;

      const auto run_powered = [&](const std::string& name,
                                   const PT16PoweredTable* table, Result& r) {
        r.name = name;
        const Parse parser(*index, table, &ref, direction.forward);
        r.needed_bytes = r.own_bytes + n;  // + the reference (the escape)
        run_threads(r, [&](const std::vector<Symbol>& input) {
          std::vector<std::tuple<std::uint64_t, std::uint64_t>> phrases;
          const std::string_view view(reinterpret_cast<const char*>(input.data()),
                                      input.size());
          Parse::Stats stats;
          parser.parse(view, phrases, stats);
          return phrases;
        });
      };

      if (wanted(direction.escape_name)) {
        Result r;
        r.load_ms = index_ms + reverse_ms;
        r.own_bytes = index_bytes;
        run_powered(direction.escape_name, nullptr, r);
      }
      if (wanted(direction.pt16_name)) {
        Result r;
        std::unique_ptr<PT16PoweredTable> table;
        r.build_ms = time_ms([&] {
          table = std::make_unique<PT16PoweredTable>(
              PT16PoweredTable::build(ref, index->sa_, nullptr));
        });
        r.own_bytes = index_bytes + table->bytes();
        run_powered(direction.pt16_name, table.get(), r);
      }
#endif
    }

    // ---------- Totals ----------

    // "vs first": against the first parser run with the same thread count.
    const auto base_wall = [&](const unsigned threads) {
      for (const Result& r : results) {
        if (r.threads == threads) return r.wall_ms;
      }
      return 0.0;
    };
    std::cerr << "[5] summary (up to " << args.threads << " of " << available
              << " available threads, "
              << inputs.files.size() << " files, " << mb(inputs.total_bytes)
              << " MB)\n"
              << "    parser                   threads   wall ms   vs first      MB/s   parallel"
                 "      phrases    own MB  needs MB\n";
    for (const Result& r : results) {
      std::cerr << "    " << std::left << std::setw(24) << r.name << std::right
                << std::setw(8) << r.threads
                << std::fixed << std::setprecision(1) << std::setw(10)
                << r.wall_ms << std::setprecision(2) << std::setw(10)
                << (r.wall_ms == 0.0 ? 0.0 : base_wall(r.threads) / r.wall_ms)
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
        csv << r.name << ',' << r.threads << ',' << available << ','
            << inputs.files.size() << ','
            << inputs.total_bytes << ',' << std::fixed << std::setprecision(3)
            << r.build_ms
            << ',' << r.load_ms << ',' << r.wall_ms << ',' << r.work_ms << ','
            << r.longest_ms << ',' << r.phrases << ',' << r.own_bytes << ','
            << r.needed_bytes << '\n';
      }
    }

    if (!args.per_file.empty()) {
      std::ofstream csv(args.per_file);
      if (!csv) throw std::runtime_error("cannot create " + args.per_file);
      csv << "parser,threads,file,input_bytes,parse_ms,phrases\n";
      for (const Result& r : results) {
        for (std::size_t f = 0; f < r.file_ms.size(); ++f) {
          const std::size_t bytes = inputs.loaded.empty()
                                        ? fs::file_size(inputs.files[f])
                                        : inputs.loaded[f].size();
          csv << r.name << ',' << r.threads << ',' << inputs.files[f] << ','
              << bytes << ',' << std::fixed << std::setprecision(4) << r.file_ms[f]
              << ','
              << r.file_phrases[f] << '\n';
        }
      }
    }

    std::cout << "threads=" << args.threads << '\n'
              << "threads_available=" << available << '\n'
              << "files=" << inputs.files.size() << '\n'
              << "input_bytes=" << inputs.total_bytes << '\n';
    for (const Result& r : results) {
      // With several thread counts, the keys name theirs (_t4_wall_ms, ...).
      const std::string key =
          args.threads_list.size() > 1
              ? r.name + "_t" + std::to_string(r.threads)
              : r.name;
      std::cout << key << "_wall_ms=" << r.wall_ms << '\n'
                << key << "_work_ms=" << r.work_ms << '\n'
                << key << "_build_ms=" << r.build_ms << '\n'
                << key << "_load_ms=" << r.load_ms << '\n'
                << key << "_phrases=" << r.phrases << '\n'
                << key << "_own_bytes=" << r.own_bytes << '\n'
                << key << "_needed_bytes=" << r.needed_bytes << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
