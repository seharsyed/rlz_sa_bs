// Checks the powered PT16 table (pt16_powered.hpp) against the powered
// index's own backward search, and measures what it would give the parse.
//
//   [3] entries: in increasing reversed-key order, covering all rows; for
//       every entry (count() on an even sample of --samples) its position
//       holds its 16-mer, looking the 16-mer up returns the entry (a range:
//       its row interval, whose first row has that position), and powered's
//       count() (first metacharacter + three rank steps) gives its size.
//   [4] random 16-mers (--samples, fixed seed): a hit exactly when count() >
//       0, with the same size; on a miss the returned (length, position) is
//       the longest suffix that occurs: the reference has it at that
//       position, and count() of one character more is 0.
//   [5] with --filenames: each input is parsed by powered (parse_tuples); at
//       every phrase start (its right end i, with i >= 16) the lookup of
//       seq[i-16, i) must hit exactly when the phrase has length >= 16, and
//       on a miss return the phrase's length. Also times the table lookups
//       against count() on those same 16-mers (a smoke test off the server).
//
// Build (x86-64, GCC), from PT16mer_RLZ/:
//
//   g++ -std=c++2a -O3 -march=native -DNDEBUG \
//       -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 \
//       powered/pt16_powered_check.cpp -o pt16_powered_check
//
// Run:
//
//   ./pt16_powered_check --reference REF --index REF_four.bwt
//       [--table REF_four.pt16] [--filenames LIST] [--max-files N]
//       [--samples N]
//
// Without --table the table is built in memory.

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <tuple>

#include "../../powered_rlz/include/types.hpp"  // bbwt::non_rle
#include "pt16_powered.hpp"

namespace {

struct CheckArgs {
  std::string reference, index, table, filenames;
  std::size_t max_files = 0;
  std::size_t samples = 1'000'000;
};

CheckArgs parse_check_args(int argc, char** argv) {
  CheckArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--reference") {
      args.reference = require_value(i, argc, argv);
    } else if (option == "--index") {
      args.index = require_value(i, argc, argv);
    } else if (option == "--table") {
      args.table = require_value(i, argc, argv);
    } else if (option == "--filenames") {
      args.filenames = require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = std::stoull(require_value(i, argc, argv));
    } else if (option == "--samples") {
      args.samples = std::stoull(require_value(i, argc, argv));
    } else if (option == "--help" || option == "-h") {
      std::cout << "Usage: " << argv[0]
                << " --reference REF --index REF_four.bwt\n"
                   "  [--table PATH]      table file (else built in memory)\n"
                   "  [--filenames LIST]  inputs for the phrase-start check\n"
                   "  [--max-files N]     only the first N inputs\n"
                   "  [--samples N]       entries / random 16-mers checked "
                   "(default 1000000)\n";
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }
  if (args.reference.empty() || args.index.empty()) {
    throw std::runtime_error("--reference and --index are required");
  }
  return args;
}

// The 16 characters of the cyclic reference at `position`.
std::string cyclic_16mer(const std::vector<unsigned char>& reference,
                         const std::uint64_t position) {
  std::string kmer(KMER_LENGTH, 'A');
  for (std::uint32_t j = 0; j < KMER_LENGTH; ++j) {
    kmer[j] = static_cast<char>(reference[(position + j) % reference.size()]);
  }
  return kmer;
}

// powered's own backward search for a pattern (any length >= 1): its first,
// possibly partial, metacharacter is the last ((length - 1) % 4) + 1
// characters, then one rank step per four characters. Returns the interval
// size.
template <typename Index>
std::uint64_t powered_count(const Index& index, const std::string_view pattern) {
  const auto offset = static_cast<std::uint8_t>((pattern.size() - 1) % 4 + 1);
  return index.count(pattern, offset, 4);
}

// Whether the cyclic reference has `pattern` at `position`.
bool occurs_at(const std::vector<unsigned char>& reference,
               const std::uint64_t position, const std::string_view pattern) {
  for (std::size_t j = 0; j < pattern.size(); ++j) {
    if (reference[(position + j) % reference.size()] !=
        static_cast<unsigned char>(pattern[j])) {
      return false;
    }
  }
  return true;
}

std::size_t failures = 0;

void report_failure(const std::string& what) {
  if (failures < 10) std::cerr << "    FAIL: " << what << '\n';
  ++failures;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const CheckArgs args = parse_check_args(argc, argv);

    std::vector<unsigned char> reference =
        load_reference<unsigned char>(args.reference);
    std::cerr << "[1] reference: " << reference.size() << " characters\n";

    std::unique_ptr<bbwt::non_rle<>> index =
        std::make_unique<bbwt::non_rle<>>(args.index);
    const auto& gca = index->gca_;
    std::cerr << "[2] powered index: " << index->size() << " rows\n";

    PT16PoweredTable table;
    if (args.table.empty()) {
      table = PT16PoweredTable::build(reference, gca, nullptr);
      std::cerr << "    table built in memory: ";
    } else {
      table = PT16PoweredTable::load(args.table);
      std::cerr << "    table loaded: ";
    }
    std::cerr << table.entries() << " entries\n";
    if (table.rows() != index->size()) {
      throw std::runtime_error("table has " + std::to_string(table.rows()) +
                               " rows, index " +
                               std::to_string(index->size()));
    }

    // ---------- [3] entries against powered's count ----------
    {
      const std::uint64_t entries = table.entries();
      const std::uint64_t stride = std::max<std::uint64_t>(
          1, entries / std::max<std::size_t>(1, args.samples));
      std::uint64_t checked = 0, e = 0, singletons = 0, rows = 0;
      std::uint32_t previous_key = 0;
      const std::size_t before = failures;

      table.for_each_entry([&](const PT16PoweredTable::EntryInfo& info) {
        if (e != 0 && info.reversed_key <= previous_key) {
          report_failure("entry " + std::to_string(e) +
                         " not in increasing reversed-key order");
        }
        previous_key = info.reversed_key;

        const std::string kmer = cyclic_16mer(reference, info.position);
        const PoweredStart s =
            table.lookup(reinterpret_cast<const unsigned char*>(kmer.data()));
        const std::uint64_t size = info.singleton ? 1 : info.b - info.a;
        rows += size;
        singletons += info.singleton ? 1 : 0;

        if (reverse_16mer_key(encode_cyclic_16mer(reference, info.position)) !=
            info.reversed_key) {
          report_failure("entry " + std::to_string(e) +
                         ": its position does not hold its 16-mer");
        } else if (!s.found || s.singleton != info.singleton ||
                   s.position != info.position ||
                   (!info.singleton && (s.a != info.a || s.b != info.b))) {
          report_failure("entry " + std::to_string(e) + " (" + kmer +
                         ") not returned by its own lookup");
        } else if (!info.singleton && gca[info.a] != info.position) {
          report_failure("entry " + std::to_string(e) +
                         ": range position is not its first row's");
        }
        if (e % stride == 0) {
          const std::uint64_t count = powered_count(*index, kmer);
          if (count != size) {
            report_failure(kmer + ": table size " + std::to_string(size) +
                           ", powered count " + std::to_string(count));
          }
          ++checked;
        }
        ++e;
      });
      if (e != entries || rows != table.rows()) {
        report_failure("entries cover " + std::to_string(rows) + " rows of " +
                       std::to_string(table.rows()));
      }
      std::cerr << "[3] entries: " << e << " (" << singletons
                << " singletons), covering " << rows << " rows, " << checked
                << " checked against powered count: "
                << (failures == before ? "OK" : "FAILED") << '\n';
    }

    // ---------- [4] random 16-mers ----------
    {
      std::mt19937_64 random(12345);
      const char dna[4] = {'A', 'C', 'G', 'T'};
      std::string kmer(KMER_LENGTH, 'A');
      std::uint64_t found = 0;
      std::array<std::uint64_t, KMER_LENGTH> miss_lengths{};
      const std::size_t before = failures;

      for (std::size_t s = 0; s < args.samples; ++s) {
        for (auto& c : kmer) c = dna[random() & 3U];
        const std::string_view view(kmer);
        const PoweredStart r =
            table.lookup(reinterpret_cast<const unsigned char*>(kmer.data()));
        const std::uint64_t count = powered_count(*index, view);

        if (r.found) {
          ++found;
          const std::uint64_t size = r.singleton ? 1 : r.b - r.a;
          if (size != count || !occurs_at(reference, r.position, view)) {
            report_failure("random " + kmer + ": table size " +
                           std::to_string(size) + ", powered count " +
                           std::to_string(count));
          }
          continue;
        }

        const std::uint32_t l = r.match_length;
        ++miss_lengths[std::min<std::uint32_t>(l, KMER_LENGTH - 1)];
        if (count != 0) {
          report_failure("random " + kmer + ": table miss, powered count " +
                         std::to_string(count));
        } else if (l == 0 || l >= KMER_LENGTH ||
                   !occurs_at(reference, r.position,
                              view.substr(KMER_LENGTH - l))) {
          report_failure("random " + kmer + ": miss length " +
                         std::to_string(l) + " at " +
                         std::to_string(r.position) +
                         " does not occur there");
        } else if (powered_count(*index, view.substr(KMER_LENGTH - l - 1)) !=
                   0) {
          report_failure("random " + kmer + ": miss length " +
                         std::to_string(l) + " is not the longest suffix");
        }
      }
      std::cerr << "[4] random 16-mers: " << args.samples << " checked, "
                << found << " found: "
                << (failures == before ? "OK" : "FAILED")
                << "\n    miss lengths:";
      for (std::uint32_t l = 0; l < KMER_LENGTH; ++l) {
        if (miss_lengths[l] != 0) std::cerr << ' ' << l << ':' << miss_lengths[l];
      }
      std::cerr << '\n';
    }

    // ---------- [5] phrase starts of powered's parse ----------
    if (!args.filenames.empty()) {
      std::vector<std::string> files = load_input_list(args.filenames);
      if (args.max_files != 0 && files.size() > args.max_files) {
        files.resize(args.max_files);
      }

      const std::uint64_t literal = std::uint64_t{1} << 63;
      std::uint64_t phrases = 0, starts = 0, hits = 0, left_end = 0;
      double lookup_ms = 0, count_ms = 0;
      std::uint64_t sink = 0;
      const std::size_t before = failures;

      for (const std::string& file : files) {
        const std::vector<unsigned char> input =
            load_input<unsigned char>(file);
        const std::string_view view(reinterpret_cast<const char*>(input.data()),
                                    input.size());

        std::vector<std::tuple<std::uint64_t, std::uint64_t>> parse;
        index->parse_tuples(view, parse, 4);

        std::vector<std::uint64_t> positions;  // phrase starts with i >= 16
        std::uint64_t i = input.size();
        for (const auto& phrase : parse) {
          const std::uint64_t length = std::get<0>(phrase) & ~literal;
          if (length == 0) continue;  // powered's zero-length phrase
          if (i >= KMER_LENGTH) {
            const PoweredStart s =
                table.lookup(input.data() + i - KMER_LENGTH);
            const bool expect_hit = length >= KMER_LENGTH;
            if (s.found != expect_hit ||
                (!s.found && s.match_length != length)) {
              report_failure(file + ": phrase ending at " + std::to_string(i) +
                             " has length " + std::to_string(length) +
                             ", the lookup " +
                             (s.found ? "hits"
                                      : "misses with length " +
                                            std::to_string(s.match_length)));
            }
            hits += s.found;
            ++starts;
            positions.push_back(i - KMER_LENGTH);
          } else {
            ++left_end;
          }
          i -= std::min(length, i);
          ++phrases;
        }

        lookup_ms += time_ms([&] {
          for (const std::uint64_t p : positions) {
            const PoweredStart s = table.lookup(input.data() + p);
            sink += s.a + s.position + s.found;
          }
        });
        count_ms += time_ms([&] {
          for (const std::uint64_t p : positions) {
            sink += powered_count(
                *index, std::string_view(view.data() + p, KMER_LENGTH));
          }
        });
      }

      std::cerr << "[5] phrase starts over " << files.size() << " inputs: "
                << phrases << " phrases, " << starts << " with a full 16-mer ("
                << left_end << " at the left end)\n"
                << "    hits (phrase >= 16) " << hits << ", misses resolved "
                << starts - hits << " (" << std::fixed << std::setprecision(1)
                << (starts == 0 ? 0.0 : 100.0 * (starts - hits) / starts)
                << "%): " << (failures == before ? "OK" : "FAILED") << '\n'
                << "    time over those " << starts << " 16-mers: table "
                << std::setprecision(2) << lookup_ms << " ms, powered count "
                << count_ms << " ms (x"
                << (lookup_ms == 0 ? 0.0 : count_ms / lookup_ms) << ")\n"
                << "    (sink " << (sink & 1) << ")\n";

      std::cout << "phrases=" << phrases << '\n'
                << "phrase_starts=" << starts << '\n'
                << "phrase_start_hits=" << hits << '\n'
                << "lookup_ms=" << lookup_ms << '\n'
                << "powered_count_ms=" << count_ms << '\n';
    }

    std::cerr << (failures == 0 ? "ALL OK" : "FAILURES: " + std::to_string(failures))
              << '\n';
    std::cout << "check_ok=" << (failures == 0 ? 1 : 0) << '\n';
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
