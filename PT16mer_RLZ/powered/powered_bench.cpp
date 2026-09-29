// Powered RLZ (../../powered_rlz) inside our benchmark conventions: a first,
// single-threaded experiment.
//
// Differences from powered_rlz's own rlz_parser, which reads its file list
// and every input itself and times that too:
//
//   - the powered index (the .bwt written by powered_rlz's `transform -b`)
//     is given as an argument and loaded once (timed on its own);
//   - the input list is our usual one, and each input is loaded with our
//     loader, untimed -- only the parse itself is timed, as in our other
//     benchmarks;
//   - one thread, one file at a time;
//   - each parse is decoded against the reference (untimed) and must
//     reproduce the input exactly.
//
// Powered RLZ parses right to left (each phrase is the longest suffix of
// the remaining prefix, by backward search), so its phrases are not the
// same as a left-to-right greedy parse; the decode check is what proves
// the parse correct. Positions are into the cyclic reference.
//
// The inputs and the reference must be plain ACGT (no N, no newline):
// powered RLZ maps any other byte to one of A/C/G/T. Inputs with other
// bytes are reported.
//
// Build (x86-64 only: the powered_rlz headers use x86 intrinsics). The
// block sizes must be the ones powered_rlz's make_bwt / transform were
// built with (its Makefile's defaults are below):
//
//   g++ -std=c++2a -O3 -march=native -DNDEBUG \
//       -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 \
//       powered/powered_bench.cpp -o powered_bench
//
// Run:
//
//   ./powered_bench --index REF_four.bwt --reference REF --filenames LIST
//                   [--code-size 4] [--results CSV] [--max-files N]

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "../../powered_rlz/include/types.hpp"  // bbwt::non_rle
#include "../pt16_utils.hpp"                    // loaders, time_ms

using Symbol = unsigned char;
using Phrase = std::tuple<std::uint64_t, std::uint64_t>;  // (length, position)

namespace {

struct PoweredArgs {
  std::string index;       // powered_rlz index (.bwt from `transform -b`)
  std::string reference;   // the plain reference, for the decode check
  std::string filenames;   // our input list
  std::string results;     // optional CSV
  std::size_t max_files = 0;
  // Characters per (meta)symbol: 4 for the powered index, 1 for an index
  // over the original alphabet (rlz_parser's --original).
  unsigned code_size = 4;
};

std::string require_value(int& i, int argc, char** argv) {
  if (i + 1 >= argc) {
    throw std::runtime_error(std::string("missing value for ") + argv[i]);
  }
  return argv[++i];
}

PoweredArgs parse_powered_args(int argc, char** argv) {
  PoweredArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--index") {
      args.index = require_value(i, argc, argv);
    } else if (option == "--reference") {
      args.reference = require_value(i, argc, argv);
    } else if (option == "--filenames") {
      args.filenames = require_value(i, argc, argv);
    } else if (option == "--results") {
      args.results = require_value(i, argc, argv);
    } else if (option == "--max-files") {
      args.max_files = std::stoull(require_value(i, argc, argv));
    } else if (option == "--code-size") {
      args.code_size = static_cast<unsigned>(
          std::stoul(require_value(i, argc, argv)));
    } else if (option == "--help" || option == "-h") {
      std::cout << "Usage: " << argv[0]
                << " --index REF_four.bwt --reference REF --filenames LIST\n"
                   "  [--code-size 4]  4: powered index, 1: original alphabet\n"
                   "  [--results CSV]  per-file CSV\n"
                   "  [--max-files N]  only the first N input files\n";
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("unknown argument: " + option);
    }
  }
  if (args.index.empty() || args.reference.empty() || args.filenames.empty()) {
    throw std::runtime_error("index, reference and filenames are required");
  }
  if (args.code_size != 1 && args.code_size != 4) {
    throw std::runtime_error("--code-size must be 1 or 4");
  }
  return args;
}

// Rebuilds the input from the phrases (kept right to left, as the parser
// produces them) against the cyclic reference, as powered_rlz's
// decompress does, and compares it with `input`. Returns the first
// differing position, or input.size() if they are equal.
std::size_t decode_check(const std::vector<Phrase>& phrases,
                         const std::vector<Symbol>& reference,
                         const std::vector<Symbol>& input,
                         std::size_t& decoded_length) {
  const std::uint64_t literal = std::uint64_t{1} << 63;
  const std::size_t m = reference.size();
  std::size_t at = 0;
  std::size_t first_difference = input.size();

  for (auto it = phrases.rbegin(); it != phrases.rend(); ++it) {
    std::uint64_t length = std::get<0>(*it);
    const std::uint64_t position = std::get<1>(*it);

    // A run of a character absent from the reference (decompress's
    // format): `length` copies of the byte `position`.
    const bool is_literal = (length & literal) != 0;
    if (is_literal) length -= literal;

    for (std::uint64_t j = 0; j < length; ++j, ++at) {
      const Symbol c = is_literal ? static_cast<Symbol>(position)
                                  : reference[(position + j) % m];
      if (first_difference == input.size() &&
          (at >= input.size() || input[at] != c)) {
        first_difference = at;
      }
    }
  }

  decoded_length = at;
  if (first_difference == input.size() && at != input.size()) {
    first_difference = std::min(at, input.size());
  }
  return first_difference;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::cerr << "========================================\n"
              << "powered RLZ benchmark (single-threaded)\n"
              << "========================================\n";

    const PoweredArgs args = parse_powered_args(argc, argv);

    std::vector<Symbol> reference;
    const double reference_ms =
        time_ms([&] { reference = load_reference<Symbol>(args.reference); });
    std::cerr << "[1] reference: " << reference.size() << " bytes in "
              << reference_ms << " ms\n";

    std::size_t reference_non_acgt = 0;
    for (const Symbol c : reference) reference_non_acgt += is_acgt(c) ? 0 : 1;
    if (reference_non_acgt != 0) {
      std::cerr << "    WARNING: the reference has " << reference_non_acgt
                << " non-ACGT bytes; powered RLZ expects plain ACGT\n";
    }

    std::vector<std::string> files = load_input_list(args.filenames);
    if (args.max_files != 0 && files.size() > args.max_files) {
      files.resize(args.max_files);
    }
    std::cerr << "[2] input files: " << files.size() << '\n';

    // The powered index, loaded once.
    std::cerr << "[3] loading the powered index " << args.index << " ...\n";
    std::unique_ptr<bbwt::non_rle<>> index;
    const double index_ms = time_ms(
        [&] { index = std::make_unique<bbwt::non_rle<>>(args.index); });
    std::cerr << "    loaded in " << index_ms << " ms (" << index->size()
              << " BWT symbols)\n";

    std::ofstream csv;
    if (!args.results.empty()) {
      csv.open(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "file,input_bytes,non_acgt,parse_ms,phrases,decoded_ok\n";
    }

    std::cerr << "[4] files\n";

    std::size_t total_bytes = 0, total_phrases = 0, files_ok = 0;
    double total_parse_ms = 0.0;
    std::vector<Phrase> phrases;

    for (std::size_t f = 0; f < files.size(); ++f) {
      // Loading is not timed (as in our other benchmarks).
      const std::vector<Symbol> input = load_input<Symbol>(files[f]);

      std::size_t non_acgt = 0;
      for (const Symbol c : input) non_acgt += is_acgt(c) ? 0 : 1;

      // Only the parse is timed.
      // A fresh, empty phrase list per file, as rlz_parser has (freed here,
      // untimed).
      std::vector<Phrase>().swap(phrases);
      const std::string_view view(reinterpret_cast<const char*>(input.data()),
                                  input.size());
      const double parse_ms = time_ms([&] {
        index->parse_tuples(view, phrases,
                            static_cast<std::uint8_t>(args.code_size));
      });

      std::size_t decoded_length = 0;
      const std::size_t difference =
          decode_check(phrases, reference, input, decoded_length);
      const bool ok = difference == input.size();

      std::cerr << "[" << f + 1 << "/" << files.size() << "] " << files[f]
                << "  (" << input.size() << " bytes";
      if (non_acgt != 0) std::cerr << ", " << non_acgt << " non-ACGT";
      std::cerr << ")\n    parse " << std::fixed << std::setprecision(2)
                << parse_ms << " ms, " << phrases.size() << " phrases, "
                << "average length "
                << (phrases.empty() ? 0.0
                                    : static_cast<double>(input.size()) /
                                          static_cast<double>(phrases.size()))
                << ", decode " << (ok ? "OK" : "DIFFERS");
      if (!ok) {
        std::cerr << " (first at " << difference << ", decoded "
                  << decoded_length << " of " << input.size() << " bytes)";
      }
      std::cerr << '\n';

      // On a mismatch: the first phrases in text order, and the input next
      // to what the first phrase decodes to, to see what kind of mismatch
      // it is (order, position offset, strand, ...).
      if (!ok && !phrases.empty()) {
        const auto show = [&](const std::vector<Symbol>& text, std::size_t at,
                              std::size_t count, bool cyclic) {
          std::string s;
          for (std::size_t k = 0; k < count; ++k) {
            const std::size_t p = cyclic ? (at + k) % text.size() : at + k;
            if (p >= text.size()) break;
            s.push_back(static_cast<char>(text[p]));
          }
          return s;
        };
        std::cerr << "    first phrases in text order (length, position):";
        std::size_t shown = 0;
        for (auto it = phrases.rbegin(); it != phrases.rend() && shown < 5;
             ++it, ++shown) {
          std::cerr << " (" << std::get<0>(*it) << ", " << std::get<1>(*it)
                    << ")";
        }
        const std::uint64_t first_position = std::get<1>(phrases.back());
        std::cerr << "\n    input[0..40)                : "
                  << show(input, 0, 40, false)
                  << "\n    reference[first position ..]: "
                  << show(reference, first_position, 40, true)
                  << "\n    last phrase (vector front)  : ("
                  << std::get<0>(phrases.front()) << ", "
                  << std::get<1>(phrases.front()) << "), input end: "
                  << show(input, input.size() >= 40 ? input.size() - 40 : 0,
                          40, false)
                  << '\n';
      }

      if (csv.is_open()) {
        csv << files[f] << ',' << input.size() << ',' << non_acgt << ','
            << std::setprecision(3) << parse_ms << ',' << phrases.size() << ','
            << (ok ? "YES" : "NO") << '\n';
      }

      total_bytes += input.size();
      total_phrases += phrases.size();
      total_parse_ms += parse_ms;
      files_ok += ok ? 1 : 0;
    }

    std::cerr << "[5] totals over " << files.size() << " files ("
              << total_bytes << " bytes)\n"
              << "    index load " << index_ms << " ms\n"
              << "    parse      " << total_parse_ms << " ms, " << total_phrases
              << " phrases, decode OK " << files_ok << "/" << files.size()
              << " files\n";

    std::cout << "files=" << files.size() << '\n'
              << "input_bytes=" << total_bytes << '\n'
              << "index_load_ms=" << index_ms << '\n'
              << "powered_parse_ms=" << total_parse_ms << '\n'
              << "powered_phrases=" << total_phrases << '\n'
              << "decode_ok=" << (files_ok == files.size() ? "YES" : "NO")
              << '\n';

    return files_ok == files.size() ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
