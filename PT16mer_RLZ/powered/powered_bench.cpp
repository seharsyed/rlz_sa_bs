// Powered RLZ (../../RLZ_powered) inside our benchmark conventions: a first,
// single-threaded experiment.
//
// Differences from RLZ_powered's own rlz_parser, which reads its file list
// and every input itself and times that too:
//
//   - the powered index (the .bwt written by RLZ_powered's `make_bwt`)
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
// Build (x86-64 only: the RLZ_powered headers use x86 intrinsics). The
// block sizes must be the ones RLZ_powered's make_bwt was
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
//                   [--pt16-table REF_four.pt16 | --pt16] [--escape] [--quiet]
//
// Variants (powered_pt16_parse.hpp), each run on every file after powered,
// timed the same way:
//   powered-escape       with --escape: powered, but a phrase whose interval
//                        is a single row is extended by character comparison;
//   powered-pt16-escape  with --pt16-table (a table from pt16_powered_build)
//                        or --pt16 (the table built in memory): a PT16 lookup
//                        (sassy layout) at each phrase start, and the escape
//                        (a singleton table entry has no row to rank from).
// They must decode. Their phrase lengths are compared with powered's for
// information: the variants take only matches inside the reference (and a
// greedy tail), where powered's cyclic parse may run over the reference's end;
// rlz_suite checks the left-to-right variants against sa-binary-search.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "../../RLZ_powered/include/types.hpp"  // bbwt::non_rle
#include "../pt16_utils.hpp"                    // loaders, time_ms
#include "powered_pt16_parse.hpp"
#include "pt16_powered.hpp"

using Symbol = unsigned char;
using Phrase = std::tuple<std::uint64_t, std::uint64_t>;  // (length, position)

namespace {

struct PoweredArgs {
  std::string index;       // RLZ_powered index (.bwt from make_bwt)
  std::string reference;   // the plain reference, for the decode check
  std::string filenames;   // our input list
  std::string results;     // optional CSV
  std::size_t max_files = 0;
  bool quiet = false;  // per-file lines only for files that fail
  // Characters per (meta)symbol: 4 for the powered index, 1 for an index
  // over the original alphabet (rlz_parser's --original).
  unsigned code_size = 4;
  // powered-pt16: a table file, or build the table in memory.
  std::string pt16_table;
  bool pt16 = false;
  bool escape = false;  // the escape-on-singleton variants
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
    } else if (option == "--quiet") {
      args.quiet = true;
    } else if (option == "--max-files") {
      args.max_files = std::stoull(require_value(i, argc, argv));
    } else if (option == "--code-size") {
      args.code_size = static_cast<unsigned>(
          std::stoul(require_value(i, argc, argv)));
    } else if (option == "--pt16-table") {
      args.pt16_table = require_value(i, argc, argv);
      args.pt16 = true;
    } else if (option == "--pt16") {
      args.pt16 = true;
    } else if (option == "--escape") {
      args.escape = true;
    } else if (option == "--help" || option == "-h") {
      std::cout << "Usage: " << argv[0]
                << " --index REF_four.bwt --reference REF --filenames LIST\n"
                   "  [--code-size 4]  4: powered index, 1: original alphabet\n"
                   "  [--results CSV]  per-file CSV\n"
                   "  [--max-files N]  only the first N input files\n"
                   "  [--quiet]        per-file lines only for failing files\n"
                   "  [--pt16-table P] also run powered-pt16-escape with this "
                   "table\n"
                   "  [--pt16]         also run powered-pt16-escape, table built "
                   "in memory\n"
                   "  [--escape]       also run powered-escape\n";
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
  if ((args.pt16 || args.escape) && args.code_size != 4) {
    throw std::runtime_error("the powered variants need --code-size 4");
  }
  return args;
}

// Rebuilds the input from the phrases (kept right to left, as the parser
// produces them) against the cyclic reference, as RLZ_powered's
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

    // The variants: the table (loaded or built) and their parsers.
    using Index = bbwt::non_rle<>;
    using Parser = PoweredPT16Parser<Index>;
    std::unique_ptr<PT16PoweredTable> table;
    double table_ms = 0.0;
    if (args.pt16) {
      table = std::make_unique<PT16PoweredTable>();
      if (!args.pt16_table.empty()) {
        table_ms = time_ms(
            [&] { *table = PT16PoweredTable::load(args.pt16_table); });
        std::cerr << "    PT16 table loaded in " << table_ms << " ms: ";
      } else {
        table_ms = time_ms([&] {
          *table = PT16PoweredTable::build(reference, index->sa_, nullptr);
        });
        std::cerr << "    PT16 table built in " << table_ms << " ms: ";
      }
      std::cerr << table->entries() << " entries, "
                << table->bytes() / (1024.0 * 1024.0) << " MB\n";
    }

    struct Variant {
      std::string name;
      std::unique_ptr<Parser> parser;
      double total_ms = 0.0;
      std::size_t files_ok = 0;
      std::size_t lengths_differ = 0;  // files where they differ from powered's
      Parser::Stats total;
    };
    std::vector<Variant> variants;
    const auto add_variant = [&](std::string name, const PT16PoweredTable* t,
                                 const std::vector<Symbol>* r) {
      variants.push_back({std::move(name),
                          std::make_unique<Parser>(*index, t, r), 0.0, 0, {}});
    };
    if (args.escape) add_variant("powered-escape", nullptr, &reference);
    if (args.pt16) add_variant("powered-pt16-escape", table.get(), &reference);

    std::ofstream csv;
    if (!args.results.empty()) {
      csv.open(args.results);
      if (!csv) throw std::runtime_error("cannot create " + args.results);
      csv << "file,input_bytes,non_acgt,parse_ms,phrases,decoded_ok";
      for (const Variant& v : variants) {
        csv << ',' << v.name << "_ms," << v.name << "_ok";
      }
      csv << '\n';
    }

    if (!args.quiet) std::cerr << "[4] files\n";

    std::size_t total_bytes = 0, total_phrases = 0, files_ok = 0;
    double total_parse_ms = 0.0;
    std::vector<Phrase> phrases;

    std::vector<Phrase> variant_phrases;
    std::size_t total_powered_wrapping = 0;

    for (std::size_t f = 0; f < files.size(); ++f) {
      // Loading is not timed (as in our other benchmarks).
      const std::vector<Symbol> input = load_input<Symbol>(files[f]);

      std::size_t non_acgt = 0;
      for (const Symbol c : input) non_acgt += is_acgt(c) ? 0 : 1;

      // This file's lines; with --quiet printed only if the file fails.
      std::ostringstream file_log;

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

      file_log << "[" << f + 1 << "/" << files.size() << "] " << files[f]
                << "  (" << input.size() << " bytes";
      if (non_acgt != 0) file_log << ", " << non_acgt << " non-ACGT";
      file_log << ")\n    parse " << std::fixed << std::setprecision(2)
                << parse_ms << " ms, " << phrases.size() << " phrases, "
                << "average length "
                << (phrases.empty() ? 0.0
                                    : static_cast<double>(input.size()) /
                                          static_cast<double>(phrases.size()))
                << ", decode " << (ok ? "OK" : "DIFFERS");
      if (!ok) {
        file_log << " (first at " << difference << ", decoded "
                  << decoded_length << " of " << input.size() << " bytes)";
      }
      file_log << '\n';

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
        file_log << "    first phrases in text order (length, position):";
        std::size_t shown = 0;
        for (auto it = phrases.rbegin(); it != phrases.rend() && shown < 5;
             ++it, ++shown) {
          file_log << " (" << std::get<0>(*it) << ", " << std::get<1>(*it)
                    << ")";
        }
        const std::uint64_t first_position = std::get<1>(phrases.back());
        file_log << "\n    input[0..40)                : "
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

      // powered's phrases that run over the reference's end (its index is
      // cyclic), for the length comparison below.
      std::size_t powered_wrapping = 0;
      for (const Phrase& phrase : phrases) {
        const std::uint64_t length = std::get<0>(phrase) & ~(std::uint64_t{1} << 63);
        powered_wrapping += std::get<1>(phrase) + length > reference.size() ? 1 : 0;
      }
      total_powered_wrapping += powered_wrapping;

      // The variants on the same input, timed the same way.
      bool variants_ok = true;
      std::vector<std::pair<double, bool>> variant_results;
      for (Variant& v : variants) {
        Parser::Stats stats;
        std::vector<Phrase>().swap(variant_phrases);
        const double ms =
            time_ms([&] { v.parser->parse(view, variant_phrases, stats); });
        // The check is the decode. The lengths are compared with powered's
        // (zero-length phrases aside, which powered emits when a phrase ends
        // at the input start) for information only: the variants take only
        // matches inside the reference and a greedy tail, powered's cyclic
        // parse may run over the reference's end and drops a partial
        // extension before its last short phrase. (rlz_suite checks the
        // left-to-right variants against sa-binary-search exactly.)
        std::size_t k = 0, m = 0;
        bool lengths_equal = true;
        while (lengths_equal) {
          while (k < phrases.size() && std::get<0>(phrases[k]) == 0) ++k;
          while (m < variant_phrases.size() &&
                 std::get<0>(variant_phrases[m]) == 0) {
            ++m;
          }
          if (k == phrases.size() || m == variant_phrases.size()) {
            lengths_equal =
                k == phrases.size() && m == variant_phrases.size();
            break;
          }
          lengths_equal =
              std::get<0>(phrases[k]) == std::get<0>(variant_phrases[m]);
          if (lengths_equal) {
            ++k;
            ++m;
          }
        }
        std::size_t variant_decoded = 0;
        const bool decodes =
            decode_check(variant_phrases, reference, input, variant_decoded) ==
            input.size();
        const bool variant_ok = decodes;

        file_log << "    " << v.name << ' ' << ms << " ms (x"
                  << (ms == 0.0 ? 0.0 : parse_ms / ms) << " vs powered), lengths "
                  << (lengths_equal ? "identical to powered's"
                                    : "differ from powered's")
                  << ", decode "
                  << (decodes ? "OK" : "DIFFERS");
        if (stats.lookups != 0) {
          file_log << ", 16-mer hits " << stats.hits << ", misses "
                    << stats.misses;
        }
        if (stats.escapes != 0) file_log << ", escapes " << stats.escapes;
        file_log << '\n';
        if (!lengths_equal) {
          file_log << "    first length difference at phrase " << k
                    << " (from the right) of " << phrases.size() << " / "
                    << variant_phrases.size() << "; powered runs over the "
                    << "reference's end in " << powered_wrapping << " phrases\n";
        }

        v.total_ms += ms;
        v.files_ok += variant_ok ? 1 : 0;
        v.lengths_differ += lengths_equal ? 0 : 1;
        v.total.lookups += stats.lookups;
        v.total.hits += stats.hits;
        v.total.misses += stats.misses;
        v.total.escapes += stats.escapes;
        variants_ok = variants_ok && variant_ok;
        variant_results.emplace_back(ms, variant_ok);
      }

      if (csv.is_open()) {
        csv << files[f] << ',' << input.size() << ',' << non_acgt << ','
            << std::setprecision(3) << parse_ms << ',' << phrases.size() << ','
            << (ok ? "YES" : "NO");
        for (const auto& [ms, variant_ok] : variant_results) {
          csv << ',' << ms << ',' << (variant_ok ? "YES" : "NO");
        }
        csv << '\n';
      }

      total_bytes += input.size();
      total_phrases += phrases.size();
      total_parse_ms += parse_ms;
      files_ok += ok ? 1 : 0;

      if (!args.quiet || !ok || !variants_ok) {
        std::cerr << file_log.str();
      }
    }

    std::cerr << std::fixed << std::setprecision(2) << "[5] totals over " << files.size() << " files ("
              << total_bytes << " bytes)\n"
              << "    index load " << index_ms << " ms\n"
              << "    parse      " << total_parse_ms << " ms, " << total_phrases
              << " phrases, decode OK " << files_ok << "/" << files.size()
              << " files\n";
    bool all_variants_ok = true;
    for (const Variant& v : variants) {
      std::cerr << "    " << v.name << ' ' << v.total_ms << " ms (x"
                << (v.total_ms == 0.0 ? 0.0 : total_parse_ms / v.total_ms)
                << " vs powered), decode OK " << v.files_ok << "/"
                << files.size() << " files, lengths differ from powered's in "
                << v.lengths_differ << " (powered runs over the reference's end "
                << "in " << total_powered_wrapping << " phrases)";
      if (v.total.lookups != 0) {
        std::cerr << ", 16-mer hits " << v.total.hits << ", misses "
                  << v.total.misses << " of " << v.total.lookups << " lookups";
      }
      if (v.total.escapes != 0) std::cerr << ", escapes " << v.total.escapes;
      std::cerr << '\n';
      all_variants_ok = all_variants_ok && v.files_ok == files.size();
    }
    if (table) {
      std::cerr << "    PT16 table " << (args.pt16_table.empty() ? "build " : "load ")
                << table_ms << " ms, " << table->bytes() / (1024.0 * 1024.0)
                << " MB\n";
    }

    std::cout << "files=" << files.size() << '\n'
              << "input_bytes=" << total_bytes << '\n'
              << "index_load_ms=" << index_ms << '\n'
              << "powered_parse_ms=" << total_parse_ms << '\n'
              << "powered_phrases=" << total_phrases << '\n'
              << "decode_ok=" << (files_ok == files.size() ? "YES" : "NO")
              << '\n';
    if (table) {
      std::cout << "pt16_table_ms=" << table_ms << '\n'
                << "pt16_table_bytes=" << table->bytes() << '\n';
    }
    for (const Variant& v : variants) {
      std::string key = v.name;
      std::replace(key.begin(), key.end(), '-', '_');
      std::cout << key << "_parse_ms=" << v.total_ms << '\n'
                << key << "_ok=" << (v.files_ok == files.size() ? "YES" : "NO")
                << '\n';
      if (v.total.lookups != 0) {
        std::cout << key << "_lookups=" << v.total.lookups << '\n'
                  << key << "_hits=" << v.total.hits << '\n'
                  << key << "_misses=" << v.total.misses << '\n';
      }
      if (v.total.escapes != 0) {
        std::cout << key << "_escapes=" << v.total.escapes << '\n';
      }
    }

    return files_ok == files.size() && all_variants_ok ? EXIT_SUCCESS : 4;
  } catch (const std::exception& error) {
    std::cerr << "\nERROR: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
