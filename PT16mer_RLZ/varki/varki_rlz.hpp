#pragma once

// Varki's RLZ (../../RLZ-Varki: FM-index of the reversed reference, sdsl-lite)
// under our benchmark conventions: the index is built once from the reference
// in memory, and an input already in memory is parsed into an in-memory
// phrase list -- the part we time, as for every other parser.
//
// Varki's RLZ_CHAR::compress / parse (include/rlz_algo_char.h) read each
// input file from disk one character at a time (ifstream::get / peek) inside
// the parse loop, and may split a file into chunks parsed by different
// threads (-t), which changes the parse at the chunk boundaries. Here there is
// no chunking at all: every input is one chunk, and parallelism is only
// across files (rlz_parallel), as for the other parsers. parse() below is
// RLZ_CHAR::parse's loop for one chunk (max_len not set), reading from memory
// instead of the file; it uses Varki's own FM_Wrapper::backward_match and
// get_suffix_array_value on Varki's index type (rlz_fm_index_t), and Varki's
// position formula. Varki's per-call chrono timers around backward_match and
// the SA lookup are instrumentation and are left out. The build is
// load_reverse_reference + construct_im + calculate_occs.
//
// The parse: the input is read forwards; each character is one backward
// step on the reversed reference, which extends the match forwards in the
// original. So each phrase is the longest prefix of the rest of the input
// that occurs in the reference: the greedy left-to-right parse, with the
// same phrase lengths as sa-binary-search (positions may be another
// occurrence; a length-1 phrase is a reference position, not a literal).
// Every input character must occur in the reference (Varki exits otherwise).
//
// Phrases are (position, length), in input order, into the plain (not
// cyclic) reference: phrase k is reference[position, position + length).
//
// This header has no sdsl types: the implementation (varki_rlz.cpp, which
// also compiles Varki's fm_wrapper.cpp) is built apart as C++17, as Varki is
// -- sdsl-lite does not compile as C++20 (std::result_of) -- and linked in.
// From PT16mer_RLZ/, with sdsl built by varki/build_sdsl.sh:
//   S=../RLZ-Varki/build/sdsl
//   g++ -std=c++17 -O3 -DNDEBUG -I$S/include -I../RLZ-Varki/include \
//       -c varki/varki_rlz.cpp -o varki_rlz.o
//   g++ -std=c++2a ... -DWITH_VARKI PROGRAM.cpp varki_rlz.o \
//       -L$S/lib -lsdsl -ldivsufsort -ldivsufsort64
// (build_if_needed in rlz_datasets.sh does this.)

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

class VarkiRLZ {
 public:
  using Phrase = std::pair<std::uint64_t, std::uint64_t>;  // (position, length)

  // Builds the index of the reversed reference, as RLZ_CHAR does
  // (load_reverse_reference, construct_im(fm_index, ref_content, 1),
  // calculate_occs).
  explicit VarkiRLZ(const std::vector<unsigned char>& reference);
  ~VarkiRLZ();
  VarkiRLZ(const VarkiRLZ&) = delete;
  VarkiRLZ& operator=(const VarkiRLZ&) = delete;

  // RLZ_CHAR::parse for one chunk covering all of text[0, n), from memory.
  // Read-only: threads may parse different inputs at once.
  void parse(const unsigned char* text, std::size_t n,
             std::vector<Phrase>& phrases) const;

  // The index's size (sdsl's own count) plus the F column.
  std::size_t bytes() const;

 private:
  struct Index;  // varki_rlz.cpp: Varki's rlz_fm_index_t and its F column
  std::unique_ptr<Index> index_;
};
