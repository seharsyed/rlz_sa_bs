// Implementation of varki_rlz.hpp: Varki's index and parse loop (see the
// header). Compiled as C++17 with sdsl-lite's headers and RLZ-Varki/include
// on the include path.

#include "varki_rlz.hpp"

#include <string>
#include <tuple>

#include <sdsl/suffix_arrays.hpp>

#include "fm_wrapper.h"                        // FM_Wrapper, rlz_fm_index_t
#include "../../RLZ-Varki/src/fm_wrapper.cpp"  // its definitions

struct VarkiRLZ::Index {
  rlz_fm_index_t fm_index;
  std::vector<std::size_t> occs = std::vector<std::size_t>(256, 0);
};

VarkiRLZ::VarkiRLZ(const std::vector<unsigned char>& reference)
    : index_(std::make_unique<Index>()) {
  // load_reverse_reference + construct_im(fm_index, ref_content, 1)
  const std::string reversed(reference.rbegin(), reference.rend());
  sdsl::construct_im(index_->fm_index, reversed, 1);

  // calculate_occs: the compressed F column (characters before each code).
  std::vector<std::size_t>& occs = index_->occs;
  for (const char c : reversed) ++occs[static_cast<unsigned char>(c)];
  std::size_t running_total = 0;
  for (std::size_t i = 0; i < 256; ++i) {
    const std::size_t frequency = occs[i];
    occs[i] = running_total;
    running_total += frequency;
  }
}

VarkiRLZ::~VarkiRLZ() = default;

void VarkiRLZ::parse(const unsigned char* text, const std::size_t n,
                     std::vector<Phrase>& phrases) const {
  const rlz_fm_index_t& fm_index = index_->fm_index;
  const std::vector<std::size_t>& occs = index_->occs;
  FM_Wrapper fm_support;  // stateless; one per call keeps threads apart
  const std::size_t bwt_size = fm_index.bwt.size();

  std::size_t pattern_len = 0;
  std::size_t prev_left = 0, prev_right = bwt_size;
  std::size_t next_left = 0, next_right = bwt_size;

  // Varki's formula: the reversed SA value is where the match ends in the
  // original (with the virtual sentinel's offset); subtract the length.
  const auto position = [&](const std::size_t row) {
    const std::size_t sa_pos = fm_support.get_suffix_array_value(fm_index, row);
    return static_cast<std::uint64_t>(bwt_size - 1 - sa_pos - pattern_len);
  };

  std::size_t k = 0;  // the next character to match (retried on a mismatch)
  while (k < n) {
    const char next_char = static_cast<char>(text[k]);
    ++pattern_len;

    const std::tuple<std::size_t, std::size_t> next_ranges =
        fm_support.backward_match(fm_index, occs,
                                  std::make_tuple(prev_left, prev_right),
                                  next_char);
    next_left = std::get<0>(next_ranges);
    next_right = std::get<1>(next_ranges);

    if (next_left == next_right) {
      // Mismatch: emit the match so far, retry this character.
      --pattern_len;
      phrases.emplace_back(position(prev_left), pattern_len);
      prev_left = 0;
      prev_right = bwt_size;
      pattern_len = 0;
    } else if (k + 1 == n) {
      // End of the input inside a match.
      phrases.emplace_back(position(next_left), pattern_len);
      ++k;
    } else {
      prev_left = next_left;
      prev_right = next_right;
      ++k;
    }
  }
}

std::size_t VarkiRLZ::bytes() const {
  return sdsl::size_in_bytes(index_->fm_index) +
         index_->occs.size() * sizeof(std::size_t);
}
