#pragma once

#include <cstdint>

#include "block_rlbwt.hpp"
#include "byte_alphabet.hpp"
#include "custom_alphabet.hpp"
#include "super_block.hpp"
#include "two_byte_block.hpp"
#include "alphabet.hpp"
#include "non_rle_block.hpp"

#ifndef LARGE_BLOCK_SIZE
#define LARGE_BLOCK_SIZE 16384
#endif

#ifndef SMALL_BLOCK_SIZE
#define SMALL_BLOCK_SIZE 512
#endif

namespace bbwt {

template <uint32_t block_size = SMALL_BLOCK_SIZE>
using non_rle_build = block_rlbwt<
super_block<non_rle_block<block_size, byte_alphabet<uint32_t>>>,
byte_alphabet<uint64_t>>;

template <uint32_t block_size = SMALL_BLOCK_SIZE>
using non_rle = block_rlbwt<
    super_block<non_rle_block<block_size, byte_alphabet<uint32_t>>>,
    byte_alphabet<uint64_t>>;

template <uint32_t block_size = SMALL_BLOCK_SIZE>
using two_byte_build = block_rlbwt<
    super_block<two_byte_block<block_size, custom_alphabet<uint32_t>>>,
    custom_alphabet<uint64_t>>;

template <uint32_t block_size = SMALL_BLOCK_SIZE>
using two_byte = block_rlbwt<
    super_block<two_byte_block<block_size, alphabet<uint32_t>>>,
    alphabet<uint64_t>>;
}  // namespace bbwt
