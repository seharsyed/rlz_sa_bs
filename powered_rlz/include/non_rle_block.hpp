//
// Library containing blocks employing no run length compression.
// The computation of the rank is inspired by Quad rank
// But extended to support the computation even for larger alphabets.
// Koerkamp, R. Groot. "QuadRank: Engineering a High Throughput Rank." arXiv preprint arXiv:2602.04103 (2026)
//

#pragma once

#include <immintrin.h>
#include <cstdint>
#include <iostream>
#include <ostream>

namespace bbwt {
template <uint32_t block_size, class alphabet_type_, bool avx = false>
class non_rle_block {
  public:
    typedef alphabet_type_ alphabet_type;
  private:
    static_assert(block_size <= ~uint32_t(0) >> 1);

    uint64_t MASK;

   public:
    static const constexpr bool has_members = false;
    static const constexpr uint32_t cap = block_size;
    static const constexpr uint32_t scratch_blocks = 2;
    static const constexpr uint32_t min_size = 2;
    static const constexpr uint32_t padding_bytes = avx ? 32 : 0;
    static const constexpr uint32_t max_size = block_size;

    static constexpr uint64_t scratch_size(uint32_t i) {
        if (i == 0) {
            return 8;
        } else {
            return max_size;
        }
    }

    non_rle_block() {
        MASK = 1ULL << 63;
    }

    non_rle_block(const non_rle_block& other) = delete;
    non_rle_block(non_rle_block&& other) = delete;
    non_rle_block& operator=(non_rle_block&& other) = delete;
    non_rle_block& operator=(const non_rle_block&) = delete;

    uint32_t append(uint8_t head, uint32_t length, uint8_t** scratch) {
        uint8_t W = alphabet_type::width;

        uint64_t* offset = reinterpret_cast<uint64_t*>(scratch[0]);
        offset[0] >>= 3; //because the offset point to first 8-bit space, not 64-bit

        if (MASK < (1ULL << 63)) {
            offset[0] -= W;
        }

        uint64_t* data = reinterpret_cast<uint64_t*>(scratch[1]);

        for (uint32_t i = 0; i < length; i++) {

            uint8_t is_set = 1 << (W-1);
            for (uint8_t j = 0; j < W; j++) {
                if (head & is_set) {
                    data[offset[0] + j] |= MASK;
                }
                is_set >>= 1;
            }

            if (MASK == 1) {
                MASK <<= 63;
                offset[0] += W;
            } else {
                MASK >>= 1;
            }
        }

        if (MASK < (1ULL << 63)) offset[0] += alphabet_type::width;

        offset[0] <<= 3;
        return offset[0];
    }

    uint8_t at(uint64_t location) const {
        uint64_t mask = 1ULL << (63 - (location & 63));
        const uint64_t* data = reinterpret_cast<const uint64_t*>(this);

        uint8_t result = 0;
        for (uint8_t j = 0; j < alphabet_type::width; j++) {
            result <<= 1;
            result |= (data[(location >> 6) * alphabet_type::width +j] & mask) > 0;
        }

        return result;
    }

    uint32_t rank(uint8_t c, uint32_t location) const {
        uint32_t result = 0;
        const uint64_t* data = reinterpret_cast<const uint64_t*>(this);
        uint8_t W = alphabet_type::width;
        uint8_t mask;
        uint64_t r;

        uint64_t div = (location >> 6) * W;
        for (uint64_t i = 0; i < div; i += W) {
            mask = 1 << (W-1);

            if (c&mask) {
                r = data[i];
            } else {
                r = ~data[i];
            }

            for (uint8_t j = 1; j < W; j++) {
                mask >>= 1;
                if (c&mask) {
                    r &= data[i+j];
                } else {
                    r &= ~data[i+j];
                }
            }
            result += __builtin_popcountll(r);
        }

        if (location & 63) {
            mask = 1 << (W-1);

            if (c&mask) {
                r = data[div];
            } else {
                r = ~data[div];
            }

            for (uint8_t j = 1; j < W; j++) {
                mask >>= 1;
                if (c&mask) {
                    r &= data[div+j];
                } else {
                    r &= ~data[div+j];
                }
            }

            r >>= (64 - (location & 63));
            result += __builtin_popcountll(r);
        }

        return result;
    }

    std::pair<uint64_t, uint64_t> extension(uint64_t a, uint64_t b, uint8_t c) const {
        const uint64_t* data = reinterpret_cast<const uint64_t*>(this);

        uint64_t start = (a >> 6) << 3;
        uint64_t end = (b >> 6) << 3;

        uint8_t extension_size = 0;
        uint64_t dex = a;

        uint64_t mask = ~0ULL >> (a&63);

        for (uint64_t i = start; i <= end; i += 8) {
            if (i > start) {
                mask = ~0ULL;
            }

            if (i == end) {
                mask &= ~0ULL << (63 - (b&63));
            }

            uint8_t char_mask = 0b00000001;

            for (int j = 3; j > 0; j--) {
                if (c&char_mask) {
                    mask &= data[i+(j<<1)+1];
                } else {
                    mask &= ~data[i+(j<<1)+1];
                }
                char_mask <<= 1;
                if (c&char_mask) {
                    mask &= data[i+(j<<1)];
                } else {
                    mask &= ~data[i+(j<<1)];
                }
                char_mask <<= 1;

                if (mask == 0)
                    break;

                if (extension_size < (4-j)) {
                    extension_size++;
                    dex = std::__countl_zero(mask) | (i<<3);
                }
            }
        }

        return {extension_size, dex};
    }

    uint64_t commit(uint8_t** scratch) {
        uint64_t bytes = reinterpret_cast<uint64_t*>(scratch[0])[0];
        uint8_t* data = reinterpret_cast<uint8_t*>(this);
        uint8_t* src = reinterpret_cast<uint8_t*>(scratch[1]);
        std::memcpy(data, src, bytes);
        return bytes;
    }

    void clear() {}
    static void write_statics(std::fstream&) {return; }
    static uint64_t load_statics(std::fstream&) {return 0; }

    void print(uint32_t sb) const {
        const uint64_t* data = reinterpret_cast<const uint64_t*>(this);

        uint64_t pos = 1ULL << 63;
        uint32_t block_num = 0;
        for (uint32_t i = 0; i < sb; i++) {
            uint8_t current = 0;
            for (uint8_t j = 0; j < alphabet_type::width; j++) {
                if (data[block_num + j] & pos)
                    current |= 1 << (alphabet_type::width-j -1);
            }
            std::cerr << "symbol " << (int)current << std::endl;
            if (pos == 1) {
                block_num += alphabet_type::width;
                pos = 1ULL << 63;
            } else {
                pos >>= 1;
            }
        }
    }
};
}  // namespace bbwt