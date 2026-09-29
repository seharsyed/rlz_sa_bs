#pragma once

#include <cstring>
#include <iostream>
#include <string>
#include <bitset>
#include <cstdint>

namespace bbwt {
template <class block_type_>
class super_block {
   public:
    typedef block_type_ block_type;
    typedef typename block_type::alphabet_type alphabet_type;
    static_assert((uint64_t(1) << 32) % block_type::cap == 0);
    static const constexpr uint64_t blocks =
        (uint64_t(1) << 32) / block_type::cap;
    static const constexpr uint32_t cap = block_type::cap;
   private:
    uint64_t offsets_[blocks];
   public:
    

    super_block() = delete;
    super_block(const super_block& other) = delete;
    super_block(super_block&& other) = delete;
    super_block& operator=(const super_block& other) = delete;
    super_block& operator=(super_block&& other) = delete;

    uint8_t at(uint64_t i) const {
        uint32_t block_i = i / cap;
        //std::cerr << " block " << block_i << std::endl;
        const block_type* block =
            reinterpret_cast<const block_type*>(data() + offsets_[block_i]);
        return block->at(i % cap);
    }

    uint32_t rank(uint8_t c, uint32_t i) const {
        uint32_t block_i = i / cap;
        //std::cerr << "rank(" << int(c) << ", " << i << ")" << std::endl;
        __builtin_prefetch(data() + offsets_[block_i]);
        const alphabet_type* alpha = reinterpret_cast<const alphabet_type*>(
            data() + offsets_[block_i] - alphabet_type::size());
        uint32_t res = alpha->p_sum(c);
        //std::cerr << res << " from previous blocks " << std::endl;
        const block_type* block =
            reinterpret_cast<const block_type*>(data() + offsets_[block_i]);
        res += block->rank(c, i % cap);
        return res;
    }

    std::pair<uint64_t, uint64_t> extension(uint64_t a, uint64_t b, uint8_t c) {
        uint64_t block_a = a / cap;
        uint64_t block_b = b / cap;

        std::pair<uint64_t, uint64_t> res;
        const block_type* start_block = reinterpret_cast<const block_type*>(data() + offsets_[block_a]);
        if (block_a == block_b) {
            res = start_block->extension(a % cap, b % cap, c);
            return {res.first, res.second + block_a*cap};
        }

        res = start_block->extension(a % cap, cap-1, c);
        res = {res.first, res.second + block_a*cap};
        block_a++;
        while (block_a < block_b) {
            const block_type* block = reinterpret_cast<const block_type*>(data() + offsets_[block_a]);
            auto r = block->extension(0, cap-1, c);

            if (r.first > res.first) res = {r.first, r.second+block_a*cap};
            block_a++;
        }
        const block_type* end_block = reinterpret_cast<const block_type*>(data() + offsets_[block_b]);
        auto r = end_block->extension(0, b % cap, c);
        if (r.first > res.first) res = {r.first, r.second+block_b*cap};

        return res;
    }

    template <class dtype>
    void print_block(uint32_t idx, uint32_t n_bytes) const {
        dtype* dp = reinterpret_cast<dtype*>(data() + offsets_[idx]);
        for (uint32_t i = 0; i < n_bytes; i++) {
            std::cerr << std::bitset<sizeof(dtype) * 8>(dp[i]) << std::endl;
        }
    }

    alphabet_type* get_psums(uint32_t i) const {
        return reinterpret_cast<alphabet_type*>(data() + offsets_[i] - alphabet_type::size());
    }

    void print(uint64_t s) const {
        bool done = false;
        for (uint32_t i = 0; i < blocks; i++) {
            std::cerr << "sub-block " << i << ": " << std::endl;
            uint64_t sb = cap * (i + 1);
            if (sb > s) {
                sb = s % cap;
                done = true;
            }
            sb = cap;
            uint32_t block_i = (sb - 1) / cap;
            const alphabet_type* alpha = reinterpret_cast<const alphabet_type*>(
                data() + offsets_[block_i] - alphabet_type::size());
            alpha->print();
            const block_type* block =
                reinterpret_cast<const block_type*>(data() + offsets_[block_i]);
            block->print(sb);
            if (done) break;
        }
    }

    static uint64_t write_statics(std::fstream&) {return 0; }
    static uint64_t load_statics(std::fstream&) {return 0; }

   private:
    const uint8_t* data() const {
        return reinterpret_cast<const uint8_t*>(this) + sizeof(super_block);
    }
};
}  // namespace bbwt
