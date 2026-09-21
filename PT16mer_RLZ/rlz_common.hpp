#pragma once

// Neutral RLZ building blocks: plain suffix-array search and greedy longest-
// match, plus file reading.
//
// They are shared by the baseline parser (parser.hpp) and by every PT16 and
// cached variant, which call them by their rlz:: names for whatever their own
// index cannot answer. Because both sides run the same code here, agreement
// between the baseline and a variant does not independently check these
// functions.

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace rlz {

template<typename T>
std::vector<T> read_file(const char* filename) {
    std::ifstream ifs(filename, std::ios::binary);

    const auto begin = ifs.tellg();
    ifs.seekg(0, std::ios::end);
    const auto end = ifs.tellg();
    const std::size_t len = (end - begin) / sizeof(T);
    ifs.seekg(0);

    std::vector<T> v(len, 0);

    for (std::size_t i = 0; i < len; ++i) {
        ifs.read(reinterpret_cast<char*>(v.data() + i), sizeof(T));
    }

    ifs.close();

    return v;
}

// SA only needs .size() and operator[], so this and the two binary searches
// below run equally well over a full std::vector<T2> suffix array or over a
// std::span<const T2> naming just a slice of one (for example, one 16-mer's
// own occurrences: see PT16SassyLookup::find_longest_matching_factor).
template<typename T1, typename SA>
inline T1 safe_ref_symbol(const std::vector<T1>& ref,
                          const SA& sa,
                          const std::size_t sa_index,
                          const std::size_t offset) {
    if (sa_index >= sa.size()) {
        throw std::out_of_range("safe_ref_symbol: sa_index out of range");
    }

    const std::size_t suffix_pos = static_cast<std::size_t>(sa[sa_index]);

    if (suffix_pos >= ref.size() || offset >= ref.size() - suffix_pos) {
        return T1{};
    }

    return ref[suffix_pos + offset];
}

template<typename T1, typename SA>
inline std::optional<std::int64_t> binarySearchLB(const std::vector<T1>& ref, const SA& sa,
                                                 const std::int64_t lo, const std::int64_t hi,
                                                 const std::int64_t offset, const T1 c) {
    std::int64_t low = lo;
    std::int64_t high = hi;

    while (low <= high) {
        const std::int64_t mid = low + (high - low) / 2;
        const auto midVal = safe_ref_symbol(ref, sa, mid, offset);

        if (midVal < c) {
            low = mid + 1;
        } else if (midVal > c) {
            high = mid - 1;
        } else { //midVal == c
            if (mid == lo) {
                return mid; // leftmost occ of key found
            }

            const auto midValLeft = safe_ref_symbol(ref, sa, mid - 1, offset);
            if (midValLeft == midVal) {
                high = mid - 1; //discard mid and the ones to the right of mid
            } else { //midValLeft must be less than midVal == c
                return mid; //leftmost occ of key found
            }
        }
    }

    return {}; // key not found.
}

template<typename T1, typename SA>
inline std::optional<std::int64_t> binarySearchRB(const std::vector<T1>& ref, const SA& sa,
                                                 const std::int64_t lo, const std::int64_t hi,
                                                 const std::int64_t offset, const T1 c) {
    std::int64_t low = lo;
    std::int64_t high = hi;

    while (low <= high) {
        const std::int64_t mid = low + (high - low) / 2;
        const auto midVal = safe_ref_symbol(ref, sa, mid, offset);

        if (midVal < c) {
            low = mid + 1;
        } else if (midVal > c) {
            high = mid - 1;
        } else { //midVal == c
            if (mid == hi) {
                return mid; // rightmost occ of key found
            }

            const auto midValRight = safe_ref_symbol(ref, sa, mid + 1, offset);
            if (midValRight == midVal) {
                low = mid + 1; //discard mid and the ones to the left of mid
            } else { //midValRight must be greater than midVal == c
                return mid; //rightmost occ of key found
            }
        }
    }
    return {}; // key not found.
}

template<typename T1, typename T2>
std::tuple<std::size_t, std::size_t> computeLZFactorAt(const std::vector<T1>& input,
                                                       const std::vector<T1>& ref,
                                                       const std::vector<T2>& sa,
                                                       const std::size_t input_pos) {
    std::size_t offset = 0;
    std::size_t j = input_pos;

    std::size_t match = 0;
    std::size_t nlb = 0;
    std::size_t nrb = ref.size() - 1;

    while (j < input.size()) {
        if (nlb == nrb) {
            if (safe_ref_symbol(ref, sa, nlb, offset) != input[j]) {
                break;
            }
        }
        else {
            if (const auto opt = binarySearchLB(ref, sa, nlb, nrb, offset, input.at(j))) {
                nlb = opt.value();
            } else {
                break;
            }

            if (const auto opt = binarySearchRB(ref, sa, nlb, nrb, offset, input.at(j))) {
                nrb = opt.value();
            } else {
                break;
            }
        }

        match = sa[nlb];
        ++j;
        ++offset;
    }

    return {match, offset};
}

}  // namespace rlz
