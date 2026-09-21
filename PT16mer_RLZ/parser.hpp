#pragma once

// The baseline RLZ parser that the PT16 and cached variants are compared
// against. The functions it is built from are in rlz_common.hpp.

#include <cstddef>
#include <tuple>
#include <vector>

#include "rlz_common.hpp"

template<typename T1, typename T2>
std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> lzFactorize(const std::vector<T1>& input,
                                                                           const std::vector<T1>& ref,
                                                                           const std::vector<T2>& sa) {
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> spl_vec;
    std::size_t i = 0;

    while (i < input.size()) {
        auto [pos, len] = rlz::computeLZFactorAt(input, ref, sa, i);

        if (len <= 1) {
            pos = static_cast<std::size_t>(input.at(i));
            len = 1;
        }

        spl_vec.push_back({i, pos, len});

        i += len;
    }

    return spl_vec;
}
