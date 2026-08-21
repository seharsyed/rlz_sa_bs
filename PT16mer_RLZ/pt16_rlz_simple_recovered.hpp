#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "parser.hpp"


template <typename T1, typename T2>
class PT16RLZParser {
public:
    using input_type = std::vector<T1>;
    using reference_type = std::vector<T1>;
    using suffix_array_type = std::vector<T2>;

    using factor_type = std::tuple<std::size_t, std::size_t>;
    using phrase_type = std::tuple<std::size_t, std::size_t, std::size_t>;
    using phrase_vector_type = std::vector<phrase_type>;

    static constexpr std::uint32_t kmer_length = 16;
    static constexpr std::uint32_t bucket_size = 65536;
    static constexpr std::uint32_t low_bits = 16;
    static constexpr std::uint32_t low_mask = bucket_size - 1;
    static constexpr std::uint32_t number_of_buckets = 65536;

    struct Stats {
        std::size_t hits = 0;
        std::size_t misses = 0;
        std::size_t entries = 0;
        std::size_t approx_bytes = 0;
    };

private:
    struct LookupResult {
        bool found = false;
        std::uint32_t sa_start = 0;
        std::uint32_t sa_end = 0;
    };

    const reference_type* ref_ = nullptr;
    const suffix_array_type* sa_ = nullptr;

    std::array<std::uint8_t, 256> alphatab_{};

    // ---------- PT16 H/L representation ----------

    std::vector<std::uint32_t> H_;
    std::vector<std::uint16_t> L_;

    // Only the SA starting position of each PT16 interval is stored.
    std::vector<std::uint32_t> interval_starts_;

    mutable Stats stats_;

public:
    PT16RLZParser(
        const reference_type& ref,
        const suffix_array_type& sa,
        const std::string& pt16_path
    ) : ref_(&ref), sa_(&sa) {
        initialise_alphatab();
        load_hl(pt16_path);
    }


    factor_type computeLZFactorAt(const input_type& input, const std::size_t input_pos) {
        // Fewer than 16 characters remain: use ordinary RLZ.
        if (input.size() - input_pos < kmer_length) {
            return ::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);
        }

        // Pack the next 16 characters and look them up in PT16.
        const std::uint32_t key = pack_16mer(input, input_pos);
        const LookupResult result = lookup(key);

        // PT16 miss: use ordinary RLZ.
        if (!result.found) {
            return ::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);
        }

        std::size_t offset = kmer_length;
        std::size_t j = input_pos + kmer_length;
        std::size_t nlb = result.sa_start;
        std::size_t nrb = result.sa_end;

        // Range case: narrow the SA interval from character 17 onward.
        while (nlb < nrb && j < input.size()) {
            const auto lb = ::binarySearchLB<T1, T2>(*ref_, *sa_, nlb, nrb, offset, input[j]);

            if (!lb) {
                break;
            }

            const auto rb = ::binarySearchRB<T1, T2>(
                *ref_,
                *sa_,
                static_cast<std::size_t>(lb.value()),
                nrb,
                offset,
                input[j]
            );

            if (!rb) {
                break;
            }

            nlb = static_cast<std::size_t>(lb.value());
            nrb = static_cast<std::size_t>(rb.value());

            ++j;
            ++offset;
        }

        std::size_t match = static_cast<std::size_t>((*sa_)[nlb]);

        // Singleton case: extend directly from character 17 onward.
        if (nlb == nrb) {
            while (
                j < input.size() &&
                match + offset < ref_->size() &&
                (*ref_)[match + offset] == input[j]
            ) {
                ++j;
                ++offset;
            }
        }

        return {match, offset};
    }


    phrase_vector_type lzFactorize(const input_type& input) {
        phrase_vector_type spl_vec;
        std::size_t i = 0;

        while (i < input.size()) {
            auto [pos, len] = computeLZFactorAt(input, i);

            if (len <= 1) {
                pos = static_cast<std::size_t>(input[i]);
                len = 1;
            }

            spl_vec.push_back({i, pos, len});
            i += len;
        }

        return spl_vec;
    }


    const Stats& stats() const {
        return stats_;
    }


private:
    // ---------- Build alphatab once ----------

    void initialise_alphatab() {
        alphatab_[static_cast<unsigned char>('A')] = 0; // 00
        alphatab_[static_cast<unsigned char>('C')] = 1; // 01
        alphatab_[static_cast<unsigned char>('G')] = 2; // 10
        alphatab_[static_cast<unsigned char>('T')] = 3; // 11
    }


    template <typename T>
    static void read_value(std::ifstream& input, T& value) {
        input.read(reinterpret_cast<char*>(&value), sizeof(T));

        if (!input) {
            throw std::runtime_error("Failed while reading PT16 H/L table");
        }
    }


    void load_hl(const std::string& path) {
        std::ifstream input(path, std::ios::binary);

        if (!input) {
            throw std::runtime_error("Cannot open PT16 H/L table: " + path);
        }

        char magic[8]{};
        input.read(magic, sizeof(magic));

        const char expected_magic[8] = {'P', 'T', '1', '6', 'H', 'L', '0', '1'};

        if (!input || std::memcmp(magic, expected_magic, sizeof(magic)) != 0) {
            throw std::runtime_error("Invalid PT16 H/L table");
        }

        std::uint64_t entry_count;
        read_value(input, entry_count);

        H_.resize(number_of_buckets + 1);
        L_.resize(entry_count);
        interval_starts_.resize(entry_count);

        input.read(
            reinterpret_cast<char*>(H_.data()),
            static_cast<std::streamsize>(H_.size() * sizeof(std::uint32_t))
        );

        input.read(
            reinterpret_cast<char*>(L_.data()),
            static_cast<std::streamsize>(L_.size() * sizeof(std::uint16_t))
        );

        input.read(
            reinterpret_cast<char*>(interval_starts_.data()),
            static_cast<std::streamsize>(interval_starts_.size() * sizeof(std::uint32_t))
        );

        if (!input) {
            throw std::runtime_error("Failed while loading PT16 H/L table");
        }

        if (H_.back() != L_.size()) {
            throw std::runtime_error("PT16 H directory does not end at m");
        }

        stats_.entries = L_.size();

        stats_.approx_bytes =
            H_.size() * sizeof(std::uint32_t) +
            L_.size() * sizeof(std::uint16_t) +
            interval_starts_.size() * sizeof(std::uint32_t);
    }


    std::uint32_t pack_16mer(const input_type& input, const std::size_t position) const {
        std::uint32_t key = 0;

        for (std::uint32_t j = 0; j < kmer_length; ++j) {
            const std::uint8_t code = alphatab_[static_cast<unsigned char>(input[position + j])];
            key = (key << 2U) | code;
        }

        return key;
    }


    std::uint32_t interval_end(const std::size_t position) const {
        std::size_t sa_end;

        if (position + 1 < interval_starts_.size()) {
            sa_end = static_cast<std::size_t>(interval_starts_[position + 1]) - 1;
        } else {
            sa_end = sa_->size() - 1;
        }

        /*
        build_entries() skips suffixes shorter than 16. Such suffixes can
        occur between two valid PT16 intervals, so remove them from the
        calculated end of the current interval.
        */
        while (
            static_cast<std::size_t>((*sa_)[sa_end]) + kmer_length >
            ref_->size()
        ) {
            --sa_end;
        }

        return static_cast<std::uint32_t>(sa_end);
    }


    LookupResult lookup(const std::uint32_t key) const {
    const std::uint32_t bucket = key >> low_bits;
    const std::uint16_t low = static_cast<std::uint16_t>(key & low_mask);

    const std::uint32_t begin = H_[bucket];
    const std::uint32_t end = H_[bucket + 1];

    if (begin == end) {
        ++stats_.misses;
        return {};
    }

    for (std::uint32_t position = begin; position < end; ++position) {
        if (L_[position] == low) {
            ++stats_.hits;

            const std::uint32_t sa_start = interval_starts_[position];
            const std::uint32_t sa_end = interval_end(position);

            return {true, sa_start, sa_end};
        }
    }

    ++stats_.misses;
    return {};
}
};