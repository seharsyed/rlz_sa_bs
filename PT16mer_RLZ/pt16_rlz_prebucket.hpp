#pragma once

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
        std::size_t singleton_hits = 0;
        std::size_t range_hits = 0;
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

    std::vector<std::uint32_t> H_;
    std::vector<std::uint16_t> L_;
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


    factor_type computeLZFactorAt(
        const input_type& input,
        const std::size_t input_pos,
        const std::vector<LookupResult>& lookup_results
    ) {
        // Fewer than 16 characters remain, so use ordinary RLZ.
        if (input.size() - input_pos < kmer_length) {
            return ::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);
        }

        // Use the PT16 result already computed during pre-bucketing.
        const LookupResult& result = lookup_results[input_pos];

        if (!result.found) {
            ++stats_.misses;
            return ::computeLZFactorAt<T1, T2>(input, *ref_, *sa_, input_pos);
        }

        ++stats_.hits;

        if (result.sa_start == result.sa_end) {
            ++stats_.singleton_hits;
        } else {
            ++stats_.range_hits;
        }

        std::size_t offset = kmer_length;
        std::size_t j = input_pos + kmer_length;
        std::size_t nlb = result.sa_start;
        std::size_t nrb = result.sa_end;

        // Range intervals are narrowed from character 17 onward.
        while (nlb < nrb && j < input.size()) {
            const auto lb = ::binarySearchLB<T1, T2>(
                *ref_,
                *sa_,
                nlb,
                nrb,
                offset,
                input[j]
            );

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

        // A singleton interval is extended by direct comparison.
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
        // Group positions of S by the upper 16 bits of their packed 16-mer.
        std::array<std::vector<std::size_t>, number_of_buckets> position_buckets;

        // Keep the packed key so the 16-mer is not packed again during
        // the grouped PT16 access.
        std::vector<std::uint32_t> keys(input.size());

        for (std::size_t i = 0; i + kmer_length <= input.size(); ++i) {
            const std::uint32_t key = pack_16mer(input, i);
            const std::uint32_t bucket = key >> low_bits;

            keys[i] = key;
            position_buckets[bucket].push_back(i);
        }

        // Store the PT16 result for every valid position of S.
        std::vector<LookupResult> lookup_results(input.size());

        // Process positions belonging to the same H/L bucket together.
        for (std::uint32_t bucket = 0; bucket < number_of_buckets; ++bucket) {
            for (const std::size_t i : position_buckets[bucket]) {
                lookup_results[i] = prebucketing_s(keys[i]);
            }
        }

        // Greedy RLZ parsing remains left-to-right.
        phrase_vector_type spl_vec;
        std::size_t i = 0;

        while (i < input.size()) {
            auto [pos, len] = computeLZFactorAt(input, i, lookup_results);

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
    void initialise_alphatab() {
        alphatab_[static_cast<unsigned char>('A')] = 0;
        alphatab_[static_cast<unsigned char>('C')] = 1;
        alphatab_[static_cast<unsigned char>('G')] = 2;
        alphatab_[static_cast<unsigned char>('T')] = 3;
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

        const char expected_magic[8] = {
            'P', 'T', '1', '6', 'H', 'L', '0', '1'
        };

        if (
            !input ||
            std::memcmp(
                magic,
                expected_magic,
                sizeof(magic)
            ) != 0
        ) {
            throw std::runtime_error("Invalid PT16 H/L table");
        }

        std::uint64_t entry_count;
        read_value(input, entry_count);

        H_.resize(number_of_buckets + 1);
        L_.resize(entry_count);
        interval_starts_.resize(entry_count);

        input.read(
            reinterpret_cast<char*>(H_.data()),
            static_cast<std::streamsize>(
                H_.size() * sizeof(std::uint32_t)
            )
        );

        input.read(
            reinterpret_cast<char*>(L_.data()),
            static_cast<std::streamsize>(
                L_.size() * sizeof(std::uint16_t)
            )
        );

        input.read(
            reinterpret_cast<char*>(interval_starts_.data()),
            static_cast<std::streamsize>(
                interval_starts_.size() *
                sizeof(std::uint32_t)
            )
        );

        if (!input) {
            throw std::runtime_error(
                "Failed while loading PT16 H/L table"
            );
        }

        if (H_.back() != L_.size()) {
            throw std::runtime_error(
                "PT16 H directory does not end at m"
            );
        }

        stats_.entries = L_.size();

        stats_.approx_bytes =
            H_.size() * sizeof(std::uint32_t) +
            L_.size() * sizeof(std::uint16_t) +
            interval_starts_.size() *
                sizeof(std::uint32_t);
    }


    std::uint32_t pack_16mer(
        const input_type& input,
        const std::size_t position
    ) const {
        std::uint32_t key = 0;

        for (
            std::uint32_t j = 0;
            j < kmer_length;
            ++j
        ) {
            const std::uint8_t code =
                alphatab_[
                    static_cast<unsigned char>(
                        input[position + j]
                    )
                ];

            key = (key << 2U) | code;
        }

        return key;
    }


    std::uint32_t interval_end(
        const std::size_t position
    ) const {
        std::size_t sa_end;

        if (position + 1 < interval_starts_.size()) {
            sa_end =
                static_cast<std::size_t>(
                    interval_starts_[position + 1]
                ) - 1;
        } else {
            sa_end = sa_->size() - 1;
        }

        // Suffixes shorter than 16 are not PT16 entries.
        while (
            static_cast<std::size_t>(
                (*sa_)[sa_end]
            ) + kmer_length >
            ref_->size()
        ) {
            --sa_end;
        }

        return static_cast<std::uint32_t>(sa_end);
    }


    // Existing PT16 lookup used by the original parser.
    LookupResult lookup(
        const std::uint32_t key
    ) const {
        const std::uint32_t bucket =
            key >> low_bits;

        const std::uint16_t low =
            static_cast<std::uint16_t>(
                key & low_mask
            );

        const std::uint32_t begin =
            H_[bucket];

        const std::uint32_t end =
            H_[bucket + 1];

        if (begin == end) {
            ++stats_.misses;
            return {};
        }

        for (
            std::uint32_t position = begin;
            position < end;
            ++position
        ) {
            if (L_[position] == low) {
                ++stats_.hits;

                const std::uint32_t sa_start =
                    interval_starts_[position];

                const std::uint32_t sa_end =
                    interval_end(position);

                if (sa_start == sa_end) {
                    ++stats_.singleton_hits;
                } else {
                    ++stats_.range_hits;
                }

                return {
                    true,
                    sa_start,
                    sa_end
                };
            }
        }

        ++stats_.misses;
        return {};
    }


    // PT16 search used while pre-bucketing S.
    // Statistics are counted later only at actual RLZ phrase starts.
    LookupResult prebucketing_s(
        const std::uint32_t key
    ) const {
        const std::uint32_t bucket =
            key >> low_bits;

        const std::uint16_t low =
            static_cast<std::uint16_t>(
                key & low_mask
            );

        const std::uint32_t begin =
            H_[bucket];

        const std::uint32_t end =
            H_[bucket + 1];

        if (begin == end) {
            return {};
        }

        for (
            std::uint32_t position = begin;
            position < end;
            ++position
        ) {
            if (L_[position] == low) {
                const std::uint32_t sa_start =
                    interval_starts_[position];

                const std::uint32_t sa_end =
                    interval_end(position);

                return {
                    true,
                    sa_start,
                    sa_end
                };
            }
        }

        return {};
    }
};