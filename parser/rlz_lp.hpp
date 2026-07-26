#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "parser.hpp"

template <typename T1, typename T2>
class RLZLPParser {
public:
    using input_type = std::vector<T1>;
    using reference_type = std::vector<T1>;
    using suffix_array_type = std::vector<T2>;
    using factor_type = std::tuple<std::size_t, std::size_t>;
    using phrase_type = std::tuple<std::size_t, std::size_t, std::size_t>;
    using phrase_vector_type = std::vector<phrase_type>;

    struct CacheInfo {
        std::size_t hits = 0;
        std::size_t misses = 0;
        std::size_t current_size = 0;
        std::size_t table_slots = 0;
        double hit_rate = 0.0;
        double load_factor = 0.0;
        std::size_t approx_bytes = 0;
        bool table_full = false;
    };

private:
    struct LookupKey {
        std::size_t lb;
        std::size_t rb;
        std::size_t offset;
        T1 symbol;

        bool operator==(const LookupKey& other) const noexcept {
            return lb == other.lb &&
                   rb == other.rb &&
                   offset == other.offset &&
                   symbol == other.symbol;
        }
    };

    struct Interval {
        std::size_t new_lb;
        std::size_t new_rb;
    };

    struct Entry {
        bool occupied = false;
        LookupKey key{};
        Interval interval{};
    };

    const reference_type* ref_ = nullptr;
    const suffix_array_type* sa_ = nullptr;

    std::size_t div_p_ = 64;
    std::vector<Entry> table_;

    std::size_t hits_ = 0;
    std::size_t misses_ = 0;
    std::size_t entries_ = 0;
    bool table_full_ = false;

public:
    explicit RLZLPParser(
        const reference_type& ref,
        const suffix_array_type& sa,
        const std::size_t div_p = 64
    )
        : ref_(&ref),
          sa_(&sa),
          div_p_(div_p == 0 ? 64 : div_p) {

        if (ref.empty()) {
            throw std::invalid_argument(
                "RLZLPParser: reference is empty"
            );
        }

        if (sa.empty()) {
            throw std::invalid_argument(
                "RLZLPParser: suffix array is empty"
            );
        }

        initialize_table();
    }

    factor_type computeLZFactorAt(
        const input_type& input,
        const std::size_t input_pos
    ) {
        std::size_t offset = 0;
        std::size_t j = input_pos;
        std::size_t match = 0;

        std::size_t nlb = 0;
        std::size_t nrb = ref_->size() - 1;

        while (j < input.size()) {
            if (nlb == nrb) {
                if ((*ref_)[
                        static_cast<std::size_t>((*sa_)[nlb]) + offset
                    ] != input[j]) {
                    break;
                }
            } else {
                Interval cached_interval{};
                const T1 symbol = input.at(j);

                if (lookup(
                        nlb,
                        nrb,
                        offset,
                        symbol,
                        cached_interval
                    )) {
                    nlb = cached_interval.new_lb;
                    nrb = cached_interval.new_rb;
                } else {
                    const auto opt_lb =
                        ::binarySearchLB<T1, T2>(
                            *ref_,
                            *sa_,
                            static_cast<std::int64_t>(nlb),
                            static_cast<std::int64_t>(nrb),
                            static_cast<std::int64_t>(offset),
                            symbol
                        );

                    if (!opt_lb) {
                        break;
                    }

                    const std::size_t new_lb =
                        static_cast<std::size_t>(*opt_lb);

                    const auto opt_rb =
                        ::binarySearchRB<T1, T2>(
                            *ref_,
                            *sa_,
                            static_cast<std::int64_t>(new_lb),
                            static_cast<std::int64_t>(nrb),
                            static_cast<std::int64_t>(offset),
                            symbol
                        );

                    if (!opt_rb) {
                        break;
                    }

                    const std::size_t new_rb =
                        static_cast<std::size_t>(*opt_rb);

                    insert(
                        nlb,
                        nrb,
                        offset,
                        symbol,
                        new_lb,
                        new_rb
                    );

                    nlb = new_lb;
                    nrb = new_rb;
                }
            }

            match = static_cast<std::size_t>((*sa_)[nlb]);

            ++j;
            ++offset;
        }

        return {match, offset};
    }

    phrase_vector_type lzFactorize(
        const input_type& input
    ) {
        phrase_vector_type spl_vec;
        std::size_t i = 0;

        while (i < input.size()) {
            auto [pos, len] =
                computeLZFactorAt(input, i);

            if (len <= 1) {
                pos = static_cast<std::size_t>(
                    input.at(i)
                );
                len = 1;
            }

            spl_vec.push_back({i, pos, len});
            i += len;
        }

        return spl_vec;
    }

    void clear_cache() {
        initialize_table();

        hits_ = 0;
        misses_ = 0;
    }

    CacheInfo cache_info() const {
        CacheInfo info;

        info.hits = hits_;
        info.misses = misses_;
        info.current_size = entries_;
        info.table_slots = table_.size();
        info.table_full = table_full_;

        const std::size_t total =
            hits_ + misses_;

        info.hit_rate =
            total == 0
                ? 0.0
                : static_cast<double>(hits_) /
                      static_cast<double>(total);

        info.load_factor =
            table_.empty()
                ? 0.0
                : static_cast<double>(entries_) /
                      static_cast<double>(table_.size());

        info.approx_bytes =
            table_.size() * sizeof(Entry);

        return info;
    }

private:
    void initialize_table() {
        const std::size_t slots =
            1 + ((sa_->size() - 1) / div_p_);

        table_.assign(slots, Entry{});
        entries_ = 0;
        table_full_ = false;
    }

    bool lookup(
        const std::size_t lb,
        const std::size_t rb,
        const std::size_t offset,
        const T1 symbol,
        Interval& cached_interval
    ) {
        const LookupKey key{
            lb,
            rb,
            offset,
            symbol
        };

        std::size_t pos = lb / div_p_;

        assert(pos < table_.size());

        for (
            std::size_t scanned = 0;
            scanned < table_.size();
            ++scanned
        ) {
            const Entry& entry = table_[pos];

            if (!entry.occupied) {
                ++misses_;
                return false;
            }

            if (entry.key == key) {
                cached_interval = entry.interval;
                ++hits_;
                return true;
            }

            ++pos;

            if (pos == table_.size()) {
                pos = 0;
            }
        }

        ++misses_;
        return false;
    }

    void insert(
        const std::size_t lb,
        const std::size_t rb,
        const std::size_t offset,
        const T1 symbol,
        const std::size_t new_lb,
        const std::size_t new_rb
    ) {
        if (table_full_) {
            return;
        }

        std::size_t pos = lb / div_p_;

        assert(pos < table_.size());

        for (
            std::size_t scanned = 0;
            scanned < table_.size();
            ++scanned
        ) {
            Entry& entry = table_[pos];

            if (!entry.occupied) {
                entry.occupied = true;

                entry.key = LookupKey{
                    lb,
                    rb,
                    offset,
                    symbol
                };

                entry.interval = Interval{
                    new_lb,
                    new_rb
                };

                ++entries_;

                if (entries_ == table_.size()) {
                    table_full_ = true;
                }

                return;
            }

            ++pos;

            if (pos == table_.size()) {
                pos = 0;
            }
        }

        table_full_ = true;
    }
};