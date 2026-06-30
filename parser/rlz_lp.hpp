#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
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
        std::size_t max_probe_cluster = 0;
        double hit_rate = 0.0;
        double load_factor = 0.0;
        std::size_t approx_bytes = 0;
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

    struct Slot {
        bool occupied = false;
        LookupKey key{};
        Interval interval{};
    };

    const reference_type* ref_ = nullptr;
    const suffix_array_type* sa_ = nullptr;

    std::size_t div_p_ = 64;
    std::vector<Slot> table_;

    std::size_t hits_ = 0;
    std::size_t misses_ = 0;
    std::size_t entries_ = 0;

    static constexpr double max_load_factor_ = 0.70;

public:
    explicit RLZLPParser(
        const reference_type& ref,
        const suffix_array_type& sa,
        std::size_t div_p = 64
    ) : ref_(&ref),
        sa_(&sa),
        div_p_(div_p == 0 ? 64 : div_p) {
        if (ref.empty()) {
            throw std::invalid_argument("RLZLPParser: reference is empty");
        }
        if (sa.empty()) {
            throw std::invalid_argument("RLZLPParser: suffix array is empty");
        }
        rebuild_table();
    }

    template <typename T>
    static std::vector<T> read_file(const char* filename) {
        return ::read_file<T>(filename);
    }

    std::optional<std::int64_t> binarySearchLB(
        const reference_type& ref,
        const suffix_array_type& sa,
        const std::int64_t lo,
        const std::int64_t hi,
        const std::int64_t offset,
        const T1 c
    ) const {
        return ::binarySearchLB<T1, T2>(ref, sa, lo, hi, offset, c);
    }

    std::optional<std::int64_t> binarySearchRB(
        const reference_type& ref,
        const suffix_array_type& sa,
        const std::int64_t lo,
        const std::int64_t hi,
        const std::int64_t offset,
        const T1 c
    ) const {
        return ::binarySearchRB<T1, T2>(ref, sa, lo, hi, offset, c);
    }

    factor_type computeLZFactorAt(
        const input_type& input,
        const reference_type& ref,
        const suffix_array_type& sa,
        const std::size_t input_pos
    ) {
        assert_bound_pair(ref, sa);
        return computeLZFactorAt(input, input_pos);
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
                if ((*ref_)[static_cast<std::size_t>((*sa_)[nlb]) + offset] != input[j]) {
                    break;
                }
            } else {
                Interval cached_interval{};
                const T1 c = input.at(j);

                if (lookup(nlb, nrb, offset, c, cached_interval)) {
                    nlb = cached_interval.new_lb;
                    nrb = cached_interval.new_rb;
                } else {
                    const auto opt_lb = ::binarySearchLB<T1, T2>(
                        *ref_, *sa_,
                        static_cast<std::int64_t>(nlb),
                        static_cast<std::int64_t>(nrb),
                        static_cast<std::int64_t>(offset),
                        c
                    );

                    if (!opt_lb) {
                        break;
                    }

                    const std::size_t new_lb = static_cast<std::size_t>(*opt_lb);

                    const auto opt_rb = ::binarySearchRB<T1, T2>(
                        *ref_, *sa_,
                        static_cast<std::int64_t>(new_lb),
                        static_cast<std::int64_t>(nrb),
                        static_cast<std::int64_t>(offset),
                        c
                    );

                    if (!opt_rb) {
                        break;
                    }

                    const std::size_t new_rb = static_cast<std::size_t>(*opt_rb);

                    if (new_lb > new_rb) {
                        break;
                    }

                    insert(nlb, nrb, offset, c, new_lb, new_rb);

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
        const input_type& input,
        const reference_type& ref,
        const suffix_array_type& sa
    ) {
        assert_bound_pair(ref, sa);
        return lzFactorize(input);
    }

    phrase_vector_type lzFactorize(const input_type& input) {
        phrase_vector_type spl_vec;
        std::size_t i = 0;

        while (i < input.size()) {
            auto [pos, len] = computeLZFactorAt(input, i);

            if (len <= 1) {
                pos = static_cast<std::size_t>(input.at(i));
                len = 1;
            }

            spl_vec.push_back({i, pos, len});
            i += len;
        }

        return spl_vec;
    }

    void clear_cache() {
        rebuild_table();
        hits_ = 0;
        misses_ = 0;
        entries_ = 0;
    }

    CacheInfo cache_info() const {
        CacheInfo info;

        info.hits = hits_;
        info.misses = misses_;
        info.current_size = entries_;
        info.table_slots = table_.size();
        info.max_probe_cluster = max_probe_cluster();

        const std::size_t total = hits_ + misses_;
        info.hit_rate = total == 0
            ? 0.0
            : static_cast<double>(hits_) / static_cast<double>(total);

        info.load_factor = table_.empty()
            ? 0.0
            : static_cast<double>(entries_) / static_cast<double>(table_.size());

        info.approx_bytes = table_.size() * sizeof(Slot);

        return info;
    }

private:
    void rebuild_table() {
        const std::size_t slots = std::max<std::size_t>(
            1,
            sa_->size() / div_p_
        );

        table_.assign(slots, Slot{});
    }

    void assert_bound_pair(const reference_type& ref,
                           const suffix_array_type& sa) const {
        if (&ref != ref_ || &sa != sa_) {
            throw std::invalid_argument(
                "RLZLPParser: this cache is bound to a different reference/suffix-array pair"
            );
        }
    }

    static std::size_t mix_hash(std::size_t x) noexcept {
        x += static_cast<std::size_t>(0x9e3779b97f4a7c15ull);
        x = (x ^ (x >> 30)) * static_cast<std::size_t>(0xbf58476d1ce4e5b9ull);
        x = (x ^ (x >> 27)) * static_cast<std::size_t>(0x94d049bb133111ebull);
        return x ^ (x >> 31);
    }

    std::size_t key_hash(const LookupKey& key) const noexcept {
        std::size_t h = mix_hash(key.lb / div_p_);
        h ^= mix_hash(key.lb + static_cast<std::size_t>(0x9e3779b97f4a7c15ull) + (h << 6) + (h >> 2));
        h ^= mix_hash(key.rb + static_cast<std::size_t>(0x9e3779b97f4a7c15ull) + (h << 6) + (h >> 2));
        h ^= mix_hash(key.offset + static_cast<std::size_t>(0x9e3779b97f4a7c15ull) + (h << 6) + (h >> 2));
        h ^= mix_hash(static_cast<std::size_t>(key.symbol) + static_cast<std::size_t>(0x9e3779b97f4a7c15ull) + (h << 6) + (h >> 2));
        return h;
    }

    std::size_t start_slot(const LookupKey& key) const {
        return key_hash(key) % table_.size();
    }

    LookupKey make_key(const std::size_t lb,
                       const std::size_t rb,
                       const std::size_t offset,
                       const T1 symbol) const {
        return LookupKey{lb, rb, offset, symbol};
    }

    bool lookup(const std::size_t lb,
                const std::size_t rb,
                const std::size_t offset,
                const T1 symbol,
                Interval& cached_interval) {
        const LookupKey key = make_key(lb, rb, offset, symbol);
        std::size_t pos = start_slot(key);

        for (std::size_t probe = 0; probe < table_.size(); ++probe) {
            const Slot& slot = table_[pos];

            if (!slot.occupied) {
                ++misses_;
                return false;
            }

            if (slot.key == key) {
                cached_interval = slot.interval;
                ++hits_;
                return true;
            }

            pos = (pos + 1 == table_.size()) ? 0 : pos + 1;
        }

        ++misses_;
        return false;
    }

    void insert(const std::size_t lb,
                const std::size_t rb,
                const std::size_t offset,
                const T1 symbol,
                const std::size_t new_lb,
                const std::size_t new_rb) {
        const LookupKey key = make_key(lb, rb, offset, symbol);
        maybe_resize();

        std::size_t pos = start_slot(key);

        for (std::size_t probe = 0; probe < table_.size(); ++probe) {
            Slot& slot = table_[pos];

            if (!slot.occupied) {
                slot.occupied = true;
                slot.key = key;
                slot.interval = Interval{new_lb, new_rb};
                ++entries_;
                return;
            }

            if (slot.key == key) {
                slot.interval = Interval{new_lb, new_rb};
                return;
            }

            pos = (pos + 1 == table_.size()) ? 0 : pos + 1;
        }

        throw std::runtime_error("RLZLPParser: linear probing table is full after resize");
    }

    void maybe_resize() {
        if (table_.empty()) {
            table_.assign(1, Slot{});
            return;
        }

        const double next_load =
            static_cast<double>(entries_ + 1) / static_cast<double>(table_.size());

        if (next_load <= max_load_factor_) {
            return;
        }

        std::vector<Slot> old_table = std::move(table_);
        table_.assign(old_table.size() * 2, Slot{});
        entries_ = 0;

        for (const Slot& slot : old_table) {
            if (slot.occupied) {
                reinsert(slot.key, slot.interval);
            }
        }
    }

    void reinsert(const LookupKey& key, const Interval& interval) {
        std::size_t pos = start_slot(key);

        for (std::size_t probe = 0; probe < table_.size(); ++probe) {
            Slot& slot = table_[pos];

            if (!slot.occupied) {
                slot.occupied = true;
                slot.key = key;
                slot.interval = interval;
                ++entries_;
                return;
            }

            pos = (pos + 1 == table_.size()) ? 0 : pos + 1;
        }

        throw std::runtime_error("RLZLPParser: reinsert failed during resize");
    }

    std::size_t max_probe_cluster() const {
        std::size_t best = 0;
        std::size_t run = 0;

        for (const Slot& slot : table_) {
            if (slot.occupied) {
                ++run;
                best = std::max(best, run);
            } else {
                run = 0;
            }
        }

        return best;
    }
};