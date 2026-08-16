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

#include "parser.hpp"   // original RLZ parser

/**
 * CachedRLZParser is an add-on wrapper for the RLZ parser.
 *
 * The file binds one reference vector and one suffix array, then mirrors
 * lzFactorize/computeLZFactorAt while caching the expensive interval-refinement step:
 *
 *  (lb, rb, offset, next_symbol) -> (new_lb, new_rb)
 *
 **/

template <typename T1, typename T2>
class CachedRLZParser {
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
        std::size_t bucket_count = 0;
        std::size_t max_bucket_size = 0;
        double hit_rate = 0.0;
        double load_factor = 0.0;
        std::size_t approx_bytes = 0;
    };

private:
    struct LookupKey {
        std::size_t lb_remainder;
        std::size_t rb;
        std::size_t offset;
        T1 symbol;

        bool operator==(const LookupKey& other) const noexcept {
            return lb_remainder == other.lb_remainder &&
                   rb == other.rb &&
                   offset == other.offset &&
                   symbol == other.symbol;
        }
    };

    struct Interval {
        std::size_t new_lb;
        std::size_t new_rb;
    };

    using Entry = std::pair<LookupKey, Interval>;
    using Bucket = std::vector<Entry>;

    const reference_type* ref_ = nullptr;
    const suffix_array_type* sa_ = nullptr;

    // Used to map lb to a bucket:
    // bucket_index = lb / bucket_divisor_.
    // This is not the number of entries allowed in a bucket.
    std::size_t bucket_divisor_ = 256;

    std::vector<Bucket> table_;

    std::size_t hits_ = 0;
    std::size_t misses_ = 0;
    std::size_t entries_ = 0;

    // Per-file trace of bucket indices where cache hits occur.
    // This is cleared before each input file, but the cache itself stays warm.
    std::vector<std::size_t> bucket_hit_sequence_;

public:
    /**
     * Bind this cache to exactly one (reference, suffix-array) pair.
     */
    explicit CachedRLZParser(
        const reference_type& ref,
        const suffix_array_type& sa,
        std::size_t bucket_divisor = 256
    ) : ref_(&ref),
        sa_(&sa),
        bucket_divisor_(bucket_divisor == 0 ? 256 : bucket_divisor) {
        if (ref.empty()) {
            throw std::invalid_argument("CachedRLZParser: reference is empty");
        }
        if (sa.empty()) {
            throw std::invalid_argument("CachedRLZParser: suffix array is empty");
        }
        rebuild_table();
    }

    /**
     * File loading.
     */
    template <typename T>
    static std::vector<T> read_file(const char* filename) {
        return ::read_file<T>(filename);
    }

    /**
     * binarySearchLB: The useful cached unit for RLZ is the combined interval
     * transition, because computeLZFactorAt always needs both LB and RB before
     * continuing.
     */
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

    /**
     * binarySearchRB is not cached alone for the same reason as binarySearchLB.
     * The wrapper caches (lb, rb, offset, c) -> (new_lb, new_rb).
     */
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

    /**
     * Cached equivalent of RLZ computeLZFactorAt(input, ref, sa, pos).
     * It caches only interval refinements; failed refinements are not cached.
     */
    factor_type computeLZFactorAt(
        const input_type& input,
        const reference_type& ref,
        const suffix_array_type& sa,
        const std::size_t input_pos
    ) {
        assert_bound_pair(ref, sa);
        return computeLZFactorAt(input, input_pos);
    }

    /**
     * Convenience overload using the reference and suffix array bound at construction.
     */
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
                // Not cached: there is no interval refinement, only one suffix check.
                if (safe_ref_symbol(*ref_, *sa_, nlb, offset) != input[j]) {
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

    /**
     * Cached equivalent of lzFactorize(input, ref, sa).
     * The whole factor list is not cached; only the internal suffix-array
     * interval transitions are cached.
     */
    phrase_vector_type lzFactorize(
        const input_type& input,
        const reference_type& ref,
        const suffix_array_type& sa
    ) {
        assert_bound_pair(ref, sa);
        return lzFactorize(input);
    }

    /**
     * Convenience overload using the reference and suffix array bound at construction.
     */
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

    /**
     * Clear cache entries and statistics. Use before a cold-cache timing run.
     */
    void clear_cache() {
        rebuild_table();
        hits_ = 0;
        misses_ = 0;
        entries_ = 0;
        bucket_hit_sequence_.clear();
    }

    /**
     * Clear only the per-file bucket-hit sequence.
     * This must not clear the cache because experiments use a warm cache across files.
     */
    void clear_bucket_hit_sequence() {
        bucket_hit_sequence_.clear();
    }

    /**
     * Return the sequence of bucket indices hit while parsing the current file.
     */
    const std::vector<std::size_t>& bucket_hit_sequence() const {
        return bucket_hit_sequence_;
    }

    /**
     * Report hit count, miss count, current dynamic cache size, and basic table stats.
     */
    CacheInfo cache_info() const {
        CacheInfo info;

        info.hits = hits_;
        info.misses = misses_;
        info.current_size = entries_;
        info.bucket_count = table_.size();
        info.max_bucket_size = max_bucket_size();

        const std::size_t total = hits_ + misses_;
        info.hit_rate =
            total == 0
                ? 0.0
                : static_cast<double>(hits_) / static_cast<double>(total);

        info.load_factor =
            table_.empty()
                ? 0.0
                : static_cast<double>(entries_) / static_cast<double>(table_.size());

        info.approx_bytes =
            table_.size() * sizeof(Bucket) + entries_ * sizeof(Entry);

        return info;
    }

private:
    void rebuild_table() {
        const std::size_t buckets = std::max<std::size_t>(
            1,
            (ref_->size() + bucket_divisor_ - 1) / bucket_divisor_
        );

        table_.clear();
        table_.resize(buckets);
    }

    void assert_bound_pair(const reference_type& ref,
                           const suffix_array_type& sa) const {
        if (&ref != ref_ || &sa != sa_) {
            throw std::invalid_argument(
                "CachedRLZParser: this cache is bound to a different reference/suffix-array pair"
            );
        }
    }

    std::size_t bucket_index(const std::size_t lb) const {
        return lb / bucket_divisor_;
    }

    LookupKey make_key(const std::size_t lb,
                       const std::size_t rb,
                       const std::size_t offset,
                       const T1 symbol) const {
        return LookupKey{lb % bucket_divisor_, rb, offset, symbol};
    }

    bool lookup(const std::size_t lb,
                const std::size_t rb,
                const std::size_t offset,
                const T1 symbol,
                Interval& cached_interval) {
        const std::size_t b = bucket_index(lb);
        assert(b < table_.size());

        const LookupKey key = make_key(lb, rb, offset, symbol);

        for (const auto& [entry_key, entry_interval] : table_[b]) {
            if (entry_key == key) {
                cached_interval = entry_interval;
                ++hits_;

                // Record the bucket index only when the cache actually hits.
                bucket_hit_sequence_.push_back(b);

                return true;
            }
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
        const std::size_t b = bucket_index(lb);
        assert(b < table_.size());

        table_[b].emplace_back(
            make_key(lb, rb, offset, symbol),
            Interval{new_lb, new_rb}
        );

        ++entries_;
    }

    std::size_t max_bucket_size() const {
        std::size_t m = 0;

        for (const auto& bucket : table_) {
            m = std::max(m, bucket.size());
        }

        return m;
    }
};