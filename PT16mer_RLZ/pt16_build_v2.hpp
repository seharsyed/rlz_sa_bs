#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <bit>
#include <unordered_map>


constexpr std::uint32_t KMER_LENGTH = 16;
constexpr std::uint32_t BUCKET_SIZE = 65536;
constexpr std::uint32_t LOW_BITS = 16;
constexpr std::uint32_t LOW_MASK = BUCKET_SIZE - 1;
constexpr std::uint32_t NUMBER_OF_BUCKETS = 65536;
constexpr std::uint32_t EMPTY_BUCKET_FLAG = 1U << 31;
constexpr std::uint16_t LARGE_OFFSET_FLAG = 65535;

// ---------- Build alphatab once ----------

static std::array<std::uint8_t, 256> build_alphatab() {
    std::array<std::uint8_t, 256> alphatab{};

    alphatab[static_cast<unsigned char>('A')] = 0; // 00
    alphatab[static_cast<unsigned char>('C')] = 1; // 01
    alphatab[static_cast<unsigned char>('G')] = 2; // 10
    alphatab[static_cast<unsigned char>('T')] = 3; // 11

    return alphatab;
}

static const std::array<std::uint8_t, 256> alphatab = build_alphatab();


// ---------- PT16 H/L representation ----------

// Each lower-level entry stores the low 16-mer value and its relative SA offset.
struct LowerEntry {
    std::uint16_t low;
    std::uint16_t sa_offset;
};

struct ScanResult {
    // Lower-level entries in lexicographic order.
    std::vector<LowerEntry> L;

    // Starting SA position corresponding to each bucket.
    std::vector<std::uint32_t> H_sa;

    // Full relative offsets for entries marked with the 65535 escape value.
    std::unordered_map<std::uint32_t, std::uint32_t> large_offsets;

    // Number of PT16 entries belonging to each bucket.
    std::array<std::uint32_t, NUMBER_OF_BUCKETS> count{};
};


// Packs one 16-mer from the reference into a 32-bit key.

static std::uint32_t encode_16mer(const std::vector<unsigned char>& reference, const std::uint32_t position) {
    std::uint32_t key = 0;

    for (std::uint32_t j = 0; j < KMER_LENGTH; ++j) {
        const std::uint8_t code = alphatab[static_cast<unsigned char>(reference[position + j])];
        key = (key << 2U) | code;
    }

    return key;
}


/*
Stores one completed PT16 entry: its low 16-bit value,
SA interval starting position, and H bucket count.
*/

static void store_entry(
    ScanResult& result,
    const std::uint32_t key,
    const std::uint32_t sa_start
) {
    const std::uint32_t bucket = key >> LOW_BITS;
    const std::uint16_t low =
        static_cast<std::uint16_t>(key & LOW_MASK);

    // The first interval in the bucket gives the bucket's SA start.
    if (result.count[bucket] == 0) {
        result.H_sa[bucket] = sa_start;
    }

    const std::uint32_t relative_offset =
        sa_start - result.H_sa[bucket];

    const std::uint32_t position =
        static_cast<std::uint32_t>(result.L.size());

    // Store the offset directly when it fits in the 16-bit field.
    if (relative_offset < LARGE_OFFSET_FLAG) {
        result.L.push_back({
            low,
            static_cast<std::uint16_t>(relative_offset)
        });
    } else {
        // 65535 indicates that the full offset is stored separately.
        result.L.push_back({
            low,
            LARGE_OFFSET_FLAG
        });

        result.large_offsets.emplace(
            position,
            relative_offset
        );
    }

    ++result.count[bucket];
}

/*
Scans the suffix array, groups identical 16-mers into SA intervals,
and stores only the starting SA position of each interval.
*/

static ScanResult build_entries(const std::vector<unsigned char>& reference, const std::vector<std::uint32_t>& suffix_array) {
    ScanResult result;

    // One SA starting position for each H bucket.
    result.H_sa.resize(NUMBER_OF_BUCKETS, 0);

    const std::size_t reserve_size = reference.size() - KMER_LENGTH + 1;

    result.L.reserve(reserve_size);
    

    bool interval_open = false;

    std::uint32_t current_key = 0;
    std::uint32_t current_start = 0;

    for (std::uint32_t sa_index = 0; sa_index < suffix_array.size(); ++sa_index) {
        const std::uint32_t position = suffix_array[sa_index];

        // Suffixes shorter than 16 symbols cannot form a PT16 entry.
        if (static_cast<std::size_t>(position) + KMER_LENGTH > reference.size()) {
            continue;
        }

        const std::uint32_t key = encode_16mer(reference, position);

        // Start the first interval.
        if (!interval_open) {
            interval_open = true;
            current_key = key;
            current_start = sa_index;
            continue;
        }

        // Same 16-mer: remain inside the current interval.
        if (key == current_key) {
            continue;
        }

        // New 16-mer: store the starting position of the completed interval.
        store_entry(result, current_key, current_start);

        current_key = key;
        current_start = sa_index;
    }

    // Store the final interval.
    if (interval_open) {
        store_entry(result, current_key, current_start);
    }

    return result;
}


// Number of common leading bits between two 16-bit bucket prefixes.

static std::uint32_t lcp_bits_16(
    const std::uint32_t first,
    const std::uint32_t second
) {
    const std::uint32_t difference = (first ^ second) << 16U;
    return std::countl_zero(difference);
}


/*
Builds the normal H prefix sums, then marks each empty bucket.
For an empty bucket b:
    MSB = 1
    lower 31 bits = non-empty bucket x with maximum LCP with b.
*/

static std::vector<std::uint32_t> build_H(const ScanResult& scan) {
    std::vector<std::uint32_t> H(static_cast<std::size_t>(NUMBER_OF_BUCKETS) + 1, 0);
    std::vector<std::uint32_t> non_empty_buckets;

    H[0] = 0;

    // First construct the ordinary H directory.
    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
        H[bucket + 1] = H[bucket] + scan.count[bucket];

        if (scan.count[bucket] != 0) {
            non_empty_buckets.push_back(bucket);
        }
    }

    // Encode every empty bucket.
    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
        if (scan.count[bucket] != 0) {
            continue;
        }

        // The best LCP candidate must be the nearest non-empty
        // bucket on the left or right in prefix order.
        const auto next = std::lower_bound(
            non_empty_buckets.begin(),
            non_empty_buckets.end(),
            bucket
        );

        std::uint32_t x;

        if (next == non_empty_buckets.begin()) {
            x = *next;
        } else if (next == non_empty_buckets.end()) {
            x = non_empty_buckets.back();
        } else {
            const std::uint32_t left = *(next - 1);
            const std::uint32_t right = *next;

            if (lcp_bits_16(bucket, left) >= lcp_bits_16(bucket, right)) {
                x = left;
            } else {
                x = right;
            }
        }

        H[bucket] = EMPTY_BUCKET_FLAG | x;
    }

    return H;
}


// ---------- Write PT16 H/L representation ----------

template <typename T>
static void write_value(std::ofstream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}


template <typename T>
static void write_vector(std::ofstream& output, const std::vector<T>& values) {
    if (!values.empty()) {
        output.write(
            reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T))
        );
    }
}

static void write_hl_table(
    const std::string& output_path,
    const ScanResult& scan,
    const std::vector<std::uint32_t>& H
) {
    std::ofstream output(output_path, std::ios::binary);

    if (!output) {
        throw std::runtime_error(
            "Cannot create PT16 H/L table: " + output_path
        );
    }

    const char magic[8] = {
        'P', 'T', '1', '6', 'H', 'L', '0', '1'
    };

    const std::uint64_t entry_count = scan.L.size();

    output.write(magic, sizeof(magic));
    write_value(output, entry_count);

    // Write the H-to-L directory.
    write_vector(output, H);

    // Write the starting SA position for each bucket.
    write_vector(output, scan.H_sa);

    // Write the low 16-mer value and relative SA offset pairs.
    write_vector(output, scan.L);

    // Write full offsets for entries marked with the 65535 escape value.
    for (
        std::uint32_t position = 0;
        position < scan.L.size();
        ++position
    ) {
        if (scan.L[position].sa_offset == LARGE_OFFSET_FLAG) {
            write_value(
                output,
                scan.large_offsets.at(position)
            );
        }
    }

    if (!output) {
        throw std::runtime_error(
            "Failed while writing PT16 H/L table."
        );
    }
}


// ---------- Public PT16 preprocessing function ----------

void build_pt16_table(
    const std::vector<unsigned char>& reference,
    const std::vector<std::uint32_t>& suffix_array,
    const std::string& output_path
) {
    const ScanResult scan = build_entries(reference, suffix_array);
    const std::vector<std::uint32_t> H = build_H(scan);

    write_hl_table(output_path, scan, H);
}