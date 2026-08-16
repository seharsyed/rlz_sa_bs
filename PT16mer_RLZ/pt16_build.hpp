#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr std::uint32_t KMER_LENGTH = 16;
constexpr std::uint32_t BUCKET_SIZE = 65536;
constexpr std::uint32_t LOW_BITS = 16;
constexpr std::uint32_t LOW_MASK = BUCKET_SIZE - 1;
constexpr std::uint32_t NUMBER_OF_BUCKETS = 65536;


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

struct ScanResult {
    std::vector<std::uint16_t> L;
    std::vector<std::uint32_t> interval_starts;
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

static void store_entry(ScanResult& result, const std::uint32_t key, const std::uint32_t sa_start) {
    result.L.push_back(static_cast<std::uint16_t>(key & LOW_MASK));
    result.interval_starts.push_back(sa_start);

    const std::uint32_t bucket = key >> LOW_BITS;
    ++result.count[bucket];
}


/*
Scans the suffix array, groups identical 16-mers into SA intervals,
and stores only the starting SA position of each interval.
*/

static ScanResult build_entries(const std::vector<unsigned char>& reference, const std::vector<std::uint32_t>& suffix_array) {
    ScanResult result;

    const std::size_t reserve_size = reference.size() - KMER_LENGTH + 1;

    result.L.reserve(reserve_size);
    result.interval_starts.reserve(reserve_size);

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


// Builds H from the number of PT16 entries assigned to each bucket.

static std::vector<std::uint32_t> build_H(const ScanResult& scan) {
    std::vector<std::uint32_t> H(static_cast<std::size_t>(NUMBER_OF_BUCKETS) + 1, 0);

    H[0] = 0;

    for (std::uint32_t bucket = 0; bucket < NUMBER_OF_BUCKETS; ++bucket) {
        H[bucket + 1] = H[bucket] + scan.count[bucket];
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


static void write_hl_table(const std::string& output_path, const ScanResult& scan, const std::vector<std::uint32_t>& H) {
    std::ofstream output(output_path, std::ios::binary);

    if (!output) {
        throw std::runtime_error("Cannot create PT16 H/L table: " + output_path);
    }

    const char magic[8] = {'P', 'T', '1', '6', 'H', 'L', '0', '1'};
    const std::uint64_t entry_count = scan.L.size();

    output.write(magic, sizeof(magic));
    write_value(output, entry_count);

    write_vector(output, H);
    write_vector(output, scan.L);
    write_vector(output, scan.interval_starts);

    if (!output) {
        throw std::runtime_error("Failed while writing PT16 H/L table");
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