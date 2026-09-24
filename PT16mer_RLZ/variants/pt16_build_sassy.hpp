#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "pt16_sassy_format.hpp"  // PackedShortSuffix, sassy_encode_*
#include "../pt16_utils.hpp"         // PT16 constants and shared helpers

// ---------- PT16 H/L representation ----------

struct SassyScanResult {
  // Lower-level entries in lexicographic order.
  std::vector<std::uint64_t> L;

  // Starting SA position corresponding to each bucket.
  std::vector<std::uint32_t> H_sa;

  // Full relative offsets for entries marked with the 65535 escape value.
  std::vector<std::uint32_t> sampled_sa;

  // Number of PT16 entries belonging to each bucket.
  std::array<std::uint32_t, NUMBER_OF_BUCKETS> count{};
};

/*
Stores one completed PT16 entry: its low 16-bit value, and either the text
position (one occurrence) or a slice of the sampled SA (several occurrences).
The occurrences are the suffixes suffix_array[sa_start, sa_end). This interval
holds only suffixes that start with the 16-mer, so it contains no suffix
shorter than 16.

A range whose count does not fit the entry's 16-bit count field escapes: the
real count is pushed onto sampled_sa as a plain uint32_t, immediately before
its positions (see the "Sassy L entry" section of pt16_utils.hpp).
*/

static void store_entry(const std::vector<std::uint32_t>& suffix_array,
                        SassyScanResult& result, const std::uint32_t key,
                        const std::uint32_t sa_start,
                        const std::uint32_t sa_end) {
  const std::uint32_t bucket = key >> LOW_BITS;
  const std::uint16_t low = static_cast<std::uint16_t>(key & LOW_MASK);

  // The first interval in the bucket gives the bucket's SA start.
  if (result.count[bucket] == 0) {
    result.H_sa[bucket] = static_cast<std::uint32_t>(result.sampled_sa.size());
  }

  if (sa_end - sa_start == 1) {
    result.L.push_back(sassy_encode_singleton(low, suffix_array[sa_start]));
  } else {
    const std::uint32_t offset = static_cast<std::uint32_t>(
        result.sampled_sa.size() - result.H_sa[bucket]);

    const std::uint32_t occurrences = sa_end - sa_start;

    std::uint16_t count_field;

    if (occurrences < SASSY_COUNT_ESCAPE) {
      count_field = static_cast<std::uint16_t>(occurrences);
    } else {
      // The real count does not fit in 16 bits: record it as a plain
      // uint32_t right before the positions themselves.
      count_field = SASSY_COUNT_ESCAPE;
      result.sampled_sa.push_back(occurrences);
    }

    for (std::uint32_t i = sa_start; i < sa_end; ++i) {
      result.sampled_sa.push_back(suffix_array[i]);
    }

    result.L.push_back(sassy_encode_range(low, offset, count_field));
  }

  ++result.count[bucket];
}

/*
Scans the suffix array, groups identical 16-mers into SA intervals, and
stores each interval as a table entry (see store_entry).
*/

static SassyScanResult build_sassy_entries(
    const std::vector<unsigned char>& reference,
    const std::vector<std::uint32_t>& suffix_array) {
  SassyScanResult result;

  // One SA starting position for each H bucket.
  result.H_sa.resize(NUMBER_OF_BUCKETS, 0);

  const std::size_t reserve_size = reference.size() - KMER_LENGTH + 1;

  result.L.reserve(reserve_size);

  bool interval_open = false;

  std::uint32_t current_key = 0;
  std::uint32_t current_start = 0;

  // One past the last suffix of the open interval. It is advanced only by
  // suffixes of the interval itself, so suffixes shorter than 16 that follow
  // it in the suffix array are never counted. They cannot lie inside the
  // interval, because they do not start with its 16-mer.
  std::uint32_t current_end = 0;

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
      current_end = sa_index + 1;
      continue;
    }

    // Same 16-mer: remain inside the current interval.
    if (key == current_key) {
      current_end = sa_index + 1;
      continue;
    }

    // New 16-mer: store the completed interval.
    store_entry(suffix_array, result, current_key, current_start, current_end);

    current_key = key;
    current_start = sa_index;
    current_end = sa_index + 1;
  }

  // Store the final interval.
  if (interval_open) {
    store_entry(suffix_array, result, current_key, current_start, current_end);
  }

  return result;
}

// ---------- Write PT16 H/L representation ----------

static void write_hl_table(
    const std::string& output_path, const SassyScanResult& scan,
    const std::vector<std::uint32_t>& H,
    const std::vector<PackedShortSuffix>& short_suffixes) {
  std::ofstream output(output_path, std::ios::binary);

  if (!output) {
    throw std::runtime_error("Cannot create PT16 H/L/SAS table: " +
                             output_path);
  }

  // Layout of the file, so a reader can size every array up front:
  //   magic                      8 bytes
  //   entry_count                uint64   number of L entries
  //   sampled_count              uint64   number of sampled_sa entries
  //   short_count                uint64   number of short suffix records
  //   H                          (NUMBER_OF_BUCKETS + 1) x uint32
  //   H_sa                       NUMBER_OF_BUCKETS x uint32
  //   L                          entry_count x uint64
  //   sampled_sa                 sampled_count x uint32
  //   short suffixes             short_count x PackedShortSuffix (3 x uint32)
  //
  // The table is self-contained: it holds everything a lookup needs, so the
  // reference and the suffix array are not needed to read it.
  //
  // The magic differs from the plain PT16 table ("PT16HL01"), whose
  // lower-level entries have a different layout, and from the earlier sassy
  // layouts ("PT16SA01" without short suffixes; "PT16SA02" with a 16-bit
  // range offset that can silently overflow on a real genome), so a reader
  // of one format rejects a table of another instead of misreading it.
  const char magic[8] = {'P', 'T', '1', '6', 'S', 'A', '0', '3'};

  const std::uint64_t entry_count = scan.L.size();
  const std::uint64_t sampled_count = scan.sampled_sa.size();
  const std::uint64_t short_count = short_suffixes.size();

  output.write(magic, sizeof(magic));
  write_value(output, entry_count);
  write_value(output, sampled_count);
  write_value(output, short_count);

  // Write the H-to-L directory.
  write_vector(output, H);

  // Write the starting SA position for each bucket.
  write_vector(output, scan.H_sa);

  // Write the low 16-mer value and the associated text position in the
  // reference or pointer to the sampled suffix array.
  write_vector(output, scan.L);

  write_vector(output, scan.sampled_sa);

  // Write the suffixes shorter than 16 that still fill a bucket.
  write_vector(output, short_suffixes);

  if (!output) {
    throw std::runtime_error("Failed while writing PT16 H/L/SAS table.");
  }
}

// ---------- Public PT16 preprocessing function ----------

// Named distinctly from pt16_build.hpp/pt16_build_v2.hpp's build_pt16_table
// (and SassyScanResult from their ScanResult), so a driver that benchmarks
// more than one PT16 format can include this header alongside theirs in the
// same translation unit.
void build_pt16_sassy_table(const std::vector<unsigned char>& reference,
                            const std::vector<std::uint32_t>& suffix_array,
                            const std::string& output_path) {
  const SassyScanResult scan = build_sassy_entries(reference, suffix_array);
  const std::vector<std::uint32_t> H = build_H(scan.count);
  const std::vector<PackedShortSuffix> short_suffixes =
      build_packed_short_suffixes(reference);

  write_hl_table(output_path, scan, H, short_suffixes);
}