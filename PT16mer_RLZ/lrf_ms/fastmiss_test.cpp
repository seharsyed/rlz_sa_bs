// Standalone correctness test for PT16FastMissParser
// (pt16_rlz_v2_fastmiss.hpp) and its MS variants (pt16_fastmiss_ms.hpp).
//
// For every input position with 16+ characters left, the fast-miss lookup
// must return the same answer as PT16RLZParser::lookupKmer over the same
// table (found, count, positions, match_length, and match_position), and
// the sorted- and bucket-chain variants' lengths must equal brute-force
// matching statistics.
//
// Build and run from PT16mer_RLZ/:
//   g++ -std=c++20 -O2 lrf_ms/fastmiss_test.cpp -o fastmiss_test
//   ./fastmiss_test

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../pt16_build_v2.hpp"
#include "../pt16_rlz_v2.hpp"
#include "../pt16_rlz_v2_fastmiss.hpp"
#include "pt16_fastmiss_ms.hpp"

using Symbol = unsigned char;
using SAType = std::uint32_t;

namespace {

int failures = 0;
int cases = 0;

std::vector<SAType> build_suffix_array(const std::vector<Symbol>& text) {
  std::vector<SAType> sa(text.size());
  std::iota(sa.begin(), sa.end(), SAType{0});

  std::sort(sa.begin(), sa.end(), [&](SAType a, SAType b) {
    return std::lexicographical_compare(text.begin() + a, text.end(),
                                        text.begin() + b, text.end());
  });

  return sa;
}

std::vector<Symbol> random_dna(std::mt19937_64& rng, std::size_t n) {
  static const char letters[] = {'A', 'C', 'G', 'T'};
  std::uniform_int_distribution<int> base(0, 3);
  std::vector<Symbol> s(n);
  for (auto& c : s) c = static_cast<Symbol>(letters[base(rng)]);
  return s;
}

// Random input with real reference substrings spliced in.
std::vector<Symbol> spliced_input(std::mt19937_64& rng,
                                  const std::vector<Symbol>& reference,
                                  std::size_t input_size, std::size_t min_len,
                                  std::size_t max_len) {
  auto input = random_dna(rng, input_size);
  std::size_t at = 0;

  while (at < input_size) {
    const std::size_t wanted =
        std::uniform_int_distribution<std::size_t>(min_len, max_len)(rng);
    const std::size_t len =
        std::min({wanted, input_size - at, reference.size()});
    if (len == 0) break;

    const std::size_t p = std::uniform_int_distribution<std::size_t>(
        0, reference.size() - len)(rng);
    std::copy(reference.begin() + p, reference.begin() + p + len,
              input.begin() + at);

    at += len;
    at += std::uniform_int_distribution<std::size_t>(3, 30)(rng);
  }

  return input;
}

std::uint32_t brute_ms_length(const std::vector<Symbol>& reference,
                              const std::vector<Symbol>& input,
                              std::size_t i) {
  std::size_t best = 0;

  for (std::size_t p = 0; p < reference.size(); ++p) {
    std::size_t l = 0;
    while (i + l < input.size() && p + l < reference.size() &&
           input[i + l] == reference[p + l]) {
      ++l;
    }
    best = std::max(best, l);
  }

  return static_cast<std::uint32_t>(best);
}

void check(const std::string& name, const std::vector<Symbol>& reference,
           const std::vector<Symbol>& input) {
  ++cases;

  const std::vector<SAType> sa = build_suffix_array(reference);
  const std::string table_path =
      (fs::temp_directory_path() / "fastmiss_test.v2bin").string();

  fs::remove(table_path);
  build_pt16_table(reference, sa, table_path);

  const PT16RLZParser<Symbol, SAType> v2(reference, sa, table_path);
  const PT16FastMissParser<Symbol, SAType> fast(reference, sa, table_path);

  int mismatches = 0;

  for (std::size_t i = 0; i + KMER_LENGTH <= input.size(); ++i) {
    const KmerLookupResult expected = v2.lookupKmer(input, i);
    const KmerLookupResult actual = fast.lookupKmerByKey(
        encode_16mer(input, static_cast<std::uint32_t>(i)));

    const bool same =
        expected.found == actual.found && expected.count == actual.count &&
        expected.match_length == actual.match_length &&
        expected.match_position == actual.match_position &&
        std::equal(expected.positions.begin(), expected.positions.end(),
                   actual.positions.begin(), actual.positions.end());

    if (!same && mismatches++ < 5) {
      std::cerr << "FAIL [" << name << "] lookup at " << i
                << ": v2 found=" << expected.found
                << " count=" << expected.count
                << " len=" << expected.match_length
                << " pos=" << expected.match_position
                << " | fastmiss found=" << actual.found
                << " count=" << actual.count
                << " len=" << actual.match_length
                << " pos=" << actual.match_position << "\n";
    }
  }

  // Both full variants (sorted and bucketed), end to end, against brute
  // force.
  PT16FastMissSortedChainMS<Symbol, SAType> sorted_chain(
      reference, sa, table_path + ".fastmiss");
  PT16FastMissBucketChainMS<Symbol, SAType> bucket_chain(
      reference, sa, table_path + ".fastmiss");

  for (const MatchingStatistics& ms :
       {sorted_chain.computeMatchingStatistics(input),
        bucket_chain.computeMatchingStatistics(input)})
  for (std::size_t i = 0; i < input.size(); ++i) {
    const std::uint32_t want = brute_ms_length(reference, input, i);
    const auto [pos, len] = ms[i];

    bool genuine = pos + len <= reference.size();
    for (std::uint32_t k = 0; genuine && k < len; ++k) {
      genuine = reference[pos + k] == input[i + k];
    }

    if ((len != want || !genuine) && mismatches++ < 5) {
      std::cerr << "FAIL [" << name << "] ms at " << i << ": brute=" << want
                << " got len=" << len << " pos=" << pos
                << (genuine ? "" : " (not a real match)") << "\n";
    }
  }

  fs::remove(table_path);
  fs::remove(table_path + ".fastmiss");

  if (mismatches > 0) {
    ++failures;
    std::cerr << "  [" << name << "] " << mismatches << " mismatches\n";
  }
}

}  // namespace

int main() {
  std::mt19937_64 rng(20260923);

  // Small references: most buckets are empty, so the precomputed
  // empty-bucket answers carry most misses.
  for (const std::size_t ref_size : {20, 60, 300, 2000}) {
    const auto reference = random_dna(rng, ref_size);
    check("random ref " + std::to_string(ref_size), reference,
          spliced_input(rng, reference, 600, 4, 40));
  }

  // Larger references: most buckets non-empty; misses go through the L
  // search and the bitmap test.
  for (const std::size_t ref_size : {20000, 60000}) {
    const auto reference = random_dna(rng, ref_size);
    check("random ref " + std::to_string(ref_size), reference,
          spliced_input(rng, reference, 1500, 4, 60));
    check("noise vs ref " + std::to_string(ref_size), reference,
          random_dna(rng, 1500));
  }

  // Repetitive reference: many range hits.
  {
    const auto unit = random_dna(rng, 37);
    std::vector<Symbol> reference;
    for (int r = 0; r < 60; ++r) {
      reference.insert(reference.end(), unit.begin(), unit.end());
      if (r % 7 == 0) reference.push_back('A');
    }
    check("repetitive ref", reference,
          spliced_input(rng, reference, 1500, 10, 120));
  }

  // Short suffix under 8 characters (the old bucket-filter bug).
  {
    auto reference = random_dna(rng, 500);
    const std::vector<Symbol> tail = {'G', 'A', 'T', 'T', 'A', 'C', 'A'};
    reference.insert(reference.end(), tail.begin(), tail.end());

    std::vector<Symbol> input = random_dna(rng, 200);
    for (std::size_t at = 0; at + tail.size() < input.size(); at += 23) {
      std::copy(tail.begin(), tail.end(), input.begin() + at);
    }
    check("7-char short suffix", reference, input);
  }

  // Short suffixes of 9..15 characters: the only case that still runs the
  // full short-suffix check (flagged buckets), in both empty and
  // non-empty buckets.
  for (const std::size_t ref_size : {300, 20000}) {
    auto reference = random_dna(rng, ref_size);
    const std::vector<Symbol> tail(reference.end() - 15, reference.end());

    std::vector<Symbol> input = random_dna(rng, 800);
    for (std::size_t at = 0, k = 0; at + 16 < input.size(); at += 40, ++k) {
      // Each tail suffix of length 9..15, followed by random characters.
      const std::size_t length = 9 + (k % 7);
      std::copy(tail.end() - length, tail.end(), input.begin() + at);
    }
    check("long short suffix, ref " + std::to_string(ref_size), reference,
          input);
  }

  std::cerr << "\ncases=" << cases << " failures=" << failures << "\n";
  std::cerr << (failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED") << "\n";
  return failures == 0 ? 0 : 1;
}
