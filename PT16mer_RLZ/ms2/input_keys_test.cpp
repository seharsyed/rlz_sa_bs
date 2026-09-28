// Tests ms2::prepare_keys, bucket_order and sorted_order against brute
// force.
//
//   g++ -std=c++20 -O2 ms2/input_keys_test.cpp -o input_keys_test
//   ./input_keys_test

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "input_keys.hpp"

namespace {

int failures = 0;

void fail(const std::string& name, const std::string& what) {
  if (failures++ < 10) {
    std::cerr << "FAIL " << name << ": " << what << '\n';
  }
}

// Brute force, one position at a time.
void expected(const std::vector<unsigned char>& input,
              std::vector<std::uint64_t>& keys,
              std::vector<std::pair<std::uint32_t, char>>& run_positions,
              std::vector<ms2::ShortQuery>& short_queries,
              std::size_t& separators) {
  for (std::size_t p = 0; p < input.size(); ++p) {
    std::uint32_t length = 0;
    std::uint32_t packed = 0;

    while (length < KMER_LENGTH && p + length < input.size() &&
           is_acgt(input[p + length])) {
      packed |= static_cast<std::uint32_t>(alphatab[input[p + length]])
                << (30U - 2U * length);
      ++length;
    }

    const auto position = static_cast<std::uint32_t>(p);

    if (length == 0) {
      ++separators;
    } else if (length < KMER_LENGTH) {
      short_queries.push_back({position, packed, length});
    } else if (packed == ms2::POLY_A_KEY) {
      run_positions.push_back({position, 'A'});
    } else if (packed == ms2::POLY_T_KEY) {
      run_positions.push_back({position, 'T'});
    } else {
      keys.push_back(static_cast<std::uint64_t>(packed) << 32 | position);
    }
  }
}

void check(const std::string& name, const std::string& text) {
  const std::vector<unsigned char> input(text.begin(), text.end());

  std::vector<std::uint64_t> want_keys;
  std::vector<std::pair<std::uint32_t, char>> want_runs;
  std::vector<ms2::ShortQuery> want_short;
  std::size_t want_separators = 0;
  expected(input, want_keys, want_runs, want_short, want_separators);

  ms2::KeyedInput got;
  ms2::prepare_keys(input, got);

  if (got.keys != want_keys) {
    fail(name, "keys differ (" + std::to_string(got.keys.size()) + " vs " +
                   std::to_string(want_keys.size()) + ")");
  }

  // Runs: expand to positions; also each run must be maximal (not
  // continued by the next run of the same base).
  std::vector<std::pair<std::uint32_t, char>> got_runs;
  for (std::size_t r = 0; r < got.runs.size(); ++r) {
    const ms2::HomopolymerRun& run = got.runs[r];

    if (run.count == 0) fail(name, "empty run");

    for (std::uint32_t k = 0; k < run.count; ++k) {
      got_runs.push_back({run.first + k, run.base});
    }

    if (r + 1 < got.runs.size() && got.runs[r + 1].base == run.base &&
        got.runs[r + 1].first == run.first + run.count) {
      fail(name, "run not maximal at " + std::to_string(run.first));
    }
  }

  if (got_runs != want_runs) {
    fail(name, "runs differ (" + std::to_string(got_runs.size()) + " vs " +
                   std::to_string(want_runs.size()) + " positions)");
  }

  const auto same_short = [](const ms2::ShortQuery& a,
                             const ms2::ShortQuery& b) {
    return a.position == b.position && a.key == b.key && a.length == b.length;
  };

  if (got.short_queries.size() != want_short.size() ||
      !std::equal(got.short_queries.begin(), got.short_queries.end(),
                  want_short.begin(), same_short)) {
    fail(name, "short queries differ (" +
                   std::to_string(got.short_queries.size()) + " vs " +
                   std::to_string(want_short.size()) + ")");
  }

  if (got.separators != want_separators) {
    fail(name, "separators differ (" + std::to_string(got.separators) +
                   " vs " + std::to_string(want_separators) + ")");
  }

  if (got.keys.size() + got_runs.size() + got.short_queries.size() +
          got.separators !=
      input.size()) {
    fail(name, "positions not all accounted for");
  }

  // Orders: stable by bucket, and by key (ties in text order, which for
  // packed values is just their numeric order).
  std::vector<std::uint64_t> want_bucket = want_keys;
  std::stable_sort(want_bucket.begin(), want_bucket.end(),
                   [](const std::uint64_t a, const std::uint64_t b) {
                     return (a >> 48) < (b >> 48);
                   });

  std::vector<std::uint64_t> want_sorted = want_keys;
  std::sort(want_sorted.begin(), want_sorted.end());

  std::vector<std::uint64_t> bucket, sorted, scratch;
  ms2::bucket_order(got.keys, bucket);
  ms2::sorted_order(got.keys, sorted, scratch);

  if (bucket != want_bucket) fail(name, "bucket order differs");
  if (sorted != want_sorted) fail(name, "sorted order differs");
}

std::string random_dna(std::mt19937& rng, const std::size_t n) {
  std::string s(n, 'A');
  for (char& c : s) c = "ACGT"[rng() % 4];
  return s;
}

}  // namespace

int main() {
  std::mt19937 rng(5);
  int cases = 0;

  const auto run = [&](const std::string& name, const std::string& text) {
    ++cases;
    check(name, text);
  };

  run("empty", "");
  run("one char", "A");
  run("15 chars", "ACGTACGTACGTACG");
  run("16 chars", "ACGTACGTACGTACGT");
  run("exactly poly-A 16", std::string(16, 'A'));
  run("poly-A 100", std::string(100, 'A'));
  run("poly-T 100", std::string(100, 'T'));
  run("A then T", std::string(40, 'A') + std::string(40, 'T'));
  run("T then A", std::string(40, 'T') + std::string(40, 'A'));
  run("A run cut by N",
      random_dna(rng, 30) + std::string(25, 'A') + "N" + std::string(25, 'A') +
          random_dna(rng, 30));
  run("A 15 then N", std::string(15, 'A') + "N" + random_dna(rng, 20));
  run("A 16 at end", random_dna(rng, 40) + std::string(16, 'A'));
  run("A 17 at start", std::string(17, 'A') + random_dna(rng, 40));
  run("two A runs", random_dna(rng, 20) + std::string(30, 'A') + "C" +
                        std::string(30, 'A') + random_dna(rng, 20));
  run("only N", std::string(50, 'N'));
  run("lowercase", "acgtACGTACGTACGTACGTacgt" + random_dna(rng, 30));
  run("separators everywhere", "ACGTNACGTACGTACGTACGTACGNNNNACGTACGTACGTACGTAC");

  for (int t = 0; t < 30; ++t) {
    std::string s;
    while (s.size() < 3000) {
      const auto r = rng() % 10;
      if (r < 6) s += random_dna(rng, 1 + rng() % 200);
      else if (r < 8) s += std::string(1 + rng() % 60, "AT"[rng() % 2]);
      else s += std::string(1 + rng() % 5, "NnX"[rng() % 3]);
    }
    run("random " + std::to_string(t), s);
  }

  std::cout << "cases=" << cases << " failures=" << failures << '\n';
  if (failures != 0) return 1;
  std::cout << "ALL TESTS PASSED\n";
  return 0;
}
