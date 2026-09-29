// Tests ms2::lazy_matching_statistics against brute force, with the v2
// and sassy tables.
//
//   g++ -std=c++20 -O2 ms2/lazy_ms_test.cpp -o lazy_ms_test && ./lazy_ms_test

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../lrf_ms/probe_pipeline.hpp"  // table builders, policies
#include "lazy_ms.hpp"

namespace fs = std::filesystem;
using Symbol = unsigned char;

namespace {

int failures = 0;

void fail(const std::string& what) {
  if (failures++ < 10) std::cerr << "FAIL " << what << '\n';
}

// Only equal ACGT characters match (a separator matches nothing).
std::uint32_t brute_ms(const std::vector<Symbol>& input, const std::size_t i,
                       const std::vector<Symbol>& reference) {
  std::uint32_t best = 0;
  for (std::size_t p = 0; p < reference.size(); ++p) {
    std::uint32_t k = 0;
    while (i + k < input.size() && p + k < reference.size() &&
           is_acgt(input[i + k]) && input[i + k] == reference[p + k]) {
      ++k;
    }
    best = std::max(best, k);
  }
  return best;
}

std::vector<std::uint32_t> suffix_array(const std::vector<Symbol>& text) {
  std::vector<std::uint32_t> sa(text.size());
  std::iota(sa.begin(), sa.end(), 0);
  std::sort(sa.begin(), sa.end(), [&](std::uint32_t a, std::uint32_t b) {
    return std::lexicographical_compare(text.begin() + a, text.end(),
                                        text.begin() + b, text.end());
  });
  return sa;
}

template <typename Policy>
void check_one(const std::string& name, const typename Policy::Table& table,
               const std::vector<Symbol>& input,
               const std::vector<Symbol>& reference) {
  ms2::LazyStats stats;
  const MatchingStatistics ms =
      ms2::lazy_matching_statistics<Policy>(table, input, reference, &stats);

  if (stats.extended + stats.breaks + stats.separators != input.size()) {
    fail(name + ": positions not all accounted for");
  }

  for (std::size_t i = 0; i < input.size(); ++i) {
    const auto [position, length] = ms[i];
    const std::uint32_t want = brute_ms(input, i, reference);

    bool genuine = position + length <= reference.size();
    for (std::uint32_t k = 0; genuine && k < length; ++k) {
      genuine = input[i + k] == reference[position + k];
    }

    if (length != want || !genuine) {
      fail(name + ": at " + std::to_string(i) + " brute " +
           std::to_string(want) + " got " + std::to_string(length) + "@" +
           std::to_string(position) + (genuine ? "" : " (not a real match)"));
      return;
    }
  }
}

void check(const std::string& name, const std::string& reference_text,
           const std::vector<std::string>& inputs) {
  const std::vector<Symbol> reference(reference_text.begin(),
                                      reference_text.end());
  const std::vector<std::uint32_t> sa = suffix_array(reference);

  const fs::path dir = fs::temp_directory_path();
  const std::string v2_path = (dir / "lazy_ms_test.v2bin").string();
  const std::string sassy_path = (dir / "lazy_ms_test.sassy").string();

  build_pt16_sassy_table(reference, sa, sassy_path);
  const PT16SassyLookup sassy(sassy_path);

  // Plain v2 does not handle separators in the reference.
  const bool acgt = std::all_of(reference.begin(), reference.end(),
                                [](Symbol c) { return is_acgt(c); });

  std::unique_ptr<PT16RLZParser<Symbol, std::uint32_t>> v2;
  if (acgt) {
    build_pt16_table(reference, sa, v2_path);
    v2 = std::make_unique<PT16RLZParser<Symbol, std::uint32_t>>(reference, sa,
                                                                v2_path);
  }

  for (std::size_t k = 0; k < inputs.size(); ++k) {
    // As in the benchmarks: every non-ACGT input byte becomes one byte
    // that does not occur in the reference.
    std::vector<Symbol> input(inputs[k].begin(), inputs[k].end());
    for (Symbol& c : input) {
      if (!is_acgt(c)) c = 1;
    }

    const std::string label = name + " input " + std::to_string(k);
    check_one<msbench::SassyPolicy>(label + " sassy", sassy, input, reference);
    if (v2) check_one<msbench::V2Policy>(label + " v2", *v2, input, reference);
  }

  fs::remove(v2_path);
  fs::remove(sassy_path);
}

std::string random_dna(std::mt19937& rng, const std::size_t n) {
  std::string s(n, 'A');
  for (char& c : s) c = "ACGT"[rng() % 4];
  return s;
}

std::string mutate(std::mt19937& rng, std::string s, const int per_mille) {
  for (char& c : s) {
    if (static_cast<int>(rng() % 1000) < per_mille) c = "ACGT"[rng() % 4];
  }
  return s;
}

}  // namespace

int main() {
  std::mt19937 rng(11);
  int cases = 0;

  for (int t = 0; t < 12; ++t) {
    // A reference with repeats and homopolymer runs; with separators in
    // some cases (then only sassy runs).
    std::string reference = random_dna(rng, 400);
    const std::string unit = random_dna(rng, 60);
    for (int r = 0; r < 6; ++r) reference += mutate(rng, unit, 30);
    reference += std::string(80, "ACGT"[t % 4]);
    reference += random_dna(rng, 300);
    for (int r = 0; r < 8; ++r) reference += "CA";
    reference += random_dna(rng, 200);
    if (t % 3 == 2) {
      reference.insert(250, std::string(20, 'N'));
      reference.insert(700, "N");
    }

    std::vector<std::string> inputs;
    // A mutated copy of part of the reference (mostly extends).
    std::string copy = mutate(rng, reference.substr(100, 900), 10);
    for (char& c : copy) {
      if (c == 'N') c = "ACGT"[rng() % 4];
    }
    inputs.push_back(copy);
    // Unrelated sequence (mostly breaks).
    inputs.push_back(random_dna(rng, 400));
    // Runs and separators.
    inputs.push_back(std::string(120, 'A') + "N" + std::string(40, 'C') +
                     "NNNN" + reference.substr(500, 200) + std::string(30, 'T'));
    // Short inputs.
    inputs.push_back("");
    inputs.push_back("ACGT");
    inputs.push_back("NACGTACGTACGTACGTACGTN");

    check("case " + std::to_string(t), reference, inputs);
    ++cases;
  }

  std::cout << "cases=" << cases << " failures=" << failures << '\n';
  if (failures != 0) return 1;
  std::cout << "ALL TESTS PASSED\n";
  return 0;
}
