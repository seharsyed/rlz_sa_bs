// Tests ms2::sampled_matching_statistics against brute force, with the v2
// and sassy tables and several sample steps.
//
//   g++ -std=c++20 -O2 ms2/sampled_ms_test.cpp -o sampled_ms_test
//   ./sampled_ms_test

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../lrf_ms/probe_pipeline.hpp"  // table builders, policies
#include "sampled_ms.hpp"

namespace fs = std::filesystem;
using Symbol = unsigned char;

namespace {

int failures = 0;

void fail(const std::string& what) {
  if (failures++ < 10) std::cerr << "FAIL " << what << '\n';
}

// Only equal ACGT characters match.
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

template <typename Finger, typename Plain>
void check_one(const std::string& name, const typename Plain::Table& table,
               const std::vector<Symbol>& input,
               const std::vector<Symbol>& reference, const std::size_t step) {
  ms2::SampledMsBuffers buffers;
  ms2::SampledMsStats stats;
  const MatchingStatistics ms = ms2::sampled_matching_statistics<Finger, Plain>(
      table, input, reference, step, buffers, stats, nullptr);

  const std::size_t accounted = stats.chained + stats.separators +
                                stats.by_short + stats.by_singleton +
                                stats.by_range + stats.by_stretch;
  if (accounted != input.size()) {
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
  const std::string v2_path = (dir / "sampled_ms_test.v2bin").string();
  const std::string sassy_path = (dir / "sampled_ms_test.sassy").string();

  build_pt16_sassy_table(reference, sa, sassy_path);
  const PT16SassyLookup sassy(sassy_path);

  const bool acgt = std::all_of(reference.begin(), reference.end(),
                                [](Symbol c) { return is_acgt(c); });
  std::unique_ptr<PT16RLZParser<Symbol, std::uint32_t>> v2;
  if (acgt) {
    build_pt16_table(reference, sa, v2_path);
    v2 = std::make_unique<PT16RLZParser<Symbol, std::uint32_t>>(reference, sa,
                                                                v2_path);
  }

  for (std::size_t k = 0; k < inputs.size(); ++k) {
    std::vector<Symbol> input(inputs[k].begin(), inputs[k].end());
    for (Symbol& c : input) {
      if (!is_acgt(c)) c = 1;  // as in the benchmarks
    }

    for (const std::size_t step : {16, 8, 5}) {
      const std::string label = name + " input " + std::to_string(k) +
                                " step " + std::to_string(step);
      check_one<msbench::SassyFingerPolicy, msbench::SassyPolicy>(
          label + " sassy", sassy, input, reference, step);
      if (v2) {
        check_one<msbench::V2FingerPolicy, msbench::V2Policy>(
            label + " v2", *v2, input, reference, step);
      }
    }
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
  std::string out;
  for (char c : s) {
    const int r = static_cast<int>(rng() % 1000);
    if (r < per_mille / 2) out.push_back("ACGT"[rng() % 4]);       // substitution
    else if (r < per_mille * 3 / 4) {}                              // deletion
    else if (r < per_mille) { out.push_back(c); out.push_back("ACGT"[rng() % 4]); }  // insertion
    else out.push_back(c);
  }
  return out;
}

}  // namespace

int main() {
  std::mt19937 rng(13);
  int cases = 0;

  for (int t = 0; t < 12; ++t) {
    std::string reference = random_dna(rng, 500);
    const std::string unit = random_dna(rng, 70);
    for (int r = 0; r < 6; ++r) reference += mutate(rng, unit, 30);
    reference += std::string(80, "ACGT"[t % 4]);
    reference += random_dna(rng, 400);
    for (int r = 0; r < 10; ++r) reference += "CA";
    reference += random_dna(rng, 300);
    if (t % 3 == 2) {
      reference.insert(300, std::string(20, 'N'));
      reference.insert(800, "N");
    }

    std::vector<std::string> inputs;
    // A mutated copy (mostly stretches, some holes).
    std::string copy = mutate(rng, reference.substr(50, 1200), 8);
    for (char& c : copy) {
      if (c == 'N') c = "ACGT"[rng() % 4];
    }
    inputs.push_back(copy);
    // Unrelated sequence (mostly holes).
    inputs.push_back(random_dna(rng, 500));
    // Repeats, runs and separators.
    inputs.push_back(unit + unit + std::string(100, 'A') + "N" +
                     std::string(40, 'C') + "NNNN" + reference.substr(600, 300) +
                     std::string(30, 'T'));
    // Short inputs.
    inputs.push_back("");
    inputs.push_back("ACGT");
    inputs.push_back(reference.substr(100, 17));
    inputs.push_back("NACGTACGTACGTACGTACGTN");

    check("case " + std::to_string(t), reference, inputs);
    ++cases;
  }

  std::cout << "cases=" << cases << " failures=" << failures << '\n';
  if (failures != 0) return 1;
  std::cout << "ALL TESTS PASSED\n";
  return 0;
}
