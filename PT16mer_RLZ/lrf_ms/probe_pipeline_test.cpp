// Standalone test of the shared PT16 pipeline (probe_pipeline.hpp).
//
// ms_main chains only one variant's lookup results and checks every other
// variant's against them with compare_lookups. That is only sound if
// agreeing results really do chain to the same matching statistics, so
// this test checks it end to end, for every variant and both probe
// orders:
//
//   - compare_lookups finds every variant's results equal to the first's;
//   - chaining each variant's own results gives the same lengths, and
//     those lengths equal brute force at every position;
//   - every reported position is a real match.
//
// References are small, and some use only 1 or 2 letters, so the cases
// the variants differ on are all exercised: empty buckets, short suffixes
// at the end of the reference, range hits (repeated 16-mers), and inputs
// shorter than 16.
//
//   g++ -std=c++20 -O2 lrf_ms/probe_pipeline_test.cpp -o probe_pipeline_test
//   ./probe_pipeline_test

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "ms_tools.hpp"
#include "probe_pipeline.hpp"

namespace fs = std::filesystem;

using msbench::SAType;
using msbench::Symbol;

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

// Random text over the first `letters` of ACGT.
std::vector<Symbol> random_dna(std::mt19937_64& rng, std::size_t n,
                               int letters) {
  static const char alphabet[] = {'A', 'C', 'G', 'T'};
  std::uniform_int_distribution<int> base(0, letters - 1);
  std::vector<Symbol> s(n);
  for (auto& c : s) c = static_cast<Symbol>(alphabet[base(rng)]);
  return s;
}

std::uint32_t brute_ms_length(const std::vector<Symbol>& reference,
                              const std::vector<Symbol>& input,
                              std::size_t i) {
  std::uint32_t best = 0;

  for (std::size_t p = 0; p < reference.size(); ++p) {
    std::uint32_t length = 0;

    while (i + length < input.size() && p + length < reference.size() &&
           reference[p + length] == input[i + length]) {
      ++length;
    }

    best = std::max(best, length);
  }

  return best;
}

void check(const std::string& name, const std::vector<Symbol>& reference,
           const std::vector<std::vector<Symbol>>& inputs) {
  ++cases;

  const std::vector<SAType> sa = build_suffix_array(reference);
  const std::string table_path =
      (fs::temp_directory_path() / "probe_pipeline_test.v2bin").string();

  msbench::ProberSet set = msbench::build_probers(reference, sa, table_path);

  int mismatches = 0;
  const auto fail = [&](const std::string& what) {
    if (mismatches++ < 5) {
      std::cerr << "FAIL [" << name << "] " << what << "\n";
    }
  };

  for (const std::vector<Symbol>& input : inputs) {
    const std::size_t n = input.size();

    std::vector<std::uint32_t> keys;
    msbench::roll_keys(input, keys);
    const std::uint32_t tail_key = msbench::first_tail_key(input, keys);

    std::vector<std::uint32_t> want(n);
    for (std::size_t i = 0; i < n; ++i) {
      want[i] = brute_ms_length(reference, input, i);
    }

    std::vector<KmerLookupResult> first(n);
    bool have_first = false;

    for (const msbench::ProbeOrder probe_order : msbench::probe_orders) {
      std::vector<std::uint32_t> order;
      std::vector<std::uint32_t> scratch;

      if (probe_order == msbench::ProbeOrder::bucket) {
        msbench::bucket_order(keys, order, scratch);
      } else {
        msbench::sorted_order(keys, order, scratch);

        // The MSD build of the same order must give it exactly.
        std::vector<std::uint32_t> msd_order;
        std::vector<std::uint64_t> packed;
        msbench::sorted_order_msd(keys, msd_order, packed);

        if (msd_order != order) {
          fail("sorted_order_msd differs from sorted_order (n=" +
               std::to_string(n) + ")");
        }
      }

      for (auto& prober : set.probers) {
        if (!prober->supports(probe_order)) {
          continue;
        }

        const std::string row = prober->name() + "-" +
                                msbench::order_name(probe_order) +
                                " n=" + std::to_string(n);

        std::vector<KmerLookupResult> results(n);
        prober->probe(keys, order, tail_key, results);

        if (!have_first) {
          first = results;
          have_first = true;
        }

        const msbench::LookupComparison comparison =
            msbench::compare_lookups(first, results, input, reference);

        if (!comparison.equal) {
          fail(row + ": lookups: " + comparison.describe());
        }

        // Chaining this variant's own results must give the true lengths.
        const MatchingStatistics ms = backwardChainExtend(results);

        for (std::size_t i = 0; i < n; ++i) {
          const auto [position, length] = ms[i];

          bool genuine = position + length <= reference.size();
          for (std::uint32_t k = 0; genuine && k < length; ++k) {
            genuine = reference[position + k] == input[i + k];
          }

          if (length != want[i] || !genuine) {
            fail(row + ": ms at " + std::to_string(i) + ": brute=" +
                 std::to_string(want[i]) + " got " + std::to_string(length) +
                 "@" + std::to_string(position) +
                 (genuine ? "" : " (not a real match)"));
            break;
          }
        }
      }
    }
  }

  // The finger lookups on their own: the same result as the plain lookup
  // for keys in any order -- text order (unsorted, so it keeps
  // restarting), every key twice in a row, and sorted -- since a smaller
  // key restarts.
  {
    const PT16SassyLookup sassy(table_path + ".sassy");
    const PT16FastMissParser<Symbol, SAType> fastmiss(reference, sa,
                                                      table_path);

    for (const std::vector<Symbol>& input : inputs) {
      std::vector<std::uint32_t> keys;
      msbench::roll_keys(input, keys);

      std::vector<std::uint32_t> sorted_keys = keys;
      std::sort(sorted_keys.begin(), sorted_keys.end());

      std::vector<std::uint32_t> doubled;
      for (const std::uint32_t key : keys) {
        doubled.push_back(key);
        doubled.push_back(key);
      }

      const auto same = [](const KmerLookupResult& a,
                           const KmerLookupResult& b) {
        return a.found == b.found && a.match_length == b.match_length &&
               a.match_position == b.match_position && a.count == b.count &&
               std::equal(a.positions.begin(), a.positions.end(),
                          b.positions.begin(), b.positions.end());
      };

      for (const auto* sequence : {&keys, &doubled, &sorted_keys}) {
        PT16SassyLookup::Finger sassy_finger;
        PT16FastMissParser<Symbol, SAType>::Finger fastmiss_finger;

        for (std::size_t j = 0; j < sequence->size(); ++j) {
          const std::uint32_t key = (*sequence)[j];

          if (!same(sassy.lookup(key, sassy_finger), sassy.lookup(key))) {
            fail("sassy finger lookup differs from plain lookup at key " +
                 std::to_string(j) + " (n=" + std::to_string(input.size()) +
                 ")");
            break;
          }

          if (!same(fastmiss.lookupKmerByKey(key, fastmiss_finger),
                    fastmiss.lookupKmerByKey(key))) {
            fail("fastmiss finger lookup differs from plain lookup at key " +
                 std::to_string(j) + " (n=" + std::to_string(input.size()) +
                 ")");
            break;
          }
        }
      }
    }
  }

  fs::remove(table_path);
  fs::remove(table_path + ".sassy");

  if (mismatches > 0) {
    ++failures;
    std::cerr << "  [" << name << "] " << mismatches << " mismatches\n";
  }
}

// Inputs of every interesting length: shorter than 16, exactly 16, and
// longer; half copied from the reference (long matches, chains), some
// ending on the reference's own last characters (short suffixes).
std::vector<std::vector<Symbol>> make_inputs(std::mt19937_64& rng,
                                             const std::vector<Symbol>& reference,
                                             int letters) {
  std::vector<std::vector<Symbol>> inputs;

  for (const std::size_t n : {1, 2, 7, 8, 9, 15, 16, 17, 31, 64, 300}) {
    std::vector<Symbol> input = random_dna(rng, n, letters);

    if (rng() % 2 == 0 && reference.size() >= n) {
      const std::size_t p = rng() % (reference.size() - n + 1);
      std::copy(reference.begin() + p, reference.begin() + p + n,
                input.begin());
    }

    if (rng() % 3 == 0) {
      const std::size_t m = std::min(n, reference.size());
      std::copy(reference.end() - m, reference.end(), input.end() - m);
    }

    inputs.push_back(std::move(input));
  }

  return inputs;
}

}  // namespace

int main() {
  std::mt19937_64 rng(20260924);

  // The table builders need at least 16 characters of reference.
  for (const std::size_t ref_size : {17, 40, 300, 3000}) {
    for (const int letters : {1, 2, 4}) {
      const auto reference = random_dna(rng, ref_size, letters);
      check("ref " + std::to_string(ref_size) + " over " +
                std::to_string(letters) + " letters",
            reference, make_inputs(rng, reference, letters));
    }
  }

  // A larger reference with repeats: most buckets non-empty, and plenty of
  // range hits.
  {
    std::vector<Symbol> reference;
    const auto unit = random_dna(rng, 37, 4);
    for (int r = 0; r < 40; ++r) {
      reference.insert(reference.end(), unit.begin(), unit.end());
    }
    const auto noise = random_dna(rng, 20000, 4);
    reference.insert(reference.end(), noise.begin(), noise.end());

    check("repeats + 20000 random", reference, make_inputs(rng, reference, 4));
  }

  std::cout << "cases=" << cases << " failures=" << failures << "\n";

  if (failures == 0) {
    std::cout << "ALL TESTS PASSED\n";
    return 0;
  }

  return 1;
}
