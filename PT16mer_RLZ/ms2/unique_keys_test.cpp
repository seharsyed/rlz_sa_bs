// Tests ms2::keep_unique, the v2 finger lookup, and that all four lookup
// variants agree on the unique keys.
//
//   g++ -std=c++20 -O2 ms2/unique_keys_test.cpp -o unique_keys_test
//   ./unique_keys_test

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "unique_chain.hpp"
#include "unique_keys.hpp"

namespace fs = std::filesystem;
using Symbol = unsigned char;

namespace {

int failures = 0;

void fail(const std::string& what) {
  if (failures++ < 10) std::cerr << "FAIL " << what << '\n';
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

bool same(const KmerLookupResult& a, const KmerLookupResult& b) {
  if (a.found != b.found || a.match_length != b.match_length ||
      a.match_position != b.match_position) {
    return false;
  }
  return !a.found || (a.count == b.count &&
                      msbench::same_occurrences(msbench::occurrences(a),
                                                msbench::occurrences(b)));
}

std::string random_dna(std::mt19937& rng, const std::size_t n) {
  std::string s(n, 'A');
  for (char& c : s) c = "ACGT"[rng() % 4];
  return s;
}

void check(const std::string& name, const std::string& reference_text,
           const std::string& input_text) {
  const std::vector<Symbol> reference(reference_text.begin(),
                                      reference_text.end());
  const std::vector<Symbol> input(input_text.begin(), input_text.end());
  const std::vector<std::uint32_t> sa = suffix_array(reference);

  // ---------- keep_unique against a count of every key ----------

  ms2::KeyedInput keyed;
  ms2::prepare_keys(input, keyed);

  std::vector<std::uint64_t> sorted, scratch, unique;
  ms2::sorted_order(keyed.keys, sorted, scratch);
  const ms2::UniqueStats stats = ms2::keep_unique(sorted, unique);

  std::map<std::uint32_t, std::size_t> count;
  for (const std::uint64_t packed : keyed.keys) ++count[ms2::key_of(packed)];

  std::vector<std::uint64_t> want;
  for (const std::uint64_t packed : sorted) {
    if (count[ms2::key_of(packed)] == 1) want.push_back(packed);
  }

  if (unique != want) {
    fail(name + ": unique keys differ (" + std::to_string(unique.size()) +
         " vs " + std::to_string(want.size()) + ")");
  }
  if (stats.unique + stats.repeated_positions != keyed.keys.size()) {
    fail(name + ": unique + repeated positions != keys");
  }

  // ---------- Tables ----------

  const fs::path dir = fs::temp_directory_path();
  const std::string v2_path = (dir / "unique_keys_test.v2bin").string();
  const std::string sassy_path = (dir / "unique_keys_test.sassy").string();

  build_pt16_table(reference, sa, v2_path);
  build_pt16_sassy_table(reference, sa, sassy_path);
  const PT16RLZParser<Symbol, std::uint32_t> v2(reference, sa, v2_path);
  const PT16SassyLookup sassy(sassy_path);

  // ---------- v2 finger = v2 plain, in any order ----------

  std::vector<std::uint32_t> keys;
  for (const std::uint64_t packed : keyed.keys) keys.push_back(ms2::key_of(packed));

  std::vector<std::vector<std::uint32_t>> orders;
  orders.push_back(keys);  // text order: restarts a lot
  std::vector<std::uint32_t> sorted_keys = keys;
  std::sort(sorted_keys.begin(), sorted_keys.end());
  orders.push_back(sorted_keys);
  std::vector<std::uint32_t> doubled;
  for (const std::uint32_t key : sorted_keys) {
    doubled.push_back(key);
    doubled.push_back(key);
  }
  orders.push_back(doubled);

  for (std::size_t o = 0; o < orders.size(); ++o) {
    PT16RLZParser<Symbol, std::uint32_t>::Finger finger;
    for (const std::uint32_t key : orders[o]) {
      if (!same(v2.lookupKmerByKey(key, finger), v2.lookupKmerByKey(key))) {
        fail(name + ": v2 finger differs from plain in order " +
             std::to_string(o) + " at key " + std::to_string(key));
        break;
      }
    }
  }

  // ---------- All four variants agree on the unique keys ----------

  std::vector<KmerLookupResult> first, other;
  ms2::probe_compact<msbench::V2Policy>(v2, unique, first);

  const auto against_first = [&](const char* variant,
                                 const std::vector<KmerLookupResult>& got) {
    std::size_t at = 0;
    std::string reason;
    const std::size_t mismatches =
        ms2::compare_compact(first, got, unique, input, reference, at, reason);
    if (mismatches != 0) {
      fail(name + ": " + variant + " differs at " + std::to_string(at) + " (" +
           reason + ")");
    }
  };

  against_first("v2", first);  // its own positions must be genuine too
  ms2::probe_compact<msbench::V2FingerPolicy>(v2, unique, other);
  against_first("v2-finger", other);
  ms2::probe_compact<msbench::SassyPolicy>(sassy, unique, other);
  against_first("sassy", other);
  ms2::probe_compact<msbench::SassyFingerPolicy>(sassy, unique, other);
  against_first("sassy-finger", other);

  // A hit's count is the number of occurrences of the 16-mer.
  for (std::size_t k = 0; k < unique.size(); ++k) {
    if (!first[k].found) continue;
    const std::uint32_t position = ms2::position_of(unique[k]);
    std::uint32_t occurrences = 0;
    for (std::size_t p = 0; p + KMER_LENGTH <= reference.size(); ++p) {
      occurrences += std::equal(input.begin() + position,
                                input.begin() + position + KMER_LENGTH,
                                reference.begin() + p);
    }
    if (occurrences != first[k].count) {
      fail(name + ": count " + std::to_string(first[k].count) + " but " +
           std::to_string(occurrences) + " occurrences at " +
           std::to_string(position));
      break;
    }
  }

  // ---------- The chain: text order, then the sweep ----------

  std::vector<std::uint64_t> by_position, packed, scratch2;
  ms2::text_order(unique, by_position, packed, scratch2);

  for (std::size_t idx = 1; idx < by_position.size(); ++idx) {
    if (ms2::key_of(by_position[idx - 1]) >= ms2::key_of(by_position[idx])) {
      fail(name + ": text order not increasing");
      break;
    }
  }

  std::vector<ms2::ChainEntry> chained;
  const ms2::UniqueChainStats chain_stats =
      ms2::chain_unique(first, by_position, chained);

  if (chain_stats.exact() + chain_stats.approximate() != unique.size()) {
    fail(name + ": chain entries not all accounted for");
  }

  for (std::size_t idx = 0; idx < chained.size(); ++idx) {
    const std::uint32_t i = ms2::key_of(by_position[idx]);
    const ms2::ChainEntry& entry = chained[idx];
    const std::uint32_t truth = brute_ms(input, i, reference);

    bool genuine = entry.reference_position + entry.length <= reference.size();
    for (std::uint32_t k = 0; genuine && k < entry.length; ++k) {
      genuine = input[i + k] == reference[entry.reference_position + k];
    }

    const bool ok = ms2::is_exact(entry.status) ? entry.length == truth
                                                : entry.length <= truth;
    if (!ok || !genuine) {
      fail(name + ": chain at " + std::to_string(i) + " got " +
           std::to_string(entry.length) + (ms2::is_exact(entry.status)
                                               ? " (exact)"
                                               : " (approximate)") +
           ", brute " + std::to_string(truth) +
           (genuine ? "" : " (not a real match)"));
      break;
    }
  }

  fs::remove(v2_path);
  fs::remove(sassy_path);
}

}  // namespace

int main() {
  std::mt19937 rng(21);
  int cases = 0;

  for (int t = 0; t < 10; ++t) {
    // A reference with repeats, so there are ranges and misses.
    std::string reference = random_dna(rng, 1500);
    const std::string unit = random_dna(rng, 80);
    for (int r = 0; r < 5; ++r) reference += unit;
    reference += std::string(60, "ACGT"[t % 4]);
    reference += random_dna(rng, 1000);

    // Input: parts of the reference (hits, some repeated so not unique),
    // random sequence (misses), a repeat copied twice, runs, separators.
    std::string input = reference.substr(200, 700) + random_dna(rng, 300) +
                        unit + unit + "N" + std::string(40, 'A') +
                        reference.substr(1400, 500) + reference.substr(300, 100);

    check("case " + std::to_string(t), reference, input);
    ++cases;
  }

  // A gap that hides a difference: a 32-character reference segment with a
  // SNP in its middle appears twice in the input (once in its reference
  // context, once between random flanks), so its inner 16-mers are
  // repeated and dropped. The unique 16-mers on either side of that
  // 18-position gap both match the reference and line up, but the SNP
  // lies between them: the chain must not link across.
  for (int t = 0; t < 5; ++t) {
    const std::string reference = random_dna(rng, 3000);
    const std::size_t a = 1000;
    std::string context = reference.substr(a - 60, 60 + 32 + 60);
    std::string segment = context.substr(60, 32);
    const char old = segment[16];
    segment[16] = old == 'A' ? 'C' : 'A';
    context.replace(60, 32, segment);

    const std::string input = random_dna(rng, 100) + context +
                              random_dna(rng, 50) + segment +
                              random_dna(rng, 50);
    check("hidden SNP " + std::to_string(t), reference, input);
    ++cases;
  }

  check("tiny", "ACGTACGTACGTACGTACGTAAAA", "ACGTACGTACGTACGTACGTAAAACCCC");
  check("empty input", random_dna(rng, 200), "");
  cases += 2;

  std::cout << "cases=" << cases << " failures=" << failures << '\n';
  if (failures != 0) return 1;
  std::cout << "ALL TESTS PASSED\n";
  return 0;
}
