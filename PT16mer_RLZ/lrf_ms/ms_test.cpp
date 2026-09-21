// Standalone correctness test for LRFMS and the ms_utils checks.
//
// Builds its own suffix arrays with std::sort, independently of the
// pipeline's SA construction, and checks LRFMS against a brute-force
// matching-statistics computation, and exercises the ms_utils length
// comparison and digest helpers.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "lrf_ms.hpp"
#include "ms_utils.hpp"

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

std::vector<Symbol> random_sequence(std::mt19937_64& rng, std::size_t length,
                                    int alphabet, int first_symbol) {
  std::uniform_int_distribution<int> distribution(0, alphabet - 1);

  std::vector<Symbol> sequence(length);

  for (std::size_t i = 0; i < length; ++i) {
    sequence[i] =
        static_cast<Symbol>(first_symbol + distribution(rng));
  }

  return sequence;
}

std::vector<Symbol> repetitive_sequence(std::mt19937_64& rng,
                                        std::size_t length, int alphabet,
                                        std::size_t period) {
  const std::vector<Symbol> unit = random_sequence(rng, period, alphabet, 'a');

  std::vector<Symbol> sequence(length);

  for (std::size_t i = 0; i < length; ++i) {
    sequence[i] = unit[i % unit.size()];
  }

  // A few point mutations, so it is repetitive but not perfectly so.
  std::uniform_int_distribution<std::size_t> where(0, length - 1);
  std::uniform_int_distribution<int> what(0, alphabet - 1);

  for (int mutation = 0; mutation < 5 && length > 0; ++mutation) {
    sequence[where(rng)] = static_cast<Symbol>('a' + what(rng));
  }

  return sequence;
}

void check(const std::string& name, const std::vector<Symbol>& reference,
           const std::vector<Symbol>& input) {
  ++cases;

  const std::vector<SAType> sa = build_suffix_array(reference);

  LRFMS<Symbol, SAType> index(reference, sa);
  const MatchingStatistics ms = index.computeMatchingStatistics(input);

  const msbench::SymbolTable<Symbol> alphabet(reference);

  // 1. Structural invariants.
  const msbench::Validation invariants =
      msbench::validate_invariants(ms, input, reference, alphabet);

  if (!invariants.ok) {
    std::cerr << "FAIL [" << name << "] invariants: " << invariants << "\n";
    ++failures;
    return;
  }

  // 2. Every reported occurrence really matches.
  const msbench::Validation occurrences =
      msbench::verify_all_positions(ms, input, reference, sa, true);

  if (!occurrences.ok) {
    std::cerr << "FAIL [" << name << "] occurrences: " << occurrences << "\n";
    ++failures;
    return;
  }

  // 3. Lengths agree with brute force.
  const msbench::Validation brute =
      msbench::verify_against_brute_force(ms, input, reference);

  if (!brute.ok) {
    std::cerr << "FAIL [" << name << "] brute force: " << brute << "\n";
    std::cerr << "     reference_size=" << reference.size()
              << " input_size=" << input.size() << "\n";
    ++failures;
    return;
  }
}

void check_length_comparison() {
  ++cases;

  MatchingStatistics baseline_ms;
  baseline_ms.emplace_back(3, 7);
  baseline_ms.emplace_back(0, 6);
  baseline_ms.emplace_back(12, 5);
  baseline_ms.emplace_back(1, 4);

  const std::vector<std::uint32_t> baseline =
      msbench::extract_lengths(baseline_ms);

  if (baseline != std::vector<std::uint32_t>({7, 6, 5, 4})) {
    std::cerr << "FAIL [lengths] extract_lengths wrong\n";
    ++failures;
    return;
  }

  // Identical lengths, different positions: must compare equal, because
  // positions are not canonical.
  MatchingStatistics same_lengths;
  same_lengths.emplace_back(99, 7);
  same_lengths.emplace_back(98, 6);
  same_lengths.emplace_back(97, 5);
  same_lengths.emplace_back(96, 4);

  const auto equal = msbench::compare_lengths(baseline, same_lengths);

  if (!equal.checked || !equal.equal) {
    std::cerr << "FAIL [lengths] equal lengths reported as differing\n";
    ++failures;
    return;
  }

  // One divergence, at a known index.
  MatchingStatistics one_off = baseline_ms;
  std::get<1>(one_off[2]) = 4;

  const auto differs = msbench::compare_lengths(baseline, one_off);

  if (differs.equal || differs.first_divergence != 2 ||
      differs.baseline_value != 5 || differs.variant_value != 4 ||
      differs.divergent_positions != 1) {
    std::cerr << "FAIL [lengths] one-off divergence misreported: "
              << differs.describe() << "\n";
    ++failures;
    return;
  }

  // Wholesale divergence: the count must reflect it.
  MatchingStatistics all_different;
  all_different.emplace_back(0, 1);
  all_different.emplace_back(0, 1);
  all_different.emplace_back(0, 1);
  all_different.emplace_back(0, 1);

  const auto wholesale = msbench::compare_lengths(baseline, all_different);

  if (wholesale.equal || wholesale.first_divergence != 0 ||
      wholesale.divergent_positions != 4) {
    std::cerr << "FAIL [lengths] wholesale divergence misreported: "
              << wholesale.describe() << "\n";
    ++failures;
    return;
  }

  // Different entry counts must be caught, not indexed past.
  MatchingStatistics short_ms;
  short_ms.emplace_back(3, 7);

  const auto mismatched = msbench::compare_lengths(baseline, short_ms);

  if (mismatched.equal) {
    std::cerr << "FAIL [lengths] differing entry counts reported equal\n";
    ++failures;
  }
}

void check_length_hash_discriminates() {
  ++cases;

  MatchingStatistics a;
  a.emplace_back(1, 5);
  a.emplace_back(2, 4);

  MatchingStatistics b;
  b.emplace_back(1, 5);
  b.emplace_back(2, 3);

  MatchingStatistics c;  // same lengths, different positions
  c.emplace_back(9, 5);
  c.emplace_back(8, 4);

  const auto da = msbench::digest_ms(a);
  const auto db = msbench::digest_ms(b);
  const auto dc = msbench::digest_ms(c);

  if (da.len_hash == db.len_hash) {
    std::cerr << "FAIL [hash] different lengths hashed equal\n";
    ++failures;
  }

  if (da.len_hash != dc.len_hash) {
    std::cerr << "FAIL [hash] same lengths hashed differently\n";
    ++failures;
  }

  if (da.pos_hash == dc.pos_hash) {
    std::cerr << "FAIL [hash] different positions hashed equal\n";
    ++failures;
  }
}

}  // namespace

int main() {
  std::mt19937_64 rng(20260921);

  // ---------- Hand-built edge cases ----------

  check("single symbol reference", {'a'}, {'a', 'a', 'b'});
  check("input all absent", {'a', 'b', 'c'}, {'x', 'y', 'z'});
  check("input equals reference", {'a', 'b', 'r', 'a', 'c', 'a'},
        {'a', 'b', 'r', 'a', 'c', 'a'});
  check("all same symbol", std::vector<Symbol>(40, 'a'),
        std::vector<Symbol>(60, 'a'));
  check("banana", {'b', 'a', 'n', 'a', 'n', 'a'},
        {'a', 'n', 'a', 'n', 'a', 'b', 'a', 'n', 'a'});

  // Byte 0 is an ordinary symbol now: this is the case the original
  // sentinel substitution got wrong.
  check("alphabet includes 0", {0, 1, 0, 2, 0, 1, 1, 0},
        {0, 1, 1, 0, 2, 0, 0, 1, 3});
  check("zero-heavy", {0, 0, 1, 0, 0, 0, 2, 0}, {0, 0, 0, 1, 0, 2, 0, 0});

  // Last reference position reachable as a genuine suffix.
  check("suffix of reference", {'m', 'i', 's', 's', 'i', 's', 's', 'i'},
        {'s', 'i'});

  check("single symbol input", {'a', 'b'}, {'b'});
  check("empty input", {'a', 'b'}, {});

  // ---------- Randomised ----------

  for (int round = 0; round < 300; ++round) {
    const std::size_t reference_size = 1 + (rng() % 200);
    const std::size_t input_size = 1 + (rng() % 200);

    // Binary alphabet: maximum tie density in the binary searches.
    check("random binary " + std::to_string(round),
          random_sequence(rng, reference_size, 2, 'a'),
          random_sequence(rng, input_size, 2, 'a'));

    // DNA-like.
    check("random dna " + std::to_string(round),
          random_sequence(rng, reference_size, 4, 'A'),
          random_sequence(rng, input_size, 4, 'A'));

    // Input drawn from a larger alphabet, so some symbols are absent.
    check("absent symbols " + std::to_string(round),
          random_sequence(rng, reference_size, 3, 'a'),
          random_sequence(rng, input_size, 6, 'a'));

    // Alphabet anchored at byte 0.
    check("zero-based alphabet " + std::to_string(round),
          random_sequence(rng, reference_size, 3, 0),
          random_sequence(rng, input_size, 4, 0));
  }

  for (int round = 0; round < 120; ++round) {
    const std::size_t reference_size = 40 + (rng() % 300);
    const std::size_t input_size = 40 + (rng() % 300);
    const std::size_t period = 2 + (rng() % 12);

    const std::vector<Symbol> reference =
        repetitive_sequence(rng, reference_size, 3, period);

    // Long matches, which is what drives the LRF skip.
    check("repetitive " + std::to_string(round), reference,
          repetitive_sequence(rng, input_size, 3, period));

    // A real substring of the reference, extended with noise.
    std::vector<Symbol> derived(reference.begin(),
                                reference.begin() + reference.size() / 2);
    const std::vector<Symbol> noise = random_sequence(rng, 10, 4, 'a');
    derived.insert(derived.end(), noise.begin(), noise.end());
    derived.insert(derived.end(), reference.begin(), reference.end());

    check("derived " + std::to_string(round), reference, derived);
  }

  // ---------- ms_utils ----------

  check_length_comparison();
  check_length_hash_discriminates();

  std::cout << "\ncases=" << cases << " failures=" << failures << "\n";

  if (failures == 0) {
    std::cout << "ALL TESTS PASSED\n";
    return 0;
  }

  std::cout << "TESTS FAILED\n";
  return 1;
}
