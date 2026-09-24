#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../pt16_build_v2.hpp"
#include "../pt16_rlz_v2.hpp"
#include "../pt16_utils.hpp"
#include "../variants/pt16_build_sassy.hpp"
#include "../variants/pt16_rlz_v2_fastmiss.hpp"
#include "../variants/pt16_sassy.hpp"
#include "ms_utils.hpp"

/**
 * The shared PT16 matching-statistics pipeline.
 *
 * Every PT16 variant computes matching statistics the same way:
 *
 *   keys    roll every 16-mer key of the input            (shared)
 *   order   the order the keys are probed in:
 *             bucket  counting sort by table bucket       (shared)
 *             sorted  LSD radix sort by the full key      (shared)
 *   probe   one table lookup per key, in that order,
 *           plus the tail positions                       (per variant)
 *   chain   multi-step chain extension (ms_tools.hpp)     (shared)
 *
 * Only the probe depends on the table variant. The keys and the two orders
 * depend only on the input, and the chain only on the lookup results,
 * which every correct variant produces identically (each result is stored
 * at its own text position, whatever order it was probed in). So ms_main
 * runs keys, each order and the chain once per file, and each variant only
 * probes -- instead of every variant redoing all of it.
 *
 * A variant is a Prober (below); build_probers is the registry.
 */

namespace msbench {

using Symbol = unsigned char;
using SAType = std::uint32_t;

// ---------- Shared stages ----------

inline std::size_t kmer_positions_of(const std::size_t n) {
  return n >= KMER_LENGTH ? n - KMER_LENGTH + 1 : 0;
}

// keys[i] = the 16-mer key at input position i, for every position with 16
// characters left. Only position 0 pays a full encode; every later key is
// `(previous key << 2) | new trailing character`.
inline void roll_keys(const std::vector<Symbol>& input,
                      std::vector<std::uint32_t>& keys) {
  const std::size_t count = kmer_positions_of(input.size());
  keys.resize(count);

  if (count == 0) {
    return;
  }

  std::uint32_t key = encode_16mer(input, 0);
  keys[0] = key;

  for (std::size_t i = 1; i < count; ++i) {
    key = (key << 2U) |
          alphatab[static_cast<unsigned char>(input[i + KMER_LENGTH - 1])];
    keys[i] = key;
  }

  phase_barrier(keys.data());
}

// The padded key of the first tail position (see encode_tail): the last
// 16-mer key rolled one step past the end, or the whole input packed if it
// is shorter than 16. Each later tail position is `key << 2` of the
// previous one.
inline std::uint32_t first_tail_key(const std::vector<Symbol>& input,
                                    const std::vector<std::uint32_t>& keys) {
  if (!keys.empty()) {
    return keys.back() << 2U;
  }

  return input.empty() ? 0 : encode_tail(input, 0);
}

// One counting-sort pass: `in` reordered stably by digit_of(keys[i]).
template <typename DigitOf>
inline void counting_sort_pass(const std::vector<std::uint32_t>& keys,
                               const std::vector<std::uint32_t>& in,
                               std::vector<std::uint32_t>& out,
                               DigitOf digit_of) {
  auto count = std::make_unique<std::array<std::uint32_t, NUMBER_OF_BUCKETS>>();
  count->fill(0);

  for (const std::uint32_t i : in) {
    ++(*count)[digit_of(keys[i])];
  }

  std::uint32_t offset = 0;

  for (std::uint32_t& c : *count) {
    const std::uint32_t here = c;
    c = offset;
    offset += here;
  }

  out.resize(in.size());

  for (const std::uint32_t i : in) {
    out[(*count)[digit_of(keys[i])]++] = i;
  }
}

// Positions grouped by table bucket (the key's top LOW_BITS), in text order
// within a bucket: one counting-sort pass.
inline void bucket_order(const std::vector<std::uint32_t>& keys,
                         std::vector<std::uint32_t>& order,
                         std::vector<std::uint32_t>& identity) {
  identity.resize(keys.size());

  for (std::size_t i = 0; i < keys.size(); ++i) {
    identity[i] = static_cast<std::uint32_t>(i);
  }

  counting_sort_pass(keys, identity, order,
                     [](std::uint32_t key) { return key >> LOW_BITS; });
  phase_barrier(order.data());
}

// Positions fully sorted by their 32-bit key: two LSD counting-sort passes,
// low 16 bits then (stably) high 16 bits. Matches the table's own layout
// within a bucket too, and puts identical keys next to each other.
inline void sorted_order(const std::vector<std::uint32_t>& keys,
                         std::vector<std::uint32_t>& order,
                         std::vector<std::uint32_t>& scratch) {
  order.resize(keys.size());

  for (std::size_t i = 0; i < keys.size(); ++i) {
    order[i] = static_cast<std::uint32_t>(i);
  }

  counting_sort_pass(keys, order, scratch,
                     [](std::uint32_t key) { return key & LOW_MASK; });
  counting_sort_pass(keys, scratch, order,
                     [](std::uint32_t key) { return key >> LOW_BITS; });
  phase_barrier(order.data());
}

// ---------- Probers ----------

// The probe orders, in report order.
enum class ProbeOrder { bucket, sorted };

inline constexpr std::array<ProbeOrder, 2> probe_orders = {ProbeOrder::bucket,
                                                           ProbeOrder::sorted};

inline const char* order_name(const ProbeOrder order) {
  return order == ProbeOrder::bucket ? "bucket" : "sorted";
}


/**
 * One table variant under test: loads its table, then answers every lookup
 * of an input in a given order.
 */
class Prober {
 public:
  virtual ~Prober() = default;

  virtual const std::string& name() const = 0;

  // Wall time of loading this variant's table (and building whatever it
  // derives from it), in ms. Writing the table file is timed separately,
  // once per table format (see ProberSet::table_builds).
  virtual double build_ms() const = 0;

  // results[i] for every i in `order` (the 16-mer positions, from keys[i]),
  // then every tail position kmer_positions .. n-1 (found = false),
  // starting from `tail_key` (first_tail_key). `results` must already
  // have n entries.
  virtual void probe(const std::vector<std::uint32_t>& keys,
                     const std::vector<std::uint32_t>& order,
                     std::uint32_t tail_key,
                     std::vector<KmerLookupResult>& results) = 0;

  // Counters from the last probe() call (e.g. bucket search, misses).
  virtual Diagnostics diagnostics() const { return {}; }

  // Whether this variant is run in `order` at all (a stateful variant may
  // only be meaningful in sorted order).
  virtual bool supports(ProbeOrder) const { return true; }
};

/**
 * A Prober over one table type. The Policy adapts the table's own names:
 *
 *   using Table
 *   using State          per-probe lookup state, value-initialised at the
 *                        start of every probe (an empty struct if the
 *                        lookup is stateless)
 *   static constexpr bool sorted_only
 *                        true: only run in sorted order
 *   static KmerLookupResult lookup(const Table&, State&, std::uint32_t key)
 *   static KmerLookupResult tail(const Table&, std::uint32_t key,
 *                                std::uint32_t length)
 *   static Diagnostics counters(const Table::Stats& before,
 *                               const Table::Stats& after)
 */
template <typename Policy>
class TableProber final : public Prober {
 public:
  using Table = typename Policy::Table;

  // Loads its own table (timed as build_ms) from the table's constructor
  // arguments.
  template <typename... TableArgs>
  explicit TableProber(std::string name, TableArgs&&... args)
      : name_(std::move(name)) {
    build_ms_ = msbench::time_ms([&] {
      table_ = std::make_shared<Table>(std::forward<TableArgs>(args)...);
    });
  }

  // Shares an already loaded table (e.g. another prober's, see table()),
  // so two ways of searching the same table do not hold it twice.
  TableProber(std::string name, std::shared_ptr<Table> table)
      : name_(std::move(name)), table_(std::move(table)) {}

  const std::shared_ptr<Table>& table() const { return table_; }

  const std::string& name() const override { return name_; }

  double build_ms() const override { return build_ms_; }

  void probe(const std::vector<std::uint32_t>& keys,
             const std::vector<std::uint32_t>& order, std::uint32_t tail_key,
             std::vector<KmerLookupResult>& results) override {
    const Table& table = *table_;
    const auto before = table.stats();

    typename Policy::State state{};

    for (const std::uint32_t i : order) {
      results[i] = Policy::lookup(table, state, keys[i]);
    }

    const std::size_t n = results.size();

    for (std::size_t i = keys.size(); i < n; ++i, tail_key <<= 2U) {
      results[i] =
          Policy::tail(table, tail_key, static_cast<std::uint32_t>(n - i));
      results[i].found = false;
    }

    phase_barrier(results.data());
    diagnostics_ = Policy::counters(before, table.stats());
  }

  Diagnostics diagnostics() const override { return diagnostics_; }

  bool supports(const ProbeOrder order) const override {
    return !Policy::sorted_only || order == ProbeOrder::sorted;
  }

 private:
  std::string name_;
  double build_ms_ = 0.0;
  std::shared_ptr<Table> table_;
  Diagnostics diagnostics_;
};

// The plain v2 table (pt16_rlz_v2.hpp).
struct V2Policy {
  using Table = PT16RLZParser<Symbol, SAType>;
  struct State {};
  static constexpr bool sorted_only = false;

  static KmerLookupResult lookup(const Table& table, State&,
                                 std::uint32_t key) {
    return table.lookupKmerByKey(key);
  }

  static KmerLookupResult tail(const Table& table, std::uint32_t key,
                               std::uint32_t length) {
    return table.lookupTailByKey(key, length);
  }

  static Diagnostics counters(const Table::Stats& before,
                              const Table::Stats& after) {
    return bucket_search_diagnostics(before, after);
  }
};

// The v2 table with the cheaper miss path (variants/pt16_rlz_v2_fastmiss.hpp).
struct FastMissPolicy {
  using Table = PT16FastMissParser<Symbol, SAType>;
  struct State {};
  static constexpr bool sorted_only = false;

  static KmerLookupResult lookup(const Table& table, State&,
                                 std::uint32_t key) {
    return table.lookupKmerByKey(key);
  }

  static KmerLookupResult tail(const Table& table, std::uint32_t key,
                               std::uint32_t length) {
    return table.lookupTailByKey(key, length);
  }

  static Diagnostics counters(const Table::Stats& before,
                              const Table::Stats& after) {
    return bucket_search_diagnostics(before, after) +
           miss_diagnostics(before, after);
  }
};

// The self-contained sassy table (variants/pt16_sassy.hpp).
struct SassyPolicy {
  using Table = PT16SassyLookup;
  struct State {};
  static constexpr bool sorted_only = false;

  static KmerLookupResult lookup(const Table& table, State&,
                                 std::uint32_t key) {
    return table.lookup(key);
  }

  static KmerLookupResult tail(const Table& table, std::uint32_t key,
                               std::uint32_t length) {
    const Table::TailResult tail = table.lookup_tail(key, length);

    KmerLookupResult result;
    result.match_position = tail.match_position;
    result.match_length = tail.match_length;
    return result;
  }

  static Diagnostics counters(const Table::Stats& before,
                              const Table::Stats& after) {
    return bucket_search_diagnostics(before, after) +
           miss_diagnostics(before, after);
  }
};

// The "finger" line of a finger search: how many lookups restarted a
// bucket search, how many continued from the previous insertion point,
// and how many L entries the continuations stepped over.
template <typename Stats>
Diagnostics finger_diagnostics(const Stats& before, const Stats& after) {
  return counter_diagnostics(
      "finger",
      {{"restarts", after.finger_restarts - before.finger_restarts},
       {"continues", after.finger_continues - before.finger_continues},
       {"steps", after.finger_steps - before.finger_steps}});
}

// The fast-miss table searched with a finger
// (PT16FastMissParser::lookupKmerByKey(key, finger)): in sorted order the
// keys never decrease, so each lookup walks on from the previous one's
// insertion point in L instead of searching its bucket from scratch; only
// a new bucket restarts the search.
struct FastMissFingerPolicy {
  using Table = PT16FastMissParser<Symbol, SAType>;
  using State = Table::Finger;
  static constexpr bool sorted_only = true;

  static KmerLookupResult lookup(const Table& table, State& finger,
                                 std::uint32_t key) {
    return table.lookupKmerByKey(key, finger);
  }

  static KmerLookupResult tail(const Table& table, std::uint32_t key,
                               std::uint32_t length) {
    return table.lookupTailByKey(key, length);
  }

  static Diagnostics counters(const Table::Stats& before,
                              const Table::Stats& after) {
    return finger_diagnostics(before, after) +
           bucket_search_diagnostics(before, after) +
           miss_diagnostics(before, after);
  }
};

// The sassy table searched with a finger (PT16SassyLookup::lookup(key,
// finger)), same scheme as FastMissFingerPolicy.
struct SassyFingerPolicy {
  using Table = PT16SassyLookup;
  using State = PT16SassyLookup::Finger;
  static constexpr bool sorted_only = true;

  static KmerLookupResult lookup(const Table& table, State& finger,
                                 std::uint32_t key) {
    return table.lookup(key, finger);
  }

  static KmerLookupResult tail(const Table& table, std::uint32_t key,
                               std::uint32_t length) {
    return SassyPolicy::tail(table, key, length);
  }

  static Diagnostics counters(const Table::Stats& before,
                              const Table::Stats& after) {
    return finger_diagnostics(before, after) +
           bucket_search_diagnostics(before, after) +
           miss_diagnostics(before, after);
  }
};

// ---------- Registry ----------

struct TableBuild {
  std::string name;
  double build_ms = 0.0;
};

struct ProberSet {
  // Writing each table file, once per table format.
  std::vector<TableBuild> table_builds;

  // The variants, in report order. The first one's results are the ones
  // chained, and every other one's are checked against them.
  std::vector<std::unique_ptr<Prober>> probers;
};

/**
 * Builds each table file once and loads every variant under test.
 *
 * THIS IS THE ONLY PLACE YOU EDIT WHEN A NEW TABLE VARIANT LANDS: write a
 * Policy for it above and add one line here.
 *
 * Table files are always rebuilt (never reused from an earlier run), so a
 * stale table from another reference or format version is never read. The
 * v2 file is shared by every v2-format variant.
 */
inline ProberSet build_probers(const std::vector<Symbol>& reference,
                               const std::vector<SAType>& suffix_array,
                               const std::string& table_path) {
  ProberSet set;

  const std::string v2_path = table_path;
  const std::string sassy_path = table_path + ".sassy";

  set.table_builds.push_back({"v2 table", msbench::time_ms([&] {
                                std::filesystem::remove(v2_path);
                                build_pt16_table(reference, suffix_array,
                                                 v2_path);
                              })});

  set.table_builds.push_back({"sassy table", msbench::time_ms([&] {
                                std::filesystem::remove(sassy_path);
                                build_pt16_sassy_table(reference, suffix_array,
                                                       sassy_path);
                              })});

  set.probers.push_back(std::make_unique<TableProber<V2Policy>>(
      "pt16-v2", reference, suffix_array, v2_path));

  auto fastmiss = std::make_unique<TableProber<FastMissPolicy>>(
      "pt16-v2-fastmiss", reference, suffix_array, v2_path);
  auto fastmiss_table = fastmiss->table();
  set.probers.push_back(std::move(fastmiss));

  // Same loaded table, searched with a finger (sorted order only).
  set.probers.push_back(std::make_unique<TableProber<FastMissFingerPolicy>>(
      "pt16-v2-fastmiss-finger", std::move(fastmiss_table)));

  auto sassy = std::make_unique<TableProber<SassyPolicy>>("pt16-sassy",
                                                           sassy_path);
  auto sassy_table = sassy->table();
  set.probers.push_back(std::move(sassy));

  // Same loaded table, searched with a finger (sorted order only).
  set.probers.push_back(std::make_unique<TableProber<SassyFingerPolicy>>(
      "pt16-sassy-finger", std::move(sassy_table)));

  return set;
}

// ---------- Checking one variant's lookups ----------

// Checks one variant's `results` against the chained variant's `expected`,
// on exactly what the chain (backwardChainExtend, ms_tools.hpp) reads, so
// that chaining either list gives the same lengths:
//
//   found and match_length equal;
//   for a hit, the same count and the same occurrence set -- a singleton's
//   one occurrence is its match_position, a range's are its positions
//   (KmerLookupResult's contract: v2 also fills a one-element positions
//   span for a singleton, sassy does not, and neither is wrong);
//   every match_position a real match: reference at match_position equals
//   the input at i for match_length characters.
//
// A miss's match_position is not compared: any occurrence of the longest
// matching prefix is correct, and the variants break ties differently.
// See LookupComparison (ms_utils.hpp).

// The occurrence set of a hit, as the chain reads it.
inline std::span<const std::uint32_t> occurrences(
    const KmerLookupResult& result) {
  return result.count == 1
             ? std::span<const std::uint32_t>(&result.match_position, 1)
             : result.positions;
}

inline bool same_occurrences(std::span<const std::uint32_t> a,
                             std::span<const std::uint32_t> b) {
  if (std::equal(a.begin(), a.end(), b.begin(), b.end())) {
    return true;
  }

  // The tables only need to agree on the set, not on its order.
  std::vector<std::uint32_t> x(a.begin(), a.end());
  std::vector<std::uint32_t> y(b.begin(), b.end());
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  return x == y;
}

inline LookupComparison compare_lookups(
    const std::vector<KmerLookupResult>& expected,
    const std::vector<KmerLookupResult>& results,
    const std::vector<Symbol>& input, const std::vector<Symbol>& reference) {
  LookupComparison comparison;
  comparison.checked = true;

  const auto mismatch = [&](std::size_t i, const char* reason) {
    if (comparison.mismatches++ == 0) {
      comparison.first_mismatch = i;
      comparison.first_reason = reason;
    }
    comparison.equal = false;
  };

  if (expected.size() != results.size()) {
    mismatch(0, "different number of results");
    return comparison;
  }

  for (std::size_t i = 0; i < results.size(); ++i) {
    const KmerLookupResult& want = expected[i];
    const KmerLookupResult& got = results[i];

    if (got.found != want.found || got.match_length != want.match_length) {
      mismatch(i, "found/match_length");
      continue;
    }

    if (got.found) {
      if (got.count == 0 || (got.count > 1 && got.positions.size() != got.count)) {
        mismatch(i, "count does not match positions");
        continue;
      }

      if (got.count != want.count ||
          !same_occurrences(occurrences(got), occurrences(want))) {
        mismatch(i, "occurrences");
        continue;
      }
    }

    const std::size_t p = got.match_position;
    const std::size_t len = got.match_length;

    if (p + len > reference.size() ||
        !std::equal(input.begin() + static_cast<std::ptrdiff_t>(i),
                    input.begin() + static_cast<std::ptrdiff_t>(i + len),
                    reference.begin() + static_cast<std::ptrdiff_t>(p))) {
      mismatch(i, "match_position is not a real match");
    }
  }

  return comparison;
}

// How the 16-mer lookups classified (singleton hit / range hit / miss).
// The same for every correct variant, so taken from the chained one.
inline SearchComposition classify_lookups(
    const std::vector<KmerLookupResult>& results,
    const std::size_t kmer_positions) {
  SearchComposition composition;
  composition.available = true;

  for (std::size_t i = 0; i < kmer_positions; ++i) {
    const KmerLookupResult& result = results[i];

    if (!result.found) {
      ++composition.misses;
    } else if (result.count == 1) {
      ++composition.singleton_hits;
    } else {
      ++composition.range_hits;
    }
  }

  return composition;
}

}  // namespace msbench
