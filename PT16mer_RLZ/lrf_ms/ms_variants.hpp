#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../pt16_build_v2.hpp"
#include "../pt16_rlz_v2_interleaved.hpp"
#include "lrf_ms.hpp"
#include "ms_utils.hpp"

/**
 * The registry of matching-statistics implementations under test.
 *
 * THIS IS THE ONLY FILE YOU EDIT WHEN A NEW VARIANT LANDS.
 *
 * An implementation qualifies if it has
 *
 *   Impl(const std::vector<Symbol>&, const std::vector<SAType>&, ...)
 *   MatchingStatistics computeMatchingStatistics(const std::vector<Symbol>&)
 *
 * Extra constructor arguments (a table path, a block size, a thread
 * count) are forwarded, so a variant that needs them does not need a
 * different adapter.
 *
 * The first registered implementation is the baseline: it is what every
 * other implementation's lengths are compared against, and what the
 * speedup column is relative to.
 */

namespace msbench {

using Symbol = unsigned char;
using SAType = std::uint32_t;

// ---------- Interface ----------

class MSImplementation {
 public:
  virtual ~MSImplementation() = default;

  virtual const std::string& name() const = 0;

  // Wall time of this implementation's construction, in milliseconds.
  virtual double build_ms() const = 0;

  virtual MatchingStatistics compute(const std::vector<Symbol>& input) = 0;
};

// ---------- Adapter ----------

/**
 * Wraps a concrete index type. The constructor is what gets timed, so
 * whatever preprocessing the variant does lands in build_ms().
 *
 * The index is held by unique_ptr because these classes are typically
 * neither copyable nor movable.
 */
template <typename Index>
class MSAdapter : public MSImplementation {
 public:
  template <typename... CtorArgs>
  MSAdapter(std::string name, CtorArgs&&... args) : name_(std::move(name)) {
    build_ms_ = time_ms([&] {
      index_ = std::make_unique<Index>(std::forward<CtorArgs>(args)...);
    });
  }

  const std::string& name() const override { return name_; }

  double build_ms() const override { return build_ms_; }

  MatchingStatistics compute(const std::vector<Symbol>& input) override {
    return index_->computeMatchingStatistics(input);
  }

 private:
  std::string name_;
  double build_ms_ = 0.0;
  std::unique_ptr<Index> index_;
};

using Implementations = std::vector<std::unique_ptr<MSImplementation>>;

template <typename Index, typename... CtorArgs>
inline void add_implementation(Implementations& implementations,
                               std::string name, CtorArgs&&... args) {
  implementations.push_back(std::make_unique<MSAdapter<Index>>(
      std::move(name), std::forward<CtorArgs>(args)...));
}

// ---------- PT16 brute-force matching statistics ----------

/**
 * PT16 matching statistics: one PT16 longest-match query per input position
 * (PT16RLZParser::computeMS_brute).
 *
 * The constructor rebuilds the PT16 table at `table_path` and loads it, so
 * build_ms() covers both the table construction and the load.
 */
template <typename T1, typename T2>
class PT16MS {
 public:
  PT16MS(const std::vector<T1>& reference, const std::vector<T2>& suffix_array,
         const std::string& table_path) {
    // Always rebuild, so a stale table from another reference or another
    // version of the table format is never used.
    fs::remove(table_path);

    build_pt16_table(reference, suffix_array, table_path);

    parser_ = std::make_unique<PT16RLZParser<T1, T2>>(reference, suffix_array,
                                                      table_path);
  }

  PT16MS(const PT16MS&) = delete;
  PT16MS& operator=(const PT16MS&) = delete;

  MatchingStatistics computeMatchingStatistics(const std::vector<T1>& input) {
    return parser_->computeMS_brute(input);
  }

 private:
  std::unique_ptr<PT16RLZParser<T1, T2>> parser_;
};

// ---------- Registry ----------

/**
 * Builds every implementation under test, in order. Index 0 is the
 * baseline.
 *
 * To add a variant, include its header above and add one line here:
 *
 *   add_implementation<MyVariant<Symbol, SAType>>(
 *       implementations, "my-variant", reference, suffix_array);
 *
 * Constructor arguments beyond the reference and suffix array are
 * forwarded as given, for example:
 *
 *   add_implementation<PT16MS<Symbol, SAType>>(
 *       implementations, "pt16-ms", reference, suffix_array, table_path);
 */
inline Implementations build_implementations(
    const std::vector<Symbol>& reference,
    const std::vector<SAType>& suffix_array, const std::string& pt16_table) {
  Implementations implementations;

  add_implementation<LRFMS<Symbol, SAType>>(implementations, "lrf-ms",
                                            reference, suffix_array);

  // --- Add variants below this line ---

  add_implementation<PT16MS<Symbol, SAType>>(implementations, "pt16-brute",
                                             reference, suffix_array,
                                             pt16_table);

  return implementations;
}

}  // namespace msbench
