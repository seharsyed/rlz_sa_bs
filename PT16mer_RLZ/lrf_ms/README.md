# Matching-statistics benchmark (`lrf_ms/`)

This folder benchmarks matching-statistics (MS) implementations against a
common baseline: for every position of every input file, the longest match
against a reference and one place it occurs. It shares a home with an
older, unrelated, OpenMP-parallel tool (`matching_statistics.cpp`) that
computes MS its own way, over its own CLI; nothing below describes that
tool, and none of the parallelization advice at the end applies to it.

## Files

| File | What it is |
| --- | --- |
| `ms_main.cpp` | The benchmark's entry point (build it, run it). |
| `ms_variants.hpp` | The registry: every implementation under test, one line each. **This is the only file you edit to add a variant.** |
| `ms_utils.hpp` | Everything generic: CLI args, loading, timing, invariant/position/brute-force checking, CSV and checksum output. |
| `lrf_ms.hpp` | The baseline: `LRFMS`, the classical O(n+m) suffix-array + LCP + RMQ matching-statistics algorithm. |
| `pt16_sassy_ms.hpp` | The sassy-PT16-table-based variants (`PT16SassyMS` and its scan/bucket-scan adapters). |
| `rmq_tree.h` | The RMQ structure `LRFMS` builds over the LCP array. |
| `ms_test.cpp` | Standalone correctness test — builds its own SA, checks `LRFMS` against brute force, exercises the `ms_utils` comparison/digest helpers. Independent of `ms_main`. |

## Quick start

Same build convention as the rest of `PT16mer_RLZ/` — plain `g++`, no CMake here:

```bash
cd PT16mer_RLZ
g++ -std=c++20 -O3 lrf_ms/ms_main.cpp -o ms_main
./ms_main \
  --reference path/to/reference.cleaned \
  --suffix-array path/to/reference.cleaned.sa \
  --filenames path/to/input_list.txt \
  --results results/ms_run.csv
```

`--filenames` is a text file with one input path per line. The reference and
suffix array must match exactly: the suffix array is built over the
reference **without** a sentinel, and the loader checks `suffix_array.size()
== reference.size()` (throws otherwise). `--results` is a CSV, one row per
(file, implementation); a `--results.checksums.csv` and a summary section are
written alongside it.

Run the correctness test the same way:

```bash
g++ -std=c++20 -O2 lrf_ms/ms_test.cpp -o ms_test && ./ms_test
```

## CLI reference

```
--reference PATH --suffix-array PATH --filenames PATH --results PATH   (required)
[--checksums PATH]        default: <results>.checksums.csv
[--pt16-table PATH]       default: "<dataset>_pt16_hl.bin" next to the reference
[--dump-dir DIR]          write every implementation's raw MS output per file
[--repeats N]             timed runs per file, the minimum is reported
[--no-invariants]         skip the structural checks (on by default)
[--verify-full]           check every position's occurrence, not a sample
[--sample N --seed S]     check N random positions instead (default 10000)
[--verify-maximality]     also prove the reported length cannot be extended
[--verify-brute]          compare against an O(|input|*|reference|) brute force (off by default)
[--brute-limit N]         size limit that gates --verify-brute (default 20000)
[--stop-on-mismatch]      stop after the first file whose lengths diverge from baseline
```

`--no-invariants` and skipping `--verify-brute` (its default) are what keep a
run fast when you only care about timing, not correctness — that's the mode
used throughout this session's scan-only experiments.

## Registered implementations

Registration order is exit order; index 0 is always the baseline everything
else is compared and speedup-ratio'd against.

| Name | What it computes | Table format |
| --- | --- | --- |
| `lrf-ms` | Baseline. Classical algorithm: suffix range narrowing plus an LRF-array skip that resolves most consecutive positions in O(1), no table lookup at all. | none (ISA/LCP/LRF/RMQ, built once) |
| `pt16-v2-bucket-chain-multi` | One 16-mer lookup per position, run grouped by table bucket, then multi-step chain extension. | non-sassy H/L (`pt16_build_v2.hpp`) |
| `pt16-v2-sorted-chain-multi` | Same, with the lookups fully sorted by their 32-bit key (`sortedKmerScan`). | same |
| `pt16-v2-fastmiss-sorted-chain-multi` | As `pt16-v2-sorted-chain-multi`, with the cheaper miss path of `PT16FastMissParser`. | v2 table, own copy (`<table>.fastmiss`) |
| `pt16-v2-fastmiss-bucket-chain-multi` | As `pt16-v2-bucket-chain-multi`, with the fast-miss lookup. | same as above |
| `pt16-sassy-chain-multi` | Bucketed lookups over the self-contained sassy table, then multi-step chain extension. | sassy (`pt16_build_sassy.hpp`, `<table>.sassy`) |
| `pt16-sassy-sorted-chain-multi` | Same, with fully sorted lookups. | same |

Only variants that compute the full, exact matching statistics are
registered: every row above is proven exact against brute force
(`chain_extend_test.cpp`), so it must report `lengths=N/N equal`. The
scan-only rows (raw 16-mer lookups capped at 16) and the one-step
`pt16-sassy-chain` were dropped from the timed runs; their classes remain in
`ms_variants.hpp` / `pt16_sassy_ms.hpp` / `pt16_fastmiss_ms.hpp`. What was
found with the scan-only rows (not exhaustively, and not yet on genome-scale
real data at the time of writing):

- Bucketing consistently beats input-order scanning for both table formats,
  often enough to flip a variant from slower-than-baseline to faster.
- The plain (non-sassy) table has so far beaten the sassy table on raw
  lookup speed, despite needing an extra suffix-array read the sassy table
  avoids — most likely because its L entry is half the size (4 vs. 8
  bytes), so more entries share a cache line during a bucket search.
- Whether `H_interleaved_` (directory + SA-start packed together) beats two
  separate arrays depends on how often a lookup has to walk past several
  empty buckets to find where the current one ends — that walk only needs
  half of what interleaving fetches. This has been observed to go either
  way depending on the dataset; it isn't settled.

Treat all of the above as leads, not conclusions — see the diagnostics below
for how to check them against your own data.

## Reading the diagnostics

Nothing is printed per implementation per file (the per-file rows go to
the results CSV). Each implementation's `diagnostics()` is collected after
every `compute()` call, summed over all files, and printed once per
implementation under `[10] DIAGNOSTICS`, after the collection totals:

```
pt16-sassy-chain-multi
        phases: prebucket 0.215 ms (5.1%)  bucket 0.422 ms (9.9%)  probe 3.620 ms (85.0%)  tail 0.003 ms (0.1%)  total 4.260 ms
        bucket search: linear=149792 binary=0
```

- **`phases`** (bucketed variants; the sorted ones report `keys`,
  `radix-low`, `radix-high`, `probe`, `tail` instead): the scan's own four stages —
  `prebucket` (pack every 16-mer, tally per-bucket counts), `bucket`
  (prefix-sum + counting-sort placement), `probe` (the lookups themselves,
  now in bucket order — the phase bucketing exists to speed up), `tail`
  (the last <16 characters, never reordered). Timed coarsely — 4
  `clock::now()` pairs per whole file, not per lookup — so it doesn't
  measurably inflate the numbers it reports on.
- **`bucket search`**: how many of this call's non-empty-bucket dispatches
  searched their bucket linearly vs. with `std::lower_bound` (the switch is
  at 64 entries). Free counters (plain increments already on a taken
  branch), reported by the bucketed v2 and sassy variants (fastmiss
  reports `misses` instead). Useful for checking whether a format's average bucket size sits mostly below or
  above that threshold on your data — the answer changes which cost model
  (bandwidth-bound linear scan vs. probe-count-bound binary search)
  actually applies.

## Adding a new variant

Edit only `ms_variants.hpp`. A class qualifies if it has:

```cpp
Impl(const std::vector<Symbol>&, const std::vector<SAType>&, /* extra args */);
MatchingStatistics computeMatchingStatistics(const std::vector<Symbol>&);
```

then one line in `build_implementations`:

```cpp
add_implementation<MyVariant<Symbol, SAType>>(implementations, "my-variant",
                                              reference, suffix_array /*, extra args */);
```

The constructor is what gets timed into `build_ms()`. Optionally add
`Diagnostics diagnostics() const` (lines of phase times or counters, built
with `phase_diagnostics` / `counter_diagnostics` from `pt16_utils.hpp`) and
`MSAdapter` will forward it automatically via `if constexpr` — no interface
class to inherit from. `ms_main` sums it over every file and prints it once
per variant under `[10] DIAGNOSTICS`.

## Parallelizing over inputs

**Where**: the per-file loop in `ms_main.cpp`, section `[8] MATCHING
STATISTICS` — one iteration per entry in `--filenames`. Each file's MS
computation is independent of every other file's: different input, own
`MatchingStatistics` output, nothing computed for file *i* is needed for
file *j*. That makes the file loop the natural (and about the only sane)
place to parallelize — not inside a single file's scan, and not the table
build (`[6] BUILD`, which already runs once, up front, before any file is
touched).

### What's already safe to share across threads

Everything built once in `[6] BUILD` and never written to afterwards:
`LRFMS`'s `isa_`/`lcp_`/`lrf_`/`rmq_`, and each PT16 table's `H_`/`L_`/
`sampled_sa_`/`short_suffixes_` arrays. `LRFMS::computeMatchingStatistics`
is `const` and touches none of its own mutable state, so **it is already
safe to call concurrently from multiple threads sharing one `LRFMS`
instance** — no changes needed there.

### What is *not* safe as this code stands today

1. **`PT16RLZParser::Stats` / `PT16SassyLookup::Stats`.** Both are
   `mutable Stats stats_;`, updated with plain `++stats_.hits` (etc.)
   inside otherwise-`const` lookup methods. Two threads calling `lookup()`
   concurrently on the *same* instance race on these counters — undefined
   behavior, not just an inaccuracy. Worse, even fixing the race (e.g. by
   making the counters atomic) wouldn't make the **diagnostics** correct:
   every PT16 variant reports a per-call delta (`stats() before` vs. `after`
   around one file's scan), and if another thread's lookups land on the
   same shared instance in between, that delta silently includes their
   work too, not just this file's.
2. **`ImplementationTotals::accumulate`** — the running per-implementation
   totals (`totals[k]`) are mutated once per file per implementation with
   no synchronization at all. Concurrent files calling `accumulate` for the
   same `k` race on plain `double`/`std::size_t` fields.
3. **`CSVWriter`/`ChecksumWriter`** — both write one row per call directly
   to an `std::ofstream`; concurrent `write_row`/`write` calls need to be
   serialized (a mutex, or route every write through one dedicated writer
   thread/queue).
4. **`--stop-on-mismatch`** stops at "the first divergent file" — a
   well-defined idea only when files run in order. Under concurrency,
   "first" becomes whichever thread happens to finish first, which will
   vary run to run. Worth deciding up front whether that's acceptable or
   whether the flag should be disabled/reinterpreted ("stop launching new
   files once *any* divergence is seen so far") when running in parallel.

### Recommended shape

Give **each worker thread its own instance** of every implementation,
rather than sharing one across threads. This sidesteps every hazard above
at once — no shared mutable `stats_`, correct per-file diagnostics, nothing
to make atomic — at the cost of one extra copy of each implementation's
in-memory table per thread.

That's cheap for `LRFMS` (rebuild is fast, and it's already thread-safe to
share besides). For the PT16 wrappers it needs one small change: their
constructors currently always rebuild the table file from scratch
(`PT16SassyMS`'s constructor does `fs::remove(table_path);
build_pt16_sassy_table(...)`, `PT16MS`/`PT16ScanMS` likewise for the
non-sassy table) — fine once, wrong N times if N threads each try to
remove-and-rebuild the same file concurrently. Build the table once up
front (as today), then give each thread a lighter constructor path that
only *loads* the already-built file (`PT16SassyLookup(path)` directly, or
the equivalent `PT16RLZParser` constructor) without touching it on disk.

Sketch (illustrative, not wired into `ms_main.cpp`):

```cpp
// Built once, before any thread starts (unchanged from today).
msbench::Implementations implementations = msbench::build_implementations(...);

// One set of per-thread instances, loaded (not rebuilt) from the files
// implementations[] already built above.
struct WorkerState {
  LRFMS<Symbol, SAType> lrf_ms;               // safe to share; kept here for symmetry
  PT16SassyLookup sassy_scan;                 // PT16SassyLookup(table_path), load-only
  // ... one load-only instance per PT16 variant this thread will run ...
};

std::vector<std::future<std::vector<FileRunResult>>> futures;
std::mutex output_mutex;   // guards csv/checksums; totals[] accumulated after join

for (auto& chunk : split(files, thread_count)) {
  futures.push_back(std::async(std::launch::async, [&, chunk] {
    WorkerState state(reference, suffix_array, pt16_table);  // load-only ctors
    std::vector<FileRunResult> results;

    for (const std::string& filename : chunk) {
      const auto input = msbench::load_input<Symbol>(filename);
      // ... run every implementation via `state`, build one FileRunResult
      //     per implementation, same logic ms_main.cpp already has ...
      results.push_back(/* ... */);
    }

    return results;   // returned, not written -- avoids the writer race entirely
  }));
}

// Join in file order, so CSV/checksums/totals[] stay exactly as
// deterministic as the sequential run, just computed concurrently.
for (auto& future : futures) {
  for (FileRunResult& result : future.get()) {
    csv.write_row(result);
    checksums.write(result);
    totals[/* k */].accumulate(result);   // single-threaded here, so no race
  }
}
```

The key property that makes this simple: **do the writing (CSV, checksums,
`totals[]`) after joining, single-threaded, in file order** — every thread
only *returns* its `FileRunResult`s, it never touches shared output state
itself. That's what removes hazards 2 and 3 above without needing any
locking in the hot path; hazard 1 is removed by not sharing instances;
hazard 4 needs an explicit decision (see above) since it's a behavior
question, not a data race.

`file_index`-based chunking (contiguous ranges, as in the sketch) keeps
memory bounded and predictable — worker `t` never holds more than one
input's worth of data at a time, same as the sequential run today, just
`thread_count` times over. A work-stealing/queue-based split would balance
better if input file sizes vary a lot, at the cost of a bit more plumbing.

Peak memory scales linearly with thread count under this design (each
thread holds its own copy of every PT16 table plus whatever
`MatchingStatistics` it's currently computing) — worth checking against
`peak_RSS_MB` from a single-threaded run times your intended thread count
before running this against a real genome-scale reference.
