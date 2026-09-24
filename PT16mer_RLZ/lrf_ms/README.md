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
| `probe_pipeline.hpp` | The shared PT16 pipeline: rolling keys, the bucket and sorted probe orders, the `Prober` interface, one `Policy` per table variant, and **`build_probers`, the registry of PT16 variants**. Also `compare_lookups`, the check between variants. |
| `ms_variants.hpp` | The registry of *full* implementations (`build_implementations`) — ones that compute complete MS on their own, currently only the baseline. Also the older self-contained PT16 variant classes, kept for the tests. |
| `ms_tools.hpp` | `backwardChainExtend`: multi-step chain extension, lookup results to exact MS. |
| `chain_extend.hpp` | `chainExtend`: the one-step version (not exact; kept for comparison). |
| `ms_utils.hpp` | Everything generic: CLI args, loading, timing, invariant/position/brute-force checking, result rows and totals, CSV and checksum output. |
| `lrf_ms.hpp` | The baseline: `LRFMS`, the classical O(n+m) suffix-array + LCP + RMQ matching-statistics algorithm. |
| `pt16_sassy_ms.hpp`, `pt16_fastmiss_ms.hpp`, `sorted_kmer_scan.hpp` | The older self-contained sassy / fast-miss MS classes and the sorted scan they share; no longer run by `ms_main`, used by the tests. |
| `rmq_tree.h` | The RMQ structure `LRFMS` builds over the LCP array. |
| `ms_test.cpp` | Standalone correctness test — builds its own SA, checks `LRFMS` against brute force, exercises the `ms_utils` comparison/digest helpers. Independent of `ms_main`. |
| `probe_pipeline_test.cpp` | Standalone test of the pipeline: every variant, both orders, on small and repetitive references — lookups agree (`compare_lookups`), and chaining each variant's own results gives brute-force lengths. |
| `chain_extend_test.cpp`, `fastmiss_test.cpp` | Tests of the chain extensions and of the fast-miss parser against brute force. |

The table variants themselves live in `../variants/` (fast-miss, sassy,
interleaved v2) and `../pt16_rlz_v2.hpp` (plain v2).

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
g++ -std=c++20 -O2 lrf_ms/probe_pipeline_test.cpp -o probe_pipeline_test && ./probe_pipeline_test
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

Two kinds, run differently.

**Full implementations** (`build_implementations` in `ms_variants.hpp`)
each compute complete matching statistics on their own. Registration order
is report order; index 0 is always the baseline everything else is compared
and speedup-ratio'd against.

| Name | What it computes |
| --- | --- |
| `lrf-ms` | Baseline. Classical algorithm: suffix range narrowing plus an LRF-array skip that resolves most consecutive positions in O(1), no table lookup at all. Builds ISA/LCP/LRF/RMQ once. |

**PT16 table variants** (`build_probers` in `probe_pipeline.hpp`) share
one pipeline, run per file as:

| Stage | Runs | What it does |
| --- | --- | --- |
| `keys` | once | roll every 16-mer key of the input |
| `bucket-order` | once | positions grouped by table bucket (one counting-sort pass) |
| `sorted-order` | once | positions fully sorted by key (two LSD counting-sort passes) |
| probe | per variant, per order | one table lookup per key in that order, plus the tail positions |
| `chain` | once | `backwardChainExtend` on the first variant's lookup results |

Only the probe depends on the variant: the keys and orders depend only on
the input, and every correct variant returns the same lookup results (each
stored at its own text position, whatever order it was probed in), so the
chain runs once. Each probe row is `<variant>-<order>`, and its *pipeline*
time — keys + its order + its probe + the chain — is what its speedup is
computed from.

| Variant | Table |
| --- | --- |
| `pt16-v2` | plain v2 (`../pt16_rlz_v2.hpp`), file built by `../pt16_build_v2.hpp`. **Currently left out of `build_probers`** to keep runs short (commented out there). |
| `pt16-v2-fastmiss` | fast-miss (`../variants/pt16_rlz_v2_fastmiss.hpp`): same v2 file, reference position per entry, precomputed empty buckets, flagged short-suffix buckets |
| `pt16-v2-fastmiss-finger` | the same loaded fast-miss table, searched with a finger (`PT16FastMissParser::lookupKmerByKey(key, finger)`); **sorted order only**, same scheme as `pt16-sassy-finger` below. |
| `pt16-sassy` | sassy (`../variants/pt16_sassy.hpp`), `<table>.sassy`: self-contained, no reference or SA reads |
| `pt16-sassy-finger` | the same loaded sassy table, searched with a finger (`PT16SassyLookup::lookup(key, finger)`); **sorted order only**. Keys never decrease there, so each lookup walks on in `L` from the previous insertion point; only a new bucket restarts the search. Reports `finger: restarts / continues / steps`. |

Each table file is written once per run (the v2 file is shared by `pt16-v2`
and `pt16-v2-fastmiss`) and every variant only loads it; `[6] BUILD`
reports the writes and loads separately.

**Correctness.** The chain's output (`pt16-chain`) gets the full checks —
invariants, positions, lengths against the baseline (`lengths N/N equal`).
Every probe's lookups are compared against the chained ones with
`compare_lookups` (`lookups N/N equal`), on exactly what the chain reads:
`found`, `match_length`, `count`, and the occurrence set (a singleton's
`match_position`, a range's `positions`), plus that every `match_position`
is a real match. A miss's `match_position` is not compared: any occurrence
of the longest matching prefix is correct. `probe_pipeline_test.cpp` checks
that agreeing results really do chain to the same (brute-force) lengths.

What was found with the older per-variant scans (not exhaustively, and not
yet on genome-scale real data at the time of writing):

## Reading the diagnostics

Nothing is printed per file except failures (the per-file rows, one per
full implementation, stage, probe and chain, go to the results CSV with a
`kind` column). `[9] COLLECTION TOTALS` has the full implementations, the
shared stages with their share of the shared time, and one row per variant
and order: load time, probe time, pipeline time, MB/s, speedup, lookup
check. `[10] DIAGNOSTICS` has each variant's counters from its probes,
summed over files:

```
pt16-sassy-bucket
        bucket search: linear=149792 binary=0
        misses: total=17638 empty-bucket=508 full-short-suffix-checks=1
16-mer lookups: singleton=88.3%  range=0.0%  miss=11.7%
```

- **`bucket search`**: how many non-empty-bucket dispatches searched their
  bucket linearly vs. with `std::lower_bound` (the switch is at 64
  entries). Useful for checking whether a format's average bucket size sits
  mostly below or above that threshold on your data — the answer changes
  which cost model (bandwidth-bound linear scan vs. probe-count-bound
  binary search) actually applies.
- **`misses`** (fast-miss and sassy): `total` misses, how many landed in an
  `empty-bucket` (answered from arrays precomputed at load time), and how
  many ran the `full-short-suffix-checks` because their bucket holds a
  short suffix longer than 8 characters. Every other miss skips the short
  suffixes entirely.
- **`16-mer lookups`**: how the lookups classified, from the chained
  results (the same for every correct variant).

## Adding a new variant

**A PT16 table variant**: in `probe_pipeline.hpp`, write a `Policy` for it
(its table type and how to call its 16-mer lookup, its tail lookup and its
counters — see `V2Policy`) and add one line to `build_probers`:

```cpp
set.probers.push_back(std::make_unique<TableProber<MyPolicy>>(
    "my-variant", /* table constructor args */));
```

It is then probed in both orders, timed, and checked against the chained
results automatically. The table constructor is what gets timed as its
load; if it needs a new file format, build the file once in
`build_probers` and add it to `table_builds`.

**A full implementation** (computes complete MS on its own): edit
`ms_variants.hpp`. A class qualifies if it has

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
`Diagnostics diagnostics() const` and `MSAdapter` forwards it; `ms_main`
sums it over files and prints it under `[10] DIAGNOSTICS`.

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
