# ms2 — new matching-statistics methods

New methods, kept apart from `../lrf_ms/` (which they reuse only for the
loaders and the `LRFMS` baseline). Built up incrementally.

## Build and run

From `PT16mer_RLZ/`:

```
g++ -std=c++20 -O3 ms2/ms2_main.cpp -o ms2_main
./ms2_main --reference REF --suffix-array REF.sa --filenames LIST \
           [--results CSV] [--table PATH] [--max-files N] [--max-bytes N] \
           [--no-check] [--no-chain]

g++ -std=c++20 -O2 ms2/input_keys_test.cpp -o input_keys_test && ./input_keys_test
```

`--max-files` / `--max-bytes` cut a run down (the first N files, the first
N bytes of each) for quick experiments. `--table` sets where the PT16
table files are written (default `<reference>.ms2_pt16` and `.sassy`;
always rebuilt). `--no-chain` skips the chain, which is still far too slow
on real-sized inputs.

## What it does so far

Per input file (non-ACGT bytes are separators, as in `ms_main`):

| Step | |
| --- | --- |
| `lrf-ms` | the baseline's full matching statistics (time, average match length) |
| `keys` | one pass (`input_keys.hpp`, `prepare_keys`): packed `key << 32 \| position` values in text order for every position with a full ACGT window — **except all-A (key 0) and all-T (key `0xFFFFFFFF`) 16-mers**, recorded as runs of consecutive positions instead and not queried yet; plus the short queries (the up to 15 positions before a separator or the end, 1–15 ACGT characters each). Separator positions (non-ACGT bytes) are only counted: their match is 0, so they are never searched |
| `bucket-order` | stable counting sort of the keys by bucket (high 16 bits); phases `alloc count prefix scatter` |
| `sorted-order` | stable LSD radix sort by key (11/11/10-bit passes); phases `alloc count pass-1 pass-2 pass-3` |
| `<variant>-<order>` | the table lookups (`probe_order.hpp`) over the keys in `bucket` and in `sorted` order, into one array of lookup results by text position; short queries get tail lookups; separator and run positions stay empty. Variants: `pt16-v2` (only on an all-ACGT reference), `pt16-sassy`, `pt16-sassy-finger` (sorted only). Phases `reset keys short`. Every row's results are checked against the first row's (`compare_lookups`) |
| `chain` | the backward chain (`chain.hpp`, a copy of `backwardChainExtend` with progress marks every 5%) on the first row's results, checked against lrf-ms: equal lengths, mismatches at run positions, mismatches whose true match reaches a run (both expected until runs are looked up), and unexplained ones (a bug); plus a sampled check that reported positions really match |

Every position is exactly one of: a key, a run position, a short query,
a separator.
Both orders keep ties in text order. Each order is checked after it is
timed (same values, correctly ordered); `--no-check` skips that.

`alloc` is growing the output buffers; they are reused across files, so
only the first file (or a larger one) pays it.

Output: live per-file lines on stderr, totals at the end (including each
row's pipeline time: keys + its order + its lookups + the chain, against
lrf-ms), `key=value` lines on stdout, and a per-file CSV with `--results`.

The table side (lookup policies, table builders, `compare_lookups`) is
reused from `../lrf_ms/probe_pipeline.hpp` unchanged.
