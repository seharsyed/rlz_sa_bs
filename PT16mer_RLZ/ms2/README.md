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
on real-sized inputs. `--quiet` drops the per-file lines (a single
`files: i/N` counter remains) and prints only the totals. `--keep-runs`
(`ms2_main`, `unique_main`) looks up all-A / all-T 16-mers like any other
key instead of setting them aside as runs.

## What it does so far

Per input file (non-ACGT bytes are separators, as in `ms_main`):

| Step | |
| --- | --- |
| `lrf-ms` | the baseline's full matching statistics (time, average match length) |
| `keys` | one pass (`input_keys.hpp`, `prepare_keys`): packed `key << 32 \| position` values in text order for every position with a full ACGT window — **except all-A (key 0) and all-T (key `0xFFFFFFFF`) 16-mers**, recorded as runs of consecutive positions instead and not queried yet; plus the short queries (the up to 15 positions before a separator or the end, 1–15 ACGT characters each). Separator positions (non-ACGT bytes) are only counted: their match is 0, so they are never searched |
| `bucket-order` | stable counting sort of the keys by bucket (high 16 bits); phases `alloc count prefix scatter` |
| `sorted-order` | stable LSD radix sort by key (11/11/10-bit passes); phases `alloc count pass-1 pass-2 pass-3` |
| `<variant>-<order>` | the table lookups (`probe_order.hpp`) over the keys in `bucket` and in `sorted` order, into one array of lookup results by text position; short queries get tail lookups; separator and run positions stay empty. Variants: `pt16-v2` and `pt16-v2-finger` (only on an all-ACGT reference), `pt16-sassy`, `pt16-sassy-finger` (the finger ones on sorted order only). Phases `reset keys short`. Every row's results are checked against the first row's (`compare_lookups`) |
| `chain` | the backward chain (`chain.hpp`, a copy of `backwardChainExtend` that can report progress; `ms2_main` does not print it) on the first row's results, checked against lrf-ms: equal lengths, mismatches at run positions, mismatches whose true match reaches a run (both expected until runs are looked up), and unexplained ones (a bug); plus a sampled check that reported positions really match |

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

## Lazy matching statistics (`lazy_ms.hpp`, `lazy_main.cpp`)

Backward left-extension with PT16 lookups only where needed — the
reference point against lrf-ms. Walking the input backward, with `(q, L)`
the answer at `i + 1`:

- **extend:** `MS[i] <= MS[i+1] + 1` always, so if `ref[q-1] == input[i]`
  then `MS[i] = (q-1, L+1)` exactly, from one comparison;
- **break:** otherwise one table lookup at `i`: a tail (fewer than 16 ACGT
  characters left) or a miss is the answer; a singleton is extended from
  its position; a range is narrowed among its occurrences (binary search,
  suffix order) and then extended;
- **separator:** length 0.

No bucketing, sorting or chain: on an input close to the reference almost
every position extends. Works with any table through its lookup policy
(`lazy-v2`, `lazy-sassy`).

```
g++ -std=c++20 -O3 ms2/lazy_main.cpp -o lazy_main
./lazy_main --reference REF --suffix-array REF.sa --filenames LIST \
            [--results CSV] [--table PATH] [--max-files N] [--max-bytes N] \
            [--sample N] [--quiet]

g++ -std=c++20 -O2 ms2/lazy_ms_test.cpp -o lazy_ms_test && ./lazy_ms_test
```

Per file: lrf-ms, then each lazy variant, with its time, speedup over
lrf-ms, and whether its lengths equal lrf-ms's at every position (plus
`--sample` positions checked to really match). `--quiet` drops the
per-file lines (a `files: i/N` counter remains) and prints only the totals. Counters: positions
extended / breaks / separators; how breaks resolved (tails, misses,
singletons, ranges); and the work done at breaks (characters extended past
16, occurrences of the ranges met, characters narrowed).

## Unique 16-mer lookups (experimental: `unique_keys.hpp`, `unique_main.cpp`)

Only the input's 16-mers that occur **once in the input** are looked up.
Per file: `keys` (the key pass), `sorted-order` (LSD radix), `unique` (one
sequential pass over the sorted list keeps keys whose neighbours differ:
equal keys are adjacent only after sorting, so this is the cheapest place
to filter), `results` (the lookup-result arrays sized once per file), then
the lookups of the sorted unique keys into a compact array aligned with
the unique list, per variant: `pt16-v2`, `pt16-v2-finger`, `pt16-sassy`,
`pt16-sassy-finger`. Every variant's results are checked against the first
one's (same found/length/count/occurrence set), and the first one's
positions against the input. lrf-ms is timed for comparison.

Then the chain (`unique_chain.hpp`), on the first variant's results:

- `text-order`: the unique keys back in input order, as a radix sort of
  `position << 32 | k` (k indexes the compact results) — no input-sized
  array and no random write;
- `chain`: one right-to-left sweep with a scratchbook holding the last entry
  seen (input position P, reference position q, length L, exact or not).
  For the entry at i, d = P − i: a miss is exact; a found 16-mer with
  d ≤ 16 whose occurrences include q − d takes (q − d, L + d) — exact if
  the previous entry was exact (MS[i] ≤ MS[i+d] + d), else an approximate
  lower bound; otherwise (gap d > 16, or no occurrence lines up) it gets a
  lower bound of 16 and is marked approximate. Checked against lrf-ms at
  every unique position: exact lengths equal, approximate ones not above;
  also counts how many approximate ones happen to be right. The lookup
  composition (singletons / ranges / misses) is printed once per file.

```
g++ -std=c++20 -O3 ms2/unique_main.cpp -o unique_main
./unique_main --reference REF --suffix-array REF.sa --filenames LIST \
              [--results CSV] [--table PATH] [--max-files N] [--max-bytes N] [--quiet]

g++ -std=c++20 -O2 ms2/unique_keys_test.cpp -o unique_keys_test && ./unique_keys_test
```

The v2 table now has a finger search too (`PT16RLZParser::lookupKmerByKey(key,
finger)`, policy `V2FingerPolicy`), the same scheme as sassy's: continue
from the previous insertion point in the same bucket, and on a new bucket
walk from the bucket's start.

## Sampled method (experimental: `sampled.hpp`, run at the end of `ms2_main`)

After each file's pipeline, per table (`sampled-v2`, `sampled-sassy`),
with step L (`--sample-step`, default 16):

1. **lookups:** one lookup at every L-th position into an array of about
   n / L results (a 16-mer lookup; a tail lookup where fewer than 16 ACGT
   characters are left; nothing at a separator). `sampled-v2-textorder` /
   `sampled-sassy-textorder` look them up in text order; `sampled-v2-sorted` /
   `sampled-sassy-sorted` pack the full windows' keys as `key << 32 | j`,
   radix sort them and look them up with the finger in sorted order
   (phases `collect sort lookups`), writing each result back to its sample;
2. **link:** right to left, a hit links to the next sample when one of its
   occurrences is exactly L before that sample's chosen reference position
   (a range keeps only that occurrence; an unlinked range stays undecided
   and may be decided by the sample to its left). Misses and tails never
   link. Maximal runs of linked samples are **stretches** (one alignment,
   no lookups needed inside);
3. **midpoints:** one 16-mer lookup at the middle of every hole (unlinked
   neighbouring samples), recording whether it lines up with the left
   sample, the right one, both or neither. Nothing is inferred yet. Holes
   between two separator samples (inside a run of non-ACGT bytes) are
   only counted (`separator-holes`); the midpoints always use plain
   lookups.

Each step is timed; every row's stretches are checked (their alignment
matches the reference, lrf-ms's match at each stretch start is at least as
long, and every row finds the same links as the first).

## Sampled method end to end (experimental: `sampled_ms.hpp`, rows `sampled-ms-<table>` in `ms2_main`)

Full, exact matching statistics with lookups only where they are needed:

1. **samples:** every L-th position (default 16), keys sorted, finger
   lookups; neighbouring samples linked into chains. Each chain gives a
   strict stretch, from its first singleton sample to its last one: the
   input matches the reference along one path between them (ranges in
   between take the path's occurrence; ranges at a chain's ends are left
   outside), and its interior needs no lookups;
2. **holes:** every position outside the stretches (not a sample, not a
   separator) gets its key rolled through the hole, one character per
   position (a stretch is skipped, the rolling restarts after it); keys
   sorted, finger lookups; positions with fewer than 16 ACGT characters
   left get a tail lookup. So every position outside a stretch has a
   lookup result;
3. **matching statistics,** one pass right to left, from the lookup
   results only: a position's lookup result is chained to i + 1 when one of
   its occurrences is q − 1 (then (q − 1, L + 1)), else resolved from the
   result alone (short phrase as is, singleton extended, range narrowed).
   At a stretch's last singleton the same; then the whole interior is
   filled in one loop, each position the previous answer + 1 at path − 1 —
   exact, no lookups and no character reads.

Phases timed (`1-collect 1-sort 1-lookups 1-link 2-collect 2-sort
2-lookups 3-ms`); counters per phase; `work: lookups … lookups-skipped …`
shows how many lookups the stretches made unnecessary. Checked against
lrf-ms at every position.

```
g++ -std=c++20 -O2 ms2/sampled_ms_test.cpp -o sampled_ms_test && ./sampled_ms_test
```
