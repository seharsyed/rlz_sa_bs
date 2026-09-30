# Powered RLZ in our benchmark

`powered_bench.cpp` runs the parser of `../../powered_rlz` (RLZ by powered
backward search in an FM-index) under our benchmark conventions:

- the powered index is given as an argument and loaded once (timed apart);
- our input list is used and each input is loaded with our loader,
  untimed; **only the parse is timed**, single-threaded, one file at a time;
- every parse is decoded against the reference and must reproduce the input
  exactly (powered RLZ parses right to left, so its phrases differ from a
  left-to-right greedy parse; the decode check is what proves it correct).

x86-64 only (the powered_rlz headers use x86 intrinsics), built with GCC.

## 1. Build powered_rlz's tools and the index (once per reference)

The index needs the reference's extended BWT and its GCA (the cyclic
rotation order), made with [PFP-eBWT](https://github.com/davidecenzato/PFP-eBWT):

```
pfpebwt --GCA REF                      # writes REF.ebwt and REF.gca (and more)

cd powered_rlz
make
./make_bwt --rle -i REF.ebwt -sa REF.gca -o REF.bwt        # run again if it reports 0 dense blocks
./transform -b -i REF.bwt -sa REF.gca -o REF_four.txt      # writes REF_four.bwt: the powered index
```

The reference and the inputs must be plain ACGT: no `N`, no newline
(powered RLZ maps any other byte to one of A/C/G/T).

## 2. Build and run the benchmark

From `PT16mer_RLZ/`. The block sizes must match the ones powered_rlz was
built with (its Makefile's defaults):

```
g++ -std=c++2a -O3 -march=native -DNDEBUG \
    -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 \
    powered/powered_bench.cpp -o powered_bench

./powered_bench --index REF_four.bwt --reference REF --filenames LIST \
                [--results CSV] [--max-files N] [--code-size 4]
```

`--code-size 1` is for an index over the original alphabet (rlz_parser's
`--original`); the powered index uses 4.

Output: per file the parse time, phrase count, average phrase length and
the decode check; totals at the end; `key=value` lines on stdout
(`index_load_ms`, `powered_parse_ms`, `powered_phrases`, `decode_ok`).

## 3. PT16 table for powered (sassy layout)

`pt16_powered.hpp` mirrors the sassy PT16 table (`variants/pt16_sassy.hpp`,
`variants/pt16_build_sassy.hpp`: same H directory with empty-bucket fallback,
same 64-bit L entries via `pt16_sassy_format.hpp`, same bucket search, miss
handling and precomputed empty-bucket answers) for powered's **backward**
search. powered's phrase is the longest *suffix* of the remaining input, so
the table is keyed by the 16-mer read **from its right end**: then a miss
gives the longest common suffix, just as a forward PT16 miss gives the
longest common prefix.

A lookup of the 16-mer ending at a phrase start:

- **singleton hit**: the entry holds the 16-mer's reference position; the
  match is unique, so the parse extends it by comparing characters;
- **range hit**: the row interval `[a, b)` in the powered index (from the
  range's slice: first row and one position) -- the state after powered's
  first metacharacter and three rank steps -- and the backward search goes on;
- **miss**: the phrase is shorter than 16 and the table gives it outright
  (length = common suffix with the neighbouring key, position from that
  entry, or precomputed for an empty bucket): no backward search at all.

A singleton hit, and a miss next to a singleton, read one L entry and nothing
else (no row array, no `gca_`).

It is built from the reference and the powered index (`REF_four.bwt`), not
from a suffix array: the index's rows are the sorted rotations of the cyclic
reference and `gca_[row]` is each rotation's start. Every row has a full
(cyclic) 16-mer, so every occurrence of any string ends some table 16-mer
(no short-suffix records are needed).

```
g++ -std=c++2a -O3 -march=native -DNDEBUG \
    -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 \
    powered/pt16_powered_build.cpp -o pt16_powered_build
g++ -std=c++2a -O3 -march=native -DNDEBUG \
    -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 \
    powered/pt16_powered_check.cpp -o pt16_powered_check

./pt16_powered_build --reference REF --index REF_four.bwt   # writes REF_four.pt16
./pt16_powered_check --reference REF --index REF_four.bwt \
                     --table REF_four.pt16 --filenames LIST [--max-files N]
```

The check verifies every entry (its position holds its 16-mer, its lookup
returns it, its size equals powered's own `count()`) and random 16-mers (on a
miss: the returned suffix occurs at the returned position and one character
more does not occur), and on the inputs that at every phrase start the lookup
hits exactly when powered's phrase is >= 16 long and otherwise returns
powered's phrase length.

### Variants in the benchmark

`powered_pt16_parse.hpp` is powered's `parse_tuples` (code size 4) with an
optional table and an optional reference:

- `powered-escape` (`--escape`): powered, but once the interval is a single
  row the phrase is extended leftwards by comparing characters with the
  cyclic reference (as the SA and PT16 parsers do) instead of rank steps;
- `powered-pt16-escape` (`--pt16-table P` or `--pt16`): the table lookup at
  every phrase start, and the escape (a singleton entry has no row to rank
  from, so the table always comes with the reference).

The phrase lengths are identical to powered's; a phrase may point at another
occurrence, so each variant is also checked by decoding. powered's C array is
private, so it is rebuilt from `count()` at load; `powered_rlz` is unchanged.

```
./powered_bench --index REF_four.bwt --reference REF --filenames LIST \
                --pt16-table REF_four.pt16 --escape --quiet
```
