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
