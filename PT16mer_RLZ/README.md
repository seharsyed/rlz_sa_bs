# PT16 RLZ

This repository contains three RLZ experiment modes:

1. Baseline RLZ
2. PT16 RLZ
3. Cached RLZ

## Build

From the repository root:

```bash
mkdir -p build
cd build
cmake ..
make -j
cd ..
```

The main executables are:

```text
build/pt16_rlz
build/compare_rlz_cache
```

## Required input files

The experiments require:

- a reference sequence
- the suffix array built for that reference
- a text file containing one input-file path per line

Example:

```text
data/
├── reference.plain
├── reference.plain.sa
├── input_list.txt
└── inputs/
    ├── input1.fna
    ├── input2.fna
    └── ...
```

Create an output directory if needed:

```bash
mkdir -p results
```

## 1. Baseline RLZ

Run the baseline using `compare_rlz_cache` in baseline mode:

```bash
./build/compare_rlz_cache \
  --reference data/reference.plain \
  --suffix-array data/reference.plain.sa \
  --filenames data/input_list.txt \
  --bucket-divisor 8 \
  --mode baseline \
  > results/baseline.csv
```

## 2. PT16 RLZ

Run:

```bash
./build/pt16_rlz \
  --reference data/reference.plain \
  --suffix-array data/reference.plain.sa \
  --filenames data/input_list.txt \
  --results results/pt16.csv
```

`pt16_rlz` runs the baseline first and then PT16 so that timing and factor output can be compared directly.

If the PT16 table does not exist, it is created automatically beside the reference as:

```text
reference.plain.pt16_hl.bin
```

If the PT16 table format or construction code changes, remove the old table before rerunning:

```bash
rm data/reference.plain.pt16_hl.bin
```

The next PT16 run will rebuild it.

## 3. Cached RLZ

Run the cached parser with the required bucket divisor. For example, with divisor 8:

```bash
./build/compare_rlz_cache \
  --reference data/reference.plain \
  --suffix-array data/reference.plain.sa \
  --filenames data/input_list.txt \
  --bucket-divisor 8 \
  --mode cached \
  > results/cache_b8.csv
```

To run baseline and cached RLZ together in one experiment:

```bash
./build/compare_rlz_cache \
  --reference data/reference.plain \
  --suffix-array data/reference.plain.sa \
  --filenames data/input_list.txt \
  --bucket-divisor 8 \
  --mode both \
  > results/cache_b8.csv
```

## Typical comparison workflow

For a direct comparison of all three methods on the same dataset:

```bash
./build/pt16_rlz \
  --reference data/reference.plain \
  --suffix-array data/reference.plain.sa \
  --filenames data/input_list.txt \
  --results results/pt16.csv

./build/compare_rlz_cache \
  --reference data/reference.plain \
  --suffix-array data/reference.plain.sa \
  --filenames data/input_list.txt \
  --bucket-divisor 8 \
  --mode both \
  > results/cache_b8.csv
```

The first command records baseline and PT16 results.  
The second command records baseline and cached-RLZ results.

Use the same reference, suffix array, and input list for all methods when comparing performance.
