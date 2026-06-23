#!/usr/bin/env bash
set -euo pipefail

REF="$HOME/PhD/RLZData/ecoli_sample/reference.plain"
SA="$HOME/PhD/RLZData/ecoli_sample/reference.plain.sa"
FILES="$HOME/PhD/RLZData/ecoli_sample/input_filenames.txt"

OUTDIR="bucket_sweep_results"
mkdir -p "$OUTDIR"

echo "Running baseline once"
 /usr/bin/time -v ./build/compare_rlz_cache \
  --reference "$REF" \
  --suffix-array "$SA" \
  --filenames "$FILES" \
  --bucket-divisor 32 \
  --mode baseline \
  > "$OUTDIR/baseline_once.csv" \
  2> "$OUTDIR/baseline_once.log"

echo "Running cached bucket sweep"
for B in 8 16 32 64 128 256 512 1024 2048 4096; do
  echo "bucket=$B"

  /usr/bin/time -v ./build/compare_rlz_cache \
    --reference "$REF" \
    --suffix-array "$SA" \
    --filenames "$FILES" \
    --bucket-divisor "$B" \
    --mode cached \
    > "$OUTDIR/cached_bucket_${B}.csv" \
    2> "$OUTDIR/cached_bucket_${B}.log"
done

echo "Done"
