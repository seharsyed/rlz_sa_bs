#!/usr/bin/env bash
set -euo pipefail

BASE="/path/to/folder_containing_reference_and_input_list"
REF_PLAIN="$BASE/reference.plain"
SA="$BASE/reference.plain.sa"
INPUT_LIST="$BASE/input_list.txt"
OUTDIR="$BASE/rlz_sa_bs_results"

SA_BUILDER="/path/to/rlz_sa_bs/sa/build/sa"
RLZ_BIN="/path/to/rlz_sa_bs/parser/build/compare_rlz_cache"

mkdir -p "$OUTDIR"

echo "Building suffix array"
"$SA_BUILDER" "$REF_PLAIN" "$SA"

echo "Running baseline"

/usr/bin/time -v "$RLZ_BIN" \
  --reference "$REF_PLAIN" \
  --suffix-array "$SA" \
  --filenames "$INPUT_LIST" \
  --bucket-divisor 32 \
  --mode baseline \
  > "$OUTDIR/baseline_once.csv" \
  2> "$OUTDIR/baseline_once.log"

echo "Running cached bucket sweep"

for B in 8 16 32 64 128 256 512 1024 2048 4096; do
  /usr/bin/time -v "$RLZ_BIN" \
    --reference "$REF_PLAIN" \
    --suffix-array "$SA" \
    --filenames "$INPUT_LIST" \
    --bucket-divisor "$B" \
    --mode cached \
    > "$OUTDIR/cached_bucket_${B}.csv" \
    2> "$OUTDIR/cached_bucket_${B}.log"
done

echo "Done. Results in $OUTDIR"