#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/PhD/RLZData/bacteria"
INPUT_DIR="$BASE/input"

REF_TXT="$BASE/reference.txt"
REF_PLAIN="$BASE/reference.plain"
SA="$BASE/reference.plain.sa"
INPUT_LIST="$BASE/input_list.txt"

SA_BUILDER="$HOME/PhD/rlz/sa/build/sa"
OUTDIR="$BASE/bucket_sweep_bacteria_results"

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "Initializing bacteria experiment files"

mapfile -t files < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fna")

if [[ ${#files[@]} -lt 101 ]]; then
    echo "ERROR: need at least 101 .fna files in $INPUT_DIR"
    echo "Found: ${#files[@]}"
    exit 1
fi

REF_SOURCE="${files[0]}"

printf "%s\n" "$REF_SOURCE" > "$REF_TXT"

cp "$REF_SOURCE" "$REF_PLAIN"
printf '$' >> "$REF_PLAIN"

printf "%s\n" "${files[@]:1:100}" > "$INPUT_LIST"

echo "Reference:"
cat "$REF_TXT"

echo "Input count:"
wc -l "$INPUT_LIST"

if grep -Fxq "$REF_SOURCE" "$INPUT_LIST"; then
    echo "ERROR: reference is present in input_list.txt"
    exit 1
fi

echo "OK: reference excluded from input_list.txt"

rm -f "$SA"

echo "Building suffix array"
"$SA_BUILDER" "$REF_PLAIN" "$SA"

echo "Running bacteria baseline once"

/usr/bin/time -v ./build/compare_rlz_cache \
  --reference "$REF_PLAIN" \
  --suffix-array "$SA" \
  --filenames "$INPUT_LIST" \
  --bucket-divisor 32 \
  --mode baseline \
  > "$OUTDIR/baseline_once.csv" \
  2> "$OUTDIR/baseline_once.log"

echo "Running bacteria cached bucket sweep"

for B in 8 16 32 64 128 256 512 1024 2048 4096; do
  echo "bucket=$B"

  /usr/bin/time -v ./build/compare_rlz_cache \
    --reference "$REF_PLAIN" \
    --suffix-array "$SA" \
    --filenames "$INPUT_LIST" \
    --bucket-divisor "$B" \
    --mode cached \
    > "$OUTDIR/cached_bucket_${B}.csv" \
    2> "$OUTDIR/cached_bucket_${B}.log"
done

echo "Done."
echo "Results in: $OUTDIR"