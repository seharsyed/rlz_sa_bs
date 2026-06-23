#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/PhD/RLZData/chr19"

REF_TXT="$BASE/reference.txt"
REF="$BASE/reference.plain"
SA="$BASE/reference.plain.sa"
FILES="$BASE/input_list.txt"

SA_BUILDER="$HOME/PhD/rlz/sa/build/sa"
OUTDIR="$BASE/bucket_sweep_chr19_results"

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "Checking chr19 experiment files"

echo "Reference source:"
cat "$REF_TXT"

echo "Input count:"
wc -l "$FILES"

if grep -Fxq -f "$REF_TXT" "$FILES"; then
    echo "ERROR: reference is present in input_list.txt"
    exit 1
fi

echo "OK: reference is not in input_list.txt"

echo "Checking reference.plain sentinel:"
last_char=$(tail -c 1 "$REF")
if [[ "$last_char" != '$' ]]; then
    echo "ERROR: reference.plain does not end with sentinel $"
    exit 1
fi

echo "OK: reference.plain ends with sentinel"

echo "Checking input files exist"
while IFS= read -r f; do
    [[ -f "$f" ]] || { echo "ERROR: missing input file: $f"; exit 1; }
done < "$FILES"

echo "Removing old suffix array"
rm -f "$SA"

echo "Building suffix array"
"$SA_BUILDER" "$REF" "$SA"

echo "Running chr19 baseline once"

/usr/bin/time -v ./build/compare_rlz_cache \
  --reference "$REF" \
  --suffix-array "$SA" \
  --filenames "$FILES" \
  --bucket-divisor 32 \
  --mode baseline \
  > "$OUTDIR/baseline_once.csv" \
  2> "$OUTDIR/baseline_once.log"

echo "Running chr19 cached bucket sweep"

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

echo "Done."
echo "Results in: $OUTDIR"
