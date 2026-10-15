#!/usr/bin/env bash
# The shape of every dataset's references' PT16 tables (config.sh), with
# ../pt16_stats.cpp: singletons vs ranges and the distribution of range sizes.
# Needs only the references and their suffix arrays; one reference at a time.
#
#   ./pt16_stats.sh [DATASET ...]      (default: all of config.sh's DATASETS)
#
# Writes into results/pt16/ (config.sh's RESULTS):
#   pt16_stats.csv     per dataset and reference: reference_bytes, entries
#                      (distinct 16-mers), singletons, ranges, positions
#                      (suffixes of length >= 16), range_positions,
#                      mean_range_size, expected_range_size (the range size of
#                      a random position; a singleton's counts as 0),
#                      max_range_size, ranges of size >= 1000, 10000, 100000,
#                      1000000
#   pt16_ranges.csv    the range sizes by power of two: dataset, reference,
#                      size_low, size_high, ranges
#   <dataset>_<ref>.log  the readable summary (largest ranges with their
#                      16-mers, entries per bucket)
# A reference whose row exists is skipped (FORCE=1 computes it again).

set -euo pipefail
source "$(dirname "$0")/config.sh"
source "$PIPELINE/common.sh"

build_programs pt16_stats
PT16_STATS=$PROGRAM_PATH

OUT=$RESULTS/pt16
TABLE=$OUT/pt16_stats.csv
BINS=$OUT/pt16_ranges.csv
mkdir -p "$OUT"
[ -s "$TABLE" ] || echo "dataset,reference,reference_bytes,entries,singletons,ranges,positions,range_positions,mean_range_size,expected_range_size,max_range_size,ranges_1000,ranges_10000,ranges_100000,ranges_1000000" > "$TABLE"
[ -s "$BINS" ] || echo "dataset,reference,size_low,size_high,ranges" > "$BINS"

for ds in ${@:-$DATASETS}; do
  for dir in $(dataset_refs "$ds"); do
    name=$(basename "$dir")
    ref=$dir/$name.txt
    if grep -q "^$ds,$name," "$TABLE" && [ "${FORCE:-0}" != 1 ]; then
      echo "=== $ds / $name: done (FORCE=1 to compute again)"
      continue
    fi
    echo "=== $ds / $name"
    values=$(${NUMA_CMD[@]+"${NUMA_CMD[@]}"} "$PT16_STATS" --reference "$ref" \
               --suffix-array "$ref.sa" 2> "$OUT/${ds}_${name}.log") ||
      { echo "FAILED: see $OUT/${ds}_${name}.log" >&2; exit 1; }
    get() { printf '%s\n' "$values" | grep "^$1=" | cut -d= -f2; }

    grep -v "^$ds,$name," "$TABLE" > "$TABLE.tmp" || true
    row="$ds,$name"
    for key in reference_bytes entries singletons ranges positions range_positions \
               mean_range_size expected_range_size max_range_size ranges_at_least_1000 \
               ranges_at_least_10000 ranges_at_least_100000 ranges_at_least_1000000; do
      row="$row,$(get $key)"
    done
    echo "$row" >> "$TABLE.tmp"
    mv "$TABLE.tmp" "$TABLE"

    grep -v "^$ds,$name," "$BINS" > "$BINS.tmp" || true
    printf '%s\n' "$values" | grep '^bin_' |
      sed -E "s/^bin_([0-9]+)_([0-9]+)=([0-9]+)$/$ds,$name,\1,\2,\3/" >> "$BINS.tmp"
    mv "$BINS.tmp" "$BINS"
    echo "    $(get entries) entries: $(get singletons) singletons, $(get ranges) ranges;" \
         "expected range size $(get expected_range_size)"
  done
done
echo "tables: $TABLE, $BINS"
