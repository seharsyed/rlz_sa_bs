#!/usr/bin/env bash
# Collects run.sh's results into two flat CSV files, for plotting elsewhere
# (config.sh's layout):
#
#   results/per_file.csv   one row per dataset, reference, variant, thread
#                          count and input file:
#                          dataset,reference,parser,threads,file,input_bytes,
#                          parse_ms,phrases
#                          (plot 1: the rows with threads = 1)
#   results/summary.csv    one row per dataset, reference, variant and thread
#                          count: dataset,reference,sweep_reference, then
#                          rlz_parallel's summary columns (wall_ms: the whole
#                          collection's parse; build_ms, load_ms, own/needed
#                          bytes, ...)
#                          (plot 2: sweep_reference = 1)
#
#   ./collect.sh [DATASET ...]    (default: all of config.sh's DATASETS)

set -euo pipefail
source "$(dirname "$0")/config.sh"

PER_FILE=$RESULTS/per_file.csv
SUMMARY=$RESULTS/summary.csv
mkdir -p "$RESULTS"
: > "$PER_FILE.tmp"
: > "$SUMMARY.tmp"
rows_per_file=0 rows_summary=0

for ds in ${@:-$DATASETS}; do
  [ -d "$RESULTS/$ds" ] || continue
  first=1
  for dir in $(dataset_refs "$ds"); do
    ref=$(basename "$dir")
    out=$RESULTS/$ds/$ref
    sweep=$first
    first=0
    [ -d "$out" ] || continue
    for variant in $VARIANTS; do
      if [ -f "$out/$variant.per_file.csv" ]; then
        # the file column: the input's name, not its full path
        awk -F, -v OFS=, -v d="$ds" -v r="$ref" 'NR > 1 {
              n = split($3, p, "/"); $3 = p[n]
              print d, r, $0 }' "$out/$variant.per_file.csv" >> "$PER_FILE.tmp"
      fi
      if [ -f "$out/$variant.summary.csv" ]; then
        awk -F, -v OFS=, -v d="$ds" -v r="$ref" -v s="$sweep" 'NR > 1 {
              print d, r, s, $0 }' "$out/$variant.summary.csv" >> "$SUMMARY.tmp"
        header=$(head -1 "$out/$variant.summary.csv")
      fi
    done
  done
done

{ echo "dataset,reference,parser,threads,file,input_bytes,parse_ms,phrases"
  cat "$PER_FILE.tmp"; } > "$PER_FILE"
{ echo "dataset,reference,sweep_reference,${header:-parser,threads}"
  cat "$SUMMARY.tmp"; } > "$SUMMARY"
rm -f "$PER_FILE.tmp" "$SUMMARY.tmp"
echo "$PER_FILE: $(($(wc -l < "$PER_FILE") - 1)) rows"
echo "$SUMMARY: $(($(wc -l < "$SUMMARY") - 1)) rows"
