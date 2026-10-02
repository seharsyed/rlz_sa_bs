#!/usr/bin/env bash
# Matching-statistics and LCP statistics of every dataset and reference
# (config.sh), with lrf-ms (../ms_stats.cpp). Timing does not matter: every
# available CPU (or CPUS) is used.
#
#   ./ms_stats.sh [DATASET ...]      (default: all of config.sh's DATASETS)
#
# Writes into results/ms/ (config.sh's RESULTS):
#   <dataset>_<ref>.per_file.csv  per input: file, input_bytes, ms_mean,
#                                 ms_median, ms_max
#   ms_stats.csv                  per dataset and reference: reference_bytes,
#                                 lcp_mean (the reference's average LCP over
#                                 its n - 1 adjacent suffix pairs), lcp_max,
#                                 files, input_bytes, ms_mean (average MS
#                                 length over all input positions),
#                                 ms_mean_of_files (mean of the per-file means)
# A reference whose row exists is skipped (FORCE=1 computes it again).

set -euo pipefail
source "$(dirname "$0")/config.sh"
source "$PIPELINE/common.sh"

build_programs ms_stats
MS_STATS=$PROGRAM_PATH

OUT=$RESULTS/ms
TABLE=$OUT/ms_stats.csv
mkdir -p "$OUT"
[ -s "$TABLE" ] || echo "dataset,reference,reference_bytes,lcp_mean,lcp_max,files,input_bytes,ms_mean,ms_mean_of_files" > "$TABLE"

for ds in ${@:-$DATASETS}; do
  list=$(input_list "$ds")
  for dir in $(dataset_refs "$ds"); do
    name=$(basename "$dir")
    ref=$dir/$name.txt
    if grep -q "^$ds,$name," "$TABLE" && [ "${FORCE:-0}" != 1 ]; then
      echo "=== $ds / $name: done (FORCE=1 to compute again)"
      continue
    fi
    echo "=== $ds / $name"
    values=$(${NUMA_CMD[@]+"${NUMA_CMD[@]}"} "$MS_STATS" --reference "$ref" \
               --suffix-array "$ref.sa" --filenames "$list" \
               --per-file "$OUT/${ds}_${name}.per_file.csv" 2> "$OUT/${ds}_${name}.log") ||
      { echo "FAILED: see $OUT/${ds}_${name}.log" >&2; exit 1; }
    get() { printf '%s\n' "$values" | grep "^$1=" | cut -d= -f2; }
    grep -v "^$ds,$name," "$TABLE" > "$TABLE.tmp" || true
    echo "$ds,$name,$(get reference_bytes),$(get lcp_mean),$(get lcp_max),$(get files),$(get input_bytes),$(get ms_mean),$(get ms_mean_of_files)" >> "$TABLE.tmp"
    mv "$TABLE.tmp" "$TABLE"
    echo "    average LCP $(get lcp_mean), average MS $(get ms_mean)"
  done
done
echo "table: $TABLE"
