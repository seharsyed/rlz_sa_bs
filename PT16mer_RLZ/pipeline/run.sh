#!/usr/bin/env bash
# The timing runs (config.sh): for every dataset, reference and variant, one
# rlz_parallel run over the dataset's whole input collection -- one
# reference, one variant at a time:
#
#   the first reference: every thread count of THREADS_SWEEP (the structures
#                        are built once per run, then one parallel pass per
#                        thread count); its 1-thread pass is also plot 1's
#   the others:          1 thread
#
# Only the parse is timed (inputs are preloaded; tables and indexes are built
# or loaded before, and reported apart). The PT16 tables are rebuilt by every
# run, into the reference's folder.
#
#   ./run.sh [DATASET ...]      (default: all of config.sh's DATASETS)
#
# Writes, per dataset, reference and variant, into results/<dataset>/<ref>/:
#   <variant>.summary.csv   one row per thread count (wall time, work, ...)
#   <variant>.per_file.csv  every file's parse time and phrases, per thread count
#   <variant>.log           the run's output
# A variant whose summary exists is skipped (FORCE=1 runs it again).

set -euo pipefail
source "$(dirname "$0")/config.sh"
source "$PIPELINE/common.sh"

build_programs rlz_parallel
PARALLEL=$PROGRAM_PATH

for ds in ${@:-$DATASETS}; do
  list=$(input_list "$ds")
  first=1
  for dir in $(dataset_refs "$ds"); do
    name=$(basename "$dir")
    ref=$dir/$name.txt
    out=$RESULTS/$ds/$name
    mkdir -p "$out"
    threads=1
    [ "$first" -eq 1 ] && threads=$THREADS_SWEEP
    first=0
    echo "=== $ds / $name: threads $threads, $(wc -l < "$list") inputs"

    for variant in $VARIANTS; do
      summary=$out/$variant.summary.csv
      if [ -s "$summary" ] && [ "${FORCE:-0}" != 1 ]; then
        echo "    $variant: done (FORCE=1 to run again)"
        continue
      fi
      args=(--reference "$ref" --suffix-array "$ref.sa" --filenames "$list"
            --parsers "$variant" --threads-list "$threads"
            --table "$dir/$name.pt16" --results "$summary.tmp"
            --per-file "$out/$variant.per_file.csv" --quiet)
      [ -f "$ref.rev_four.bwt" ] && args+=(--powered-fwd-index "$ref.rev_four.bwt")

      start=$(date +%s)
      echo "# $(date), $(uptime), prefix: ${NUMA_CMD[*]:-none}" > "$out/$variant.log"
      ${NUMA_CMD[@]+"${NUMA_CMD[@]}"} "$PARALLEL" "${args[@]}" >> "$out/$variant.log" 2>&1 ||
        { echo "FAILED: $variant, see $out/$variant.log" >&2; exit 1; }
      # A variant this build does not have (varki without sdsl, powered off
      # x86-64, a missing index) gives no rows: stop rather than skip it.
      if [ "$(grep -c "^$variant," "$summary.tmp" || true)" -eq 0 ]; then
        echo "FAILED: $variant did not run, see $out/$variant.log" >&2
        rm -f "$summary.tmp"
        exit 1
      fi
      mv "$summary.tmp" "$summary"
      wall=$(awk -F, -v v="$variant" '$1 == v && $2 == 1 {print $8}' "$summary")
      echo "    $variant: $(( $(date +%s) - start )) s (1-thread wall ${wall} ms)"
    done
  done
done
echo "next: ./plot.py"
