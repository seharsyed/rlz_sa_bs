#!/usr/bin/env bash
# The correctness stage: rlz_suite on every dataset and reference (config.sh).
# Every parser's phrases are decoded against the reference and must give the
# input back, and the left-to-right parsers (all of config.sh's VARIANTS) must
# have sa-binary-search's phrases: the same count and every length. Run it
# before the timing runs (run.sh); it takes a while (single-threaded, every
# parser on every file) but catches data or build problems.
#
#   ./verify.sh [DATASET ...]      (default: all of config.sh's DATASETS)
#
# Writes results/<dataset>/<ref>/verify.{log,csv}; stops at the first failure.
# A reference that already passed (verify.ok) is skipped (FORCE=1 checks it
# again). VERIFY_FILES=N checks only the first N inputs.

set -euo pipefail
source "$(dirname "$0")/config.sh"
source "$PIPELINE/common.sh"

build_programs rlz_suite
SUITE=$PROGRAM_PATH

for ds in ${@:-$DATASETS}; do
  list=$(input_list "$ds")
  for dir in $(dataset_refs "$ds"); do
    name=$(basename "$dir")
    ref=$dir/$name.txt
    out=$RESULTS/$ds/$name
    mkdir -p "$out"
    echo "=== $ds / $name ($(wc -l < "$list") inputs${VERIFY_FILES:+, checking the first $VERIFY_FILES})"
    if [ -f "$out/verify.ok" ] && [ "${FORCE:-0}" != 1 ]; then
      echo "    passed before (FORCE=1 to check again)"
      continue
    fi
    rm -f "$out/verify.ok"
    args=(--reference "$ref" --suffix-array "$ref.sa" --filenames "$list"
          --table "$dir/$name.suite_pt16" --results "$out/verify.csv" --quiet)
    [ -n "${VERIFY_FILES:-}" ] && args+=(--max-files "$VERIFY_FILES")
    [ -f "$ref.rev_four.bwt" ] && args+=(--powered-fwd-index "$ref.rev_four.bwt")
    if ${NUMA_CMD[@]+"${NUMA_CMD[@]}"} "$SUITE" "${args[@]}" > "$out/verify.log" 2>&1; then
      sed -n '/totals over/,/peak RSS/p' "$out/verify.log" | sed 's/^/  /'
      echo "    OK"
      date > "$out/verify.ok"
    else
      sed -n '/totals over/,$p' "$out/verify.log" | sed 's/^/  /'
      echo "FAILED: see $out/verify.log" >&2
      exit 1
    fi
  done
done
