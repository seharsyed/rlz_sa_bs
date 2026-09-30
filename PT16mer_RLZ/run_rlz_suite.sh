#!/usr/bin/env bash
# Runs the RLZ parsing suite (rlz_suite: single-threaded, with checks) on a
# dataset.
#
#   ./run_rlz_suite.sh DATASET [LIST] [rlz_suite options...]
#
# DATASET  a name from rlz_datasets.sh (e.g. ecoli)
# LIST     the input list; optional when the dataset has a default one
# options  passed on, e.g. --quiet, --max-files 10, --results FILE
#
# rlz_suite is rebuilt when it is missing or older than its sources (with
# powered on x86-64, without elsewhere; see build_if_needed in rlz_datasets.sh).
#
# NUMA=N pins the run to NUMA node N (its CPU and memory; see numa_setup in
# rlz_datasets.sh), e.g.  NUMA=0 ./run_rlz_suite.sh ecoli --quiet

set -euo pipefail

cd "$(dirname "$0")"
source ./rlz_datasets.sh

if [ $# -lt 1 ]; then
  sed -n '2,/^$/p' "$0"
  echo "datasets: $(list_datasets)"
  exit 1
fi

DATASET=$1
shift
select_dataset "$DATASET"

if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then
  LIST=$1
  shift
fi
[ -n "$LIST" ] || { echo "no input list for $DATASET: give one" >&2; exit 1; }

for f in "$REF" "$SA" "$LIST"; do
  [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

build_if_needed rlz_suite

POWERED_ARGS=()
if [ -f "$POWERED" ]; then
  POWERED_ARGS=(--powered-index "$POWERED")
else
  echo "(no powered index at $POWERED: powered is skipped)" >&2
fi
if [ -f "$POWERED_FWD" ]; then
  POWERED_ARGS+=(--powered-fwd-index "$POWERED_FWD")
else
  echo "(no reversed powered index at $POWERED_FWD: the fwd variants are skipped)" >&2
fi

numa_setup

exec ${NUMA_CMD[@]+"${NUMA_CMD[@]}"} "./$BINARY" --reference "$REF" --suffix-array "$SA" --filenames "$LIST" \
  ${POWERED_ARGS[@]+"${POWERED_ARGS[@]}"} "$@"
