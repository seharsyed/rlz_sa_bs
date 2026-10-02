#!/usr/bin/env bash
# Runs rlz_parallel (RLZ parsing, parallel over input files) on a dataset.
#
#   scripts/run_rlz_parallel.sh DATASET [LIST] [rlz_parallel options...]
#
# DATASET  a name from rlz_datasets.sh (e.g. ecoli)
# LIST     the input list; optional when the dataset has a default one
# options  passed on, e.g. --threads 16, --parsers sassy,powered-pt16-escape,
#          --max-files 10, --no-preload, --results FILE
#
# rlz_parallel is rebuilt when it is missing or older than its sources (with
# powered on x86-64, without elsewhere; see build_if_needed in rlz_datasets.sh).
# Unless --results is given, the CSV goes to
# results/parallel_<dataset>_t<threads>_<numa>_<date>.csv.
#
# NUMA=N pins the run to NUMA node N (threads and memory; see numa_setup in
# rlz_datasets.sh), e.g.  NUMA=0 scripts/run_rlz_parallel.sh ecoli --threads 16

set -euo pipefail

cd "$(dirname "$0")/.."
source scripts/rlz_datasets.sh

if [ $# -lt 1 ]; then
  sed -n '2,/^$/p' "$0"
  echo "datasets: $(list_datasets)"
  exit 1
fi

DATASET=$1
shift
select_dataset "$DATASET"

# An optional LIST (any first argument that is not an option).
if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then
  LIST=$1
  shift
fi
[ -n "$LIST" ] || { echo "no input list for $DATASET: give one" >&2; exit 1; }

for f in "$REF" "$SA" "$LIST"; do
  [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

# ---------- Build if needed ----------

build_if_needed rlz_parallel

# ---------- Run ----------

numa_setup

ARGS=(--reference "$REF" --suffix-array "$SA" --filenames "$LIST")
if [ -f "$POWERED" ]; then
  ARGS+=(--powered-index "$POWERED")
else
  echo "(no powered index at $POWERED: powered is skipped)" >&2
fi
if [ -f "$POWERED_FWD" ]; then
  ARGS+=(--powered-fwd-index "$POWERED_FWD")
else
  echo "(no reversed powered index at $POWERED_FWD: the fwd variants are skipped)" >&2
fi

# A default results file, named after the dataset and the thread count.
if ! printf '%s\n' "$@" | grep -qx -- "--results"; then
  THREADS=$(nproc 2>/dev/null || echo all)
  prev=""
  for a in "$@"; do
    [ "$prev" = "--threads" ] && THREADS=$a
    prev=$a
  done
  mkdir -p results
  ARGS+=(--results "results/parallel_${DATASET}_t${THREADS}_${NUMA_LABEL}_$(date +%Y%m%d_%H%M).csv")
fi

echo "dataset $DATASET: reference $REF" >&2
exec ${NUMA_CMD[@]+"${NUMA_CMD[@]}"} "./$BINARY" "${ARGS[@]}" "$@"
