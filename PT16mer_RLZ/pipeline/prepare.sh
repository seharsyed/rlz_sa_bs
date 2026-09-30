#!/usr/bin/env bash
# Computes everything each reference needs, into its folder
# data/<dataset>/refs/<ref>/ (config.sh), keeping whatever already exists:
#
#   <ref>.txt.sa           its suffix array (4-byte entries, no sentinel), for
#                          sa-binary-search, lrf-ms and the PT16 tables; an
#                          existing one is kept if its size fits
#   <ref>.txt.rev, ...     the reversed reference, its PFP-eBWT files and its
#                          powered index <ref>.txt.rev_four.bwt (+ _data), for
#                          the powered fwd variants
#                          (../powered/prepare_reference.sh --reversed)
#
# The PT16 tables are not made here: every run rebuilds them.
#
#   ./prepare.sh [DATASET ...]    (default: all of config.sh's DATASETS)

set -euo pipefail
source "$(dirname "$0")/config.sh"

# The suffix array tool: the parallel one (psa, OpenMP) if built, else sa.
SA_TOOL=$TOOLS/bin/psa
[ -x "$SA_TOOL" ] || SA_TOOL=$TOOLS/bin/sa
[ -x "$SA_TOOL" ] || { echo "missing $TOOLS/bin/psa (or sa): run setup.sh first" >&2; exit 1; }
[ -x "$PFP/pfpebwt" ] || { echo "missing $PFP/pfpebwt: run setup.sh first" >&2; exit 1; }

file_size() { stat -c %s "$1" 2>/dev/null || wc -c < "$1"; }

for ds in ${@:-$DATASETS}; do
  echo "=== $ds"
  refs=$(dataset_refs "$ds")
  [ -n "$refs" ] || { echo "    no references in $DATA/$ds/refs (run clean.sh)"; continue; }

  for dir in $refs; do
    name=$(basename "$dir")
    ref=$dir/$name.txt
    [ -s "$ref" ] || { echo "    $name: no $ref (run clean.sh)" >&2; exit 1; }
    n=$(file_size "$ref")
    echo "--- $name ($n bases)"

    # ---------- suffix array ----------
    if [ -s "$ref.sa" ] && [ "$(file_size "$ref.sa")" -eq $((4 * n)) ]; then
      echo "    suffix array: exists"
    else
      [ -e "$ref.sa" ] && echo "    suffix array: wrong size, rebuilt"
      start=$(date +%s)
      "$SA_TOOL" "$ref" "$ref.sa.tmp"
      [ "$(file_size "$ref.sa.tmp")" -eq $((4 * n)) ] ||
        { echo "ERROR: $ref.sa.tmp has the wrong size" >&2; exit 1; }
      mv "$ref.sa.tmp" "$ref.sa"
      echo "    suffix array: built in $(( $(date +%s) - start )) s ($(basename "$SA_TOOL"))"
    fi

    # ---------- powered index of the reversed reference ----------
    if [ -s "$ref.rev_four.bwt" ] && [ -s "$ref.rev_four_data.bwt" ]; then
      echo "    reversed powered index: exists"
    else
      echo "    reversed powered index: building"
      PFP=$PFP REPO=$REPO "$PT16/powered/prepare_reference.sh" --reversed "$ref" "$dir" \
        2>&1 | sed 's/^/      /'
      [ -s "$ref.rev_four.bwt" ] ||
        { echo "ERROR: no $ref.rev_four.bwt after prepare_reference.sh" >&2; exit 1; }
    fi
  done
done
