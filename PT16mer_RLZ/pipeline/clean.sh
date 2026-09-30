#!/usr/bin/env bash
# Cleans the raw inputs and references of the datasets (config.sh):
#
#   raw/<dataset>/inputs/<file>  ->  data/<dataset>/inputs/<name>.txt
#   raw/<dataset>/refs/<file>    ->  data/<dataset>/refs/<name>/<name>.txt
#
# Cleaning (clean_reference): FASTA headers and whitespace are removed,
# lower-case (soft-masked) bases are upper-cased, and every other character
# (N, IUPAC codes, ...) becomes A; the result is plain A/C/G/T, one line, no
# newline. .gz files are decompressed on the fly. A cleaned file newer than
# its raw file is kept. Statistics per file go to data/<dataset>/cleaning/.
#
#   ./clean.sh [DATASET ...]      (default: all of config.sh's DATASETS)

set -euo pipefail
source "$(dirname "$0")/config.sh"

CLEAN=$TOOLS/bin/clean_reference
[ -x "$CLEAN" ] || { echo "missing $CLEAN: run setup.sh first" >&2; exit 1; }

# clean_one RAW OUT STATS: cleans RAW into OUT unless OUT is up to date.
clean_one() {
  local raw=$1 out=$2 stats=$3
  if [ -s "$out" ] && [ "$out" -nt "$raw" ]; then
    return 1
  fi
  mkdir -p "$(dirname "$out")" "$(dirname "$stats")"
  if [[ "$raw" == *.gz ]]; then
    local tmp
    tmp=$(mktemp "${out}.XXXX.raw")
    gzip -dc "$raw" > "$tmp"
    "$CLEAN" "$tmp" "$out" "$stats" > /dev/null
    rm -f "$tmp"
  else
    "$CLEAN" "$raw" "$out" "$stats" > /dev/null
  fi
}

# The number of replaced characters recorded in a statistics CSV.
replaced() { grep '^replacement_count,' "$1" | cut -d, -f2; }

for ds in ${@:-$DATASETS}; do
  echo "=== $ds"
  [ -d "$RAW/$ds" ] || { echo "    no $RAW/$ds: skipped"; continue; }

  # ---------- inputs ----------
  done_count=0 kept=0 total=0
  shopt -s nullglob
  for raw in "$RAW/$ds/inputs"/*; do
    [ -f "$raw" ] || continue
    total=$((total + 1))
    name=$(file_stem "$raw")
    out=$DATA/$ds/inputs/$name.txt
    stats=$DATA/$ds/cleaning/inputs/$name.csv
    if clean_one "$raw" "$out" "$stats"; then
      done_count=$((done_count + 1))
    else
      kept=$((kept + 1))
    fi
  done
  echo "    inputs: $total raw files, $done_count cleaned, $kept already clean"

  # ---------- references ----------
  for raw in "$RAW/$ds/refs"/*; do
    [ -f "$raw" ] || continue
    name=$(file_stem "$raw")
    out=$DATA/$ds/refs/$name/$name.txt
    stats=$DATA/$ds/cleaning/refs/$name.csv
    if clean_one "$raw" "$out" "$stats"; then
      echo "    reference $name: cleaned ($(stat -c %s "$out" 2>/dev/null || wc -c < "$out") bases, $(replaced "$stats") replaced by A)"
    else
      echo "    reference $name: already clean"
    fi
  done
  shopt -u nullglob

  # Files with many replacements (N runs, ...), for a look.
  if [ -d "$DATA/$ds/cleaning" ]; then
    many=$(for s in "$DATA/$ds/cleaning"/*/*.csv; do
             r=$(replaced "$s")
             if [ "${r:-0}" -gt 0 ]; then echo "$r $(basename "$s" .csv)"; fi
           done | sort -rn | head -3 | tr '\n' ';')
    if [ -n "$many" ]; then echo "    most replacements (count name): $many"; fi
  fi
done
