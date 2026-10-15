#!/usr/bin/env bash
# Single-thread runtimes of exact matching statistics with the PT16 tables,
# against lrf-ms (the classical suffix-array + LCP + RMQ algorithm), per
# dataset (config.sh), on its first reference (the plots' reference 1;
# ALL_REFS=1 for all of them), over its inputs. No chain anywhere.
#
#   ./ms_times.sh [DATASET ...]      (default: all of config.sh's DATASETS)
#
# Two programs (../ms2/README.md), each run once per reference:
#   lazy_main         lrf-ms, lazy-v2, lazy-sassy: backward left-extension, one
#                     PT16 lookup only where the extension breaks
#   ms2_main          lrf-ms again, sampled-ms-v2, sampled-ms-sassy: lookups
#     --no-chain      at every 16th position, linked into stretches, then
#                     lookups only in the holes (its other rows are lookups
#                     only, not matching statistics)
# Every variant's lengths are checked against lrf-ms at every position (not
# timed). Times exclude file IO; the tables (v2, sassy) and the lrf-ms index
# are built once per reference, from the given suffix array, and timed apart.
#
# Writes into results/ms_times/ (config.sh's RESULTS):
#   ms_times.csv                per dataset and reference: files, input_bytes,
#                               lrf_build_ms, v2_build_ms, sassy_build_ms,
#                               lrf_ms, lazy_v2_ms, lazy_sassy_ms, lazy_ok,
#                               sampled_v2_ms, sampled_sassy_ms, sampled_ok
#                               (MS times summed over the inputs)
#   <dataset>_<ref>.{lazy,ms2}.{out,log,csv}  the programs' totals, logs and
#                               per-file rows
# A reference whose row exists is skipped (FORCE=1 times it again). Each run
# is stopped after MS_TIMEOUT (default 6h).

set -euo pipefail
source "$(dirname "$0")/config.sh"
source "$PIPELINE/common.sh"

# build SOURCE NAME: ROOT/bin/NAME.<host> from PT16mer_RLZ/SOURCE, rebuilt when
# the contents of SOURCE or any header under PT16mer_RLZ change; sets BUILT.
build() {
  local source=$1 name=$2 sum
  BUILT=$ROOT/bin/$name.$(hostname -s)
  mkdir -p "$ROOT/bin"
  sum=$(cd "$PT16" && { find . -name '*.hpp' -o -name '*.h' | sort; echo "$source"; } |
          xargs cat | cksum)
  if [ ! -x "$BUILT" ] || [ "$(cat "$BUILT.sources" 2>/dev/null)" != "$sum" ]; then
    echo "(building $BUILT)" >&2
    local march=(); [ "$(uname -m)" = "x86_64" ] && march=(-march=native)
    g++ -std=c++20 -O3 ${march[@]+"${march[@]}"} -DNDEBUG "$PT16/$source" -o "$BUILT.tmp"
    mv "$BUILT.tmp" "$BUILT"
    echo "$sum" > "$BUILT.sources"
  fi
}
build ms2/lazy_main.cpp lazy_main; LAZY=$BUILT
build ms2/ms2_main.cpp ms2_main; MS2=$BUILT

OUT=$RESULTS/ms_times
TABLE=$OUT/ms_times.csv
mkdir -p "$OUT/tables"
[ -s "$TABLE" ] || echo "dataset,reference,files,input_bytes,lrf_build_ms,v2_build_ms,sassy_build_ms,lrf_ms,lazy_v2_ms,lazy_sassy_ms,lazy_ok,sampled_v2_ms,sampled_sassy_ms,sampled_ok" > "$TABLE"

# run TAG PROGRAM ARGS...: the program's stdout in TAG.out, stderr in TAG.log.
run() {
  local tag=$1; shift
  ${NUMA_CMD[@]+"${NUMA_CMD[@]}"} timeout "${MS_TIMEOUT:-6h}" "$@" > "$tag.out" 2> "$tag.log" ||
    echo "    $(basename "$tag"): exit $? (see $tag.log)" >&2
}
value() { grep "^$2=" "$1" 2>/dev/null | cut -d= -f2 || true; }

for ds in ${@:-$DATASETS}; do
  list=$(input_list "$ds")
  refs=$(dataset_refs "$ds")
  [ "${ALL_REFS:-0}" = 1 ] || refs=$(printf '%s\n' "$refs" | head -n 1)
  for dir in $refs; do
    name=$(basename "$dir")
    ref=$dir/$name.txt
    if grep -q "^$ds,$name," "$TABLE" && [ "${FORCE:-0}" != 1 ]; then
      echo "=== $ds / $name: done (FORCE=1 to time it again)"
      continue
    fi
    base=$OUT/${ds}_${name}
    common=(--reference "$ref" --suffix-array "$ref.sa" --filenames "$list"
            --table "$OUT/tables/${ds}_${name}" --quiet)

    echo "=== $ds / $name: lazy ($(date '+%F %T'))"
    run "$base.lazy" "$LAZY" "${common[@]}" --results "$base.lazy.csv"
    echo "=== $ds / $name: sampled ($(date '+%F %T'))"
    run "$base.ms2" "$MS2" "${common[@]}" --no-chain --results "$base.ms2.csv"
    rm -f "$OUT/tables/${ds}_${name}"*  # rebuilt every run

    l=$base.lazy.out m=$base.ms2.out
    files=$(value "$l" files)
    sampled_ok=NO
    [ -n "$files" ] && [ "$(value "$m" sampled-ms-v2_equal_files)" = "$files" ] &&
      [ "$(value "$m" sampled-ms-sassy_equal_files)" = "$files" ] && sampled_ok=YES
    grep -v "^$ds,$name," "$TABLE" > "$TABLE.tmp" || true
    echo "$ds,$name,$files,$(value "$l" input_bytes),$(value "$l" lrf_build_ms),$(value "$l" v2_build_ms),$(value "$l" sassy_build_ms),$(value "$l" lrf_ms),$(value "$l" lazy-v2_ms),$(value "$l" lazy-sassy_ms),$(value "$l" all_ok),$(value "$m" sampled-ms-v2_ms),$(value "$m" sampled-ms-sassy_ms),$sampled_ok" >> "$TABLE.tmp"
    mv "$TABLE.tmp" "$TABLE"
    echo "    lrf-ms $(value "$l" lrf_ms) ms; lazy-v2 $(value "$l" lazy-v2_ms), lazy-sassy" \
         "$(value "$l" lazy-sassy_ms) (exact: $(value "$l" all_ok)); sampled-ms-v2" \
         "$(value "$m" sampled-ms-v2_ms), sampled-ms-sassy $(value "$m" sampled-ms-sassy_ms)" \
         "(exact: $sampled_ok)"
  done
done
echo "table: $TABLE"
