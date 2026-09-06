#!/bin/bash

set -euo pipefail

DATE=$(date +%Y%m%d)
OUT=/home/sehar/PhD/RLZResults

ECOLI=/home/sehar/PhD/rlz_sa_bs/Data/ecoli_16mer
CHR19=/home/sehar/PhD/rlz_sa_bs/Data/chr19_16mer

mkdir -p "$OUT"

echo "========================================"
echo "Building executables"
echo "========================================"

g++ -std=c++20 -O3 compare_rlz_cache.cpp -o cache
g++ -std=c++20 -O3 pt16_main.cpp -o pt16
g++ -std=c++20 -O3 pt16_main_v2.cpp -o pt16_v2

echo
echo "Build complete."
echo

run_dataset () {
    NAME=$1
    DATA=$2

    REF="$DATA/reference_acgt.plain"
    SA="$DATA/reference_acgt.plain.sa"
    FILES="$DATA/cleaned_input_list.txt"

    echo
    echo "########################################"
    echo "DATASET: $NAME"
    echo "########################################"

    echo
    echo "[1/4] SA BASELINE"
    ./cache \
      --reference "$REF" \
      --suffix-array "$SA" \
      --filenames "$FILES" \
      --bucket-divisor 8 \
      --mode baseline \
      | tee "$OUT/${NAME}_baseline_${DATE}.txt"

    echo
    echo "[2/4] CACHE B=8"
    ./cache \
      --reference "$REF" \
      --suffix-array "$SA" \
      --filenames "$FILES" \
      --bucket-divisor 8 \
      --mode cached \
      | tee "$OUT/${NAME}_cache_b8_${DATE}.txt"

    echo
    echo "[3/4] PT16 V1"
    ./pt16 \
      --reference "$REF" \
      --suffix-array "$SA" \
      --filenames "$FILES" \
      --pt16-table "$DATA/pt16_v1_${DATE}.bin" \
      --results "$OUT/${NAME}_pt16_v1_${DATE}.csv" \
      2>&1 | tee "$OUT/${NAME}_pt16_v1_${DATE}.log"

    echo
    echo "[4/4] PT16 V2"
    ./pt16_v2 \
      --reference "$REF" \
      --suffix-array "$SA" \
      --filenames "$FILES" \
      --pt16-table "$DATA/pt16_v2_${DATE}.bin" \
      --results "$OUT/${NAME}_pt16_v2_${DATE}.csv" \
      2>&1 | tee "$OUT/${NAME}_pt16_v2_${DATE}.log"

    echo
    echo "$NAME COMPLETE"
}

run_dataset "ecoli" "$ECOLI"
run_dataset "chr19" "$CHR19"

echo
echo "========================================"
echo "ALL SA EXPERIMENTS COMPLETE"
echo "Date: $DATE"
echo "Results: $OUT"
echo "========================================"
