#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Usage
# ============================================================
#
# The number of datasets is specified first.
#
# For each dataset, provide:
#
#   <name> <reference> <suffix_array> <input_list>
#
# followed by the output directory.
#
# General form:
#
#   ./run_rlz_experiments.sh \
#       <number_of_datasets> \
#       <name1> <reference1> <suffix_array1> <input_list1> \
#       <name2> <reference2> <suffix_array2> <input_list2> \
#       ... \
#       <output_dir>
#
# Example with one dataset:
#
#   ./run_rlz_experiments.sh \
#       1 \
#       ecoli \
#       /data/ecoli/reference.plain \
#       /data/ecoli/reference.plain.sa \
#       /data/ecoli/inputs.txt \
#       ./results
#
# Example with two datasets:
#
#   ./run_rlz_experiments.sh \
#       2 \
#       ecoli \
#       /data/ecoli/reference.plain \
#       /data/ecoli/reference.plain.sa \
#       /data/ecoli/inputs.txt \
#       chr19 \
#       /data/chr19/reference.plain \
#       /data/chr19/reference.plain.sa \
#       /data/chr19/inputs.txt \
#       ./results
#
# ============================================================


# ============================================================
# Locate this script
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ============================================================
# Check minimum arguments
# ============================================================

if [[ $# -lt 5 ]]; then
    echo "Usage:"
    echo
    echo "  $0 <number_of_datasets> \\"
    echo "     <name1> <reference1> <suffix_array1> <input_list1> \\"
    echo "     [<name2> <reference2> <suffix_array2> <input_list2> ...] \\"
    echo "     <output_dir>"
    echo
    echo "Example:"
    echo "  $0 1 ecoli ref.plain ref.plain.sa inputs.txt ./results"
    exit 1
fi


# ============================================================
# Parse number of datasets
# ============================================================

NUM_DATASETS="$1"
shift

if ! [[ "$NUM_DATASETS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: Number of datasets must be a positive integer."
    exit 1
fi


# ============================================================
# Check argument count
#
# After removing NUM_DATASETS:
#
#   4 arguments per dataset
#   + 1 output directory
# ============================================================

EXPECTED_ARGS=$((NUM_DATASETS * 4 + 1))

if [[ $# -ne "$EXPECTED_ARGS" ]]; then
    echo "ERROR: Expected $EXPECTED_ARGS arguments after the dataset count."
    echo "       Received $#."
    echo
    echo "Each dataset requires:"
    echo "  <name> <reference> <suffix_array> <input_list>"
    echo
    echo "Usage:"
    echo "  $0 <number_of_datasets> \\"
    echo "     <name1> <reference1> <suffix_array1> <input_list1> \\"
    echo "     ... \\"
    echo "     <output_dir>"
    exit 1
fi


# ============================================================
# Parse datasets
# ============================================================

declare -a DATASET_NAMES
declare -a DATASET_REFS
declare -a DATASET_SAS
declare -a DATASET_FILES

for ((i=0; i<NUM_DATASETS; i++)); do
    DATASET_NAMES[$i]="$1"
    DATASET_REFS[$i]="$2"
    DATASET_SAS[$i]="$3"
    DATASET_FILES[$i]="$4"

    shift 4
done

OUT="$1"

mkdir -p "$OUT"

DATE=$(date +%Y%m%d)


# ============================================================
# Source files and executables
# ============================================================

CACHE_CPP="$SCRIPT_DIR/compare_rlz_cache.cpp"
PT16_V2_CPP="$SCRIPT_DIR/pt16_main_v2.cpp"

CACHE="$SCRIPT_DIR/cache"
PT16_V2="$SCRIPT_DIR/pt16_v2"


# ============================================================
# Validation helpers
# ============================================================

check_file() {
    local description="$1"
    local path="$2"

    if [[ ! -f "$path" ]]; then
        echo "ERROR: $description does not exist:"
        echo "  $path"
        exit 1
    fi
}

check_executable() {
    local description="$1"
    local path="$2"

    if [[ ! -x "$path" ]]; then
        echo "ERROR: $description is not found or not executable:"
        echo "  $path"
        exit 1
    fi
}


# ============================================================
# Validate source files
# ============================================================

echo "========================================"
echo "Validating RLZ experiment inputs"
echo "========================================"

check_file "cache source" "$CACHE_CPP"
check_file "PT16 V2 source" "$PT16_V2_CPP"


# ============================================================
# Validate dataset inputs
# ============================================================

for ((i=0; i<NUM_DATASETS; i++)); do
    NAME="${DATASET_NAMES[$i]}"
    REF="${DATASET_REFS[$i]}"
    SA="${DATASET_SAS[$i]}"
    FILES="${DATASET_FILES[$i]}"

    echo
    echo "Dataset: $NAME"

    check_file "Reference for $NAME" "$REF"
    check_file "Suffix array for $NAME" "$SA"
    check_file "Input list for $NAME" "$FILES"
done

echo
echo "All inputs validated successfully."


# ============================================================
# Build executables
# ============================================================

echo
echo "========================================"
echo "Building RLZ executables"
echo "========================================"

echo
echo "Building cache..."
g++ -std=c++20 -O3 "$CACHE_CPP" -o "$CACHE"

echo
echo "Building PT16 V2..."
g++ -std=c++20 -O3 "$PT16_V2_CPP" -o "$PT16_V2"

echo
echo "Build complete."

check_executable "cache executable" "$CACHE"
check_executable "PT16 V2 executable" "$PT16_V2"


# ============================================================
# Print configuration
# ============================================================

echo
echo "========================================"
echo "RLZ EXPERIMENT CONFIGURATION"
echo "========================================"

echo
echo "Number of datasets: $NUM_DATASETS"

for ((i=0; i<NUM_DATASETS; i++)); do
    echo
    echo "Dataset $((i + 1)): ${DATASET_NAMES[$i]}"
    echo "  Reference:    ${DATASET_REFS[$i]}"
    echo "  Suffix array: ${DATASET_SAS[$i]}"
    echo "  Input list:   ${DATASET_FILES[$i]}"
done

echo
echo "Source directory:"
echo "  $SCRIPT_DIR"

echo
echo "Executables:"
echo "  cache:   $CACHE"
echo "  PT16 V2: $PT16_V2"

echo
echo "Output directory:"
echo "  $OUT"


# ============================================================
# RLZ experiment runner
# ============================================================

run_dataset() {
    local NAME="$1"
    local REF="$2"
    local SA="$3"
    local FILES="$4"

    echo
    echo "########################################"
    echo "RLZ EXPERIMENT: $NAME"
    echo "########################################"

    # --------------------------------------------------------
    # 1. BASELINE
    # --------------------------------------------------------

    echo
    echo "[1/2] BASELINE"

    "$CACHE" \
      --reference "$REF" \
      --suffix-array "$SA" \
      --filenames "$FILES" \
      --bucket-divisor 8 \
      --mode baseline \
      | tee "$OUT/${NAME}_baseline_${DATE}.txt"


    # --------------------------------------------------------
    # 2. PT16 V2
    # --------------------------------------------------------

    echo
    echo "[2/2] PT16 V2"

    "$PT16_V2" \
      --reference "$REF" \
      --suffix-array "$SA" \
      --filenames "$FILES" \
      --pt16-table "$OUT/${NAME}_pt16_v2_${DATE}.bin" \
      --results "$OUT/${NAME}_pt16_v2_${DATE}.csv" \
      2>&1 | tee "$OUT/${NAME}_pt16_v2_${DATE}.log"


    echo
    echo "$NAME RLZ EXPERIMENT COMPLETE"
}


# ============================================================
# Run all requested datasets
# ============================================================

for ((i=0; i<NUM_DATASETS; i++)); do
    run_dataset \
        "${DATASET_NAMES[$i]}" \
        "${DATASET_REFS[$i]}" \
        "${DATASET_SAS[$i]}" \
        "${DATASET_FILES[$i]}"
done


# ============================================================
# Done
# ============================================================

echo
echo "========================================"
echo "ALL RLZ EXPERIMENTS COMPLETE"
echo "========================================"
echo "Date: $DATE"
echo "Datasets processed: $NUM_DATASETS"
echo "Results: $OUT"
echo "========================================"