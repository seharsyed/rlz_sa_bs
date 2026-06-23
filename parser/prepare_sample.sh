#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

RAW_DIR="/path/to/raw_fasta_files"
DATA_DIR="/path/to/output_data_folder"

INPUT_DIR="$DATA_DIR/input"
REF_TXT="$DATA_DIR/reference.txt"
REF_PLAIN="$DATA_DIR/reference.plain"
INPUT_LIST="$DATA_DIR/input_list.txt"

mkdir -p "$INPUT_DIR"

echo "Cleaning FASTA files..."
for raw in "$RAW_DIR"/*; do
  [ -f "$raw" ] || continue
  base="$(basename "$raw")"
  sed '/^>/d' "$raw" | tr -d '[:space:]' > "$INPUT_DIR/$base"
done

mapfile -t files < <(find "$INPUT_DIR" -maxdepth 1 -type f | sort)

if [ "${#files[@]}" -lt 101 ]; then
  echo "ERROR: need at least 101 cleaned files"
  exit 1
fi

REF_SOURCE="${files[0]}"

printf "%s\n" "$REF_SOURCE" > "$REF_TXT"
cp "$REF_SOURCE" "$REF_PLAIN"
printf '$' >> "$REF_PLAIN"

printf "%s\n" "${files[@]:1:100}" > "$INPUT_LIST"

echo "Done."
echo "Reference: $REF_PLAIN"
echo "Input list: $INPUT_LIST"
echo "Input count: $(wc -l < "$INPUT_LIST")"