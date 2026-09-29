#!/usr/bin/env bash
# Runs the RLZ parsing suite (rlz_suite) on the server's data.
#
#   ./run_rlz_suite.sh REFERENCE LIST [extra rlz_suite options]
#
# REFERENCE: ecoli | GCA_000179135.1 | GCA_001012175.1
# LIST:      the input list (one input path per line)
#
# The suffix arrays are in RLZ_DATA/cleaned-references/ (REF.sa). The E. coli
# reference is there too; the other references are in RLZ_DATA/cleaned_inputs/
# under the same name as their .sa file without the extension. The powered
# index is used when it exists (built with powered_rlz's tools, see
# powered/README.md); otherwise powered is skipped.
#
# Build rlz_suite first (from PT16mer_RLZ/):
#   g++ -std=c++2a -O3 -march=native -DNDEBUG -DWITH_POWERED \
#       -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 rlz_suite.cpp -o rlz_suite

set -euo pipefail

DATA=${DATA:-$HOME/proj/RLZ_DATA}
SA_DIR=$DATA/cleaned-references

if [ $# -lt 2 ]; then
  sed -n '2,20p' "$0"
  exit 1
fi

case "$1" in
  ecoli|GCA_000005845.2)
    REF=$SA_DIR/GCA_000005845.2_ASM584v2.cleaned
    SA=$REF.sa
    POWERED=$DATA/BWT_FM_files/GCA_000005845.2_ASM584v2.cleaned_four.bwt
    ;;
  GCA_000179135.1)
    REF=$DATA/cleaned_inputs/GCA_000179135.1_ASM17913v1.txt
    SA=$SA_DIR/GCA_000179135.1_ASM17913v1.txt.sa
    POWERED=$DATA/BWT_FM_files/GCA_000179135.1_ASM17913v1.txt_four.bwt
    ;;
  GCA_001012175.1)
    REF=$DATA/cleaned_inputs/GCA_001012175.1_CFSAN026787_02.0.txt
    SA=$SA_DIR/GCA_001012175.1_CFSAN026787_02.0.txt.sa
    POWERED=$DATA/BWT_FM_files/GCA_001012175.1_CFSAN026787_02.0.txt_four.bwt
    ;;
  *)
    echo "unknown reference: $1 (ecoli | GCA_000179135.1 | GCA_001012175.1)" >&2
    exit 1
    ;;
esac

LIST=$2
shift 2

for f in "$REF" "$SA" "$LIST"; do
  [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

POWERED_ARGS=()
if [ -f "$POWERED" ]; then
  POWERED_ARGS=(--powered-index "$POWERED")
else
  echo "(no powered index at $POWERED: powered is skipped)" >&2
fi

cd "$(dirname "$0")"
exec ./rlz_suite --reference "$REF" --suffix-array "$SA" --filenames "$LIST" \
  "${POWERED_ARGS[@]}" "$@"
