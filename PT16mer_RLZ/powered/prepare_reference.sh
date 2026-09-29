#!/usr/bin/env bash
# Prepares the powered index of one cleaned reference for the RLZ suite
# (rlz_suite): REF_four.bwt (+ REF_four_data.bwt), via
# PFP-eBWT -> make_bwt --rle -> transform -b.
#
# The suffix array the other parsers need (REF.sa) is not made here: use the
# one you already have. The PT16 tables and lrf-ms's structures are built by
# the suite itself.
#
#   ./prepare_reference.sh REF [OUTDIR]
#
# REF must be a cleaned reference: plain A/C/G/T, no header, no newline (our
# clean_reference / clean_inputs output). Everything is written to OUTDIR
# (default: next to REF). Steps already done (their output exists) are
# skipped, so the script can be rerun after a failure.
#
# Paths, overridable from the environment:
#   REPO   the rlz_sa_bs checkout            (default: two levels above here)
#   PFP    PFP-eBWT's build directory, with pfpebwt
#                                            (default: ~/proj/ebwt/PFP-eBWT/build)

set -euo pipefail

if [ $# -lt 1 ]; then
  sed -n '2,22p' "$0"
  exit 1
fi

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${REPO:-$(cd "$HERE/../.." && pwd)}
PFP=${PFP:-$HOME/proj/ebwt/PFP-eBWT/build}
POWERED=$REPO/powered_rlz

SRC=$(realpath "$1")
OUTDIR=$(realpath "${2:-$(dirname "$SRC")}")
mkdir -p "$OUTDIR"
NAME=$(basename "$SRC")
REF=$OUTDIR/$NAME
[ "$SRC" = "$REF" ] || ln -sf "$SRC" "$REF"   # keep all outputs together

step() { echo; echo "=== $*"; }

# ---------- 0. The reference must be plain ACGT, no newline ----------

step "0. checking $SRC"
N=$(stat -c %s "$SRC")
OTHER=$(tr -d 'ACGT' < "$SRC" | head -c 1 | wc -c)
if [ "$OTHER" -ne 0 ]; then
  echo "ERROR: $SRC contains bytes other than A/C/G/T (clean it first)" >&2
  exit 1
fi
echo "    $N bytes, plain ACGT"

# ---------- 1. FASTA framing for PFP-eBWT ----------

step "1. FASTA $REF.fa"
if [ -s "$REF.fa" ]; then
  echo "    exists"
else
  # header, the sequence, a newline ending it, then an empty second record
  # (the format of powered_rlz's example)
  { echo ">reference"; cat "$SRC"; echo; echo ">empty"; } > "$REF.fa"
  echo "    written"
fi

# ---------- 2. Extended BWT + GCA ----------

step "2. PFP-eBWT: $REF.fa.ebwt, $REF.fa.gca"
if [ -s "$REF.fa.ebwt" ] && [ -s "$REF.fa.gca" ]; then
  echo "    exists"
else
  [ -x "$PFP/pfpebwt" ] || { echo "ERROR: no pfpebwt in $PFP (set PFP=...)" >&2; exit 1; }
  (cd "$PFP" && ./pfpebwt --GCA "$REF.fa")
fi
LEN=$(grep -i "Length of the eBWT" "$REF.fa.info" | grep -o '[0-9]*$')
echo "    eBWT length $LEN (reference $N)"
[ "$LEN" -eq "$N" ] || { echo "ERROR: eBWT length differs from the reference" >&2; exit 1; }

# ---------- 3. powered index ----------
# make_bwt --rle regenerates powered_rlz's alphabet header for this eBWT and
# rebuilds the tools, but builds the index with the binary that was already
# running; so it is run twice (the second run uses the rebuilt tools).

step "3. powered index ${REF}_four.bwt"
if [ -s "${REF}_four.bwt" ] && [ -s "${REF}_four_data.bwt" ]; then
  echo "    exists (delete it to rebuild)"
else
  cd "$POWERED"
  make
  ./make_bwt --rle -i "$REF.fa.ebwt" -sa "$REF.fa.gca" -o "$REF.bwt" > /dev/null
  OUT=$(./make_bwt --rle -i "$REF.fa.ebwt" -sa "$REF.fa.gca" -o "$REF.bwt" 2>&1)
  echo "$OUT" | grep "dense blocks" | sed 's/^/    /'
  if echo "$OUT" | grep -q " 0 dense blocks"; then
    ./make_bwt --rle -i "$REF.fa.ebwt" -sa "$REF.fa.gca" -o "$REF.bwt" 2>&1 \
      | grep "dense blocks" | sed 's/^/    third run: /'
  fi
  ./transform -b -i "$REF.bwt" -sa "$REF.fa.gca" -o "${REF}_four.txt" | tail -1
fi

# ---------- Done ----------

step "done: run the suite with"
cat <<EOF
    ./rlz_suite --reference $SRC --suffix-array <its .sa file> \\
                --powered-index ${REF}_four.bwt --filenames LIST ...

(LIST should not contain this reference itself.)
EOF
