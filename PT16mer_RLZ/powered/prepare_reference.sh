#!/usr/bin/env bash
# Prepares the powered index of one cleaned reference for the RLZ suite
# (rlz_suite): REF_four.bwt (+ REF_four_data.bwt), via
# PFP-eBWT -> make_bwt (RLZ_powered's; it builds the powered index directly).
#
# The suffix array the other parsers need (REF.sa) is not made here: use the
# one you already have. The PT16 tables and lrf-ms's structures are built by
# the suite itself.
#
#   ./prepare_reference.sh [--reversed] REF [OUTDIR]
#
# --reversed: index the REVERSED reference instead, for the forward powered
# variants (powered-fwd-escape, powered-pt16-fwd-escape: the greedy
# left-to-right parse): writes REF.rev (REF reversed) to OUTDIR and builds
# REF.rev_four.bwt (+ REF.rev_four_data.bwt) from it. The parsers still take
# the original REF as --reference (and reverse it in memory themselves).
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
  sed -n '2,/^$/p' "$0"
  exit 1
fi

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${REPO:-$(cd "$HERE/../.." && pwd)}
PFP=${PFP:-$HOME/proj/ebwt/PFP-eBWT/build}
POWERED=$REPO/RLZ_powered

REVERSED=0
if [ "$1" = "--reversed" ]; then
  REVERSED=1
  shift
fi

SRC=$(realpath "$1")
OUTDIR=$(realpath "${2:-$(dirname "$SRC")}")
mkdir -p "$OUTDIR"
NAME=$(basename "$SRC")
if [ "$REVERSED" -eq 1 ]; then
  REF=$OUTDIR/$NAME.rev   # written in step 0b
else
  REF=$OUTDIR/$NAME
  [ "$SRC" = "$REF" ] || ln -sf "$SRC" "$REF"   # keep all outputs together
fi

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

# ---------- 0b. The reversed reference (--reversed) ----------

if [ "$REVERSED" -eq 1 ]; then
  step "0b. reversed reference $REF"
  if [ -s "$REF" ] && [ "$REF" -nt "$SRC" ]; then
    echo "    exists"
  else
    perl -0777 -pe '$_ = reverse $_' "$SRC" > "$REF"
  fi
  [ "$(stat -c %s "$REF")" -eq "$N" ] ||
    { echo "ERROR: $REF has another length than $SRC" >&2; exit 1; }
fi

# ---------- 1. FASTA framing for PFP-eBWT ----------

step "1. FASTA $REF.fa"
if [ -s "$REF.fa" ]; then
  echo "    exists"
else
  # header, the sequence, a newline ending it, then an empty second record
  # (the format of RLZ_powered's example)
  { echo ">reference"; cat "$REF"; echo; echo ">empty"; } > "$REF.fa"
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
# make_bwt regenerates RLZ_powered's alphabet header for this eBWT and
# rebuilds the tools, but builds the index with the binary that was already
# running; so it is run twice (the second run uses the rebuilt tools), and a
# third time if it still reports 0 dense blocks (RLZ_powered's README).

step "3. powered index ${REF}_four.bwt"
if [ -s "${REF}_four.bwt" ] && [ -s "${REF}_four_data.bwt" ]; then
  echo "    exists (delete it to rebuild)"
else
  cd "$POWERED"
  make
  ./make_bwt -i "$REF.fa.ebwt" -sa "$REF.fa.gca" -o "${REF}_four.bwt" > /dev/null 2>&1
  OUT=$(./make_bwt -i "$REF.fa.ebwt" -sa "$REF.fa.gca" -o "${REF}_four.bwt" 2>&1)
  echo "$OUT" | grep "dense blocks\|took" | sed 's/^/    /'
  if echo "$OUT" | grep -q " 0 dense blocks"; then
    ./make_bwt -i "$REF.fa.ebwt" -sa "$REF.fa.gca" -o "${REF}_four.bwt" 2>&1 \
      | grep "dense blocks" | sed 's/^/    third run: /'
  fi
fi

# ---------- Done ----------

step "done: run the suite with"
cat <<EOF
    ./rlz_suite --reference $SRC --suffix-array <its .sa file> \\
                $([ "$REVERSED" -eq 1 ] && echo --powered-fwd-index || echo --powered-index) ${REF}_four.bwt --filenames LIST ...

(LIST should not contain this reference itself.)
EOF
