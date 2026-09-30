#!/usr/bin/env bash
# Picks a dataset's inputs and references from ONE multi-FASTA file (each
# record, >name, is a genome), into config.sh's layout:
#
#   COUNT random records          -> raw/<dataset>/inputs/<name>.fa
#   REFS random other records     -> raw/<dataset>/refs/<name>.fa
#
# Only the chosen records are written (the file itself is only read). The
# draw is random but reproducible (SEED): first the inputs, then the
# references among the records not drawn as inputs, so the two never overlap.
# It is recorded in raw/<dataset>/sample.txt. An existing sample is not
# replaced (delete raw/<dataset>/inputs and refs to draw again). Then run
# clean.sh DATASET.
#
#   ./sample_fasta.sh DATASET FASTA [--count 100] [--refs 3] [--seed 1]
#
# FASTA may be gzipped. A record's name is the first word of its header;
# characters other than letters, digits, '.', '-' and '_' become '_' in the
# file name.

set -euo pipefail
source "$(dirname "$0")/config.sh"

usage() { sed -n '2,/^$/p' "$0"; exit 1; }
[ $# -ge 2 ] || usage
DS=$1 FASTA=$2
shift 2
COUNT=100 NREFS=3 SEED=1
while [ $# -gt 0 ]; do
  case "$1" in
    --count) COUNT=$2; shift 2 ;;
    --refs) NREFS=$2; shift 2 ;;
    --seed) SEED=$2; shift 2 ;;
    *) usage ;;
  esac
done
[ -f "$FASTA" ] || { echo "no such file: $FASTA" >&2; exit 1; }

if [ -n "$(ls -A "$RAW/$DS/inputs" 2>/dev/null)" ] || [ -n "$(ls -A "$RAW/$DS/refs" 2>/dev/null)" ]; then
  echo "$RAW/$DS/inputs or refs is not empty: keeping the existing sample" >&2
  echo "(delete $RAW/$DS/inputs and $RAW/$DS/refs to draw a new one)" >&2
  exit 1
fi
mkdir -p "$RAW/$DS/inputs" "$RAW/$DS/refs"

python3 - "$FASTA" "$RAW/$DS" "$COUNT" "$NREFS" "$SEED" <<'PY'
import gzip, random, re, sys, time
fasta, target, count, nrefs, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])

def open_fasta():
    return gzip.open(fasta, "rb") if fasta.endswith(".gz") else open(fasta, "rb", buffering=1 << 24)

def safe(name):
    return re.sub(rb"[^A-Za-z0-9._-]", b"_", name).decode()

# Pass 1: the records' names, in file order.
start = time.time()
names = []
with open_fasta() as f:
    for line in f:
        if line.startswith(b">"):
            names.append(line[1:].split()[0] if line[1:].split() else b"record%d" % len(names))
print(f"=== {len(names)} records in {fasta} ({time.time() - start:.0f} s)")
if len(set(names)) != len(names):
    sys.exit("ERROR: record names are not unique")
if count + nrefs > len(names):
    sys.exit(f"ERROR: {count} inputs + {nrefs} references > {len(names)} records")

# The draw: inputs first, then references among the rest.
rng = random.Random(seed)
indices = list(range(len(names)))
inputs = set(rng.sample(indices, count))
rest = [i for i in indices if i not in inputs]
refs = set(rng.sample(rest, nrefs))

# Pass 2: write the chosen records.
start = time.time()
out = None
k = -1
with open_fasta() as f:
    for line in f:
        if line.startswith(b">"):
            if out:
                out.close()
                out = None
            k += 1
            if k in inputs or k in refs:
                folder = "inputs" if k in inputs else "refs"
                out = open(f"{target}/{folder}/{safe(names[k])}.fa", "wb", buffering=1 << 24)
        if out:
            out.write(line)
if out:
    out.close()
print(f"    wrote {count} inputs and {nrefs} references ({time.time() - start:.0f} s)")

with open(f"{target}/sample.txt", "w") as s:
    s.write(f"from {fasta}: {len(names)} records, seed {seed}\n")
    s.write(f"{count} inputs, then {nrefs} references among the other records\n")
    s.write("references:\n")
    for i in sorted(refs):
        s.write(f"  {names[i].decode()}\n")
    s.write("inputs:\n")
    for i in sorted(inputs):
        s.write(f"  {names[i].decode()}\n")
print("=== references:")
for i in sorted(refs):
    print(f"    {names[i].decode()}")
PY
echo "=== recorded in $RAW/$DS/sample.txt; next: ./clean.sh $DS"
