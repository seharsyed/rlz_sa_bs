#!/usr/bin/env bash
# Picks a dataset's references and a random sample of its inputs from an
# archive of genomes (config.sh's layout):
#
#   ARCHIVE            -> raw/<dataset>/pool/       (extracted once)
#   the given refs     -> raw/<dataset>/refs/       (symbolic links into pool)
#   COUNT other files  -> raw/<dataset>/inputs/     (symbolic links into pool)
#
# The sample is random but reproducible (SEED), never contains a reference,
# and is recorded in raw/<dataset>/sample.txt. An existing sample is not
# replaced (delete raw/<dataset>/inputs and refs to draw again). Then run
# clean.sh DATASET.
#
#   ./sample.sh DATASET ARCHIVE --refs NAME1,NAME2,NAME3 [--count 100] [--seed 1]
#   ./sample.sh DATASET ARCHIVE --random-refs 3 [--count 100] [--seed 1]
#
# --random-refs N: the references are drawn at random too: first the COUNT
# inputs, then N references among the other genomes (as sample_fasta.sh).
# --glob PATTERN: the genome files are the files matching PATTERN (e.g.
# '*.txt' for plain sequence files) instead of the FASTA endings below.
#
# ARCHIVE: .zip, .tar, .tar.gz/.tgz, or a directory already holding the
# genomes. The genomes are the files ending in .fa/.fasta/.fna/.fas (optionally
# .gz) anywhere below it, except CDS, RNA and protein FASTA files (NCBI's
# cds_from_genomic.fna, rna.fna, ...). A reference NAME (e.g. an accession,
# GCA_903819205.2; trailing dots are ignored) matches a genome file whose name
# is NAME, or starts with NAME_ or NAME. (NCBI's GCA_..._<assembly>_genomic.fna),
# or that lies in a folder named NAME; it must match exactly one genome.

set -euo pipefail
source "$(dirname "$0")/config.sh"

usage() { sed -n '2,/^$/p' "$0"; exit 1; }
[ $# -ge 2 ] || usage
DS=$1 ARCHIVE=$2
shift 2
REFS="" RANDOM_REFS=0 GLOB="" COUNT=100 SEED=1
while [ $# -gt 0 ]; do
  case "$1" in
    --refs) REFS=$2; shift 2 ;;
    --random-refs) RANDOM_REFS=$2; shift 2 ;;
    --glob) GLOB=$2; shift 2 ;;
    --count) COUNT=$2; shift 2 ;;
    --seed) SEED=$2; shift 2 ;;
    *) usage ;;
  esac
done
if [ -n "$REFS" ] && [ "$RANDOM_REFS" -gt 0 ]; then
  echo "give --refs or --random-refs, not both" >&2; exit 1
fi
[ -n "$REFS" ] || [ "$RANDOM_REFS" -gt 0 ] ||
  { echo "--refs NAMES or --random-refs N is required" >&2; exit 1; }

POOL=$RAW/$DS/pool
mkdir -p "$RAW/$DS"

# ---------- 1. extract ----------
if [ -d "$ARCHIVE" ]; then
  POOL=$(cd "$ARCHIVE" && pwd)
  echo "=== genomes in $POOL"
elif [ -d "$POOL" ]; then
  echo "=== $POOL exists (delete it to extract again)"
else
  echo "=== extracting $ARCHIVE into $POOL"
  # Into a temporary folder, renamed to pool/ only when everything was
  # extracted: a failed or interrupted extraction leaves no pool/ behind.
  TMP_POOL=$POOL.partial
  rm -rf "$TMP_POOL"
  mkdir -p "$TMP_POOL"
  case "$ARCHIVE" in
    *.zip)
      # Python's zipfile: checks every file's CRC and handles archives over
      # 4 GB (some unzip versions report bad CRCs on those).
      if ! python3 - "$ARCHIVE" "$TMP_POOL" <<'PY'
import sys, zipfile
archive, target = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(archive) as z:
    members = z.infolist()
    for k, member in enumerate(members, 1):
        z.extract(member, target)  # raises BadZipFile on a CRC mismatch
        if k % 500 == 0 or k == len(members):
            print(f"    {k}/{len(members)} files", flush=True)
PY
      then
        echo "ERROR: extracting $ARCHIVE failed (a damaged archive? see above)" >&2
        rm -rf "$TMP_POOL"
        exit 1
      fi
      ;;
    *.tar.gz|*.tgz) tar -xzf "$ARCHIVE" -C "$TMP_POOL" ;;
    *.tar) tar -xf "$ARCHIVE" -C "$TMP_POOL" ;;
    *) echo "unknown archive type: $ARCHIVE" >&2; rm -rf "$TMP_POOL"; exit 1 ;;
  esac
  mv "$TMP_POOL" "$POOL"
fi

GENOMES=$RAW/$DS/genomes.txt
if [ -n "$GLOB" ]; then
  find -L "$POOL" -type f -name "$GLOB" | sort > "$GENOMES"
else
  find -L "$POOL" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fas' \
       -o -name '*.fa.gz' -o -name '*.fasta.gz' -o -name '*.fna.gz' -o -name '*.fas.gz' \) \
       ! -name '*cds_from_genomic*' ! -name '*rna_from_genomic*' ! -name 'rna.f*' \
       ! -name 'cds.f*' ! -name '*protein*' \
    | sort > "$GENOMES"
fi
echo "    $(wc -l < "$GENOMES") genome files (listed in $GENOMES)"

# ---------- 2. the references ----------
if [ -d "$RAW/$DS/inputs" ] && [ -n "$(ls -A "$RAW/$DS/inputs" 2>/dev/null)" ]; then
  echo "=== $RAW/$DS/inputs is not empty: keeping the existing sample" >&2
  echo "    (delete $RAW/$DS/inputs and $RAW/$DS/refs to draw a new one)" >&2
  exit 1
fi
mkdir -p "$RAW/$DS/refs" "$RAW/$DS/inputs"

REF_FILES=()
if [ "$RANDOM_REFS" -gt 0 ]; then
  # The inputs first, then the references among the other genomes.
  python3 - "$GENOMES" "$COUNT" "$RANDOM_REFS" "$SEED" \
      "$RAW/$DS/sample_files.txt" "$RAW/$DS/refs.txt" <<'PY'
import random, sys
genomes = [l.strip() for l in open(sys.argv[1]) if l.strip()]
count, nrefs, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
if count + nrefs > len(genomes):
    sys.exit(f"only {len(genomes)} genomes for {count} inputs and {nrefs} references")
rng = random.Random(seed)
inputs = rng.sample(genomes, count)
chosen = set(inputs)
refs = rng.sample([g for g in genomes if g not in chosen], nrefs)
open(sys.argv[5], "w").write("".join(g + "\n" for g in sorted(inputs)))
open(sys.argv[6], "w").write("".join(g + "\n" for g in sorted(refs)))
PY
  while read -r f; do REF_FILES+=("$f"); done < "$RAW/$DS/refs.txt"
  REF_NAMES=()
fi
IFS=',' read -r -a REF_NAMES <<< "$REFS"
for name in ${REF_NAMES[@]+"${REF_NAMES[@]}"}; do
  name=$(echo "$name" | tr -d '[:space:]')
  name=${name%%.}          # a trailing dot (punctuation)
  while [ "${name%.}" != "$name" ]; do name=${name%.}; done
  matches=$(while read -r f; do
              b=$(basename "$f")
              parent=$(basename "$(dirname "$f")")
              if [ "$b" = "$name" ] || [ "$(file_stem "$f")" = "$(file_stem "$name")" ] ||
                 [[ "$b" == "$name"_* ]] || [[ "$b" == "$name".* ]] ||
                 [ "$parent" = "$name" ]; then
                echo "$f"
              fi
            done < "$GENOMES")
  count=$(printf '%s' "$matches" | grep -c . || true)
  if [ "$count" -ne 1 ]; then
    echo "reference '$name' matches $count files:" >&2
    printf '%s\n' "$matches" | head -5 >&2
    rm -rf "$RAW/$DS/refs" "$RAW/$DS/inputs"
    exit 1
  fi
  REF_FILES+=("$matches")
done
echo "=== references"
for f in "${REF_FILES[@]}"; do
  ln -s "$f" "$RAW/$DS/refs/$(basename "$f")"
  echo "    $(basename "$f")"
done

# ---------- 3. the sample ----------
if [ "$RANDOM_REFS" -eq 0 ]; then
printf '%s\n' "${REF_FILES[@]}" > "$RAW/$DS/refs.txt"
python3 - "$GENOMES" "$RAW/$DS/refs.txt" "$COUNT" "$SEED" > "$RAW/$DS/sample_files.txt" <<'PY'
import random, sys
genomes = [l.strip() for l in open(sys.argv[1]) if l.strip()]
refs = {l.strip() for l in open(sys.argv[2]) if l.strip()}
candidates = [g for g in genomes if g not in refs]
count, seed = int(sys.argv[3]), int(sys.argv[4])
if count > len(candidates):
    sys.exit(f"only {len(candidates)} candidates for a sample of {count}")
for g in sorted(random.Random(seed).sample(candidates, count)):
    print(g)
PY
fi
while read -r f; do
  ln -s "$f" "$RAW/$DS/inputs/$(basename "$f")"
done < "$RAW/$DS/sample_files.txt"

{
  echo "dataset $DS, archive $ARCHIVE, $(date)"
  if [ "$RANDOM_REFS" -gt 0 ]; then
    echo "seed $SEED, $COUNT inputs from $(wc -l < "$GENOMES") genomes, then $RANDOM_REFS random references among the others"
  else
    echo "seed $SEED, $COUNT inputs from $(wc -l < "$GENOMES") genomes, references excluded"
  fi
  echo "references:"; sed 's/^/  /' "$RAW/$DS/refs.txt"
  echo "inputs:"; sed 's/^/  /' "$RAW/$DS/sample_files.txt"
} > "$RAW/$DS/sample.txt"
echo "=== $(ls "$RAW/$DS/inputs" | wc -l) inputs sampled (seed $SEED; recorded in $RAW/$DS/sample.txt)"
echo "next: ./clean.sh $DS"
