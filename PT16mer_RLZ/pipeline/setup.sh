#!/usr/bin/env bash
# Sets up a new machine for the pipeline (config.sh): checks the
# prerequisites and builds every tool. Steps whose result exists are skipped,
# so it can be rerun after a failure.
#
#   1. sdsl-lite + divsufsort for Varki     (../varki/build_sdsl.sh)
#   2. RLZ_powered's tools (make_bwt, ...)  (make in ../../RLZ_powered)
#   3. PFP-eBWT                             (cloned and built in TOOLS)
#   4. the suffix array tools (sa, and psa: parallel) and the cleaning tool,
#      in TOOLS/bin
#   5. rlz_suite and rlz_parallel           (with powered and varki)
#
#   ./setup.sh                 everything (the preprocessing machine)
#   ./setup.sh --experiment    only what the runs need: steps 0, 1 and 5 (the
#                              experiment machine: no PFP-eBWT, make_bwt or
#                              cleaning tools -- the data comes prepared)

set -euo pipefail
source "$(dirname "$0")/config.sh"

EXPERIMENT=0
[ "${1:-}" = "--experiment" ] && EXPERIMENT=1

step() { echo; echo "=== $*"; }

# ---------- 0. prerequisites ----------
step "0. prerequisites"
missing=0
for tool in g++ make cmake git gzip python3; do
  if command -v "$tool" > /dev/null; then
    echo "    $tool: $(command -v "$tool")"
  else
    echo "    $tool: MISSING" >&2
    missing=1
  fi
done
[ "$missing" -eq 0 ] || { echo "install the missing tools first" >&2; exit 1; }
echo "    g++ version: $(g++ -dumpfullversion 2>/dev/null || g++ -dumpversion) (needs C++20: GCC 10 or newer)"
[ "$(uname -m)" = "x86_64" ] || echo "    WARNING: not x86-64: powered (and so the powered variants) will not build" >&2
command -v numactl > /dev/null && echo "    numactl: yes (NUMA pinning available)" ||
  echo "    numactl: no (runs cannot be pinned; set NUMA only if it is installed)"
python3 -c "import matplotlib, pandas" 2> /dev/null && echo "    python: matplotlib and pandas found" ||
  echo "    python: matplotlib and/or pandas missing (needed for the plots: pip install matplotlib pandas)"
mkdir -p "$TOOLS/bin"

# ---------- 1. sdsl-lite for Varki ----------
step "1. sdsl-lite (Varki)"
"$PT16/varki/build_sdsl.sh"

if [ "$EXPERIMENT" -eq 0 ]; then
# ---------- 2. RLZ_powered ----------
step "2. RLZ_powered tools"
(cd "$REPO/RLZ_powered" && make)

# ---------- 3. PFP-eBWT ----------
step "3. PFP-eBWT in $TOOLS/PFP-eBWT"
if [ -x "$PFP/pfpebwt" ]; then
  echo "    exists"
else
  [ -d "$TOOLS/PFP-eBWT" ] ||
    git clone --recursive https://github.com/davidecenzato/PFP-eBWT.git "$TOOLS/PFP-eBWT"
  mkdir -p "$PFP"
  (cd "$PFP" && cmake .. && make -j"$(nproc 2>/dev/null || echo 4)")
  [ -x "$PFP/pfpebwt" ] || { echo "ERROR: no $PFP/pfpebwt after the build" >&2; exit 1; }
fi

# ---------- 4. small tools ----------
step "4. suffix array and cleaning tools in $TOOLS/bin"
g++ -std=c++20 -O3 -march=native "$REPO/sa/suffix_array.cpp" -o "$TOOLS/bin/sa"
g++ -std=c++20 -O3 -march=native -fopenmp "$REPO/sa/parallel_suffix_array.cpp" -o "$TOOLS/bin/psa"
g++ -std=c++20 -O3 "$PT16/clean_reference.cpp" -o "$TOOLS/bin/clean_reference"
echo "    built sa, psa (parallel, preferred by prepare.sh), clean_reference"
fi

# ---------- 5. our programs ----------
step "5. rlz_suite and rlz_parallel"
(
  cd "$PT16"
  source scripts/rlz_datasets.sh
  build_if_needed rlz_suite
  build_if_needed rlz_parallel
  ls -la rlz_suite."$(hostname -s)" rlz_parallel."$(hostname -s)"
)

if [ "$EXPERIMENT" -eq 1 ]; then
  step "done: next ./verify.sh, then ./run.sh, then ./collect.sh"
else
  step "done: next ./clean.sh, then ./prepare.sh"
fi
