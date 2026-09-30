# Helpers shared by verify.sh and run.sh (sourced after config.sh).

# The command prefix of every run (NUMA_CMD): numactl --cpunodebind=N
# --membind=N when config.sh's NUMA is set, taskset -c CPUS when CPUS is set;
# empty otherwise.
NUMA_CMD=()
if [ -n "${NUMA:-}" ]; then
  command -v numactl > /dev/null ||
    { echo "NUMA=$NUMA set, but numactl is not installed" >&2; exit 1; }
  NUMA_CMD=(numactl "--cpunodebind=$NUMA" "--membind=$NUMA")
  echo "(pinned to NUMA node $NUMA)"
fi
# CPU restriction (config.sh's CPUS), added to the same command prefix.
if [ -n "${CPUS:-}" ]; then
  command -v taskset > /dev/null ||
    { echo "CPUS=$CPUS set, but taskset is not installed" >&2; exit 1; }
  NUMA_CMD=(taskset -c "$CPUS" ${NUMA_CMD[@]+"${NUMA_CMD[@]}"})
  echo "(restricted to CPUs $CPUS)"
fi

# build_programs PROGRAM: builds (if needed) PT16mer_RLZ/PROGRAM.<host> with
# powered and varki (rlz_datasets.sh's build_if_needed), copies it to
# ROOT/bin (so that syncing the repository during a long run cannot remove the
# program being run) and sets PROGRAM_PATH to that copy.
build_programs() {
  local program=$1
  pushd "$PT16" > /dev/null
  # shellcheck disable=SC1091
  source ./rlz_datasets.sh
  build_if_needed "$program"
  popd > /dev/null
  mkdir -p "$ROOT/bin"
  cp -p "$PT16/$BINARY" "$ROOT/bin/$BINARY.tmp"
  mv "$ROOT/bin/$BINARY.tmp" "$ROOT/bin/$BINARY"
  PROGRAM_PATH=$ROOT/bin/$BINARY
  [ -f "$PT16/../RLZ-Varki/build/sdsl/lib/libsdsl.a" ] ||
    echo "WARNING: no sdsl for varki (run ../varki/build_sdsl.sh): varki is not built" >&2
}

# input_list DATASET: writes the list of the dataset's cleaned inputs
# (data/<dataset>/inputs/*.txt, sorted) to results/<dataset>/inputs.txt and
# prints its path.
input_list() {
  local ds=$1
  local list=$RESULTS/$ds/inputs.txt
  mkdir -p "$RESULTS/$ds"
  find "$DATA/$ds/inputs" -maxdepth 1 -type f -name '*.txt' | sort > "$list"
  [ -s "$list" ] || { echo "no inputs in $DATA/$ds/inputs" >&2; exit 1; }
  echo "$list"
}
