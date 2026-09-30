# The datasets of the RLZ benchmarks, in one place. Sourced by
# run_rlz_suite.sh and run_rlz_parallel.sh:
#
#   source rlz_datasets.sh
#   select_dataset NAME     # sets REF, SA, POWERED, LIST (LIST may be empty)
#
# REF      the cleaned reference (plain ACGT)
# SA       its suffix array (4-byte entries, no sentinel)
# POWERED  its powered index (REF_four.bwt; skipped if the file is missing)
# LIST     the default input list (one input path per line), or empty
#
# DATA is the data root; override it from the environment.

DATA=${DATA:-$HOME/proj/RLZ_DATA}

list_datasets() {
  echo "ecoli GCA_000179135.1 GCA_001012175.1"
}

select_dataset() {
  REF="" SA="" POWERED="" LIST=""
  case "$1" in
    ecoli|GCA_000005845.2)
      REF=$DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned
      SA=$REF.sa
      POWERED=$DATA/BWT_FM_files/GCA_000005845.2_ASM584v2.cleaned_four.bwt
      LIST=$DATA/ecoli-input-lists/input-list-100-lex.txt
      ;;
    GCA_000179135.1)
      REF=$DATA/cleaned_inputs/GCA_000179135.1_ASM17913v1.txt
      SA=$DATA/cleaned-references/GCA_000179135.1_ASM17913v1.txt.sa
      POWERED=$DATA/BWT_FM_files/GCA_000179135.1_ASM17913v1.txt_four.bwt
      ;;
    GCA_001012175.1)
      REF=$DATA/cleaned_inputs/GCA_001012175.1_CFSAN026787_02.0.txt
      SA=$DATA/cleaned-references/GCA_001012175.1_CFSAN026787_02.0.txt.sa
      POWERED=$DATA/BWT_FM_files/GCA_001012175.1_CFSAN026787_02.0.txt_four.bwt
      ;;
    # chr19 references: add them here once chosen, e.g.
    # chr19-HG00236)
    #   REF=$DATA/<cleaned chr19 dir>/HG00236.1.19.txt
    #   SA=$REF.sa
    #   POWERED=$DATA/BWT_FM_files/chr19/HG00236.1.19.txt_four.bwt
    #   LIST=$DATA/<chr19 list without this reference>
    #   ;;
    *)
      echo "unknown dataset: $1 (known: $(list_datasets))" >&2
      return 1
      ;;
  esac
}

# build_if_needed PROGRAM   (rlz_suite or rlz_parallel; run from PT16mer_RLZ/)
#
# Each machine gets its own binary, PROGRAM.<hostname> (set in BINARY): with a
# home directory shared between machines, one binary would be rebuilt by each
# machine in turn (and -march=native code may not run on another CPU). It is
# (re)built when missing or older than PROGRAM.cpp or any header here or in
# ../RLZ_powered/include -- with powered on x86-64, without elsewhere. A failed
# build stops the calling script. Two runs started at once on one machine can
# still race on a rebuild: start the second after the first has built.
build_if_needed() {
  local program=$1
  BINARY="$program.$(hostname -s)"
  local newer=""
  if [ ! -x "$BINARY" ]; then
    newer="(no binary for this machine)"
  else
    newer=$(find . ../RLZ_powered/include \( -name '*.hpp' -o -name "$program.cpp" \) \
              -newer "$BINARY" 2>/dev/null | head -1)
    [ -n "$newer" ] && echo "($newer changed since $BINARY was built)" >&2
  fi
  [ -z "$newer" ] && return 0

  if [ "$(uname -m)" = "x86_64" ]; then
    echo "(building $BINARY with powered)" >&2
    g++ -std=c++2a -O3 -march=native -pthread -DNDEBUG -DWITH_POWERED \
        -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384 "$program.cpp" -o "$BINARY"
  else
    echo "(building $BINARY without powered: not x86-64)" >&2
    g++ -std=c++20 -O3 -pthread "$program.cpp" -o "$BINARY"
  fi
}

# numa_setup   (reads NUMA from the environment)
#
# NUMA=N pins the run to NUMA node N: its threads run on node N's CPUs and
# all its memory (tables, index, inputs) is allocated on node N, so no access
# crosses to the other socket (numactl --cpunodebind=N --membind=N). Unset or
# empty: no pinning. Sets NUMA_CMD (the command prefix, an array) and
# NUMA_LABEL ("numaN" or "unpinned"), and prints the choice.
numa_setup() {
  NUMA_CMD=()
  NUMA_LABEL="unpinned"
  if [ -n "${NUMA:-}" ]; then
    case "$NUMA" in
      *[!0-9]*) echo "NUMA must be a node number, got '$NUMA'" >&2; exit 1 ;;
    esac
    command -v numactl >/dev/null ||
      { echo "NUMA=$NUMA given, but numactl is not installed" >&2; exit 1; }
    numactl --hardware | grep -q "^node $NUMA cpus:" ||
      { echo "no NUMA node $NUMA (see numactl --hardware)" >&2; exit 1; }
    NUMA_CMD=(numactl "--cpunodebind=$NUMA" "--membind=$NUMA")
    NUMA_LABEL="numa$NUMA"
    echo "(pinned to NUMA node $NUMA: its CPUs and its memory only)" >&2
  else
    echo "(not pinned to a NUMA node; NUMA=N pins the run to node N)" >&2
  fi
  # For the results: also on stdout, next to the program's key=value lines.
  echo "numa=$NUMA_LABEL"
}
