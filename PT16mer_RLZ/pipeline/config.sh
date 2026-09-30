# The experiment's paths and settings, in one place. Sourced by every
# pipeline script; any value can be overridden from the environment, e.g.
#   ROOT=/scratch/rlz ./clean.sh
#
# Layout under ROOT:
#
#   raw/<dataset>/inputs/        the sampled input files, as they come (FASTA or
#                                plain sequence, optionally .gz)
#   raw/<dataset>/refs/          the reference files (1-3 per dataset), likewise
#
#   data/<dataset>/inputs/       the cleaned inputs: plain A/C/G/T, one file per
#                                raw input, named <name>.txt (clean.sh)
#   data/<dataset>/refs/<ref>/   one folder per reference, with everything
#                                computed for it (prepare.sh): the cleaned
#                                reference <ref>.txt, its suffix array
#                                <ref>.txt.sa, the reversed reference and its
#                                powered index <ref>.txt.rev_four.bwt (+ _data),
#                                PFP-eBWT's intermediate files. Whatever exists
#                                is kept; the PT16 tables are rebuilt by every
#                                run (fast, and a sanity check).
#   data/<dataset>/cleaning/     one statistics CSV per cleaned file
#
#   results/<dataset>/<ref>/     the runs' CSVs and logs (run.sh)
#   tools/                       PFP-eBWT and the small tools (setup.sh)

# This machine's own settings, if any: pipeline/local.sh (not in git), e.g.
#   ROOT=/abga/work/anastasd/DCC/rlz_experiments
#   CPUS=0-31
# Variables already set in the environment take precedence over it.
_local_config="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/local.sh"
if [ -f "$_local_config" ]; then
  while IFS='=' read -r _key _value; do
    case "$_key" in ''|\#*) continue ;; esac
    [ -z "${!_key:-}" ] && eval "$_key=$_value"
  done < "$_local_config"
fi
unset _local_config _key _value

ROOT=${ROOT:-$HOME/rlz_experiments}
RAW=${RAW:-$ROOT/raw}
DATA=${DATA:-$ROOT/data}
RESULTS=${RESULTS:-$ROOT/results}
TOOLS=${TOOLS:-$ROOT/tools}

# PFP-eBWT's build directory (with pfpebwt), built by setup.sh.
PFP=${PFP:-$TOOLS/PFP-eBWT/build}

# The datasets (folder names under raw/ and data/).
DATASETS=${DATASETS:-"ecoli chr19 yeast"}

# The parsers that are compared (all give the same greedy left-to-right parse).
VARIANTS=${VARIANTS:-"sa-binary-search lrf-ms pt16 pt16-v2 sassy varki powered-fwd-escape powered-pt16-fwd-escape"}

# verify.sh: check only the first N input files of each dataset (empty: all).
VERIFY_FILES=${VERIFY_FILES:-}

# Thread counts of the parallel sweep (on each dataset's first reference).
THREADS_SWEEP=${THREADS_SWEEP:-"1,2,4,8,16,32"}

# NUMA node to pin every run to (empty: not pinned; not needed on a machine
# with one NUMA node).
NUMA=${NUMA:-}

# CPUs every run may use (taskset -c; empty: all). On a machine with two
# hardware threads per core, one per core keeps threads from sharing a core,
# e.g. CPUS=0-31 when cpu0's thread_siblings_list is "0,32"
# (/sys/devices/system/cpu/cpu0/topology/thread_siblings_list).
CPUS=${CPUS:-}

# The repository (this file is PT16mer_RLZ/pipeline/config.sh).
PIPELINE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PT16=$(cd "$PIPELINE/.." && pwd)
REPO=$(cd "$PT16/.." && pwd)

# The references of a dataset: the folders under data/<dataset>/refs, sorted
# (the first one is the thread sweep's reference).
dataset_refs() {
  local ds=$1
  [ -d "$DATA/$ds/refs" ] || return 0
  find "$DATA/$ds/refs" -mindepth 1 -maxdepth 1 -type d | sort
}

# The name of a raw file without .gz and a FASTA/text extension.
file_stem() {
  local name
  name=$(basename "$1")
  name=${name%.gz}
  for ext in .fasta .fa .fna .fas .txt .seq; do
    name=${name%"$ext"}
  done
  echo "$name"
}
