#!/usr/bin/env bash
# Records the experiment machine's specifications -- what a paper's setup
# section needs and what affects these experiments' timings -- into
# results/machine_specs.txt (config.sh). Run it on the experiment machine,
# e.g. right before run.sh. Nothing needs root; a few details (memory type
# and speed) are only readable with it and are then skipped.
#
#   ./machine_specs.sh

set -uo pipefail
source "$(dirname "$0")/config.sh"

OUT=$RESULTS/machine_specs.txt
mkdir -p "$RESULTS"

section() { echo; echo "== $*"; }
try() { "$@" 2>/dev/null || echo "(not available: $*)"; }

{
  echo "Machine specifications, $(date), host $(hostname)"

  section "CPU"
  try lscpu | grep -E "^(Model name|Architecture|Socket\(s\)|Core\(s\) per socket|Thread\(s\) per core|CPU\(s\)|CPU max MHz|CPU min MHz|NUMA node\(s\)|L1d|L1i|L2|L3|Frequency boost)"
  echo "hardware threads sharing a core with cpu0: $(cat /sys/devices/system/cpu/cpu0/topology/thread_siblings_list 2>/dev/null)"
  # How the L3 is split (which CPUs share one L3 slice): the level-3 cache of cpu0.
  for index in /sys/devices/system/cpu/cpu0/cache/index*; do
    [ "$(cat "$index/level" 2>/dev/null)" = 3 ] &&
      echo "L3 of cpu0: $(cat "$index/size"), shared by CPUs $(cat "$index/shared_cpu_list")"
  done

  section "Clock and power settings (affect timings)"
  echo "scaling governor (cpu0): $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
  echo "scaling driver (cpu0): $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_driver 2>/dev/null || echo unknown)"
  echo "boost: $(cat /sys/devices/system/cpu/cpufreq/boost 2>/dev/null || echo unknown) (1 = on)"
  echo "base frequency (cpu0): $(cat /sys/devices/system/cpu/cpu0/cpufreq/base_frequency 2>/dev/null || echo unknown) kHz"
  echo "energy/performance preference (cpu0): $(cat /sys/devices/system/cpu/cpu0/cpufreq/energy_performance_preference 2>/dev/null || echo unknown)"
  echo "frequency range (cpu0): $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_min_freq 2>/dev/null)-$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq 2>/dev/null) kHz"

  section "Memory"
  grep -E "^(MemTotal|MemAvailable|HugePages_Total|Hugepagesize)" /proc/meminfo
  echo "transparent huge pages: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null)"
  echo "THP defrag: $(cat /sys/kernel/mm/transparent_hugepage/defrag 2>/dev/null)"
  if command -v dmidecode > /dev/null && dmidecode -t memory > /dev/null 2>&1; then
    dmidecode -t memory | grep -E "^\s+(Type|Speed|Configured Memory Speed|Size):" | sort | uniq -c
  else
    echo "memory type/speed: needs root (sudo dmidecode -t memory)"
  fi

  section "NUMA"
  try numactl --hardware

  section "Operating system"
  try grep PRETTY_NAME /etc/os-release
  echo "kernel: $(uname -r)"

  section "Compilers and build"
  echo "g++: $(g++ --version | head -1)"
  echo "cmake: $(cmake --version 2>/dev/null | head -1)"
  echo "our programs: g++ -std=c++2a -O3 -march=native -pthread -DNDEBUG (+ -DWITH_POWERED -DSMALL_BLOCK_SIZE=256 -DLARGE_BLOCK_SIZE=16384; varki part -std=c++17 -O3 -march=native)"
  echo "sdsl-lite: $(cat "$REPO/RLZ-Varki/thirdparty/sdsl-lite/VERSION" 2>/dev/null || echo '?'), bundled in RLZ-Varki/thirdparty (repository commit below); built on $(cat "$REPO/RLZ-Varki/build/sdsl/.built-on" 2>/dev/null || echo 'unknown')"

  section "Code"
  echo "repository: $(git -C "$REPO" log -1 --format='%H (%cd) %s' 2>/dev/null)"
  echo "branch: $(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null)"
  echo "uncommitted changes: $(git -C "$REPO" status --porcelain 2>/dev/null | wc -l) files"

  section "Experiment settings (config.sh)"
  echo "ROOT=$ROOT"
  echo "CPUS=${CPUS:-} NUMA=${NUMA:-} THREADS_SWEEP=$THREADS_SWEEP"
  echo "VARIANTS=$VARIANTS"

  section "Load at the time of recording"
  uptime
} > "$OUT"

cat "$OUT"
echo
echo "written to $OUT"
