#!/usr/bin/env bash
# Builds sdsl-lite and divsufsort for Varki's RLZ into RLZ-Varki/build/sdsl:
# headers in include/, libsdsl.a and libdivsufsort{,64}.a in lib/. Once per
# machine; skipped if the libraries exist.
#
#   ./build_sdsl.sh
#
# The same steps as sdsl-lite's install.sh (which RLZ-Varki's CMake calls):
# cmake with the install prefix, make sdsl, make install. install.sh itself
# cannot run in this repository: it builds in sdsl-lite's own build/
# directory (with its clean.sh), which is missing here -- RLZ-Varki's
# .gitignore ignores build/, so it was never committed. So this script
# configures sdsl-lite in RLZ-Varki/build/sdsl-build instead.
#
# After this, run_rlz_suite.sh / run_rlz_parallel.sh build their programs
# with varki (build_if_needed in rlz_datasets.sh adds -DWITH_VARKI and the
# sdsl flags when RLZ-Varki/build/sdsl/lib/libsdsl.a exists).

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
VARKI=$(cd "$HERE/../../RLZ-Varki" && pwd)
SOURCE=$VARKI/thirdparty/sdsl-lite
PREFIX=$VARKI/build/sdsl
BUILD=$VARKI/build/sdsl-build

# sdsl-lite compiles with -march=native: a library built on another machine
# (copied along with the repository) may not suit this CPU, so the machine is
# recorded next to it and a library from elsewhere is rebuilt.
STAMP=$PREFIX/.built-on
MACHINE="$(uname -sm) $(hostname)"
if [ -f "$PREFIX/lib/libsdsl.a" ] && [ -f "$PREFIX/lib/libdivsufsort.a" ] &&
   [ "$(cat "$STAMP" 2>/dev/null)" = "$MACHINE" ]; then
  echo "sdsl already built in $PREFIX for this machine (delete it to rebuild)"
  exit 0
fi
if [ -e "$PREFIX" ]; then
  echo "sdsl in $PREFIX was not built on this machine: rebuilding"
  rm -rf "$PREFIX" "$BUILD"
fi

mkdir -p "$PREFIX" "$BUILD"
cd "$BUILD"
# CMAKE_POLICY_VERSION_MINIMUM: sdsl-lite asks for an old CMake version,
# which CMake 4 refuses without it (older CMake ignores it).
# Position-independent code: the programs are linked as PIE (the default of
# e.g. Ubuntu's GCC), which a non-PIC static libsdsl.a cannot be part of.
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 "$SOURCE"
make -j"$(nproc 2>/dev/null || echo 4)" sdsl
make install
echo "$MACHINE" > "$STAMP"
ls "$PREFIX/lib"
