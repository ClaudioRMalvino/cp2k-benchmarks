#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/toolchain_setup_utils.sh"

SCRATCH_USER="crm98"
SCRATCH_ROOT="${SCRATCH_ROOT:-/local/data/public}"
SCRATCH_REPO="${SCRATCH_ROOT}/${SCRATCH_USER}/original_cp2k"
JOBS="${JOBS:-32}"

export TMPDIR="${SCRATCH_ROOT}/${SCRATCH_USER}/tmp"
export PIP_CACHE_DIR="${SCRATCH_ROOT}/${SCRATCH_USER}/pip-cache"

mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

echo "SCRATCH_USER:  $SCRATCH_USER"
echo "SCRATCH_ROOT:  $SCRATCH_ROOT"
echo "SCRATCH_REPO:  $SCRATCH_REPO"
echo "JOBS:          $JOBS"
echo

if [ ! -f "$SCRATCH_REPO/tools/toolchain/install/setup" ]; then
  echo "Missing scratch toolchain setup."
  echo "Run ./install_cp2k.sh first"
  exit 1
fi

cd "$SCRATCH_REPO"

cp2k_source_toolchain_setup "$SCRATCH_REPO"

export PKG_CONFIG_PATH="$SCRATCH_REPO/tools/toolchain/install/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CMAKE_PREFIX_PATH="$SCRATCH_REPO/tools/toolchain/install:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$SCRATCH_REPO/tools/toolchain/install/lib:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$SCRATCH_REPO/tools/toolchain/install/lib:${LIBRARY_PATH:-}"
export CPATH="$SCRATCH_REPO/tools/toolchain/install/include:${CPATH:-}"

rm -rf build install

CMAKE_OPTS="${CP2K_CMAKE_OPTIONS:-}"
if [ -z "$CMAKE_OPTS" ]; then
  CMAKE_OPTS="-DCP2K_DATA_DIR=$SCRATCH_REPO/data -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF -DCP2K_USE_DFTD4=OFF"
fi

cmake -S . -B build \
  -DCMAKE_INSTALL_PREFIX=./install \
  $CMAKE_OPTS \
  -DCMAKE_BUILD_TYPE=Release \
  "-DCMAKE_Fortran_FLAGS=-march=native -funroll-loops -ftree-vectorize" \
  "-DCMAKE_C_FLAGS=-march=native" \
  "-DCMAKE_CXX_FLAGS=-march=native"
cmake --build build --target install -j "$JOBS"

echo
echo "Full CP2K install build complete."
echo "Installed under: $SCRATCH_REPO/install"
find "$SCRATCH_REPO/install" -type f -name 'cp2k*' 2>/dev/null || true
