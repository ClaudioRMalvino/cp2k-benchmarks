#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/toolchain_setup_utils.sh"

SCRATCH_USER="crm98"
SCRATCH_ROOT="${SCRATCH_ROOT:-/local/data/public}"
SCRATCH_REPO="${SCRATCH_ROOT}/${SCRATCH_USER}/cp2k-buildtree"
SRC_HOME="/home/raid/crm98/cp2k-upstream-master"
# Reuse the feature-branch toolchain — avoids a separate multi-hour toolchain install
# and guarantees DBCSR and all other libraries are present on this node.
TOOLCHAIN_REPO="${SCRATCH_ROOT}/${SCRATCH_USER}/original_cp2k"
JOBS="${JOBS:-32}"

export TMPDIR="${SCRATCH_ROOT}/${SCRATCH_USER}/tmp"
export PIP_CACHE_DIR="${SCRATCH_ROOT}/${SCRATCH_USER}/pip-cache"

mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

echo "SCRATCH_USER:  $SCRATCH_USER"
echo "SCRATCH_ROOT:  $SCRATCH_ROOT"
echo "SCRATCH_REPO:  $SCRATCH_REPO"
echo "TOOLCHAIN:     $TOOLCHAIN_REPO"
echo "JOBS:          $JOBS"
echo

# Clone source if scratch doesn't exist yet (toolchain install is separate)
if [ ! -d "$SCRATCH_REPO" ]; then
  echo "Cloning upstream-master source to scratch..."
  git clone "$SRC_HOME" "$SCRATCH_REPO"
fi

if [ ! -f "$TOOLCHAIN_REPO/tools/toolchain/install/setup" ]; then
  echo "Shared toolchain missing at $TOOLCHAIN_REPO"
  echo "Run the feature branch's ./install_cp2k.sh first"
  exit 1
fi

cd "$SCRATCH_REPO"

cp2k_source_toolchain_setup "$TOOLCHAIN_REPO"

export PKG_CONFIG_PATH="$TOOLCHAIN_REPO/tools/toolchain/install/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CMAKE_PREFIX_PATH="$TOOLCHAIN_REPO/tools/toolchain/install:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$TOOLCHAIN_REPO/tools/toolchain/install/lib:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$TOOLCHAIN_REPO/tools/toolchain/install/lib:${LIBRARY_PATH:-}"
export CPATH="$TOOLCHAIN_REPO/tools/toolchain/install/include:${CPATH:-}"

rm -rf build install

CMAKE_OPTS="${CP2K_CMAKE_OPTIONS:-}"
if [ -z "$CMAKE_OPTS" ]; then
  CMAKE_OPTS="-DCP2K_DATA_DIR=$SCRATCH_REPO/data -DCP2K_USE_MPI=ON -DCP2K_USE_DBCSR_CONFIG=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF -DCP2K_USE_DFTD4=OFF -DCP2K_USE_ELPA=OFF -DCP2K_USE_SIRIUS=OFF"
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
