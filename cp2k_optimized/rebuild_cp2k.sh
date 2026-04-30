#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/toolchain_setup_utils.sh"

SCRATCH_USER="crm98"
SCRATCH_ROOT="${SCRATCH_ROOT:-/local/data/public}"
HOME_REPO="/home/raid/crm98/cp2k"
SCRATCH_REPO="/local/data/public/crm98/original_cp2k"
BUILD_DIR="${BUILD_DIR:-$SCRATCH_REPO/build}"
INSTALL_DIR="${INSTALL_DIR:-$SCRATCH_REPO/install}"
JOBS="${JOBS:-16}"
CP2K_FORCE_CONFIGURE="${CP2K_FORCE_CONFIGURE:-0}"

export TMPDIR="${SCRATCH_ROOT}/${SCRATCH_USER}/tmp"
export PIP_CACHE_DIR="${SCRATCH_ROOT}/${SCRATCH_USER}/pip-cache"

mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

echo "SCRATCH_USER:  $SCRATCH_USER"
echo "SCRATCH_ROOT:  $SCRATCH_ROOT"
echo "SCRATCH_REPO:  $SCRATCH_REPO"
echo "HOME_REPO:     $HOME_REPO"
echo "JOBS:          $JOBS"
echo

if [ ! -d "$SCRATCH_REPO" ]; then
  echo "Scratch repo missing: $SCRATCH_REPO"
  echo "Run ./install_cp2k.sh and ./build_cp2k.sh first."
  exit 1
fi

if [ ! -f "$SCRATCH_REPO/tools/toolchain/install/setup" ]; then
  echo "Scratch toolchain missing."
  echo "Run ./install_cp2k.sh first."
  exit 1
fi

# --- Incremental rsync (checksum-based to avoid mtime issues) ---
SYNC_LOG="$(mktemp)"
trap 'rm -f "$SYNC_LOG"' EXIT

echo "Syncing source from home to scratch..."
rsync -aci --checksum --delete \
  --exclude '.git' \
  --exclude '*.mod' \
  --exclude '*.smod' \
  --exclude 'tools/toolchain/install' \
  --exclude 'tools/toolchain/build' \
  --exclude 'build' \
  --exclude 'install' \
  "$HOME_REPO/" "$SCRATCH_REPO/" | tee "$SYNC_LOG"

cd "$SCRATCH_REPO"

# --- Source the toolchain environment ---
cp2k_source_toolchain_setup "$SCRATCH_REPO"

export PKG_CONFIG_PATH="$SCRATCH_REPO/tools/toolchain/install/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CMAKE_PREFIX_PATH="$SCRATCH_REPO/tools/toolchain/install:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$SCRATCH_REPO/tools/toolchain/install/lib:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$SCRATCH_REPO/tools/toolchain/install/lib:${LIBRARY_PATH:-}"
export CPATH="$SCRATCH_REPO/tools/toolchain/install/include:${CPATH:-}"

CMAKE_OPTS="${CP2K_CMAKE_OPTIONS:-}"
if [ -z "$CMAKE_OPTS" ]; then
  CMAKE_OPTS="-DCP2K_DATA_DIR=$SCRATCH_REPO/data -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_LIBXSMM=OFF -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF -DCP2K_USE_DFTD4=OFF"
fi

# --- Decide whether to reconfigure ---
NEED_CONFIGURE=0

if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
  NEED_CONFIGURE=1
  echo "No CMakeCache.txt found, will configure."
fi

if [ "$CP2K_FORCE_CONFIGURE" = "1" ]; then
  NEED_CONFIGURE=1
  echo "CP2K_FORCE_CONFIGURE=1, will reconfigure."
fi

if grep -Eq '(^>f|\.[^ ]*t).* (CMakeLists\.txt|.*\.cmake)$' "$SYNC_LOG" 2>/dev/null; then
  NEED_CONFIGURE=1
  echo "CMake input changed during rsync, will reconfigure."
fi

mkdir -p "$BUILD_DIR" "$INSTALL_DIR"

if [ "$NEED_CONFIGURE" -eq 1 ]; then
  echo "Configuring build tree..."
  eval cmake -S . -B "\"$BUILD_DIR\"" -DCMAKE_INSTALL_PREFIX="\"$INSTALL_DIR\"" $CMAKE_OPTS
else
  echo "Skipping configure, no CMake input changes detected."
fi

echo "Starting incremental build..."
cmake --build "$BUILD_DIR" --target install -j "$JOBS"

echo
echo "CP2K rebuild complete."
echo "Installed under: $INSTALL_DIR"
echo

# --- Run NNP Regression Tests ---
echo "Starting NNP regression tests..."
WORK_BASE_DIR="${SCRATCH_ROOT}/${SCRATCH_USER}/regtesting"
mkdir -p "$WORK_BASE_DIR"

cd "$SCRATCH_REPO/tests"

./do_regtest.py \
  --mpiranks 8 \
  --ompthreads 2 \
  --restrictdir NNP/regtest-1 \
  --workbasedir "$WORK_BASE_DIR" \
  "$INSTALL_DIR/bin" \
  psmp

echo "Regression tests complete."

source /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup

echo "Shared libraries sourced from source /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup"