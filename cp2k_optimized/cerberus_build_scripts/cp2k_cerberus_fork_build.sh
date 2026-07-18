#!/usr/bin/env bash
set -euo pipefail

# Build the optimized CP2K fork (~/cp2k, DhruvSkyy/cp2k) on cerberus.
#
# Merges:
#   - cerberus machinery (cp2k_cerberus_opt_install.sh / _build.sh / _rebuild.sh):
#     rsync home -> scratch, toolchain bootstrap, incremental configure/build,
#     DBCSR config workaround, GCC + system OpenMPI/OpenBLAS.
#   - CSD3 feature flags (cp2k_CSD3_opt_build.sh): the CP2K_USE_* ON/OFF set,
#     the optional MACE (symmetrix) backend block, and SKIP_REGTEST gating.
#
# Intel-specific CSD3 bits (mpiifort/mpiicc, -xCORE-AVX512, the CXX link hack)
# do not apply on cerberus: only GCC 13 is available, and -march=native on
# this node already enables AVX512.
#
# Home quota is 5 GB, so everything heavy (source copy, toolchain, build,
# install, regtest work dirs, tmp, pip cache) lives under /data/cerberus1.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/toolchain_setup_utils.sh"

SCRATCH_USER="crm98"
SCRATCH_ROOT="${SCRATCH_ROOT:-/data/cerberus1}"
HOME_REPO="${HOME_REPO:-/home/raid/crm98/cp2k}"
SCRATCH_REPO="${SCRATCH_REPO:-${SCRATCH_ROOT}/${SCRATCH_USER}/cp2k_optimized}"
BUILD_DIR="${BUILD_DIR:-$SCRATCH_REPO/build}"
INSTALL_DIR="${INSTALL_DIR:-$SCRATCH_REPO/install}"
JOBS="${JOBS:-32}"
CP2K_FORCE_CONFIGURE="${CP2K_FORCE_CONFIGURE:-0}"

export TMPDIR="${SCRATCH_ROOT}/${SCRATCH_USER}/tmp"
export PIP_CACHE_DIR="${SCRATCH_ROOT}/${SCRATCH_USER}/pip-cache"
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

# The login environment has miniforge/conda active (/lsc/opt/miniforge3),
# which ships stale GCC-12.4 builds of mctc-lib/s-dftd3 that CMake picks up
# over the packages' own bundled subprojects (this broke the tblite toolchain
# stage). Strip conda from the environment for the whole build.
if command -v conda > /dev/null 2>&1 || [ -n "${CONDA_PREFIX:-}" ]; then
  echo "Stripping conda/miniforge from build environment."
  PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -v 'miniforge3\|condabin' | paste -sd:)"
  export PATH
  unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL CONDA_EXE CONDA_PYTHON_EXE PYTHONHOME || true
fi

echo "SCRATCH_USER:  $SCRATCH_USER"
echo "SCRATCH_ROOT:  $SCRATCH_ROOT"
echo "HOME_REPO:     $HOME_REPO"
echo "SCRATCH_REPO:  $SCRATCH_REPO"
echo "JOBS:          $JOBS"
echo

# --- Sync source: home fork -> scratch build tree ---------------------------
if [ ! -d "$SCRATCH_REPO" ]; then
  echo "Creating scratch copy of $HOME_REPO ..."
  mkdir -p "$(dirname "$SCRATCH_REPO")"
fi

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

# --- Toolchain: bootstrap once, reuse afterwards ----------------------------
# A crashed toolchain run leaves a partial install/setup behind, so completion
# is tracked with a sentinel written only after the installer succeeds.
# Re-running resumes past already-built packages.
# tblite/dftd4 are excluded: the CP2K build disables both (CP2K_USE_TBLITE=OFF,
# CP2K_USE_DFTD4=OFF), and tblite is what broke against the conda mctc-lib.
TOOLCHAIN_DONE="$SCRATCH_REPO/tools/toolchain/install/.cerberus_toolchain_complete"
if [ ! -f "$TOOLCHAIN_DONE" ]; then
  echo "Toolchain not complete - (re)running toolchain install (this takes a while)..."
  cd "$SCRATCH_REPO/tools/toolchain"
  ./install_cp2k_toolchain.sh \
    --install-all \
    --with-gcc=system \
    --with-openmpi=system \
    --with-openblas=system \
    --math-mode=openblas \
    --with-elpa=no \
    --with-sirius=no \
    --with-tblite=no \
    --with-dftd4=no
  touch "$TOOLCHAIN_DONE"
fi

cd "$SCRATCH_REPO"

cp2k_source_toolchain_setup "$SCRATCH_REPO"

export PKG_CONFIG_PATH="$SCRATCH_REPO/tools/toolchain/install/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CMAKE_PREFIX_PATH="$SCRATCH_REPO/tools/toolchain/install:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$SCRATCH_REPO/tools/toolchain/install/lib:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$SCRATCH_REPO/tools/toolchain/install/lib:${LIBRARY_PATH:-}"
export CPATH="$SCRATCH_REPO/tools/toolchain/install/include:${CPATH:-}"

# CMake 4.x case-sensitive search misses the lowercase 'dbcsr' directory, so locate DBCSRConfig.cmake explicitly.
DBCSR_CMAKE_DIR=""
for d in "$SCRATCH_REPO/tools/toolchain/install"/dbcsr-*/lib/cmake/dbcsr; do
  [ -f "$d/DBCSRConfig.cmake" ] && DBCSR_CMAKE_DIR="$d" && break
done
if [ -z "$DBCSR_CMAKE_DIR" ]; then
  echo "ERROR: DBCSRConfig.cmake not found under $SCRATCH_REPO/tools/toolchain/install"
  exit 1
fi
echo "DBCSR cmake:   $DBCSR_CMAKE_DIR"

# --- MACE (symmetrix) backend ------------------------------------------------
# Carried over from the CSD3 script. Default OFF on cerberus: the current fork
# branch (master) has no MACE support in its CMake, and no symmetrix build
# exists on this machine yet. Once symmetrix is built here and a MACE-enabled
# branch is checked out, run with CP2K_USE_MACE=ON. On a branch without MACE
# support the flag is unused (harmless CMake warning). The staged prefix goes
# to scratch, not home, because of the 5 GB home quota.
CP2K_USE_MACE=${CP2K_USE_MACE:-OFF}
SYMMETRIX_SRC=${SYMMETRIX_SRC:-$HOME/symmetrix/libsymmetrix}
SYMMETRIX_PREFIX=${SYMMETRIX_PREFIX:-${SCRATCH_ROOT}/${SCRATCH_USER}/symmetrix_cp2k_prefix}
if [[ "$CP2K_USE_MACE" == "ON" ]]; then
  if [[ ! -f "$SYMMETRIX_SRC/build/libsymmetrix.a" ]]; then
    echo "MACE: $SYMMETRIX_SRC/build/libsymmetrix.a not found - build symmetrix first" >&2
    exit 1
  fi
  echo "MACE: staging symmetrix prefix -> $SYMMETRIX_PREFIX"
  rm -rf "$SYMMETRIX_PREFIX"
  mkdir -p "$SYMMETRIX_PREFIX/lib" "$SYMMETRIX_PREFIX/include"
  cp -a "$SYMMETRIX_SRC/build/libsymmetrix.a" \
        "$SYMMETRIX_SRC/build/external/sphericart/sphericart/libsphericart.a" \
        "$SYMMETRIX_PREFIX/lib/"
  cp -a "$SYMMETRIX_SRC"/source/*.hpp "$SYMMETRIX_PREFIX/include/"
  export SYMMETRIX_ROOT="$SYMMETRIX_PREFIX"
fi

# --- CMake options: CSD3 feature-flag set + cerberus-specific flags ---------
CMAKE_OPTS="${CP2K_CMAKE_OPTIONS:-}"
if [ -z "$CMAKE_OPTS" ]; then
  CMAKE_OPTS="-DCP2K_DATA_DIR=$SCRATCH_REPO/data \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_ELPA=OFF \
    -DCP2K_USE_DEEPMD=OFF \
    -DCP2K_USE_TBLITE=OFF \
    -DCP2K_USE_MIMIC=OFF \
    -DCP2K_USE_SIRIUS=OFF \
    -DCP2K_USE_SPLA=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_ACE=OFF \
    -DCP2K_USE_MACE=${CP2K_USE_MACE} \
    -DCP2K_USE_DFTD4=OFF \
    -DCP2K_USE_LIBXS=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_DBCSR_CONFIG=ON \
    -DDBCSR_DIR=$DBCSR_CMAKE_DIR \
    -DDBCSR_USE_MPI=ON"
fi
# DBCSR_USE_MPI=ON: the toolchain builds DBCSR 2.9.1 with -DUSE_MPI=ON, but
# 2.9.1's DBCSRConfig.cmake does not export the DBCSR_USE_MPI variable that
# CP2K's CMakeLists checks (added to DBCSR after 2.9.1), so the check would
# spuriously fail without this.
# SPLA=OFF: the toolchain only builds SPLA as a SIRIUS dependency, and SIRIUS
# is disabled here; CP2K only uses SPLA for GPU-offloaded GEMM paths.

# --- Configure only when needed (fresh tree, forced, or CMake inputs changed)
NEED_CONFIGURE=0

if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
  NEED_CONFIGURE=1
  echo "No CMakeCache.txt found, will configure."
elif [ ! -f "$BUILD_DIR/Makefile" ] && [ ! -f "$BUILD_DIR/build.ninja" ]; then
  # A failed configure leaves CMakeCache.txt without generated build files.
  NEED_CONFIGURE=1
  echo "Previous configure did not complete, will reconfigure."
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
  # also remove generated build files so a failed configure can't leave a
  # stale Makefile that masks the incomplete-configure check on the next run
  rm -f "$BUILD_DIR/CMakeCache.txt" "$BUILD_DIR/Makefile" "$BUILD_DIR/build.ninja"
  rm -rf "$BUILD_DIR/CMakeFiles"
  # -march=native is the GCC equivalent of CSD3's -xCORE-AVX512 on this
  # AVX512 node; Fortran keeps CSD3's -funroll-loops -ftree-vectorize.
  cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    $CMAKE_OPTS \
    -DCMAKE_BUILD_TYPE=Release \
    "-DCMAKE_Fortran_FLAGS=-march=native -funroll-loops -ftree-vectorize" \
    "-DCMAKE_C_FLAGS=-march=native" \
    "-DCMAKE_CXX_FLAGS=-march=native"
else
  echo "Skipping configure, no CMake input changes detected."
fi

echo "Starting build..."
cmake --build "$BUILD_DIR" --target install -j "$JOBS"

echo
echo "CP2K fork build complete."
echo "Installed under: $INSTALL_DIR"
echo

# cp2k.psmp links libcp2k.so from the install tree (CSD3 avoided this with
# BUILD_SHARED_LIBS=OFF); the toolchain setup sourced above covers everything
# else. Job scripts need this same export.
export LD_LIBRARY_PATH="$INSTALL_DIR/lib:${LD_LIBRARY_PATH:-}"

# --- NNP regression tests (gate with SKIP_REGTEST=1, as on CSD3) ------------
if [[ "${SKIP_REGTEST:-0}" != "1" ]]; then
  WORK_BASE_DIR="${SCRATCH_ROOT}/${SCRATCH_USER}/regtesting"
  mkdir -p "$WORK_BASE_DIR"

  cd "$SCRATCH_REPO/tests"

  echo "Starting NNP regression tests..."
  python3 ./do_regtest.py \
    --mpiranks 8 \
    --ompthreads 2 \
    --restrictdir NNP/regtest-1 \
    --workbasedir "$WORK_BASE_DIR" \
    "$INSTALL_DIR/bin" \
    psmp

  echo "NNP regression tests complete."
else
  echo "SKIP_REGTEST=1 - skipping NNP regression tests."
fi

echo "Binary: $INSTALL_DIR/bin/cp2k.psmp"
echo
echo "Runtime environment needed in job scripts:"
echo "  source $SCRATCH_REPO/tools/toolchain/install/setup"
echo "  export LD_LIBRARY_PATH=$INSTALL_DIR/lib:\$LD_LIBRARY_PATH"
