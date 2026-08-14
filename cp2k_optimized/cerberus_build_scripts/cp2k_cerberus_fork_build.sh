#!/usr/bin/env bash
set -euo pipefail

# Build the optimized CP2K fork (~/cp2k, DhruvSkyy/cp2k) on cerberus: cerberus machinery
# (rsync home->scratch, toolchain bootstrap, incremental configure; cf. cp2k_cerberus_opt_*.sh)
# + the CSD3 feature-flag set from cp2k_CSD3_opt_build.sh. GCC 13 only — Intel CSD3 bits don't
# apply, -march=native already gives AVX512. Heavy artefacts on /data/cerberus1 (home quota 5 GB).

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

# Strip conda/miniforge: its stale GCC-12.4 mctc-lib/s-dftd3 shadow the bundled subprojects and broke the tblite toolchain stage.
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
# Sentinel written only on installer success (a crash leaves a partial install/setup); re-runs resume.
# tblite/dftd4 excluded: the CP2K build disables both, and tblite broke against the conda mctc-lib.
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

# --- MACE (symmetrix) backend (from the CSD3 script) -------------------------
# Default OFF: fork master has no MACE CMake and no symmetrix build exists here yet.
# Prefix staged to scratch, not home (5 GB quota).
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
# DBCSR_USE_MPI=ON: DBCSR 2.9.1's config file doesn't export the DBCSR_USE_MPI variable
# CP2K checks (added post-2.9.1), so the check would spuriously fail without it.
# SPLA=OFF: toolchain builds SPLA only as a SIRIUS dep (disabled); CP2K uses it only for GPU GEMM.

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
  # remove generated build files too, so a stale Makefile can't mask the incomplete-configure check
  rm -f "$BUILD_DIR/CMakeCache.txt" "$BUILD_DIR/Makefile" "$BUILD_DIR/build.ninja"
  rm -rf "$BUILD_DIR/CMakeFiles"
  # -march=native = GCC equivalent of CSD3's -xCORE-AVX512 on this AVX512 node; Fortran keeps CSD3's extra flags.
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

# cp2k.psmp links libcp2k.so from the install tree (CSD3 used BUILD_SHARED_LIBS=OFF); job scripts need this export too.
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
