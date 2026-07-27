#!/bin/bash
set -e

# Build Dhruv's optimised-HDNNP branch (CP2K PR #5295, pr/nnp-cpu-cell-list) on
# CSD3 Ice Lake. cp2k_dhruv is checked out on branch pr-nnp-cpu-cell-list.
#
# Mirrors cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_build.sh with three
# deliberate changes for an apples-to-apples benchmark vs feature/nnp-mace:
#   1. source tree    -> /home/crm98/cp2k_dhruv_pr (PR-branch worktree)
#   2. Fortran flags   -> Dhruv's thesis App. A.4 benchmark profile:
#                         -O3 -xCORE-AVX2 -fp-model=precise   (AVX-2, NOT AVX-512;
#                         his audit shows AVX-512 throttles on the 8358)
#   3. -DCP2K_USE_MACE=OFF and no symmetrix staging (his branch has no MACE)
#
# Toolchain is REUSED from ~/cp2k_master/tools/toolchain/install (ifort 2022.1.0,
# Intel MPI 2021.6.0, MKL 2022.1.0 = exactly Dhruv's versions), so no toolchain
# rebuild and no home-quota risk.

source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

SRC=/rds/user/crm98/hpc-work/cp2k_dhruv_omp_worktree
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
STAGE="$BIN_ROOT/dhruv-omp-bitexact"

# Dhruv benchmark profile (App. A.4). -qopenmp for the psmp (hybrid MPI+OMP) build.
DHRUV_FFLAGS="-O3 -xCORE-AVX2 -fp-model=precise -qopenmp"
DHRUV_CFLAGS="-O3 -xCORE-AVX2 -fp-model=precise -qopenmp"

cd "$SRC"
# CLEAN_BUILD=0 keeps the existing build/ for an incremental rebuild (e.g. after
# only toggling a CMake option); default 1 = clean from scratch.
[[ "${CLEAN_BUILD:-1}" == "1" ]] && rm -rf build || echo "incremental: keeping existing build/"

# --- DBCSR version-metadata shim ---------------------------------------------
# The reused toolchain's dbcsr-2.9.0 reports a corrupt version string
# ("2026.1.v2026.1v2026.1-dirty", major 2026), so the newer master's
# find_package(DBCSR 2.8/2.9) SameMajorVersion check rejects it. The library
# itself is correct (2.9.0). Stage a copy OUTSIDE cp2k_master with the version
# string corrected to 2.9.0 and prefer it on CMAKE_PREFIX_PATH. Idempotent.
DBCSR_SRC=/home/crm98/cp2k_master/tools/toolchain/install/dbcsr-2.9.0
DBCSR_SHIM=/rds/user/$USER/hpc-work/dbcsr-shim/dbcsr-2.9.0
if [[ ! -f "$DBCSR_SHIM/lib/cmake/dbcsr/DBCSRConfigVersion.cmake" ]]; then
  echo "DBCSR shim: staging corrected-version copy -> $DBCSR_SHIM"
  rm -rf "$DBCSR_SHIM"; mkdir -p "$(dirname "$DBCSR_SHIM")"
  cp -a "$DBCSR_SRC" "$DBCSR_SHIM"
  sed -i 's/2026\.1\.v2026\.1v2026\.1-dirty/2.9.0/g' \
      "$DBCSR_SHIM/lib/cmake/dbcsr/DBCSRConfigVersion.cmake"
fi
export CMAKE_PREFIX_PATH="$DBCSR_SHIM:${CMAKE_PREFIX_PATH:-}"

cmake -S . -B build \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_ELPA=OFF \
    -DCP2K_USE_DEEPMD=OFF \
    -DCP2K_USE_TBLITE=OFF \
    -DCP2K_USE_MIMIC=OFF \
    -DCP2K_USE_SIRIUS=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_ACE=OFF \
    -DCP2K_USE_LIBXS=OFF \
    -DCP2K_USE_GAUXC=OFF \
    -DCP2K_USE_GREENX=OFF \
    -DCP2K_USE_LIBFCI=OFF \
    -DCP2K_USE_LIBSMEAGOL=OFF \
    -DCP2K_USE_HDF5=OFF \
    -DCP2K_USE_TREXIO=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_MACE=OFF \
    -DCP2K_BLAS_VENDOR=auto \
    -DCP2K_BLAS_THREADING=sequential \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_C_COMPILER=mpiicc \
    -DCMAKE_CXX_COMPILER=mpiicpc \
    -DCMAKE_CXX_LINK_EXECUTABLE="mpiifort <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES> -nofor-main -cxxlib -L${GCC11_LIB} -lstdc++" \
    -DCMAKE_INSTALL_PREFIX="$SRC/install" \
    -DCP2K_DATA_DIR="$SRC/data" \
    -DCMAKE_C_FLAGS="$DHRUV_CFLAGS" \
    -DCMAKE_CXX_FLAGS="$DHRUV_CFLAGS" \
    -DCMAKE_Fortran_FLAGS="$DHRUV_FFLAGS"

cmake --build build -j 16
cmake --install build

# --- stage into the benchmark binary root (parity with master/chebyshev/etc.) --
echo "Staging binary -> $STAGE"
mkdir -p "$STAGE"
cp -a "$SRC/install/bin/cp2k.psmp" "$STAGE/cp2k.psmp"
rm -rf "$STAGE/lib"; cp -a "$SRC/install/lib" "$STAGE/lib"
echo "Staged: $(ls -la "$STAGE/cp2k.psmp")"

# --- Dhruv's bit-exact NNP regtest (his thesis A.5) ---------------------------
if [[ "${SKIP_REGTEST:-0}" != "1" ]]; then
   WORK_BASE_DIR="/rds/user/$USER/hpc-work/regtesting"
   INSTALL_DIR="$SRC/install/bin"
   mkdir -p "$WORK_BASE_DIR"
   cd "$SRC/tests"
   echo "Starting NNP regression tests (Dhruv branch)..."
   python3 ./do_regtest.py \
     --mpiranks 8 \
     --ompthreads 2 \
     --restrictdir NNP/regtest-1 \
     --workbasedir "$WORK_BASE_DIR" \
     "$INSTALL_DIR" \
     psmp
   echo "NNP regression tests complete."
else
   echo "SKIP_REGTEST=1 - skipping NNP regression tests."
fi

echo "=== Dhruv build done. Binary: $STAGE/cp2k.psmp ==="
