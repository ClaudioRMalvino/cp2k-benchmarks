#!/bin/bash
set -e

# Build feature/nnp-mace at Dhruv's exact AVX-2 profile (-O3 -xCORE-AVX2 -fp-model=precise -qopenmp,
# sequential BLAS) so the head-to-head vs pr/nnp-cpu-cell-list differs only in the NNP source.
# Stages cp2k.psmp -> cp2k_binaries/csd3/feature-nnp-mace/; MACE stays ON but is idle on pure-NNP
# runs, so NNP-path timing equals nnp-chebyshev. Otherwise identical to cp2k_CSD3_opt_build.sh.

source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

CP2K_USE_MACE=${CP2K_USE_MACE:-ON}
SYMMETRIX_SRC=${SYMMETRIX_SRC:-$HOME/symmetrix/libsymmetrix}
SYMMETRIX_PREFIX=${SYMMETRIX_PREFIX:-$HOME/symmetrix/cp2k_prefix}
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

AVX2_FFLAGS="-O3 -xCORE-AVX2 -fp-model=precise -qopenmp"
AVX2_CFLAGS="-O3 -xCORE-AVX2 -fp-model=precise -qopenmp"

cd ~/cp2k_optimized
[[ "${CLEAN_BUILD:-1}" == "1" ]] && rm -rf build || echo "incremental: keeping existing build/"

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
    -DCP2K_USE_MACE=${CP2K_USE_MACE} \
    -DCP2K_BLAS_THREADING=sequential \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_C_COMPILER=mpiicc \
    -DCMAKE_CXX_COMPILER=mpiicpc \
    -DCMAKE_CXX_LINK_EXECUTABLE="mpiifort <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES> -nofor-main -cxxlib -L${GCC11_LIB} -lstdc++" \
    -DCMAKE_INSTALL_PREFIX=/home/crm98/cp2k_optimized/install \
    -DCP2K_DATA_DIR=/home/crm98/cp2k_optimized/data \
    -DCMAKE_C_FLAGS="$AVX2_CFLAGS" \
    -DCMAKE_CXX_FLAGS="$AVX2_CFLAGS" \
    -DCMAKE_Fortran_FLAGS="$AVX2_FFLAGS"

cmake --build build -j 16
cmake --install build

# --- stage into the benchmark binary root (parity with dhruv-cell-list) -------
STAGE=/rds/user/$USER/hpc-work/cp2k_binaries/csd3/feature-nnp-mace
echo "Staging binary -> $STAGE"
mkdir -p "$STAGE"
cp -a /home/crm98/cp2k_optimized/install/bin/cp2k.psmp "$STAGE/cp2k.psmp"
rm -rf "$STAGE/lib"; cp -a /home/crm98/cp2k_optimized/install/lib "$STAGE/lib"
echo "Staged: $(ls -la "$STAGE/cp2k.psmp")"

# --- NNP regtest --------------------------------------------------------------
if [[ "${SKIP_REGTEST:-0}" != "1" ]]; then
   WORK_BASE_DIR="/rds/user/$USER/hpc-work/regtesting"
   INSTALL_DIR="/home/$USER/cp2k_optimized/install/bin"
   mkdir -p "$WORK_BASE_DIR"
   cd /home/$USER/cp2k_optimized/tests
   echo "Starting NNP regression tests (nnp-mace AVX-2)..."
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

echo "=== nnp-mace AVX-2 build done. Binary: $STAGE/cp2k.psmp ==="
