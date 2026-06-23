#!/bin/bash
set -e

source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

# --- MACE (symmetrix) backend ----------------------------------------------
# feature/nnp-mace adds an optional MACE manybody potential via the symmetrix
# library. Default ON; run with CP2K_USE_MACE=OFF to build without it. On a
# branch that has no MACE support (e.g. feature/nnp-chebyshev) the CMake flag
# is simply unused (CMake prints a harmless "unused variable" warning), so the
# block below is safe to keep across branches.
# FindSymmetrix expects a lib/+include/ prefix and resolves it via the
# SYMMETRIX_ROOT environment variable, so we stage one from the symmetrix
# build tree (~/symmetrix/libsymmetrix) here.
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

cd ~/cp2k_optimized
rm -rf build

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
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_C_COMPILER=mpiicc \
    -DCMAKE_CXX_COMPILER=mpiicpc \
    -DCMAKE_CXX_LINK_EXECUTABLE="mpiifort <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES> -nofor-main -cxxlib -L${GCC11_LIB} -lstdc++" \
    -DCMAKE_INSTALL_PREFIX=/home/crm98/cp2k_optimized/install \
    -DCP2K_DATA_DIR=/home/crm98/cp2k_optimized/data \
    -DCMAKE_C_FLAGS="-O2 -g -xCORE-AVX512 -qopenmp" \
    -DCMAKE_CXX_FLAGS="-O2 -g -xCORE-AVX512 -qopenmp" \
    -DCMAKE_Fortran_FLAGS="-O2 -g -xCORE-AVX512 -qopenmp -funroll-loops -ftree-vectorize"

cmake --build build -j 16
cmake --install build

if [[ "${SKIP_REGTEST:-0}" != "1" ]]; then
   WORK_BASE_DIR="/rds/user/$USER/hpc-work/regtesting"
   INSTALL_DIR="/home/$USER/cp2k_optimized/install/bin"
   mkdir -p "$WORK_BASE_DIR"
   cd /home/$USER/cp2k_optimized/tests

   echo "Starting NNP regression tests..."
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
