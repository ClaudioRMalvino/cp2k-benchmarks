#!/bin/bash
set -e

source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

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
   echo "SKIP_REGTEST=1 — skipping NNP regression tests."
fi
