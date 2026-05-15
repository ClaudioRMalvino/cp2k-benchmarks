#!/bin/bash
set -e

source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

cd ~/cp2k_optimized

cmake --build build -j 16
cmake --install build

WORK_BASE_DIR="/rds/user/$USER/hpc-work/cp2k_optimized_regtesting"
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
