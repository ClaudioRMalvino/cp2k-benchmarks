#!/bin/bash
set -e

source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

# --- MACE (symmetrix) backend ----------------------------------------------
# Rebuild reuses the symmetrix prefix staged by the full build; re-stage it
# here if it is missing so an incremental rebuild (or a triggered reconfigure)
# can still find the headers/libs. Ignored on branches without MACE support.
CP2K_USE_MACE=${CP2K_USE_MACE:-ON}
SYMMETRIX_SRC=${SYMMETRIX_SRC:-$HOME/symmetrix/libsymmetrix}
SYMMETRIX_PREFIX=${SYMMETRIX_PREFIX:-$HOME/symmetrix/cp2k_prefix}
if [[ "$CP2K_USE_MACE" == "ON" ]]; then
  if [[ ! -d "$SYMMETRIX_PREFIX/include" && -f "$SYMMETRIX_SRC/build/libsymmetrix.a" ]]; then
    mkdir -p "$SYMMETRIX_PREFIX/lib" "$SYMMETRIX_PREFIX/include"
    cp -a "$SYMMETRIX_SRC/build/libsymmetrix.a" \
          "$SYMMETRIX_SRC/build/external/sphericart/sphericart/libsphericart.a" \
          "$SYMMETRIX_PREFIX/lib/"
    cp -a "$SYMMETRIX_SRC"/source/*.hpp "$SYMMETRIX_PREFIX/include/"
  fi
  [[ -d "$SYMMETRIX_PREFIX" ]] && export SYMMETRIX_ROOT="$SYMMETRIX_PREFIX"
fi

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
