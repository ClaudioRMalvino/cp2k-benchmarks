#!/usr/bin/env bash
set -euo pipefail

set +u
# --- Load the Toolchain Environment ---
source /home/raid/crm98/cp2k_binaries/phy-cerberus/setup
set -u

# --- Pass the branch as an argument (defaults to master) ---
# Accepted values: master, feature-nnp-verlet-cells, feature-nnp-native-spline
TARGET_BRANCH=${1:-master}
TIMESTAMP=$(date +%d-%m_%H-%M)

# Per-branch binary cache populated by benchmark_slurm.sh.  Each branch has
# its own cp2k.psmp + lib/, so LD_LIBRARY_PATH cannot accidentally pick up
# another branch's libcp2k.so.
BIN_ROOT=/home/raid/crm98/cp2k_binaries/phy-cerberus

case "$TARGET_BRANCH" in
  feature-nnp-verlet-cells)
      CP2K_EXE="$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-verlet-cells/lib"
      LABEL="feature-nnp-verlet-cells"
      PROJECT_ROOT="/home/raid/crm98/cp2k"
      BENCHMARK_ROOT="/home/raid/crm98/cp2k-benchmarks/cp2k_optimized/NNP"
      OUTDIR_PARENT="cp2k_feature_verlet_cells"
      ;;
  feature-nnp-native-spline)
      CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-native-spline/lib"
      LABEL="feature-nnp-native-spline"
      PROJECT_ROOT="/home/raid/crm98/cp2k"
      BENCHMARK_ROOT="/home/raid/crm98/cp2k-benchmarks/cp2k_optimized/NNP"
      OUTDIR_PARENT="cp2k_feature_native_spline"
      ;;
  master|*)
      CP2K_EXE="$BIN_ROOT/master/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/master/lib"
      LABEL="upstream-master"
      PROJECT_ROOT="/home/raid/crm98/cp2k-upstream-master"
      BENCHMARK_ROOT="/home/raid/crm98/cp2k-benchmarks/cp2k_master/NNP"
      OUTDIR_PARENT="cp2k_master"
      ;;
esac

OUTDIR="/local/data/public/crm98/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_size_scaling_${LABEL}_${TIMESTAMP}"

# Ensure this branch's libcp2k.so.2026.1 is loaded, not another branch's.
export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"

# --- Variables ---
MPI_RANKS=36
export OMP_NUM_THREADS=1
STEPS=100

mkdir -p "$OUTDIR"

CSV_FILE="${OUTDIR}/results_size_scaling_${LABEL}_${TIMESTAMP}.csv"


cat <<EOF >"$CSV_FILE"
# branch:    $LABEL
# exe:       $CP2K_EXE
# MPI ranks: $MPI_RANKS
# steps:     $STEPS
# n_molecules,walltime_s
EOF

echo "Starting NNP size scaling benchmark on branch: $LABEL"
printf "%-12s %-15s %-15s\n" "Molecules" "Repeated Tiles" "Walltime (s)"
echo "---------------------------------------------------"

declare -A MULT
MULT[64]="1 1 1"
MULT[256]="2 2 1"
MULT[512]="2 2 2"
MULT[1024]="4 2 2"
MULT[2048]="4 4 2"
MULT[4096]="4 4 4"

for size in 64 256 512 1024 2048 4096; do
  rundir="${OUTDIR}/N_${size}"
  mkdir -p "$rundir"

  export current_size=$size
  export mx=$(echo "${MULT[$size]}" | awk '{print $1}')
  export my=$(echo "${MULT[$size]}" | awk '{print $2}')
  export mz=$(echo "${MULT[$size]}" | awk '{print $3}')
  export target_file="${rundir}/run.inp"
  export base_inp=$BASE_INP

python3 - <<'PYEOF'
import os, re

base = os.environ['base_inp']
target = os.environ['target_file']
size = os.environ['current_size']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']

with open(base, "r") as f:
    txt = f.read()

if size != "64":
    txt = re.sub(r'(&CELL\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.IGNORECASE)
    if "&TOPOLOGY" in txt:
        txt = re.sub(r'(&TOPOLOGY\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.IGNORECASE)
    else:
        txt = re.sub(r'(&SUBSYS\s*\n)', rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n', txt, flags=re.IGNORECASE)

with open(target, "w") as f:
    f.write(txt)
PYEOF

  ln -sfn "$NNP_DATA" "${rundir}/NNP"

  # Execute CP2K inside the run directory
  (cd "$rundir" && mpiexec -n "$MPI_RANKS" --bind-to core "$CP2K_EXE" -i run.inp > cp2k_${size}.out 2>&1) || true

  if grep -q "PROGRAM ENDED" "${rundir}/cp2k_${size}.out" 2>/dev/null; then
    wt=$(grep -E "^ CP2K +[0-9]" "${rundir}/cp2k_${size}.out" | awk '{print $NF}' | tail -1)
    printf "%-12s %-15s %-15s\n" "$size" "${MULT[$size]}" "$wt"
    echo "$size,$wt" >> "$CSV_FILE"
  else
    printf "%-12s %-15s %-15s\n" "$size" "${MULT[$size]}" "FAILED"
  fi

done
