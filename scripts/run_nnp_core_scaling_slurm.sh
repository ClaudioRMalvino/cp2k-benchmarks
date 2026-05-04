#!/usr/bin/env bash
set -euo pipefail

# Per-branch binary cache populated by benchmark_slurm.sh.  Each branch has
# its own cp2k.psmp + lib/, so LD_LIBRARY_PATH cannot accidentally pick up
# another branch's libcp2k.so.  Lives on scratch (/local/data/public) — too
# big for /home (libcp2k.so.* alone is hundreds of MB per branch).
BIN_ROOT=/local/data/public/crm98/cp2k_binaries/phy-cerberus

set +u
# --- Load the Toolchain Environment ---
source "$BIN_ROOT/setup"
set -u

# --- Pass the branch as an argument (defaults to master) ---
# Accepted values: master, feature-nnp-verlet-cells, feature-nnp-native-spline
TARGET_BRANCH=${1:-master}
TIMESTAMP=$(date +%d-%m_%H-%M)

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

OUTDIR="/local/data/public/crm98/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_core_scaling_${LABEL}_${TIMESTAMP}"

# Ensure this branch's libcp2k.so.2026.1 is loaded, not another branch's.
export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"

# --- Variables ---
N_MOLECULES=${N_MOLECULES:-64}
MPI_LIST=${MPI_LIST:-"1 2 4 8 16 32"}
OMP_LIST=${OMP_LIST:-"1"}
STEPS=${STEPS:-100}

mkdir -p "$OUTDIR"

CSV_FILE="${OUTDIR}/results_core_scaling_${LABEL}_${TIMESTAMP}.csv"

cat <<EOF >"$CSV_FILE"
# branch:      $LABEL
# exe:         $CP2K_EXE
# N_molecules: $N_MOLECULES
# steps:       $STEPS
# mpi_ranks,omp_threads,total_cores,walltime_s,speedup
EOF

echo "Starting core scaling benchmark (System size: $N_MOLECULES) on branch: $LABEL"
printf "%-10s %-10s %-12s %-14s %-10s %-10s\n" "MPI_Ranks" "OMP_Threads" "Total_Cores" "Walltime(s)" "Speedup" "Efficiency"
echo "------------------------------------------------------------------------"

BASELINE_TIME=""
BASELINE_CORES=""

for omp in $OMP_LIST; do
  export OMP_NUM_THREADS=$omp
  for mpi in $MPI_LIST; do
    total=$((mpi * omp))
    rundir="${OUTDIR}/mpi${mpi}_omp${omp}"
    mkdir -p "$rundir"

    cp "$BASE_INP" "${rundir}/run.inp"
    ln -sfn "$NNP_DATA" "${rundir}/NNP"

    printf "%-10s %-10s %-12s " "$mpi" "$omp" "$total"

    (cd "$rundir" && mpiexec -n "$mpi" --bind-to core "$CP2K_EXE" -i run.inp >cp2k.out 2>&1) || true

    if grep -q "PROGRAM ENDED" "${rundir}/cp2k.out" 2>/dev/null; then
      wt=$(grep -E "^ CP2K +[0-9]" "${rundir}/cp2k.out" | awk '{print $NF}' | tail -1)

      if [[ -z "$BASELINE_TIME" ]]; then
        BASELINE_TIME=$wt
        BASELINE_CORES=$total
      fi

      speedup=$(awk "BEGIN {printf \"%.2f\", $BASELINE_TIME / $wt}")
      efficiency=$(awk "BEGIN {printf \"%.1f%%\", 100 * $BASELINE_TIME / $wt / ($total / $BASELINE_CORES)}")

      printf "%-14s %-10s %-10s\n" "$wt" "${speedup}x" "$efficiency"
      echo "$mpi,$omp,$total,$wt,$speedup" >>"$CSV_FILE"
    else
      printf "%-14s %-10s %-10s\n" "FAILED" "N/A" "N/A"
    fi
  done
done
