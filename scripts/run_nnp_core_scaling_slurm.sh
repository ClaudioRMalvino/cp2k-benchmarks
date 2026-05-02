#!/usr/bin/env bash
set -euo pipefail

set +u
# --- Load the Toolchain Environment ---
source /home/raid/crm98/cp2k_binaries/phy-cerberus/setup
set -u

# --- Pass the branch as an argument (defaults to master) ---
TARGET_BRANCH=${1:-master} 
TIMESTAMP=$(date +%d-%m_%H-%M)

if [ "$TARGET_BRANCH" == "optimized" ]; then
    CP2K_EXE="/home/raid/crm98/cp2k_binaries/phy-cerberus/cp2k_feature_verlet_cells.psmp"
    LABEL="feature-nnp-verlet-cells"
    PROJECT_ROOT="/home/raid/crm98/cp2k"
    BENCHMARK_ROOT="/home/raid/crm98/cp2k-benchmarks/cp2k_optimized/NNP"
    OUTDIR="/local/data/public/crm98/cp2k-benchmarks/results/cp2k_optimized/NNP/NNP_core_scaling_${LABEL}_${TIMESTAMP}"
    INSTALL_LIB="/local/data/public/crm98/original_cp2k/install/lib"
else
    CP2K_EXE="/home/raid/crm98/cp2k_binaries/phy-cerberus/cp2k_master.psmp"
    LABEL="upstream-master"
    PROJECT_ROOT="/home/raid/crm98/cp2k-upstream-master"
    BENCHMARK_ROOT="/home/raid/crm98/cp2k-benchmarks/cp2k_master/NNP"
    OUTDIR="/local/data/public/crm98/cp2k-benchmarks/results/cp2k_master/NNP/NNP_core_scaling_${LABEL}_${TIMESTAMP}"
    INSTALL_LIB="/local/data/public/crm98/cp2k-buildtree/install/lib"
fi
# Ensure this branch's libcp2k.so.2026.1 is loaded, not the other branch's
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
    
    # Copy your pre-configured file into the run directory and name it run.inp
    cp "$BASE_INP" "${rundir}/run.inp"
    
    # Symlink the required NNP data directory
    ln -sfn "$NNP_DATA" "${rundir}/NNP"

    printf "%-10s %-10s %-12s " "$mpi" "$omp" "$total"

    # Execute CP2K
    (cd "$rundir" && mpiexec -n "$mpi" --bind-to core "$CP2K_EXE" -i run.inp >cp2k.out 2>&1) || true

    # 1. Check if the run finished successfully
    if grep -q "PROGRAM ENDED" "${rundir}/cp2k.out" 2>/dev/null; then
      
      # 2. Extract the Walltime using grep and awk
      wt=$(grep -E "^ CP2K +[0-9]" "${rundir}/cp2k.out" | awk '{print $NF}' | tail -1)

      # 3. Save the 1-core baseline time to do the math against
      if [[ -z "$BASELINE_TIME" ]]; then
        BASELINE_TIME=$wt
        BASELINE_CORES=$total
      fi

      # 4. Calculate Speedup and Efficiency
      speedup=$(awk "BEGIN {printf \"%.2f\", $BASELINE_TIME / $wt}")
      efficiency=$(awk "BEGIN {printf \"%.1f%%\", 100 * $BASELINE_TIME / $wt / ($total / $BASELINE_CORES)}")

      # 5. Print to the terminal log
      printf "%-14s %-10s %-10s\n" "$wt" "${speedup}x" "$efficiency"
      
      # 6. SAVE TO THE CSV
      echo "$mpi,$omp,$total,$wt,$speedup" >>"$CSV_FILE"
      
    else
      printf "%-14s %-10s %-10s\n" "FAILED" "N/A" "N/A"
    fi
  done
done
