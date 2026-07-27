#!/usr/bin/env bash
# Head-to-head AXIS 3: OpenMP THREAD scaling — Dhruv's linked-cell NNP
# (pr/nnp-cpu-cell-list) vs my feature/nnp-mace. Single MPI rank, thread count
# swept 1 -> 76, fixed system size. This isolates the OpenMP implementations
# (his centre-level threading vs mine) with NO MPI in the mix, so it directly
# answers "is his OMP slower than ours?".
#
#   * MPI ranks  : 1 (fixed)
#   * OMP threads: {1,2,4,8,16,19,32,38,76}  (38 = one socket, 76 = both -> the
#                  >38 points expose cross-NUMA OMP degradation)
#   * size       : N = 1024 H2O (matches the pure-MPI axis for comparability)
#   * 100 MD steps, 5 timed reps + 1 discarded warm-up, t/step from qs_mol_dyn_low
#
# Both binaries built at the SAME AVX-2 profile + sequential MKL, so the only
# difference is the NNP source.
#
#SBATCH -J NNP_dhruv_vs_mace_omp
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=10:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_dhruv_vs_mace_omp_%j.out

set -uo pipefail
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
SCALING=/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling
OMP_LIST_R1="1 2 4 8 16 19 32 38 76"

echo "=== BINARY VERIFICATION (md5) ==="
for b in dhruv-cell-list feature-nnp-mace; do
   exe="$BIN_ROOT/$b/cp2k.psmp"
   [[ -x "$exe" ]] || { echo "!! missing $exe"; exit 1; }
   md5sum "$exe"
done
cd "$SCALING"

# Per-branch PROJECT_ROOT (data/NNP) + OUTDIR_PARENT so results are labelled right.
run_omp () {
   local label=$1 proot=$2 oparent=$3
   echo ""; echo "---- OMP thread scaling: $label ----"
   TARGET_LABEL="$label" \
   PROJECT_ROOT_OVERRIDE="$proot" \
   OUTDIR_PARENT_OVERRIDE="$oparent" \
   MPI_RANKS=1 N_MOLECULES=1024 STEPS=100 N_REPS=5 OMP_LIST="$OMP_LIST_R1" \
      ./run_nnp_omp_thread_scaling_slurm.sh
}

echo ""
echo "############ AXIS 3: OMP THREAD SCALING (1 MPI rank, N=1024) ############"
echo "threads: $OMP_LIST_R1 | steps=100 | reps=5(+1 warm-up)"
run_omp dhruv-cell-list  /home/crm98/cp2k_dhruv      cp2k_dhruv_cell_list
run_omp feature-nnp-mace /home/crm98/cp2k_optimized  cp2k_feature_mace

# Mirror CSVs from scratch to the in-repo results tree (CSVs only, no run dirs).
SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_R=/home/crm98/cp2k-benchmarks/results
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs "$SCRATCH/" "$HOME_R/"
echo ""; echo "=== done; Dhruv-vs-mace OMP thread-scaling CSVs synced to $HOME_R ==="
