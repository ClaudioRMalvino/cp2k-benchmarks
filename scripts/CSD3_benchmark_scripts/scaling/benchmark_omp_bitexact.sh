#!/usr/bin/env bash
# Single-binary benchmark of the bit-exact atom-level OpenMP graft
# (dhruv-omp-bitexact) on Dhruv's NNP kernel. Runs ONLY the new binary, on the
# SAME two axes / parameters as the earlier dhruv-vs-mace jobs (31038227 +
# 31038352), so the CSVs overlay directly onto the existing baselines:
#   * dhruv-cell-list  (baseline, his fine-grained inner-loop OMP)  -- already measured
#   * feature-nnp-mace (reference outer-loop OMP)                   -- already measured
# We do NOT re-run those (unchanged binaries) -> saves node hours.
#
# Axis 1 (headline): OMP THREAD scaling, 1 MPI rank, threads {1..76}, N=1024.
#   -> does the bit-exact outer-loop graft scale where his went flat (~0.8x)?
# Axis 2: PURE-MPI strong scaling, N=1024, cores {1..76}, 1 OMP.
#   -> did the per-thread arc-clone / fold refactor cost anything MPI-only?
# Both: 100 MD steps, 5 timed reps + 1 warm-up, t/step from qs_mol_dyn_low.
# No thread pinning -> identical conditions to the dhruv-cell-list 0.8x run.
#
#SBATCH -J NNP_omp_bitexact
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=08:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_omp_bitexact_%j.out

set -uo pipefail
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
SCALING=/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling
LABEL=dhruv-omp-bitexact
PROOT=/rds/user/$USER/hpc-work/cp2k_dhruv_omp_worktree
OPARENT=cp2k_dhruv_omp_bitexact
LIST_R1="1 2 4 8 16 19 32 38 76"   # report-1 NUMA-aware set, matches the baselines

EXE="$BIN_ROOT/$LABEL/cp2k.psmp"
echo "=== BINARY VERIFICATION (md5) ==="
[[ -x "$EXE" ]] || { echo "!! missing $EXE"; exit 1; }
md5sum "$EXE"
cd "$SCALING"

echo ""
echo "############ AXIS 1: OMP THREAD SCALING (1 MPI rank, N=1024) ############"
echo "threads: $LIST_R1 | steps=100 | reps=5(+1 warm-up)"
TARGET_LABEL="$LABEL" \
PROJECT_ROOT_OVERRIDE="$PROOT" \
OUTDIR_PARENT_OVERRIDE="$OPARENT" \
MPI_RANKS=1 N_MOLECULES=1024 STEPS=100 N_REPS=5 OMP_LIST="$LIST_R1" \
   ./run_nnp_omp_thread_scaling_slurm.sh

echo ""
echo "############ AXIS 2: PURE-MPI STRONG SCALING (N=1024, 1 OMP) ############"
echo "cores: $LIST_R1 | steps=100 | reps=5(+1 warm-up)"
CP2K_EXE_OVERRIDE="$EXE" \
INSTALL_LIB_OVERRIDE="$BIN_ROOT/$LABEL/lib" \
LABEL_OVERRIDE="$LABEL" \
PROJECT_ROOT_OVERRIDE="$PROOT" \
OUTDIR_PARENT_OVERRIDE="$OPARENT" \
N_MOLECULES=1024 STEPS=100 N_REPS=5 CORE_LIST="$LIST_R1" \
   ./run_nnp_core_scaling_slurm.sh "$LABEL"

# Mirror CSVs from scratch to the in-repo results tree (CSVs only, no run dirs).
SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_R=/home/crm98/cp2k-benchmarks/results
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs "$SCRATCH/" "$HOME_R/"
echo ""; echo "=== done; dhruv-omp-bitexact CSVs synced to $HOME_R/$OPARENT ==="
