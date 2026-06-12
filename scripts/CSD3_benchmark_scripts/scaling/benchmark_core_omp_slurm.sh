#!/usr/bin/env bash
# Full Report-2 single-node battery -- the same measurement sets Report 1
# ran for native-spline / native-spline-omp, EVERY ONE AT STEPS=100 (Report
# 1 mixed 50- and 100-step runs; Report 2 is uniform):
#   1. pure-MPI size scaling   (N=64..4096, 76 ranks)        chebyshev
#   2. pure-MPI core scaling   (N=1024, cores 1..76)         chebyshev
#   3. pure-MPI core scaling   (N=1024, cores 1..76)         MASTER
#      -- re-measured because the Report 1 master curve used STEPS=50;
#      reusing it would break the uniform-steps rule.  One-off ~2 h.
#   4. OMP thread scaling      (N=64, 1 rank, 1..16 threads) chebyshev
#   5. OMP x size 2D sweep     (1 rank, OMP 1..8, N=64..2048) chebyshev
# Master size scaling from 10-06 was already STEPS=100 and stays reused.
#
# Walltime budget (conservative, pre-rework chebyshev speeds): size ~1 h +
# cheby core ~1.5 h + master core ~2 h + OMP thread ~0.75 h + OMP x size
# ~4 h (the N=2048 single-thread row alone is ~1 h at 100 steps) => ~9.5 h;
# 11 h walltime for headroom (SL3 cap is 12 h).
#
# Submit with:  sbatch benchmark_core_omp_slurm.sh

#SBATCH -J NNP_suite_100
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=11:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_suite_100_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl
module list

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3

echo "=== BINARY VERIFICATION ==="
for b in feature-nnp-chebyshev master; do
   exe="$BIN_ROOT/$b/cp2k.psmp"
   [[ -x "$exe" ]] || { echo "!! missing $exe -- snapshot the binary first"; exit 1; }
   md5sum "$exe"
   [[ -f "$BIN_ROOT/$b/PROVENANCE.txt" ]] && \
      sed 's/^/    /' "$BIN_ROOT/$b/PROVENANCE.txt"
done
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

echo "=== SIZE SCALING (chebyshev, pure MPI) ==="
./run_nnp_size_scaling_slurm.sh feature-nnp-chebyshev

echo ""
echo "=== STRONG (CORE) SCALING (chebyshev) ==="
N_REPS=3 ./run_nnp_core_scaling_slurm.sh feature-nnp-chebyshev

echo ""
echo "=== STRONG (CORE) SCALING (master, 100-step baseline) ==="
N_REPS=3 ./run_nnp_core_scaling_slurm.sh master

echo ""
echo "=== OMP THREAD SCALING (chebyshev, 1 rank) ==="
OMP_LIST="1 2 4 8 16" ./run_nnp_omp_thread_scaling_slurm.sh

echo ""
echo "=== OMP x SIZE SWEEP (chebyshev, 1 rank) ==="
./run_nnp_omp_size_scaling_slurm.sh

SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS (full results remain on $SCRATCH_RESULTS)"
