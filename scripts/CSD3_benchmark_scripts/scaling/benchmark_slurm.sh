#!/usr/bin/env bash
# CSD3 / Peta4-IceLake NNP benchmark driver (Report 2 scope):
# upstream master vs feature/nnp-chebyshev.
#
# Runs, on one exclusive node:
#   1. size scaling   (full node, N = 64..4096)        both binaries
#   2. strong scaling (N = 1024, cores 1..76)          both binaries
#   3. OMP thread scaling (1 rank, threads 1..16)      chebyshev only
#      (master's NNP path is pure-MPI; threads would idle)
#
# Binaries are NOT built here: they come from the provenance-stamped
# snapshot cache (see PROVENANCE.txt next to each cp2k.psmp).  Rebuild +
# re-snapshot deliberately stays a manual, user-driven step so a queued
# benchmark can never silently measure a different binary.

#SBATCH -J NNP_scaling
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
# 8 h: size scaling ~1.5 h + core scaling ~3.5 h (the N=1024 1-core points
# alone are ~80 min) + OMP scaling ~1 h.  3 h was not enough (job 30352037).
#SBATCH --time=08:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_scaling_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

# rhel8/default-icl alone gives compilers + MPI but NOT Intel MKL; cp2k.psmp
# loads libmkl_intel_thread.so.2 dynamically at startup, so the MKL module
# is mandatory or every srun exits 127.
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl
module list

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3

echo "=== BINARY VERIFICATION ==="
for label in master feature-nnp-chebyshev; do
   exe="$BIN_ROOT/$label/cp2k.psmp"
   [[ -x "$exe" ]] || { echo "!! missing $exe -- snapshot the binary first"; exit 1; }
   md5sum "$exe"
   [[ -f "$BIN_ROOT/$label/PROVENANCE.txt" ]] && sed 's/^/    /' "$BIN_ROOT/$label/PROVENANCE.txt"
done
if cmp -s "$BIN_ROOT/master/cp2k.psmp" "$BIN_ROOT/feature-nnp-chebyshev/cp2k.psmp"; then
   echo "!! the two binaries are identical -- aborting"; exit 1
fi
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

echo "=== SIZE SCALING ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh feature-nnp-chebyshev

echo ""
echo "=== STRONG (CORE) SCALING ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh feature-nnp-chebyshev

echo ""
echo "=== OMP THREAD SCALING (chebyshev, 1 rank) ==="
OMP_LIST="1 2 4 8 16" ./run_nnp_omp_thread_scaling_slurm.sh

SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS (full results remain on $SCRATCH_RESULTS)"
