#!/usr/bin/env bash
# Extensive 3-way comparison driver: master vs native-spline vs chebyshev.
# Captures the three axes the pros/cons table needs:
#   1. hybrid MPI x OMP decomposition at fixed 76 cores
#      -> does the branch convert threads to speed, or idle them?
#   2. aggregate node memory per decomposition (peak/proc x mpi)
#      -> cost of many ranks (replication) vs few-ranks/many-threads
#   3. OMP thread scaling at 1 rank (chebyshev ladder)
# Pure-MPI size/strong scaling already exists; native-spline core is re-run
# fresh to remove the May-vs-June cross-month caveat.
#
# Note: low-MPI decompositions on master/native-spline are pathologically slow
# (NNP path is serial per rank: 1 rank x N=4096 ~ 260 s/step). Those branches
# gain nothing from threads, so at N=4096 they are run ONLY at high-MPI points;
# the full hybrid ladder runs for all branches at N=1024 (where 1x76 ~ 11 min)
# and for chebyshev at N=4096.
#
#SBATCH -J NNP_hybrid_mem
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=05:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_hybrid_mem_%j.out

set -uo pipefail
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
SCALING=/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling
BRANCHES="master feature-nnp-native-spline feature-nnp-chebyshev"
FULL_LADDER="76x1 38x2 19x4 4x19 2x38 1x76"
HIMPI_ONLY="76x1 38x2 19x4"

echo "=== BINARY VERIFICATION ==="
for b in $BRANCHES; do
   exe="$BIN_ROOT/$b/cp2k.psmp"
   [[ -x "$exe" ]] || { echo "!! missing $exe"; exit 1; }
   md5sum "$exe"
   [[ -f "$BIN_ROOT/$b/PROVENANCE.txt" ]] && sed 's/^/    /' "$BIN_ROOT/$b/PROVENANCE.txt"
done
cd "$SCALING"

echo ""; echo "############ HYBRID DECOMP SWEEP — N=1024 (full ladder, all branches) ############"
for b in $BRANCHES; do
   echo "---- $b ----"
   N_MOLECULES=1024 STEPS=20 N_REPS=2 DECOMP_LIST="$FULL_LADDER" ./run_nnp_decomp_sweep_slurm.sh "$b"
done

echo ""; echo "############ HYBRID DECOMP SWEEP — N=4096 ############"
echo "---- feature-nnp-chebyshev (full ladder) ----"
N_MOLECULES=4096 STEPS=20 N_REPS=2 DECOMP_LIST="$FULL_LADDER" ./run_nnp_decomp_sweep_slurm.sh feature-nnp-chebyshev
for b in master feature-nnp-native-spline; do
   echo "---- $b (high-MPI only) ----"
   N_MOLECULES=4096 STEPS=15 N_REPS=2 DECOMP_LIST="$HIMPI_ONLY" ./run_nnp_decomp_sweep_slurm.sh "$b"
done

echo ""; echo "############ chebyshev OMP THREAD LADDER (1 rank, N=1024) ############"
for th in 1 2 4 8 16; do
   TOTAL_CORES=$th N_MOLECULES=1024 STEPS=20 N_REPS=2 DECOMP_LIST="1x${th}" \
      ./run_nnp_decomp_sweep_slurm.sh feature-nnp-chebyshev
done

echo ""; echo "############ FRESH PURE-MPI STRONG SCALING (native-spline, June bridge) ############"
N_MOLECULES=1024 STEPS=50 N_REPS=3 ./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline

SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_R=/home/crm98/cp2k-benchmarks/results
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs "$SCRATCH/" "$HOME_R/"
echo ""; echo "=== done; CSVs synced to $HOME_R ==="
