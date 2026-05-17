#!/usr/bin/env bash
# Strong (core) scaling sweep only, at N = 1024 H2O. Reads cached binaries
# from $BIN_ROOT; does not rebuild. STEPS dropped from 100 to 50 to keep the
# 1-core master point tractable (master is O(N^2) so the 1-core run at
# N=1024 is ~256x slower per step than at N=64).

#SBATCH -J NNP_strong_3branch_N1024
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_strong_N1024_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

# MKL module is mandatory or cp2k.psmp fails to find
# libmkl_intel_thread.so.2 at startup.
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl
module list

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3

for b in master feature-nnp-native-spline feature-nnp-native-spline-omp; do
   if [[ ! -x "$BIN_ROOT/$b/cp2k.psmp" ]]; then
      echo "!! missing binary: $BIN_ROOT/$b/cp2k.psmp"; exit 1
   fi
done

echo "=== BINARY VERIFICATION ==="
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

export N_MOLECULES=${N_MOLECULES:-1024}
export STEPS=${STEPS:-50}
export CORE_LIST=${CORE_LIST:-"1 2 4 8 16 19 32 38 76"}  # NUMA-safe; 19/38/76 = divisors of 76, rest <38 fit one socket

echo "=== STRONG (CORE) SCALING at N=$N_MOLECULES H2O, STEPS=$STEPS, CORES=$CORE_LIST ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline-omp

SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS (full results remain on $SCRATCH_RESULTS)"
