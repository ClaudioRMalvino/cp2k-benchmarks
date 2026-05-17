#!/usr/bin/env bash
# Size-scaling sweep only for all three branches (N = 64 -> 4096 H2O at
# full-node parallelism). Reads cached binaries from $BIN_ROOT; does not
# rebuild.

#SBATCH -J NNP_size_rerun
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_size_rerun_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

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
      echo "!! missing binary: $BIN_ROOT/$b/cp2k.psmp"
      echo "   run rebuild_all_branches.sh first."
      exit 1
   fi
done

echo "=== BINARY VERIFICATION ==="
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

echo "=== SIZE SCALING (full-node, N = 64 -> 4096 H2O) ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline-omp

SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo
echo "CSVs synced to $HOME_RESULTS"
