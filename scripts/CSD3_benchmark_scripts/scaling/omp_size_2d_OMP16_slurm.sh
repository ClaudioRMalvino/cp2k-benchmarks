#!/usr/bin/env bash
# Re-run of the OMP 2D size sweep (OMP_LIST x SIZE_LIST) for
# feature/nnp-native-spline-omp. Invokes run_nnp_omp_size_scaling_slurm.sh
# only; thread scaling, full-node size sweep, and the build cache are not
# re-run.

#SBATCH -J NNP_omp_size_2d_OMP16
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_omp_size_2d_OMP16_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl
module list

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
if [[ ! -x "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp" ]]; then
   echo "!! missing binary: $BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
   echo "   run benchmark_slurm.sh (or set REBUILD=1) to populate the cache first"
   exit 1
fi

md5sum "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

# Tell the OMP size script that no MPI:OMP decomposition driver is chaining
# it (matches the omp_benchmark_slurm.sh pattern - prevents double-run).
export DECOMP_VIA_DRIVER=1
export OMP_LIST="16"
export SIZE_LIST="64 256 512 1024 2048"

echo "=== OMP 2D size sweep (OMP=$OMP_LIST x N in $SIZE_LIST) ==="
./run_nnp_omp_size_scaling_slurm.sh

SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS"
