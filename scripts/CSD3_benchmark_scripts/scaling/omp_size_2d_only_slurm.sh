#!/usr/bin/env bash
#! CSD3 / Peta4-IceLake -- 2D OMP-x-size sweep ONLY, for feature/nnp-native-spline-omp.
#!
#! Companion to omp_benchmark_slurm.sh.  Re-runs just the OMP 2-D size sweep
#! (the OMP_LIST x SIZE_LIST grid) because the previous run was killed by
#! SLURM mid-grid -- the parent job (29351927) had a 2 h budget for three
#! sub-sweeps and timed out after the thread sweep + the first ~2.5 rows of
#! the 2-D grid, leaving OMP=4 N>=1024 and the entire OMP=8 row missing.
#!
#! Thread scaling, the full-node size sweep, and the build cache are NOT
#! re-run -- this script only invokes run_nnp_omp_size_scaling_slurm.sh.
#!
#! Walltime budget (estimated from the previous run's per-point walltimes):
#!   OMP=1 row : ~20 min     OMP=2 row : ~12 min
#!   OMP=4 row :  ~8 min     OMP=8 row :  ~5 min
#! Total ~45 min for the default OMP_LIST="1 2 4 8" x SIZE_LIST="64..2048";
#! 1 h walltime gives a comfortable margin.  Override with longer if extending
#! either list, e.g.:   OMP_LIST="1 2 4 8 16" SIZE_LIST="64 256 512 1024 2048 4096" sbatch ...

#SBATCH -J NNP_omp_size_2d
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_omp_size_2d_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

# --- environment -----------------------------------------------------------
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
   echo "   run benchmark_slurm.sh (or set REBUILD=1) to populate the cache first."
   exit 1
fi

md5sum "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

# Tell the OMP size script that no MPI:OMP decomposition driver is chaining
# it (matches the omp_benchmark_slurm.sh pattern -- prevents double-run).
export DECOMP_VIA_DRIVER=1

echo "=== OMP 2D size sweep (omp x N) ==="
./run_nnp_omp_size_scaling_slurm.sh

# --- sync CSVs back to home for plotting -----------------------------------
SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS"
