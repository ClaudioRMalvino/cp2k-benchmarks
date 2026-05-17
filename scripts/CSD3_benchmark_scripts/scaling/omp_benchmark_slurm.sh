#!/usr/bin/env bash
# CSD3 OpenMP-scaling benchmark for feature/nnp-native-spline-omp:
# thread-scaling sweep, size-scaling sweep, and the MPI:OMP decomposition
# sweep, all on a single node.

#SBATCH -J NNP_omp_scaling
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_omp_scaling_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

. /etc/profile.d/modules.sh
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
mkdir -p "$BIN_ROOT/feature-nnp-native-spline-omp/lib"
REBUILD=${REBUILD:-0}

CACHE="$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
if [[ "$REBUILD" -eq 1 || ! -x "$CACHE" ]]; then
   echo "==> Building feature-nnp-native-spline-omp"
   git -C /home/crm98/cp2k_optimized stash || true
   git -C /home/crm98/cp2k_optimized checkout feature/nnp-native-spline-omp
   bash /home/crm98/cp2k-benchmarks/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_rebuild.sh
   git -C /home/crm98/cp2k_optimized stash pop || true
   if [[ -x /home/crm98/cp2k_optimized/install/bin/cp2k.psmp ]]; then
      cp /home/crm98/cp2k_optimized/install/bin/cp2k.psmp "$CACHE"
   else
      echo "!! build did not produce cp2k.psmp - aborting"; exit 1
   fi
else
   echo "==> Using cached binary: $CACHE"
fi

cp /home/crm98/cp2k_master/tools/toolchain/install/setup "$BIN_ROOT/setup"

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

# This driver runs the MPI:OMP decomposition sweep itself (below); tell
# run_nnp_omp_size_scaling_slurm.sh so it does not also chain it.
export DECOMP_VIA_DRIVER=1

echo "=== OMP THREAD SCALING ==="
./run_nnp_omp_thread_scaling_slurm.sh

echo ""
echo "=== OMP SIZE SCALING ==="
./run_nnp_omp_size_scaling_slurm.sh

echo ""
echo "=== MPI:OMP DECOMPOSITION SWEEP ==="
./run_nnp_decomposition_sweep_slurm.sh

SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS (full results remain on $SCRATCH_RESULTS)"
