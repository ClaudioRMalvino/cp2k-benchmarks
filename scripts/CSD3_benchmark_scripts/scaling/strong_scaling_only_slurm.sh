#!/usr/bin/env bash
#! CSD3 / Peta4-IceLake -- strong (core) scaling sweep ONLY, at N = 1024 H2O.
#!
#! Companion to benchmark_slurm.sh.  Re-runs just the strong scaling for the
#! three branches at a production-typical system size (1024 H2O = 3072 atoms),
#! where master's O(N^2) neighbour search is the visible bottleneck.  The
#! earlier sweep used N = 64 H2O = 192 atoms (cerberus default), at which
#! master and feature/nnp-native-spline differ by only ~14% per core because
#! the algorithmic gap is too small to show up over communication noise.
#!
#! Size scaling, OMP thread scaling, OMP 2D size sweep, and the build cache
#! are NOT re-run -- this script reads the cached binaries from BIN_ROOT.
#!
#! STEPS is dropped from 100 to 50 to keep the 1-core master point tractable:
#! master is O(N^2) so the 1-core run at N=1024 is ~256x slower per step than
#! at N=64.  STEPS=50 trims walltime in half without meaningfully harming the
#! per-step timing precision (still 5 reps of 50 steps = 250 MD steps each).

#SBATCH -J NNP_strong_N1024
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_strong_N1024_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

# ---------------------------------------------------------------------------
# Runtime environment -- same as benchmark_slurm.sh (MKL module is mandatory
# or cp2k.psmp fails to find libmkl_intel_thread.so.2 at startup).
# ---------------------------------------------------------------------------
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl
module list

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3

# Bail early if any binary is missing -- this script does not rebuild.
for b in master feature-nnp-native-spline feature-nnp-native-spline-omp; do
   if [[ ! -x "$BIN_ROOT/$b/cp2k.psmp" ]]; then
      echo "!! missing binary: $BIN_ROOT/$b/cp2k.psmp"
      echo "   run benchmark_slurm.sh (or set REBUILD=1) to populate the cache first."
      exit 1
   fi
done

echo "=== BINARY VERIFICATION ==="
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

# ---------------------------------------------------------------------------
# Strong scaling sweep.  Defaults to N=1024 / STEPS=50 (the canonical run);
# override via the submission environment to sweep additional system sizes:
#   N_MOLECULES=256  STEPS=100 sbatch --time=01:30:00 --job-name=NNP_strong_N256  strong_scaling_only_slurm.sh
#   sbatch strong_scaling_only_slurm.sh                                                                       # default N=1024
#   N_MOLECULES=2048 STEPS=50  sbatch --time=04:00:00 --job-name=NNP_strong_N2048 strong_scaling_only_slurm.sh
# The plotting script auto-discovers the latest CSV per (branch, N) by mtime.
# ---------------------------------------------------------------------------
export N_MOLECULES=${N_MOLECULES:-1024}
export STEPS=${STEPS:-50}

echo "=== STRONG (CORE) SCALING at N=$N_MOLECULES H2O, STEPS=$STEPS ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline-omp

# ---------------------------------------------------------------------------
# Sync the new CSVs from /rds scratch to home for plotting.
# ---------------------------------------------------------------------------
SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS (full results remain on $SCRATCH_RESULTS)"
