#!/usr/bin/env bash
#! CSD3 / Peta4-IceLake -- size scaling sweep ONLY for all three branches.
#!
#! Companion to benchmark_slurm.sh / strong_scaling_only_slurm.sh.  Re-runs
#! just the size scaling (N = 64 → 4096 H2O) for master, native-spline, and
#! native-spline-omp at full-node parallelism.  Needed any time the cached
#! binaries change (e.g. after rebuilding with -g and -xCORE-AVX512).
#!
#! The binaries are assumed to already be in $BIN_ROOT/<branch>/cp2k.psmp;
#! this script does NOT rebuild.  If a binary is missing it bails early.
#!
#! Walltime budget (estimated from previous full-node N = 4096 runs):
#!   * master            ~7 min total sweep (dominated by N=4096 at 329 s)
#!   * native-spline     ~1 min total sweep
#!   * native-spline-omp ~1 min total sweep
#! Plus 5 timed reps per size + 1 warm-up; the 6 sizes × 5 reps + overhead
#! keeps the whole job under ~1.5 h.  Bumping to 2 h gives a margin if the
#! new -xCORE-AVX512 build performs differently.

#SBATCH -J NNP_size_rerun
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_size_rerun_%j.out

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

echo "=== SIZE SCALING (full-node, N = 64 → 4096 H2O) ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline-omp

# --- sync CSVs back to home for plotting -----------------------------------
SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_RESULTS=/home/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo
echo "CSVs synced to $HOME_RESULTS"
