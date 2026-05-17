#!/usr/bin/env bash
# CSD3 / Peta4-IceLake NNP benchmark driver: build three branches into a
# per-branch binary cache on /rds scratch, then run size- and strong-scaling
# sweeps for each.

#SBATCH -J NNP_scaling
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
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
mkdir -p "$BIN_ROOT/master/lib"
mkdir -p "$BIN_ROOT/feature-nnp-native-spline/lib"
mkdir -p "$BIN_ROOT/feature-nnp-native-spline-omp/lib"

REBUILD=${REBUILD:-0}

build_branch() {
   local branch_label=$1
   local repo=$2
   local git_branch=$3
   local build_script=$4

   local cache_dir="$BIN_ROOT/$branch_label"
   if [[ "$REBUILD" -eq 0 && -x "$cache_dir/cp2k.psmp" ]]; then
      echo "==> Cached binary present for $branch_label, skipping build"
      return
   fi

   echo "==> Building $branch_label"
   if [[ -n "$git_branch" ]]; then
      git -C "$repo" stash || true
      git -C "$repo" checkout "$git_branch"
   fi
   bash "$build_script"
   local install_dir="$repo/install"
   # libcp2k.so* is not produced (BUILD_SHARED_LIBS=OFF in the build configs);
   # guard the copy so a missing .so does not abort the driver.
   if [[ -x "$install_dir/bin/cp2k.psmp" ]]; then
      cp "$install_dir/bin/cp2k.psmp" "$cache_dir/cp2k.psmp"
   else
      echo "!! $branch_label build did not produce cp2k.psmp - cache left untouched"
   fi
   cp -P "$install_dir/lib"/libcp2k.so* "$cache_dir/lib/" 2>/dev/null || true
   if [[ -n "$git_branch" ]]; then
      git -C "$repo" stash pop || true
   fi
}

CP2K_OPT_REPO=/home/crm98/cp2k_optimized
CP2K_MASTER_REPO=/home/crm98/cp2k_master
OPT_REBUILD=/home/crm98/cp2k-benchmarks/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_rebuild.sh
MASTER_BUILD=/home/crm98/cp2k-benchmarks/cp2k_master/CSD3_build_scripts/cp2k_CSD3_master_build.sh

build_branch "master"                          "$CP2K_MASTER_REPO" ""                              "$MASTER_BUILD"
build_branch "feature-nnp-native-spline"        "$CP2K_OPT_REPO"    "feature/nnp-native-spline"     "$OPT_REBUILD"
build_branch "feature-nnp-native-spline-omp"    "$CP2K_OPT_REPO"    "feature/nnp-native-spline-omp" "$OPT_REBUILD"

cp /home/crm98/cp2k_master/tools/toolchain/install/setup "$BIN_ROOT/setup"

# matching md5sums would mean a build was silently skipped or the wrong
# source tree was used
echo "=== BINARY VERIFICATION ==="
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/

echo "=== SIZE SCALING ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline-omp

echo ""
echo "=== STRONG (CORE) SCALING ==="
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
