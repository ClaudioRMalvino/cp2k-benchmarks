#!/usr/bin/env bash
#! Replicates Figure 2 of Iannuzzi et al., J. Phys. Chem. B 2026, 130, 1237.
#!   * Fig 2a: wall-time per MD step vs. core count for 64 / 512 / 4096 H2O
#!             on master, feature/nnp-native-spline, feature/nnp-native-spline-omp.
#!   * Fig 2b: O-O / O-H / H-H RDFs from a 64-H2O trajectory per branch.

#SBATCH -J NNP_ACS_PUB_REPLICATION
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=72
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_paper_fig2_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

# Per-branch binary cache: each branch has its own cp2k.psmp + lib/ so
# LD_LIBRARY_PATH cannot pick up another branch's libcp2k.so. Set REBUILD=1 to
# force rebuilds from the home-clones; otherwise the cache is reused.
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
   cp     "$install_dir/bin/cp2k.psmp" "$cache_dir/cp2k.psmp"
   cp -P  "$install_dir/lib"/libcp2k.so* "$cache_dir/lib/"
   if [[ -n "$git_branch" ]]; then
      git -C "$repo" stash pop || true
   fi
}

CP2K_OPT_REPO=/home/crm98/cp2k_optimized
CP2K_MASTER_REPO=/home/crm98/cp2k_master
OPT_BUILD=/home/crm98/cp2k-benchmarks/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_build.sh
OPT_REBUILD=/home/crm98/cp2k-benchmarks/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_rebuild.sh
MASTER_BUILD=/home/crm98/cp2k-benchmarks/cp2k_master/CSD3_build_scripts/cp2k_CSD3_master_build.sh

build_branch "master"                          "$CP2K_MASTER_REPO" ""                                "$MASTER_BUILD"
build_branch "feature-nnp-native-spline"       "$CP2K_OPT_REPO"    "feature/nnp-native-spline"       "$OPT_REBUILD"
build_branch "feature-nnp-native-spline-omp"   "$CP2K_OPT_REPO"    "feature/nnp-native-spline-omp"   "$OPT_REBUILD"

cp /home/crm98/cp2k_master/tools/toolchain/install/setup "$BIN_ROOT/setup"

# Matching md5sums would mean a build was silently skipped or LD_LIBRARY_PATH
# crossed branches.
echo "=== BINARY VERIFICATION ==="
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
md5sum "$BIN_ROOT/master/lib"/libcp2k.so* \
       "$BIN_ROOT/feature-nnp-native-spline/lib"/libcp2k.so* \
       "$BIN_ROOT/feature-nnp-native-spline-omp/lib"/libcp2k.so*
echo ""

cd /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/nnp_acs_pub/

echo "=== FIGURE 2a SCALING ==="
./run_nnp_acs_pub_scaling.sh master
./run_nnp_acs_pub_scaling.sh feature-nnp-native-spline
./run_nnp_acs_pub_scaling.sh feature-nnp-native-spline-omp

echo ""
echo "=== FIGURE 2b RDF TRAJECTORIES ==="
./run_nnp_acs_pub_rdf.sh master
./run_nnp_acs_pub_rdf.sh feature-nnp-native-spline
./run_nnp_acs_pub_rdf.sh feature-nnp-native-spline-omp

echo ""
echo "=== COMPUTING RDFs ==="
SCRATCH_RESULTS=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/paper_fig2
for branch in master feature-nnp-native-spline feature-nnp-native-spline-omp; do
   traj=$(ls -t "$SCRATCH_RESULTS"/fig2b_rdf/"$branch"/H2O-64_NNP_MD-pos-1.xyz 2>/dev/null | head -1)
   if [[ -f "$traj" ]]; then
      echo "RDF: $branch"
      python3 compute_rdf.py "$traj" \
              --box 12.42 \
              --out "$SCRATCH_RESULTS/fig2b_rdf/${branch}/rdf.csv"
   else
      echo "RDF: $branch - trajectory missing, skipping"
   fi
done

HOME_RESULTS=/home/crm98/cp2k-benchmarks/results/paper_fig2
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo ""
echo "CSVs synced to $HOME_RESULTS"
echo "Full trajectories remain at $SCRATCH_RESULTS"
