#!/usr/bin/env bash
#SBATCH --job-name=NNP_Benchmarking
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --nodelist=phy-cerberus4
#SBATCH --ntasks=36
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:30:00
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/NNP_Benchmarking_%j.out
#SBATCH --mail-type=ALL

mkdir -p /home/raid/crm98/cp2k-benchmarks/logs/

# Per-branch binary cache: each branch has its own bin + lib so the run scripts
# can point LD_LIBRARY_PATH at $BIN_ROOT/<branch>/lib without picking up another
# branch's libcp2k.so from the scratch install tree (which gets overwritten on
# the next rebuild). Must live on /local/data/public, not /home: libcp2k.so.*
# is hundreds of MB per branch and exceeds the home quota.
BIN_ROOT=/local/data/public/crm98/cp2k_binaries/phy-cerberus
mkdir -p "$BIN_ROOT/master/lib"
mkdir -p "$BIN_ROOT/feature-nnp-verlet-cells/lib"
mkdir -p "$BIN_ROOT/feature-nnp-native-spline/lib"
mkdir -p "$BIN_ROOT/feature-nnp-native-spline-omp/lib"

CP2K_REPO=/home/raid/crm98/cp2k

# Master is built from a dedicated upstream clone into 'cp2k-buildtree' so it
# doesn't collide with the feature builds, which use a different scratch tree.
cd /home/raid/crm98/cp2k-benchmarks/cp2k_master/cerberus_build_scripts/
./cp2k_cerberus_master_build.sh -j 32
MASTER_INSTALL=/local/data/public/crm98/cp2k-buildtree/install
cp     "$MASTER_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/master/cp2k.psmp"
# -P preserves the versioned symlinks instead of dereferencing them.
cp -P "$MASTER_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/master/lib/"

cp /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup \
   "$BIN_ROOT/setup"

# rebuild_cp2k.sh rsyncs /home/raid/crm98/cp2k -> /local/data/public/crm98/original_cp2k,
# so all feature branches share the same scratch source AND install/{bin,lib}.
git -C "$CP2K_REPO" stash
git -C "$CP2K_REPO" checkout feature/nnp-verlet-cells
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/cerberus_build_scripts/
CP2K_FORCE_CONFIGURE=1 ./cp2k_cerberus_opt_rebuild.sh
FEATURE_INSTALL=/local/data/public/crm98/original_cp2k/install
cp     "$FEATURE_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp"
cp -P "$FEATURE_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/feature-nnp-verlet-cells/lib/"
git -C "$CP2K_REPO" stash pop

git -C "$CP2K_REPO" stash
git -C "$CP2K_REPO" checkout feature/nnp-native-spline
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/cerberus_build_scripts/
CP2K_FORCE_CONFIGURE=1 ./cp2k_cerberus_opt_rebuild.sh
cp     "$FEATURE_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
cp -P "$FEATURE_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/feature-nnp-native-spline/lib/"
git -C "$CP2K_REPO" stash pop

git -C "$CP2K_REPO" stash
git -C "$CP2K_REPO" checkout feature/nnp-native-spline-omp
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/cerberus_build_scripts/
CP2K_FORCE_CONFIGURE=1 ./cp2k_cerberus_opt_rebuild.sh
cp     "$FEATURE_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
cp -P "$FEATURE_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/feature-nnp-native-spline-omp/lib/"
git -C "$CP2K_REPO" stash pop

# Matching hashes would mean a rebuild was skipped or the wrong source tree
# was used.
echo "=== BINARY VERIFICATION ==="
echo "-- cp2k.psmp md5sums --"
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo "-- libcp2k.so md5sums --"
md5sum "$BIN_ROOT/master/lib"/libcp2k.so \
       "$BIN_ROOT/feature-nnp-verlet-cells/lib"/libcp2k.so \
       "$BIN_ROOT/feature-nnp-native-spline/lib"/libcp2k.so \
       "$BIN_ROOT/feature-nnp-native-spline-omp/lib"/libcp2k.so
echo "-- LD_LIBRARY_PATH sanity (ldd on verlet-cells binary) --"
source "$BIN_ROOT/setup"
export LD_LIBRARY_PATH="$BIN_ROOT/feature-nnp-verlet-cells/lib:${LD_LIBRARY_PATH:-}"
ldd "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp" | grep libcp2k
echo ""

cd /home/raid/crm98/cp2k-benchmarks/scripts/

echo "=== SIZE SCALING ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh feature-nnp-verlet-cells
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline-omp

echo ""
echo "=== CORE SCALING ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh feature-nnp-verlet-cells
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline-omp

SCRATCH_RESULTS=/local/data/public/crm98/cp2k-benchmarks/results
HOME_RESULTS=/home/raid/crm98/cp2k-benchmarks/results
mkdir -p "$HOME_RESULTS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH_RESULTS/" "$HOME_RESULTS/"
echo "CSVs synced to $HOME_RESULTS (full results remain on $SCRATCH_RESULTS)"
