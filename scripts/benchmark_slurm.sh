#!/usr/bin/env bash
#SBATCH --job-name=Master_NNP_Scaling
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=36
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/Master_Scaling_%j.out
#SBATCH --mail-type=ALL

mkdir -p /home/raid/crm98/cp2k-benchmarks/logs/

# ---------------------------------------------------------------------------
# Per-branch binary cache.  Each branch gets its own bin + lib directory so
# we never accidentally run binary-A against libcp2k.so-of-B: the run scripts
# point LD_LIBRARY_PATH at $BIN_ROOT/<branch>/lib, not at the scratch install
# tree (which gets overwritten by the next rebuild).
# ---------------------------------------------------------------------------
BIN_ROOT=~/cp2k_binaries/phy-cerberus
mkdir -p "$BIN_ROOT/master/lib"
mkdir -p "$BIN_ROOT/feature-nnp-verlet-cells/lib"
mkdir -p "$BIN_ROOT/feature-nnp-native-spline/lib"

CP2K_REPO=/home/raid/crm98/cp2k

# ---------------------------------------------------------------------------
# 1) Master — built from the dedicated upstream clone into 'cp2k-buildtree'.
#    No collision risk with the feature builds (different scratch tree).
# ---------------------------------------------------------------------------
cd /home/raid/crm98/cp2k-benchmarks/cp2k_master/
./build_cp2k.sh
MASTER_INSTALL=/local/data/public/crm98/cp2k-buildtree/install
cp "$MASTER_INSTALL/bin/cp2k.psmp"   "$BIN_ROOT/master/cp2k.psmp"
cp -r "$MASTER_INSTALL/lib/."        "$BIN_ROOT/master/lib/"

# Toolchain setup is shared across every build (same toolchain dir).  Save once.
cp /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup \
   "$BIN_ROOT/setup"

# ---------------------------------------------------------------------------
# 2) Feature branch: nnp-verlet-cells.  IMPORTANT: rebuild_cp2k.sh rsyncs
#    /home/raid/crm98/cp2k -> /local/data/public/crm98/original_cp2k, so
#    'feature' and 'native-spline' share the same scratch source AND the
#    same install/{bin,lib}.  We MUST snapshot bin+lib into $BIN_ROOT before
#    moving on to the next branch.
# ---------------------------------------------------------------------------
git -C "$CP2K_REPO" checkout feature/nnp-verlet-cells
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/
CP2K_FORCE_CONFIGURE=1 ./rebuild_cp2k.sh
FEATURE_INSTALL=/local/data/public/crm98/original_cp2k/install
cp "$FEATURE_INSTALL/bin/cp2k.psmp" "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp"
cp -r "$FEATURE_INSTALL/lib/."      "$BIN_ROOT/feature-nnp-verlet-cells/lib/"

# ---------------------------------------------------------------------------
# 3) Feature branch: nnp-native-spline.  Overwrites $FEATURE_INSTALL — that's
#    fine because verlet-cells is already snapshotted in $BIN_ROOT above.
# ---------------------------------------------------------------------------
git -C "$CP2K_REPO" checkout feature/nnp-native-spline
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/
CP2K_FORCE_CONFIGURE=1 ./rebuild_cp2k.sh
cp "$FEATURE_INSTALL/bin/cp2k.psmp" "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
cp -r "$FEATURE_INSTALL/lib/."      "$BIN_ROOT/feature-nnp-native-spline/lib/"

# ---------------------------------------------------------------------------
# Benchmarks.  The run scripts read from $BIN_ROOT/<branch>/{cp2k.psmp,lib},
# so the home-repo's currently-checked-out branch is irrelevant from here on.
# ---------------------------------------------------------------------------
cd /home/raid/crm98/cp2k-benchmarks/scripts/

echo "=== SIZE SCALING ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh feature-nnp-verlet-cells
./run_nnp_size_scaling_slurm.sh feature-nnp-native-spline

echo ""
echo "=== CORE SCALING ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh feature-nnp-verlet-cells
./run_nnp_core_scaling_slurm.sh feature-nnp-native-spline
