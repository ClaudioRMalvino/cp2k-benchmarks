#!/usr/bin/env bash
#SBATCH --job-name=NNP_Benchmarking
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --nodelist=phy-cerberus4
#SBATCH --ntasks=36
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/NNP_Benchmarking_%j.out
#SBATCH --mail-type=ALL

mkdir -p /home/raid/crm98/cp2k-benchmarks/logs/

# ---------------------------------------------------------------------------
# Per-branch binary cache.  Each branch gets its own bin + lib directory so
# we never accidentally run binary-A against libcp2k.so-of-B: the run scripts
# point LD_LIBRARY_PATH at $BIN_ROOT/<branch>/lib, not at the scratch install
# tree (which gets overwritten by the next rebuild).
#
# Lives on scratch (/local/data/public), NOT on /home — libcp2k.so.* alone is
# hundreds of MB, three copies will blow out a typical home quota.  Job 10330
# died with "Disk quota exceeded" mid-rebuild because BIN_ROOT was on /home.
# ---------------------------------------------------------------------------
BIN_ROOT=/local/data/public/crm98/cp2k_binaries/phy-cerberus
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
cp     "$MASTER_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/master/cp2k.psmp"
# Copy only the runtime shared library + its versioned symlinks.  -P preserves
# the symlinks instead of dereferencing them.  We deliberately skip install/lib's
# cmake/ and pkgconfig/ subdirs — those are build-time metadata for downstream
# CMake consumers, not required for LD_LIBRARY_PATH at runtime.
cp -P "$MASTER_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/master/lib/"

# Toolchain setup is shared across every build (same toolchain dir).  Save once.
cp /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup \
   "$BIN_ROOT/setup"

# ---------------------------------------------------------------------------
# 2) Feature branch: nnp-verlet-cells.  rebuild_cp2k.sh rsyncs
#    /home/raid/crm98/cp2k -> /local/data/public/crm98/original_cp2k, so
#    'feature' and 'native-spline' share the same scratch source AND the
#    same install/{bin,lib}.  We MUST snapshot bin+lib into $BIN_ROOT before
#    moving on to the next branch.
# ---------------------------------------------------------------------------
git -C "$CP2K_REPO" stash
git -C "$CP2K_REPO" checkout feature/nnp-verlet-cells
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/
CP2K_FORCE_CONFIGURE=1 ./rebuild_cp2k.sh
FEATURE_INSTALL=/local/data/public/crm98/original_cp2k/install
cp     "$FEATURE_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp"
cp -P "$FEATURE_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/feature-nnp-verlet-cells/lib/"
git -C "$CP2K_REPO" stash pop

# ---------------------------------------------------------------------------
# 3) Feature branch: nnp-native-spline.  Overwrites $FEATURE_INSTALL — that's
#    fine because verlet-cells is already snapshotted in $BIN_ROOT above.
# ---------------------------------------------------------------------------
git -C "$CP2K_REPO" stash
git -C "$CP2K_REPO" checkout feature/nnp-native-spline
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/
CP2K_FORCE_CONFIGURE=1 ./rebuild_cp2k.sh
cp     "$FEATURE_INSTALL/bin/cp2k.psmp"     "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
cp -P "$FEATURE_INSTALL/lib"/libcp2k.so*    "$BIN_ROOT/feature-nnp-native-spline/lib/"
git -C "$CP2K_REPO" stash pop

# ---------------------------------------------------------------------------
# Sanity check: confirm the three cp2k.psmp binaries and their libcp2k.so
# are genuinely distinct.  Matching hashes would mean a rebuild was skipped
# or the wrong source tree was used.
# ---------------------------------------------------------------------------
echo "=== BINARY VERIFICATION ==="
echo "-- cp2k.psmp md5sums --"
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
echo "-- libcp2k.so md5sums --"
md5sum "$BIN_ROOT/master/lib"/libcp2k.so \
       "$BIN_ROOT/feature-nnp-verlet-cells/lib"/libcp2k.so \
       "$BIN_ROOT/feature-nnp-native-spline/lib"/libcp2k.so
echo "-- LD_LIBRARY_PATH sanity (ldd on verlet-cells binary) --"
source "$BIN_ROOT/setup"
export LD_LIBRARY_PATH="$BIN_ROOT/feature-nnp-verlet-cells/lib:${LD_LIBRARY_PATH:-}"
ldd "$BIN_ROOT/feature-nnp-verlet-cells/cp2k.psmp" | grep libcp2k
echo ""

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


#----------------------------------------------------------------------------
# Plots the results of the above benchmarks and saves the graphs as PNG files.
#----------------------------------------------------------------------------
python3 /home/raid/crm98/cp2k-benchmarks/plots/plot_scaling.py