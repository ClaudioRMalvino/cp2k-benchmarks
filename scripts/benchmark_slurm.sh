#!/usr/bin/env bash
#SBATCH --job-name=Master_NNP_Scaling
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:25:00
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/Master_Scaling_%j.out
#SBATCH --mail-type=ALL

mkdir -p /home/raid/crm98/cp2k-benchmarks/logs/

export OMP_NUM_THREADS=2

# --- Build upstream/master (clean build, always reconfigures via rm -rf build) ---
cd /home/raid/crm98/cp2k-benchmarks/cp2k_master/
./build_cp2k.sh
cp /local/data/public/crm98/cp2k-buildtree/install/bin/cp2k.psmp \
   ~/cp2k_binaries/phy-cerberus/cp2k_master.psmp

# --- Build feature branch (incremental, force reconfigure to pick up flag changes) ---
cd /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/
CP2K_FORCE_CONFIGURE=1 ./rebuild_cp2k.sh
cp /local/data/public/crm98/original_cp2k/install/bin/cp2k.psmp \
   ~/cp2k_binaries/phy-cerberus/cp2k_feature_verlet_cells.psmp
cp /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup \
   ~/cp2k_binaries/phy-cerberus/setup

cd /home/raid/crm98/cp2k-benchmarks/scripts/

echo "=== SIZE SCALING ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh optimized

echo ""
echo "=== CORE SCALING ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh optimized
