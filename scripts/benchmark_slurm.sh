#!/usr/bin/env bash
#SBATCH --job-name=Master_NNP_Scaling
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:15:00
#SBATCH --output=/home/raid/crm98/cp2k_benchmarks/logs/Master_Scaling_%j.out
#SBATCH --mail-type=ALL

mkdir -p /home/raid/crm98/cp2k_benchmarks/logs/

# Rebuilds and recompiles CP2k with the latest changes to the optimized branch
cd /home/raid/crm98/cp2k_benchmarks/cp2k_optimized/
./rebuild_cp2k.sh 

# Copies over the necessary binaries and tool setup for sourcing in the benchmark scripts
cp /local/data/public/crm98/original_cp2k/install/bin/cp2k.psmp ~/cp2k_binaries/phy-cerberus/cp2k_feature_verlet_cells.psmp
cp /local/data/public/crm98/original_cp2k/tools/toolchain/install/setup ~/cp2k_binaries/phy-cerberus/setup

cd /home/raid/crm98/cp2k_benchmarks/scripts/

echo "=== SIZE SCALING ==="
./run_nnp_size_scaling_slurm.sh master
./run_nnp_size_scaling_slurm.sh optimized

echo ""
echo "=== CORE SCALING ==="
./run_nnp_core_scaling_slurm.sh master
./run_nnp_core_scaling_slurm.sh optimized