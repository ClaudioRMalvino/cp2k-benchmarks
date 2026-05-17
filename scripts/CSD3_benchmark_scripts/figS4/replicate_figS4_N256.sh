#!/usr/bin/env bash
#SBATCH -J figS4_N256
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_N256_%j.out
# Budget: 32 cores × 3 h = 96 CPU-hr

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/replicate_N256_$(date +%H%M%S)
/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/figS4/replicate_figS4_pipeline.sh \
    256 2 2 1 "$CP2K_EXE" "$RUNDIR"
