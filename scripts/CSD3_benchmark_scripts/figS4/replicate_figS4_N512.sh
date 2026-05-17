#!/usr/bin/env bash
#SBATCH -J figS4_N512
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=05:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_N512_%j.out
# Budget: 64 cores × 5 h = 320 CPU-hr

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/replicate_N512_$(date +%H%M%S)
/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/figS4/replicate_figS4_pipeline.sh \
    512 2 2 2 "$CP2K_EXE" "$RUNDIR"
