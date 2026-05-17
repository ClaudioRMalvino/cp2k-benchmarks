#!/usr/bin/env bash
#SBATCH -J figS4_N128
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_N128_%j.out
# Budget: 16 cores × 1.5 h = 24 CPU-hr  (well under the original 32-rank prod)

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/replicate_N128_$(date +%H%M%S)
/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/figS4/replicate_figS4_pipeline.sh \
    128 2 1 1 "$CP2K_EXE" "$RUNDIR"
