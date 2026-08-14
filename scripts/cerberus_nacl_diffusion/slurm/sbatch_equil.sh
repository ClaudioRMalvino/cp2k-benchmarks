#!/usr/bin/env bash
#SBATCH --job-name=nacl_equil
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/nacl_equil_%A_%a.out

# Submit from athena (CSC head node). csc-mphil = phy-cerberus4/5/6: 48 cores +
# 248 GB each, 6 h max; requesting all 48 cores gives exclusive access.
# One array task = one box size. Idempotent: finished cells skipped, killed ones
# resume from checkpoint (only cell222 should hit the 6 h limit; resubmit).
# Concentration series: sbatch --export=ALL,N_PAIRS=2 sbatch_equil.sh
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

CELL_LIST=(111 211 221 222)
export CELLS="${CELL_LIST[$SLURM_ARRAY_TASK_ID]}"
export TOTAL_RANKS="$SLURM_NTASKS"
export N_PAIRS="${N_PAIRS:-1}"
export MODEL="${MODEL:-revPBE-D3}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/run_equil.sh
