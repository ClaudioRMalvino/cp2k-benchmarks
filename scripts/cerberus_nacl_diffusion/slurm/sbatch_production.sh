#!/usr/bin/env bash
#SBATCH --job-name=nacl_prod
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --array=0-19
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/nacl_prod_%A_%a.out

# One array task = one (box size, NVE segment): index = cell*5 + (seg-1).
# Chain after equil: sbatch --dependency=afterok:<equil id> sbatch_production.sh.
# Idempotent: finished segments skipped; killed ones resume from 5 ps checkpoints
# (cell222 segments need ~9 h = two 6 h passes).
# Concentration series: sbatch --export=ALL,N_PAIRS=2 sbatch_production.sh
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

CELL_LIST=(111 211 221 222)
idx="$SLURM_ARRAY_TASK_ID"
export CELLS="${CELL_LIST[$(( idx / 5 ))]}"
export SEGMENTS="$(( idx % 5 + 1 ))"
export TOTAL_RANKS="$SLURM_NTASKS"
export N_PAIRS="${N_PAIRS:-1}"
export MODEL="${MODEL:-revPBE-D3}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/run_production.sh
