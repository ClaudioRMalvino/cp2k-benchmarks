#!/usr/bin/env bash
#SBATCH --job-name=conc_prod
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --array=0-19
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/conc_prod_%A_%a.out

# Concentration series production at FIXED box cell211.
# array index = np_slot*5 + (seg-1), np_slot -> NP_LIST=(0 2 3 4).
# Chain after equil:  sbatch --dependency=afterok:<conc_equil id> this_script
# cell211 segments are ~2 h each (< 6 h) so no resume juggling expected.
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

NP_LIST=(0 2 3 4)
idx="$SLURM_ARRAY_TASK_ID"
export N_PAIRS="${NP_LIST[$(( idx / 5 ))]}"
export SEGMENTS="$(( idx % 5 + 1 ))"
export CELLS="211"
export TOTAL_RANKS="$SLURM_NTASKS"
export MODEL="${MODEL:-revPBE-D3}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/run_production.sh
