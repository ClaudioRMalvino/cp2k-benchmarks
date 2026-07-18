#!/usr/bin/env bash
#SBATCH --job-name=mp2cube_pr
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --array=0-4
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/mp2cube_prod_%A_%a.out

# MP2 cubic anchor production: one array task per NVE segment (5 segments).
# Each segment (~9 h at 1500 atoms) exceeds the 6 h wall, so a task will be
# killed once and must be resubmitted; it resumes from its 5 ps checkpoint
# (RESTART_COUNTERS fix means it stops at the right total, no 2x overrun).
# Chain after equil:  sbatch --dependency=afterok:<equil id> this_script
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

export MODEL=MP2
export CELLS="cube2"
export SEGMENTS="$(( SLURM_ARRAY_TASK_ID + 1 ))"
export TOTAL_RANKS="$SLURM_NTASKS"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/run_production.sh
