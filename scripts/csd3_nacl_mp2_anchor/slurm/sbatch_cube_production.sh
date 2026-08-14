#!/usr/bin/env bash
#SBATCH -J nacl_mp2_prod
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/nacl_mp2_prod_%A_%a.out

# NVE production: one array task per 100 ps segment (5 segments = replicas for
# error bars), each from its stage-1 snapshot.
#     eq=$(sbatch --parsable sbatch_cube_equil.sh cube2)
#     sbatch --dependency=afterok:$eq sbatch_cube_production.sh cube2
# Checkpoint-safe: killed tasks resume from 5 ps checkpoints (RESTART_COUNTERS
# fix, no 2x overrun); completed segments are skipped on resubmit.
set -euo pipefail

export MODEL=MP2
export CELLS="${1:-cube2}"
export SEGMENTS="$(( SLURM_ARRAY_TASK_ID + 1 ))"
export TOTAL_RANKS="$SLURM_NTASKS"
# optional 2nd arg = Fist ranks of the MIXED split (from the split-sweep job)
[ -n "${2:-}" ] && export FIST_RANKS_OVERRIDE="$2"

exec /home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor/run_production.sh
