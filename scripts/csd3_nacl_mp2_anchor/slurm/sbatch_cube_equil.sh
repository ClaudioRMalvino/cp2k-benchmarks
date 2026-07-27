#!/usr/bin/env bash
#SBATCH -J nacl_mp2_equil
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/nacl_mp2_equil_%j.out

# NVT equilibration of one cubic anchor cell (210k steps = 30 ps settle +
# 5 snapshots 15 ps apart). Cell as argument:
#     sbatch sbatch_cube_equil.sh cube2      # 24.84 A, 1500 atoms
#     sbatch sbatch_cube_equil.sh cube3      # 37.26 A, 5064 atoms
# Checkpoint-safe: if the 12 h wall kills it, resubmit the same command and
# it continues from the last 15 ps checkpoint (finished work is skipped).
set -euo pipefail

export MODEL=MP2
export CELLS="${1:-cube2}"
export TOTAL_RANKS="$SLURM_NTASKS"
# optional 2nd arg = Fist ranks of the MIXED split (from the split-sweep job)
[ -n "${2:-}" ] && export FIST_RANKS_OVERRIDE="$2"

exec /home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor/run_equil.sh
