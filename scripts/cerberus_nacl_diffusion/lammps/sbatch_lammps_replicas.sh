#!/usr/bin/env bash
#SBATCH --job-name=madrid_rep
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=06:00:00
#SBATCH --array=0-19
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/madrid_rep_%A_%a.out

# 4 cubic boxes x 5 independent seeds, all at 1.0 mol/kg. Each task is one
# replica. Ion statistics scale with (n_ions x seeds x time origins); the seed
# spread gives the error bar on D(L). Safe to resubmit: finished replicas skip.
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

SIZES=(24.84 31.05 37.26 43.47)
SEEDS=(4321 8765 1597 2718 3141)

i=$SLURM_ARRAY_TASK_ID
export L="${SIZES[$(( i % 4 ))]}"
export SEED="${SEEDS[$(( i / 4 ))]}"
export RANKS="$SLURM_NTASKS"
export EQUIL_PS="${EQUIL_PS:-150}"
export PROD_PS="${PROD_PS:-2000}"
export MOLALITY="${MOLALITY:-1.0}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/lammps/run_lammps_replica.sh
