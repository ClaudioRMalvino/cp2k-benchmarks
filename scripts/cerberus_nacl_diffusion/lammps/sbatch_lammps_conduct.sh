#!/usr/bin/env bash
#SBATCH --job-name=madrid_conduct
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=06:00:00
#SBATCH --array=0-49
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/madrid_conduct_%A_%a.out

# Conductivity vs concentration: 2 box sizes x 5 molalities x 5 seeds = 50.
# Two sizes give a finite-size check on the collective (Onsager) quantities
# themselves (Celebi warns they are size-dependent; Gullbrekken found them
# empirically flat -- we test it for Madrid-2019).
# index i: size = i % 2, molality = (i / 2) % 5, seed = i / 10
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

SIZES=(31.05 43.47)
MOLS=(0.25 0.5 1.0 2.0 4.0)
SEEDS=(4321 8765 1597 2718 3141)

i=$SLURM_ARRAY_TASK_ID
export L="${SIZES[$(( i % 2 ))]}"
export MOLALITY="${MOLS[$(( (i / 2) % 5 ))]}"
export SEED="${SEEDS[$(( i / 10 ))]}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/lammps/run_lammps_conduct.sh
