#!/usr/bin/env bash
#SBATCH --job-name=madrid_fs
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/madrid_fs_%A_%a.out

# Madrid-2019 finite-size sweep as a SLURM array: one cubic box per task, all
# at 1.0 mol/kg. Classical MD is cheap so each box finishes well under 6 h at
# 48 cores; array runs them in parallel. The driver skips any box already
# marked DONE (e.g. L24.84 completed earlier), so this is safe to resubmit.
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

SIZES_LIST=(24.84 31.05 37.26 43.47)
export SIZES="${SIZES_LIST[$SLURM_ARRAY_TASK_ID]}"
export RANKS="$SLURM_NTASKS"
export EQUIL_PS="${EQUIL_PS:-150}"
export PROD_PS="${PROD_PS:-2000}"
export MOLALITY="${MOLALITY:-1.0}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/lammps/run_lammps_sweep.sh
