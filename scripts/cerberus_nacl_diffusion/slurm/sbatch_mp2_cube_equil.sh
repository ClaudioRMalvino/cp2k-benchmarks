#!/usr/bin/env bash
#SBATCH --job-name=mp2cube_eq
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/mp2cube_equil_%j.out

# MP2 cubic 1.0 mol/kg anchor, 24.84 A box (cube2 = 1500 atoms, 9 NaCl pairs).
# ~9 h equilibration at 48 cores > 6 h wall, so it will hit the limit once:
# just resubmit this script and it continues from the last checkpoint.
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

export MODEL=MP2
export CELLS="cube2"
export TOTAL_RANKS="$SLURM_NTASKS"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/run_equil.sh
