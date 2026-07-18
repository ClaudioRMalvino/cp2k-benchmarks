#!/usr/bin/env bash
#SBATCH --job-name=conc_equil
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/conc_equil_%A_%a.out

# Concentration series (D vs molality) at FIXED box size cell211.
# One array task = one NaCl loading. N_PAIRS is pairs per 62-water base cell:
#   np0 -> pure water (0 mol/kg, D(c=0) reference)
#   np2 -> 1.85 mol/kg,  np3 -> 2.87 mol/kg,  np4 -> 3.97 mol/kg
# np1 (0.90 mol/kg) is already done by the finite-size run - reuse it.
# Idempotent: finished loadings are skipped; walltime-killed ones resume.
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

NP_LIST=(0 2 3 4)
export N_PAIRS="${NP_LIST[$SLURM_ARRAY_TASK_ID]}"
export CELLS="211"
export TOTAL_RANKS="$SLURM_NTASKS"
export MODEL="${MODEL:-revPBE-D3}"

exec /home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/run_equil.sh
