#!/usr/bin/env bash
#SBATCH --job-name=madrid_vacf2
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --array=0-19
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/madrid_vacf2_%A_%a.out

# Corrected Green-Kubo VACF campaign: restart from the REPLICA final.data files
# (correct NPT-equilibrated density), NOT the old July 10 sweep configs, which
# ran ~2.5% under-dense and taint the original vacf/ results. 4 boxes x 5 seeds
# gives seed error bars on the VACF D for the MSD-vs-VACF estimator cross-check.
# index i: box = i % 4, seed = i / 4
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

LD=/home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/lammps
PY=/data/cerberus1/crm98/pyenv/bin/python3
LMP=/lsc/opt/lammps-2025/bin/lmp
REP=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/replicas
ROOT=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/vacf_replica

SIZES=(24.84 31.05 37.26 43.47)
SEEDS=(1597 2718 3141 4321 8765)
i=$SLURM_ARRAY_TASK_ID
L="${SIZES[$(( i % 4 ))]}"
SEED="${SEEDS[$(( i / 4 ))]}"

tag="L${L}_s${SEED}_vacf"
W="$ROOT/L${L}_s${SEED}"; mkdir -p "$W"; cd "$W"

if [ -f "$W/${tag}.vacf" ]; then echo "$tag: done, skipping"; exit 0; fi

rc=0
mpirun -n "$SLURM_NTASKS" "$LMP" -in "$LD/in.madrid_vacf" \
  -var data "$REP/L${L}_s${SEED}/L${L}_s${SEED}.final.data" \
  -var settings "$LD/madrid.settings" \
  -var seed "$SEED" -var retherm_ps 20 -var prod_ps 40 \
  -var tag "$tag" > "$W/vacf.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || ! grep -q "DONE $tag" "$W/vacf.log"; then
  echo "$tag FAILED (rc=$rc):" >&2; tail -20 "$W/vacf.log" >&2; exit 1
fi

"$PY" "$LD/vacf_analyze.py" --dump "$W/${tag}.vel.lammpstrj" \
  --dt-fs 4.0 --tmax-ps 5.0 --out "$W/${tag}.vacf"
rm -f "$W/${tag}.vel.lammpstrj"
echo "VACF COMPLETE: $tag"
