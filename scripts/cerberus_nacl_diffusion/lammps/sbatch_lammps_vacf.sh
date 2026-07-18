#!/usr/bin/env bash
#SBATCH --job-name=madrid_vacf
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --array=0-3
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/madrid_vacf_%A_%a.out

# Green-Kubo VACF benchmark: one short NVE run per box, restarted from the
# finite-size sweep's equilibrated final.data. D from the VACF integral is then
# compared against D from the MSD slope — the estimator cross-check.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

LD=/home/raid/crm98/cp2k-benchmarks/scripts/cerberus_nacl_diffusion/lammps
PY=/data/cerberus1/crm98/pyenv/bin/python3
LMP=/lsc/opt/lammps-2025/bin/lmp
FS=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/finite_size
ROOT=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/vacf

SIZES=(24.84 31.05 37.26 43.47)
L="${SIZES[$SLURM_ARRAY_TASK_ID]}"
tag="L${L}_vacf"
W="$ROOT/L${L}"; mkdir -p "$W"; cd "$W"

if [ -f "$W/${tag}.vacf" ]; then echo "$tag: done, skipping"; exit 0; fi

rc=0
mpirun -n "$SLURM_NTASKS" "$LMP" -in "$LD/in.madrid_vacf" \
  -var data "$FS/L${L}/L${L}.final.data" -var settings "$LD/madrid.settings" \
  -var retherm_ps 20 -var prod_ps 40 -var tag "$tag" > "$W/vacf.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || ! grep -q "DONE $tag" "$W/vacf.log"; then
  echo "$tag FAILED (rc=$rc):" >&2; tail -20 "$W/vacf.log" >&2; exit 1
fi

"$PY" "$LD/vacf_analyze.py" --dump "$W/${tag}.vel.lammpstrj" \
  --dt-fs 4.0 --tmax-ps 5.0 --out "$W/${tag}.vacf"
rm -f "$W/${tag}.vel.lammpstrj"
echo "VACF COMPLETE: $tag"
