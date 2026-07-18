#!/usr/bin/env bash
set -euo pipefail
# One (box, molality, seed) conductivity run for the Madrid-2019 concentration
# series. Env: L (box side, A), MOLALITY (mol/kg), SEED. Produces
# <tag>.onsager (self + collective multi-origin MSDs) and <tag>.rdf, then
# deletes the position dump (disk at 97% -- dumps must be transient).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1

LD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/data/cerberus1/crm98/pyenv/bin/python3
LMP=/lsc/opt/lammps-2025/bin/lmp
ROOT=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/conductivity
RANKS="${RANKS:-24}"
EQUIL_PS="${EQUIL_PS:-200}"
PROD_PS="${PROD_PS:-2000}"
L="${L:?set L}"
MOLALITY="${MOLALITY:?set MOLALITY}"
SEED="${SEED:?set SEED}"

tag="L${L}_m${MOLALITY}_s${SEED}"
W="$ROOT/m${MOLALITY}/L${L}_s${SEED}"
mkdir -p "$W"; cd "$W"

if [ -f "$W/${tag}.onsager" ]; then echo "$tag: done, skipping"; exit 0; fi

"$PY" "$LD/make_lammps_data.py" --side "$L" --target-molality "$MOLALITY" \
  --out "$W/data.nacl"

rc=0
mpirun -n "$RANKS" "$LMP" -in "$LD/in.madrid_conduct" \
  -var data data.nacl -var settings "$LD/madrid.settings" \
  -var equil_ps "$EQUIL_PS" -var prod_ps "$PROD_PS" \
  -var seed "$SEED" -var tag "$tag" > "$W/prod.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || ! grep -q "DONE $tag" "$W/prod.log"; then
  echo "$tag FAILED (mpirun rc=$rc). Last lines of $W/prod.log:" >&2
  tail -20 "$W/prod.log" >&2
  exit 1
fi

"$PY" "$LD/onsager_multiorigin.py" --dump "$W/${tag}.pos.lammpstrj" \
  --dt-ps 1.0 --out "$W/${tag}.onsager"
rm -f "$W/${tag}.pos.lammpstrj"
echo "CONDUCT COMPLETE: $tag"
