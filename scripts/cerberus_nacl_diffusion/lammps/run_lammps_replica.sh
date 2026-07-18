#!/usr/bin/env bash
set -euo pipefail
# One (box size, seed) replica of the Madrid-2019 finite-size study.
# Seeds are independent samples: their spread across replicas IS the error bar
# on D(L), which a single trajectory cannot give you.

# lmp links OpenBLAS, libgomp and threaded FFTW, each defaulting to one thread
# per core. Unset, every MPI rank spawns $(nproc) threads and the node dies
# allocating stacks before step 0. One thread per rank; MPI does the work.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1

LD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/data/cerberus1/crm98/pyenv/bin/python3
LMP=/lsc/opt/lammps-2025/bin/lmp
ROOT=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/replicas
RANKS="${RANKS:-24}"
EQUIL_PS="${EQUIL_PS:-150}"
PROD_PS="${PROD_PS:-2000}"
MOLALITY="${MOLALITY:-1.0}"
L="${L:?set L (cubic side in A)}"
SEED="${SEED:?set SEED}"

tag="L${L}_s${SEED}"
W="$ROOT/$tag"; mkdir -p "$W"; cd "$W"
if grep -q "DONE $tag" "$W/prod.log" 2>/dev/null; then
  echo "$tag: done, skipping"; exit 0
fi

"$PY" "$LD/make_lammps_data.py" --side "$L" --target-molality "$MOLALITY" \
  --seed "$SEED" --out "$W/data.nacl"

echo "=== $tag: $EQUIL_PS ps NPT + $PROD_PS ps NVT on $RANKS ranks ==="
rc=0
mpirun -n "$RANKS" "$LMP" -in "$LD/in.madrid_nacl" \
  -var data data.nacl -var settings "$LD/madrid.settings" \
  -var equil_ps "$EQUIL_PS" -var prod_ps "$PROD_PS" \
  -var seed "$SEED" -var tag "$tag" > "$W/prod.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || ! grep -q "DONE $tag" "$W/prod.log"; then
  echo "$tag FAILED (mpirun rc=$rc). Last lines of $W/prod.log:" >&2
  tail -20 "$W/prod.log" >&2
  exit 1
fi

# Collapse the position dump into a multi-origin MSD immediately, then drop the
# dump: it is ~300 MB-1 GB per replica and we only ever need the MSD from it.
echo "=== $tag: multi-origin MSD ==="
"$PY" "$LD/msd_multiorigin.py" --dump "$W/${tag}.pos.lammpstrj" \
  --dt-ps 1.0 --out "$W/${tag}.msdmo"
rm -f "$W/${tag}.pos.lammpstrj"
echo "REPLICA COMPLETE: $tag"
