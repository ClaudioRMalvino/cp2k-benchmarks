#!/usr/bin/env bash
set -euo pipefail
# Madrid-2019 finite-size sweep: cubic NaCl(aq) boxes at 1.0 mol/kg over a range
# of sizes, for the Yeh-Hummer D-vs-1/L study. Classical MD is cheap, so we use
# big boxes (into the linear YH regime, matching Dhruv's 21-40 A) and long runs.
# Cubic sides chosen to overlap the CP2K anchors (24.84 = cube2, 37.26 = cube3).

# lmp links OpenBLAS, libgomp and threaded FFTW, each of which defaults to one
# thread per core. Unset, every MPI rank spawns $(nproc) threads and the node
# dies allocating their stacks before step 0. One thread per rank; MPI does the
# parallelism.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1

LD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/data/cerberus1/crm98/pyenv/bin/python3
LMP=/lsc/opt/lammps-2025/bin/lmp
ROOT=/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/finite_size
RANKS="${RANKS:-8}"
EQUIL_PS="${EQUIL_PS:-200}"
PROD_PS="${PROD_PS:-1500}"
SIZES="${SIZES:-24.84 31.05 37.26 43.47}"
MOLALITY="${MOLALITY:-1.0}"

mkdir -p "$ROOT"
for L in $SIZES; do
  W="$ROOT/L${L}"; mkdir -p "$W"; cd "$W"
  tag="L${L}"
  if grep -q "DONE $tag" "$W/prod.log" 2>/dev/null; then
    echo "$tag: done, skipping"; continue
  fi
  "$PY" "$LD/make_lammps_data.py" --side "$L" --target-molality "$MOLALITY" \
    --out "$W/data.nacl"
  echo "=== $tag: $EQUIL_PS ps NPT + $PROD_PS ps NVT on $RANKS ranks ==="
  rc=0
  mpirun -n "$RANKS" "$LMP" -in "$LD/in.madrid_nacl" \
    -var data data.nacl -var settings "$LD/madrid.settings" \
    -var equil_ps "$EQUIL_PS" -var prod_ps "$PROD_PS" \
    -var tag "$tag" > "$W/prod.log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ] || ! grep -q "DONE $tag" "$W/prod.log"; then
    echo "$tag FAILED (mpirun rc=$rc). Last lines of $W/prod.log:" >&2
    tail -20 "$W/prod.log" >&2
    exit 1
  fi
done
echo "SWEEP COMPLETE: sizes=$SIZES"
