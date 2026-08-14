#!/usr/bin/env bash
#SBATCH -J rpa_smoke
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=00:40:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_smoke_%j.out

# Transport campaign, run 0: smoke of the campaign engine (Dhruv's PR #5295
# dhruv-cell-list binary) with the RPA committee (supervisor decision 2026-07-30:
# better dynamics than MP2, identical cost) on the 1 mol/kg 37.26 A cube3 deck.
# SL2 fast queue; gates every other campaign job (afterok). Checks: MIXED deck at
# 5064 atoms, NONZERO stress (GK viscosity needs it), committee/ener/traj streams,
# check_run.py, s/step projections, step-0 E_pot vs master_ref (same coordinates).
set -euo pipefail
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport
export BIN_LABEL=dhruv-cell-list
source "$ADIR/env_csd3.sh"

MODEL=RPA
ANCHOR_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor
MODEL_DIR="$ANCHOR_ROOT/models/$MODEL/final_model"
COORD="$ANCHOR_ROOT/runs/$MODEL/cubic_1M/cube_n3.xyz"
RUNDIR="$ANCHOR_ROOT/runs/$MODEL/smoke"
STEPS=200
RANKS=$SLURM_NTASKS
FIST_RANKS=8   # cube3 split-sweep winner from the anchor campaign
proj="NaCl_${MODEL}_smoke"

# provenance without exec'ing the binary (direct exec aborts in MPI_Init; launch only via srun)
echo "node: $(hostname)  binary: $BIN  model: $MODEL"
cat "$BIN_ROOT/$BIN_LABEL/PROVENANCE.txt" 2>/dev/null || true; echo
mkdir -p "$RUNDIR"; rm -f "$RUNDIR"/smoke.out "$RUNDIR/$proj"*

# energy print every 10 steps so check_run.py has enough records at 200 steps
sed -e "s|__PROJECT__|$proj|" \
    -e "s|__STEPS__|$STEPS|" \
    -e "s|__SNAPSTEPS__|100000|" \
    -e "s|__MODEL__|$MODEL_DIR|" \
    -e "s|__COORD__|$COORD|" \
    -e "s|__ABCX__|37.2600|" -e "s|__ABCY__|37.2600|" -e "s|__ABCZ__|37.2600|" \
    -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
    -e "s|__GMAXX__|120|" -e "s|__GMAXY__|120|" -e "s|__GMAXZ__|120|" \
    -e "s|__NNPRANKS__|$(( RANKS - FIST_RANKS ))|" \
    -e "s|__FISTRANKS__|$FIST_RANKS|" \
    -e "s|MD 100$|MD 10|" \
    "$ADIR/nacl_diffusion_template.inp" > "$RUNDIR/smoke.inp"

cd "$RUNDIR"
t0=$SECONDS
srun --ntasks="$RANKS" --cpus-per-task=1 --hint=nomultithread \
  "$BIN" -i smoke.inp -o smoke.out
wall=$(( SECONDS - t0 ))

fail=0
if grep -q "PROGRAM ENDED" smoke.out; then echo "OK  cp2k finished"; else
  echo "FAIL: no PROGRAM ENDED in smoke.out"; fail=1; fi

st="$proj-1.stress"
if [ -f "$st" ] && awk 'NR>1 {for(i=3;i<=NF;i++) if ($i+0!=0) {found=1; exit}} END{exit !found}' "$st"
then echo "OK  stress tensor nonzero ($st)"; else
  echo "FAIL: stress file missing or all-zero (GK viscosity would break)"; fail=1; fi

ce=$(ls *nnp-energies-1.data 2>/dev/null | head -1 || true)
if [ -n "$ce" ]; then echo "OK  committee energies written ($ce)"; else
  echo "FAIL: no *nnp-energies-1.data"; fail=1; fi

if [ -f "$proj-pos-1.xyz" ] && [ -f "$proj-1.ener" ]; then
  echo "OK  trajectory + ener streams present"; else
  echo "FAIL: missing pos/ener stream"; fail=1; fi

# step-0 potential energy (col 5 of .ener) for the cross-binary check
e0=$(awk 'NR==2 {print $5}' "$proj-1.ener")
echo "step-0 potential energy [Ha]: $e0   (compare against master_ref)"

tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' smoke.out || true)
if [ -n "${tstep:-}" ]; then
  ps=$(awk -v t="$tstep" -v s="$STEPS" 'BEGIN{printf "%.4f", t/s}')
  echo "t/step = $ps s  ($RANKS ranks, cube3/5064 atoms, $BIN_LABEL, $MODEL)"
  awk -v p="$ps" 'BEGIN{
    printf "  equil (210k steps): %5.1f h;  prod segment (320k steps / 160 ps): %5.1f h\n",
           p*210000/3600, p*320000/3600 }'
  bash "$TDIR/log_timing.sh" smoke cubic_1M "$BIN_LABEL" "$RANKS" "$STEPS" "$wall" "$ps" "cube3 5064 atoms RPA"
fi

set +u; source ~/.fortran_env/bin/activate; set -u
python "$ADIR/check_run.py" "$RUNDIR" || fail=1

[ $fail -eq 0 ] && echo "SMOKE TEST PASSED" || { echo "SMOKE TEST FAILED"; exit 1; }
