#!/usr/bin/env bash
#SBATCH -J nacl_mp2_smoke
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/nacl_mp2_smoke_%j.out

# Run this FIRST, once, before any production submission (SL2 = fast queue,
# tiny cost). Verifies on a real icelake compute node that:
#   1. the dhruv-cell-list (PR #5295) binary runs the MIXED force eval
#      (C-NNP committee of 8 + Fist SPME Coulomb baseline) on cube_n2,
#   2. the per-step stress file is NONZERO (Green-Kubo viscosity needs it),
#   3. the committee-energies file is written (check_run.py needs it),
#   4. prints t/step -> walltime projections for equil/production.
set -euo pipefail
SDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
source "$SDIR/env_csd3.sh"

ANCHOR_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor
MODEL_DIR="$ANCHOR_ROOT/models/MP2/final_model"
CUBE_DIR="$ANCHOR_ROOT/runs/MP2/cubic_1M"
RUNDIR="$ANCHOR_ROOT/runs/MP2/smoke"
STEPS=200
RANKS=$SLURM_NTASKS

echo "node: $(hostname)"; lscpu | grep -E 'Model name|Socket|Core\(s\)'; echo
mkdir -p "$RUNDIR"; rm -f "$RUNDIR"/smoke.out "$RUNDIR"/NaCl_MP2_smoke*

# cube2: 24.84 A cubic, 1500 atoms, GMAX 80
sed -e "s|__PROJECT__|NaCl_MP2_smoke|" \
    -e "s|__STEPS__|$STEPS|" \
    -e "s|__SNAPSTEPS__|100000|" \
    -e "s|__MODEL__|$MODEL_DIR|" \
    -e "s|__COORD__|$CUBE_DIR/cube_n2.xyz|" \
    -e "s|__ABCX__|24.8400|" -e "s|__ABCY__|24.8400|" -e "s|__ABCZ__|24.8400|" \
    -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
    -e "s|__GMAXX__|80|" -e "s|__GMAXY__|80|" -e "s|__GMAXZ__|80|" \
    -e "s|__NNPRANKS__|$(( RANKS - 1 ))|" \
    -e "s|__FISTRANKS__|1|" \
    "$SDIR/nacl_diffusion_template.inp" > "$RUNDIR/smoke.inp"

cd "$RUNDIR"
/usr/bin/time -p srun --ntasks="$RANKS" --cpus-per-task=1 --hint=nomultithread \
  "$BIN" -i smoke.inp -o smoke.out

fail=0
if grep -q "PROGRAM ENDED" smoke.out; then echo "OK  cp2k finished"; else
  echo "FAIL: no PROGRAM ENDED in smoke.out"; fail=1; fi

st=$(ls NaCl_MP2_smoke-1.stress 2>/dev/null || true)
if [ -n "$st" ] && awk 'NR>1 {for(i=3;i<=NF;i++) if ($i+0!=0) {found=1; exit}} END{exit !found}' "$st"
then echo "OK  stress tensor nonzero ($st)"; else
  echo "FAIL: stress file missing or all-zero (GK viscosity would break)"; fail=1; fi

ce=$(ls *nnp-energies-1.data 2>/dev/null | head -1 || true)
if [ -n "$ce" ]; then echo "OK  committee energies written ($ce)"; else
  echo "FAIL: no *nnp-energies-1.data (check_run.py sigma canary needs it)"; fail=1; fi

# t/step from the self-timing report -> projections
tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' smoke.out || true)
if [ -n "${tstep:-}" ]; then
  awk -v t="$tstep" -v s="$STEPS" 'BEGIN{
    p=t/s;
    printf "t/step = %.4f s  (%d ranks, cube2/1500 atoms)\n", p, '"$RANKS"';
    printf "  equil  (210k steps): %5.1f h\n", p*210000/3600;
    printf "  prod   (200k steps): %5.1f h /segment\n", p*200000/3600;
    printf "  cube3 estimate (x3.4 atoms): equil %5.1f h, prod %5.1f h/seg\n",
           3.4*p*210000/3600, 3.4*p*200000/3600 }'
fi

[ $fail -eq 0 ] && echo "SMOKE TEST PASSED" || { echo "SMOKE TEST FAILED"; exit 1; }
