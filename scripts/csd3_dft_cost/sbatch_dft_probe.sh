#!/usr/bin/env bash
#SBATCH -J dft_probe
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/dft_probe_%j.out

# Diagnostic probe for the 5064-atom (37.26 A) rung of the DFT ladder, which
# never starts its SCF. Evidence: 33128702 (76x1, 6 h 25), 33170265 (19x4, 1 h),
# 33170266 (19x4 + RS_GRID REPLICATED, 1 h) - zero SCF iterations each; gdb on
# 33128702 put every rank in calculate_rho_core -> transfer_rs2pw_distributed ->
# mp_waitany (rs->pw transfer of the core density, pre-SCF). Not memory (80 of
# 502 GB), not the FULL_ALL preconditioner (never reached); slab-vs-halo ruled out
# (three layouts + replicated grid stall identically; diffuse Gaussians live on
# the COARSE multigrids).
# This probe: srun under `timeout` so the report below always runs; PRINT_LEVEL
# MEDIUM prints the grid tables before the stall; two gdb backtraces (8/20 min)
# distinguish stuck from grinding. Accuracy settings (cutoff/EPS_*/XC/basis/D3)
# untouched - layout, FFT planning and print level change no computed number.
#   sbatch --ntasks=19 --cpus-per-task=4 sbatch_dft_probe.sh 19 4 base 2700
set -euo pipefail
RANKS=${1:?usage: sbatch_dft_probe.sh <ranks> <threads> <base|replicated> [run_seconds]}
THREADS=${2:?}
VARIANT=${3:-base}
RUNSECS=${4:-2700}

TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_dft_cost
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
ROOT=/rds/user/$USER/hpc-work/dft_cost_ladder
export BIN_LABEL=master
export OMP_NUM_THREADS=$THREADS
source "$ADIR/env_csd3.sh"
export OMP_NUM_THREADS=$THREADS

DATA=/home/crm98/cp2k_master/data
COORD="$ROOT/cube_n3_equil.xyz"
RUNDIR="$ROOT/probe_${RANKS}x${THREADS}_${VARIANT}"
proj="nacl_revpbe_d3_probe"

echo "node: $(hostname)  layout: ${RANKS} ranks x ${THREADS} threads = $((RANKS*THREADS)) cores  variant: $VARIANT  run_seconds: $RUNSECS"
[ -f "$COORD" ] || { echo "FATAL: missing $COORD"; exit 6; }
mkdir -p "$RUNDIR"; rm -f "$RUNDIR/$proj"* "$RUNDIR/probe.out"

sed -e "s|__PROJECT__|$proj|" \
    -e "s|__NREP__|1 1 1|g" \
    -e "s|__STEPS__|2|" \
    -e "s|__DATA__|$DATA|g" \
    -e "s|__ABC__|37.2600|g" \
    -e "s|__COORD__|$COORD|" \
    "$TDIR/dft_ladder.inp.template" > "$RUNDIR/probe.inp"

sed -i 's|^  FFTW_PLAN_TYPE MEASURE$|  FFTW_PLAN_TYPE ESTIMATE|' "$RUNDIR/probe.inp"
sed -i 's|^  PRINT_LEVEL LOW$|  PRINT_LEVEL MEDIUM|' "$RUNDIR/probe.inp"
grep -q "PRINT_LEVEL MEDIUM" "$RUNDIR/probe.inp" || { echo "FATAL: PRINT_LEVEL patch did not apply"; exit 9; }

if [ "$VARIANT" = replicated ]; then
  sed -i 's|^      CUTOFF 1200$|      CUTOFF 1200\n      \&RS_GRID\n        DISTRIBUTION_TYPE REPLICATED\n      \&END RS_GRID|' "$RUNDIR/probe.inp"
  grep -q "DISTRIBUTION_TYPE REPLICATED" "$RUNDIR/probe.inp" || { echo "FATAL: RS_GRID patch did not apply"; exit 7; }
fi

cd "$RUNDIR"
t0=$SECONDS
timeout "$RUNSECS" srun --ntasks="$RANKS" --cpus-per-task="$THREADS" --hint=nomultithread \
     "$BIN" -i probe.inp -o probe.out &
srun_pid=$!

# sample a live rank twice: same frame at both = stuck, different = grinding
sample() {
  local tag=$1
  local p
  p=$(pgrep -u "$USER" -x cp2k.psmp | head -1 || true)
  echo; echo "=== STACK SAMPLE ($tag, t+$((SECONDS - t0))s, pid=${p:-none}) ==="
  [ -n "${p:-}" ] || { echo "no cp2k.psmp rank found"; return; }
  timeout 120 gdb -p "$p" -batch -ex "bt 14" 2>/dev/null | grep -E "^#" | head -14 || echo "(gdb unavailable)"
  echo "--- probe.out size: $(stat -c %s probe.out 2>/dev/null || echo 0) bytes"
}
( sleep 480;  sample "8 min"  ) &
( sleep 1200; sample "20 min" ) &

wait "$srun_pid" || echo "(srun stopped: walled out or timed out - expected)"
wall=$((SECONDS - t0))

echo; echo "=== PROBE RESULT (${RANKS}x${THREADS} $VARIANT, wall ${wall}s) ==="
echo "SCF iterations printed : $(grep -cE 'OT (DIIS|SD) ' probe.out 2>/dev/null || echo 0)"
echo "SCF converged banners  : $(grep -c 'SCF run converged' probe.out 2>/dev/null || echo 0)"
echo; echo "--- grid geometry as CP2K actually built it ---"
grep -E "^ (GRID|MGRID|QS)\|" probe.out 2>/dev/null | head -30 || echo "(no grid tables printed)"
echo; echo "--- realspace grid distribution ---"
grep -iE "RS_GRID|distribut" probe.out 2>/dev/null | head -20 || true
if [ -f "$proj-1.ener" ]; then
  echo; echo "--- .ener (last col = UsedTime per step) ---"; cat "$proj-1.ener"
else
  echo; echo "no .ener - did not complete a single MD step"
fi
echo; echo "--- last 10 non-blank lines of probe.out ---"
grep -v '^ *$' probe.out 2>/dev/null | tail -10 || true
