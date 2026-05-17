#!/usr/bin/env bash
# Validation of the Verlet-skin fix (default_skin 1.0 → 3.0, cap 0.1*max_cut
# → 0.5*max_cut in nnp_acsf.F::nnp_prepare_cell_list_cache).
#
# Runs the SAME 200-step NVT from the SAME GEO_OPT-minimised restart used in
# the earlier master-vs-feature comparison, but with the new skinfix binary.
# Pass criterion (correctness): no single-step Pot jumps > 0.10 ha, ConsQty
# drift |end-0| < 0.01 ha — same bar master cleared (0.0011 ha drift, 0 jumps).
#
# Also reports MD-loop time per step so we can immediately see the perf hit
# vs the original cached feature-nnp-native-spline binary (which timed
# ~0.015 s/step at 32 ranks on N=128 in your existing benchmarks).

#SBATCH -J figS4_diag_skinfix
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/diag_skinfix_%j.out
# Budget: 4 cores × 10 min = 0.67 CPU-hr

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

CP2K_FIX=$BIN_ROOT/feature-nnp-native-spline-omp-skinfix/cp2k.psmp
[[ ! -x "$CP2K_FIX" ]] && { echo "*** skinfix binary missing: $CP2K_FIX"; exit 1; }
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline-omp-skinfix/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

FEATURE_DIR=/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N128_223052
RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_skinfix_$(date +%H%M%S)
mkdir -p "$RUNDIR"
ln -sfn /home/crm98/cp2k-benchmarks/potentials "$RUNDIR/NNP"
cp "$FEATURE_DIR/relax_N128-1.restart" "$RUNDIR/relax_input.restart"

python3 <<PYEOF
import re
src = open("$FEATURE_DIR/02_nvt.inp").read()
src = re.sub(r'PROJECT\s+\S+',  'PROJECT  skinfix_N128', src)
src = re.sub(r'STEPS\s+\d+',    'STEPS 200', src, count=1)
src = re.sub(r'RESTART_FILE_NAME\s+\S+',
             'RESTART_FILE_NAME relax_input.restart', src)
open("$RUNDIR/run.inp", "w").write(src)
PYEOF

cd "$RUNDIR"
echo "=== running NVT 200 steps with SKINFIX binary (omp branch, default_skin=3.0) ==="
t0=$SECONDS
srun --ntasks=4 --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_FIX" -i run.inp > cp2k.out 2>&1 || true
wall=$(( SECONDS - t0 ))
echo "wall: ${wall} s"

ener=$(ls *.ener 2>/dev/null | head -1)
[[ -z "$ener" ]] && { echo "*** no .ener — skinfix run failed:"; tail -30 cp2k.out; exit 1; }

echo ""
echo "=== trajectory (steps 0/4/8/17/29/60/100/150/200) ==="
awk 'NR>1 && ($1==0||$1==4||$1==8||$1==17||$1==29||$1==60||$1==100||$1==150||$1==200) {
  printf "  step=%-3d  T=%7.1f  Pot=%10.5f  ConsQty=%10.5f\n", $1, $4, $5, $6
}' "$ener"

echo ""
echo "=== VERDICT ==="
python3 <<PYEOF
import re

def trace(path):
    out = []
    for ln in open(path):
        if ln.startswith('#') or not ln.strip(): continue
        p = ln.split()
        if len(p) < 6: continue
        try: out.append((int(p[0]), float(p[3]), float(p[4]), float(p[5])))
        except: pass
    return out

r = trace("$ener")
big_pot_jumps = sum(1 for i in range(1, len(r)) if abs(r[i][2] - r[i-1][2]) > 0.10)
big_cq_jumps  = sum(1 for i in range(1, len(r)) if abs(r[i][3] - r[i-1][3]) > 0.10)
maxT = max(x[1] for x in r)
finalT = r[-1][1]
cq_drift = abs(r[-1][3] - r[0][3])
print(f"  steps                          : {len(r)}")
print(f"  T (final / max)                : {finalT:.1f} K / {maxT:.1f} K")
print(f"  ConsQty drift (|end-0|)        : {cq_drift:.4f} hartree")
print(f"  single-step Pot jumps > 0.10   : {big_pot_jumps}")
print(f"  single-step ConsQty jumps>0.10 : {big_cq_jumps}")
print()
print(f"  reference (original feature)   : 3+ jumps, ConsQty drift > 0.4 ha in 30 steps")
print(f"  reference (master, smooth)     : 0 jumps, ConsQty drift = 0.0011 ha")
print()
if big_pot_jumps == 0 and cq_drift < 0.01:
    print("  ==>  SKINFIX WORKS.  Matches master's stability.  Bug confirmed in")
    print("       Verlet-skin sizing; fix is the 2-line edit in nnp_acsf.F.")
elif big_pot_jumps < 3:
    print(f"  ==>  PARTIAL.  Fewer jumps than original ({big_pot_jumps}) but not zero.")
    print(f"       Skin may need to be even larger, or another mechanism is in play.")
else:
    print("  ==>  SKINFIX DID NOT HELP.  Jumps still present at same frequency.")
    print("       The bug isn't (only) about skin sizing.")
PYEOF

# Performance: extract MD-loop time-per-step from the printout
echo ""
echo "=== PERFORMANCE: time/step (skinfix binary, 4 ranks, 200 NVT steps) ==="
md_loop=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' cp2k.out | tail -1)
if [[ -n "${md_loop:-}" ]]; then
    tps=$(awk -v t="$md_loop" 'BEGIN{printf "%.5f", t/200}')
    echo "  qs_mol_dyn_low total : ${md_loop} s"
    echo "  time/step (this run) : ${tps} s/step  (at 4 ranks, N=128, on icelake)"
    echo "  for comparison: your prior benchmarks at 32 ranks gave ~0.015 s/step for N=128"
    echo "  (NOT directly comparable — different rank counts; full perf re-test would need"
    echo "   matching 32-rank runs at multiple N values)"
else
    echo "  (couldn't parse qs_mol_dyn_low timing from cp2k.out — see file directly)"
fi
echo ""
echo "rundir: $RUNDIR"
