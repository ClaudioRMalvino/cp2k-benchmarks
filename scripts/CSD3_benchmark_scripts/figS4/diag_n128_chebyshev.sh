#!/usr/bin/env bash
# Fig. S4 gate test for feature/nnp-chebyshev.
#
# The May 2026 diagnostics showed single-step PES jumps (>0.1 Ha) and
# kinetic runaway on the feature branch at N=128 (24.84x12.42x12.42
# replicated box) while upstream master was smooth, and while N=64 cubic
# was stable on BOTH.  The cell-list code has since been heavily revised
# (minimax branch + June optimisation pass).  This re-runs the exact
# failing configuration with the CURRENT install-tree binaries:
#
#   smooth on both  -> blocker gone, proceed to the full Fig. S4 chain
#   jumps on branch -> the neighbour-list bug survives; debug before S4
#
# Submit with:  sbatch diag_n128_chebyshev.sh

#SBATCH -J figS4_diag_chebyN128
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/diag_chebyN128_%j.out
# Budget: 8 cores x <=30 min = <=4 CPU-hr; typically ~0.5.

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
set -euo pipefail

MASTER=/home/crm98/cp2k_master/install/bin/cp2k.psmp
BRANCH=/home/crm98/cp2k_optimized/install/bin/cp2k.psmp
[[ -x "$MASTER" ]] || { echo "missing $MASTER"; exit 1; }
[[ -x "$BRANCH" ]] || { echo "missing $BRANCH"; exit 1; }
export OMP_NUM_THREADS=1

BENCH=/home/crm98/cp2k-benchmarks
FEATURE_DIR=/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N128_223052
RESTART="$FEATURE_DIR/relax_N128-1.restart"
[[ -f "$RESTART" ]] || { echo "missing GEO_OPT restart $RESTART"; exit 1; }
[[ -f "$BENCH/potentials/RPBE-vdW-2016/input.nn" ]] || \
    { echo "missing potentials/RPBE-vdW-2016 (restore from ~/.snapshot)"; exit 1; }

RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_cheby_$(date +%H%M%S)
mkdir -p "$RUNDIR"
cd "$RUNDIR"
cp "$RESTART" relax_input.restart

# 500-step NVT input derived from the original failing run's 02_nvt.inp.
python3 <<PYEOF
import re
src = open("$FEATURE_DIR/02_nvt.inp").read()
src = re.sub(r'PROJECT\s+\S+', 'PROJECT  chebyN128', src)
src = re.sub(r'STEPS\s+\d+',   'STEPS 500', src, count=1)
src = re.sub(r'RESTART_FILE_NAME\s+\S+',
             'RESTART_FILE_NAME relax_input.restart', src)
open("run.inp", "w").write(src)
PYEOF

for side in master branch; do
    case $side in
        master) exe=$MASTER ;;
        branch) exe=$BRANCH ;;
    esac
    mkdir -p $side
    cp run.inp relax_input.restart $side/
    ln -sfn "$BENCH/potentials" "$side/NNP"
    echo "=== [$side] 500 NVT steps on the N=128 failing configuration ==="
    ( cd $side && srun --ntasks=8 --cpus-per-task=1 --hint=nomultithread \
          "$exe" -i run.inp > cp2k.out 2>&1 ) || true
    if ! ls $side/*.ener >/dev/null 2>&1; then
        echo "*** [$side] no .ener file -- run failed. cp2k.out tail:"
        tail -30 $side/cp2k.out
    fi
done

echo
echo "================ VERDICT ================"
python3 <<'PYEOF'
import glob, os

def trace(path):
    out = []
    with open(path) as f:
        for ln in f:
            if ln.startswith('#') or not ln.strip(): continue
            p = ln.split()
            if len(p) < 6: continue
            try: out.append((int(p[0]), float(p[3]), float(p[4]), float(p[5])))
            except ValueError: pass
    return out

results = {}
for side in ("master", "branch"):
    eners = glob.glob(f"{side}/*.ener")
    if not eners:
        print(f"  {side}: NO OUTPUT")
        continue
    m = trace(eners[0])
    pot_jumps = sum(1 for i in range(1, len(m)) if abs(m[i][2]-m[i-1][2]) > 0.10)
    cq_jumps  = sum(1 for i in range(1, len(m)) if abs(m[i][3]-m[i-1][3]) > 0.10)
    maxT   = max(r[1] for r in m)
    drift  = abs(m[-1][3]-m[0][3])
    results[side] = (pot_jumps, cq_jumps, maxT, drift, len(m))
    print(f"  {side:7s} steps={len(m):4d}  T_max={maxT:7.1f} K  "
          f"ConsQty drift={drift:.4f} Ha  Pot jumps>0.1Ha={pot_jumps}  CQ jumps={cq_jumps}")

print()
if "branch" in results and "master" in results:
    bp, bc, bT, bd, _ = results["branch"]
    mp, mc, mT, md, _ = results["master"]
    if bp == 0 and bc == 0 and bT < 450.0 and bd < 5*max(md, 1e-3):
        print("  ==> BRANCH IS SMOOTH on the previously-failing N=128 configuration.")
        print("      The May neighbour-list discontinuity is GONE.  Fig. S4 chain is unblocked.")
    elif bp > 0 or bc > 0:
        print(f"  ==> BRANCH STILL JUMPS ({bp} Pot / {bc} ConsQty jumps; master: {mp}/{mc}).")
        print("      The cell-list bug survives -- debug before launching Fig. S4 production.")
    else:
        print("  ==> No discrete jumps, but check drift/T against master above.")
PYEOF
echo
echo "rundir: $RUNDIR"
