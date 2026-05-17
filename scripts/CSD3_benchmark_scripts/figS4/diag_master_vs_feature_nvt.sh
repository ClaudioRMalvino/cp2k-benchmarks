#!/usr/bin/env bash
# Critical test: does upstream master show the same step-wise PES jumps
# that the feature/nnp-native-spline binary showed, when both run NVT from
# the SAME GEO_OPT-minimised restart of the N=128 Morawietz setup?
#
# If master is smooth (no jumps in Pot, no ConsQty drift): the discontinuity
# is in the feature branch's neighbour-list / Verlet-skin cache, written by
# this project — fixable.
#
# If master shows the same jumps: it's deeper than the Verlet code (in
# the shared NNP machinery), upstream concern.
#
# Tiny job: 4 cores, 200 NVT steps from existing restart.

#SBATCH -J figS4_diag_masterNVT
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/diag_masterNVT_%j.out
# Budget: 4 cores × 10 min = 0.67 CPU-hr

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

CP2K_MASTER=$BIN_ROOT/master/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/master/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

FEATURE_DIR=/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N128_223052
RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_masterNVT_$(date +%H%M%S)
mkdir -p "$RUNDIR"
ln -sfn /home/crm98/cp2k-benchmarks/potentials "$RUNDIR/NNP"

# Reuse the feature branch's GEO_OPT restart so both binaries start from
# identical positions.
cp "$FEATURE_DIR/relax_N128-1.restart" "$RUNDIR/relax_input.restart"

# Generate an NVT input identical to feature's 02_nvt.inp but:
#   - PROJECT renamed (avoid file collisions)
#   - STEPS reduced to 200 (we only need to see the step-by-step Pot trace)
#   - EXT_RESTART pointing to our copied restart
python3 <<PYEOF
import re
src = open("$FEATURE_DIR/02_nvt.inp").read()
src = re.sub(r'PROJECT\s+\S+',  'PROJECT  masterNVT_N128', src)
src = re.sub(r'STEPS\s+\d+',    'STEPS 200', src, count=1)
src = re.sub(r'RESTART_FILE_NAME\s+\S+',
             'RESTART_FILE_NAME relax_input.restart', src)
open("$RUNDIR/run.inp", "w").write(src)
PYEOF

cd "$RUNDIR"
echo "=== running NVT (200 steps) with MASTER binary on feature's GEO_OPT minimum ==="
echo "    rundir: $RUNDIR"
srun --ntasks=4 --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_MASTER" -i run.inp > cp2k.out 2>&1 || true

echo ""
ener=$(ls *.ener 2>/dev/null | head -1)
if [[ -z "$ener" ]]; then
    echo "*** no .ener — master failed to start. cp2k.out tail:"
    tail -30 cp2k.out
    exit 1
fi

echo "=== MASTER trajectory (step-by-step, first 32 + last 5) ==="
awk 'NR>1 && NR<=33 {printf "  step=%-3s  T=%7.1f  Pot=%10.5f  ConsQty=%10.5f\n", $1, $4, $5, $6}' "$ener"
echo "  ..."
tail -5 "$ener" | awk '{printf "  step=%-3s  T=%7.1f  Pot=%10.5f  ConsQty=%10.5f\n", $1, $4, $5, $6}'

echo ""
echo "=== VERDICT ==="
python3 <<PYEOF
import re

def trace(path):
    out = []
    with open(path) as f:
        for ln in f:
            if ln.startswith('#') or not ln.strip(): continue
            parts = ln.split()
            if len(parts) < 6: continue
            try: out.append((int(parts[0]), float(parts[3]), float(parts[4]), float(parts[5])))
            except: pass
    return out

m = trace("$ener")
# Look for any single-step Pot jump > 0.1 hartree (= 2.7 eV)
big_pot_jumps = sum(1 for i in range(1, len(m)) if abs(m[i][2] - m[i-1][2]) > 0.10)
big_cq_jumps  = sum(1 for i in range(1, len(m)) if abs(m[i][3] - m[i-1][3]) > 0.10)
maxT = max(r[1] for r in m)
finalT = m[-1][1]
cq_drift = abs(m[-1][3] - m[0][3])

print(f"  steps recorded         : {len(m)}")
print(f"  T (final / max)         : {finalT:.1f} K / {maxT:.1f} K")
print(f"  ConsQty drift (|end−0|) : {cq_drift:.4f} hartree")
print(f"  single-step Pot jumps  > 0.10 ha : {big_pot_jumps}")
print(f"  single-step ConsQty jumps > 0.10 ha : {big_cq_jumps}")
print()

# Comparison reference: feature showed Pot jump 0.23 ha at step 4, 0.19 at 17, 0.32 at 29 (3 jumps in 30 steps)
if big_pot_jumps == 0 and cq_drift < 0.01:
    print("  -> MASTER IS SMOOTH. The discontinuities are in feature/nnp-native-spline's")
    print("     neighbour-list / Verlet-skin cache (your code). FIXABLE.")
elif big_pot_jumps > 0:
    print(f"  -> MASTER ALSO SHOWS {big_pot_jumps} JUMPS. Bug is in shared NNP machinery, not Verlet cache.")
    print("     Upstream concern.")
else:
    print("  -> MASTER SHOWS NO JUMPS but ConsQty drift is high — possibly different failure mode")
PYEOF
