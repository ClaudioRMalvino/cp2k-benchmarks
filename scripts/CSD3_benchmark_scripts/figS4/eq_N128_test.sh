#!/usr/bin/env bash
# Equilibration-rescue diagnostic for N=128 (Morawietz NNP).
#
# FD test established the heating is real dynamics from a non-equilibrium
# starting config, NOT a bad force calc. Three standard remedies, run in
# parallel so we can see which (if any) drains the excess PE without the
# system running away:
#
#   tight  : TIMECON 10 fs (was 100), MD for 2 ps. Minimal change — if this
#            works, every production input just needs that one-line edit.
#   anneal : T_target = 1 K, MD for 2 ps. CSVR will drain KE aggressively,
#            letting the system slide to its local PES minimum slowly. We
#            then check whether T stays bounded.
#   geoopt : RUN_TYPE GEO_OPT (BFGS, 200 iters, force tol 1e-3 ha/bohr) —
#            removes the ~9 kcal/mol/H2O excess before MD even starts.
#            Then no thermostat-rescue is needed.

#SBATCH -J figS4_eq_N128_test
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/eq_N128_test_%j.out

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"

set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

BENCH=/home/crm98/cp2k-benchmarks
rundir=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/eq_N128_test_$(date +%H%M%S)
mkdir -p "$rundir"
ln -sfn "$BENCH/potentials" "$rundir/NNP"

for f in eq_N128_tight eq_N128_anneal eq_N128_geoopt; do
   cp "$BENCH/${f}.inp" "$rundir/${f}.inp"
done

cd "$rundir"
for f in eq_N128_tight eq_N128_anneal eq_N128_geoopt; do
   echo "=== running $f ==="
   srun --ntasks=32 --cpus-per-task=1 --hint=nomultithread \
        "$CP2K_EXE" -i "${f}.inp" > "${f}.out" 2>&1 || true
done

echo ""
echo "==============================================================="
echo "  VERDICTS"
echo "==============================================================="
python3 <<PYEOF
import re, os, glob

def parse_ener(path):
    """Return list of (step, T, Pot, ConsQty) from CP2K .ener file."""
    out = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'): continue
            parts = ln.split()
            if len(parts) < 6: continue
            try:
                step = int(parts[0]); T = float(parts[3])
                pot  = float(parts[4]); cq = float(parts[5])
                out.append((step, T, pot, cq))
            except ValueError:
                continue
    return out

# Find .ener files (MD runs produce them) for tight + anneal
for label, prefix, note in [
    ("TIGHT (TIMECON=10 fs)",  "H2O-128_tight_thermostat", "want T to stay bounded; baseline standalone hit 1261 K in 200 steps"),
    ("ANNEAL (T_target=1K)",   "H2O-128_anneal",           "want T to fall toward 1 K"),
]:
    print(f"\n--- {label} ---")
    print(f"  ({note})")
    eners = sorted(glob.glob(f"{prefix}*-1.ener"))
    if not eners:
        print(f"  no .ener file found for {prefix} — likely crashed early")
        continue
    rows = parse_ener(eners[0])
    if not rows:
        print(f"  .ener file empty")
        continue
    sampler = [rows[0]] + [r for r in rows if r[0] in (10, 50, 100, 500, 1000, 2000, 4000)]
    seen = set()
    for r in sampler:
        if r[0] in seen: continue
        seen.add(r[0])
        s, T, pot, cq = r
        print(f"    step={s:5d}  T={T:8.1f} K  Pot={pot:12.5f}  ConsQty={cq:12.5f}")
    maxT = max(r[1] for r in rows)
    minT = min(r[1] for r in rows)
    finalT = rows[-1][1]
    print(f"  Tmin={minT:.1f}, Tmax={maxT:.1f}, Tfinal={finalT:.1f} K  ({len(rows)} steps recorded)")
    if maxT < 500:
        print(f"  -> THERMOSTAT HELD  (max T = {maxT:.1f} K, well below explosion regime)")
    elif maxT < 1000:
        print(f"  -> PARTIAL HELP  (max T = {maxT:.1f} K, less catastrophic than baseline)")
    else:
        print(f"  -> STILL DIVERGED  (max T = {maxT:.1f} K — thermostat couldn't keep up)")

# GEO_OPT diagnostic
print(f"\n--- GEO_OPT (BFGS minimisation, 200 iter, force tol 1e-3) ---")
gopt_outs = sorted(glob.glob("eq_N128_geoopt.out"))
if not gopt_outs:
    print("  no geoopt output found")
else:
    txt = open(gopt_outs[0]).read()
    # CP2K writes "GEOMETRY OPTIMIZATION COMPLETED" on success
    converged = "GEOMETRY OPTIMIZATION COMPLETED" in txt
    n_iter = len(re.findall(r"OPTIMIZATION STEP:\s+\d+", txt))
    # Initial and final energies
    energies = [float(m.group(1)) for m in re.finditer(
        r"Total FORCE_EVAL.*?\[hartree\]\s+(-?\d+\.\d+(?:[eE][+-]?\d+)?)", txt)]
    e_init  = energies[0]  if energies else None
    e_final = energies[-1] if energies else None
    print(f"  converged   : {converged}")
    print(f"  iterations  : {n_iter}")
    if e_init is not None:
        print(f"  E_initial   = {e_init:.6f} hartree")
        print(f"  E_final     = {e_final:.6f} hartree")
        print(f"  Delta E     = {e_final - e_init: .4e} hartree  (= {(e_final-e_init)*627.51:.2f} kcal/mol total)")
        print(f"               = {(e_final - e_init)/128*627.51:.3f} kcal/mol per H2O drained")
        if abs(e_final - e_init) > 0.1:
            print(f"  -> GEO_OPT successfully drained substantial PE — use this minimised config as MD start")

print()
print(f"rundir: {os.getcwd()}")
PYEOF
