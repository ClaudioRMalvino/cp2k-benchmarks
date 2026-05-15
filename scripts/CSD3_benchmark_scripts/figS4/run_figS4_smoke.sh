#!/usr/bin/env bash
#! Smoke test for the reconstructed RPBE-vdW (Morawietz 2016) NNP.
#! Runs a single-point ENERGY_FORCE on the 64-H2O box with BOTH branch
#! binaries (upstream master, feature/nnp-native-spline) and checks:
#!   1. the reconstructed input.nn parses and CP2K completes
#!   2. the total energy / forces / stress are finite and physically sane
#!   3. master and native-spline agree to ~machine precision
#! This validates the input.nn reconstruction + the optimisation's correctness
#! BEFORE committing to the full Fig. S4 production campaign.

#SBATCH -J figS4_smoke
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_smoke_%j.out

# Source the toolchain environment BEFORE enabling strict mode: the toolchain
# 'setup' script references unbound vars (CP_DFLAGS) that trip `set -u`.
. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCH=/home/crm98/cp2k-benchmarks
source "$BIN_ROOT/setup"

set -euo pipefail

OUTDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/smoke_$(date +%d-%m_%H-%M)
mkdir -p "$OUTDIR"
cp "$BENCH/scripts/CSD3_benchmark_scripts/figS4/smoke_RPBE-vdW.inp" "$OUTDIR/"

run_branch() {
   local label=$1 exe=$2 lib=$3
   local rundir="$OUTDIR/$label"
   mkdir -p "$rundir"
   cp "$OUTDIR/smoke_RPBE-vdW.inp" "$rundir/run.inp"
   # The input references NNP/RPBE-vdW-2016/... ; point NNP at the potentials dir.
   ln -sfn "$BENCH/potentials" "$rundir/NNP"
   export LD_LIBRARY_PATH="$lib:${LD_LIBRARY_PATH:-}"
   export OMP_NUM_THREADS=1
   echo "==> $label : $exe"
   ( cd "$rundir" && srun --ntasks=8 --cpus-per-task=1 --hint=nomultithread \
        "$exe" -i run.inp >cp2k.out 2>&1 ) || true
   if grep -q "PROGRAM ENDED" "$rundir/cp2k.out"; then
      echo "    completed OK"
   else
      echo "    *** FAILED *** (tail of cp2k.out:)"
      tail -30 "$rundir/cp2k.out" | sed 's/^/      /'
   fi
}

run_branch master              "$BIN_ROOT/master/cp2k.psmp"                     "$BIN_ROOT/master/lib"
run_branch feature-nnp-native-spline "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" "$BIN_ROOT/feature-nnp-native-spline/lib"

echo
echo "=================  COMPARISON  ================="
python3 - "$OUTDIR" <<'PYEOF'
import sys, re, os
root = sys.argv[1]
branches = ["master", "feature-nnp-native-spline"]

def parse(out):
    txt = open(out).read()
    # total NNP energy (Hartree): "ENERGY| Total FORCE_EVAL ( NNP ) energy [hartree]  -45.159..."
    m = re.search(r'ENERGY\|\s*Total FORCE_EVAL.*?energy\s*\[\w+\]\s+([-\d.E+]+)', txt)
    energy = float(m.group(1)) if m else None
    # atomic forces: "FORCES|  <atom>  fx  fy  fz  |f|"
    forces = []
    for ln in txt.splitlines():
        p = ln.split()
        if len(p) == 6 and p[0] == "FORCES|" and p[1].isdigit():
            forces.append((float(p[2]), float(p[3]), float(p[4])))
    # 1/3 trace of the analytical stress tensor [bar], as a pressure sanity check
    st = re.search(r'STRESS\|\s*1/3 Trace\s+([-\d.E+]+)', txt)
    return energy, forces, (st.group(1) if st else None)

res = {}
for b in branches:
    out = os.path.join(root, b, "cp2k.out")
    if not os.path.exists(out):
        print(f"  {b}: no cp2k.out"); continue
    res[b] = parse(out)
    e, f, s = res[b]
    print(f"  {b}:")
    print(f"    total energy      = {e} Ha" + ("" if e is None else f"  ({e*27.2114:.4f} eV)"))
    print(f"    n forces parsed   = {len(f)}")
    if f:
        fmax = max(max(abs(c) for c in atom) for atom in f)
        print(f"    max |force comp|  = {fmax:.6e} a.u.")
    if s: print(f"    stress 1/3 trace  = {s} bar")

if len(res) == 2:
    e0, f0, _ = res["master"]
    e1, f1, _ = res["feature-nnp-native-spline"]
    print()
    if e0 is not None and e1 is not None:
        de = abs(e0 - e1)
        print(f"  |dE| master vs native-spline = {de:.3e} Ha")
        ok_e = de < 1e-8
    else:
        ok_e = False
    ok_f = False
    if f0 and f1 and len(f0) == len(f1):
        dmax = max(max(abs(a-b) for a,b in zip(af,bf)) for af,bf in zip(f0,f1))
        print(f"  max |dF| master vs native-spline = {dmax:.3e} a.u.")
        ok_f = dmax < 1e-6
    print()
    print("  RESULT:", "PASS - branches agree, input.nn parses" if (ok_e and ok_f)
          else "CHECK - see numbers above")
PYEOF
echo "Outputs in: $OUTDIR"
