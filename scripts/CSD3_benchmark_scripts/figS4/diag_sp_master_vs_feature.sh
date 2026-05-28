#!/usr/bin/env bash
# Single-point ENERGY_FORCE: master vs feature/nnp-native-spline on the same
# coordinates.
#
# Two configurations × two binaries = 4 evals (each ~1 s of compute):
#   sp_init   : initial config (H2O-128_explicit_RPBE-vdW_NNP.inp coords)
#   sp_relax  : feature's GEO_OPT minimum (relax_N128-1.restart)
#
# Reveals:
#   - whether master and feature agree on total E for identical positions
#   - if not: is the offset uniform (∝ N_atoms) or concentrated on specific atoms
#   - whether per-atom forces differ in proportion to per-atom energies

#SBATCH -J figS4_diag_sp
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/diag_sp_%j.out
# Budget: 4 cores × 10 min = 0.67 CPU-hr

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

MASTER=$BIN_ROOT/master/cp2k.psmp
FEATURE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
[[ ! -x "$MASTER"  ]] && { echo "missing $MASTER";  exit 1; }
[[ ! -x "$FEATURE" ]] && { echo "missing $FEATURE"; exit 1; }
export OMP_NUM_THREADS=1

BENCH=/home/crm98/cp2k-benchmarks
RESTART=/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N128_223052/relax_N128-1.restart
[[ ! -f "$RESTART" ]] && { echo "GEO_OPT restart missing: $RESTART"; exit 1; }

RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_sp_$(date +%H%M%S)
mkdir -p "$RUNDIR"; cd "$RUNDIR"
ln -sfn "$BENCH/potentials" NNP
cp "$RESTART" relax_input.restart

# Materialise the two input templates, filling in the restart path for sp_relax
cp "$BENCH/diffusion/sp_init.inp"  sp_init.inp
sed "s|__RESTART__|relax_input.restart|" "$BENCH/diffusion/sp_relax.inp" > sp_relax.inp

# Run the four combinations.  Subdir per run so .ener/.xyz/.restart don't collide.
for cfg in init relax; do
   for binlabel in master feature; do
      sub="${cfg}_${binlabel}"
      mkdir -p "$sub"
      cp "sp_${cfg}.inp"          "$sub/run.inp"
      cp relax_input.restart      "$sub/" 2>/dev/null || true
      ln -sfn "$BENCH/potentials" "$sub/NNP"
      case "$binlabel" in
         master)  CP2K_EXE=$MASTER;  LIB=$BIN_ROOT/master/lib  ;;
         feature) CP2K_EXE=$FEATURE; LIB=$BIN_ROOT/feature-nnp-native-spline/lib ;;
      esac
      export LD_LIBRARY_PATH="$LIB:${LD_LIBRARY_PATH:-}"
      echo "=== running $sub ==="
      ( cd "$sub" && srun --ntasks=4 --cpus-per-task=1 --hint=nomultithread \
              "$CP2K_EXE" -i run.inp > cp2k.out 2>&1 ) || true
   done
done

echo ""
echo "==============================================================="
echo "  PARSED RESULTS"
echo "==============================================================="
python3 <<'PYEOF'
import re, os, glob, statistics

def parse_total_E(path):
    """Total NNP-method energy (hartree)."""
    txt = open(path).read()
    m = re.search(r'ENERGY\|\s+Total FORCE_EVAL.*?\[hartree\]\s+(-?\d+\.\d+(?:[eE][+-]?\d+)?)', txt)
    return float(m.group(1)) if m else None

def parse_forces(path):
    """Parse per-atom forces from the FORCES| block. Returns dict {atom_idx: (Fx,Fy,Fz)}."""
    out = {}
    with open(path) as f:
        for ln in f:
            m = re.match(r'^\s*FORCES\|\s+(\d+)\s+(-?\d+\.\d+E[+-]\d+)\s+(-?\d+\.\d+E[+-]\d+)\s+(-?\d+\.\d+E[+-]\d+)', ln)
            if m:
                out[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return out

def parse_per_atom_E(path):
    """Try to recover per-atom NNP energies from &ENERGIES HIGH printout.
       CP2K may format these as 'NNP|  Atom XXX  energy = -0.XXXX' or similar.
       If absent, return empty dict."""
    out = {}
    with open(path) as f:
        for ln in f:
            # try a few common formats
            m = re.match(r'^\s*NNP\|\s+(?:atom\s+|Atom\s+)?(\d+)\s+(?:energy\s*[:=]?\s*)?(-?\d+\.\d+(?:[eE][+-]?\d+)?)', ln)
            if m:
                out[int(m.group(1))] = float(m.group(2))
    return out

def fnorm(vec):
    return (vec[0]**2 + vec[1]**2 + vec[2]**2) ** 0.5

for cfg in ['init', 'relax']:
    print(f"\n--- configuration: {cfg} ---")
    Em = parse_total_E(f"{cfg}_master/cp2k.out")
    Ef = parse_total_E(f"{cfg}_feature/cp2k.out")
    if Em is None or Ef is None:
        print(f"  PARSE FAIL  E_master={Em}  E_feature={Ef}")
        print(f"  tail of master/cp2k.out:"); os.system(f"tail -10 {cfg}_master/cp2k.out | sed 's/^/    /'")
        continue
    dE = Ef - Em
    print(f"  E_master         = {Em:.10f} ha")
    print(f"  E_feature        = {Ef:.10f} ha")
    print(f"  ΔE (feature-master) = {dE:+.6e} ha  =  {dE*27.2114:+.6e} eV")
    n_h2o = 128 if cfg == 'init' else 128
    print(f"  per H2O          = {dE/n_h2o*627.51:+.4f} kcal/mol/H2O")

    Fm = parse_forces(f"{cfg}_master/cp2k.out")
    Ff = parse_forces(f"{cfg}_feature/cp2k.out")
    if not Fm or not Ff or set(Fm) != set(Ff):
        print(f"  forces missing or atom-set mismatch (master {len(Fm)}, feature {len(Ff)})")
        continue
    diffs = []
    rel_diffs = []
    for k in Fm:
        d = (Ff[k][0]-Fm[k][0], Ff[k][1]-Fm[k][1], Ff[k][2]-Fm[k][2])
        nd = fnorm(d)
        nf = fnorm(Fm[k])
        diffs.append(nd)
        if nf > 1e-12: rel_diffs.append(nd/nf)
    diffs.sort(reverse=True)
    rel_diffs.sort(reverse=True)
    print(f"  per-atom |ΔF|  median = {statistics.median(diffs):.3e} ha/bohr")
    print(f"                 max    = {diffs[0]:.3e} ha/bohr (atom {[k for k in Fm if fnorm((Ff[k][0]-Fm[k][0],Ff[k][1]-Fm[k][1],Ff[k][2]-Fm[k][2]))==diffs[0]][0]})")
    print(f"                 95th % = {diffs[int(0.05*len(diffs))]:.3e} ha/bohr")
    print(f"  per-atom rel |ΔF|/|F| median = {statistics.median(rel_diffs):.3e}")
    print(f"                              max    = {rel_diffs[0]:.3e}")

    Em_at = parse_per_atom_E(f"{cfg}_master/cp2k.out")
    Ef_at = parse_per_atom_E(f"{cfg}_feature/cp2k.out")
    if Em_at and Ef_at and set(Em_at) == set(Ef_at):
        atomE_diffs = sorted([abs(Ef_at[k] - Em_at[k]) for k in Em_at], reverse=True)
        print(f"  per-atom ΔE     median = {statistics.median(atomE_diffs):.3e} ha")
        print(f"                  max    = {atomE_diffs[0]:.3e} ha")
        # if pattern is "uniform across all atoms" vs "concentrated", report
        spread = atomE_diffs[0] / max(statistics.median(atomE_diffs), 1e-30)
        print(f"                  max/median = {spread:.2f}  (high = concentrated, ~1 = uniform)")
    else:
        print(f"  (per-atom NNP energies not in expected print format — diag works without them)")

print()
print("Interpretation:")
print("  uniform per-atom ΔE  → systematic baseline offset (e.g. wrong atom_energies constant)")
print("  concentrated per-atom ΔE → specific atomic environments hit a bug")
print("  large rel|ΔF|/|F| → forces also wrong, not just energy")
print("  small rel|ΔF|/|F| but big ΔE → energy bookkeeping bug, forces consistent")
PYEOF
echo ""
echo "rundir: $(pwd)"
