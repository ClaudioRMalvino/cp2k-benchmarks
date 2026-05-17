#!/usr/bin/env bash
# Two diagnostics in one job, on N=64 (stable) and N=128 (unstable):
#
#  TEST 1 — finite-difference energy-force consistency.
#    For each system: single-point ENERGY_FORCE at original positions, then
#    shift the first O atom's x by delta = 1e-4 A and re-evaluate. Compare
#       FD       :  -(E_shift - E_orig) / delta
#       reported :   F_x of that atom from the original run
#    If they agree to ~5+ digits, the NNP force IS the energy gradient
#    (energy-force consistent). If they disagree, the implementation
#    returns a force that is NOT -grad(E) — definitive code bug.
#
#  TEST 2 — extrapolation incidence.
#    Built-in NNP extrapolation tracker fires when ANY symmetry-function
#    value exceeds the [min, max] range from the training set (per
#    scaling.data, threshold 1e-4). If N=128 extrapolates from step 0 but
#    N=64 does not, the NNP is being used outside its trained regime
#    from the very start of the simulation — extrapolation cascade, not
#    necessarily a CP2K bug.
#
# Both tests use ENERGY_FORCE (single-shot, no MD), so the job is short.

#SBATCH -J figS4_diag_fd
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:15:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/diag_fd_extrapol_%j.out

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
rundir=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_fd_extrapol_$(date +%H%M%S)
mkdir -p "$rundir"
ln -sfn "$BENCH/potentials" "$rundir/NNP"

# Stage inputs
for f in fd_N64_orig fd_N64_shift fd_N128_orig fd_N128_shift; do
   cp "$BENCH/${f}.inp" "$rundir/${f}.inp"
done

# Run all four single-point evaluations (each takes seconds)
cd "$rundir"
for f in fd_N64_orig fd_N64_shift fd_N128_orig fd_N128_shift; do
   echo "=== running $f ==="
   srun --ntasks=32 --cpus-per-task=1 --hint=nomultithread \
        "$CP2K_EXE" -i "${f}.inp" > "${f}.out" 2>&1 || true
done

echo ""
echo "==============================================================="
echo "  PARSED RESULTS"
echo "==============================================================="
python3 <<PYEOF
import re, os, glob

DELTA = 1.0e-4    # Angstrom (must match what fd_N*_shift.inp was generated with)

def parse_energy(path):
    """Total NNP energy in Hartree from the ENERGY_FORCE output."""
    with open(path) as f:
        txt = f.read()
    # CP2K writes:  ENERGY| Total FORCE_EVAL ( NNP ) energy [hartree]   -42.123456...
    m = re.search(r'ENERGY\|\s+Total FORCE_EVAL.*?\[hartree\]\s+(-?\d+\.\d+([eE][+-]?\d+)?)', txt)
    if not m:
        return None
    return float(m.group(1))

def parse_first_atom_fx(path):
    """F_x of atom #1 (the perturbed O) from the printed FORCES table."""
    with open(path) as f:
        lines = f.readlines()
    # Find a block like:
    # ATOMIC FORCES in [a.u.]
    # # Atom   Kind   Element          X              Y              Z
    #      1      1      O      0.123E+00  ...
    for i, ln in enumerate(lines):
        if 'ATOMIC FORCES' in ln:
            for j in range(i+1, min(i+6, len(lines))):
                m = re.match(r'\s*1\s+\d+\s+\S+\s+(-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+', lines[j])
                if m:
                    return float(m.group(1))
    return None

def count_extrapolations(path):
    """Number of 'NNP extrapolation point N =' lines emitted, plus the
    final counter value if printed."""
    n_hits = 0
    last_N = None
    with open(path) as f:
        for ln in f:
            if 'NNP extrapolation point N =' in ln:
                n_hits += 1
                try:
                    last_N = int(ln.split('=')[-1].strip())
                except ValueError:
                    pass
    return n_hits, last_N

# Bohr_to_Angstrom; CP2K forces are in [a.u.] = hartree/bohr.
# Energies in hartree, displacements in Angstrom → need to convert.
BOHR_PER_A = 1.0 / 0.529177210903
H_PER_eV   = 1.0 / 27.211386245988

for label, orig, shift in [
    ('N=64  (cubic 12.42, stable)',          'fd_N64_orig.out',  'fd_N64_shift.out'),
    ('N=128 (24.84 x 12.42 x 12.42, fails)', 'fd_N128_orig.out', 'fd_N128_shift.out'),
]:
    print(f"\n--- {label} ---")
    E_orig  = parse_energy(orig)
    E_shift = parse_energy(shift)
    Fx_rep  = parse_first_atom_fx(orig)
    n_ext_orig, last_orig   = count_extrapolations(orig)
    n_ext_shift, last_shift = count_extrapolations(shift)

    if None in (E_orig, E_shift, Fx_rep):
        print(f"  PARSE FAILED  E_orig={E_orig}  E_shift={E_shift}  Fx_reported={Fx_rep}")
        print(f"  (check {orig}/{shift} manually)")
        continue

    # F_x (atomic units, hartree/bohr) is reported by CP2K.
    # FD: -(E_shift - E_orig) [hartree] / (DELTA [A] * BOHR_PER_A [bohr/A])
    fd_force = -(E_shift - E_orig) / (DELTA * BOHR_PER_A)   # hartree/bohr
    rel_err  = abs(fd_force - Fx_rep) / max(abs(Fx_rep), 1e-30)

    print(f"  E_orig             = {E_orig:.10f} hartree")
    print(f"  E_shift            = {E_shift:.10f} hartree")
    print(f"  Delta E            = {E_shift - E_orig: .6e} hartree   (over dx = {DELTA:.1e} A)")
    print(f"  -dE/dx (FD)        = {fd_force: .10e} hartree/bohr")
    print(f"  F_x reported       = {Fx_rep: .10e} hartree/bohr")
    print(f"  |FD - F| / |F|     = {rel_err:.3e}")
    if rel_err < 1e-3:
        print(f"  -> CONSISTENT  (F is the gradient of E to within FD step error)")
    elif rel_err < 1e-1:
        print(f"  -> MARGINAL  (could be DELTA too coarse — re-run with smaller delta to confirm)")
    else:
        print(f"  -> INCONSISTENT  (force is NOT -dE/dx — definitive energy-force bug)")

    print(f"  extrapolation frames (orig):  {n_ext_orig}  (counter ended at N={last_orig})")
    print(f"  extrapolation frames (shift): {n_ext_shift} (counter ended at N={last_shift})")
    if n_ext_orig == 0:
        print(f"  -> initial config IS within training-set symfunc range")
    else:
        print(f"  -> initial config is OUTSIDE training-set symfunc range  (NNP extrapolating from step 0)")

print()
print("rundir:", os.getcwd())
PYEOF
