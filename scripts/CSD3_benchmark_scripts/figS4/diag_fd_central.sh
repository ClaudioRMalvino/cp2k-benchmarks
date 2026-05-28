#!/usr/bin/env bash
# Follow-up to diag_fd_extrapol: central-difference FD test.
#
# Forward FD with delta=1e-4 A gave ~1.5% rel error on BOTH N=64 (stable) and
# N=128 (unstable), with bit-identical absolute error magnitudes. That is the
# fingerprint of O(delta) truncation, not energy-force inconsistency. Central
# FD has O(delta^2) truncation so the same delta yields ~3e-5 rel error if
# the implementation is consistent; if either system has a real F != -dE/dx
# bug, central FD will still show large error there.
#
# Only the new -delta evaluations are run here; +delta and orig are reused.

#SBATCH -J figS4_diag_fd_central
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/diag_fd_central_%j.out

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
PREV=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_fd_extrapol_212925
rundir=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/diag_fd_central_$(date +%H%M%S)
mkdir -p "$rundir"
ln -sfn "$BENCH/potentials" "$rundir/NNP"

# Copy the -delta inputs + reuse the previously-generated orig/+delta outputs
for f in fd_N64_minus fd_N128_minus; do
   cp "$BENCH/${f}.inp" "$rundir/${f}.inp"
done
for f in fd_N64_orig fd_N64_shift fd_N128_orig fd_N128_shift; do
   cp "$PREV/${f}.out" "$rundir/${f}.out"
done

# Run only the new -delta evals
cd "$rundir"
for f in fd_N64_minus fd_N128_minus; do
   echo "=== running $f ==="
   srun --ntasks=4 --cpus-per-task=1 --hint=nomultithread \
        "$CP2K_EXE" -i "${f}.inp" > "${f}.out" 2>&1 || true
done

echo ""
echo "==============================================================="
echo "  CENTRAL-DIFFERENCE FD vs CP2K FORCE"
echo "==============================================================="
python3 <<PYEOF
import re, os

DELTA  = 1.0e-4
BOHR_A = 1.0 / 0.529177210903

def parse_E(path):
    txt = open(path).read()
    m = re.search(r'ENERGY\|\s+Total FORCE_EVAL.*?\[hartree\]\s+(-?\d+\.\d+(?:[eE][+-]?\d+)?)', txt)
    return float(m.group(1)) if m else None

def parse_F1x(path):
    txt = open(path).read()
    m = re.search(r'^\s*FORCES\|\s+1\s+(-?\d+\.\d+E[+-]\d+)', txt, flags=re.M)
    return float(m.group(1)) if m else None

for label, prefix in [
    ("N=64  (cubic 12.42 — STABLE)",       "fd_N64"),
    ("N=128 (24.84x12.42x12.42 — FAILS)",  "fd_N128"),
]:
    Eo = parse_E(f"{prefix}_orig.out")
    Ep = parse_E(f"{prefix}_shift.out")    # E(+delta)
    Em = parse_E(f"{prefix}_minus.out")    # E(-delta)
    Fx = parse_F1x(f"{prefix}_orig.out")

    fd_fwd  = -(Ep - Eo) / (DELTA * BOHR_A)             # O(delta)
    fd_ctr  = -(Ep - Em) / (2 * DELTA * BOHR_A)         # O(delta^2)
    err_fwd = abs(fd_fwd - Fx); rel_fwd = err_fwd / max(abs(Fx), 1e-30)
    err_ctr = abs(fd_ctr - Fx); rel_ctr = err_ctr / max(abs(Fx), 1e-30)

    print(f"--- {label} ---")
    print(f"  E_orig       = {Eo:.12f} ha")
    print(f"  E(+delta)    = {Ep:.12f} ha")
    print(f"  E(-delta)    = {Em:.12f} ha")
    print(f"  Reported F_x = {Fx: .10e} ha/bohr")
    print(f"  Forward FD   = {fd_fwd: .10e} ha/bohr  |err|/|F| = {rel_fwd:.3e}")
    print(f"  Central FD   = {fd_ctr: .10e} ha/bohr  |err|/|F| = {rel_ctr:.3e}")
    print()
    if rel_ctr < 1e-4:
        v = "CONSISTENT  (central-FD agrees with F to 4+ digits — F IS -grad E)"
    elif rel_ctr < 1e-2:
        v = "MARGINAL    (some inconsistency at this delta — try delta = 1e-5)"
    else:
        v = "INCONSISTENT (central FD still disagrees — definitive code bug)"
    print(f"  verdict      : {v}")
    print()
PYEOF
