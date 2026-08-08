#!/usr/bin/env bash
#SBATCH -J dft_probe
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/dft_probe_%j.out

# Short probe for the 5064-atom (37.26 A) rung of the on-the-fly DFT ladder.
#
# THE PROBLEM IT DIAGNOSES: job 33128702 ran 6 h 25 min on 76 pure-MPI ranks
# and never printed a single SCF iteration. A gdb backtrace put every rank in
#     qs_forces -> qs_energies_init_hamiltonians -> qs_env_update_s_mstruct
#     -> calculate_rho_core_c1d_gs -> transfer_rs2pw_distributed -> mp_waitany
# i.e. stuck in the realspace->planewave halo exchange for the core density,
# before the SCF even starts. The node was not short of memory (80 GB of
# 502 GB used), so this is communication, not capacity.
#
# THE HYPOTHESIS: at 1200 Ry the fine realspace grid is roughly 776^3. Split
# over 76 ranks that is ~10 planes of slab per rank, while the collocation
# halo (set by the Gaussian radii, which EPS_DEFAULT 1e-12 / EPS_PGF_ORB 1e-14
# make generous) is far wider than 10 planes. Every rank then exchanges many
# times its own slab with many neighbours and the transfer degenerates. Fewer,
# fatter slabs should fix it - hence MPI x OpenMP instead of pure MPI.
#
# Variants:
#   base        - hybrid ranks x threads, deck otherwise untouched
#   replicated  - additionally force &RS_GRID DISTRIBUTION_TYPE REPLICATED, so
#                 the halo exchange becomes a plain reduction over full grids
#                 (~3.7 GB per rank at 19 ranks = ~70 GB, comfortable)
#
# Accuracy settings (cutoff, EPS_*, XC, basis, D3) are NOT touched by any
# variant: those define the level of theory the NNP was fitted to. Rank/thread
# layout and grid distribution are solver bookkeeping and change no number.
#
#   sbatch --ntasks=19 --cpus-per-task=4 sbatch_dft_probe.sh 19 4 base
set -euo pipefail
RANKS=${1:?usage: sbatch_dft_probe.sh <ranks> <threads> <base|replicated>}
THREADS=${2:?}
VARIANT=${3:-base}

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

echo "node: $(hostname)  layout: ${RANKS} ranks x ${THREADS} threads = $((RANKS*THREADS)) cores  variant: $VARIANT"
[ -f "$COORD" ] || { echo "FATAL: missing $COORD"; exit 6; }
mkdir -p "$RUNDIR"; rm -f "$RUNDIR/$proj"* "$RUNDIR/probe.out"

sed -e "s|__PROJECT__|$proj|" \
    -e "s|__NREP__|1 1 1|g" \
    -e "s|__STEPS__|2|" \
    -e "s|__DATA__|$DATA|g" \
    -e "s|__ABC__|37.2600|g" \
    -e "s|__COORD__|$COORD|" \
    "$TDIR/dft_ladder.inp.template" > "$RUNDIR/probe.inp"

# MEASURE plans every FFT size by timing trials; at ~776^3 that is a large
# one-off charge landing inside the very step we are trying to time.
sed -i 's|^  FFTW_PLAN_TYPE MEASURE$|  FFTW_PLAN_TYPE ESTIMATE|' "$RUNDIR/probe.inp"

if [ "$VARIANT" = replicated ]; then
  sed -i 's|^      CUTOFF 1200$|      CUTOFF 1200\n      \&RS_GRID\n        DISTRIBUTION_TYPE REPLICATED\n      \&END RS_GRID|' "$RUNDIR/probe.inp"
  grep -q "DISTRIBUTION_TYPE REPLICATED" "$RUNDIR/probe.inp" || { echo "FATAL: RS_GRID patch did not apply"; exit 7; }
fi

cd "$RUNDIR"
t0=$SECONDS
srun --ntasks="$RANKS" --cpus-per-task="$THREADS" --hint=nomultithread \
     "$BIN" -i probe.inp -o probe.out || echo "(srun exited non-zero - probe walls out by design)"
wall=$((SECONDS - t0))

echo "=== PROBE RESULT (${RANKS}x${THREADS} $VARIANT, wall ${wall}s) ==="
scf_iters=$(grep -cE "OT (DIIS|SD) " probe.out 2>/dev/null || echo 0)
echo "SCF iterations printed : $scf_iters"
echo "SCF converged banners  : $(grep -c 'SCF run converged' probe.out 2>/dev/null || echo 0)"
if [ -f "$proj-1.ener" ]; then
  echo "--- .ener (last col = UsedTime per step) ---"; cat "$proj-1.ener"
else
  echo "no .ener - did not complete a single MD step"
fi
echo "--- tail ---"; tail -6 probe.out 2>/dev/null || true
