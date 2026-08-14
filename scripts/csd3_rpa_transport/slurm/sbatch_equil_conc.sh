#!/usr/bin/env bash
#SBATCH -J rpa_equil
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_equil_%j.out

# MP2 transport campaign: NVT (CSVR) equilibration of one concentration's cube3
# cell on the PR #5295 engine; 30 ps equil + 5 snapshots 15 ps apart (210k steps).
# Resubmit-idempotent (continues from checkpoint). Usage: sbatch <this> <conc_dir> [seed];
# seed REQUIRED for the 1 m redo (independence from the anchor equil on the same cell).
set -euo pipefail
CONC=${1:?usage: sbatch sbatch_equil_conc.sh <conc_dir e.g. cubic_2m> [seed]}
SEED=${2:-}
# CELL selects the box (default cube3); cube2 (24.84 A, 1500 atoms) = 1 m Yeh-Hummer
# partner. Fist rank split per anchor split-sweep winners (cube2 -> 4, cube3 -> 8).
CELL="${CELL:-cube3}"
case "$CELL" in cube2) FRANKS=4 ;; *) FRANKS=8 ;; esac
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport
export BIN_LABEL=dhruv-cell-list

MODEL=RPA
echo "node: $(hostname)  model: $MODEL  conc: $CONC  cell: $CELL  binary: $BIN_LABEL  seed: ${SEED:-default}"
t0=$SECONDS
MODEL="$MODEL" CONC_DIR="$CONC" CELLS="$CELL" TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=$FRANKS \
  SEED_OVERRIDE="$SEED" bash "$ADIR/run_equil.sh"
wall=$(( SECONDS - t0 ))

out=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs/$MODEL/$CONC/equil/$CELL/equil.out
tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' "$out" || true)
ps=""
if [ -n "${tstep:-}" ]; then
  # s/step only meaningful when this invocation ran the full 210k steps
  ps=$(awk -v t="$tstep" 'BEGIN{printf "%.4f", t/210000}')
fi
bash "$TDIR/log_timing.sh" equil "$CONC" "$BIN_LABEL" "$SLURM_NTASKS" 210000 "$wall" "${ps:-NA}" "$CELL $MODEL NVT+snapshots; s_per_step valid only for single-window runs"
echo "EQUIL DONE ($CONC), wall ${wall}s"
