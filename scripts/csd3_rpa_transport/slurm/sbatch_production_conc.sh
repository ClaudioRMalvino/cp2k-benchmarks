#!/usr/bin/env bash
#SBATCH -J rpa_prod
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_prod_%A_%a.out

# Transport campaign: one 80 ps NVE production segment per array task, from the
# concentration's equil snapshots. 80 ps = 160k steps x 0.2554 s/step ~ 11.4 h ->
# one 12 h window. 5x80 design (2026-07-30): RPA decorrelates ~2x faster than MP2,
# so 5x80 matches the kappa bars validated on the MP2 anchor at 5x160.
# Resume sweep = insurance: finished segments skipped, killed ones continue from
# 5 ps checkpoints (RESTART_COUNTERS fix).
#   p=$(sbatch --parsable --array=1-5 sbatch_production_conc.sh cubic_1M production_transport)
#   sbatch --array=1-5 --dependency=afterany:$p sbatch_production_conc.sh cubic_1M production_transport
# args: <conc_dir> [prod_subdir]; 1 m pilot -> production_transport (anchor segments untouched).
set -euo pipefail
CONC=${1:?usage: sbatch --array=1-5 sbatch_production_conc.sh <conc_dir> [prod_subdir]}
PSUB=${2:-production}
SEG=${SLURM_ARRAY_TASK_ID:?submit with --array=1-5}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport
export BIN_LABEL=dhruv-cell-list
# CELL env selects the box (default cube3); cube2 = 1 m Yeh-Hummer partner
CELL="${CELL:-cube3}"
case "$CELL" in cube2) FRANKS=4 ;; *) FRANKS=8 ;; esac

MODEL=RPA
echo "node: $(hostname)  model: $MODEL  conc: $CONC  cell: $CELL  segment: $SEG  ->  $PSUB/  binary: $BIN_LABEL"
t0=$SECONDS
MODEL="$MODEL" CONC_DIR="$CONC" PROD_SUBDIR="$PSUB" CELLS="$CELL" SEGMENTS="$SEG" \
  TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=$FRANKS PROD_PS=80 \
  bash "$ADIR/run_production.sh"
wall=$(( SECONDS - t0 ))

out=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs/$MODEL/$CONC/$PSUB/$CELL/seg$SEG/prod.out
tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' "$out" || true)
ps=""
[ -n "${tstep:-}" ] && ps=$(awk -v t="$tstep" 'BEGIN{printf "%.4f", t/160000}')
bash "$TDIR/log_timing.sh" prod_seg$SEG "$CONC" "$BIN_LABEL" "$SLURM_NTASKS" 160000 "$wall" "${ps:-NA}" "$CELL $MODEL 80ps NVE; s_per_step valid only for single-window runs"
echo "PRODUCTION SEGMENT DONE ($CONC seg$SEG), wall ${wall}s"
