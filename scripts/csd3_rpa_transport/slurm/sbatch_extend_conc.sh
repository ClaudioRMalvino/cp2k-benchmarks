#!/usr/bin/env bash
#SBATCH -J rpa_ext
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_ext_%A_%a.out

# ROUND 2, CPU TRACK: extend a completed 80 ps NVE segment to 160 ps on one
# icelake node. The CPU twin of sbatch_gpu_extend_conc.sh - same EXTEND=1
# checkpoint logic, same 320k-step target, only the partition, the binary and
# RUN_ROOT differ (runs/ rather than runs_gpu/).
#
# WHY: the 5x80 ps design was sized against the MP2 anchor's kappa error bars,
# and the diffusion/conductivity-NE/viscosity numbers it produced are solid.
# kappa_Onsager is not: it scatters 3.0-33.7 S/m across the 4 m CPU segments,
# which is the signature of a correlation function that has not converged in
# 80 ps. The 4 m Delta_NE > 0 pairing signal has the same problem - it rests on
# the GPU track's tighter bars with the CPU track merely consistent. Both need
# 160 ps.
#
# This does NOT replace the queued ampere round-2 jobs (32644668/70/72,
# 32648051) - it runs alongside them. Extending the CPU segments on CPU keeps
# each hardware track's trajectory generated end-to-end on one architecture,
# so the 80 ps hardware cross-check stays clean and upgrades to 160 ps if both
# tracks land. Continuing a GPU trajectory on CPU would have forfeited that.
#
# Every CPU segment is checkpointed at exactly step 160000, so all 15 need the
# full 160000 additional steps: ~10.9 h at the measured 0.214-0.259 s/step,
# inside the 12 h QOS wall with ~1 h of margin. The resume twin covers a slow
# node via the capped-STEPS checkpoint logic in run_production.sh.
#
#   p=$(sbatch --parsable --array=1-5 sbatch_extend_conc.sh cubic_4m)
#   sbatch --array=1-5 --dependency=afterany:$p sbatch_extend_conc.sh cubic_4m
set -euo pipefail
CONC=${1:?usage: sbatch --array=1-5 sbatch_extend_conc.sh <conc_dir> [prod_subdir]}
PSUB=${2:-production}
SEG=${SLURM_ARRAY_TASK_ID:?submit with --array=1-5}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport

export BIN_LABEL=dhruv-cell-list          # the binary round 1 was produced with
RUN_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs
CELL=cube3
FRANKS=8                                  # cube3 Fist/NNP rank split, as round 1

rundir="$RUN_ROOT/RPA/$CONC/$PSUB/$CELL/seg$SEG"
restart="$rundir/NaCl_RPA_${CELL}_seg${SEG}-1.restart"
[ -f "$restart" ] || { echo "FATAL: no checkpoint to extend: $restart" >&2; exit 98; }
start0=$(awk '$1=="STEP_START_VAL"{print $2; exit}' "$restart")
added=$(( 320000 - ${start0:-0} )); [ "$added" -lt 0 ] && added=0

echo "node: $(hostname)  conc: $CONC  segment: $SEG  extend 80->160 ps (from step ${start0:-?}, +$added steps)"
t0=$SECONDS
MODEL=RPA CONC_DIR="$CONC" PROD_SUBDIR="$PSUB" CELLS="$CELL" SEGMENTS="$SEG" \
  RUN_ROOT="$RUN_ROOT" TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=$FRANKS \
  PROD_PS=160 EXTEND=1 \
  bash "$ADIR/run_production.sh"
wall=$(( SECONDS - t0 ))

ps=NA
[ "$added" -gt 0 ] && ps=$(awk -v w="$wall" -v n="$added" 'BEGIN{printf "%.4f", w/n}')
bash "$TDIR/log_timing.sh" cpu_ext_seg$SEG "$CONC" "$BIN_LABEL" "$SLURM_NTASKS" "$added" "$wall" "$ps" \
  "$CELL RPA CPU extend to 160ps NVE; steps col = steps actually added"
echo "CPU EXTEND SEGMENT DONE ($CONC seg$SEG), wall ${wall}s"
