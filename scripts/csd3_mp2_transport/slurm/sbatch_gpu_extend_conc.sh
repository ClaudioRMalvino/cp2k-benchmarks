#!/usr/bin/env bash
#SBATCH -J rpa_gpu_ext
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_ext_%A_%a.out

# ROUND 2 P1: extend a completed 80 ps NVE segment to 160 ps on one A100.
# EXTEND=1 makes run_production.sh continue from the segment's final
# checkpoint (fail-closed if it is missing - never re-runs from the equil
# snapshot); the capped-STEPS logic computes the remaining budget, so the
# 1M overrun segments (125/155 ps banked) top up with only 35/5 ps. CP2K
# appends pos/stress to the round-1 files; the analyzers monotonic-slice
# the checkpoint overlap - REFRESH the reduced_traj_e*.npz caches before
# re-analysis.
#
#   p=$(sbatch --parsable -A nikiforakis-csc-funds-sl3-gpu --array=1-5 sbatch_gpu_extend_conc.sh cubic_1M)
#   sbatch -A nikiforakis-csc-funds-sl3-gpu --array=1-5 --dependency=afterany:$p sbatch_gpu_extend_conc.sh cubic_1M
set -euo pipefail
CONC=${1:?usage: sbatch --array=1-5 sbatch_gpu_extend_conc.sh <conc_dir>}
SEG=${SLURM_ARRAY_TASK_ID:?submit with --array=1-5}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_mp2_transport

export ENV_SCRIPT="$TDIR/env_csd3_gpu.sh"
export USE_GPU_OVERRIDE=ON
RUN_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs_gpu

rundir="$RUN_ROOT/RPA/$CONC/production/cube3/seg$SEG"
restart="$rundir/NaCl_RPA_cube3_seg${SEG}-1.restart"
start0=$(awk '$1=="STEP_START_VAL"{print $2; exit}' "$restart")
added=$(( 320000 - ${start0:-0} )); [ "$added" -lt 0 ] && added=0

echo "node: $(hostname)  conc: $CONC  segment: $SEG  extend 80->160 ps (from step ${start0:-?}, +$added steps)"
nvidia-smi -L 2>/dev/null || true
t0=$SECONDS
MODEL=RPA CONC_DIR="$CONC" CELLS=cube3 SEGMENTS="$SEG" RUN_ROOT="$RUN_ROOT" \
  TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=3 PROD_PS=160 EXTEND=1 \
  bash "$ADIR/run_production.sh"
wall=$(( SECONDS - t0 ))

ps=NA
[ "$added" -gt 0 ] && ps=$(awk -v w="$wall" -v n="$added" 'BEGIN{printf "%.4f", w/n}')
bash "$TDIR/log_timing.sh" gpu_ext_seg$SEG "$CONC" gpu-collab-branch "$SLURM_NTASKS" "$added" "$wall" "$ps" "cube3 RPA GPU extend to 160ps NVE; steps col = steps actually added"
echo "GPU EXTEND SEGMENT DONE ($CONC seg$SEG), wall ${wall}s"
