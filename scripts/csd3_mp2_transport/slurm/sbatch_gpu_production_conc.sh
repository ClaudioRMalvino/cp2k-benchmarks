#!/usr/bin/env bash
#SBATCH -J rpa_gpu_prod
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_prod_%A_%a.out

# GPU-track production: one 80 ps NVE segment per array task on one A100
# (USE_GPU ON, 1 NNP device rank + 3 CPU SPME ranks), from the GPU-track
# equil snapshots in runs_gpu/. Same run_production.sh, same checkpoints,
# same RESTART_COUNTERS resume logic.
#
#   p=$(sbatch --parsable -A <SL3-GPU> --array=1-5 sbatch_gpu_production_conc.sh cubic_1M)
#   sbatch -A <SL3-GPU> --array=1-5 --dependency=afterany:$p sbatch_gpu_production_conc.sh cubic_1M
# Adjust --time at submit once the smoke's measured s/step is known.
set -euo pipefail
CONC=${1:?usage: sbatch --array=1-5 sbatch_gpu_production_conc.sh <conc_dir>}
SEG=${SLURM_ARRAY_TASK_ID:?submit with --array=1-5}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_mp2_transport

export ENV_SCRIPT="$TDIR/env_csd3_gpu.sh"
export USE_GPU_OVERRIDE=ON
RUN_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs_gpu

echo "node: $(hostname)  conc: $CONC  segment: $SEG  track: GPU (1x A100 + 3 SPME ranks)"
nvidia-smi -L 2>/dev/null || true
t0=$SECONDS
MODEL=RPA CONC_DIR="$CONC" CELLS=cube3 SEGMENTS="$SEG" RUN_ROOT="$RUN_ROOT" \
  TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=3 PROD_PS=80 \
  bash "$ADIR/run_production.sh"
wall=$(( SECONDS - t0 ))

out="$RUN_ROOT/RPA/$CONC/production/cube3/seg$SEG/prod.out"
tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' "$out" || true)
ps=""
[ -n "${tstep:-}" ] && ps=$(awk -v t="$tstep" 'BEGIN{printf "%.4f", t/160000}')
bash "$TDIR/log_timing.sh" gpu_prod_seg$SEG "$CONC" gpu-collab-branch "$SLURM_NTASKS" 160000 "$wall" "${ps:-NA}" "cube3 RPA GPU 80ps NVE; s_per_step valid only single-window"
echo "GPU PRODUCTION SEGMENT DONE ($CONC seg$SEG), wall ${wall}s"
