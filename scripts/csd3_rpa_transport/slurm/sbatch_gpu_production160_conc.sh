#!/usr/bin/env bash
#SBATCH -J rpa_gpu_p160
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_p160_%A_%a.out

# ROUND 2 P2 stage 2: one FRESH 160 ps NVE segment per array task on one
# A100, from the round-2 NVT snapshots (snapshot_6..10). ~15.5 h at healthy
# pace vs the 12 h wall, so the resume twin is expected to finish the last
# ~40 ps via the capped-STEPS checkpoint logic (that is the plan, not an
# accident - 12 h primaries backfill far better than 17 h one-shots).
#
#   p=$(sbatch --parsable -A nikiforakis-csc-funds-sl3-gpu --array=6-10 \
#         --dependency=afterok:<nvt_job> sbatch_gpu_production160_conc.sh cubic_4m)
#   sbatch -A nikiforakis-csc-funds-sl3-gpu --array=6-10 --dependency=afterany:$p \
#         sbatch_gpu_production160_conc.sh cubic_4m
set -euo pipefail
CONC=${1:?usage: sbatch --array=6-10 sbatch_gpu_production160_conc.sh <conc_dir>}
SEG=${SLURM_ARRAY_TASK_ID:?submit with --array}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport

export ENV_SCRIPT="$TDIR/env_csd3_gpu.sh"
export USE_GPU_OVERRIDE=ON
RUN_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs_gpu

echo "node: $(hostname)  conc: $CONC  segment: $SEG  track: GPU 160 ps NVE (round 2)"
nvidia-smi -L 2>/dev/null || true
t0=$SECONDS
MODEL=RPA CONC_DIR="$CONC" CELLS=cube3 SEGMENTS="$SEG" RUN_ROOT="$RUN_ROOT" \
  TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=3 PROD_PS=160 \
  bash "$ADIR/run_production.sh"
wall=$(( SECONDS - t0 ))

out="$RUN_ROOT/RPA/$CONC/production/cube3/seg$SEG/prod.out"
tstep=$(awk '/qs_mol_dyn_low/ {v=$NF} END{print v}' "$out" || true)
bash "$TDIR/log_timing.sh" gpu_p160_seg$SEG "$CONC" gpu-collab-branch "$SLURM_NTASKS" 320000 "$wall" "${tstep:-NA}" "cube3 RPA GPU 160ps NVE round-2 segment; col12 = qs_mol_dyn_low wall of last window"
echo "GPU 160ps SEGMENT DONE ($CONC seg$SEG), wall ${wall}s"
