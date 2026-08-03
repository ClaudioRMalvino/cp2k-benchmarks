#!/usr/bin/env bash
#SBATCH -J rpa_gpu_equil
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_equil_%j.out

# GPU-track equilibration of one concentration's cube3 cell: same protocol
# (30 ps CSVR + 5 snapshots @15 ps = 210k steps) via the same run_equil.sh,
# but on one A100 (USE_GPU ON, 1 NNP rank on the device + 3 CPU SPME ranks)
# and into the separate runs_gpu/ tree so CPU and GPU runs never collide.
# Account at submit: sbatch -A <SL3-GPU acct> sbatch_gpu_equil_conc.sh cubic_1M
# Adjust --time at submit once the smoke's measured s/step is known.
set -euo pipefail
CONC=${1:?usage: sbatch sbatch_gpu_equil_conc.sh <conc_dir e.g. cubic_1M>}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport

export ENV_SCRIPT="$TDIR/env_csd3_gpu.sh"
export USE_GPU_OVERRIDE=ON
RUN_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs_gpu

echo "node: $(hostname)  conc: $CONC  track: GPU (1x A100 + 3 SPME ranks)"
nvidia-smi -L 2>/dev/null || true
t0=$SECONDS
MODEL=RPA CONC_DIR="$CONC" CELLS=cube3 RUN_ROOT="$RUN_ROOT" \
  TOTAL_RANKS=$SLURM_NTASKS FIST_RANKS_OVERRIDE=3 \
  bash "$ADIR/run_equil.sh"
wall=$(( SECONDS - t0 ))

out="$RUN_ROOT/RPA/$CONC/equil/cube3/equil.out"
tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' "$out" || true)
ps=""
[ -n "${tstep:-}" ] && ps=$(awk -v t="$tstep" 'BEGIN{printf "%.4f", t/210000}')
bash "$TDIR/log_timing.sh" gpu_equil "$CONC" gpu-collab-branch "$SLURM_NTASKS" 210000 "$wall" "${ps:-NA}" "cube3 RPA GPU NVT+snapshots; s_per_step valid only single-window"
echo "GPU EQUIL DONE ($CONC), wall ${wall}s"
