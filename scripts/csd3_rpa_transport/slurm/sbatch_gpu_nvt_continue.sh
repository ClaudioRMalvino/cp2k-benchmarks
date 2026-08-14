#!/usr/bin/env bash
#SBATCH -J rpa_gpu_nvt_ext
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_nvt_ext_%j.out

# ROUND 2 P2 stage 1: continue a finished NVT equil 95 ps (step 210000 -> 400000),
# RESTART_HISTORY retargeted 30000 -> 40000 to give 5 new NVE starts 20 ps apart
# (snapshot_6..10 at steps 240000..400000; snapshot_5 = 210000 seeded round-1 seg5).
# &RESTART checkpoint cadence stays 30000. Idempotent on resume.
#   sbatch -A mphil-nikiforakis-crm98-sl2-gpu sbatch_gpu_nvt_continue.sh cubic_4m
set -euo pipefail
CONC=${1:?usage: sbatch sbatch_gpu_nvt_continue.sh <conc_dir>}
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport
source "$TDIR/env_csd3_gpu.sh"

RUN_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs_gpu
rundir="$RUN_ROOT/RPA/$CONC/equil/cube3"
proj="NaCl_RPA_cube3_equil"
restart="$rundir/${proj}-1.restart"
NEW_TOTAL=400000          # 210000 round-1 NVT + 190000 continuation (95 ps)
SNAP_STEPS_NEW=40000      # 20 ps snapshot spacing

[ -f "$restart" ] || { echo "no $restart - nothing to continue" >&2; exit 98; }

sed -i "/&RESTART_HISTORY/,/&END RESTART_HISTORY/ s/^\( *\)MD 30000$/\1MD $SNAP_STEPS_NEW/" "$restart"
awk '/&RESTART_HISTORY/,/&END RESTART_HISTORY/' "$restart" | grep -q "MD $SNAP_STEPS_NEW" \
  || { echo "RESTART_HISTORY retarget failed - aborting" >&2; exit 99; }

sed -i -E 's/^( *RESTART_COUNTERS +)[FT.]+/\1T/' "$restart"
start=$(awk '$1=="STEP_START_VAL"{print $2; exit}' "$restart")
remaining=$(( NEW_TOTAL - ${start:-0} )); [ "$remaining" -lt 0 ] && remaining=0
sed -i "0,/^\( *\)STEPS  *[0-9][0-9]*$/s//\1STEPS $remaining/" "$restart"
got=$(awk '$1=="STEPS"{print $2; exit}' "$restart")
[ "$got" = "$remaining" ] || { echo "STEPS cap failed ($got != $remaining) - aborting" >&2; exit 99; }
echo "NVT continuation $CONC: from step ${start:-0}, $remaining steps remain (target $NEW_TOTAL)"

echo "node: $(hostname)"
nvidia-smi -L 2>/dev/null || true
t0=$SECONDS
ended0=$(grep -c "PROGRAM ENDED" "$rundir/equil.out" 2>/dev/null || :)
( cd "$rundir" && \
  srun --ntasks="$SLURM_NTASKS" --cpus-per-task=1 --hint=nomultithread \
       "$BIN" -i "${proj}-1.restart" -o equil.out )
ended1=$(grep -c "PROGRAM ENDED" "$rundir/equil.out" 2>/dev/null || :)
[ "${ended1:-0}" -gt "${ended0:-0}" ] || { echo "did not finish - resubmit to continue" >&2; exit 1; }
wall=$(( SECONDS - t0 ))

k=6
for s in 240000 280000 320000 360000 400000; do
  f="$rundir/${proj}-1_${s}.restart"
  [ -f "$f" ] || { echo "missing history file $f" >&2; exit 1; }
  cp "$f" "$rundir/snapshot_$k.restart"
  k=$(( k + 1 ))
done
echo "staged snapshot_6..10 (steps 240000..400000, 20 ps apart)"
bash "$TDIR/log_timing.sh" gpu_nvt_ext "$CONC" gpu-collab-branch "$SLURM_NTASKS" "$remaining" "$wall" NA "cube3 RPA GPU NVT continuation 210k->400k, snapshot_6..10 for round-2 segments"
echo "NVT CONTINUATION DONE ($CONC), wall ${wall}s"
