#!/usr/bin/env bash
set -euo pipefail

# STAGE 1 - NVT (CSVR) equilibration of the cubic 1.0 mol/kg NaCl(aq) anchor
# cells on CSD3. Port of the cerberus run_equil.sh: 30 ps equilibration +
# 5 restart snapshots 15 ps apart (shared NVE starting points).
# Cube cells only:
#   cube2 = 24.84 A, 1500 atoms,  9 NaCl pairs
#   cube3 = 37.26 A, 5064 atoms, 30 NaCl pairs
# Must run inside an sbatch allocation (uses srun) - see slurm/ wrappers.
#
# Env:  MODEL=MP2  CELLS="cube2"  TOTAL_RANKS=76  SKIP_DONE=1
#       STEPS_OVERRIDE (smoke test only)
#       CONC_DIR=cubic_1M (concentration subdir under runs/$MODEL/; the
#       transport campaign uses cubic_2m / cubic_4m with their own cube_n3.xyz)
#       SEED_OVERRIDE (optional int): non-default &GLOBAL SEED, so a redo of
#       an already-equilibrated cell is statistically independent of the
#       original run instead of bit-identical to it

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="$SCRIPT_DIR/nacl_diffusion_template.inp"
# ENV_SCRIPT override: the GPU track sources the ampere env instead
source "${ENV_SCRIPT:-$SCRIPT_DIR/env_csd3.sh}"

ANCHOR_ROOT="${ANCHOR_ROOT:-/rds/user/$USER/hpc-work/nacl_mp2_anchor}"
NACL_REPO="${NACL_REPO:-$ANCHOR_ROOT}"
RUN_ROOT="${RUN_ROOT:-$ANCHOR_ROOT/runs}"

MODEL="${MODEL:-MP2}"
CELLS="${CELLS:-cube2}"
TOTAL_RANKS="${TOTAL_RANKS:-76}"

# figS4 stage-1 parameters (0.5 fs timestep)
EQUIL_PS=30
N_SNAPSHOTS=5
SNAP_SPACING_PS=15
SNAP_STEPS=$(( SNAP_SPACING_PS * 2000 ))                       # 30000
TOTAL_STEPS=$(( (EQUIL_PS + N_SNAPSHOTS * SNAP_SPACING_PS) * 2000 ))  # 210000

# NNP/Fist rank split of the MIXED groups. 1 Fist rank (cerberus default)
# bottlenecks on CSD3's 76 fast ranks -> tune with the split-sweep job and
# override here. Performance only; the physics is split-independent.
FIST_RANKS="${FIST_RANKS_OVERRIDE:-1}"
NNP_RANKS=$(( TOTAL_RANKS - FIST_RANKS ))

MODEL_DIR="$NACL_REPO/models/$MODEL/final_model"
[ -d "$MODEL_DIR" ] || { echo "no such model dir: $MODEL_DIR" >&2; exit 1; }
CUBE_DIR="$RUN_ROOT/$MODEL/${CONC_DIR:-cubic_1M}"

for cell in $CELLS; do
  case "$cell" in
    cube*) ;;
    *) echo "only cube cells are staged on CSD3 (got '$cell'); tiled digit" \
           "cells live in the cerberus campaign" >&2; exit 1 ;;
  esac
  # pre-built full cubic cell at 1.0 mol/kg (side 12.42*N, mult 1x1x1)
  n=${cell#cube}
  mx=1; my=1; mz=1
  abc=$(awk "BEGIN{printf \"%.4f\", 12.42*$n}")
  gx=$(( 40 * n ))
  coord_cell="$CUBE_DIR/cube_n${n}.xyz"
  [ -f "$coord_cell" ] || { echo "missing $coord_cell" >&2; exit 1; }
  rundir="$CUBE_DIR/equil/$cell"
  proj="NaCl_${MODEL}_${cell}_equil"

  if [ "${SKIP_DONE:-1}" = "1" ] && [ -f "$rundir/snapshot_$N_SNAPSHOTS.restart" ]; then
    echo "$cell: snapshots exist, skipping (SKIP_DONE=0 to force)"
    continue
  fi
  mkdir -p "$rundir"

  sed -e "s|__PROJECT__|$proj|" \
      -e "s|__STEPS__|${STEPS_OVERRIDE:-$TOTAL_STEPS}|" \
      -e "s|__SNAPSTEPS__|$SNAP_STEPS|" \
      -e "s|__MODEL__|$MODEL_DIR|" \
      -e "s|__COORD__|$coord_cell|" \
      -e "s|__ABCX__|$abc|" -e "s|__ABCY__|$abc|" -e "s|__ABCZ__|$abc|" \
      -e "s|__MX__|$mx|" -e "s|__MY__|$my|" -e "s|__MZ__|$mz|" \
      -e "s|__GMAXX__|$gx|" \
      -e "s|__GMAXY__|$gx|" \
      -e "s|__GMAXZ__|$gx|" \
      -e "s|__NNPRANKS__|$NNP_RANKS|" \
      -e "s|__FISTRANKS__|$FIST_RANKS|" \
      "$TEMPLATE" > "$rundir/equil.inp"

  if [ -n "${SEED_OVERRIDE:-}" ]; then
    sed -i "s|^&GLOBAL$|\&GLOBAL\n  SEED $SEED_OVERRIDE|" "$rundir/equil.inp"
  fi

  # GPU track: inject the USE_GPU keyword into &NNP (ON aborts without a device)
  if [ -n "${USE_GPU_OVERRIDE:-}" ]; then
    sed -i "s|^  &NNP$|  \&NNP\n    USE_GPU $USE_GPU_OVERRIDE|" "$rundir/equil.inp"
  fi

  # 12 h walltime survival: if a previous attempt was killed mid-run, continue
  # from its last checkpoint instead of starting over.
  if grep -q "PROGRAM ENDED" "$rundir/equil.out" 2>/dev/null; then
    echo "$cell: MD already finished, staging snapshots"
  else
    INPUT=equil.inp
    if [ -f "$rundir/${proj}-1.restart" ]; then
      INPUT="${proj}-1.restart"
      # CP2K's MD STEPS is per-invocation (STEP_START_VAL only offsets the
      # numbering), so a resume running the restart verbatim would do a full
      # TOTAL_STEPS again, blow the walltime, and never reach the snapshot
      # staging below. Cap STEPS to the remaining budget.
      start=$(awk '$1=="STEP_START_VAL"{print $2; exit}' "$rundir/$INPUT")
      remaining=$(( ${STEPS_OVERRIDE:-$TOTAL_STEPS} - ${start:-0} ))
      [ "$remaining" -lt 0 ] && remaining=0
      sed -i "0,/^\( *\)STEPS  *[0-9][0-9]*$/s//\1STEPS $remaining/" "$rundir/$INPUT"
      # fail-closed: abort (spending nothing) rather than launch MD with an
      # uncapped STEPS if the rewrite did not take
      got=$(awk '$1=="STEPS"{print $2; exit}' "$rundir/$INPUT")
      [ "$got" = "$remaining" ] || { echo "$cell: STEPS cap failed ($got != $remaining) - aborting" >&2; exit 99; }
      echo "$cell: continuing from checkpoint $INPUT (step ${start:-0}, $remaining steps remain)"
    fi
    echo "=== equil $cell (${abc} A cubic), ${STEPS_OVERRIDE:-$TOTAL_STEPS} steps, $TOTAL_RANKS ranks ==="
    ( cd "$rundir" && \
      srun --ntasks="$TOTAL_RANKS" --cpus-per-task=1 --hint=nomultithread \
           "$BIN" -i "$INPUT" -o equil.out )
    if ! grep -q "PROGRAM ENDED" "$rundir/equil.out"; then
      echo "$cell: run did not finish (walltime?) - resubmit to continue from checkpoint" >&2
      exit 1
    fi
  fi

  # last N_SNAPSHOTS restart-history files -> snapshot_k.restart
  mapfile -t snaps < <(ls -v "$rundir/${proj}-1_"*.restart 2>/dev/null | tail -n "$N_SNAPSHOTS")
  if [ "${#snaps[@]}" -ne "$N_SNAPSHOTS" ]; then
    echo "$cell: expected $N_SNAPSHOTS restart snapshots, found ${#snaps[@]}" >&2
    exit 1
  fi
  k=1
  for s in "${snaps[@]}"; do
    cp "$s" "$rundir/snapshot_$k.restart"
    k=$(( k + 1 ))
  done
  echo "$cell: staged $N_SNAPSHOTS snapshots"
done

echo "Equilibration complete for: $CELLS (model $MODEL)"
