#!/usr/bin/env bash
set -euo pipefail

# STAGE 2 - NVE production for the cubic NaCl(aq) MP2 anchor on CSD3.
# Port of the cerberus run_production.sh: 100 ps NVE per segment, 5
# independent segments started from the stage-1 snapshots (positions +
# velocities only). Stress every step, positions every 10 steps - exactly
# what compute_viscosity.py / compute_diffusion.py expect.
# Must run inside an sbatch allocation (uses srun) - see slurm/ wrappers.
#
# Env:  MODEL=MP2  CELLS="cube2"  SEGMENTS="1 2 3 4 5"  TOTAL_RANKS=76
#       PROD_PS=100
#       CONC_DIR=cubic_1M (concentration subdir; transport campaign uses
#       cubic_2m / cubic_4m)
#       PROD_SUBDIR=production (output subdir; the 1 m transport pilot uses
#       production_transport so the completed anchor segments stay untouched)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="$SCRIPT_DIR/nacl_diffusion_template.inp"
# ENV_SCRIPT override: the GPU track sources the ampere env instead
source "${ENV_SCRIPT:-$SCRIPT_DIR/env_csd3.sh}"

ANCHOR_ROOT="${ANCHOR_ROOT:-/rds/user/$USER/hpc-work/nacl_mp2_anchor}"
NACL_REPO="${NACL_REPO:-$ANCHOR_ROOT}"
RUN_ROOT="${RUN_ROOT:-$ANCHOR_ROOT/runs}"

MODEL="${MODEL:-MP2}"
CELLS="${CELLS:-cube2}"
SEGMENTS="${SEGMENTS:-1 2 3 4 5}"
TOTAL_RANKS="${TOTAL_RANKS:-76}"
PROD_PS="${PROD_PS:-100}"

PROD_STEPS=$(( PROD_PS * 2000 ))   # 0.5 fs timestep
RESTART_EVERY=10000                # checkpoint every 5 ps (12 h walltime survival)
# see run_equil.sh: Fist ranks tuned via split-sweep; physics-neutral
FIST_RANKS="${FIST_RANKS_OVERRIDE:-1}"
NNP_RANKS=$(( TOTAL_RANKS - FIST_RANKS ))

MODEL_DIR="$NACL_REPO/models/$MODEL/final_model"
[ -d "$MODEL_DIR" ] || { echo "no such model dir: $MODEL_DIR" >&2; exit 1; }
CUBE_DIR="$RUN_ROOT/$MODEL/${CONC_DIR:-cubic_1M}"

for cell in $CELLS; do
  case "$cell" in
    cube*) ;;
    *) echo "only cube cells are staged on CSD3 (got '$cell')" >&2; exit 1 ;;
  esac
  n=${cell#cube}
  mx=1; my=1; mz=1
  abc=$(awk "BEGIN{printf \"%.4f\", 12.42*$n}")
  gx=$(( 40 * n ))
  coord_cell="$CUBE_DIR/cube_n${n}.xyz"
  equildir="$CUBE_DIR/equil/$cell"
  pbase="NaCl_${MODEL}_${cell}"
  prod_parent="$CUBE_DIR/${PROD_SUBDIR:-production}/$cell"

  for seg in $SEGMENTS; do
    snap="$equildir/snapshot_$seg.restart"
    [ -f "$snap" ] || { echo "missing $snap - run run_equil.sh first" >&2; exit 1; }

    rundir="$prod_parent/seg$seg"
    proj="${pbase}_seg${seg}"
    # EXTEND=1 (round-2): lengthen an already-finished segment to the new
    # PROD_PS by continuing from its final checkpoint - so no skip-if-done,
    # and the resume branch below must fire even though prod.out says ENDED.
    if [ "${EXTEND:-0}" != "1" ] && [ "${SKIP_DONE:-1}" = "1" ] && \
       [ -f "$rundir/${proj}-1.stress" ] && \
       grep -q "PROGRAM ENDED" "$rundir/prod.out" 2>/dev/null; then
      echo "$cell seg$seg: already complete, skipping"
      continue
    fi
    mkdir -p "$rundir"

    sed -e "s|__PROJECT__|$proj|" \
        -e "s|__STEPS__|$PROD_STEPS|" \
        -e "s|__SNAPSTEPS__|$RESTART_EVERY|" \
        -e "s|__MODEL__|$MODEL_DIR|" \
        -e "s|__COORD__|$coord_cell|" \
        -e "s|__ABCX__|$abc|" -e "s|__ABCY__|$abc|" -e "s|__ABCZ__|$abc|" \
        -e "s|__MX__|$mx|" -e "s|__MY__|$my|" -e "s|__MZ__|$mz|" \
        -e "s|__GMAXX__|$gx|" \
        -e "s|__GMAXY__|$gx|" \
        -e "s|__GMAXZ__|$gx|" \
        -e "s|__NNPRANKS__|$NNP_RANKS|" \
        -e "s|__FISTRANKS__|$FIST_RANKS|" \
        -e "s|ENSEMBLE NVT|ENSEMBLE NVE|" \
        "$TEMPLATE" \
      | sed -e '/&THERMOSTAT/,/&END THERMOSTAT/d' \
            -e 's|^\( *\)&RESTART_HISTORY$|\1\&RESTART_HISTORY OFF|' > "$rundir/prod.inp"

    # GPU track: inject the USE_GPU keyword into &NNP
    if [ -n "${USE_GPU_OVERRIDE:-}" ]; then
      sed -i "s|^  &NNP$|  \&NNP\n    USE_GPU $USE_GPU_OVERRIDE|" "$rundir/prod.inp"
    fi

    # positions + velocities from the snapshot; everything else fresh
    cat >> "$rundir/prod.inp" <<EOF

&EXT_RESTART
  RESTART_FILE_NAME $snap
  RESTART_DEFAULT .FALSE.
  RESTART_POS .TRUE.
  RESTART_VEL .TRUE.
&END EXT_RESTART
EOF

    # 12 h walltime survival: continue a killed run from its last checkpoint.
    # EXTEND=1 additionally treats a FINISHED segment's checkpoint as the
    # input; fail closed if there is nothing to continue from, because the
    # fallthrough would re-run the whole segment from the equil snapshot.
    INPUT=prod.inp
    if [ "${EXTEND:-0}" = "1" ] && [ ! -f "$rundir/${proj}-1.restart" ]; then
      echo "$cell seg$seg: EXTEND=1 but no ${proj}-1.restart - aborting" >&2
      exit 98
    fi
    if [ -f "$rundir/${proj}-1.restart" ] && \
       { [ "${EXTEND:-0}" = "1" ] || \
         ! grep -q "PROGRAM ENDED" "$rundir/prod.out" 2>/dev/null; }; then
      INPUT="${proj}-1.restart"
      # Counters ON keeps the step numbering continuous across the resume;
      # the initial launch keeps them OFF (production starts at step 0 from
      # the equil snapshot). Numbering is ALL counters fix: CP2K's MD STEPS
      # is per-invocation even with STEP_START_VAL set, so STEPS must also
      # be capped to the remaining budget or the resume runs a full extra
      # segment past the target.
      sed -i -E 's/^( *RESTART_COUNTERS +)[FT.]+/\1T/' "$rundir/$INPUT"
      start=$(awk '$1=="STEP_START_VAL"{print $2; exit}' "$rundir/$INPUT")
      remaining=$(( PROD_STEPS - ${start:-0} ))
      [ "$remaining" -lt 0 ] && remaining=0
      sed -i "0,/^\( *\)STEPS  *[0-9][0-9]*$/s//\1STEPS $remaining/" "$rundir/$INPUT"
      # fail-closed: abort (spending nothing) rather than launch MD with an
      # uncapped STEPS if the rewrite did not take
      got=$(awk '$1=="STEPS"{print $2; exit}' "$rundir/$INPUT")
      [ "$got" = "$remaining" ] || { echo "$cell seg$seg: STEPS cap failed ($got != $remaining) - aborting" >&2; exit 99; }
      echo "$cell seg$seg: continuing from checkpoint $INPUT (step ${start:-0}, $remaining steps remain)"
    fi
    echo "=== production $cell seg$seg: $PROD_STEPS steps ($PROD_PS ps) ==="
    # CP2K appends to prod.out on resume/extend, so success = the count of
    # ENDED banners going up, not mere presence (round-1's banner is already
    # in the file when extending)
    ended0=$(grep -c "PROGRAM ENDED" "$rundir/prod.out" 2>/dev/null || :)
    ( cd "$rundir" && \
      srun --ntasks="$TOTAL_RANKS" --cpus-per-task=1 --hint=nomultithread \
           "$BIN" -i "$INPUT" -o prod.out )
    ended1=$(grep -c "PROGRAM ENDED" "$rundir/prod.out" 2>/dev/null || :)
    if [ "${ended1:-0}" -le "${ended0:-0}" ]; then
      echo "$cell seg$seg: did not finish (walltime?) - resubmit to continue" >&2
      exit 1
    fi
  done
done

echo "Production complete: model=$MODEL cells=$CELLS segments=$SEGMENTS"
