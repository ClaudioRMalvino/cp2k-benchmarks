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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="$SCRIPT_DIR/nacl_diffusion_template.inp"
source "$SCRIPT_DIR/env_csd3.sh"

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
CUBE_DIR="$RUN_ROOT/$MODEL/cubic_1M"

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
  prod_parent="$CUBE_DIR/production/$cell"

  for seg in $SEGMENTS; do
    snap="$equildir/snapshot_$seg.restart"
    [ -f "$snap" ] || { echo "missing $snap - run run_equil.sh first" >&2; exit 1; }

    rundir="$prod_parent/seg$seg"
    proj="${pbase}_seg${seg}"
    if [ "${SKIP_DONE:-1}" = "1" ] && [ -f "$rundir/${proj}-1.stress" ] && \
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

    # positions + velocities from the snapshot; everything else fresh
    cat >> "$rundir/prod.inp" <<EOF

&EXT_RESTART
  RESTART_FILE_NAME $snap
  RESTART_DEFAULT .FALSE.
  RESTART_POS .TRUE.
  RESTART_VEL .TRUE.
&END EXT_RESTART
EOF

    # 12 h walltime survival: continue a killed run from its last checkpoint
    INPUT=prod.inp
    if [ -f "$rundir/${proj}-1.restart" ] && \
       ! grep -q "PROGRAM ENDED" "$rundir/prod.out" 2>/dev/null; then
      INPUT="${proj}-1.restart"
      # Restore the step counter on continuation so the run stops at the total
      # STEPS instead of doing a whole fresh segment. The initial launch keeps
      # counters OFF (production must start at step 0 from the equil snapshot);
      # a continuation must turn them ON, otherwise CP2K treats STEPS as a
      # per-invocation count and overruns to 2x the segment length.
      sed -i -E 's/^( *RESTART_COUNTERS +)[FT.]+/\1T/' "$rundir/$INPUT"
      echo "$cell seg$seg: continuing from checkpoint $INPUT (counters restored)"
    fi
    echo "=== production $cell seg$seg: $PROD_STEPS steps ($PROD_PS ps) ==="
    ( cd "$rundir" && \
      srun --ntasks="$TOTAL_RANKS" --cpus-per-task=1 --hint=nomultithread \
           "$BIN" -i "$INPUT" -o prod.out )
    if ! grep -q "PROGRAM ENDED" "$rundir/prod.out"; then
      echo "$cell seg$seg: did not finish (walltime?) - resubmit to continue" >&2
      exit 1
    fi
  done
done

echo "Production complete: model=$MODEL cells=$CELLS segments=$SEGMENTS"
