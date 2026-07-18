#!/usr/bin/env bash
set -euo pipefail

# STAGE 2 - NVE production for the NaCl(aq) diffusion study on cerberus.
# Mirrors scripts/CSD3_benchmark_scripts/figS4/run_nnp_figS4_production_*.sh:
# 100 ps NVE per segment, 5 independent segments started from the stage-1
# snapshots (positions + velocities only). Stress every step, positions
# every 10 steps - exactly what compute_viscosity.py / compute_diffusion.py
# expect.
#
# Usage:   ./run_production.sh
# Env:     MODEL=revPBE-D3   CELLS="111 211 221 222"   SEGMENTS="1 2 3 4 5"
#          TOTAL_RANKS=48    PROD_PS=100

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="$SCRIPT_DIR/nacl_diffusion_template.inp"

SCRATCH_ROOT="${SCRATCH_ROOT:-/data/cerberus1}"
SCRATCH_USER="crm98"
CP2K_ROOT="${CP2K_ROOT:-${SCRATCH_ROOT}/${SCRATCH_USER}/cp2k_optimized}"
BIN="${BIN:-$CP2K_ROOT/install/bin/cp2k.psmp}"
NACL_REPO="${NACL_REPO:-${SCRATCH_ROOT}/${SCRATCH_USER}/nacl-water}"
RUN_ROOT="${RUN_ROOT:-${SCRATCH_ROOT}/${SCRATCH_USER}/nacl_diffusion}"

MODEL="${MODEL:-revPBE-D3}"
CELLS="${CELLS:-111 211 221 222}"
SEGMENTS="${SEGMENTS:-1 2 3 4 5}"
N_PAIRS="${N_PAIRS:-1}"            # NaCl pairs per base cell (0 = pure water)
TOTAL_RANKS="${TOTAL_RANKS:-32}"
PROD_PS="${PROD_PS:-100}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

PROD_STEPS=$(( PROD_PS * 2000 ))   # 0.5 fs timestep
RESTART_EVERY=10000                # checkpoint every 5 ps (6 h walltime survival)
NNP_RANKS=$(( TOTAL_RANKS - 1 ))
FIST_RANKS=1

set +u
source "$CP2K_ROOT/tools/toolchain/install/setup"
set -u
export LD_LIBRARY_PATH="$CP2K_ROOT/install/lib:${LD_LIBRARY_PATH:-}"

MODEL_DIR="$NACL_REPO/models/$MODEL/final_model"
WORK="$RUN_ROOT/$MODEL/np${N_PAIRS}"
CUBE_DIR="$RUN_ROOT/$MODEL/cubic_1M"
COORD="$WORK/base_cell.xyz"

for cell in $CELLS; do
  if [[ "$cell" == cube* ]]; then
    n=${cell#cube}
    mx=1; my=1; mz=1
    abc=$(awk "BEGIN{printf \"%.4f\", 12.42*$n}")
    abcx=$abc; abcy=$abc; abcz=$abc
    gx=$(( 40 * n )); gy=$gx; gz=$gx
    coord_cell="$CUBE_DIR/cube_n${n}.xyz"
    equildir="$CUBE_DIR/equil/$cell"
    pbase="NaCl_${MODEL}_${cell}"
    prod_parent="$CUBE_DIR/production/$cell"
  else
    mx=${cell:0:1}; my=${cell:1:1}; mz=${cell:2:1}
    abcx=12.42; abcy=12.42; abcz=12.42
    gx=$(( 40 * mx )); gy=$(( 40 * my )); gz=$(( 40 * mz ))
    coord_cell="$COORD"
    equildir="$WORK/equil/cell$cell"
    pbase="NaCl_${MODEL}_np${N_PAIRS}_cell${cell}"
    prod_parent="$WORK/production/cell$cell"
  fi

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
        -e "s|__ABCX__|$abcx|" -e "s|__ABCY__|$abcy|" -e "s|__ABCZ__|$abcz|" \
        -e "s|__MX__|$mx|" -e "s|__MY__|$my|" -e "s|__MZ__|$mz|" \
        -e "s|__GMAXX__|$gx|" \
        -e "s|__GMAXY__|$gy|" \
        -e "s|__GMAXZ__|$gz|" \
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

    # 6 h walltime survival: continue a killed run from its last checkpoint
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
      echo "cell$cell seg$seg: continuing from checkpoint $INPUT (counters restored)"
    fi
    echo "=== production cell$cell seg$seg: $PROD_STEPS steps ($PROD_PS ps) ==="
    ( cd "$rundir" && \
      mpirun -n "$TOTAL_RANKS" --bind-to core "$BIN" -i "$INPUT" -o prod.out )
    if ! grep -q "PROGRAM ENDED" "$rundir/prod.out"; then
      echo "cell$cell seg$seg: did not finish (walltime?) - rerun to continue" >&2
      exit 1
    fi
  done
done

echo "Production complete: model=$MODEL cells=$CELLS segments=$SEGMENTS"
