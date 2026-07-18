#!/usr/bin/env bash
set -euo pipefail

# STAGE 1 - NVT equilibration of the NaCl(aq) cells on cerberus.
# Mirrors scripts/CSD3_benchmark_scripts/figS4/run_nnp_figS4_equil.sh
# (30 ps equilibration + 5 restart snapshots 15 ps apart, shared NVE
# starting points) but runs directly with mpirun on the cerberus node.
#
# Usage:   ./run_equil.sh
# Env:     MODEL=revPBE-D3   CELLS="111 211 221 222"   TOTAL_RANKS=48
#          SKIP_DONE=1 (skip cells whose snapshots already exist)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="$SCRIPT_DIR/nacl_diffusion_template.inp"

SCRATCH_ROOT="${SCRATCH_ROOT:-/data/cerberus1}"
SCRATCH_USER="crm98"
CP2K_ROOT="${CP2K_ROOT:-${SCRATCH_ROOT}/${SCRATCH_USER}/cp2k_optimized}"
BIN="${BIN:-$CP2K_ROOT/install/bin/cp2k.psmp}"
NACL_REPO="${NACL_REPO:-${SCRATCH_ROOT}/${SCRATCH_USER}/nacl-water}"
RUN_ROOT="${RUN_ROOT:-${SCRATCH_ROOT}/${SCRATCH_USER}/nacl_diffusion}"

MODEL="${MODEL:-revPBE-D3}"
CELLS="${CELLS:-111 211 221 222}"     # tiling of the 12.42 A base cell
N_PAIRS="${N_PAIRS:-1}"               # NaCl pairs per base cell (0 = pure water)
TOTAL_RANKS="${TOTAL_RANKS:-32}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# figS4 stage-1 parameters (0.5 fs timestep)
EQUIL_PS=30
N_SNAPSHOTS=5
SNAP_SPACING_PS=15
SNAP_STEPS=$(( SNAP_SPACING_PS * 2000 ))                       # 30000
TOTAL_STEPS=$(( (EQUIL_PS + N_SNAPSHOTS * SNAP_SPACING_PS) * 2000 ))  # 210000

NNP_RANKS=$(( TOTAL_RANKS - 1 ))
FIST_RANKS=1

# runtime environment for cp2k.psmp
set +u
source "$CP2K_ROOT/tools/toolchain/install/setup"
set -u
export LD_LIBRARY_PATH="$CP2K_ROOT/install/lib:${LD_LIBRARY_PATH:-}"

MODEL_DIR="$NACL_REPO/models/$MODEL/final_model"
[ -d "$MODEL_DIR" ] || { echo "no such model dir: $MODEL_DIR" >&2; exit 1; }

WORK="$RUN_ROOT/$MODEL/np${N_PAIRS}"
mkdir -p "$WORK"
CUBE_DIR="$RUN_ROOT/$MODEL/cubic_1M"
# tiled base cell is only needed for the digit-style cells (111/211/...);
# cubeN cells use a pre-built full cubic cell, so skip the base build (and the
# make_base_cell water frame, which some models e.g. MP2 do not ship)
COORD="$WORK/base_cell.xyz"
if [[ " $CELLS " == *[0-9][0-9][0-9]* ]] && [ ! -f "$COORD" ]; then
  PYTHON="${PYTHON:-${SCRATCH_ROOT}/${SCRATCH_USER}/pyenv/bin/python3}"
  "$PYTHON" "$SCRIPT_DIR/make_base_cell.py" \
    --datasets "$NACL_REPO/data-sets" --model "$MODEL" \
    --n-pairs "$N_PAIRS" --out "$COORD"
fi

for cell in $CELLS; do
  if [[ "$cell" == cube* ]]; then
    # pre-built full cubic cell at 1.0 mol/kg (side 12.42*N, mult 1x1x1)
    n=${cell#cube}
    mx=1; my=1; mz=1
    abc=$(awk "BEGIN{printf \"%.4f\", 12.42*$n}")
    abcx=$abc; abcy=$abc; abcz=$abc
    gx=$(( 40 * n )); gy=$gx; gz=$gx
    coord_cell="$CUBE_DIR/cube_n${n}.xyz"
    [ -f "$coord_cell" ] || { echo "missing $coord_cell (build it with make_base_cell --cubic-mult $n --target-molality 1.0)" >&2; exit 1; }
    rundir="$CUBE_DIR/equil/$cell"
    proj="NaCl_${MODEL}_${cell}_equil"
  else
    mx=${cell:0:1}; my=${cell:1:1}; mz=${cell:2:1}
    abcx=12.42; abcy=12.42; abcz=12.42
    gx=$(( 40 * mx )); gy=$(( 40 * my )); gz=$(( 40 * mz ))
    coord_cell="$COORD"
    rundir="$WORK/equil/cell$cell"
    proj="NaCl_${MODEL}_np${N_PAIRS}_cell${cell}_equil"
  fi
  if [ "${SKIP_DONE:-1}" = "1" ] && [ -f "$rundir/snapshot_$N_SNAPSHOTS.restart" ]; then
    echo "$cell: snapshots exist, skipping (SKIP_DONE=0 to force)"
    continue
  fi
  mkdir -p "$rundir"

  sed -e "s|__PROJECT__|$proj|" \
      -e "s|__STEPS__|$TOTAL_STEPS|" \
      -e "s|__SNAPSTEPS__|$SNAP_STEPS|" \
      -e "s|__MODEL__|$MODEL_DIR|" \
      -e "s|__COORD__|$coord_cell|" \
      -e "s|__ABCX__|$abcx|" -e "s|__ABCY__|$abcy|" -e "s|__ABCZ__|$abcz|" \
      -e "s|__MX__|$mx|" -e "s|__MY__|$my|" -e "s|__MZ__|$mz|" \
      -e "s|__GMAXX__|$gx|" \
      -e "s|__GMAXY__|$gy|" \
      -e "s|__GMAXZ__|$gz|" \
      -e "s|__NNPRANKS__|$NNP_RANKS|" \
      -e "s|__FISTRANKS__|$FIST_RANKS|" \
      "$TEMPLATE" > "$rundir/equil.inp"

  # 6 h walltime survival: if a previous attempt was killed mid-run, continue
  # from its last checkpoint instead of starting over.
  if grep -q "PROGRAM ENDED" "$rundir/equil.out" 2>/dev/null; then
    echo "cell$cell: MD already finished, staging snapshots"
  else
    INPUT=equil.inp
    if [ -f "$rundir/${proj}-1.restart" ]; then
      INPUT="${proj}-1.restart"
      echo "cell$cell: continuing from checkpoint $INPUT"
    fi
    echo "=== equil cell$cell ($mx x $my x $mz), $TOTAL_STEPS steps, $TOTAL_RANKS ranks ==="
    ( cd "$rundir" && \
      mpirun -n "$TOTAL_RANKS" --bind-to core "$BIN" -i "$INPUT" -o equil.out )
    if ! grep -q "PROGRAM ENDED" "$rundir/equil.out"; then
      echo "cell$cell: run did not finish (walltime?) - rerun to continue from checkpoint" >&2
      exit 1
    fi
  fi

  # last N_SNAPSHOTS restart-history files -> snapshot_k.restart
  mapfile -t snaps < <(ls -v "$rundir/${proj}-1_"*.restart 2>/dev/null | tail -n "$N_SNAPSHOTS")
  if [ "${#snaps[@]}" -ne "$N_SNAPSHOTS" ]; then
    echo "cell$cell: expected $N_SNAPSHOTS restart snapshots, found ${#snaps[@]}" >&2
    exit 1
  fi
  k=1
  for s in "${snaps[@]}"; do
    cp "$s" "$rundir/snapshot_$k.restart"
    k=$(( k + 1 ))
  done
  echo "cell$cell: staged $N_SNAPSHOTS snapshots"
done

echo "Equilibration complete for: $CELLS (model $MODEL)"
