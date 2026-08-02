#!/usr/bin/env bash
#SBATCH -J rpa_vacf
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=05:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_vacf_%j.out

# MP2 transport campaign: the ONE short velocity-dump run per concentration
# (VACF spot-check of the MSD-derived diffusivities). 20 ps NVE from
# snapshot_1 with velocities every 4 steps (2 fs). Not part of the MSD/GK
# statistics - purely a cross-method check.
#
#   sbatch sbatch_vacf.sh cubic_1M
set -euo pipefail
CONC=${1:?usage: sbatch sbatch_vacf.sh <conc_dir>}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_mp2_transport
export BIN_LABEL=dhruv-cell-list
source "$ADIR/env_csd3.sh"

MODEL=RPA
ANCHOR_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor
MODEL_DIR="$ANCHOR_ROOT/models/$MODEL/final_model"
CUBE_DIR="$ANCHOR_ROOT/runs/$MODEL/$CONC"
COORD="$CUBE_DIR/cube_n3.xyz"
SNAP="$CUBE_DIR/equil/cube3/snapshot_1.restart"
[ -f "$SNAP" ] || { echo "missing $SNAP - run equil first" >&2; exit 1; }
RUNDIR="$CUBE_DIR/vacf/cube3"
STEPS=40000   # 20 ps at 0.5 fs
RANKS=$SLURM_NTASKS
FIST_RANKS=8
proj="NaCl_${MODEL}_vacf"

echo "node: $(hostname)  conc: $CONC  binary: $BIN"
mkdir -p "$RUNDIR"; rm -f "$RUNDIR"/vacf.out "$RUNDIR/$proj"*

# production-style NVE deck + a MOTION-level VELOCITIES print every 4 steps
sed -e "s|__PROJECT__|$proj|" \
    -e "s|__STEPS__|$STEPS|" \
    -e "s|__SNAPSTEPS__|100000|" \
    -e "s|__MODEL__|$MODEL_DIR|" \
    -e "s|__COORD__|$COORD|" \
    -e "s|__ABCX__|37.2600|" -e "s|__ABCY__|37.2600|" -e "s|__ABCZ__|37.2600|" \
    -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
    -e "s|__GMAXX__|120|" -e "s|__GMAXY__|120|" -e "s|__GMAXZ__|120|" \
    -e "s|__NNPRANKS__|$(( RANKS - FIST_RANKS ))|" \
    -e "s|__FISTRANKS__|$FIST_RANKS|" \
    -e "s|ENSEMBLE NVT|ENSEMBLE NVE|" \
    "$ADIR/nacl_diffusion_template.inp" \
  | sed -e '/&THERMOSTAT/,/&END THERMOSTAT/d' \
        -e 's|^\( *\)&RESTART_HISTORY$|\1\&RESTART_HISTORY OFF|' \
        -e 's|^  &PRINT$|  \&PRINT\n    \&VELOCITIES\n      FORMAT XYZ\n      \&EACH\n        MD 4\n      \&END EACH\n    \&END VELOCITIES|' \
  > "$RUNDIR/vacf.inp"

cat >> "$RUNDIR/vacf.inp" <<EOF

&EXT_RESTART
  RESTART_FILE_NAME $SNAP
  RESTART_DEFAULT .FALSE.
  RESTART_POS .TRUE.
  RESTART_VEL .TRUE.
&END EXT_RESTART
EOF

cd "$RUNDIR"
t0=$SECONDS
srun --ntasks="$RANKS" --cpus-per-task=1 --hint=nomultithread \
  "$BIN" -i vacf.inp -o vacf.out
wall=$(( SECONDS - t0 ))

grep -q "PROGRAM ENDED" vacf.out && echo "OK  cp2k finished" || { echo "FAIL"; exit 1; }
[ -f "$proj-vel-1.xyz" ] && echo "OK  velocity stream present" || { echo "FAIL: no vel file"; exit 1; }
bash "$TDIR/log_timing.sh" vacf "$CONC" "$BIN_LABEL" "$RANKS" "$STEPS" "$wall" \
  "$(awk '/qs_mol_dyn_low/ {printf "%.4f", $NF/'"$STEPS"'; exit}' vacf.out)" "cube3 $MODEL 20ps NVE + velocities"
echo "VACF RUN DONE ($CONC)"
