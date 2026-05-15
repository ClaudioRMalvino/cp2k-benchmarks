#!/usr/bin/env bash
# ============================================================================
# Fig. S4 replication -- STAGE 1: NVT equilibration + snapshot generation.
#
# SLURM ARRAY JOB: one task per system size (64/128/256/512/1024 H2O).
# Each task runs ONE NVT trajectory at rho ~ 1 g/cm3, 300 K, and saves
# N_SNAPSHOTS decorrelated restart files.  Those snapshots are the shared,
# IDENTICAL starting points for the NVE production segments of BOTH branches
# (upstream master and feature/nnp-native-spline) -- so any difference in the
# resulting viscosity / diffusion is purely numerical, not a difference in
# initial conditions.
#
# Run once (not per branch); the native-spline binary is used purely as the
# config generator since equilibration physics is branch-independent.
#
# Trajectory layout per size:
#   [0, EQUIL_PS)                       -- discarded equilibration
#   restart snapshot every SNAP_SPACING_PS ps, last N_SNAPSHOTS kept
# ============================================================================
#SBATCH -J figS4_equil
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
# Slowest task is N=1024 equilibration (~1.5 h on one node); 3 h is ample.
#SBATCH --time=03:00:00
#SBATCH --array=0-4
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_equil_%A_%a.out

# Source the toolchain environment BEFORE enabling strict mode: the toolchain
# 'setup' script references unbound vars (CP_DFLAGS) that trip `set -u`.
. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCH=/home/crm98/cp2k-benchmarks
RESULTS_SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4
source "$BIN_ROOT/setup"

set -euo pipefail

# --- tunables (override via --export on sbatch) ----------------------------
EQUIL_PS=${EQUIL_PS:-20}                 # discarded equilibration window
N_SNAPSHOTS=${N_SNAPSHOTS:-5}            # independent NVE starting points / size
SNAP_SPACING_PS=${SNAP_SPACING_PS:-10}   # spacing between kept snapshots
TIMESTEP_FS=0.5

# array index -> system size
SIZE_LIST=(64 128 256 512 1024)
size=${SIZE_LIST[$SLURM_ARRAY_TASK_ID]}

# size -> MULTIPLE_UNIT_CELL multipliers (base cell = 64 H2O) and MPI ranks
declare -A MULT RANKS
MULT[64]="1 1 1";  MULT[128]="2 1 1"; MULT[256]="2 2 1"
MULT[512]="2 2 2"; MULT[1024]="4 2 2"
RANKS[64]=32; RANKS[128]=64; RANKS[256]=76; RANKS[512]=76; RANKS[1024]=76

SNAP_INT_STEPS=$(python3 -c "print(int(${SNAP_SPACING_PS}/${TIMESTEP_FS}*1000))")
TOTAL_STEPS=$(python3 -c "print(int((${EQUIL_PS}+${N_SNAPSHOTS}*${SNAP_SPACING_PS})/${TIMESTEP_FS}*1000))")

CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
export LD_LIBRARY_PATH="$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

mx=$(awk '{print $1}' <<<"${MULT[$size]}")
my=$(awk '{print $2}' <<<"${MULT[$size]}")
mz=$(awk '{print $3}' <<<"${MULT[$size]}")
mpi=${RANKS[$size]}
rundir="$RESULTS_SCRATCH/equil/N${size}"
mkdir -p "$rundir"
ln -sfn "$BENCH/potentials" "$rundir/NNP"

echo "=== Fig.S4 STAGE 1: NVT equilibration  N=${size} H2O ==="
echo "cells ${MULT[$size]} | $mpi ranks | $TOTAL_STEPS steps | snapshot every $SNAP_INT_STEPS steps"

# --- build the NVT input ----------------------------------------------------
export base_inp="$BENCH/H2O-64_RPBE-vdW_NNP.inp" target="$rundir/equil.inp" \
       size mx my mz TOTAL_STEPS SNAP_INT_STEPS
python3 - <<'PYEOF'
import os, re
txt = open(os.environ['base_inp']).read()
size = os.environ['size']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
total, snap = os.environ['TOTAL_STEPS'], os.environ['SNAP_INT_STEPS']

# Every edit is anchored to line-start (re.M) so that comment prose in the
# base input can never be matched -- only real, indented CP2K directives are.
txt = txt.replace('PROJECT H2O_RPBE-vdW_figS4', f'PROJECT H2O-{size}_equil')
txt = re.sub(r'^([ \t]*)STEPS[ \t]+\d+', rf'\g<1>STEPS {total}', txt, count=1, flags=re.M)
# Equilibration only needs the restart snapshots.  &STRESS / &TRAJECTORY are
# print keys: deleting the section is not a reliable disable (a low-print-
# level key reverts to its on-default at PRINT_LEVEL LOW), so set the section
# parameter to OFF explicitly.
txt = re.sub(r'^([ \t]*)&STRESS\b[ \t]*$', r'\g<1>&STRESS OFF', txt, count=1, flags=re.M)
txt = re.sub(r'^([ \t]*)&TRAJECTORY\b[ \t]*$', r'\g<1>&TRAJECTORY OFF', txt, count=1, flags=re.M)
# dump a restart-history snapshot every SNAP_INT_STEPS steps
txt = re.sub(r'(^[ \t]*&RESTART\b(?!_HISTORY).*?&EACH[ \t]*\n[ \t]*MD[ \t]+)\d+',
             rf'\g<1>{snap}', txt, count=1, flags=re.S | re.M)
txt = re.sub(r'(^[ \t]*&RESTART_HISTORY\b.*?&EACH[ \t]*\n[ \t]*MD[ \t]+)\d+',
             rf'\g<1>{snap}', txt, count=1, flags=re.S | re.M)
# multi-cell replication for sizes > 64
if size != '64':
    txt = re.sub(r'^([ \t]*&CELL[ \t]*\n)',
                 rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, count=1, flags=re.M)
    txt = re.sub(r'^([ \t]*&SUBSYS[ \t]*\n)',
                 rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n',
                 txt, count=1, flags=re.M)
open(os.environ['target'], 'w').write(txt)
PYEOF

# --- run --------------------------------------------------------------------
( cd "$rundir" && srun --ntasks=$mpi --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i equil.inp >equil.out 2>&1 )

if ! grep -q "PROGRAM ENDED" "$rundir/equil.out"; then
   echo "*** equilibration FAILED for N=${size} -- tail of equil.out:"
   tail -25 "$rundir/equil.out" | sed 's/^/   /'
   exit 1
fi

# --- collect the last N_SNAPSHOTS restart-history files --------------------
# RESTART_HISTORY writes <PROJECT>-1_<step>.restart ; sort -V orders them by
# the embedded step number regardless of digit count (8000 before 30000).
mapfile -t snaps < <(ls -1 "$rundir"/H2O-${size}_equil-1_*.restart 2>/dev/null | sort -V)
if (( ${#snaps[@]} < N_SNAPSHOTS )); then
   echo "*** only ${#snaps[@]} restart snapshots found (need $N_SNAPSHOTS)"
   exit 1
fi
k=1
for s in "${snaps[@]: -$N_SNAPSHOTS}"; do
   cp "$s" "$rundir/snapshot_${k}.restart"
   echo "snapshot_${k}.restart  <-  $(basename "$s")"
   k=$((k+1))
done
echo "=== STAGE 1 done for N=${size}.  Snapshots in $rundir/ ==="
