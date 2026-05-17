#!/usr/bin/env bash
# Fig. S4 STAGE 1: NVT equilibration array job, one task per system size.
# Generates N_SNAPSHOTS shared NVE starting points used identically by both
# branches in STAGE 2.
#SBATCH -J figS4_equil
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=04:00:00
#SBATCH --array=1,3
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_equil_%A_%a.out

# Source toolchain env BEFORE strict mode: 'setup' references unbound
# CP_DFLAGS that would trip `set -u`.
. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCH=/home/crm98/cp2k-benchmarks
RESULTS_SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4
source "$BIN_ROOT/setup"

set -euo pipefail

EQUIL_PS=${EQUIL_PS:-30}
N_SNAPSHOTS=${N_SNAPSHOTS:-5}
SNAP_SPACING_PS=${SNAP_SPACING_PS:-15}
TIMESTEP_FS=0.5

SIZE_LIST=(64 128 256 512 1024)
size=${SIZE_LIST[$SLURM_ARRAY_TASK_ID]}

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

export base_inp="$BENCH/H2O-64_RPBE-vdW_NNP.inp" target="$rundir/equil.inp" \
       size mx my mz TOTAL_STEPS SNAP_INT_STEPS
python3 - <<'PYEOF'
import os, re
txt = open(os.environ['base_inp']).read()
size = os.environ['size']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
total, snap = os.environ['TOTAL_STEPS'], os.environ['SNAP_INT_STEPS']

# Anchor every edit at line-start so comment prose in the base input cannot
# be matched - only real, indented CP2K directives are.
txt = txt.replace('PROJECT H2O_RPBE-vdW_figS4', f'PROJECT H2O-{size}_equil')
txt = re.sub(r'^([ \t]*)STEPS[ \t]+\d+', rf'\g<1>STEPS {total}', txt, count=1, flags=re.M)
# &STRESS / &TRAJECTORY are low-print-level print keys: deleting the section
# is not a reliable disable (reverts to on-default at PRINT_LEVEL LOW), so
# set the section parameter to OFF explicitly.
txt = re.sub(r'^([ \t]*)&STRESS\b[ \t]*$', r'\g<1>&STRESS OFF', txt, count=1, flags=re.M)
txt = re.sub(r'^([ \t]*)&TRAJECTORY\b[ \t]*$', r'\g<1>&TRAJECTORY OFF', txt, count=1, flags=re.M)
txt = re.sub(r'(^[ \t]*&RESTART\b(?!_HISTORY).*?&EACH[ \t]*\n[ \t]*MD[ \t]+)\d+',
             rf'\g<1>{snap}', txt, count=1, flags=re.S | re.M)
txt = re.sub(r'(^[ \t]*&RESTART_HISTORY\b.*?&EACH[ \t]*\n[ \t]*MD[ \t]+)\d+',
             rf'\g<1>{snap}', txt, count=1, flags=re.S | re.M)
if size != '64':
    txt = re.sub(r'^([ \t]*&CELL[ \t]*\n)',
                 rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, count=1, flags=re.M)
    txt = re.sub(r'^([ \t]*&SUBSYS[ \t]*\n)',
                 rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n',
                 txt, count=1, flags=re.M)
open(os.environ['target'], 'w').write(txt)
PYEOF

( cd "$rundir" && srun --ntasks=$mpi --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i equil.inp >equil.out 2>&1 )

if ! grep -q "PROGRAM ENDED" "$rundir/equil.out"; then
   echo "*** equilibration FAILED for N=${size} - tail of equil.out:"
   tail -25 "$rundir/equil.out" | sed 's/^/   /'
   exit 1
fi

# sort -V orders by embedded step number regardless of digit count
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
