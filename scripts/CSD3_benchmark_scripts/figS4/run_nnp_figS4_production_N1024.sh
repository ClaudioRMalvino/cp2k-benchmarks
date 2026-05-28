#!/usr/bin/env bash
# Fig. S4 STAGE 2: NVE production array job (50 tasks = 2 branches x 5 sizes
# x 5 segments). Each task runs PROD_PS ps from a shared equilibration
# snapshot, writing stress every step and positions every 10 steps.
#SBATCH -J figS4_prod_N1024
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=12:00:00
#SBATCH --array=45-49%5
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_prod_N1024_%A_%a.out

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

PROD_PS=${PROD_PS:-100}
TIMESTEP_FS=0.5
PROD_STEPS=$(python3 -c "print(int(${PROD_PS}/${TIMESTEP_FS}*1000))")

idx=$SLURM_ARRAY_TASK_ID
BRANCH_LIST=(master feature-nnp-native-spline)
SIZE_LIST=(64 128 256 512 1024)
branch=${BRANCH_LIST[$(( idx / 25 ))]}
size=${SIZE_LIST[$(( (idx % 25) / 5 ))]}
seg=$(( (idx % 25) % 5 + 1 ))

case "$branch" in
  master)                    CP2K_EXE="$BIN_ROOT/master/cp2k.psmp";                   LIB="$BIN_ROOT/master/lib";                   LABEL=upstream-master ;;
  feature-nnp-native-spline) CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"; LIB="$BIN_ROOT/feature-nnp-native-spline/lib"; LABEL=feature-nnp-native-spline ;;
esac
export LD_LIBRARY_PATH="$LIB:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

# Same rank count for both branches at a given size keeps timing comparable.
# Smaller sizes use fewer ranks to avoid over-decomposition.
declare -A MULT RANKS
MULT[64]="1 1 1";  MULT[128]="2 1 1"; MULT[256]="2 2 1"
MULT[512]="2 2 2"; MULT[1024]="4 2 2"
RANKS[64]=16; RANKS[128]=32; RANKS[256]=64; RANKS[512]=76; RANKS[1024]=76
mx=$(awk '{print $1}' <<<"${MULT[$size]}")
my=$(awk '{print $2}' <<<"${MULT[$size]}")
mz=$(awk '{print $3}' <<<"${MULT[$size]}")
mpi=${RANKS[$size]}

SNAPSHOT="$RESULTS_SCRATCH/equil/N${size}/snapshot_${seg}.restart"
rundir="$RESULTS_SCRATCH/production/${branch}/N${size}/seg${seg}"
mkdir -p "$rundir"
ln -sfn "$BENCH/potentials" "$rundir/NNP"

echo "=== Fig.S4 STAGE 2: NVE production ==="
echo "branch=$LABEL  size=$size  segment=$seg  cells=${MULT[$size]}  ranks=$mpi  steps=$PROD_STEPS"
echo "snapshot: $SNAPSHOT"
if [[ ! -f "$SNAPSHOT" ]]; then
   echo "*** snapshot missing - did STAGE 1 (equil) complete for N=${size}?"
   exit 1
fi

export base_inp="$BENCH/H2O-64_RPBE-vdW_NNP.inp" target="$rundir/run.inp" \
       size mx my mz PROD_STEPS SNAPSHOT branch seg
python3 - <<'PYEOF'
import os, re
txt  = open(os.environ['base_inp']).read()
size = os.environ['size']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
steps = os.environ['PROD_STEPS']
snap  = os.environ['SNAPSHOT']
branch, seg = os.environ['branch'], os.environ['seg']

# Anchor every edit at line-start so comment prose in the base input cannot
# be matched - only real, indented CP2K directives are.
txt = txt.replace('PROJECT H2O_RPBE-vdW_figS4', f'PROJECT H2O-{size}_{branch}_seg{seg}')
txt = re.sub(r'^([ \t]*)ENSEMBLE[ \t]+NVT', r'\g<1>ENSEMBLE NVE', txt, count=1, flags=re.M)
txt = re.sub(r'^([ \t]*)STEPS[ \t]+\d+', rf'\g<1>STEPS {steps}', txt, count=1, flags=re.M)
txt = re.sub(r'^[ \t]*&THERMOSTAT\b.*?^[ \t]*&END THERMOSTAT[ \t]*\n', '',
             txt, flags=re.S | re.M)
# &RESTART_HISTORY is a low-print-level key: deleting the section reverts to
# its on-default at PRINT_LEVEL LOW. Set the section parameter to OFF to
# disable at any print level.
txt = re.sub(r'^([ \t]*)&RESTART_HISTORY\b[ \t]*$', r'\g<1>&RESTART_HISTORY OFF',
             txt, count=1, flags=re.M)
txt = re.sub(r'(^[ \t]*&RESTART\b(?!_HISTORY).*?&EACH[ \t]*\n[ \t]*MD[ \t]+)\d+',
             rf'\g<1>{steps}', txt, count=1, flags=re.S | re.M)
if size != '64':
    txt = re.sub(r'^([ \t]*&CELL[ \t]*\n)',
                 rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, count=1, flags=re.M)
    txt = re.sub(r'^([ \t]*&SUBSYS[ \t]*\n)',
                 rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n',
                 txt, count=1, flags=re.M)
# RESTART_DEFAULT .FALSE.: pull only coords + velocities from the snapshot,
# everything else (ensemble, cell, counters) comes from this fresh NVE input.
txt += (f"\n&EXT_RESTART\n"
        f"  RESTART_FILE_NAME {snap}\n"
        f"  RESTART_DEFAULT .FALSE.\n"
        f"  RESTART_POS .TRUE.\n"
        f"  RESTART_VEL .TRUE.\n"
        f"&END EXT_RESTART\n")
open(os.environ['target'], 'w').write(txt)
PYEOF

t0=$SECONDS
( cd "$rundir" && srun --ntasks=$mpi --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i run.inp >cp2k.out 2>&1 ) || true
wall=$(( SECONDS - t0 ))

CSV="$rundir/timing.csv"
if grep -q "PROGRAM ENDED" "$rundir/cp2k.out"; then
   total_wt=$(grep -E "^ CP2K +[0-9]" "$rundir/cp2k.out" | awk '{print $NF}' | tail -1)
   md_loop=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' "$rundir/cp2k.out" | tail -1)
   tps=$(awk -v t="${md_loop:-0}" -v s="$PROD_STEPS" 'BEGIN{ if(s>0) printf "%.6f", t/s; else print "NA" }')
   echo "branch,label,n_molecules,segment,mpi_ranks,prod_steps,time_per_step_s,md_loop_s,total_walltime_s" >"$CSV"
   echo "$branch,$LABEL,$size,$seg,$mpi,$PROD_STEPS,$tps,${md_loop:-NA},${total_wt:-NA}" >>"$CSV"
   echo "OK  N=$size seg=$seg $LABEL  time/step=${tps}s  total=${total_wt}s  (job wall ${wall}s)"
else
   echo "branch,label,n_molecules,segment,mpi_ranks,prod_steps,time_per_step_s,md_loop_s,total_walltime_s" >"$CSV"
   echo "$branch,$LABEL,$size,$seg,$mpi,$PROD_STEPS,FAILED,FAILED,FAILED" >>"$CSV"
   echo "*** FAILED  N=$size seg=$seg $LABEL - tail of cp2k.out:"
   tail -25 "$rundir/cp2k.out" | sed 's/^/   /'
   exit 1
fi
