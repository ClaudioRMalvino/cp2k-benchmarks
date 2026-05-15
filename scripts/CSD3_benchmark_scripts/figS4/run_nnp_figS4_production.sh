#!/usr/bin/env bash
# ============================================================================
# Fig. S4 replication -- STAGE 2: NVE production.
#
# SLURM ARRAY JOB: 50 tasks = 2 branches x 5 sizes x 5 segments.
#   array index -> (branch, size, segment)
#       branch  = idx / 25            0 = upstream master
#                                     1 = feature/nnp-native-spline
#       size    = (idx % 25) / 5      [64,128,256,512,1024]
#       segment = (idx % 25) % 5 + 1  [1..5]
#
# Each task runs ONE NVE trajectory (PROD_PS ps, 0.5 fs) started from the
# shared equilibration snapshot for that (size, segment) -- so master and
# native-spline see IDENTICAL initial conditions and any divergence in eta / D
# is purely numerical.  Stress tensor is printed every step (Green-Kubo
# viscosity), positions every 10 steps (MSD diffusion).  Per-step and total
# wall-time go to timing.csv -> the same runs give both the eta/D agreement
# check AND the master-vs-native-spline performance comparison.
# ============================================================================
#SBATCH -J figS4_prod
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
# Slowest task is master @ N=1024, 100 ps (~4-5 h estimated on one node);
# 8 h cap covers it with margin.  Tighten once the validation run gives a
# real CSD3 per-step number -- most tasks (native-spline, small N) are far
# shorter, so a tight cap is what makes the 50-task array backfill well.
#SBATCH --time=08:00:00
#SBATCH --array=0-49%12
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_prod_%A_%a.out

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

# --- tunables --------------------------------------------------------------
PROD_PS=${PROD_PS:-100}            # NVE production length per segment
TIMESTEP_FS=0.5
PROD_STEPS=$(python3 -c "print(int(${PROD_PS}/${TIMESTEP_FS}*1000))")

# --- decode array index -> (branch, size, segment) -------------------------
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

# size -> MULTIPLE_UNIT_CELL multipliers and MPI rank count.  The SAME rank
# count is used for both branches at a given size, so the master-vs-native-
# spline timing is apples-to-apples.  Smaller sizes use fewer ranks to avoid
# over-decomposition (<1 atom/rank would skew the per-step timing).
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
   echo "*** snapshot missing -- did STAGE 1 (equil) complete for N=${size}?"
   exit 1
fi

# --- build the NVE input ---------------------------------------------------
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

# Every edit is anchored to line-start (re.M) so that comment prose in the
# base input can never be matched -- only real, indented CP2K directives are.
txt = txt.replace('PROJECT H2O_RPBE-vdW_figS4', f'PROJECT H2O-{size}_{branch}_seg{seg}')
txt = re.sub(r'^([ \t]*)ENSEMBLE[ \t]+NVT', r'\g<1>ENSEMBLE NVE', txt, count=1, flags=re.M)
txt = re.sub(r'^([ \t]*)STEPS[ \t]+\d+', rf'\g<1>STEPS {steps}', txt, count=1, flags=re.M)
# NVE: remove the thermostat block entirely
txt = re.sub(r'^[ \t]*&THERMOSTAT\b.*?^[ \t]*&END THERMOSTAT[ \t]*\n', '',
             txt, flags=re.S | re.M)
# Production doesn't need restart-history snapshots.  Deleting the section
# does NOT disable it: &RESTART_HISTORY is a low-print-level key, so at
# PRINT_LEVEL LOW it falls back to its default (a restart every 500 steps).
# Set the section parameter to OFF -- a true disable at any print level.
txt = re.sub(r'^([ \t]*)&RESTART_HISTORY\b[ \t]*$', r'\g<1>&RESTART_HISTORY OFF',
             txt, count=1, flags=re.M)
txt = re.sub(r'(^[ \t]*&RESTART\b(?!_HISTORY).*?&EACH[ \t]*\n[ \t]*MD[ \t]+)\d+',
             rf'\g<1>{steps}', txt, count=1, flags=re.S | re.M)
# multi-cell replication for sizes > 64 (must match the equilibration cell)
if size != '64':
    txt = re.sub(r'^([ \t]*&CELL[ \t]*\n)',
                 rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, count=1, flags=re.M)
    txt = re.sub(r'^([ \t]*&SUBSYS[ \t]*\n)',
                 rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n',
                 txt, count=1, flags=re.M)
# pull ONLY coordinates + velocities from the equilibration snapshot;
# RESTART_DEFAULT .FALSE. keeps everything else (ensemble, cell, counters,
# thermostat state) coming from this fresh NVE input.
txt += (f"\n&EXT_RESTART\n"
        f"  RESTART_FILE_NAME {snap}\n"
        f"  RESTART_DEFAULT .FALSE.\n"
        f"  RESTART_POS .TRUE.\n"
        f"  RESTART_VEL .TRUE.\n"
        f"&END EXT_RESTART\n")
open(os.environ['target'], 'w').write(txt)
PYEOF

# --- run --------------------------------------------------------------------
t0=$SECONDS
( cd "$rundir" && srun --ntasks=$mpi --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i run.inp >cp2k.out 2>&1 ) || true
wall=$(( SECONDS - t0 ))

CSV="$rundir/timing.csv"
if grep -q "PROGRAM ENDED" "$rundir/cp2k.out"; then
   # total wall-time from CP2K's timing table; per-step from the MD loop timer
   total_wt=$(grep -E "^ CP2K +[0-9]" "$rundir/cp2k.out" | awk '{print $NF}' | tail -1)
   md_loop=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' "$rundir/cp2k.out" | tail -1)
   tps=$(awk -v t="${md_loop:-0}" -v s="$PROD_STEPS" 'BEGIN{ if(s>0) printf "%.6f", t/s; else print "NA" }')
   echo "branch,label,n_molecules,segment,mpi_ranks,prod_steps,time_per_step_s,md_loop_s,total_walltime_s" >"$CSV"
   echo "$branch,$LABEL,$size,$seg,$mpi,$PROD_STEPS,$tps,${md_loop:-NA},${total_wt:-NA}" >>"$CSV"
   echo "OK  N=$size seg=$seg $LABEL  time/step=${tps}s  total=${total_wt}s  (job wall ${wall}s)"
else
   echo "branch,label,n_molecules,segment,mpi_ranks,prod_steps,time_per_step_s,md_loop_s,total_walltime_s" >"$CSV"
   echo "$branch,$LABEL,$size,$seg,$mpi,$PROD_STEPS,FAILED,FAILED,FAILED" >>"$CSV"
   echo "*** FAILED  N=$size seg=$seg $LABEL -- tail of cp2k.out:"
   tail -25 "$rundir/cp2k.out" | sed 's/^/   /'
   exit 1
fi
