#!/usr/bin/env bash
set -euo pipefail

# CSD3 / Peta4-IceLake size-scaling sweep — reproduction of the cerberus
# run_nnp_size_scaling_slurm.sh.  Fixes parallelism at one full node and
# sweeps the H2O system size 64 -> 4096 via MULTIPLE_UNIT_CELL.
#
# Each size point is measured N_REPS+1 times: the first run is a discarded
# warm-up (cold caches / first MPI init), the remaining N_REPS are timed.
# Two CSVs are written into the run directory:
#   results_size_scaling_<label>_<ts>_raw.csv  every individual repeat
#   results_size_scaling_<label>_<ts>.csv      mean / std / min per size
# A single snapshot can land on a slow node or a noisy moment; the repeats
# give an honest error bar on the master-vs-native-spline comparison.
#
# Differences from cerberus: srun launcher, intel-oneapi-mkl module load,
# /rds scratch paths, per-step trajectory/energy/force printing stripped.

# Cores per Peta4-IceLake node — change here if your allocation differs.
CORES_PER_NODE=76

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks

# Self-sufficient module environment (idempotent if the driver already loaded
# them) so this script also works when invoked standalone.  Strict mode is
# relaxed across the sourcing: cp2k_CSD3_env.sh pulls in the toolchain 'setup',
# which references unbound vars (CP_DFLAGS) that would otherwise trip `set -u`.
. /etc/profile.d/modules.sh
set +u
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
source "$BIN_ROOT/setup"
set -u

TARGET_BRANCH=${1:-master}
TIMESTAMP=$(date +%d-%m_%H-%M)

case "$TARGET_BRANCH" in
  feature-nnp-native-spline)
      CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-native-spline/lib"
      LABEL="feature-nnp-native-spline"
      PROJECT_ROOT="/home/crm98/cp2k_optimized"
      OUTDIR_PARENT="cp2k_feature_native_spline"
      OMP_THREADS=1
      ;;
  feature-nnp-native-spline-omp)
      CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-native-spline-omp/lib"
      LABEL="feature-nnp-native-spline-omp"
      PROJECT_ROOT="/home/crm98/cp2k_optimized"
      OUTDIR_PARENT="cp2k_feature_native_spline_omp"
      OMP_THREADS=2
      ;;
  master|*)
      CP2K_EXE="$BIN_ROOT/master/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/master/lib"
      LABEL="upstream-master"
      PROJECT_ROOT="/home/crm98/cp2k_master"
      OUTDIR_PARENT="cp2k_master"
      OMP_THREADS=1
      ;;
esac

# One full node of parallelism.  OMP branch uses OMP=2, so MPI ranks are
# halved to keep total cores = CORES_PER_NODE.
MPI_RANKS=$(( CORES_PER_NODE / OMP_THREADS ))

# Number of TIMED repeats per size (a warm-up run is always done first and
# discarded).  Override with e.g. N_REPS=3 for a quicker, coarser sweep.
N_REPS=${N_REPS:-5}

OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_size_scaling_${LABEL}_${TIMESTAMP}"
mkdir -p "$OUTDIR"

export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=$OMP_THREADS

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"
STEPS=100

CSV_FILE="${OUTDIR}/results_size_scaling_${LABEL}_${TIMESTAMP}.csv"
RAW_CSV="${OUTDIR}/results_size_scaling_${LABEL}_${TIMESTAMP}_raw.csv"
for f in "$CSV_FILE" "$RAW_CSV"; do
cat <<EOF >"$f"
# branch:      $LABEL
# exe:         $CP2K_EXE
# MPI ranks:   $MPI_RANKS
# OMP threads: $OMP_THREADS
# total cores: $CORES_PER_NODE
# steps:       $STEPS
# reps:        $N_REPS timed (+1 discarded warm-up) per size
EOF
done
echo "# n_molecules,rep,time_per_step_s,walltime_s" >>"$RAW_CSV"
echo "# n_molecules,n_reps,time_per_step_mean_s,time_per_step_std_s,time_per_step_min_s,walltime_mean_s,walltime_std_s,walltime_min_s" >>"$CSV_FILE"

# --- helpers ----------------------------------------------------------------
# stats <num>... -> "mean sample_std min n"  (n=0 and NAs if no valid input)
stats() {
   printf '%s\n' "$@" | awk '
      /^[0-9.eE+-]+$/ { x[++n]=$1; s+=$1; if (n==1 || $1<mn) mn=$1 }
      END { if (n==0) { print "NA NA NA 0"; exit }
            m=s/n; v=0; for (i=1;i<=n;i++) v+=(x[i]-m)^2
            sd=(n>1)?sqrt(v/(n-1)):0
            printf "%.6f %.6f %.6f %d\n", m, sd, mn, n }'
}

# bench_point <rundir> <mpi> <omp> : run a warm-up + N_REPS timed reps of
# run.inp in <rundir>, append per-rep rows to $RAW_CSV (prefixed with the
# caller-set $RAW_PREFIX), and echo:
#   "tps_mean tps_std tps_min wt_mean wt_std wt_min n_ok"
bench_point() {
   local rundir=$1 mpi=$2 omp=$3
   local r out wt md tps tps_list="" wt_list=""
   export OMP_NUM_THREADS=$omp
   for r in $(seq 0 "$N_REPS"); do                      # r=0 is the warm-up
      out="$rundir/cp2k_rep${r}.out"
      ( cd "$rundir" && srun --ntasks="$mpi" --cpus-per-task="$omp" --hint=nomultithread \
           "$CP2K_EXE" -i run.inp >"$out" 2>&1 ) || true
      if ! grep -q "PROGRAM ENDED" "$out" 2>/dev/null; then
         [[ $r -gt 0 ]] && echo "${RAW_PREFIX},${r},FAILED,FAILED" >>"$RAW_CSV"
         continue
      fi
      wt=$(grep -E "^ CP2K +[0-9]" "$out" 2>/dev/null | awk '{print $NF}' | tail -1 || true)
      md=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' "$out" 2>/dev/null | tail -1 || true)
      tps=$(awk -v t="${md:-0}" -v s="$STEPS" 'BEGIN{ if (s>0) printf "%.6f", t/s; else print "NA" }')
      [[ $r -eq 0 ]] && continue                        # discard the warm-up
      tps_list+="$tps "
      wt_list+="$wt "
      echo "${RAW_PREFIX},${r},${tps},${wt}" >>"$RAW_CSV"
   done
   local tm tsd tmin tn wm wsd wmin wn
   read -r tm tsd tmin tn <<<"$(stats $tps_list)" || true
   read -r wm wsd wmin wn <<<"$(stats $wt_list)"  || true
   echo "$tm $tsd $tmin $wm $wsd $wmin $tn"
}

echo "Starting NNP size scaling on branch: $LABEL  ($MPI_RANKS MPI x $OMP_THREADS OMP, $N_REPS reps)"
echo "------------------------------------------------------------------------"

declare -A MULT
MULT[64]="1 1 1"
MULT[256]="2 2 1"
MULT[512]="2 2 2"
MULT[1024]="4 2 2"
MULT[2048]="4 4 2"
MULT[4096]="4 4 4"

for size in 64 256 512 1024 2048 4096; do
   rundir="${OUTDIR}/N_${size}"
   mkdir -p "$rundir"

   export current_size=$size
   export mx=$(echo "${MULT[$size]}" | awk '{print $1}')
   export my=$(echo "${MULT[$size]}" | awk '{print $2}')
   export mz=$(echo "${MULT[$size]}" | awk '{print $3}')
   export target_file="${rundir}/run.inp"
   export base_inp=$BASE_INP
   export STEPS

python3 - <<'PYEOF'
import os, re

base = os.environ['base_inp']
target = os.environ['target_file']
size = os.environ['current_size']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
steps = os.environ['STEPS']

with open(base, "r") as f:
    txt = f.read()

# Fix the step count, and strip per-step trajectory/energy/force printing —
# a size-scaling benchmark should measure compute, not file I/O.
txt = re.sub(r'STEPS\s+\d+', f'STEPS {steps}', txt, count=1)
# Performance benchmark: switch OFF all per-step file I/O so the timing
# reflects compute, not disk.  Deleting a print-key section does NOT disable
# it -- &TRAJECTORY / &RESTART_HISTORY are low-print-level keys that revert to
# their (on) defaults at PRINT_LEVEL LOW -- so set the section parameter to
# OFF explicitly (also covers &RESTART, which is on by default at any level).
txt = re.sub(r'&TRAJECTORY\b[\s\S]*?&END TRAJECTORY',
             '&TRAJECTORY OFF\n    &END TRAJECTORY\n'
             '    &RESTART OFF\n    &END RESTART\n'
             '    &RESTART_HISTORY OFF\n    &END RESTART_HISTORY', txt, count=1)
txt = re.sub(r'&ENERGIES\b', '&ENERGIES OFF', txt, count=1)
txt = re.sub(r'&FORCES\b', '&FORCES OFF', txt, count=1)

if size != "64":
    txt = re.sub(r'(&CELL\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.IGNORECASE)
    if "&TOPOLOGY" in txt:
        txt = re.sub(r'(&TOPOLOGY\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.IGNORECASE)
    else:
        txt = re.sub(r'(&SUBSYS\s*\n)', rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n', txt, flags=re.IGNORECASE)

with open(target, "w") as f:
    f.write(txt)
PYEOF

   ln -sfn "$NNP_DATA" "${rundir}/NNP"

   RAW_PREFIX="$size"
   read -r TM TS TMIN WM WS WMIN NOK \
        <<<"$(bench_point "$rundir" "$MPI_RANKS" "$OMP_THREADS")" || true

   if [[ "${NOK:-0}" == "0" ]]; then
      printf "%-10s tiles %-9s  *** ALL %d REPS FAILED ***\n" "$size" "${MULT[$size]}" "$N_REPS"
      echo "$size,0,NA,NA,NA,NA,NA,NA" >>"$CSV_FILE"
   else
      printf "%-10s tiles %-9s  reps=%-3s  t/step=%s+/-%s s  wall=%s+/-%s s\n" \
         "$size" "${MULT[$size]}" "$NOK" "$TM" "$TS" "$WM" "$WS"
      echo "$size,$NOK,$TM,$TS,$TMIN,$WM,$WS,$WMIN" >>"$CSV_FILE"
   fi
done

echo "Wrote $CSV_FILE"
echo "      $RAW_CSV"
