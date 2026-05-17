#!/usr/bin/env bash
set -euo pipefail

# MPI:OMP decomposition sweep: fix system size and total cores, sweep the
# MPI:OMP split across every factor pair, running
# feature/nnp-native-spline-omp at each. Pure-MPI references (master and
# native-spline at TOTAL_CORES x 1) run alongside. TOTAL_CORES defaults to
# 72 because 2^3*3^2 factors into many MPI:OMP pairs (76 would only give
# 1,2,4,19,38,76). OMP_PROC_BIND/PLACES are set so a hybrid split is not
# penalised by threads being scattered across sockets/NUMA domains.

CORES_PER_NODE=76
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks

# Strict mode relaxed across the sourcing: cp2k_CSD3_env.sh pulls in the
# toolchain 'setup' which references unbound CP_DFLAGS - trips `set -u`.
. /etc/profile.d/modules.sh
set +u
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
source "$BIN_ROOT/setup"
set -u
ORIG_LD="${LD_LIBRARY_PATH:-}"

N_MOLECULES=${N_MOLECULES:-512}
TOTAL_CORES=${TOTAL_CORES:-72}
STEPS=${STEPS:-100}
N_REPS=${N_REPS:-5}
TIMESTAMP=$(date +%d-%m_%H-%M)

[[ $TOTAL_CORES -le $CORES_PER_NODE ]] || { echo "TOTAL_CORES=$TOTAL_CORES > $CORES_PER_NODE"; exit 1; }

# Pin OMP threads to each rank's cores; without this the runtime can scatter
# them and make hybrid look bad for the wrong reason.
export OMP_PROC_BIND=close
export OMP_PLACES=cores

OUTDIR_PARENT="cp2k_feature_native_spline_omp"
OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_decomposition_sweep_${TIMESTAMP}"
mkdir -p "$OUTDIR"

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"

declare -A EXE LIB NNPDATA
EXE[master]="$BIN_ROOT/master/cp2k.psmp"
EXE[feature-nnp-native-spline]="$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
EXE[feature-nnp-native-spline-omp]="$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
LIB[master]="$BIN_ROOT/master/lib"
LIB[feature-nnp-native-spline]="$BIN_ROOT/feature-nnp-native-spline/lib"
LIB[feature-nnp-native-spline-omp]="$BIN_ROOT/feature-nnp-native-spline-omp/lib"
NNPDATA[master]="/home/crm98/cp2k_master/data/NNP"
NNPDATA[feature-nnp-native-spline]="/home/crm98/cp2k_optimized/data/NNP"
NNPDATA[feature-nnp-native-spline-omp]="/home/crm98/cp2k_optimized/data/NNP"

CSV_FILE="${OUTDIR}/results_decomposition_sweep_${TIMESTAMP}.csv"
RAW_CSV="${OUTDIR}/results_decomposition_sweep_${TIMESTAMP}_raw.csv"
for f in "$CSV_FILE" "$RAW_CSV"; do
cat <<EOF >"$f"
# decomposition sweep: fixed N=$N_MOLECULES H2O, fixed total cores=$TOTAL_CORES
# steps:    $STEPS
# reps:     $N_REPS timed (+1 discarded warm-up) per point
# binding:  OMP_PROC_BIND=close  OMP_PLACES=cores
# speedup_vs_pure_mpi is relative to feature-nnp-native-spline @ ${TOTAL_CORES}x1
EOF
done
echo "# branch,mpi_ranks,omp_threads,total_cores,rep,time_per_step_s,walltime_s" >>"$RAW_CSV"
echo "# branch,mpi_ranks,omp_threads,total_cores,n_reps,time_per_step_mean_s,time_per_step_std_s,time_per_step_min_s,walltime_mean_s,walltime_std_s,walltime_min_s,speedup_vs_pure_mpi" >>"$CSV_FILE"

stats() {
   printf '%s\n' "$@" | awk '
      /^[0-9.eE+-]+$/ { x[++n]=$1; s+=$1; if (n==1 || $1<mn) mn=$1 }
      END { if (n==0) { print "NA NA NA 0"; exit }
            m=s/n; v=0; for (i=1;i<=n;i++) v+=(x[i]-m)^2
            sd=(n>1)?sqrt(v/(n-1)):0
            printf "%.6f %.6f %.6f %d\n", m, sd, mn, n }'
}

bench_point() {
   local rundir=$1 mpi=$2 omp=$3
   local r out wt md tps tps_list="" wt_list=""
   export OMP_NUM_THREADS=$omp
   for r in $(seq 0 "$N_REPS"); do
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
      [[ $r -eq 0 ]] && continue
      tps_list+="$tps "
      wt_list+="$wt "
      echo "${RAW_PREFIX},${r},${tps},${wt}" >>"$RAW_CSV"
   done
   local tm tsd tmin tn wm wsd wmin wn
   read -r tm tsd tmin tn <<<"$(stats $tps_list)" || true
   read -r wm wsd wmin wn <<<"$(stats $wt_list)"  || true
   echo "$tm $tsd $tmin $wm $wsd $wmin $tn"
}

declare -A MULT
MULT[64]="1 1 1";   MULT[256]="2 2 1";  MULT[512]="2 2 2"
MULT[1024]="4 2 2"; MULT[2048]="4 4 2"; MULT[4096]="4 4 4"
[[ -n "${MULT[$N_MOLECULES]:-}" ]] || { echo "no MULTIPLE_UNIT_CELL for N=$N_MOLECULES"; exit 1; }
export current_size=$N_MOLECULES STEPS base_inp=$BASE_INP target_file="$OUTDIR/run.inp"
export mx=$(awk '{print $1}' <<<"${MULT[$N_MOLECULES]}")
export my=$(awk '{print $2}' <<<"${MULT[$N_MOLECULES]}")
export mz=$(awk '{print $3}' <<<"${MULT[$N_MOLECULES]}")
python3 - <<'PYEOF'
import os, re
txt = open(os.environ['base_inp']).read()
size = os.environ['current_size']; steps = os.environ['STEPS']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
txt = re.sub(r'STEPS\s+\d+', f'STEPS {steps}', txt, count=1)
# &TRAJECTORY / &RESTART_HISTORY are low-print-level keys: deleting the
# section reverts to on-default at PRINT_LEVEL LOW. Set the section param
# to OFF (also covers &RESTART, on by default at any level).
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
open(os.environ['target_file'], "w").write(txt)
PYEOF

REF_TPS=""
run_split() {
   local label=$1 mpi=$2 omp=$3
   local rundir="$OUTDIR/${label}__${mpi}mpi_x_${omp}omp"
   mkdir -p "$rundir"
   cp "$OUTDIR/run.inp" "$rundir/run.inp"
   ln -sfn "${NNPDATA[$label]}" "$rundir/NNP"
   CP2K_EXE="${EXE[$label]}"
   export LD_LIBRARY_PATH="${LIB[$label]}:$ORIG_LD"
   RAW_PREFIX="$label,$mpi,$omp,$TOTAL_CORES"

   local TM TS TMIN WM WS WMIN NOK
   read -r TM TS TMIN WM WS WMIN NOK \
        <<<"$(bench_point "$rundir" "$mpi" "$omp")" || true

   if [[ "${NOK:-0}" == "0" ]]; then
      printf "  %-28s %3smpi x %-3somp : *** ALL %d REPS FAILED ***\n" "$label" "$mpi" "$omp" "$N_REPS"
      echo "$label,$mpi,$omp,$TOTAL_CORES,0,NA,NA,NA,NA,NA,NA,NA" >>"$CSV_FILE"
      return
   fi
   [[ -z "$REF_TPS" ]] && REF_TPS=$TM
   local sp
   sp=$(awk -v r="$REF_TPS" -v t="$TM" 'BEGIN{ printf "%.3f", r/t }')
   printf "  %-28s %3smpi x %-3somp : reps=%-2s  t/step=%s+/-%s s  speedup_vs_pureMPI=%sx\n" \
      "$label" "$mpi" "$omp" "$NOK" "$TM" "$TS" "$sp"
   echo "$label,$mpi,$omp,$TOTAL_CORES,$NOK,$TM,$TS,$TMIN,$WM,$WS,$WMIN,$sp" >>"$CSV_FILE"
}

echo "=== MPI:OMP decomposition sweep - N=$N_MOLECULES H2O, $TOTAL_CORES cores, $N_REPS reps ==="
echo "--------------------------------------------------------------------------------"
# Pure-MPI reference runs FIRST so REF_TPS is set for the speedup column.
run_split feature-nnp-native-spline      "$TOTAL_CORES" 1
run_split master                         "$TOTAL_CORES" 1
for ((omp = 1; omp <= TOTAL_CORES; omp++)); do
   (( TOTAL_CORES % omp == 0 )) || continue
   mpi=$(( TOTAL_CORES / omp ))
   run_split feature-nnp-native-spline-omp "$mpi" "$omp"
done

echo "--------------------------------------------------------------------------------"
echo "Wrote $CSV_FILE"
echo "      $RAW_CSV"
