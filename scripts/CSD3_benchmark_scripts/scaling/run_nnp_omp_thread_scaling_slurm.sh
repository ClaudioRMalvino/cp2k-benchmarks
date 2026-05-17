#!/usr/bin/env bash
set -euo pipefail

# OMP thread scaling for feature/nnp-native-spline-omp at fixed system size
# and MPI rank count; reports time-per-MD-step, speedup and efficiency
# relative to OMP=1. time-per-step comes from the qs_mol_dyn_low row, so it
# is independent of SCF-init/shutdown overhead.

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

LABEL="feature-nnp-native-spline-omp"
CP2K_EXE="$BIN_ROOT/$LABEL/cp2k.psmp"
INSTALL_LIB="$BIN_ROOT/$LABEL/lib"
PROJECT_ROOT="/home/crm98/cp2k_optimized"
OUTDIR_PARENT="cp2k_feature_native_spline_omp"

MPI_RANKS=${MPI_RANKS:-1}
N_MOLECULES=${N_MOLECULES:-64}
STEPS=${STEPS:-100}
N_REPS=${N_REPS:-5}
OMP_LIST="${OMP_LIST:-1 2 4 8 16 32 ${CORES_PER_NODE}}"
TIMESTAMP=$(date +%d-%m_%H-%M)

OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_omp_thread_scaling_${LABEL}_${TIMESTAMP}"
mkdir -p "$OUTDIR"

export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"

declare -A MULT
MULT[64]="1 1 1"
MULT[256]="2 2 1"
MULT[512]="2 2 2"
MULT[1024]="4 2 2"
MULT[2048]="4 4 2"
MULT[4096]="4 4 4"
mx=$(echo "${MULT[$N_MOLECULES]}" | awk '{print $1}')
my=$(echo "${MULT[$N_MOLECULES]}" | awk '{print $2}')
mz=$(echo "${MULT[$N_MOLECULES]}" | awk '{print $3}')

CSV_FILE="${OUTDIR}/results_omp_thread_scaling_${LABEL}_${TIMESTAMP}.csv"
RAW_CSV="${OUTDIR}/results_omp_thread_scaling_${LABEL}_${TIMESTAMP}_raw.csv"
for f in "$CSV_FILE" "$RAW_CSV"; do
cat <<EOF >"$f"
# branch:      $LABEL
# exe:         $CP2K_EXE
# MPI ranks:   $MPI_RANKS  (fixed)
# N_molecules: $N_MOLECULES  (fixed)
# steps:       $STEPS
# reps:        $N_REPS timed (+1 discarded warm-up) per OMP count
EOF
done
echo "# omp_threads,mpi_ranks,total_cores,rep,time_per_step_s,walltime_s" >>"$RAW_CSV"
echo "# omp_threads,mpi_ranks,total_cores,n_reps,time_per_step_mean_s,time_per_step_std_s,time_per_step_min_s,walltime_mean_s,walltime_std_s,walltime_min_s,speedup,parallel_efficiency" >>"$CSV_FILE"

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

echo "OMP thread scaling - $LABEL  (N=$N_MOLECULES H2O, $MPI_RANKS MPI rank(s), $N_REPS reps)"
echo "-----------------------------------------------------------------------"

BASELINE_TPS=""
for omp in $OMP_LIST; do
   total_cores=$(( MPI_RANKS * omp ))
   [[ $total_cores -gt $CORES_PER_NODE ]] && continue

   rundir="${OUTDIR}/omp_${omp}"
   mkdir -p "$rundir"

   export current_size=$N_MOLECULES mx my mz \
          target_file="${rundir}/run.inp" base_inp=$BASE_INP STEPS
python3 - <<'PYEOF'
import os, re

base = os.environ['base_inp']
target = os.environ['target_file']
size = os.environ['current_size']
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
steps = os.environ['STEPS']

with open(base, "r") as f:
    txt = f.read()

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

with open(target, "w") as f:
    f.write(txt)
PYEOF

   ln -sfn "$NNP_DATA" "${rundir}/NNP"

   RAW_PREFIX="$omp,$MPI_RANKS,$total_cores"
   read -r TM TS TMIN WM WS WMIN NOK \
        <<<"$(bench_point "$rundir" "$MPI_RANKS" "$omp")" || true

   if [[ "${NOK:-0}" == "0" ]]; then
      printf "OMP=%-3s (%-3s cores) : *** ALL %d REPS FAILED ***\n" "$omp" "$total_cores" "$N_REPS"
      echo "$omp,$MPI_RANKS,$total_cores,0,NA,NA,NA,NA,NA,NA,NA,NA" >>"$CSV_FILE"
      continue
   fi

   if [[ -z "$BASELINE_TPS" ]]; then BASELINE_TPS=$TM; fi
   speedup=$(awk -v b="$BASELINE_TPS" -v t="$TM" 'BEGIN{ printf "%.3f", b/t }')
   eff=$(awk -v sp="$speedup" -v o="$omp" 'BEGIN{ printf "%.1f", 100*sp/o }')

   printf "OMP=%-3s (%-3s cores) : reps=%-2s  t/step=%s+/-%s s  speedup=%sx  eff=%s%%\n" \
      "$omp" "$total_cores" "$NOK" "$TM" "$TS" "$speedup" "$eff"
   echo "$omp,$MPI_RANKS,$total_cores,$NOK,$TM,$TS,$TMIN,$WM,$WS,$WMIN,$speedup,$eff" >>"$CSV_FILE"
done

echo "Wrote $CSV_FILE"
echo "      $RAW_CSV"
