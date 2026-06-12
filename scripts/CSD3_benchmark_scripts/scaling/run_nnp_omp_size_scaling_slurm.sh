#!/usr/bin/env bash
set -euo pipefail

# 2D OMP-x-size sweep for feature/nnp-native-spline-omp. For each
# OMP_NUM_THREADS in OMP_LIST, sweep the H2O system size and record
# time-per-MD-step (qs_mol_dyn_low row, startup-independent). MPI_RANKS is
# fixed (default 1) so only the per-rank OMP layer changes.

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

LABEL="${TARGET_LABEL:-feature-nnp-chebyshev}"
CP2K_EXE="$BIN_ROOT/$LABEL/cp2k.psmp"
INSTALL_LIB="$BIN_ROOT/$LABEL/lib"
PROJECT_ROOT="/home/crm98/cp2k_optimized"
OUTDIR_PARENT="cp2k_feature_chebyshev"

MPI_RANKS=${MPI_RANKS:-1}
STEPS=${STEPS:-100}   # 100 MD steps: uniform across every Report-2 measurement
N_REPS=${N_REPS:-5}
OMP_LIST="${OMP_LIST:-1 2 4 8}"
SIZE_LIST="${SIZE_LIST:-64 256 512 1024 2048}"
TIMESTAMP=$(date +%d-%m_%H-%M)

OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_omp_size_scaling_${LABEL}_${TIMESTAMP}"
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

CSV_FILE="${OUTDIR}/results_omp_size_scaling_${LABEL}_${TIMESTAMP}.csv"
RAW_CSV="${OUTDIR}/results_omp_size_scaling_${LABEL}_${TIMESTAMP}_raw.csv"
for f in "$CSV_FILE" "$RAW_CSV"; do
cat <<EOF >"$f"
# branch:    $LABEL
# exe:       $CP2K_EXE
# MPI ranks: $MPI_RANKS  (fixed)
# steps:     $STEPS
# reps:      $N_REPS timed (+1 discarded warm-up) per point
EOF
done
echo "# omp_threads,n_molecules,mpi_ranks,total_cores,rep,time_per_step_s,walltime_s" >>"$RAW_CSV"
echo "# omp_threads,n_molecules,mpi_ranks,total_cores,n_reps,time_per_step_mean_s,time_per_step_std_s,time_per_step_min_s,walltime_mean_s,walltime_std_s,walltime_min_s" >>"$CSV_FILE"

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

echo "OMP size scaling - $LABEL  ($MPI_RANKS MPI rank(s), OMP in: $OMP_LIST, $N_REPS reps)"
echo "-----------------------------------------------------------------------"

for omp in $OMP_LIST; do
   total_cores=$(( MPI_RANKS * omp ))
   if [[ $total_cores -gt $CORES_PER_NODE ]]; then
      echo "OMP=$omp skipped: $MPI_RANKS x $omp = $total_cores > $CORES_PER_NODE cores"
      continue
   fi

   for size in $SIZE_LIST; do
      mx=$(echo "${MULT[$size]}" | awk '{print $1}')
      my=$(echo "${MULT[$size]}" | awk '{print $2}')
      mz=$(echo "${MULT[$size]}" | awk '{print $3}')

      rundir="${OUTDIR}/omp_${omp}/N_${size}"
      mkdir -p "$rundir"

      export current_size=$size mx my mz \
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

      RAW_PREFIX="$omp,$size,$MPI_RANKS,$total_cores"
      read -r TM TS TMIN WM WS WMIN NOK \
           <<<"$(bench_point "$rundir" "$MPI_RANKS" "$omp")" || true

      if [[ "${NOK:-0}" == "0" ]]; then
         printf "OMP=%-3s N=%-6s : *** ALL %d REPS FAILED ***\n" "$omp" "$size" "$N_REPS"
         echo "$omp,$size,$MPI_RANKS,$total_cores,0,NA,NA,NA,NA,NA,NA" >>"$CSV_FILE"
      else
         printf "OMP=%-3s N=%-6s : reps=%-2s  t/step=%s+/-%s s  wall=%s+/-%s s\n" \
            "$omp" "$size" "$NOK" "$TM" "$TS" "$WM" "$WS"
         echo "$omp,$size,$MPI_RANKS,$total_cores,$NOK,$TM,$TS,$TMIN,$WM,$WS,$WMIN" >>"$CSV_FILE"
      fi
   done
done

echo "Wrote $CSV_FILE"
echo "      $RAW_CSV"

# Transitional shim: if reached from an older driver submitted before the
# decomposition sweep existed (DECOMP_VIA_DRIVER unset), run the sweep here
# so an already-queued job picks it up with no cancel/resubmit.
if [[ "${DECOMP_VIA_DRIVER:-0}" != 1 ]]; then
   echo
   echo "=== MPI:OMP DECOMPOSITION SWEEP (chained) ==="
   ./run_nnp_decomposition_sweep_slurm.sh || true
fi
