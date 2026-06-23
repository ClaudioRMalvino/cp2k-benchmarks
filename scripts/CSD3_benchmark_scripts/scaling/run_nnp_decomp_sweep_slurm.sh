#!/usr/bin/env bash
set -euo pipefail

# Hybrid MPI x OMP decomposition sweep at FIXED total cores and FIXED system
# size.  For each (mpi, omp) factorisation of $TOTAL_CORES it records the mean
# qs_mol_dyn_low time per MD step AND the CP2K-reported per-process peak memory,
# then derives aggregate node memory = peak_per_proc * mpi.
#
# Purpose: quantify (a) whether a branch can convert OMP threads into speed
# (chebyshev) or merely idles them (master / native-spline, whose NNP path is
# pure-MPI), and (b) the aggregate-memory cost of needing many MPI ranks
# (each replicates the model + tables) vs few ranks x many threads.
#
#   run_nnp_decomp_sweep_slurm.sh <branch>
# env: N_MOLECULES (1024), STEPS (100), N_REPS (3), TOTAL_CORES (76),
#      DECOMP_LIST ("76x1 38x2 19x4 4x19 2x38 1x76")

CORES_PER_NODE=76
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks

. /etc/profile.d/modules.sh
set +u
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
source "$BIN_ROOT/setup"
set -u

TARGET_BRANCH=${1:-master}
TIMESTAMP=$(date +%d-%m_%H-%M-%S)${SLURM_JOB_ID:+_${SLURM_JOB_ID}}

case "$TARGET_BRANCH" in
  feature-nnp-native-spline)
      CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-native-spline/lib"
      LABEL="feature-nnp-native-spline"; PROJECT_ROOT="/home/crm98/cp2k_optimized"
      OUTDIR_PARENT="cp2k_feature_native_spline" ;;
  feature-nnp-native-spline-omp)
      CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-native-spline-omp/lib"
      LABEL="feature-nnp-native-spline-omp"; PROJECT_ROOT="/home/crm98/cp2k_optimized"
      OUTDIR_PARENT="cp2k_feature_native_spline_omp" ;;
  feature-nnp-chebyshev)
      CP2K_EXE="$BIN_ROOT/feature-nnp-chebyshev/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-chebyshev/lib"
      LABEL="feature-nnp-chebyshev"; PROJECT_ROOT="/home/crm98/cp2k_optimized"
      OUTDIR_PARENT="cp2k_feature_chebyshev" ;;
  master|*)
      CP2K_EXE="$BIN_ROOT/master/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/master/lib"
      LABEL="upstream-master"; PROJECT_ROOT="/home/crm98/cp2k_master"
      OUTDIR_PARENT="cp2k_master" ;;
esac

N_MOLECULES=${N_MOLECULES:-1024}
STEPS=${STEPS:-100}
N_REPS=${N_REPS:-3}
TOTAL_CORES=${TOTAL_CORES:-76}
DECOMP_LIST=${DECOMP_LIST:-"76x1 38x2 19x4 4x19 2x38 1x76"}

OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_decomp_sweep_${LABEL}_${TIMESTAMP}"
mkdir -p "$OUTDIR"
export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"
BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"

CSV_FILE="${OUTDIR}/results_decomp_sweep_${LABEL}_${TIMESTAMP}.csv"
RAW_CSV="${OUTDIR}/results_decomp_sweep_${LABEL}_${TIMESTAMP}_raw.csv"
cat <<EOF >"$CSV_FILE"
# branch:       $LABEL
# exe:          $CP2K_EXE
# N_molecules:  $N_MOLECULES
# total_cores:  $TOTAL_CORES
# steps:        $STEPS
# reps:         $N_REPS timed (+1 warm-up) per decomposition
EOF
cp "$CSV_FILE" "$RAW_CSV"
echo "# mpi_ranks,omp_threads,total_cores,n_reps,tps_mean_s,tps_std_s,tps_min_s,peakmem_per_proc_mib,aggregate_mem_mib" >>"$CSV_FILE"
echo "# mpi_ranks,omp_threads,total_cores,rep,tps_s,peakmem_per_proc_mib" >>"$RAW_CSV"

mean_std_min() {  # stdin: list of numbers -> "mean std min n"
   awk '/^[0-9.eE+-]+$/{x[++n]=$1; s+=$1; if(n==1||$1<m)m=$1}
        END{ if(n==0){print "NA NA NA 0"; exit}
             mu=s/n; v=0; for(i=1;i<=n;i++)v+=(x[i]-mu)^2;
             sd=(n>1)?sqrt(v/(n-1)):0; printf "%.6f %.6f %.6f %d\n",mu,sd,m,n }'
}

declare -A MULT
MULT[64]="1 1 1"; MULT[256]="2 2 1"; MULT[512]="2 2 2"
MULT[1024]="4 2 2"; MULT[2048]="4 4 2"; MULT[4096]="4 4 4"
read -r mx my mz <<<"${MULT[$N_MOLECULES]}"

echo "Decomposition sweep: $LABEL  N=$N_MOLECULES  total_cores=$TOTAL_CORES  reps=$N_REPS"
echo "------------------------------------------------------------------------"

# Preflight: catch a node-local SLURM/SPANK fault (e.g. lua.so failing to
# dlopen liblua-5.3.so) BEFORE burning the whole campaign on silent reps=0.
# Job 30882875 died this way on a bad node and recorded zeros everywhere.
preflight="${OUTDIR}/srun_preflight.log"
srun --ntasks=1 hostname >"$preflight" 2>&1 || true
if grep -qiE 'Plug-in initialization failed|spank:|Dlopen of plugin' "$preflight"; then
   echo "!! FATAL: srun cannot launch on this node -- SLURM SPANK plugin fault:" >&2
   sed 's/^/   /' "$preflight" >&2
   echo "!! Aborting before recording bogus zeros. Resubmit (likely a bad node)." >&2
   exit 2
fi

for decomp in $DECOMP_LIST; do
   mpi=${decomp%x*}; omp=${decomp#*x}
   [[ $((mpi * omp)) -ne $TOTAL_CORES ]] && { echo "skip $decomp (!= $TOTAL_CORES cores)"; continue; }

   rundir="${OUTDIR}/decomp_${mpi}x${omp}"; mkdir -p "$rundir"
   export current_size=$N_MOLECULES mx my mz target_file="${rundir}/run.inp" base_inp=$BASE_INP STEPS
python3 - <<'PYEOF'
import os, re
txt = open(os.environ['base_inp']).read()
size, mx, my, mz, steps = (os.environ[k] for k in ('current_size','mx','my','mz','STEPS'))
txt = re.sub(r'STEPS\s+\d+', f'STEPS {steps}', txt, count=1)
txt = re.sub(r'&TRAJECTORY\b[\s\S]*?&END TRAJECTORY',
             '&TRAJECTORY OFF\n    &END TRAJECTORY\n    &RESTART OFF\n    &END RESTART\n'
             '    &RESTART_HISTORY OFF\n    &END RESTART_HISTORY', txt, count=1)
txt = re.sub(r'&ENERGIES\b', '&ENERGIES OFF', txt, count=1)
txt = re.sub(r'&FORCES\b', '&FORCES OFF', txt, count=1)
if size != "64":
    txt = re.sub(r'(&CELL\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.I)
    if "&TOPOLOGY" in txt:
        txt = re.sub(r'(&TOPOLOGY\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.I)
    else:
        txt = re.sub(r'(&SUBSYS\s*\n)', rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n', txt, flags=re.I)
open(os.environ['target_file'],'w').write(txt)
PYEOF
   ln -sfn "$NNP_DATA" "${rundir}/NNP"

   export OMP_NUM_THREADS=$omp OMP_PROC_BIND=spread OMP_PLACES=cores
   tps_list=""; mem_list=""
   for r in $(seq 0 "$N_REPS"); do
      out="$rundir/cp2k_rep${r}.out"
      ( cd "$rundir" && srun --ntasks="$mpi" --cpus-per-task="$omp" --hint=nomultithread \
           "$CP2K_EXE" -i run.inp >"$out" 2>&1 ) || true
      grep -q "PROGRAM ENDED" "$out" 2>/dev/null || continue
      md=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' "$out" | tail -1)
      tps=$(awk -v t="${md:-0}" -v s="$STEPS" 'BEGIN{ if(s>0) printf "%.6f", t/s; else print "NA" }')
      mem=$(grep 'Estimated peak process memory' "$out" | awk '{print $NF}' | tail -1)
      [[ $r -eq 0 ]] && continue
      tps_list+="$tps "; mem_list+="${mem:-NA} "
      echo "$mpi,$omp,$TOTAL_CORES,$r,$tps,${mem:-NA}" >>"$RAW_CSV"
   done
   read -r tm ts tmin tn <<<"$(printf '%s\n' $tps_list | mean_std_min)"
   read -r pm  _  _  _   <<<"$(printf '%s\n' $mem_list | mean_std_min)"
   peak=$(awk -v x="$pm" 'BEGIN{printf "%.0f", x}')
   agg=$(awk -v p="$pm" -v n="$mpi" 'BEGIN{printf "%.0f", p*n}')
   printf "%2s MPI x %2s OMP : reps=%s  t/step=%s+/-%s s  peak/proc=%s MiB  aggregate=%s MiB\n" \
      "$mpi" "$omp" "$tn" "$tm" "$ts" "$peak" "$agg"
   echo "$mpi,$omp,$TOTAL_CORES,$tn,$tm,$ts,$tmin,$peak,$agg" >>"$CSV_FILE"
done
echo "Wrote $CSV_FILE"
