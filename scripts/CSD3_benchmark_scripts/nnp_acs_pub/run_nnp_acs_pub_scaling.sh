#!/usr/bin/env bash
set -euo pipefail

# Per-branch sweep over (system_size, total_cores) for Figure 2a:
# time-per-MD-step vs. core count, one curve per system size.

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks

set +u
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

OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/paper_fig2/fig2a_scaling/${OUTDIR_PARENT}/${LABEL}_${TIMESTAMP}"
mkdir -p "$OUTDIR"

export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=$OMP_THREADS

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"
STEPS=200

CSV_FILE="${OUTDIR}/fig2a_${LABEL}_${TIMESTAMP}.csv"
cat <<EOF >"$CSV_FILE"
# Figure 2a wall-time per MD step
# branch: $LABEL
# exe:    $CP2K_EXE
# steps:  $STEPS
# omp:    $OMP_THREADS
n_molecules,mpi_ranks,omp_threads,total_cores,time_per_step_s,total_walltime_s
EOF

declare -A MULT
MULT[64]="1 1 1"
MULT[512]="2 2 2"
MULT[4096]="4 4 4"

# OMP branch uses OMP=2; pure-MPI branches use total_cores ranks.
CORE_LIST="2 4 8 16 32 64 72 144 288"
if [[ "$OMP_THREADS" -eq 1 ]]; then
   CORE_LIST="1 $CORE_LIST"
fi

echo "Branch: $LABEL  (OMP=$OMP_THREADS)"
printf "%-12s %-10s %-12s %-15s\n" "Molecules" "Cores" "Time/step(s)" "Walltime(s)"
echo "-------------------------------------------------------"

for size in 64 512 4096; do
   mx=$(echo "${MULT[$size]}" | awk '{print $1}')
   my=$(echo "${MULT[$size]}" | awk '{print $2}')
   mz=$(echo "${MULT[$size]}" | awk '{print $3}')

   for total_cores in $CORE_LIST; do
      mpi=$(( total_cores / OMP_THREADS ))
      [[ $mpi -lt 1 ]] && continue
      [[ $((mpi * OMP_THREADS)) -ne $total_cores ]] && continue

      rundir="${OUTDIR}/N${size}_c${total_cores}"
      mkdir -p "$rundir"

      export current_size=$size mx my mz target_file="${rundir}/run.inp" base_inp=$BASE_INP STEPS
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

txt = re.sub(r'&TRAJECTORY[\s\S]*?&END TRAJECTORY', '', txt)
txt = re.sub(r'&ENERGIES[\s\S]*?&END ENERGIES', '', txt)
txt = re.sub(r'&FORCES[\s\S]*?&END FORCES', '', txt)

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

      (cd "$rundir" && \
         srun --ntasks=$mpi --cpus-per-task=$OMP_THREADS --hint=nomultithread \
              "$CP2K_EXE" -i run.inp >cp2k.out 2>&1) || true

      if grep -q "PROGRAM ENDED" "${rundir}/cp2k.out" 2>/dev/null; then
         total_wt=$(grep -E "^ CP2K +[0-9]" "${rundir}/cp2k.out" | awk '{print $NF}' | tail -1)
         # Per-step time = qs_mol_dyn_low total / steps, which drops the SCF
         # init/shutdown overhead.
         md_loop=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' "${rundir}/cp2k.out" | tail -1)
         time_per_step=$(awk -v t="$md_loop" -v s="$STEPS" 'BEGIN{printf "%.6f", t/s}')
         printf "%-12s %-10s %-12s %-15s\n" "$size" "$total_cores" "$time_per_step" "$total_wt"
         echo "$size,$mpi,$OMP_THREADS,$total_cores,$time_per_step,$total_wt" >>"$CSV_FILE"
      else
         printf "%-12s %-10s %-12s %-15s\n" "$size" "$total_cores" "FAILED" "FAILED"
         echo "$size,$mpi,$OMP_THREADS,$total_cores,FAILED,FAILED" >>"$CSV_FILE"
      fi
   done
done

echo "Wrote $CSV_FILE"
