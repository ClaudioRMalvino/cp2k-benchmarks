#!/usr/bin/env bash
# Login-node perf-stat cache-miss benchmark.  Same two-pass design as the
# SLURM version, but runs cp2k directly without srun.  Total wall time is
# a few minutes for STEPS=20 across both branches — well inside what CSD3
# permits for short login-node work, and costs nothing against the SL3
# budget.

set +e

. /etc/profile.d/modules.sh
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
set +e  # re-disable in case the env script turned it on

which perf
perf --version
echo "perf_event_paranoid = $(cat /proc/sys/kernel/perf_event_paranoid)"

# cp2k.psmp links MPI and OpenMP.  On the login node there is no SLURM PMI,
# so we use Intel MPI's Hydra launcher (mpirun -n 1) and force fabrics to
# shared memory to avoid touching any IB/UCX startup.  OMP_NUM_THREADS=1 keeps
# the cache counters clean (no thread-private working sets overlapping).
export I_MPI_FABRICS=shm
export OMP_NUM_THREADS=1

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks
OUT_ROOT=/rds/user/$USER/hpc-work/cp2k-benchmarks/perf_cache
TIMESTAMP=$(date +%d-%m_%H-%M)
RUNDIR=$OUT_ROOT/run_$TIMESTAMP
mkdir -p "$RUNDIR"

STEPS=20
echo "==> Preparing dataset: H2O-64, STEPS=$STEPS"
export base_inp="$BENCHMARK_ROOT/H2O-64_NNP_MD.inp" \
       target_file="$RUNDIR/run.inp" STEPS
python3 - <<'PYEOF'
import os, re
txt = open(os.environ['base_inp']).read()
steps = os.environ['STEPS']
txt = re.sub(r'STEPS\s+\d+', f'STEPS {steps}', txt, count=1)
txt = re.sub(r'&TRAJECTORY\b[\s\S]*?&END TRAJECTORY',
             '&TRAJECTORY OFF\n    &END TRAJECTORY\n'
             '    &RESTART OFF\n    &END RESTART\n'
             '    &RESTART_HISTORY OFF\n    &END RESTART_HISTORY', txt, count=1)
txt = re.sub(r'&ENERGIES\b', '&ENERGIES OFF', txt, count=1)
txt = re.sub(r'&FORCES\b',   '&FORCES OFF',   txt, count=1)
open(os.environ['target_file'], 'w').write(txt)
PYEOF
cp -rL /home/crm98/cp2k_optimized/data/NNP "$RUNDIR/NNP"

EVENTS_LOADS="cycles,instructions,\
mem_inst_retired.all_loads,\
mem_load_retired.l1_hit,mem_load_retired.l1_miss,\
mem_load_retired.l2_hit,mem_load_retired.l2_miss,\
mem_load_retired.l3_hit,mem_load_retired.l3_miss"

EVENTS_BACKING="cycles,instructions,\
l2_rqsts.all_demand_data_rd,l2_rqsts.all_demand_miss,\
l1d.replacement,\
dTLB-load-misses,dTLB-loads"

run_branch() {
   local label=$1
   local exe="$BIN_ROOT/$label/cp2k.psmp"
   local lib="$BIN_ROOT/$label/lib"
   local bdir="$RUNDIR/$label"
   mkdir -p "$bdir"

   echo
   echo "==> PERF CACHE : $label"
   if [[ ! -x "$exe" ]]; then
      echo "!! missing binary $exe"; return 1
   fi
   md5sum "$exe"

   export LD_LIBRARY_PATH="$lib:${LD_LIBRARY_PATH:-}"

   cd "$bdir"
   ln -sfn "$RUNDIR/NNP" NNP
   cp "$RUNDIR/run.inp" run.inp

   echo "    pass 1: mem_load_retired hierarchy"
   perf stat -e "$EVENTS_LOADS" \
             -o "$bdir/perf_loads.txt" \
             -- mpirun -n 1 "$exe" -i run.inp \
        > "$bdir/cp2k_loads.out" 2> "$bdir/cp2k_loads.err"
   echo "    pass 1 exit: $?"

   echo "    pass 2: l2_rqsts + l1d.replacement + dTLB"
   perf stat -e "$EVENTS_BACKING" \
             -o "$bdir/perf_backing.txt" \
             -- mpirun -n 1 "$exe" -i run.inp \
        > "$bdir/cp2k_backing.out" 2> "$bdir/cp2k_backing.err"
   echo "    pass 2 exit: $?"

   if [[ -f "$bdir/perf_loads.txt" ]]; then
      echo "    --- perf_loads.txt ---"
      cat "$bdir/perf_loads.txt"
   fi
   if [[ -f "$bdir/perf_backing.txt" ]]; then
      echo "    --- perf_backing.txt ---"
      cat "$bdir/perf_backing.txt"
   fi

   unset LD_LIBRARY_PATH
}

HOME_OUT=$BENCHMARK_ROOT/results/perf_cache/$TIMESTAMP

copy_back() {
   mkdir -p "$HOME_OUT"
   for b in master feature-nnp-native-spline; do
      if [[ -d "$RUNDIR/$b" ]]; then
         mkdir -p "$HOME_OUT/$b"
         cp "$RUNDIR/$b/perf_loads.txt"    "$HOME_OUT/$b/" 2>/dev/null
         cp "$RUNDIR/$b/perf_backing.txt"  "$HOME_OUT/$b/" 2>/dev/null
         cp "$RUNDIR/$b/cp2k_loads.out"    "$HOME_OUT/$b/" 2>/dev/null
         cp "$RUNDIR/$b/cp2k_loads.err"    "$HOME_OUT/$b/" 2>/dev/null
         cp "$RUNDIR/$b/cp2k_backing.out"  "$HOME_OUT/$b/" 2>/dev/null
         cp "$RUNDIR/$b/cp2k_backing.err"  "$HOME_OUT/$b/" 2>/dev/null
      fi
   done
   echo
   echo "Done. Reports in $HOME_OUT and raw cp2k output on $RUNDIR"
}
trap copy_back EXIT

ORIG_LD=${LD_LIBRARY_PATH:-}
run_branch master
export LD_LIBRARY_PATH="$ORIG_LD"
run_branch feature-nnp-native-spline
export LD_LIBRARY_PATH="$ORIG_LD"
