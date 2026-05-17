#!/usr/bin/env bash
# Cachegrind uses dynamic binary translation (not ptrace), so the CSD3
# ptrace mitigation that broke MAQAO does not affect it.

#SBATCH -J NNP_cachegrind
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_cachegrind_%j.out

set +e
mkdir -p /home/crm98/cp2k-benchmarks/logs

. /etc/profile.d/modules.sh
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

which valgrind
valgrind --version

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks
OUT_ROOT=/rds/user/$USER/hpc-work/cp2k-benchmarks/cachegrind
TIMESTAMP=$(date +%d-%m_%H-%M)
RUNDIR=$OUT_ROOT/run_$TIMESTAMP
mkdir -p "$RUNDIR"

STEPS=5
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

run_branch() {
   local label=$1
   local exe="$BIN_ROOT/$label/cp2k.psmp"
   local lib="$BIN_ROOT/$label/lib"
   local bdir="$RUNDIR/$label"
   mkdir -p "$bdir"

   echo
   echo "==> CACHEGRIND : $label"
   if [[ ! -x "$exe" ]]; then
      echo "!! missing binary $exe"; return 1
   fi
   md5sum "$exe"

   export LD_LIBRARY_PATH="$lib:${LD_LIBRARY_PATH:-}"
   local cg_out="$bdir/cachegrind.out"

   cd "$bdir"
   ln -sfn "$RUNDIR/NNP" NNP
   cp "$RUNDIR/run.inp" run.inp

   echo "    starting valgrind..."
   /usr/bin/time -v \
      srun --ntasks=1 --cpus-per-task=1 --hint=nomultithread \
        valgrind --tool=cachegrind \
                 --cache-sim=yes \
                 --branch-sim=no \
                 --cachegrind-out-file="$cg_out" \
                 "$exe" -i run.inp > "$bdir/cp2k.out" 2> "$bdir/valgrind.err"
   local rc=$?
   echo "    valgrind exit code: $rc"

   if [[ -f "$cg_out" ]]; then
      echo "    cachegrind output: $cg_out ($(stat -c %s "$cg_out") bytes)"
      cg_annotate "$cg_out" > "$bdir/cg_annotate.report.txt" 2>&1
      echo "    cg_annotate report: $bdir/cg_annotate.report.txt"
      head -60 "$bdir/cg_annotate.report.txt"
      grep -nE 'nnp_|cp2k.psmp' "$bdir/cg_annotate.report.txt" | head -30 > "$bdir/cg_nnp_rows.txt" || true
   else
      echo "    !! no cachegrind output produced"
      tail -30 "$bdir/valgrind.err"
   fi
   unset LD_LIBRARY_PATH
}

ORIG_LD=${LD_LIBRARY_PATH:-}
run_branch master
export LD_LIBRARY_PATH="$ORIG_LD"
run_branch feature-nnp-native-spline
export LD_LIBRARY_PATH="$ORIG_LD"

HOME_OUT=$BENCHMARK_ROOT/results/cachegrind/$TIMESTAMP
mkdir -p "$HOME_OUT"
for b in master feature-nnp-native-spline; do
   if [[ -d "$RUNDIR/$b" ]]; then
      mkdir -p "$HOME_OUT/$b"
      cp "$RUNDIR/$b/cg_annotate.report.txt" "$HOME_OUT/$b/" 2>/dev/null
      cp "$RUNDIR/$b/cg_nnp_rows.txt"        "$HOME_OUT/$b/" 2>/dev/null
      # Raw cachegrind.out kept on /rds only due to size.
   fi
done
echo
echo "Done. Reports in $HOME_OUT and raw cachegrind output on $RUNDIR"
