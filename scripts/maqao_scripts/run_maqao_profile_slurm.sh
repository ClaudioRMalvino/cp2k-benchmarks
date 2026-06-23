#!/usr/bin/env bash
# CSD3 / Peta4-IceLake MAQAO ONE View profiling of the CP2K NNP code path.
# Override the MD length via: STEPS=1000 sbatch run_maqao_profile_slurm.sh

#SBATCH -J maqao_cp2k
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/maqao_cp2k_%j.out

mkdir -p /home/crm98/cp2k-benchmarks/logs/

# cp2k_CSD3_env.sh loads the Intel MKL / MPI modules + toolchain 'setup' that
# the cp2k.psmp binaries link against (missing these caused earlier runs to
# fail with 'libmkl_intel_thread.so.2: cannot open ...').
. /etc/profile.d/modules.sh
set +u
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
set +e   # cp2k_CSD3_env.sh enables `set -e`; MAQAO exit codes handled inline

# LProf attaches to the srun-spawned ranks with ptrace.  The zero-day
# hardening is being rolled back node-by-node: the June 10 probe worked on
# cpu-q-206 (ptrace_scope 0) but job 30352216 died on cpu-q-547 with
# "ptrace cannot attach: Operation not permitted".  Gate on a ptrace
# self-test so a blocked node costs seconds, not the whole walltime;
# on failure just resubmit to land on a different node.
echo "node: $(hostname)   ptrace_scope: $(cat /proc/sys/kernel/yama/ptrace_scope 2>/dev/null)   perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null)"
if ! srun --ntasks=1 python3 - <<'PYEOF'
import ctypes, os, sys, time
libc = ctypes.CDLL("libc.so.6", use_errno=True)
pid = os.fork()
if pid == 0:
    time.sleep(10)
    os._exit(0)
rc = libc.ptrace(16, pid, None, None)          # PTRACE_ATTACH
err = ctypes.get_errno()
if rc == 0:
    os.waitpid(pid, 0)
    libc.ptrace(17, pid, None, None)           # PTRACE_DETACH
    print("ptrace self-test: OK")
    sys.exit(0)
print(f"ptrace self-test: BLOCKED (errno {err}: {os.strerror(err)})")
sys.exit(1)
PYEOF
then
   echo "!!! ptrace is blocked on this compute node -- MAQAO LProf cannot work here."
   echo "!!! Resubmit the job; a different node may have the mitigation rolled back."
   exit 1
fi

MAQAO=/home/crm98/maqao.x86_64.2026.0.0-b/bin/maqao
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks
MAQAO_DIR=$BENCHMARK_ROOT/scripts/maqao_scripts
DATASET=$MAQAO_DIR/dataset
XP_ROOT=/rds/user/$USER/hpc-work/cp2k-benchmarks/maqao
HOME_REPORTS=$BENCHMARK_ROOT/results/maqao
mkdir -p "$XP_ROOT" "$HOME_REPORTS"

ORIG_LD=${LD_LIBRARY_PATH:-}
STEPS=${STEPS:-2000}

# Shared dataset across binaries. Per-step trajectory/energy/force I/O is
# switched OFF so the profile reflects compute (the ACSF kernels), not disk.
echo "==> Building profiling dataset ($STEPS MD steps) in $DATASET"
mkdir -p "$DATASET"
export base_inp="$BENCHMARK_ROOT/H2O-64_NNP_MD.inp" \
       target_file="$DATASET/run.inp" STEPS
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

rm -rf "$DATASET/NNP"
cp -rL /home/crm98/cp2k_optimized/data/NNP "$DATASET/NNP"
echo "    dataset ready: $(ls "$DATASET" | tr '\n' ' ')"

# LD_LIBRARY_PATH reset per-branch so each binary picks up its own lib/.
run_oneview () {
   local label=$1 cfg=$2
   echo
   echo "=== MAQAO ONE View : $label"
   export LD_LIBRARY_PATH="$BIN_ROOT/$label/lib:$ORIG_LD"
   if "$MAQAO" oneview --create-report=one \
            --config="$MAQAO_DIR/$cfg" \
            --experiment-path="$XP_ROOT/xp_$label" \
            --replace; then
      echo "--- $label : report written to $XP_ROOT/xp_$label/RESULTS"
   else
      echo "!!! $label : ONE View FAILED -- see output above"
   fi
}

# Report 2 scope: three profiles, mirroring Report 1's three-branch layout
# but for ONE branch under different parallelisation:
#   1. master            16x1  (pure-MPI baseline)
#   2. chebyshev         16x1  (pure MPI)
#   3. chebyshev (omp)    8x2  (OMP threads = 2, same 16 cores)
# The OMP chebyshev run is intentionally NEVER diffed against master:
# different parallelisation makes the per-function delta meaningless.  The
# two comparisons below are each like-for-like (see compare section).
run_oneview master                 cp2k_master-config.json
run_oneview feature-nnp-chebyshev  cp2k_chebyshev-config.json
XP_CHEBY_OMP="$XP_ROOT/xp_feature-nnp-chebyshev-omp"
echo
echo "=== MAQAO ONE View : feature-nnp-chebyshev (hybrid 8x2)"
export LD_LIBRARY_PATH="$BIN_ROOT/feature-nnp-chebyshev/lib:$ORIG_LD"
if "$MAQAO" oneview --create-report=one \
         --config="$MAQAO_DIR/cp2k_chebyshev_omp-config.json" \
         --experiment-path="$XP_CHEBY_OMP" \
         --replace; then
   echo "--- feature-nnp-chebyshev (8x2) : report written to $XP_CHEBY_OMP/RESULTS"
else
   echo "!!! feature-nnp-chebyshev (8x2) : ONE View FAILED -- see output above"
fi
export LD_LIBRARY_PATH="$ORIG_LD"

# Two like-for-like comparisons.  NEITHER diffs the OMP run against master.
#   1. BRANCH    : master 16x1 vs chebyshev 16x1 -- both pure MPI (OMP=1),
#                  identical NNP symbol names, so MAQAO lines the kernels up
#                  directly.  This is the "is chebyshev faster than master"
#                  story (the one that showed 1.42x slower pre-rework).
#   2. THREADING : chebyshev 16x1 (pure MPI) vs chebyshev 8x2 (OMP=2) -- ONE
#                  branch, same 16 cores, two parallelisation strategies.
#                  Shows what the OMP layer costs/saves vs pure MPI.
echo
echo "=== MAQAO compare 1/2 (BRANCH): master 16x1  vs  chebyshev 16x1 (both pure MPI) ==="
if "$MAQAO" oneview --compare-reports \
         --inputs="$XP_ROOT/xp_master,$XP_ROOT/xp_feature-nnp-chebyshev" \
         --include-detailed \
         --experiment-path="$XP_ROOT/xp_compare_branch" \
         --replace; then
   echo "--- branch compare written to $XP_ROOT/xp_compare_branch/RESULTS"
else
   echo "!!! branch compare failed -- the individual reports above are still intact"
fi

echo
echo "=== MAQAO compare 2/2 (THREADING): chebyshev pure-MPI 16x1  vs  OMP 8x2 ==="
if "$MAQAO" oneview --compare-reports \
         --inputs="$XP_ROOT/xp_feature-nnp-chebyshev,$XP_CHEBY_OMP" \
         --include-detailed \
         --experiment-path="$XP_ROOT/xp_compare_threading" \
         --replace; then
   echo "--- threading compare written to $XP_ROOT/xp_compare_threading/RESULTS"
else
   echo "!!! threading compare failed -- the individual reports above are still intact"
fi

# Raw LProf/CQA data stays on /rds scratch; only HTML reports go to /home.
echo
echo "==> Syncing HTML reports to $HOME_REPORTS"
for d in xp_master xp_feature-nnp-chebyshev xp_feature-nnp-chebyshev-omp \
         xp_compare_branch xp_compare_threading; do
   if [[ -d "$XP_ROOT/$d/RESULTS" ]]; then
      mkdir -p "$HOME_REPORTS/$d"
      rsync -a "$XP_ROOT/$d/RESULTS/" "$HOME_REPORTS/$d/"
      echo "    $HOME_REPORTS/$d/index.html"
   fi
done

echo
echo "Done. View the reports in a browser, e.g.:"
echo "  firefox $HOME_REPORTS/xp_compare_threading/index.html   # MPI vs OMP (chebyshev)"
echo "  firefox $HOME_REPORTS/xp_compare_branch/index.html      # chebyshev vs master (pure MPI)"
echo "  firefox $HOME_REPORTS/xp_feature-nnp-chebyshev-omp/index.html"
echo "Raw experiment directories remain on $XP_ROOT"
