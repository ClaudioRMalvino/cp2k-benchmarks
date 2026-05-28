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

run_oneview master                        cp2k_master-config.json
run_oneview feature-nnp-native-spline      cp2k_optimized-config.json
run_oneview feature-nnp-native-spline-omp  cp2k_optimized_omp-config.json
export LD_LIBRARY_PATH="$ORIG_LD"

# Both binaries are pure-MPI (OMP=1) with identical NNP source symbol names,
# so MAQAO's function/loop name matching lines the kernels up directly.
echo
echo "=== MAQAO ONE View compare : master  vs  feature-nnp-native-spline"
if "$MAQAO" oneview --compare-reports \
         --inputs="$XP_ROOT/xp_master,$XP_ROOT/xp_feature-nnp-native-spline" \
         --include-detailed \
         --experiment-path="$XP_ROOT/xp_compare" \
         --replace; then
   echo "--- compare report written to $XP_ROOT/xp_compare/RESULTS"
else
   echo "!!! compare step failed -- the individual reports above are still intact"
fi

# Raw LProf/CQA data stays on /rds scratch; only HTML reports go to /home.
echo
echo "==> Syncing HTML reports to $HOME_REPORTS"
for d in xp_master xp_feature-nnp-native-spline xp_feature-nnp-native-spline-omp xp_compare; do
   if [[ -d "$XP_ROOT/$d/RESULTS" ]]; then
      mkdir -p "$HOME_REPORTS/$d"
      rsync -a "$XP_ROOT/$d/RESULTS/" "$HOME_REPORTS/$d/"
      echo "    $HOME_REPORTS/$d/index.html"
   fi
done

echo
echo "Done. View the reports in a browser, e.g.:"
echo "  firefox $HOME_REPORTS/xp_compare/index.html"
echo "  firefox $HOME_REPORTS/xp_master/index.html"
echo "Raw experiment directories remain on $XP_ROOT"
