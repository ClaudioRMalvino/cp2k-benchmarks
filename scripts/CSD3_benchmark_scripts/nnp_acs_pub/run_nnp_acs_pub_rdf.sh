#!/usr/bin/env bash
set -euo pipefail

# Per-branch 64-H2O NVT trajectory at 300 K for RDF computation (Figure 2b).
# 10,000 steps × 0.5 fs = 5 ps simulation time, trajectory printed every 10
# steps → 1000 frames.  Single-node run on 72 cores.  RDFs are post-processed
# by compute_rdf.py from the SBATCH driver.

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
      OMP_THREADS=1
      MPI_RANKS=72
      ;;
  feature-nnp-native-spline-omp)
      CP2K_EXE="$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/feature-nnp-native-spline-omp/lib"
      LABEL="feature-nnp-native-spline-omp"
      PROJECT_ROOT="/home/crm98/cp2k_optimized"
      OMP_THREADS=2
      MPI_RANKS=36
      ;;
  master|*)
      CP2K_EXE="$BIN_ROOT/master/cp2k.psmp"
      INSTALL_LIB="$BIN_ROOT/master/lib"
      LABEL="upstream-master"
      PROJECT_ROOT="/home/crm98/cp2k_master"
      OMP_THREADS=1
      MPI_RANKS=72
      ;;
esac

# All three branches write into the same per-branch dir under fig2b_rdf/<label>/
# so the SBATCH driver can pick the trajectory up afterward without timestamp
# matching.
OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/paper_fig2/fig2b_rdf/${TARGET_BRANCH}"
mkdir -p "$OUTDIR"

export LD_LIBRARY_PATH="$INSTALL_LIB:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=$OMP_THREADS

BASE_INP="${BENCHMARK_ROOT}/H2O-64_NNP_MD.inp"
NNP_DATA="${PROJECT_ROOT}/data/NNP"
STEPS=10000
TRAJ_EVERY=10

# Generate input: 10k steps, trajectory every 10 steps, no force/energy spam
export base_inp=$BASE_INP target_file="${OUTDIR}/run.inp" STEPS TRAJ_EVERY
python3 - <<'PYEOF'
import os, re
base = os.environ['base_inp']
target = os.environ['target_file']
steps = os.environ['STEPS']
each = os.environ['TRAJ_EVERY']

with open(base, "r") as f:
    txt = f.read()

txt = re.sub(r'STEPS\s+\d+', f'STEPS {steps}', txt, count=1)
# trajectory every TRAJ_EVERY steps
txt = re.sub(r'(&TRAJECTORY\s*\n\s*&EACH\s*\n\s*MD\s+)\d+', rf'\g<1>{each}', txt)
# silence the per-step force/energy printing inside &NNP — pure RDF run
txt = re.sub(r'&ENERGIES[\s\S]*?&END ENERGIES', '', txt)
txt = re.sub(r'&FORCES[\s\S]*?&END FORCES', '', txt)

with open(target, "w") as f:
    f.write(txt)
PYEOF

ln -sfn "$NNP_DATA" "${OUTDIR}/NNP"

echo "Branch: $LABEL"
echo "MPI ranks: $MPI_RANKS, OMP threads: $OMP_THREADS, total cores: $((MPI_RANKS * OMP_THREADS))"
echo "Steps: $STEPS, trajectory every: $TRAJ_EVERY"

(cd "$OUTDIR" && \
   srun --ntasks=$MPI_RANKS --cpus-per-task=$OMP_THREADS --hint=nomultithread \
        "$CP2K_EXE" -i run.inp >cp2k.out 2>&1)

if grep -q "PROGRAM ENDED" "${OUTDIR}/cp2k.out" 2>/dev/null; then
   total_wt=$(grep -E "^ CP2K +[0-9]" "${OUTDIR}/cp2k.out" | awk '{print $NF}' | tail -1)
   echo "Done.  Total wall-time: ${total_wt}s.  Trajectory: ${OUTDIR}/H2O-64_NNP_MD-pos-1.xyz"
else
   echo "FAILED — see ${OUTDIR}/cp2k.out"
   exit 1
fi
