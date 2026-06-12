#!/usr/bin/env bash
# Single-core 100-step NNP run of master and feature/nnp-native-spline.
# Goal: capture CP2K's hierarchical timing block (printed at the end of
# cp2k.out) for a side-by-side comparison.  Runs on the login node, no
# SLURM, no perf, no MAQAO — just the two binaries.

set +e
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
set +e

export OMP_NUM_THREADS=1

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks
TIMESTAMP=$(date +%d-%m_%H-%M)
OUT_DIR=$BENCHMARK_ROOT/results/timing/$TIMESTAMP
mkdir -p "$OUT_DIR"

STEPS=100

for label in master feature-nnp-native-spline; do
   workdir=$OUT_DIR/$label
   mkdir -p "$workdir"
   ln -sfn /home/crm98/cp2k_optimized/data/NNP "$workdir/NNP"

   # Copy the H2O-64 deck and patch STEPS + silence trajectory/forces/energies/restart
   # output so the timing report is dominated by the NNP cost, not file I/O.
   sed -e "s/STEPS  *[0-9]*/STEPS $STEPS/" \
       -e 's/&TRAJECTORY/\&TRAJECTORY OFF/' \
       -e 's/&FORCES/\&FORCES OFF/'         \
       -e 's/&ENERGIES/\&ENERGIES OFF/'     \
       -e 's/&RESTART/\&RESTART OFF/'       \
       "$BENCHMARK_ROOT/H2O-64_NNP_MD.inp" > "$workdir/run.inp"

   echo "==> timing run: $label"
   ( cd "$workdir"
     LD_LIBRARY_PATH="$BIN_ROOT/$label/lib:$LD_LIBRARY_PATH" \
        mpirun -n 1 "$BIN_ROOT/$label/cp2k.psmp" -i run.inp > cp2k.out 2>&1 )
   echo "    done: $workdir/cp2k.out"
done

echo
echo "Timing tables in $OUT_DIR/<branch>/cp2k.out (look for 'T I M I N G' block)"
