#!/usr/bin/env bash
# Single-core 100-step NNP run of master vs feature/nnp-native-spline on the
# login node (no SLURM/perf/MAQAO); captures CP2K's end-of-run timing block.

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

   # Patch STEPS and silence trajectory/forces/energies/restart so the timing is NNP-dominated, not I/O.
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
