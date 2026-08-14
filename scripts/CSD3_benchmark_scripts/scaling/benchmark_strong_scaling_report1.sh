#!/usr/bin/env bash
# Report-1-methodology (thesis 4.1.1, fig 4.4) pure-MPI strong scaling: master,
# native-spline, chebyshev as a matched set for thesis_figures.py (fig1+fig2).
# N=1024 H2O, NUMA-aware cores (19=half socket, 38=socket, 76=node; 64 straddles
# NUMA, omitted), 100 steps @ 0.5 fs, 5 reps + 1 warm-up, cost = qs_mol_dyn_low/steps; native-spline-omp excluded (pure-MPI comparison).
#SBATCH -J NNP_strong_r1
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=08:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_strong_r1_%j.out

set -uo pipefail
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
SCALING=/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling
BRANCHES="master feature-nnp-native-spline feature-nnp-chebyshev"
CORE_LIST_R1="1 2 4 8 16 19 32 38 76"   # report-1 NUMA-aware set

echo "=== BINARY VERIFICATION (md5) ==="
for b in $BRANCHES; do
   exe="$BIN_ROOT/$b/cp2k.psmp"
   [[ -x "$exe" ]] || { echo "!! missing $exe"; exit 1; }
   md5sum "$exe"
done
cd "$SCALING"

echo ""
echo "############ REPORT-1 PURE-MPI STRONG SCALING ############"
echo "N=1024 H2O | steps=100 | reps=5(+1 warm-up) | cores: $CORE_LIST_R1"
for b in $BRANCHES; do
   echo ""; echo "---- $b ----"
   N_MOLECULES=1024 STEPS=100 N_REPS=5 CORE_LIST="$CORE_LIST_R1" \
      ./run_nnp_core_scaling_slurm.sh "$b"
done

# Mirror CSVs from scratch to the in-repo results tree (CSVs only, no run dirs).
SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_R=/home/crm98/cp2k-benchmarks/results
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs "$SCRATCH/" "$HOME_R/"
echo ""; echo "=== done; report-1 strong-scaling CSVs synced to $HOME_R ==="
