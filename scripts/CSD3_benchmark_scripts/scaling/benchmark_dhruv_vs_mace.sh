#!/usr/bin/env bash
# Head-to-head: Dhruv's linked-cell NNP (pr/nnp-cpu-cell-list) vs feature/nnp-mace,
# both built -O3 -xCORE-AVX2 -fp-model=precise + sequential MKL on the same
# toolchain, so the only difference is the NNP source. Report-1 axes (Ch.4): size
# scaling (76 MPI x 1 OMP, N=64..4096) + pure-MPI strong scaling (N=1024); 100 steps, 5 reps + 1 warm-up, qs_mol_dyn_low.
#SBATCH -J NNP_dhruv_vs_mace
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=10:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_dhruv_vs_mace_%j.out

set -uo pipefail
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
SCALING=/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling
BRANCHES="dhruv-cell-list feature-nnp-mace"
CORE_LIST_R1="1 2 4 8 16 19 32 38 76"   # report-1 NUMA-aware set

echo "=== BINARY VERIFICATION (md5) ==="
for b in $BRANCHES; do
   exe="$BIN_ROOT/$b/cp2k.psmp"
   [[ -x "$exe" ]] || { echo "!! missing $exe"; exit 1; }
   md5sum "$exe"
done
cd "$SCALING"

echo ""
echo "############ AXIS 1: SIZE SCALING (76 MPI x 1 OMP, N=64..4096) ############"
for b in $BRANCHES; do
   echo ""; echo "---- $b ----"
   ./run_nnp_size_scaling_slurm.sh "$b"
done

echo ""
echo "############ AXIS 2: PURE-MPI STRONG SCALING (N=1024) ############"
echo "steps=100 | reps=5(+1 warm-up) | cores: $CORE_LIST_R1"
for b in $BRANCHES; do
   echo ""; echo "---- $b ----"
   N_MOLECULES=1024 STEPS=100 N_REPS=5 CORE_LIST="$CORE_LIST_R1" \
      ./run_nnp_core_scaling_slurm.sh "$b"
done

# Mirror CSVs from scratch to the in-repo results tree (CSVs only, no run dirs).
SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_R=/home/crm98/cp2k-benchmarks/results
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs "$SCRATCH/" "$HOME_R/"
echo ""; echo "=== done; Dhruv-vs-mace CSVs synced to $HOME_R ==="
