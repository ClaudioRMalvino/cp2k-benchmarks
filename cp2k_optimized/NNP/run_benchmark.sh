#!/usr/bin/env bash

TIMESTAMP=$(date +%d-%m_%H-%M)
MPI_RANKS=${1:-1}
export OMP_NUM_THREADS=1

mpirun -np $MPI_RANKS /home/raid/crm98/cp2k_binaries/phy-cerberus/cp2k_feature_verlet_cells.psmp -i /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/NNP/H2O-64_NNP_MD.inp > /home/raid/crm98/cp2k-benchmarks/cp2k_optimized/NNP/cerberus1_benchmarks/cerberus-${TIMESTAMP}-proc${MPI_RANKS}.out 2>&1
