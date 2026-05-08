#!/bin/bash
set -e


module load rhel8/default-icl 2>/dev/null
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw 2>/dev/null
module load gcc/11 2>/dev/null
module load python/3.11.0-icl 2>/dev/null

export PATH=/usr/local/software/spack/spack-views/._rocky8-icelake-20220710/vdcqmdfi2lxsqy6qm3npptgktlrkt37k/intel-oneapi-mpi-2021.6.0/intel-2021.6.0/guxuvcpmykplbrr2e3af2yd7njqhau5e/mpi/2021.6.0/bin:$PATH

source ~/cp2k_master/tools/toolchain/install/setup

export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort

GCC11_LIB=$(dirname $(gcc -print-file-name=libstdc++.so))

module list
