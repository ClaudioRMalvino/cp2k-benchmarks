#!/usr/bin/env bash
# GPU (ampere / A100) runtime environment for the transport-campaign run
# scripts, selected with ENV_SCRIPT=<this file> so run_equil.sh /
# run_production.sh work unchanged on the GPU track. Mirrors the module +
# library recipe of cp2k-mlp's run_water_diffusion.sh.
. /etc/profile.d/modules.sh 2>/dev/null || true
set +u
module purge
module load rhel8/default-amp || { echo "cannot load rhel8/default-amp" >&2; exit 1; }
set -u

HPC_WORK="${HPC_WORK:-/rds/user/$USER/hpc-work}"
BIN="${BIN:-$HPC_WORK/cp2k-mlp/install-collab-gpu-water-diffusion-gpu/bin/cp2k.psmp}"
[ -x "$BIN" ] || { echo "no GPU binary at $BIN" >&2; exit 1; }

_inst="$(cd "$(dirname "$BIN")/.." && pwd)"
for _d in lib lib64; do
  [ -d "$_inst/$_d" ] && export LD_LIBRARY_PATH="$_inst/$_d:${LD_LIBRARY_PATH:-}"
done
_ch="${CUDA_HOME:-}"
if [ -z "$_ch" ] && command -v nvcc >/dev/null 2>&1; then
  _ch="$(dirname "$(dirname "$(command -v nvcc)")")"
fi
[ -n "$_ch" ] && [ -d "$_ch/lib64" ] && export LD_LIBRARY_PATH="$_ch/lib64:${LD_LIBRARY_PATH:-}"

# The binary links MKL 2020.4 through its GNU (gf/gnu_thread) bindings plus
# scalapack/blacs_openmpi (paths from the build log, job 32427234); the amp
# module stack does not export these, so add them explicitly. libgomp comes
# from the spack gcc-9.4 runtime the build linked.
_mkl=/usr/local/Cluster-Apps/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mkl
for _d in "$_mkl/lib/intel64" "$_mkl/lib/intel64_lin" \
  /usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen3/gcc-11.2.0/gcc-9.4.0-72sgv5z3s2sluudvrz5eefoxh6oalq6v/lib64; do
  [ -d "$_d" ] && export LD_LIBRARY_PATH="$_d:${LD_LIBRARY_PATH:-}"
done

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR=/rds/user/$USER/hpc-work/tmp
mkdir -p "$TMPDIR"
