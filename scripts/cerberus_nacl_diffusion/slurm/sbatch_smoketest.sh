#!/usr/bin/env bash
#SBATCH --job-name=nacl_smoketest
#SBATCH --partition=csc-mphil
#SBATCH --clusters=CSC
#SBATCH --account=crm98
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --output=/home/raid/crm98/cp2k-benchmarks/logs/nacl_smoketest_%j.out

# Run this FIRST, once, before any real SLURM jobs. It verifies on the
# allocated compute node that:
#   1. the scratch filesystem with the binary + models is mounted,
#   2. the cerberus1-compiled (-march=native Ice Lake) binary runs on this
#      node's CPU without SIGILL,
#   3. a 200-step MD of the real NaCl deck completes and prints its speed.
set -euo pipefail
mkdir -p /home/raid/crm98/cp2k-benchmarks/logs

echo "node: $(hostname)"
lscpu | grep -E 'Model name|Socket|Core\(s\)'
echo

CP2K_ROOT=/data/cerberus1/crm98/cp2k_optimized
BENCH=/data/cerberus1/crm98/nacl_diffusion/bench
if [ ! -x "$CP2K_ROOT/install/bin/cp2k.psmp" ]; then
  echo "FAIL: $CP2K_ROOT/install/bin/cp2k.psmp not reachable from this node"
  echo "(is /data/cerberus1 mounted here? 'mount | grep cerberus':)"
  mount | grep -i cerberus || true
  exit 1
fi

set +u
source "$CP2K_ROOT/tools/toolchain/install/setup"
set -u
export LD_LIBRARY_PATH="$CP2K_ROOT/install/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

cd "$BENCH"
rm -f smoke.out
/usr/bin/time -p mpirun -n "$SLURM_NTASKS" --bind-to core \
  "$CP2K_ROOT/install/bin/cp2k.psmp" -i bench.inp -o smoke.out
grep -q "PROGRAM ENDED" smoke.out && echo "SMOKE TEST PASSED"
grep -m1 'qs_mol_dyn_low' smoke.out || true
