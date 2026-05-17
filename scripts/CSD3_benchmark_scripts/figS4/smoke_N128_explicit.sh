#!/usr/bin/env bash
# 2-minute diagnostic smoke test for the N=128 setup, using an EXPLICIT
# 384-atom .inp (no MULTIPLE_UNIT_CELL). If T stays near 300 K, the bug
# is specifically MULTIPLE_UNIT_CELL + NNP interaction.

#SBATCH -J figS4_smoke_N128
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:15:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/smoke_N128_explicit_%j.out

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"

set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

rundir=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/smoke_N128_explicit_$(date +%H%M%S)
mkdir -p "$rundir"
cp /home/crm98/cp2k-benchmarks/H2O-128_explicit_RPBE-vdW_NNP.inp "$rundir/run.inp"
ln -sfn /home/crm98/cp2k-benchmarks/potentials "$rundir/NNP"

echo "=== smoke N=128 explicit (no MULTIPLE_UNIT_CELL) ==="
cd "$rundir"
srun --ntasks=32 --cpus-per-task=1 --hint=nomultithread "$CP2K_EXE" -i run.inp >cp2k.out 2>&1 || true

echo ""
echo "=== verdict: temperature column over the 200 steps ==="
ener=$(ls *.ener 2>/dev/null | head -1)
[ -z "$ener" ] && { echo "no .ener produced!"; tail -30 cp2k.out; exit 1; }
awk 'NR==2 || NR==10 || NR==50 || NR==100 || NR==200 {printf "  step=%-4s  t=%.3f ps  T=%.1f K  ConsQty=%.4f\n", $1, $2/1000, $4, $6}' "$ener"
maxT=$(awk 'NR>1 {if ($4>m) m=$4} END{print m}' "$ener")
echo ""
if (( $(echo "$maxT > 1000" | bc -l) )); then
   echo "*** STILL BLOWING UP — bug is NOT just MULTIPLE_UNIT_CELL ***"
else
   echo "STABLE — bug IS MULTIPLE_UNIT_CELL (explicit-coord workaround needed)"
fi
echo "rundir: $rundir"
