#!/usr/bin/env bash
# Diagnostic smoke test for the MIXED FORCE_EVAL workaround on N=128.
# Same physics as smoke_N128_explicit.sh (Morawietz RPBE-vdW NNP, 384 atoms,
# explicit 24.84 x 12.42 x 12.42 cell) but wrapped in METHOD MIXED + a
# zero-charge Fist sibling, matching the dispatch pattern Schran et al.
# (J. Phys. Chem. Lett. 2024, 15, 6081) use for NaCl/water.
#
# Hypothesis: if T stays near 300 K and ConsQty is conserved, the
# instability lives in CP2K's standalone METHOD-NNP top-level dispatch,
# NOT in the NNP force calculation itself. MIXED routes the same NNP
# forces through a different integrator path and may bypass the bug.

#SBATCH -J figS4_smoke_N128_MIXED
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:15:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/smoke_N128_MIXED_%j.out

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"

set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

rundir=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/smoke_N128_MIXED_$(date +%H%M%S)
mkdir -p "$rundir"
cp /home/crm98/cp2k-benchmarks/H2O-128_explicit_RPBE-vdW_NNP_MIXED.inp "$rundir/run.inp"
ln -sfn /home/crm98/cp2k-benchmarks/potentials "$rundir/NNP"

echo "=== smoke N=128 MIXED (NNP + zero-charge Fist, Schran-style dispatch) ==="
echo "rundir: $rundir"
cd "$rundir"
srun --ntasks=32 --cpus-per-task=1 --hint=nomultithread "$CP2K_EXE" -i run.inp >cp2k.out 2>&1 || true

echo ""
echo "=== verdict: temperature column over the 200 steps ==="
ener=$(ls *.ener 2>/dev/null | head -1)
[ -z "$ener" ] && { echo "no .ener produced! tail of cp2k.out:"; tail -40 cp2k.out; exit 1; }
awk 'NR==2 || NR==10 || NR==50 || NR==100 || NR==200 {printf "  step=%-4s  t=%.3f ps  T=%.1f K  ConsQty=%.4f\n", $1, $2/1000, $4, $6}' "$ener"
maxT=$(awk 'NR>1 {if ($4>m) m=$4} END{print m}' "$ener")
echo ""

# Compare against the known-failing standalone case:
#   smoke_N128_explicit at step=198 saw T=1261 K, ConsQty drift from -87.98 to -84.54.
# Threshold of 1000 K is generous; if MIXED is the right fix we expect T to
# stay near 300 K and ConsQty to drift no more than O(10^-3).
if (( $(echo "$maxT > 1000" | bc -l) )); then
   echo "*** MIXED DID NOT HELP — max T = ${maxT} K — bug is in NNP force calc itself ***"
   echo "    (or in something the Fist sibling doesn't compensate for)"
else
   echo "MIXED IS STABLE — max T = ${maxT} K — bug is in standalone METHOD-NNP dispatch."
   echo "Workaround for production: wrap NNP in MIXED+zero-charge-Fist as in this input."
fi
echo ""
echo "rundir: $rundir"
echo "(for comparison: standalone smoke at logs/smoke_N128_explicit_29487913.out hit 1261 K by step 198)"
