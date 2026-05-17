#!/usr/bin/env bash
#SBATCH -J figS4_smoke_N64w
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/smoke_N64_wrapped_%j.out

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
source "$BIN_ROOT/setup"
set -euo pipefail

CP2K_EXE=$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp
export LD_LIBRARY_PATH=$BIN_ROOT/feature-nnp-native-spline/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

rundir=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4/smoke_N64_wrapped_$(date +%H%M%S)
mkdir -p "$rundir"
# Use the molecule-wrapped base input directly (no MULTIPLE_UNIT_CELL, just N=64)
sed 's/STEPS\s\+[0-9]\+/STEPS 200/; s/PROJECT H2O_RPBE-vdW_figS4/PROJECT H2O-64_wrapped_smoke/' \
    /home/crm98/cp2k-benchmarks/H2O-64_RPBE-vdW_NNP.inp > "$rundir/run.inp"
ln -sfn /home/crm98/cp2k-benchmarks/potentials "$rundir/NNP"

echo "=== smoke N=64 WRAPPED coords (no replication) ==="
cd "$rundir"
srun --ntasks=16 --cpus-per-task=1 --hint=nomultithread "$CP2K_EXE" -i run.inp >cp2k.out 2>&1 || true

ener=$(ls *.ener 2>/dev/null | head -1)
[ -z "$ener" ] && { echo "no .ener!"; tail -30 cp2k.out; exit 1; }
echo ""
echo "=== T column over 200 steps ==="
awk 'NR==2 || NR==10 || NR==50 || NR==100 || NR==200 {printf "  step=%-4s  T=%.1f K  ConsQty=%.4f\n", $1, $4, $6}' "$ener"
maxT=$(awk 'NR>1 {if ($4>m) m=$4} END{print m}' "$ener")
echo ""
if (( $(echo "$maxT > 500" | bc -l) )); then
   echo "*** WRAPPED N=64 ALSO BLOWS UP — wrap broke physics even at N=64 ***"
   echo "    -> need to restore original UN-wrapped coords"
else
   echo "WRAPPED N=64 IS STABLE — wrap is fine at N=64, problem emerges on doubling"
fi
echo "rundir: $rundir"
