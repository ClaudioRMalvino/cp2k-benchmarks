#!/usr/bin/env bash
# Cheap physics sanity check for feature/nnp-chebyshev vs master.
#
#   Test 1 (seconds):  single-point energy + forces, master vs branch.
#                      Differences should be ~1e-7 a.u. or below (the spline
#                      truncation error of MASTER). Anything >1e-5 = bug.
#   Test 2 (minutes):  2 ps NVE on 64 waters with EACH binary.
#                      Conserved-quantity drift is the sharpest probe of
#                      energy/force consistency (the clenshaw_d derivatives).
#                      A kinetic explosion like the Morawietz run shows up
#                      here within a few hundred steps.
#
# Submit with:  sbatch run_check.sh

#SBATCH -J nnp_physics_check
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/physics_check_%j.out
# Budget: 8 cores x <=30 min = <=4 CPU-hr worst case; typically well under 1.

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
set -euo pipefail

cd /home/crm98/cp2k-benchmarks/validation/physics_check

MASTER=/home/crm98/cp2k_master/install/bin/cp2k.psmp
BRANCH=/home/crm98/cp2k_optimized/install/bin/cp2k.psmp
[[ ! -x "$MASTER" ]] && { echo "missing $MASTER"; exit 1; }
[[ ! -x "$BRANCH" ]] && { echo "missing $BRANCH"; exit 1; }
export OMP_NUM_THREADS=1

run_side () {
    local side=$1 exe=$2
    mkdir -p out_$side
    # CP2K appends to existing .out files; start clean so the parser
    # never mixes runs.
    rm -f out_$side/sp_check.out out_$side/nve_drift.out out_$side/nve_drift-1.ener
    cp coord_H2O-64.inc out_$side/
    echo ">>> [$side] single-point energy+forces"
    ( cd out_$side && srun --ntasks=8 --cpus-per-task=1 --hint=nomultithread \
          "$exe" -i ../sp_check.inp -o sp_check.out )
    echo ">>> [$side] 2 ps NVE drift run"
    ( cd out_$side && srun --ntasks=8 --cpus-per-task=1 --hint=nomultithread \
          "$exe" -i ../nve_drift.inp -o nve_drift.out )
}

run_side master "$MASTER"
run_side branch "$BRANCH"

echo
echo "================ RESULTS ================"

echo "--- Test 1: single-point master vs branch ---"
for side in master branch; do
    grep "ENERGY| Total FORCE_EVAL" out_$side/sp_check.out | awk -v s=$side '{printf "  %-7s E_tot = %.12f Ha\n", s, $NF}'
done
extract_forces () {
    # CP2K prints "FORCES|   <atom>  fx  fy  fz  |f|" lines
    awk '/^ FORCES\|[ ]+[0-9]/{print $3, $4, $5}' "$1"
}
extract_forces out_master/sp_check.out > f_master.dat
extract_forces out_branch/sp_check.out > f_branch.dat
paste f_master.dat f_branch.dat | awk '
    {for(i=1;i<=3;i++){d=$(i)-$(i+3); if(d<0)d=-d; if(d>m)m=d}}
    END{printf "  max |dF| component = %.3e Ha/bohr  (expect ~1e-7, flag if >1e-5)\n", m}'

echo "--- Test 2: NVE conserved-quantity drift (192 atoms, 2 ps) ---"
for side in master branch; do
    awk -v s=$side 'NR==2{t0=$2;c0=$6} NR>1{t=$2;c=$6; if($4>tmax)tmax=$4}
        END{printf "  %-7s drift = %+.3e Ha total  (%+.3e Ha/ps/atom),  T_max = %.1f K\n",
            s, c-c0, (c-c0)/((t-t0)/1000.0)/192.0, tmax}' out_$side/nve_drift-1.ener
done

echo
echo "--- Extrapolation warnings (NNP leaving training range) ---"
for side in master branch; do
    n=$(grep -ci extrapolation out_$side/nve_drift.out || true)
    echo "  $side: $n"
done
echo
echo "Pass criteria: both drifts similar magnitude, |drift| < ~1e-6 Ha/ps/atom,"
echo "T stays in the 200-400 K band (no runaway heating)."
