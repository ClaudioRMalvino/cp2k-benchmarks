#!/usr/bin/env bash
# Fig. S4 replication driver (Morawietz et al., PNAS 113 (2016) 8368).
# Login-node orchestrator: verifies inputs and submits three chained SLURM
# jobs (equil array -> prod array -> analysis).
set -euo pipefail

BENCH=/home/crm98/cp2k-benchmarks
SCRIPTS="$BENCH/scripts/CSD3_benchmark_scripts/figS4"
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
REBUILD=${REBUILD:-0}

EQUIL_PS=${EQUIL_PS:-20}
N_SNAPSHOTS=${N_SNAPSHOTS:-5}
SNAP_SPACING_PS=${SNAP_SPACING_PS:-10}
PROD_PS=${PROD_PS:-100}

echo "=== Fig. S4 replication - submission ==="
echo "equil=${EQUIL_PS}ps  snapshots=${N_SNAPSHOTS}x${SNAP_SPACING_PS}ps  production=${PROD_PS}ps/segment"

POT="$BENCH/potentials/RPBE-vdW-2016"
for fql in input.nn scaling.data weights.001.data weights.008.data; do
   [[ -f "$POT/$fql" ]] || { echo "ERROR: missing $POT/$fql"; exit 1; }
done
echo "  potential : $POT  OK"

build_branch() {
   local label=$1 repo=$2 gitref=$3 script=$4
   local cache="$BIN_ROOT/$label"
   if [[ "$REBUILD" -eq 0 && -x "$cache/cp2k.psmp" ]]; then
      echo "  binary    : $label  cached  OK"; return
   fi
   echo "  binary    : building $label ..."
   mkdir -p "$cache/lib"
   [[ -n "$gitref" ]] && { git -C "$repo" stash || true; git -C "$repo" checkout "$gitref"; }
   bash "$script"
   cp    "$repo/install/bin/cp2k.psmp"   "$cache/cp2k.psmp"
   cp -P "$repo/install/lib"/libcp2k.so* "$cache/lib/" 2>/dev/null || true
   [[ -n "$gitref" ]] && { git -C "$repo" stash pop || true; }
}
build_branch master                    /home/crm98/cp2k_master    "" \
              "$BENCH/cp2k_master/CSD3_build_scripts/cp2k_CSD3_master_build.sh"
build_branch feature-nnp-native-spline /home/crm98/cp2k_optimized "feature/nnp-native-spline" \
              "$BENCH/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_rebuild.sh"

# identical md5 means a build was skipped or LD_LIBRARY_PATH crossed branches
echo "  --- binary verification ---"
md5sum "$BIN_ROOT/master/cp2k.psmp" "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp"
if [[ "$(md5sum < "$BIN_ROOT/master/cp2k.psmp")" == \
      "$(md5sum < "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp")" ]]; then
   echo "ERROR: the two branch binaries are identical - aborting"; exit 1
fi
[[ -f "$BIN_ROOT/setup" ]] || cp /home/crm98/cp2k_master/tools/toolchain/install/setup "$BIN_ROOT/setup"

mkdir -p "$BENCH/logs"
EXPORTS="ALL,EQUIL_PS=$EQUIL_PS,N_SNAPSHOTS=$N_SNAPSHOTS,SNAP_SPACING_PS=$SNAP_SPACING_PS,PROD_PS=$PROD_PS"

echo
echo "=== submitting SLURM chain ==="
EQUIL_JID=$(sbatch --parsable --export="$EXPORTS" "$SCRIPTS/run_nnp_figS4_equil.sh")
echo "  STAGE 1  figS4_equil      : $EQUIL_JID   (array 0-4, 5 sizes)"

PROD_JID=$(sbatch --parsable --dependency=afterok:$EQUIL_JID --export="$EXPORTS" \
                  "$SCRIPTS/run_nnp_figS4_production.sh")
echo "  STAGE 2  figS4_prod       : $PROD_JID   (array 0-49, 2 branches x 5 sizes x 5 segments, afterok:$EQUIL_JID)"

ANALYSIS_JID=$(sbatch --parsable --dependency=afterok:$PROD_JID \
                      "$SCRIPTS/run_nnp_figS4_analysis.sh")
echo "  STAGE 3  figS4_analysis   : $ANALYSIS_JID   (afterok:$PROD_JID)"

echo
echo "monitor:   squeue -u $USER -j $EQUIL_JID,$PROD_JID,$ANALYSIS_JID"
echo "logs:      $BENCH/logs/figS4_{equil,prod,analysis}_*.out"
echo "results:   $BENCH/results/figS4/analysis/figS4_summary.csv  (after STAGE 3)"
echo "figures:   $BENCH/plots/figS4_{replication,performance,accuracy}.png"
echo "$EQUIL_JID $PROD_JID $ANALYSIS_JID" > /tmp/figS4_chain_jobids
