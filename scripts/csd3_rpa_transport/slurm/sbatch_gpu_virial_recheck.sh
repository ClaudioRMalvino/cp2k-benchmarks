#!/usr/bin/env bash
#SBATCH -J rpa_gpu_vircheck
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_vircheck_%j.out

# Pair-active GPU virial validation: the smoke's ON-vs-OFF stress check ran
# from the builder geometry, where ion-ion angular groups have EMPTY pair
# lists (the very condition that exposed the max_n_ang2 bug) - so those
# kernel instances were skipped in both paths. This job repeats the same
# 20-step USE_GPU ON vs OFF comparison from equil snapshot_1, where ion
# pairs exist and every kernel instance is live. Gates GPU production.
# Expected agreement ~1e-9 relative; hard-fail above 1e-6.
set -uo pipefail
CONC=${1:-cubic_1M}
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport
source "$TDIR/env_csd3_gpu.sh"

ANCHOR_ROOT="$HPC_WORK/nacl_mp2_anchor"
MODEL_DIR="$ANCHOR_ROOT/models/RPA/final_model"
CUBE_DIR="$ANCHOR_ROOT/runs_gpu/RPA/$CONC"
COORD="$CUBE_DIR/cube_n3.xyz"
SNAP="$CUBE_DIR/equil/cube3/snapshot_1.restart"
[ -f "$SNAP" ] || { echo "FAIL: missing $SNAP"; exit 1; }
BASE="$CUBE_DIR/virial_recheck"

fail=0
for mode in ON OFF; do
  d="$BASE/$(echo $mode | tr A-Z a-z)"; mkdir -p "$d"; cd "$d"; rm -f v.out NaCl_RPA_v*
  sed -e "s|__PROJECT__|NaCl_RPA_v|" \
      -e "s|__STEPS__|20|" \
      -e "s|__SNAPSTEPS__|100000|" \
      -e "s|__MODEL__|$MODEL_DIR|" \
      -e "s|__COORD__|$COORD|" \
      -e "s|__ABCX__|37.2600|" -e "s|__ABCY__|37.2600|" -e "s|__ABCZ__|37.2600|" \
      -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
      -e "s|__GMAXX__|120|" -e "s|__GMAXY__|120|" -e "s|__GMAXZ__|120|" \
      -e "s|__NNPRANKS__|1|" -e "s|__FISTRANKS__|3|" \
      -e "s|ENSEMBLE NVT|ENSEMBLE NVE|" \
      "$ADIR/nacl_diffusion_template.inp" \
    | sed -e '/&THERMOSTAT/,/&END THERMOSTAT/d' \
          -e 's|^\( *\)&RESTART_HISTORY$|\1\&RESTART_HISTORY OFF|' \
          -e "s|^  &NNP$|  \&NNP\n    USE_GPU $mode|" > v.inp
  cat >> v.inp <<EOF

&EXT_RESTART
  RESTART_FILE_NAME $SNAP
  RESTART_DEFAULT .FALSE.
  RESTART_POS .TRUE.
  RESTART_VEL .TRUE.
&END EXT_RESTART
EOF
  srun --ntasks="$SLURM_NTASKS" "$BIN" -i v.inp -o v.out
  grep -q "PROGRAM ENDED" v.out || { echo "FAIL: $mode run died"; fail=1; }
done
[ $fail -eq 0 ] || { echo "PAIR-ACTIVE VIRIAL CHECK FAILED (run error)"; exit 1; }

echo "--- pair-active GPU-vs-CPU virial (20 steps from $SNAP) ---"
paste "$BASE/on/NaCl_RPA_v-1.stress" "$BASE/off/NaCl_RPA_v-1.stress" \
  | awk 'NR>1 { n=(NF/2); for(i=3;i<=n;i++){ d=$i-$(i+n); if(d<0)d=-d; a=($i<0?-$i:$i);
                if(d>md)md=d; if(a>ma)ma=a } }
         END{ printf "max|dS| = %.3e bar   max|S| = %.3e bar   rel = %.3e\n", md, ma, md/ma;
              exit (md/ma > 1e-6) ? 1 : 0 }' || { echo "PAIR-ACTIVE VIRIAL CHECK FAILED (rel > 1e-6)"; exit 1; }
paste "$BASE/on/NaCl_RPA_v-1.ener" "$BASE/off/NaCl_RPA_v-1.ener" \
  | awk '!/^#/ { n=(NF/2); d=$5-$(5+n); if(d<0)d=-d; if(d>md)md=d } END{ printf "max|dE_pot| = %.3e Ha\n", md }'
echo "PAIR-ACTIVE VIRIAL CHECK PASSED"
