#!/usr/bin/env bash
#SBATCH -J nacl_mp2_split
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=00:45:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/nacl_mp2_split_%j.out

# GROUP_PARTITION (NNP/Fist ranks) sweep: the 75/1 cerberus split leaves the single
# Coulomb rank as the bottleneck (smoke 31514146: fist 0.134 vs nnp 0.057 s/step).
# Measures 4 splits on both boxes; rank split is performance-only, physics-neutral.
set -euo pipefail
SDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
source "$SDIR/env_csd3.sh"

ANCHOR_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor
MODEL_DIR="$ANCHOR_ROOT/models/MP2/final_model"
CUBE_DIR="$ANCHOR_ROOT/runs/MP2/cubic_1M"
SWEEP="$ANCHOR_ROOT/runs/MP2/split_sweep"
RANKS=$SLURM_NTASKS

declare -A STEPS=( [2]=200 [3]=100 )
SPLITS="1 4 8 16"          # Fist ranks; NNP gets the rest

echo "cell fist_ranks nnp_ranks steps t_per_step_s"
for n in 2 3; do
  abc=$(awk "BEGIN{printf \"%.4f\", 12.42*$n}")
  gx=$(( 40 * n ))
  for nf in $SPLITS; do
    nn=$(( RANKS - nf ))
    d="$SWEEP/cube${n}_f${nf}"
    mkdir -p "$d"
    sed -e "s|__PROJECT__|sweep_c${n}_f${nf}|" \
        -e "s|__STEPS__|${STEPS[$n]}|" \
        -e "s|__SNAPSTEPS__|1000000|" \
        -e "s|__MODEL__|$MODEL_DIR|" \
        -e "s|__COORD__|$CUBE_DIR/cube_n${n}.xyz|" \
        -e "s|__ABCX__|$abc|" -e "s|__ABCY__|$abc|" -e "s|__ABCZ__|$abc|" \
        -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
        -e "s|__GMAXX__|$gx|" -e "s|__GMAXY__|$gx|" -e "s|__GMAXZ__|$gx|" \
        -e "s|__NNPRANKS__|$nn|" \
        -e "s|__FISTRANKS__|$nf|" \
        "$SDIR/nacl_diffusion_template.inp" > "$d/sweep.inp"
    ( cd "$d" && srun --ntasks="$RANKS" --cpus-per-task=1 --hint=nomultithread \
        "$BIN" -i sweep.inp -o sweep.out ) || { echo "cube$n f$nf: RUN FAILED"; continue; }
    t=$(awk -v s="${STEPS[$n]}" '/qs_mol_dyn_low/ {printf "%.4f", $NF/s; exit}' "$d/sweep.out")
    echo "cube$n $nf $nn ${STEPS[$n]} ${t:-parse-error}"
  done
done
echo
echo "Pick the split with the lowest t/step per cell, then pass it to the"
echo "equil/production wrappers via FIST_RANKS_OVERRIDE (see run scripts)."
