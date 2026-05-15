#!/usr/bin/env bash
# ============================================================================
# Fig. S4 replication -- STAGE 3: analysis + plotting.
#
# Runs after the production array job.  Aggregates the 50 NVE segments into
# mean +/- SEM viscosity / diffusion / per-step-walltime per (branch, size),
# produces the 4-panel Fig. S4 figure plus performance and accuracy plots,
# and syncs all CSVs + PNGs back to ~/cp2k-benchmarks so they are reachable
# from the login node.
# ============================================================================
#SBATCH -J figS4_analysis
#SBATCH -A MPHIL-NIKIFORAKIS-CRM98-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=4
# Pure post-processing: aggregate 50 segments (FFT ACF + MSD) and plot.
# Real cost is ~5-10 min; 30 min cap leaves generous backfill-friendly margin.
#SBATCH --time=00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/figS4_analysis_%j.out

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

set -euo pipefail

BENCH=/home/crm98/cp2k-benchmarks
SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/figS4
SCRIPTS="$BENCH/scripts/CSD3_benchmark_scripts/figS4"

python3 -c "import numpy" 2>/dev/null || {
   echo "numpy not in module python -- installing to --user"; pip install --user numpy; }

# --- aggregate the 50 segments -> panel-ready CSVs (on scratch) -------------
echo "=== aggregating production segments ==="
python3 "$SCRIPTS/aggregate_figS4.py" \
   --prod-root "$SCRATCH/production" \
   --out-dir   "$SCRATCH/analysis"

# --- sync CSVs to home so plotting / inspection works from the login node --
HOME_ANALYSIS="$BENCH/results/figS4/analysis"
mkdir -p "$HOME_ANALYSIS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH/analysis/" "$HOME_ANALYSIS/"
# also bring back the per-segment timing/viscosity/diffusion CSVs for the record
rsync -a --include='*/' --include='timing.csv' --include='viscosity.csv' \
      --include='diffusion.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH/production/" "$BENCH/results/figS4/production/"

# --- produce the figures ---------------------------------------------------
echo "=== plotting ==="
python3 "$BENCH/plots/plot_figS4.py" \
   --analysis-dir "$HOME_ANALYSIS" \
   --plot-dir     "$BENCH/plots"

echo
echo "=== Fig. S4 replication complete ==="
echo "  summary : $HOME_ANALYSIS/figS4_summary.csv"
echo "  figures : $BENCH/plots/figS4_replication.png"
echo "            $BENCH/plots/figS4_performance.png"
echo "            $BENCH/plots/figS4_accuracy.png"
column -t -s, "$HOME_ANALYSIS/figS4_summary.csv" 2>/dev/null || cat "$HOME_ANALYSIS/figS4_summary.csv"
