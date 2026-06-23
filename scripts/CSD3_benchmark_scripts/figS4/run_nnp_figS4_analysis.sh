#!/usr/bin/env bash
# Fig. S4 STAGE 3: aggregate 50 NVE segments and plot.
#SBATCH -J figS4_analysis
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=8
# 1 h (was 30 min): aggregate_figS4.py parses trajectories SERIALLY; the 1.9 GB
# N512 (and the 3.8 GB N1024 once it lands) FFT-MSD passes need the headroom.
# ntasks=8 reserves ~27 GB on icelake so the largest trajectory loads safely.
#SBATCH --time=01:00:00
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
   echo "numpy not in module python - installing to --user"; pip install --user numpy; }

echo "=== aggregating production segments (chebyshev only; see equivalence table) ==="
python3 "$SCRIPTS/aggregate_figS4.py" \
   --prod-root "$SCRATCH/production" \
   --out-dir   "$SCRATCH/analysis" \
   --branches  feature-nnp-chebyshev

HOME_ANALYSIS="$BENCH/results/figS4/analysis"
mkdir -p "$HOME_ANALYSIS"
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH/analysis/" "$HOME_ANALYSIS/"
rsync -a --include='*/' --include='timing.csv' --include='viscosity.csv' \
      --include='diffusion.csv' --exclude='*' --prune-empty-dirs \
      "$SCRATCH/production/" "$BENCH/results/figS4/production/"

echo "=== plotting ==="
# The cluster module python lacks matplotlib/packaging; ~/.fortran_env has them.
# Subshell so activating the venv does not disturb the numpy/module python above.
(
   source "$HOME/.fortran_env/bin/activate"
   python3 "$BENCH/plots/plot_figS4.py" \
      --analysis-dir "$HOME_ANALYSIS" \
      --plot-dir     "$BENCH/plots"
)

echo
echo "=== Fig. S4 replication complete ==="
echo "  summary : $HOME_ANALYSIS/figS4_summary.csv"
echo "  figures : $BENCH/plots/figS4_replication.png"
echo "            $BENCH/plots/figS4_performance.png"
echo "            $BENCH/plots/figS4_accuracy.png"
column -t -s, "$HOME_ANALYSIS/figS4_summary.csv" 2>/dev/null || cat "$HOME_ANALYSIS/figS4_summary.csv"
