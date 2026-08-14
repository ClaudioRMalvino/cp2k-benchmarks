#!/usr/bin/env bash
#SBATCH -J equil_cube1
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=03:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/equil_cube1_%j.out

# Equilibrate the 188-atom cube1 cell (12.42 A, 1 NaCl + 62 H2O, 0.895 mol/kg)
# with the campaign's C-NNP + SPME model so DFT ladder rungs 1-4 start from a
# real liquid, not make_base_cell.py's strained substituted lattice (rung 1
# timed with +-24% spread, rung 2 failed 21 outer-SCF iterations; timing a
# strained frame OVERSTATES the AIMD cost, biasing the comparison pro-MLIP).
# Rung 5 already uses the production frame.
# Caveat (checked): ACSF cutoff 12 bohr = 6.350 A > half-side 6.210 A, so some
# neighbours double-count periodic images - well defined, and irrelevant since
# the DFT rungs use revPBE-D3, not the NNP.
# 150000 steps = 75 ps = 5x snapshot spacing -> the 5 restart-history files
# run_equil.sh needs; rungs 1-4 use snapshot_5.
#   sbatch sbatch_equil_cube1.sh
set -euo pipefail
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor

echo "node: $(hostname)  ranks: $SLURM_NTASKS"
MODEL=RPA CONC_DIR=cubic_1M CELLS=cube1 TOTAL_RANKS=$SLURM_NTASKS \
  FIST_RANKS_OVERRIDE=4 STEPS_OVERRIDE=150000 SKIP_DONE=0 \
  bash "$ADIR/run_equil.sh"

# hand the final snapshot to the DFT ladder as a bare @INCLUDE coordinate list
SNAP=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs/RPA/cubic_1M/equil/cube1/snapshot_5.restart
OUT=/rds/user/$USER/hpc-work/dft_cost_ladder/cube_n1_equil.xyz
python3 - "$SNAP" "$OUT" <<'PY'
import re, sys
src, dst = sys.argv[1], sys.argv[2]
text = open(src).read()
# &COORD ... &END COORD inside the restart holds "El x y z" in Angstrom
m = re.search(r'&COORD\n(.*?)\n\s*&END COORD', text, re.S)
if not m:
    sys.exit("no &COORD block in %s" % src)
rows = [ln.split() for ln in m.group(1).splitlines() if ln.split()]
rows = [r for r in rows if r[0] in ("O", "H", "Na", "Cl")]
with open(dst, "w") as fh:
    for el, x, y, z in (r[:4] for r in rows):
        fh.write("%-3s %18s %18s %18s\n" % (el, x, y, z))
print("wrote %d atoms -> %s" % (len(rows), dst))
PY

n=$(grep -cE "^ *(O|H|Na|Cl) " "$OUT")
[ "$n" -eq 188 ] || { echo "FATAL: expected 188 atoms in $OUT, got $n"; exit 8; }
echo "cube1 equilibrated; $OUT ready for DFT ladder rungs 1-4"
