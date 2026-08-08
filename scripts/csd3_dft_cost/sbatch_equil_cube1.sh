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
# with the campaign's own C-NNP + SPME model, so the DFT cost ladder's rungs
# 1-4 start from a real liquid rather than the builder's substituted lattice.
#
# WHY THIS EXISTS: make_base_cell.py places ions on former water-oxygen sites
# with no solvation shell. That configuration is electronically strained -
# rung 1 (n=188) converged but with a +-24% step-time spread, and rung 2
# (n=376) blew through 21 outer-SCF iterations without converging. Timing a
# strained frame OVERSTATES the AIMD cost, which biases the cost comparison in
# the MLIP's favour. The headline rung 5 already avoids this by using the
# production frame; this gives rungs 1-4 the same footing.
#
# CAVEAT (checked, expected to be fine): the ACSF cutoff is 12 bohr = 6.350 A,
# and this box's half-side is 6.210 A, so a few neighbour shells reach past the
# minimum image and some atoms see two periodic images of the same neighbour.
# The symmetry functions are sums over all neighbours inside r_c, so this is
# well defined - it just means the small cell carries a little spurious self-
# correlation. That is irrelevant here: the config only has to be a plausible
# equilibrated liquid for the SCF to converge on, and the DFT rungs use the
# revPBE-D3 electronic structure, not the NNP. If CP2K rejects the cell
# outright this job fails in the first seconds, which is the cheap answer.
#
# 150000 steps = 75 ps at 0.5 fs. That is 5x the 30000-step snapshot spacing,
# which is exactly the 5 restart-history files run_equil.sh requires; rungs
# 1-4 then use snapshot_5 (the most equilibrated).
#
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
