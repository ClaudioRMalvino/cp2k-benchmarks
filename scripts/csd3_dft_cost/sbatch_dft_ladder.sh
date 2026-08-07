#!/usr/bin/env bash
#SBATCH -J dft_ladder
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/dft_ladder_%j.out

# On-the-fly DFT cost ladder: one rung per invocation.
#
#   usage: sbatch --time=<wall> sbatch_dft_ladder.sh <rung 1-4>
#
# Measures the steady-state Born-Oppenheimer AIMD cost per MD step for
# liquid water at revPBE-D3 (the O'Neill training-set level), at four system
# sizes built by replicating one 64-water box. The fitted cost-vs-N exponent
# extrapolates to the 5064-atom production cell, giving the "what this
# campaign would have cost with on-the-fly DFT" number that sits above the
# measured master / optimised-CPU / GPU NNP rungs in the cost ladder.
#
# Run on icelake with the same 1-node / 76-rank layout and the same pristine
# master binary (757bb76a80) as the master NNP reference, so DFT and NNP
# numbers are directly comparable.
set -euo pipefail

# The ladder is NaCl(aq) throughout - the campaign's own cell lineage, not a
# water proxy. Rungs 1-4 replicate cube_n1 (1 NaCl + 62 H2O, 12.42 A,
# 0.895 mol/kg) so molality is constant across the fitted rungs; rung 5 is
# the PRODUCTION cube_n3 configuration itself (30 NaCl + 1668 H2O, 37.26 A,
# 1.000 mol/kg) - the same coordinates the RPA campaign started from, so the
# headline number is measured on the system we actually study rather than
# extrapolated to it. Rung 5 gets fewer MD steps because each costs hours.
ROOT=/rds/user/$USER/hpc-work/dft_cost_ladder
PROD_CELL=/rds/user/$USER/hpc-work/nacl_mp2_anchor/runs/RPA/cubic_1M/cube_n3.xyz

RUNG="${1:?usage: sbatch_dft_ladder.sh <rung 1-5>}"
case "$RUNG" in
  1) NREP="1 1 1"; NAT=188  ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1.xyz"; MOL=0.895 ;;
  2) NREP="2 1 1"; NAT=376  ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1.xyz"; MOL=0.895 ;;
  3) NREP="2 2 1"; NAT=752  ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1.xyz"; MOL=0.895 ;;
  4) NREP="2 2 2"; NAT=1504 ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1.xyz"; MOL=0.895 ;;
  5) NREP="1 1 1"; NAT=5064 ; STEPS=4 ; ABC=37.2600; COORD="$PROD_CELL";        MOL=1.000 ;;
  *) echo "rung must be 1-5" >&2; exit 2 ;;
esac
TIMED=$((STEPS - 1))   # step 1 runs from the ATOMIC guess and is discarded

TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_dft_cost
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
export BIN_LABEL=master
source "$ADIR/env_csd3.sh"

DATA=/home/crm98/cp2k_master/data
RUNDIR="$ROOT/n${NAT}"
CSV=/home/crm98/cp2k-benchmarks/results/dft_cost/dft_ladder_timings.csv
proj="nacl_revpbe_d3_n${NAT}"

echo "node: $(hostname)  rung: $RUNG  atoms: $NAT  cell: $ABC A x $NREP  molality: $MOL  ranks: $SLURM_NTASKS"
echo "coord: $COORD"
cat "$BIN_ROOT/$BIN_LABEL/PROVENANCE.txt" 2>/dev/null || true; echo
[ -f "$COORD" ] || { echo "FATAL: coordinate file missing: $COORD"; exit 6; }

mkdir -p "$RUNDIR" "$(dirname "$CSV")"
rm -f "$RUNDIR/$proj"* "$RUNDIR/dft.out"

sed -e "s|__PROJECT__|$proj|" \
    -e "s|__NREP__|$NREP|g" \
    -e "s|__STEPS__|$STEPS|" \
    -e "s|__DATA__|$DATA|g" \
    -e "s|__ABC__|$ABC|g" \
    -e "s|__COORD__|$COORD|" \
    "$TDIR/dft_ladder.inp.template" > "$RUNDIR/dft.inp"

# the deck must describe the system we think it does
nat_in=$(grep -cE "^ *(O|H|Na|Cl) " "$COORD")
nrep_prod=$(awk -v n="$NREP" 'BEGIN{split(n,r," "); print r[1]*r[2]*r[3]}')
[ $((nat_in * nrep_prod)) -eq "$NAT" ] || {
  echo "FATAL: $COORD has $nat_in atoms x $nrep_prod replicas != $NAT"; exit 7; }

cd "$RUNDIR"
t0=$SECONDS
srun --cpu-bind=cores "$BIN" -i dft.inp -o dft.out
wall=$((SECONDS - t0))

ener="$RUNDIR/$proj-1.ener"
[ -f "$ener" ] || { echo "FATAL: no .ener written"; exit 3; }

# steady-state s/step = mean UsedTime over the last TIMED MD steps (step 0
# and step 1, which pays for the ATOMIC-guess SCF, fall outside the tail).
# NOTE: the .ener file is written incrementally, so if a long rung hits its
# wall limit the per-step timings still survive on disk and can be reduced by
# hand - a TIMEOUT loses the CSV row, not the measurement.
read -r nsteps sps sd <<<"$(tail -n $TIMED "$ener" | awk '{s+=$NF; q+=$NF*$NF; n++} END{m=(n? s/n:0); printf "%d %.4f %.4f", n, m, (n>1? sqrt((q/n-m*m)*n/(n-1)) : 0)}')"
scf=$(grep -c "SCF run converged" dft.out || true)
last=$(tail -1 "$ener" | awk '{print $1}')

echo "atoms=$NAT  timed_steps=$nsteps  s_per_step=$sps (sd $sd)  last_step=$last  scf_converged=$scf  wall=${wall}s"
[ "$nsteps" -eq "$TIMED" ] || { echo "FATAL: expected $TIMED timed steps, got $nsteps"; exit 4; }
[ "$last" -eq "$STEPS" ] || { echo "FATAL: MD stopped at step $last of $STEPS"; exit 5; }

[ -f "$CSV" ] || echo "date,job_id,natoms,system,molality,replication,cell_abc,method,basis,cutoff_ry,binary,partition,nodes,ranks,omp,timed_steps,s_per_step,sd_s,wall_s,scf_converged,notes" > "$CSV"
cell=$(awk -v n="$NREP" -v u="$ABC" 'BEGIN{split(n,r," "); printf "%.2fx%.2fx%.2f", u*r[1], u*r[2], u*r[3]}')
echo "$(date +%F),${SLURM_JOB_ID:-login},$NAT,NaCl(aq),$MOL,\"$NREP\",$cell,revPBE-D3,TZV2P-GTH,1200,$BIN_LABEL,${SLURM_JOB_PARTITION:-},${SLURM_NNODES:-1},$SLURM_NTASKS,${OMP_NUM_THREADS:-1},$nsteps,$sps,$sd,$wall,$scf,on-the-fly BOMD; steady-state mean over last $TIMED steps" >> "$CSV"
echo "logged to $CSV"
