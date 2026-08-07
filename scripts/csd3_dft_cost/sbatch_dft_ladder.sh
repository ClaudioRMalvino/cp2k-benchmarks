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

RUNG="${1:?usage: sbatch_dft_ladder.sh <rung 1-4>}"
case "$RUNG" in
  1) NREP="1 1 1"; NAT=192  ;;
  2) NREP="2 1 1"; NAT=384  ;;
  3) NREP="2 2 1"; NAT=768  ;;
  4) NREP="2 2 2"; NAT=1536 ;;
  *) echo "rung must be 1-4" >&2; exit 2 ;;
esac

TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_dft_cost
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
export BIN_LABEL=master
source "$ADIR/env_csd3.sh"

DATA=/home/crm98/cp2k_master/data
ROOT=/rds/user/$USER/hpc-work/dft_cost_ladder
RUNDIR="$ROOT/n${NAT}"
COORD="$ROOT/unit64.xyz"
CSV=/home/crm98/cp2k-benchmarks/results/dft_cost/dft_ladder_timings.csv
STEPS=11          # step 1 runs from the ATOMIC guess and is discarded;
                  # steps 2-11 are the 10 steady-state AIMD steps we time
proj="water_revpbe_d3_n${NAT}"

echo "node: $(hostname)  rung: $RUNG  atoms: $NAT  replication: $NREP  ranks: $SLURM_NTASKS"
cat "$BIN_ROOT/$BIN_LABEL/PROVENANCE.txt" 2>/dev/null || true; echo

mkdir -p "$RUNDIR" "$(dirname "$CSV")"
rm -f "$RUNDIR/$proj"* "$RUNDIR/dft.out"

sed -e "s|__PROJECT__|$proj|" \
    -e "s|__NREP__|$NREP|g" \
    -e "s|__STEPS__|$STEPS|" \
    -e "s|__DATA__|$DATA|g" \
    -e "s|__COORD__|$COORD|" \
    "$TDIR/dft_ladder.inp.template" > "$RUNDIR/dft.inp"

cd "$RUNDIR"
t0=$SECONDS
srun --cpu-bind=cores "$BIN" -i dft.inp -o dft.out
wall=$((SECONDS - t0))

ener="$RUNDIR/$proj-1.ener"
[ -f "$ener" ] || { echo "FATAL: no .ener written"; exit 3; }

# steady-state s/step = mean UsedTime over the last 10 MD steps (step 1,
# which pays for the ATOMIC-guess SCF, is excluded by taking the tail)
read -r nsteps sps <<<"$(tail -n 10 "$ener" | awk '{s+=$NF; n++} END{printf "%d %.4f", n, (n? s/n : 0)}')"
scf=$(grep -c "SCF run converged" dft.out || true)
last=$(tail -1 "$ener" | awk '{print $1}')

echo "atoms=$NAT  timed_steps=$nsteps  s_per_step=$sps  last_step=$last  scf_converged=$scf  wall=${wall}s"
[ "$nsteps" -eq 10 ] || { echo "FATAL: expected 10 timed steps, got $nsteps"; exit 4; }
[ "$last" -eq "$STEPS" ] || { echo "FATAL: MD stopped at step $last of $STEPS"; exit 5; }

[ -f "$CSV" ] || echo "date,job_id,natoms,nwaters,replication,method,basis,cutoff_ry,binary,partition,nodes,ranks,omp,timed_steps,s_per_step,wall_s,scf_converged,notes" > "$CSV"
echo "$(date +%F),${SLURM_JOB_ID:-login},$NAT,$((NAT/3)),\"$NREP\",revPBE-D3,TZV2P-GTH,1200,$BIN_LABEL,${SLURM_JOB_PARTITION:-},${SLURM_NNODES:-1},$SLURM_NTASKS,${OMP_NUM_THREADS:-1},$nsteps,$sps,$wall,$scf,on-the-fly BOMD; steady-state mean over last 10 steps" >> "$CSV"
echo "logged to $CSV"
