#!/usr/bin/env bash
#SBATCH -J dft_ladder
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/dft_ladder_%j.out

# On-the-fly revPBE-D3 AIMD cost per MD step, one rung per invocation.
# usage: sbatch --time=<wall> sbatch_dft_ladder.sh <rung 1-5>
#
# Rung 4 needs -p icelake-himem (job 33197088: OOM at ~283 GB vs 250 GB on
# plain icelake; himem is the same 76-core silicon, only more DRAM, so the
# 76x1 layout is preserved). Rung 3 is re-run on himem as a timing control;
# the CSV records the partition per row.
# Layout is fixed at 1 node x 76 MPI ranks - the layout the NNP reference
# rungs were measured on, on the same pristine master binary (757bb76a80).
# Rung 5 (5064 atoms) stalls in transfer_rs2pw_distributed before the first
# SCF at every layout tried (jobs 33128702, 33170265/6); rungs 1-4 carry the
# cost-vs-N extrapolation.
# FFTW_PLAN_TYPE ESTIMATE (CP2K default): MEASURE's one-off planning cost
# would land inside the few steps being timed. No accuracy setting touched.
set -euo pipefail

# Rungs 1-4 replicate an NNP-equilibrated cube_n1 (1 NaCl + 62 H2O, 12.42 A);
# rung 5 is the equilibrated production cube_n3 frame the RPA campaign started
# from. Unequilibrated frames cost extra SCF iterations (+-24% spread at rung
# 1, rung 2 non-convergent) and would overstate the AIMD cost.
ROOT=/rds/user/$USER/hpc-work/dft_cost_ladder

RUNG="${1:?usage: sbatch_dft_ladder.sh <rung 1-5>}"
case "$RUNG" in
  1) NREP="1 1 1"; NAT=188  ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1_equil.xyz"; MOL=0.895 ;;
  2) NREP="2 1 1"; NAT=376  ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1_equil.xyz"; MOL=0.895 ;;
  3) NREP="2 2 1"; NAT=752  ; STEPS=11; ABC=12.4200; COORD="$ROOT/cube_n1_equil.xyz"; MOL=0.895 ;;
  4) NREP="2 2 2"; NAT=1504 ; STEPS=6 ; ABC=12.4200; COORD="$ROOT/cube_n1_equil.xyz"; MOL=0.895 ;;
  5) NREP="1 1 1"; NAT=5064 ; STEPS=3 ; ABC=37.2600; COORD="$ROOT/cube_n3_equil.xyz"; MOL=1.000 ;;
  *) echo "rung must be 1-5" >&2; exit 2 ;;
esac
TIMED=$((STEPS - 1))   # step 1 runs from the ATOMIC guess and is discarded

TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_dft_cost
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
export BIN_LABEL=master
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
source "$ADIR/env_csd3.sh"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

DATA=/home/crm98/cp2k_master/data
RUNDIR="$ROOT/n${NAT}"
CSV=/home/crm98/cp2k-benchmarks/results/dft_cost/dft_ladder_timings.csv
proj="nacl_revpbe_d3_n${NAT}"

echo "node: $(hostname)  rung: $RUNG  atoms: $NAT  cell: $ABC A x $NREP  molality: $MOL" \
     " layout: $SLURM_NTASKS ranks x $OMP_NUM_THREADS threads = $((SLURM_NTASKS * OMP_NUM_THREADS)) cores  rs_grid: ${RS_GRID_MODE:-auto}"
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

sed -i 's|^  FFTW_PLAN_TYPE MEASURE$|  FFTW_PLAN_TYPE ESTIMATE|' "$RUNDIR/dft.inp"

# RS_GRID_MODE=replicated: full realspace grid per rank, rs->pw becomes a
# reduction instead of a halo exchange (set per rung by sbatch_dft_probe.sh).
if [ "${RS_GRID_MODE:-auto}" = replicated ]; then
  sed -i 's|^      CUTOFF 1200$|      CUTOFF 1200\n      \&RS_GRID\n        DISTRIBUTION_TYPE REPLICATED\n      \&END RS_GRID|' "$RUNDIR/dft.inp"
  grep -q "DISTRIBUTION_TYPE REPLICATED" "$RUNDIR/dft.inp" || { echo "FATAL: RS_GRID patch did not apply"; exit 8; }
fi

nat_in=$(grep -cE "^ *(O|H|Na|Cl) " "$COORD")
nrep_prod=$(awk -v n="$NREP" 'BEGIN{split(n,r," "); print r[1]*r[2]*r[3]}')
[ $((nat_in * nrep_prod)) -eq "$NAT" ] || {
  echo "FATAL: $COORD has $nat_in atoms x $nrep_prod replicas != $NAT"; exit 7; }

cd "$RUNDIR"
t0=$SECONDS
srun --ntasks="$SLURM_NTASKS" --cpus-per-task="${SLURM_CPUS_PER_TASK:-4}" \
     --hint=nomultithread "$BIN" -i dft.inp -o dft.out
wall=$((SECONDS - t0))

ener="$RUNDIR/$proj-1.ener"
[ -f "$ener" ] || { echo "FATAL: no .ener written"; exit 3; }

# steady-state s/step = mean over the last TIMED steps; step 1 pays for the
# ATOMIC-guess SCF and is excluded.
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
