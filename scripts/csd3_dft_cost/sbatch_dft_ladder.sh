#!/usr/bin/env bash
#SBATCH -J dft_ladder
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/dft_ladder_%j.out

# On-the-fly DFT cost ladder: one rung per invocation.
#
#   usage: sbatch --time=<wall> sbatch_dft_ladder.sh <rung 1-4>
#
#   RUNG 4 MUST BE SUBMITTED WITH -p icelake-himem:
#     sbatch -p icelake-himem --time=04:00:00 sbatch_dft_ladder.sh 4
#
# RUNG 4 DOES NOT FIT IN A PLAIN icelake NODE. Job 33197088 converged the
# step-0 wavefunction cleanly (7 SCF steps, outer loop converged in 2,
# E = -9065.61728603745905 Ha) and was then OOM-killed on rank 35 during the
# first force evaluation; the job hung on the survivors until the 5 h wall
# and wrote no .ener at all. sacct puts MaxRSS at 3,911,944 K per rank, so
# 76 ranks need ~283 GB against the 250 GB an icelake node has - it is short
# by about 13%, not by an order of magnitude. icelake-himem is the same Ice
# Lake silicon with the same 76 cores and 502 GB, so the 76x1 pure-MPI layout
# the whole ladder depends on is preserved exactly and only the DRAM ceiling
# moves. Hybrid 19x4 would also have fitted, but it changes the layout the
# NNP reference rungs were measured on, which is precisely what this ladder
# is not allowed to do.
#
# Because rung 4 now sits on different nodes from rungs 1-3, rung 3 is re-run
# on himem as a control - if it reproduces its icelake 250.2436 s/step then
# the partitions are timing-equivalent and rung 4 can be fitted alongside
# rungs 1-3 without an asterisk. The CSV records SLURM_JOB_PARTITION per row,
# so both rung-3 measurements coexist and are told apart by that column.
#
# Measures the steady-state Born-Oppenheimer AIMD cost per MD step for
# liquid water at revPBE-D3 (the O'Neill training-set level), at four system
# sizes built by replicating one 64-water box. The fitted cost-vs-N exponent
# extrapolates to the 5064-atom production cell, giving the "what this
# campaign would have cost with on-the-fly DFT" number that sits above the
# measured master / optimised-CPU / GPU NNP rungs in the cost ladder.
#
# Run on icelake with the same 1-node / 76-core footprint and the same pristine
# master binary (757bb76a80) as the master NNP reference, so DFT and NNP
# numbers are directly comparable.
#
# LAYOUT: 76 pure-MPI ranks on one node, deliberately - that is exactly the
# layout the NNP reference rungs were measured on (master 4.0473 s/step,
# dhruv-cell-list 0.2554 s/step, both 1 node / 76 ranks), and the ladder only
# means anything if the DFT and NNP numbers describe the same machine
# footprint. Consistency with what we are comparing against beats internal
# consistency with rung 5.
#
# RUNG 5 IS NOT SUBMITTED THROUGH THIS PATH. The 5064-atom cell hangs before
# the first SCF iteration - job 33128702 burned 6 h 25 min at 76x1 with every
# rank inside transfer_rs2pw_distributed collocating the core density, and the
# 1 h probes 33170265 / 33170266 reproduced it exactly at 19 ranks x 4 threads
# and again with a replicated realspace grid. Three layouts, identical stall,
# so it is not a slab/halo ratio. It is tracked separately; rungs 1-4 carry the
# extrapolation, which is what the cost argument needs.
#
# FFTW plans with ESTIMATE rather than MEASURE: MEASURE times trial plans for
# every FFT size, and on the larger rungs that one-off charge lands inside the
# handful of steps being timed. ESTIMATE is CP2K's own default. Neither
# setting changes a computed number, and no accuracy setting is touched.
set -euo pipefail

# The ladder is NaCl(aq) throughout - the campaign's own cell lineage, not a
# water proxy. Rungs 1-4 replicate cube_n1 (1 NaCl + 62 H2O, 12.42 A,
# 0.895 mol/kg) so molality is constant across the fitted rungs; rung 5 is
# the PRODUCTION cube_n3 configuration itself (30 NaCl + 1668 H2O, 37.26 A,
# 1.000 mol/kg) - the same coordinates the RPA campaign started from, so the
# headline number is measured on the system we actually study rather than
# extrapolated to it. Rung 5 gets fewer MD steps because each costs hours.
# CONFIGURATIONS MUST BE EQUILIBRATED. The builder's substituted cells place
# ions on former water-oxygen sites with no solvation shell; that is a strained
# electronic structure which costs extra SCF iterations (rung 1 measured a
# +-24% step-time spread) or fails to converge outright (rung 2 aborted after
# 21 outer-SCF iterations). Timing a strained frame would OVERSTATE the AIMD
# cost, biasing the comparison in the MLIP's favour - not acceptable.
# So: rung 5 uses the equilibrated frame production actually started from
# (t=0 of seg1, i.e. equil snapshot_1 after 30 ps), and rungs 1-4 use an
# NNP-equilibrated cube_n1.
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

# MEASURE times trial plans for every FFT size; on the large rungs that one-off
# charge lands inside the steps we are trying to time. ESTIMATE is CP2K's own
# default and costs nothing up front.
sed -i 's|^  FFTW_PLAN_TYPE MEASURE$|  FFTW_PLAN_TYPE ESTIMATE|' "$RUNDIR/dft.inp"

# RS_GRID_MODE=replicated: each rank holds the full realspace grid and the
# rs->pw step becomes a plain reduction instead of a halo exchange. Set by the
# probe result, per rung - see sbatch_dft_probe.sh.
if [ "${RS_GRID_MODE:-auto}" = replicated ]; then
  sed -i 's|^      CUTOFF 1200$|      CUTOFF 1200\n      \&RS_GRID\n        DISTRIBUTION_TYPE REPLICATED\n      \&END RS_GRID|' "$RUNDIR/dft.inp"
  grep -q "DISTRIBUTION_TYPE REPLICATED" "$RUNDIR/dft.inp" || { echo "FATAL: RS_GRID patch did not apply"; exit 8; }
fi

# the deck must describe the system we think it does
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
