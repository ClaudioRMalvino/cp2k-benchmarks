#!/usr/bin/env bash
#SBATCH -J rpa_master_ref
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
#SBATCH --time=02:30:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_master_ref_%j.out

# Transport campaign, run 2: MASTER-BRANCH reference timing at 1 mol/kg.
# Pristine upstream 2026.1 (commit 757bb76a80), same RPA cube3 deck, same
# node type and rank split as production. 1000 MD steps for a stable s/step;
# no physics needed from this run. This single number x the campaign step
# count = the "what this campaign would have cost on master" extrapolation,
# and its step-0 potential energy is the cross-binary consistency check
# against the smoke test (PR #5295 engine, same coordinates).
set -euo pipefail
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_mp2_transport
export BIN_LABEL=master
source "$ADIR/env_csd3.sh"

MODEL=RPA
ANCHOR_ROOT=/rds/user/$USER/hpc-work/nacl_mp2_anchor
MODEL_DIR="$ANCHOR_ROOT/models/$MODEL/final_model"
COORD="$ANCHOR_ROOT/runs/$MODEL/cubic_1M/cube_n3.xyz"
RUNDIR="$ANCHOR_ROOT/runs/$MODEL/master_ref"
STEPS=1000
RANKS=$SLURM_NTASKS
FIST_RANKS=8
proj="NaCl_${MODEL}_master_ref"

# provenance without exec'ing the binary (direct exec aborts in MPI_Init -
# the PMI gotcha; the binary must only ever be launched via srun)
echo "node: $(hostname)  binary: $BIN  model: $MODEL"
cat "$BIN_ROOT/$BIN_LABEL/PROVENANCE.txt" 2>/dev/null || true; echo
mkdir -p "$RUNDIR"; rm -f "$RUNDIR"/ref.out "$RUNDIR/$proj"*

sed -e "s|__PROJECT__|$proj|" \
    -e "s|__STEPS__|$STEPS|" \
    -e "s|__SNAPSTEPS__|100000|" \
    -e "s|__MODEL__|$MODEL_DIR|" \
    -e "s|__COORD__|$COORD|" \
    -e "s|__ABCX__|37.2600|" -e "s|__ABCY__|37.2600|" -e "s|__ABCZ__|37.2600|" \
    -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
    -e "s|__GMAXX__|120|" -e "s|__GMAXY__|120|" -e "s|__GMAXZ__|120|" \
    -e "s|__NNPRANKS__|$(( RANKS - FIST_RANKS ))|" \
    -e "s|__FISTRANKS__|$FIST_RANKS|" \
    "$ADIR/nacl_diffusion_template.inp" > "$RUNDIR/ref.inp"

cd "$RUNDIR"
t0=$SECONDS
srun --ntasks="$RANKS" --cpus-per-task=1 --hint=nomultithread \
  "$BIN" -i ref.inp -o ref.out
wall=$(( SECONDS - t0 ))

grep -q "PROGRAM ENDED" ref.out && echo "OK  cp2k finished" || { echo "FAIL"; exit 1; }

e0=$(awk 'NR==2 {print $5}' "$proj-1.ener")
echo "step-0 potential energy [Ha]: $e0   (compare against smoke)"

tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' ref.out)
ps=$(awk -v t="$tstep" -v s="$STEPS" 'BEGIN{printf "%.4f", t/s}')
echo "MASTER t/step = $ps s  ($RANKS ranks, cube3/5064 atoms, $MODEL)"
awk -v p="$ps" 'BEGIN{
  printf "  implied per 160 ps segment: %6.1f h;  15-segment campaign: %7.0f node-h\n",
         p*320000/3600, 15*p*320000/3600 }'
bash "$TDIR/log_timing.sh" master_ref cubic_1M master "$RANKS" "$STEPS" "$wall" "$ps" "cube3 RPA deck, pristine 2026.1"
echo "MASTER REFERENCE DONE"
