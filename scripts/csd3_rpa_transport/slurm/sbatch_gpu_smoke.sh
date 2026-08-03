#!/usr/bin/env bash
#SBATCH -J rpa_gpu_smoke
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_gpu_smoke_%j.out

# GPU-track smoke for the transport campaign: OUR 1 mol/kg cube3 deck
# (5064 atoms, 8-net RPA committee MIXED with FIST/SPME, per-step stress)
# on one A100 with Dhruv's collab/gpu-water-diffusion binary. Verifies the
# three things his 1008-atom water protocol does not (campaign brief):
#   1. MIXED force eval runs on this branch at our size,
#   2. the stress tensor is produced on the GPU path (no stress = no eta),
#   3. our deck/cell work; record s/step for the cost ladder.
# Plus the runbook-mandated virial validation: a 20-step USE_GPU ON vs OFF
# pair through the SAME binary, diffing stress and energies (the GPU virial
# path is new; expected agreement ~1e-9 relative on energies).
# Account is passed literally at submit: sbatch -A <gpu-acct> ... (lowercase).
set -uo pipefail
ADIR=/home/crm98/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor
TDIR=/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport

# ampere runtime env (modules, BIN, MKL/CUDA/install lib paths) — single
# source of truth shared with the GPU equil/production wrappers
source "$TDIR/env_csd3_gpu.sh"

echo "node: $(hostname)"; nvidia-smi -L 2>/dev/null || echo "(no nvidia-smi)"
"$BIN" --version | tr ' ' '\n' | grep -x nnp_gpu \
  && echo "OK  binary carries nnp_gpu" || { echo "FAIL: nnp_gpu flag missing"; exit 1; }

ANCHOR_ROOT="$HPC_WORK/nacl_mp2_anchor"
MODEL_DIR="$ANCHOR_ROOT/models/RPA/final_model"
COORD="$ANCHOR_ROOT/runs/RPA/cubic_1M/cube_n3.xyz"
RUNDIR="$ANCHOR_ROOT/runs/RPA/gpu_smoke"
NNP_RANKS=1   # one rank drives the GPU; 3 CPU ranks carry SPME
FIST_RANKS=$(( SLURM_NTASKS - NNP_RANKS ))
export CP2K_DATA_DIR=$HOME/cp2k_optimized/data

mkdeck () {  # $1=outfile $2=steps $3=USE_GPU value $4=project
  sed -e "s|__PROJECT__|$4|" \
      -e "s|__STEPS__|$2|" \
      -e "s|__SNAPSTEPS__|100000|" \
      -e "s|__MODEL__|$MODEL_DIR|" \
      -e "s|__COORD__|$COORD|" \
      -e "s|__ABCX__|37.2600|" -e "s|__ABCY__|37.2600|" -e "s|__ABCZ__|37.2600|" \
      -e "s|__MX__|1|" -e "s|__MY__|1|" -e "s|__MZ__|1|" \
      -e "s|__GMAXX__|120|" -e "s|__GMAXY__|120|" -e "s|__GMAXZ__|120|" \
      -e "s|__NNPRANKS__|$NNP_RANKS|" \
      -e "s|__FISTRANKS__|$FIST_RANKS|" \
      -e "s|MD 100$|MD 10|" \
      -e "s|^  &NNP$|  \&NNP\n    USE_GPU $3|" \
      "$ADIR/nacl_diffusion_template.inp" > "$1"
}

fail=0

# --- main smoke: 200 steps on the GPU ---------------------------------------
mkdir -p "$RUNDIR/main"; cd "$RUNDIR/main"; rm -f smoke.out NaCl_RPA_gpu*
mkdeck smoke.inp 200 ON NaCl_RPA_gpu
t0=$SECONDS
srun --ntasks="$SLURM_NTASKS" "$BIN" -i smoke.inp -o smoke.out
wall=$(( SECONDS - t0 ))

grep -q "PROGRAM ENDED" smoke.out && echo "OK  cp2k finished" || { echo "FAIL: no PROGRAM ENDED"; fail=1; }
grep -i "USE_GPU\|GPU" smoke.out | head -3
st=NaCl_RPA_gpu-1.stress
if [ -f "$st" ] && awk 'NR>1 {for(i=3;i<=NF;i++) if ($i+0!=0) {found=1; exit}} END{exit !found}' "$st"
then echo "OK  stress tensor nonzero on GPU path"; else
  echo "FAIL: stress missing/all-zero on GPU path (eta would be blocked)"; fail=1; fi
ce=$(ls *nnp-energies-1.data 2>/dev/null | head -1 || true)
[ -n "$ce" ] && echo "OK  committee energies written" || { echo "FAIL: no committee energies"; fail=1; }
[ -f NaCl_RPA_gpu-pos-1.xyz ] && [ -f NaCl_RPA_gpu-1.ener ] \
  && echo "OK  trajectory + ener streams" || { echo "FAIL: missing streams"; fail=1; }
e0=$(awk 'NR==2 {print $5}' NaCl_RPA_gpu-1.ener)
echo "step-0 potential energy [Ha]: $e0"
echo "   (CPU dhruv binary + master both gave: -30777.866694782)"
tstep=$(awk '/qs_mol_dyn_low/ {print $NF; exit}' smoke.out || true)
if [ -n "${tstep:-}" ]; then
  ps=$(awk -v t="$tstep" 'BEGIN{printf "%.4f", t/200}')
  echo "GPU t/step = $ps s  (1x A100 + $FIST_RANKS SPME ranks, cube3/5064 atoms)"
  awk -v p="$ps" 'BEGIN{printf "  80 ps segment: %5.1f h  (CPU node: 11.4 h)\n", p*160000/3600}'
  bash "$TDIR/log_timing.sh" gpu_smoke cubic_1M gpu-collab-branch "$SLURM_NTASKS" 200 "$wall" "$ps" "cube3 RPA, 1xA100+3 SPME ranks"
fi

# --- virial validation: 20-step USE_GPU ON vs OFF, same binary --------------
for mode in on off; do
  d="$RUNDIR/virial_$mode"; mkdir -p "$d"; cd "$d"; rm -f v.out NaCl_RPA_v*
  mkdeck v.inp 20 "$(echo $mode | tr a-z A-Z)" NaCl_RPA_v
  srun --ntasks="$SLURM_NTASKS" "$BIN" -i v.inp -o v.out
  grep -q "PROGRAM ENDED" v.out || { echo "FAIL: virial $mode run died"; fail=1; }
done
if [ $fail -eq 0 ]; then
  echo "--- GPU-vs-CPU virial comparison (20 steps, same binary) ---"
  paste "$RUNDIR/virial_on/NaCl_RPA_v-1.stress" "$RUNDIR/virial_off/NaCl_RPA_v-1.stress" \
    | awk 'NR>1 { n=(NF/2); for(i=3;i<=n;i++){ d=$i-$(i+n); if(d<0)d=-d; a=($i<0?-$i:$i);
                  if(d>md)md=d; if(a>ma)ma=a } } END{ printf "max|dS| = %.3e bar   max|S| = %.3e bar   rel = %.3e\n", md, ma, md/ma }'
  paste "$RUNDIR/virial_on/NaCl_RPA_v-1.ener" "$RUNDIR/virial_off/NaCl_RPA_v-1.ener" \
    | awk '!/^#/ { n=(NF/2); d=$5-$(5+n); if(d<0)d=-d; if(d>md)md=d } END{ printf "max|dE_pot| over 20 steps = %.3e Ha\n", md }'
  echo "(runbook expectation: ~1e-9 level; large deviations mean the GPU virial path is NOT trustworthy for eta)"
fi

[ $fail -eq 0 ] && echo "GPU SMOKE PASSED (inspect virial numbers above)" || { echo "GPU SMOKE FAILED"; exit 1; }
