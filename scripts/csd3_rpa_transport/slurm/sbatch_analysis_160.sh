#!/usr/bin/env bash
#SBATCH -J rpa_ana160
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:30:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/rpa_ana160_%j.out

# Round-2 transport analysis over the 5 x 160 ps CPU campaign (usage: sbatch).
# Reads runs/RPA/{cubic_1M,cubic_2m,cubic_4m}/production/cube3/seg{1..5};
# writes *_5x160_* products to results/rpa_transport/ (round-1 *_5x80_* untouched).
# Batch not login-node: phase A streams 15 x 10.8 GB pos files (162 GB total).
# 16 cores = widest phase (A, 15 single-threaded procs; icelake bills reserved cores).
# Pilot 2026-08-10: 84 s + 624 MB RSS per segment; himem chosen for queue time, not memory.

# No set -e: one concentration failing must not throw away the other two.
set -uo pipefail

BENCH=/home/crm98/cp2k-benchmarks
S="$BENCH/scripts/csd3_rpa_transport"
RR=/rds/user/crm98/hpc-work/nacl_mp2_anchor/runs/RPA
RES="$BENCH/results/rpa_transport"
BOX=37.26          # all three cells are the same cube3 box; only N differs
EVERY=20           # 20 x 5 fs dumps = 0.1 ps sampling, matches round 1
LOGD="$BENCH/logs/rpa_ana160_${SLURM_JOB_ID:-manual}"
mkdir -p "$LOGD" "$RES"

# single-threaded workers; parallelism is across segments/concentrations
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

source ~/.fortran_env/bin/activate

echo "node      : $(hostname)"
echo "started   : $(date)"
echo "logs      : $LOGD"
echo "results   : $RES"
echo

CONCS="cubic_1M cubic_2m cubic_4m"
# label stems follow the round-1 convention (1 m has no molality token in kappa/eta labels)
stem_kappa() { case "$1" in cubic_1M) echo rpa_cpu_cube3    ;; cubic_2m) echo rpa_cpu_cube3_2m ;; cubic_4m) echo rpa_cpu_cube3_4m ;; esac; }
stem_pair()  { case "$1" in cubic_1M) echo rpa_cpu_cube3_1M ;; cubic_2m) echo rpa_cpu_cube3_2m ;; cubic_4m) echo rpa_cpu_cube3_4m ;; esac; }

segdirs() { for s in 1 2 3 4 5; do echo "$RR/$1/production/cube3/seg$s"; done; }

# ---------------------------------------------------------------- preflight
# every segment must be at step 320000 = 160 ps
fail=0
for c in $CONCS; do
  for s in 1 2 3 4 5; do
    d="$RR/$c/production/cube3/seg$s"
    e="$d/NaCl_RPA_cube3_seg${s}-1.ener"
    [ -f "$e" ] || { echo "PREFLIGHT FAIL: missing $e"; fail=1; continue; }
    last=$(awk 'END{print $1}' "$e")
    [ "$last" = "320000" ] || { echo "PREFLIGHT FAIL: $c seg$s last step $last != 320000"; fail=1; }
    [ -f "$d/NaCl_RPA_cube3_seg${s}-pos-1.xyz" ] || { echo "PREFLIGHT FAIL: no pos file in $d"; fail=1; }
    [ -f "$d/NaCl_RPA_cube3_seg${s}-1.stress" ]  || { echo "PREFLIGHT FAIL: no stress file in $d"; fail=1; }
  done
done
[ "$fail" -eq 0 ] || { echo "preflight failed - nothing submitted for analysis"; exit 3; }
echo "preflight OK: 15/15 segments at step 320000 with pos + stress present"
echo

# ------------------------------------------------- phase A: rebuild caches
cat > "$LOGD/rebuild_cache.py" <<'PY'
import os, sys, time
sys.path.insert(0, "/home/crm98/cp2k-benchmarks/scripts/csd3_rpa_transport")
from analyze_onsager_cp2k import reduce_segment

sd, every = sys.argv[1], int(sys.argv[2])
tag = f"{os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(sd))))}/{os.path.basename(sd)}"
t0 = time.time()
# Unconditional --refresh: every cache is rebuilt from the whole pos file,
# start to finish. The extension appended to the same trajectories, so there
# is no "new part" to splice on - the 160 ps analysis simply re-reads all of
# it. steps 0..320000 in the line below is the check that that happened.
steps, codes, _ = reduce_segment(sd, every, refresh=True)
ps = (steps[-1] - steps[0]) * 0.0005          # 0.5 fs per MD step
print(f"{tag}: {len(steps)} frames  {len(codes)} atoms  "
      f"steps {steps[0]}..{steps[-1]}  {ps:.1f} ps  [{time.time()-t0:.0f} s]",
      flush=True)
PY

echo "=== phase A: rebuilding 15 reduced-trajectory caches ==="
tA=$SECONDS
for c in $CONCS; do
  for d in $(segdirs "$c"); do
    ( python "$LOGD/rebuild_cache.py" "$d" "$EVERY" \
        > "$LOGD/cacheA_${c}_$(basename "$d").log" 2>&1 \
        || echo "CACHE FAIL: $d" ) &
  done
done
wait
grep -h . "$LOGD"/cacheA_*.log
nok=$(grep -l "frames" "$LOGD"/cacheA_*.log 2>/dev/null | wc -l)
echo "phase A done in $((SECONDS-tA)) s: $nok/15 caches rebuilt"
[ "$nok" -eq 15 ] || { echo "FATAL: not all caches rebuilt - aborting before analysis"; exit 4; }
echo

# ------------------------------- phases B, C, D run concurrently from here
echo "=== phases B (Onsager) / C (viscosity) / D (pairing) ==="
tB=$SECONDS

# ---- B: Onsager. Windows serialised per concentration (shared .dat file).
for c in $CONCS; do
(
  k=$(stem_kappa "$c")
  for w in "10 30" "10 40" "10 60"; do
    set -- $w
    lab="${k}_5x160_w$1$2"
    python "$S/analyze_onsager_cp2k.py" \
      --segdirs $(segdirs "$c") --box-a "$BOX" --temp 300.0 \
      --every "$EVERY" --tmin "$1" --tmax "$2" --traj-max-ps 160.0 \
      --label "$lab" > "$LOGD/onsager_${lab}.log" 2>&1 \
      || echo "ONSAGER FAIL: $lab"
  done
) &
done

# ---- C: viscosity. *_final: acf 40 / plateau 10-18 ps (what 160 ps buys);
#         *_p0210: acf 20 / plateau 2-10 ps (matches round 1, like-for-like).
for c in $CONCS; do
(
  k=$(stem_kappa "$c")
  python "$S/analyze_viscosity_cp2k.py" \
    --segdirs $(segdirs "$c") --box-a "$BOX" --temp 300.0 --dt-fs 0.5 \
    --acf-max-ps 40 --plateau-lo 10 --plateau-hi 18 \
    --label "${k}_5x160_final" > "$LOGD/eta_${k}_5x160_final.log" 2>&1 \
    || echo "ETA FAIL: ${k}_5x160_final"
  python "$S/analyze_viscosity_cp2k.py" \
    --segdirs $(segdirs "$c") --box-a "$BOX" --temp 300.0 --dt-fs 0.5 \
    --acf-max-ps 20 --plateau-lo 2 --plateau-hi 10 \
    --label "${k}_5x160_p0210" > "$LOGD/eta_${k}_5x160_p0210.log" 2>&1 \
    || echo "ETA FAIL: ${k}_5x160_p0210"
) &
done

# ---- D: pairing / PMF. Reads the phase-A caches only.
for c in $CONCS; do
(
  p=$(stem_pair "$c")
  lab="${p}_5x160"
  python "$S/analyze_ion_pairing_cp2k.py" \
    --segdirs $(segdirs "$c") --box-a "$BOX" --temp 298.15 \
    --every "$EVERY" --tmax-ps 160.0 \
    --label "$lab" > "$LOGD/pairing_${lab}.log" 2>&1 \
    || echo "PAIRING FAIL: $lab"
) &
done

wait
echo "phases B/C/D done in $((SECONDS-tB)) s"
echo

# ------------------------------------------------------------- transcripts
for f in "$LOGD"/onsager_*.log "$LOGD"/eta_*.log "$LOGD"/pairing_*.log; do
  [ -f "$f" ] || continue
  echo "################################################## $(basename "$f")"
  cat "$f"
  echo
done

# --------------------------------------------------------------- manifest
MAN="$RES/MANIFEST_5x160.txt"
{
  echo "# round-2 transport products, 5 x 160 ps CPU campaign"
  echo "# job ${SLURM_JOB_ID:-manual}, generated $(date -Iseconds)"
  echo "# source: $RR/{cubic_1M,cubic_2m,cubic_4m}/production/cube3/seg{1..5} @ step 320000"
  echo
  ls -la "$RES"/*5x160* 2>/dev/null
} > "$MAN"

echo "=== products written ==="
ls -la "$RES"/*5x160* 2>/dev/null || echo "NONE - check the transcripts above"
echo
echo "manifest: $MAN"
echo "finished: $(date)"
