#!/usr/bin/env bash
# Watchdog for the three Fig.S4 replication jobs (N=128, N=256, N=512).
# Polls each job's .ener file once a minute. If max T exceeds KILL_T (2000K),
# scancels the job to stop wasting compute on a runaway. Otherwise just
# reports the current state. Exits when all three jobs are no longer running.
#
# Designed to be launched with nohup so it survives session disconnect:
#   nohup bash watchdog_figS4.sh > logs/watchdog.log 2>&1 &

set -u

# Job IDs and the result-dir prefix they each populate.
declare -A JOB2GLOB
JOB2GLOB[29496515]="/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N128_*"
JOB2GLOB[29496518]="/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N256_*"
JOB2GLOB[29496519]="/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/replicate_N512_*"

KILL_T=2000.0      # K — only kill on unambiguous runaway
WARN_T=800.0       # K — log a warning above this
POLL_S=60          # seconds between checks (cheap; worst-case ~1 CPU-hr
                   # wasted on N=512 if a runaway slips past one tick)
MAX_HOURS=6        # hard cap so the watchdog can't loop forever

ts() { date "+%F %H:%M:%S"; }
log() { echo "[$(ts)] $*"; }

# Pick newest matching rundir for each job (in case there are stale ones).
latest_dir()  { ls -td $1 2>/dev/null | head -1; }
latest_ener() { ls -t $1/*.ener 2>/dev/null | head -1; }

# Quick liveness check via squeue (avoid sacct rate limiting).
job_state() {
    local s
    s=$(squeue -j "$1" -h -o '%T' 2>/dev/null)
    [[ -z "$s" ]] && echo "DONE" || echo "$s"
}

log "watchdog start. KILL_T=${KILL_T}K  POLL=${POLL_S}s  jobs: ${!JOB2GLOB[*]}"

deadline=$(( SECONDS + MAX_HOURS*3600 ))
while [ "$SECONDS" -lt "$deadline" ]; do
    all_done=1
    for jid in "${!JOB2GLOB[@]}"; do
        state=$(job_state "$jid")
        if [[ "$state" == "DONE" ]]; then
            continue
        fi
        all_done=0
        if [[ "$state" != "RUNNING" ]]; then
            log "  job $jid : $state"
            continue
        fi
        rundir=$(latest_dir "${JOB2GLOB[$jid]}")
        if [[ -z "$rundir" ]]; then
            log "  job $jid : RUNNING (no rundir yet)"
            continue
        fi
        ener=$(latest_ener "$rundir")
        if [[ -z "$ener" ]]; then
            log "  job $jid : RUNNING (no .ener yet; stage is GEO_OPT or NVT just starting)"
            continue
        fi
        # Parse: col 1 = step, col 4 = T(K), col 6 = ConsQty
        read maxT lastT lastStep <<< "$(awk '
            NR>1 && NF>=6 { if ($4>m) m=$4; lT=$4; lS=$1 }
            END { printf "%.1f %.1f %s", m, lT, lS }' "$ener")"
        msg="  job $jid : RUNNING  step=$lastStep  T_now=${lastT}K  T_max=${maxT}K"
        if (( $(echo "$maxT > $KILL_T" | bc -l 2>/dev/null) )); then
            log "$msg  *** RUNAWAY (max T > ${KILL_T} K) — scancel ***"
            scancel "$jid"
            continue
        fi
        if (( $(echo "$maxT > $WARN_T" | bc -l 2>/dev/null) )); then
            log "$msg  (above WARN_T=${WARN_T}K, watching closely)"
        else
            log "$msg"
        fi
    done
    [ "$all_done" -eq 1 ] && { log "all jobs finished — watchdog exits"; break; }
    sleep "$POLL_S"
done

if [ "$SECONDS" -ge "$deadline" ]; then
    log "watchdog hit ${MAX_HOURS}h cap — exiting without killing anything"
fi
