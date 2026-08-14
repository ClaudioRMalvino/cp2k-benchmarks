#!/usr/bin/env bash
# Runaway watchdog for Fig. S4 jobs (login-node, free).
# Usage: nohup ./watchdog_figS4.sh JOBID=/run/dir [...] > watchdog.log 2>&1 &
# Polls every POLL_S; scancels a job when latest T > KILL_T (default 1500 K) or
# a single-step Pot jump >0.1 Ha appears. Exits when none left or after MAX_HOURS.

set -u
KILL_T=${KILL_T:-1500}
POLL_S=${POLL_S:-60}
MAX_HOURS=${MAX_HOURS:-24}

(( $# > 0 )) || { echo "usage: $0 JOBID=/run/dir [...]"; exit 1; }

declare -A DIR
for arg in "$@"; do
    DIR[${arg%%=*}]=${arg#*=}
done

ts() { date "+%F %H:%M:%S"; }
deadline=$(( $(date +%s) + MAX_HOURS*3600 ))

echo "[$(ts)] watchdog: ${!DIR[*]} (KILL_T=${KILL_T} K, poll ${POLL_S}s)"
while (( $(date +%s) < deadline )); do
    alive=0
    for jid in "${!DIR[@]}"; do
        state=$(squeue -h -j "$jid" -o %T 2>/dev/null || true)
        [[ -z "$state" ]] && continue
        alive=1
        ener=$(ls -t "${DIR[$jid]}"/*.ener 2>/dev/null | head -1)
        [[ -z "$ener" ]] && { echo "[$(ts)] $jid $state (no .ener yet)"; continue; }
        read -r T jump <<< "$(awk 'NR>1 { T=$4; if (NR>2) { d=$5-p; if (d<0) d=-d; if (d>0.1) j++ } p=$5 }
                                  END { printf "%.1f %d", T+0, j+0 }' "$ener")"
        if (( $(echo "$T > $KILL_T" | bc -l) )) || (( jump > 0 )); then
            echo "[$(ts)] $jid RUNAWAY (T=$T K, jumps=$jump) -- scancel"
            scancel "$jid"
        else
            echo "[$(ts)] $jid $state  T=$T K  jumps=$jump  ($ener)"
        fi
    done
    (( alive == 0 )) && { echo "[$(ts)] no monitored jobs left -- exiting"; exit 0; }
    sleep "$POLL_S"
done
echo "[$(ts)] MAX_HOURS reached -- exiting"
