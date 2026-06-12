#!/usr/bin/env bash
# Probe: is MAQAO LProf usable on icelake compute nodes again after the
# zero-day mitigation?  Runs a 20-step N=64 NNP MD under maqao lprof and
# reports whether sampling worked.  Cost: 16 cores x ~5 min.
#
# Submit with:  sbatch probe_maqao_available.sh

#SBATCH -J maqao_probe
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:15:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/maqao_probe_%j.out

. /etc/profile.d/modules.sh
module purge
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
set +e

MAQAO=/home/crm98/maqao.x86_64.2026.0.0-b/bin/maqao
EXE=/home/crm98/cp2k_optimized/install/bin/cp2k.psmp
BENCH=/home/crm98/cp2k-benchmarks
RUNDIR=/rds/user/$USER/hpc-work/cp2k-benchmarks/maqao/probe_$(date +%H%M%S)
mkdir -p "$RUNDIR" && cd "$RUNDIR"

echo "=== compute-node security knobs ==="
echo "  perf_event_paranoid : $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null)"
echo "  ptrace_scope        : $(cat /proc/sys/kernel/yama/ptrace_scope 2>/dev/null)"

# Tiny dataset: 20-step N=64 from the standard benchmark input.
sed -e 's/STEPS 100/STEPS 20/' "$BENCH/H2O-64_NNP_MD.inp" > probe.inp
ln -sfn /home/crm98/cp2k_optimized/data/NNP NNP

# MAQAO must run ONCE and own the MPI launch (one srun inside maqao), not
# the other way round: srun-ning N maqao instances races on the experiment
# directory and hangs the MPI world.
try_lprof () {
    local label=$1 xp=$2; shift 2
    echo "=== maqao lprof ($label) ==="
    rm -rf "$xp"
    "$MAQAO" lprof --mpi-command="srun --ntasks=16 --cpus-per-task=1 --hint=nomultithread" \
        xp="$xp" "$@" -- "$EXE" -i probe.inp -o probe.out
    local rc=$?
    echo "    exit code: $rc"
    if [[ -d "$xp" ]] && "$MAQAO" lprof -df xp="$xp" 2>/dev/null | grep -qi "cp2k\|qs_\|nnp"; then
        echo "    -> $label: SAMPLING WORKED (function-level profile present)"
        return 0
    fi
    echo "    -> $label: no usable profile"
    return 1
}

echo
echo "=== VERDICT ==="
if try_lprof "default engine" probe_xp; then
    echo "  MAQAO LPROF WORKS on compute nodes — full ONE View profiling is back on the menu."
elif try_lprof "no-perf engine (time-based)" probe_xp_noperf engine=no-perf; then
    echo "  MAQAO works in TIME-BASED mode only (kernel counters still blocked)."
    echo "  Hotspot shares (Fig 4.7 style) are fine; use 'perf stat' for the counter table (Fig 4.8)."
else
    echo "  MAQAO LPROF STILL BLOCKED in both modes.  Fall back to 'perf stat' counters"
    echo "  + CP2K's internal timing report."
fi
echo "rundir: $RUNDIR"
