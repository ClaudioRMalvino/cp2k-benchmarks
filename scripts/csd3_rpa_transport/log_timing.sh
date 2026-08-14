#!/usr/bin/env bash
# Append one row to the transport-campaign timings CSV (cost-proof bookkeeping
# behind the "campaign cost on master" extrapolation).
# usage: log_timing.sh <stage> <conc> <binary_label> <ranks> <steps> <wall_s> <s_per_step> [notes]
CSV=${TIMINGS_CSV:-/home/crm98/cp2k-benchmarks/results/rpa_transport/timings.csv}
mkdir -p "$(dirname "$CSV")"
[ -f "$CSV" ] || echo "date,job_id,stage,conc,binary,partition,nodes,ranks,omp,steps,wall_s,s_per_step,notes" > "$CSV"
echo "$(date +%F),${SLURM_JOB_ID:-login}${SLURM_ARRAY_TASK_ID:+_$SLURM_ARRAY_TASK_ID},$1,$2,$3,${SLURM_JOB_PARTITION:-},${SLURM_NNODES:-},$4,${OMP_NUM_THREADS:-1},$5,$6,$7,${8:-}" >> "$CSV"
