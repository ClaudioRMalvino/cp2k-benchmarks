#!/usr/bin/env bash
# Multi-node strong scaling: master vs feature/nnp-chebyshev across
# 1 / 2 / 4 Peta4-IceLake nodes (76 ranks per node, pure MPI), at two
# fixed system sizes.  Measures the cross-node MPI_Allreduce wall that
# Report 1 SS5.2 predicted but never measured past one socket.
#
# One allocation of MAX nodes; each sweep point sruns a SUBSET of the
# allocated nodes, so node count varies without resubmitting.  Same
# methodology as run_nnp_core_scaling_slurm.sh: N_REPS timed reps + 1
# discarded warm-up, time/step from qs_mol_dyn_low, CSV per branch.
#
# Submit with:  sbatch run_nnp_multinode_scaling_slurm.sh
# Budget: 4 nodes x <=2 h = <=2432 CPU-h worst case; realistically ~300.

#SBATCH -J NNP_multinode
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=76
#SBATCH --time=02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_multinode_%j.out

CORES_PER_NODE=76
BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BENCHMARK_ROOT=/home/crm98/cp2k-benchmarks

. /etc/profile.d/modules.sh
module purge
set +u
source "$BENCHMARK_ROOT/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh"
source "$BIN_ROOT/setup"
set -u
set +e

TIMESTAMP=$(date +%d-%m_%H-%M-%S)${SLURM_JOB_ID:+_${SLURM_JOB_ID}}
NODE_LIST="${NODE_LIST:-1 2 4}"
SIZE_LIST="${SIZE_LIST:-1024 4096}"
STEPS=${STEPS:-100}   # 100 MD steps: uniform across every Report-2 measurement
N_REPS=${N_REPS:-3}
BASE_INP="$BENCHMARK_ROOT/H2O-64_NNP_MD.inp"

declare -A MULT
MULT[1024]="4 2 2"
MULT[2048]="4 4 2"
MULT[4096]="4 4 4"

stats() {
   printf '%s\n' "$@" | awk '
      /^[0-9.eE+-]+$/ { x[++n]=$1; s+=$1; if (n==1 || $1<mn) mn=$1 }
      END { if (n==0) { print "NA NA NA 0"; exit }
            m=s/n; v=0; for (i=1;i<=n;i++) v+=(x[i]-m)^2
            sd=(n>1)?sqrt(v/(n-1)):0
            printf "%.6f %.6f %.6f %d\n", m, sd, mn, n }'
}

make_input() {   # $1=target  $2=size
   export target_file="$1" current_size="$2" base_inp="$BASE_INP" STEPS
   export mx=$(awk '{print $1}' <<<"${MULT[$2]}") \
          my=$(awk '{print $2}' <<<"${MULT[$2]}") \
          mz=$(awk '{print $3}' <<<"${MULT[$2]}")
python3 - <<'PYEOF'
import os, re
txt = open(os.environ['base_inp']).read()
txt = re.sub(r'STEPS\s+\d+', f"STEPS {os.environ['STEPS']}", txt, count=1)
txt = re.sub(r'&TRAJECTORY\b[\s\S]*?&END TRAJECTORY',
             '&TRAJECTORY OFF\n    &END TRAJECTORY\n'
             '    &RESTART OFF\n    &END RESTART\n'
             '    &RESTART_HISTORY OFF\n    &END RESTART_HISTORY', txt, count=1)
txt = re.sub(r'&ENERGIES\b', '&ENERGIES OFF', txt, count=1)
txt = re.sub(r'&FORCES\b', '&FORCES OFF', txt, count=1)
mx, my, mz = os.environ['mx'], os.environ['my'], os.environ['mz']
txt = re.sub(r'(&CELL\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.IGNORECASE)
if "&TOPOLOGY" in txt:
    txt = re.sub(r'(&TOPOLOGY\s*\n)', rf'\g<1>      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n', txt, flags=re.IGNORECASE)
else:
    txt = re.sub(r'(&SUBSYS\s*\n)', rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL {mx} {my} {mz}\n    &END TOPOLOGY\n', txt, flags=re.IGNORECASE)
open(os.environ['target_file'], 'w').write(txt)
PYEOF
}

# BRANCH_LIST override lets a rerun measure one branch only, e.g.
#   BRANCH_LIST=feature-nnp-chebyshev sbatch run_nnp_multinode_scaling_slurm.sh
# NOTE: the 10-06 multinode data (job 30351703, both branches) was taken at
# STEPS=50; Report 2 standardises on STEPS=100, so the default both-branch
# run is the right one unless you are only refreshing chebyshev.
for branch in ${BRANCH_LIST:-master feature-nnp-chebyshev}; do
   case "$branch" in
      master) LABEL=upstream-master; PROJECT_ROOT=/home/crm98/cp2k_master; OUTDIR_PARENT=cp2k_master ;;
      *)      LABEL=$branch; PROJECT_ROOT=/home/crm98/cp2k_optimized; OUTDIR_PARENT=cp2k_feature_chebyshev ;;
   esac
   CP2K_EXE="$BIN_ROOT/$branch/cp2k.psmp"
   [[ -x "$CP2K_EXE" ]] || { echo "!! missing $CP2K_EXE"; continue; }
   export LD_LIBRARY_PATH="$BIN_ROOT/$branch/lib:${LD_LIBRARY_PATH:-}"
   export OMP_NUM_THREADS=1

   OUTDIR="/rds/user/$USER/hpc-work/cp2k-benchmarks/results/${OUTDIR_PARENT}/NNP/NNP_multinode_${LABEL}_${TIMESTAMP}"
   mkdir -p "$OUTDIR"
   CSV_FILE="$OUTDIR/results_multinode_${LABEL}_${TIMESTAMP}.csv"
   {
      echo "# branch: $LABEL"
      echo "# exe:    $CP2K_EXE ($(md5sum < "$CP2K_EXE" | awk '{print $1}'))"
      echo "# steps:  $STEPS   reps: $N_REPS (+1 warm-up)"
      echo "# n_molecules,nodes,mpi_ranks,n_reps,time_per_step_mean_s,time_per_step_std_s,time_per_step_min_s,speedup_vs_1node,parallel_efficiency_pct"
   } >"$CSV_FILE"

   for size in $SIZE_LIST; do
      BASE_TPS=""
      for nodes in $NODE_LIST; do
         mpi=$(( nodes * CORES_PER_NODE ))
         rundir="$OUTDIR/N${size}_nodes${nodes}"
         mkdir -p "$rundir"
         make_input "$rundir/run.inp" "$size"
         ln -sfn "$PROJECT_ROOT/data/NNP" "$rundir/NNP"

         tps_list=""
         for r in $(seq 0 "$N_REPS"); do
            out="$rundir/cp2k_rep${r}.out"
            ( cd "$rundir" && srun --nodes="$nodes" --ntasks="$mpi" \
                 --ntasks-per-node="$CORES_PER_NODE" --cpus-per-task=1 \
                 --hint=nomultithread --exact \
                 "$CP2K_EXE" -i run.inp >"$out" 2>&1 )
            grep -q "PROGRAM ENDED" "$out" 2>/dev/null || continue
            md=$(awk '/^ qs_mol_dyn_low/ {print $(NF-1)}' "$out" | tail -1)
            tps=$(awk -v t="${md:-0}" -v s="$STEPS" 'BEGIN{printf "%.6f", t/s}')
            [[ $r -eq 0 ]] && continue
            tps_list+="$tps "
         done
         read -r TM TS TMIN NOK <<<"$(stats $tps_list)"
         if [[ "${NOK:-0}" == "0" ]]; then
            printf "%-22s N=%-5s nodes=%s : *** ALL REPS FAILED ***\n" "$LABEL" "$size" "$nodes"
            echo "$size,$nodes,$mpi,0,NA,NA,NA,NA,NA" >>"$CSV_FILE"
            continue
         fi
         [[ -z "$BASE_TPS" ]] && BASE_TPS=$TM
         spd=$(awk -v b="$BASE_TPS" -v t="$TM" 'BEGIN{printf "%.3f", b/t}')
         eff=$(awk -v sp="$spd" -v n="$nodes" 'BEGIN{printf "%.1f", 100*sp/n}')
         printf "%-22s N=%-5s nodes=%s (%3s ranks) : t/step=%s+/-%s s  speedup=%sx  eff=%s%%\n" \
            "$LABEL" "$size" "$nodes" "$mpi" "$TM" "$TS" "$spd" "$eff"
         echo "$size,$nodes,$mpi,$NOK,$TM,$TS,$TMIN,$spd,$eff" >>"$CSV_FILE"
      done
   done
   echo "Wrote $CSV_FILE"
done

rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs \
      "/rds/user/$USER/hpc-work/cp2k-benchmarks/results/" "$BENCHMARK_ROOT/results/"
echo "CSVs synced to $BENCHMARK_ROOT/results"
