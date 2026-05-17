#!/usr/bin/env bash
# Single-pass Morawietz Fig. S4 replication pipeline for ONE system size.
#
# Stages (chained via CP2K EXT_RESTART):
#   1. GEO_OPT  — BFGS, up to 500 iter, drains the ~9 kcal/mol/H2O excess PE
#                 that we identified in the FD diagnostic.
#   2. NVT      — CSVR, T=300K, TIMECON=100 fs, 30 ps (60000 × 0.5 fs steps).
#   3. NVE      — no thermostat, 30 ps, positions every 10 steps for MSD.
#   4. MSD/D    — compute_msd.py extracts D_PBC.
#
# Called as:  ./replicate_figS4_pipeline.sh  <N_molecules>  <mx> <my> <mz>  <cp2k_exe>  <rundir>
#
# This script is INVOKED inside a SLURM job (see replicate_figS4_N*.sh) which
# provides the appropriate ntasks / wall time per size.

set -euo pipefail

N_H2O=$1
MX=$2; MY=$3; MZ=$4
CP2K_EXE=$5
RUNDIR=$6
NTASKS=${SLURM_NTASKS:-16}

# Production trajectory parameters (kept short to fit single-pass overnight)
GEO_MAX_ITER=500
NVT_STEPS=60000     # 30 ps at 0.5 fs
NVE_STEPS=60000     # 30 ps at 0.5 fs
TRAJ_EACH=10        # write positions every 10 steps -> 6000 frames

BENCH=/home/crm98/cp2k-benchmarks
mkdir -p "$RUNDIR"
ln -sfn "$BENCH/potentials" "$RUNDIR/NNP"
cd "$RUNDIR"

# ---------------------------------------------------------------------------
# Generate the three stage inputs from the N=64 template, swapping
# MULTIPLE_UNIT_CELL, RUN_TYPE, ENSEMBLE, and EXT_RESTART blocks as needed.
# ---------------------------------------------------------------------------
python3 <<PYEOF
import re

base = open("$BENCH/H2O-64_RPBE-vdW_NNP.inp").read()

# Common: add MULTIPLE_UNIT_CELL into &SUBSYS &CELL and &TOPOLOGY (need both for tiling)
def with_mult(src, n_h2o):
    if "$MX $MY $MZ" == "1 1 1":
        return src   # no replication needed
    src = re.sub(r'^([ \t]*&CELL[ \t]*\n)',
                 rf'\g<1>      MULTIPLE_UNIT_CELL $MX $MY $MZ\n', src, count=1, flags=re.M)
    src = re.sub(r'^([ \t]*&SUBSYS[ \t]*\n)',
                 rf'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL $MX $MY $MZ\n    &END TOPOLOGY\n',
                 src, count=1, flags=re.M)
    return src

# Stage 1: GEO_OPT
g = with_mult(base, $N_H2O)
g = re.sub(r'PROJECT\s+\S+',  f'PROJECT  relax_N$N_H2O', g)
g = re.sub(r'RUN_TYPE\s+\S+', f'RUN_TYPE GEO_OPT', g)
g = re.sub(r'&MD\b.*?&END MD\s*\n',
           f'&GEO_OPT\n'
           f'    OPTIMIZER BFGS\n'
           f'    MAX_ITER $GEO_MAX_ITER\n'
           f'    MAX_FORCE 1.0E-3\n'
           f'  &END GEO_OPT\n',
           g, flags=re.S, count=1)
g = re.sub(r'\s*&THERMOSTAT.*?&END THERMOSTAT\s*\n', '\n', g, flags=re.S)
open("01_geoopt.inp", "w").write(g)

# Stage 2: NVT, restart positions from GEO_OPT
nvt = with_mult(base, $N_H2O)
nvt = re.sub(r'PROJECT\s+\S+',  f'PROJECT  equil_N$N_H2O', nvt)
nvt = re.sub(r'RUN_TYPE\s+\S+', f'RUN_TYPE MD', nvt)
nvt = re.sub(r'STEPS\s+\d+',    f'STEPS $NVT_STEPS', nvt, count=1)
nvt += ("\n&EXT_RESTART\n"
        "  RESTART_FILE_NAME relax_N$N_H2O-1.restart\n"
        "  RESTART_DEFAULT  .FALSE.\n"
        "  RESTART_POS      .TRUE.\n"
        "  RESTART_COUNTERS .FALSE.\n"
        "&END EXT_RESTART\n")
open("02_nvt.inp", "w").write(nvt)

# Stage 3: NVE production, restart both pos and vel from NVT
nve = with_mult(base, $N_H2O)
nve = re.sub(r'PROJECT\s+\S+',  f'PROJECT  prod_N$N_H2O', nve)
nve = re.sub(r'RUN_TYPE\s+\S+', f'RUN_TYPE MD', nve)
nve = re.sub(r'ENSEMBLE\s+NVT', f'ENSEMBLE NVE', nve)
nve = re.sub(r'STEPS\s+\d+',    f'STEPS $NVE_STEPS', nve, count=1)
nve = re.sub(r'\s*&THERMOSTAT.*?&END THERMOSTAT\s*\n', '\n', nve, flags=re.S)
# Enable trajectory dump every 10 steps
nve = re.sub(r'&TRAJECTORY\s+OFF', '&TRAJECTORY ON', nve, count=1)
nve = re.sub(r'(&TRAJECTORY ON.*?MD\s+)\d+', rf'\g<1>$TRAJ_EACH', nve, count=1, flags=re.S)
nve += ("\n&EXT_RESTART\n"
        "  RESTART_FILE_NAME equil_N$N_H2O-1.restart\n"
        "  RESTART_DEFAULT  .FALSE.\n"
        "  RESTART_POS      .TRUE.\n"
        "  RESTART_VEL      .TRUE.\n"
        "  RESTART_COUNTERS .FALSE.\n"
        "&END EXT_RESTART\n")
open("03_nve.inp", "w").write(nve)

print(f"wrote 01_geoopt.inp 02_nvt.inp 03_nve.inp  (N=$N_H2O, MULT=$MX×$MY×$MZ)")
PYEOF

# ---------------------------------------------------------------------------
# Stage 1: GEO_OPT
# ---------------------------------------------------------------------------
echo ""
echo "=== [N=$N_H2O] STAGE 1: GEO_OPT (BFGS, max $GEO_MAX_ITER iter) ==="
t0=$SECONDS
srun --ntasks=$NTASKS --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i 01_geoopt.inp > 01_geoopt.out 2>&1 || true
echo "  wall: $(( SECONDS - t0 )) s"
if [[ ! -f relax_N$N_H2O-1.restart ]]; then
    echo "*** GEO_OPT did not produce a restart file — aborting"
    tail -30 01_geoopt.out
    exit 1
fi
grep -E "OPTIMIZATION STEP:|GEOMETRY OPTIMIZATION COMPLETED|Total energy" 01_geoopt.out | tail -5

# ---------------------------------------------------------------------------
# Stage 2: NVT equilibration
# ---------------------------------------------------------------------------
echo ""
echo "=== [N=$N_H2O] STAGE 2: NVT 30 ps (CSVR, T=300K, TIMECON=100 fs) ==="
t0=$SECONDS
srun --ntasks=$NTASKS --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i 02_nvt.inp > 02_nvt.out 2>&1 || true
echo "  wall: $(( SECONDS - t0 )) s"
if [[ ! -f equil_N$N_H2O-1.restart ]]; then
    echo "*** NVT did not produce a restart file — aborting"
    tail -30 02_nvt.out
    exit 1
fi
awk 'NR>1' equil_N$N_H2O-1.ener | awk 'END {printf "  NVT final: step=%s  T=%.1f K  ConsQty=%.5f\n", $1, $4, $6}'

# ---------------------------------------------------------------------------
# Stage 3: NVE production
# ---------------------------------------------------------------------------
echo ""
echo "=== [N=$N_H2O] STAGE 3: NVE 30 ps (positions every $TRAJ_EACH steps) ==="
t0=$SECONDS
srun --ntasks=$NTASKS --cpus-per-task=1 --hint=nomultithread \
     "$CP2K_EXE" -i 03_nve.inp > 03_nve.out 2>&1 || true
echo "  wall: $(( SECONDS - t0 )) s"
TRAJ=$(ls prod_N${N_H2O}-pos-1.xyz 2>/dev/null | head -1)
if [[ -z "$TRAJ" ]]; then
    echo "*** NVE did not produce a positions trajectory — aborting"
    tail -30 03_nve.out
    exit 1
fi
awk 'NR>1' prod_N$N_H2O-1.ener | awk 'END {printf "  NVE final: step=%s  T=%.1f K  ConsQty=%.5f (should be conserved)\n", $1, $4, $6}'

# ---------------------------------------------------------------------------
# Stage 4: MSD → D_PBC
# ---------------------------------------------------------------------------
echo ""
echo "=== [N=$N_H2O] STAGE 4: MSD / D_PBC ==="
# Cell dims for unwrap.  N=64 cubic 12.42; MULT scales each axis.
PY_GET_BOX=$(python3 -c "
n=$N_H2O; mx,my,mz=$MX,$MY,$MZ
Lx,Ly,Lz = 12.42*mx, 12.42*my, 12.42*mz
print(f'{Lx} {Ly} {Lz}')")
read LX LY LZ <<<"$PY_GET_BOX"
# dt_fs = TIMESTEP * TRAJ_EACH = 0.5 * 10 = 5 fs between trajectory frames
python3 "$BENCH/scripts/CSD3_benchmark_scripts/figS4/compute_msd.py" \
    "$TRAJ" 5.0 "$LX" "$LY" "$LZ"

echo ""
echo "rundir: $RUNDIR"
