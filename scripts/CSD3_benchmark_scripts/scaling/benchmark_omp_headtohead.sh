#!/usr/bin/env bash
# HEAD-TO-HEAD: Dhruv's pristine PR head (8be1dfa50f) vs the SAME head with
# the bit-exact atom-level OpenMP graft. Both binaries come from the SAME
# configured CMake build tree (identical flags -O3 -xCORE-AVX2
# -fp-model=precise -qopenmp, identical toolchain) and run on the SAME node
# in the SAME job, so the ONLY variable is the graft.
#
#   Leg 1 (accuracy): step-0/1 NNP forces on the N=1024 benchmark system,
#          baseline serial vs graft OMP=1 and vs graft OMP=76.
#          Claim to support: bit-identical (max abs diff 0.0).
#   Leg 2: OMP thread scaling, 1 MPI rank, threads {1..76}, N=1024, both.
#   Leg 3: pure-MPI strong scaling, OMP=1, cores {1..76}, N=1024, both.
#   100 MD steps, 5 timed reps + 1 warm-up, t/step from qs_mol_dyn_low.
#
#SBATCH -J NNP_omp_h2h
#SBATCH -A NIKIFORAKIS-CSC-FUNDS-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=76
#SBATCH --time=10:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/home/crm98/cp2k-benchmarks/logs/NNP_omp_h2h_%j.out

set -uo pipefail
. /etc/profile.d/modules.sh
set +u
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
source /rds/user/$USER/hpc-work/cp2k_binaries/csd3/setup
set -u

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
SCALING=/home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling
PROOT=/rds/user/$USER/hpc-work/cp2k_dhruv_omp_worktree
OPARENT=cp2k_omp_headtohead
LIST_R1="1 2 4 8 16 19 32 38 76"       # report-1 NUMA-aware set

LABEL_BASE=dhruv-head-8be1dfa          # pristine PR head
LABEL_OMP=dhruv-omp-bitexact-v2        # head + graft
EXE_BASE="$BIN_ROOT/$LABEL_BASE/cp2k.psmp"
EXE_OMP="$BIN_ROOT/$LABEL_OMP/cp2k.psmp"

echo "=== BINARY PROVENANCE ==="
for e in "$EXE_BASE" "$EXE_OMP"; do
   [[ -x "$e" ]] || { echo "!! missing $e"; exit 1; }
   md5sum "$e"
done
export LD_LIBRARY_PATH="$BIN_ROOT/$LABEL_OMP/lib:${LD_LIBRARY_PATH:-}"

############################################################################
echo ""
echo "######## LEG 1: FORCE ACCURACY, N=1024 (3072 atoms), 1 MD step ########"
ACC=/rds/user/$USER/hpc-work/cp2k-benchmarks/results/$OPARENT/accuracy_${SLURM_JOB_ID}
mkdir -p "$ACC"

prep_input() {  # $1 = rundir
   python3 - "/home/crm98/cp2k-benchmarks/H2O-64_NNP_MD.inp" "$1/run.inp" <<'PY'
import sys, re
src, dst = sys.argv[1], sys.argv[2]
t = open(src).read()
t = re.sub(r'STEPS\s+\d+', 'STEPS 1', t, count=1)
t = re.sub(r'&TRAJECTORY\b[\s\S]*?&END TRAJECTORY',
           '&TRAJECTORY OFF\n    &END TRAJECTORY\n'
           '    &RESTART OFF\n    &END RESTART\n'
           '    &RESTART_HISTORY OFF\n    &END RESTART_HISTORY', t, count=1)
# keep &FORCES ON (that is what we compare); N=1024 = 4x2x2 of the 64 cell
t = re.sub(r'(&CELL\s*\n)', r'\g<1>      MULTIPLE_UNIT_CELL 4 2 2\n', t, flags=re.I)
if "&TOPOLOGY" in t:
    t = re.sub(r'(&TOPOLOGY\s*\n)', r'\g<1>      MULTIPLE_UNIT_CELL 4 2 2\n', t, flags=re.I)
else:
    t = re.sub(r'(&SUBSYS\s*\n)',
               r'\g<1>    &TOPOLOGY\n      MULTIPLE_UNIT_CELL 4 2 2\n    &END TOPOLOGY\n',
               t, flags=re.I)
open(dst, 'w').write(t)
PY
}

acc_run() {  # $1 = name, $2 = exe, $3 = omp
   local d="$ACC/$1"; mkdir -p "$d"
   prep_input "$d"
   ln -sfn "$PROOT/data/NNP" "$d/NNP"
   ( cd "$d" && OMP_NUM_THREADS=$3 srun --ntasks=1 --cpus-per-task=$3 \
        --hint=nomultithread "$2" -i run.inp > run.out 2>&1 ) || true
   grep -q "PROGRAM ENDED" "$d/run.out" || echo "!! accuracy run $1 FAILED"
}

acc_run base_omp1  "$EXE_BASE" 1
acc_run graft_omp1 "$EXE_OMP"  1
acc_run graft_omp76 "$EXE_OMP" 76

python3 - "$ACC" <<'PY'
import sys, glob, os
acc = sys.argv[1]
def forces(name):
    d = os.path.join(acc, name)
    f = sorted(glob.glob(os.path.join(d, "*forces*.xyz")) +
               glob.glob(os.path.join(d, "*-frc-*.xyz")))
    if not f:
        print(f"!! no forces file in {d}"); sys.exit(1)
    frames, cur = [], None
    for ln in open(f[0]):
        p = ln.split()
        if len(p) == 1 and p[0].isdigit():
            cur = []; frames.append(cur); continue
        if cur is not None and len(p) >= 4:
            try: cur.append(tuple(float(x) for x in p[-3:]))
            except ValueError: pass
    return frames, os.path.basename(f[0])
def energy(name):
    for ln in open(os.path.join(acc, name, "run.out")):
        if "INITIAL POTENTIAL ENERGY" in ln or "POTENTIAL ENERGY" in ln.upper():
            try: return float(ln.split()[-1])
            except ValueError: pass
    return None
ref_frames, src = forces("base_omp1")
print(f"reference: base_omp1 {src}, {len(ref_frames)} frame(s), {len(ref_frames[0])} atoms")
for name in ("graft_omp1", "graft_omp76"):
    fr, _ = forces(name)
    mabs = 0.0
    for fa, fb in zip(ref_frames, fr):
        assert len(fa) == len(fb) > 0, f"atom-count mismatch vs {name}"
        for va, vb in zip(fa, fb):
            for u, v in zip(va, vb):
                mabs = max(mabs, abs(u - v))
    print(f"FORCES base_omp1 vs {name:11s}: max_abs={mabs:.3e}  bit-identical? {mabs==0.0}")
ea = energy("base_omp1")
for name in ("graft_omp1", "graft_omp76"):
    eb = energy(name)
    if ea is not None and eb is not None:
        print(f"ENERGY base_omp1 vs {name:11s}: diff={abs(ea-eb):.3e} Ha")
PY

############################################################################
echo ""
echo "########## LEG 2: OMP THREAD SCALING (1 MPI rank, N=1024) ##########"
cd "$SCALING"
for L in "$LABEL_BASE" "$LABEL_OMP"; do
   echo ""; echo "--- $L ---"
   TARGET_LABEL="$L" \
   PROJECT_ROOT_OVERRIDE="$PROOT" \
   OUTDIR_PARENT_OVERRIDE="$OPARENT" \
   MPI_RANKS=1 N_MOLECULES=1024 STEPS=100 N_REPS=5 OMP_LIST="$LIST_R1" \
      ./run_nnp_omp_thread_scaling_slurm.sh
done

############################################################################
echo ""
echo "########## LEG 3: PURE-MPI STRONG SCALING (OMP=1, N=1024) ##########"
for L in "$LABEL_BASE" "$LABEL_OMP"; do
   echo ""; echo "--- $L ---"
   CP2K_EXE_OVERRIDE="$BIN_ROOT/$L/cp2k.psmp" \
   INSTALL_LIB_OVERRIDE="$BIN_ROOT/$L/lib" \
   LABEL_OVERRIDE="$L" \
   PROJECT_ROOT_OVERRIDE="$PROOT" \
   OUTDIR_PARENT_OVERRIDE="$OPARENT" \
   N_MOLECULES=1024 STEPS=100 N_REPS=5 CORE_LIST="$LIST_R1" \
      ./run_nnp_core_scaling_slurm.sh "$L"
done

# Mirror CSVs (only) from scratch into the in-repo results tree.
SCRATCH=/rds/user/$USER/hpc-work/cp2k-benchmarks/results
HOME_R=/home/crm98/cp2k-benchmarks/results
rsync -a --include='*/' --include='*.csv' --exclude='*' --prune-empty-dirs "$SCRATCH/" "$HOME_R/"
echo ""; echo "=== done; head-to-head CSVs synced to $HOME_R/$OPARENT ==="
