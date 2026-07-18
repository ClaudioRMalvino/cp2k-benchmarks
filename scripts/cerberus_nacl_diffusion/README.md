# NaCl(aq) diffusion with finite-size (Yeh-Hummer) correction — cerberus

Project: compute diffusion coefficients of a **nonhomogeneous liquid**
(NaCl solution: water + Na+ + Cl-) with the optimized CP2K NNP code, apply
the **Yeh-Hummer finite-size correction**, and produce the same
D-vs-1/L extrapolation graph that the group already made for pure water
(the "Morawietz figS4 replication", `plots/plot_figS4.py`).

## Translation of the supervisor's comments

1. *"Homogeneous liquids isotropic Yeh-Hummer correction"* — for a periodic
   cubic box, the measured diffusion coefficient is suppressed by
   hydrodynamic self-interaction across the boundaries. Yeh & Hummer (2004):
   `D_0 = D_PBC + kB*T*xi / (6*pi*eta*L)`, xi = 2.837297, eta = shear
   viscosity, L = box edge. This was derived for a homogeneous isotropic
   liquid, and the group's water study already validated it
   (`scripts/CSD3_benchmark_scripts/figS4/`).
2. *"replicate diffusion finite correction for NaCl nonhomogeneous liquids"* —
   repeat that finite-size study for the NaCl solution. The YH slope
   `-kB*T*xi/(6*pi*eta)` is species-independent (it only depends on the
   solvent's hydrodynamics), so the test is: do D(1/L) for water, Na+ and
   Cl- all extrapolate consistently with the same YH slope?
3. *"use Dhruv's models to do the NaCl diffusion"* — the committee NNP
   models staged in this repo (`potentials/oneill-nacl-water-2024`, also at
   `/data/cerberus1/crm98/nacl-water`): O'Neill et al. JPCL 2024, 8-member
   C-NNP committees at 6 levels of theory, run through Dhruv's optimized
   CP2K NNP implementation.

## Pipeline (mirrors the water figS4 pipeline)

| stage | script | what it does |
|---|---|---|
| 0 | `make_base_cell.py` | Builds the 12.42 A base cell: 1 Na + 1 Cl + 62 H2O (~0.9 mol/kg) by substituting 2 waters in a training-set water frame (called automatically by stage 1). |
| 1 | `run_equil.sh` | NVT (CSVR) equilibration per box size: 30 ps + 5 snapshots 15 ps apart. Box sizes via MULTIPLE_UNIT_CELL tiling: cell111/211/221/222 = 188/376/752/1504 atoms. |
| 2 | `run_production.sh` | 5 x 100 ps NVE segments per size from the snapshots (pos+vel only). Writes XYZ trajectory every 10 steps and stress every step. |
| 2.5 | `check_run.py` | Physics sanity check on any run dir: temperature, conserved-energy drift (per atom/ps), committee disagreement. Run it on the first finished equil and on every production segment BEFORE aggregating: `pyenv python3 check_run.py .../production/cell*/seg*` |
| 3 | `aggregate_nacl.py` | Green-Kubo viscosity per segment + MSD diffusion per species (O, Na, Cl) + YH correction, averaged over segments. Reuses the figS4 `compute_viscosity.py`/`compute_diffusion.py` unchanged. |
| 4 | `plot_nacl_diffusion.py` | The Morawietz-style D vs 1/L extrapolation figure. |
| 5 | `compute_rdf.py` | O-O / Na-O / Cl-O g(r) from production trajectories, for validation against O'Neill Fig 1: `compute_rdf.py .../cell221/seg*/*-pos-1.xyz --cell-ang 24.84 24.84 12.42 --plot` |

All heavy data goes to `/data/cerberus1/crm98/nacl_diffusion/<model>/`.

## How to run

```bash
# 0) build the optimized binary first (see "CP2K branch" below), then:
cd ~/cp2k-benchmarks/scripts/cerberus_nacl_diffusion

# 1) smoke test one tiny run before committing days of compute:
CELLS="111" ./run_equil.sh          # ~1.2 h at 32 ranks on an idle node

# check the smoke test: stress must be nonzero (viscosity depends on it)
head -3 /data/cerberus1/crm98/nacl_diffusion/revPBE-D3/equil/cell111/*-1.stress

# 2) full equilibration + production (run in a persistent session, e.g. tmux)
./run_equil.sh
./run_production.sh

# 3) analysis + plot (system python has no numpy - use the scratch venv)
PY=/data/cerberus1/crm98/pyenv/bin/python3
$PY aggregate_nacl.py \
  --prod-root /data/cerberus1/crm98/nacl_diffusion/revPBE-D3/np1/production \
  --out-dir   /data/cerberus1/crm98/nacl_diffusion/revPBE-D3/np1/analysis
$PY plot_nacl_diffusion.py \
  --summary /data/cerberus1/crm98/nacl_diffusion/revPBE-D3/np1/analysis/nacl_summary.csv \
  --out nacl_diffusion_replication.png
```

Knobs (env vars): `MODEL` (revPBE-D3 default; MP2, RPA, r2SCAN, revPBE0-D3,
optB88-vdW available), `CELLS`, `SEGMENTS`, `TOTAL_RANKS` (default 32 = physical cores; hyperthreads hurt),
`PROD_PS`, `N_PAIRS` (NaCl pairs per base cell; 0 = pure-water reference).

Run layout: `/data/cerberus1/crm98/nacl_diffusion/<model>/np<N_PAIRS>/{equil,production}/cell<mult>/...`

## SLURM (6 h walltime) — `slurm/`

The run scripts are checkpoint-safe: equilibration writes a restart every
15 ps, production every 5 ps, and both **continue from the last checkpoint**
if rerun after a walltime kill (finished work is skipped). So the SLURM
strategy is simply "one array task per unit of work, resubmit until done":

```bash
sbatch slurm/sbatch_smoketest.sh      # ONCE first: node CPU/mount/binary check
eq=$(sbatch --parsable slurm/sbatch_equil.sh)
sbatch --dependency=afterok:$eq slurm/sbatch_production.sh
# cell222 tasks exceed 6 h -> resubmit the same script; they resume.
```

Cluster facts (from CSC docs, July 2026): submit host is **athena**
(`ssh athena`, home dir is the same); `csc-mphil` partition =
phy-cerberus4/5/6 (48 cores, 248 GB each, **6 h max**); nodes are shared
up to 4 jobs unless you request all 48 cores (the wrappers do).
There is also an `lsc` partition (phy-cerberus7/8, **36 h max**) — if
account crm98 is allowed there, one `--partition=lsc --time=36:00:00`
job per cell/segment removes all resubmission juggling; test with a tiny
job first.

Filesystems on compute nodes: home and `/data/athena` are mounted
everywhere; each node's disk is exported as `/data/<nodename>`. Our data
lives on cerberus1's disk (`/data/cerberus1`) — the smoke test verifies
the compute nodes can see it. NOTE: the binary is compiled `-march=native`
on cerberus1 (Ice Lake); the smoke test also catches CPU incompatibility
(SIGILL) on the allocated node. Alternatively the run scripts work
directly with `mpirun` on any idle cerberus1/2/3 node (no SLURM there).

## Concentration series (D vs molality, paper SI Fig S10 analogue)

At fixed cell size (211 recommended), vary `N_PAIRS`: 0, 1, 2, 3, 4
-> 0, 0.90, 1.85, 2.87, 3.97 mol/kg:

```bash
for np in 0 1 2 3 4; do
  N_PAIRS=$np CELLS=211 ./run_equil.sh && N_PAIRS=$np CELLS=211 ./run_production.sh
done
# aggregate each np<k>/production separately, then:
PY=/data/cerberus1/crm98/pyenv/bin/python3
$PY plot_conc.py --cell 211 0.0=<np0 summary> 0.90=<np1 summary> ...
```

Expected physics: water D *decreases* with NaCl concentration
(experimentally ~linearly; ions bind their hydration shells and raise the
viscosity), and ion D also decreases (Na+/Cl- attract, pair, and drag).
An *increase* would be anomalous — worth double-checking with the
supervisor what trend he expects.

## Hardware / runtime (measured 2026-07-02)

cerberus1/2/3 are identical: 2x Xeon Silver 4314 (Ice Lake), 32 physical
cores, 62 GB RAM, no SLURM (plain mpirun). `/data/cerberus1` and home are
shared across nodes (NFS), so the scripts run unchanged on any node —
**run on an idle node** (`ssh cerberus2`, check `uptime` first; cerberus1
and cerberus3 were fully loaded by other users at time of writing).

Measured with the optimized branch, cell111 (188 atoms, committee of 8 +
Coulomb baseline, 32 ranks): **0.020 s/step**. Assuming ~linear scaling
with atom count:

| cell | atoms | equil (210k steps) | 5x100 ps production |
|---|---|---|---|
| 111 | 188  | ~1.2 h | ~5.5 h |
| 211 | 376  | ~2.4 h | ~11 h |
| 221 | 752  | ~4.7 h | ~22 h |
| 222 | 1504 | ~9.4 h | ~44 h |

Full single-model campaign: ~4-4.5 days on one idle node. Refine after the
first cell211 run using the printed walltimes.

## CP2K branch — IMPORTANT

`~/cp2k` `master` is **vanilla upstream** — it does NOT contain the NNP
optimizations. They live on `origin/pr/nnp-cpu-cell-list` (cell-list
neighbour cache, spline ACSFs, OpenMP, sparse force scatter). To build the
optimized binary:

```bash
cd ~/cp2k && git checkout -b nnp-cpu-cell-list origin/pr/nnp-cpu-cell-list
~/cp2k-benchmarks/cp2k_optimized/cerberus_build_scripts/cp2k_cerberus_fork_build.sh
```

## Decisions taken (confirm with supervisor)

1. **Concentration**: 1 NaCl pair per 62 waters (~0.9 mol/kg). The paper's
   PMF system was 1.0 mol/kg (6 pairs/332 waters); training frames are
   ~2.1 mol/kg (2 pairs/52 waters, `--mode extract`). Change with
   `make_base_cell.py --n-pairs N`.
2. **Timestep/H mass**: 0.5 fs with natural hydrogen, consistent with the
   water figS4 study (and hence directly comparable to it). O'Neill et al.
   used 1.0 fs with H mass 2.0 (deuterated) — toggle in the template.
3. **Ensemble/density**: NVT/NVE at the training density (12.42 A cell),
   like the water study. No NpT re-equilibration per model.
4. **Ion statistics**: 1 pair/cell means ion MSDs average over few atoms
   (1 ion at cell111 up to 8 at cell222) x 5 segments — expect large error
   bars on D_Na/D_Cl at small sizes. Remedies: more segments
   (`SEGMENTS="1 ... 10"` after extending equilibration), longer PROD_PS,
   or more pairs per cell.
5. **Viscosity from the MIXED force eval**: STRESS_TENSOR ANALYTICAL is set
   on all three force evals; verify the smoke-test stress file is nonzero
   before production (if the NNP branch lacks the analytical virial, the
   GK viscosity route fails and eta must come from literature/experiment
   instead).
