# Madrid-2019 NaCl(aq) transport pipeline (LAMMPS, cerberus/CSC)

All simulation data lands under `/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/`
(home has a 5 GB quota — never write trajectories there). Submit sbatch scripts
from athena (the only host with a SLURM client); logs go to `~/cp2k-benchmarks/logs/`.

## Three campaigns — what to run, in order

### 1. Self-diffusion + Yeh–Hummer finite-size (replicas)
```
sbatch sbatch_lammps_replicas.sh        # 4 box sizes x 5 seeds, 1 m NaCl
python3 analyze_replicas.py             # -> replicas/analysis/replica_D.npz
```
Chain: `sbatch_lammps_replicas.sh` → `run_lammps_replica.sh` →
`make_lammps_data.py` (builds the box) + `in.madrid_nacl` (NPT→NVT + MSD production)
→ `msd_multiorigin.py` (multi-origin Einstein D per species).
`analyze_replicas.py` pools seeds and fits D(1/L) (Yeh–Hummer extrapolation).

### 2. Conductivity / Onsager coefficients (concentration series)
```
sbatch sbatch_lammps_conduct.sh         # 5 molalities x 2 boxes x 5 seeds = 50 runs
python3 kappa_analysis.py               # -> conductivity/analysis/kappa_vs_c.npz
```
Chain: `sbatch_lammps_conduct.sh` → `run_lammps_conduct.sh` →
`make_lammps_data.py` + `in.madrid_conduct` → `onsager_multiorigin.py`
(collective displacement correlations C_ij). `kappa_analysis.py` computes
κ_Onsager (z=1 and z=0.85), κ_Nernst–Einstein, and t_Na with seed error bars.

### 3. VACF cross-check (Green–Kubo vs Einstein)
```
sbatch sbatch_lammps_vacf2.sh           # restarts from replica final.data, NVE
python3 msd_vs_vacf.py                  # agreement table (ratio + pull sigma)
```
Chain: `sbatch_lammps_vacf2.sh` → `in.madrid_vacf` (20 ps NVT retherm → 40 ps NVE,
velocity dump every 4 fs) → `vacf_analyze.py` (FFT VACF + running integral D(t);
the dump is deleted after analysis). `msd_vs_vacf.py` compares the VACF plateau D
against the MSD D from campaign 1 — the two routes are mathematically equivalent,
so agreement validates both pipelines.

### 4. Paper figures
```
python3 make_figures.py                 # -> <data>/figures/fig1..fig6 (.png + .pdf)
python3 pair_association.py             # K_A / paired fraction from the 1 m PMF
```
Reads only the analysis products of campaigns 1–3 plus the conductivity runs'
`.rdf` files (Na–Cl PMF). Experimental κ anchors (m = 1, 2, 4 at 298.15 K) from
the compilation in Blazquez et al., J. Chem. Theory Comput. 2023.

fig6 decomposes κ_Ons into κ_NaNa + κ_ClCl − 2 κ_NaCl (per-run slopes saved in
kappa_vs_c.npz by kappa_analysis.py) and splits the diagonal terms into
self + distinct: the Na–Cl cross term is zero within noise at all molalities,
so the NE deviation comes from like-ion (mostly Na–Na) anticorrelation.
`pair_association.py` shows this is quantitatively consistent with the
under-paired PMF: ~2.5% contact pairs at 1 m → predicted pairing contribution
to Δ_NE of a few %, below our seed noise.

## Shared inputs
- `madrid.settings` — Madrid-2019 force field: pair/bond/angle coeffs incl.
  fitted cross terms, SHAKE, groups. Included by every `in.*` file.
- `make_lammps_data.py` — builds an NaCl(aq) box at a given molality and L.
- `in.madrid_nacl` / `in.madrid_conduct` / `in.madrid_vacf` — the three protocols.

## Gotchas an assessor should know
- `in.madrid_vacf` must declare bond/angle/pair styles *before* `read_data` and
  re-include `madrid.settings` *after*: `write_data` stores only i–i pair coeffs,
  so the fitted Madrid cross terms would otherwise be silently lost.
- VACF production is NVE on purpose — a thermostat distorts exactly the
  short-time velocity correlations the Green–Kubo integral measures.
- `vacf_TAINTED_underdense/` and `finite_size/` in the data tree are superseded
  (an early sweep ran ~2.5% under-dense); nothing in this directory reads them.
