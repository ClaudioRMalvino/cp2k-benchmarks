# cp2k-benchmarks

Build, benchmark, and analysis harness for the MPhil thesis work:
**Report 1** (Part I: optimising the CP2K Neural Network Potential
module) and **Report 2** (Part II: NaCl(aq) transport with the
optimised implementation).

The repository holds:

- per-cluster build / install / rebuild scripts for the **master** and
  **optimised** CP2K trees on CSD3 (Cambridge HPC) and Cerberus (lab
  workstation),
- SLURM driver scripts for the size, strong-scaling, and OMP scaling
  sweeps, plus MAQAO and `perf` profiling drivers (Report 1),
- the complete Report 2 simulation and analysis pipelines: the
  Madrid-2019 force-field campaigns (LAMMPS), the MP2 and RPA
  committee-potential runs (CP2K), the on-the-fly DFT cost ladder,
  and the cost/carbon accounting,
- the input decks and the external neural-network potentials used by
  every run (Schran et al. bulk water; O'Neill et al. NaCl(aq)
  committee models),
- the derived data (`results/`) and Python scripts that regenerate
  **every figure and table in both reports** on any machine — no
  cluster access needed.

## Reproducing the results — what an assessor needs to know

Every number and figure in both reports regenerates from data bundled
in this repository. Reproduction has two tiers:

1. **Analysis and figures (any machine, minutes).** `results/` holds
   the derived data of every campaign (scaling CSVs, pooled MSD /
   Onsager / viscosity / RDF outputs, profiling reports). The Python
   scripts below turn them into exactly the figures shipped in the
   reports. Set up once:

   ```bash
   python3 -m venv .venv && source .venv/bin/activate
   pip install numpy matplotlib pandas scipy
   ```

2. **Full re-simulation (HPC, weeks).** The SLURM drivers, input
   decks, build scripts and per-pipeline READMEs regenerate the raw
   data from scratch on a CSD3-like machine. Raw MD trajectories are
   not shipped (terabytes); everything derived from them is.

The two sections below give the per-report commands.

---

## Layout

```
cp2k-benchmarks/
├── H2O-64_NNP_MD.inp          headline benchmark input
├── potentials/                external reference NN potentials (see potentials/README.md)
│   ├── morawietz-water-2016/  bulk-H2O committee NNPs (was NNPs/)
│   ├── RPBE-vdW-2016/         runnable RPBE-vdW potential (with input.nn)
│   └── oneill-nacl-water-2024/ NaCl-in-water committee models (was nacl/)
├── diffusion/                 fig S4 inputs (MSD / diffusion / viscosity)
│
├── cp2k_master/               build scripts for the upstream CP2K tree
│   ├── cerberus_build_scripts/
│   └── CSD3_build_scripts/
├── cp2k_optimized/            build scripts for the optimised tree
│   ├── cerberus_build_scripts/
│   └── CSD3_build_scripts/
│
├── scripts/
│   ├── CSD3_benchmark_scripts/        — Report 1 —
│   │   ├── cp2k_CSD3_env.sh           module + toolchain setup
│   │   ├── rebuild_all_branches.sh    full 3-branch rebuild + cache
│   │   ├── scaling/                   size/core/OMP scaling SLURM drivers
│   │   ├── perf_cache/                perf stat hardware-counter runs
│   │   ├── figS4/                     fig S4 (diffusion / viscosity)
│   │   └── nnp_acs_pub/               RDF reproduction of the cnnp paper
│   ├── cerberus_benchmark_scripts/    size + core scaling on Cerberus
│   ├── maqao_scripts/                 MAQAO profile + CQA configs
│   ├── benchmark_figures/             Report-1 + Chebyshev figure scripts
│   ├── cerberus_nacl_diffusion/       — Report 2 — Madrid-2019 LAMMPS
│   │   └── lammps/                    campaigns 1-3 + make_figures.py (see its README)
│   ├── csd3_nacl_mp2{,_anchor}/       MP2 C-NNP runs + two-cell anchor analysis
│   ├── csd3_rpa_transport/            RPA series runs, Onsager/viscosity/pairing
│   │                                  analysis, make_series_figures.py (see README)
│   └── csd3_dft_cost/                 on-the-fly revPBE-D3 cost ladder (see README)
│
├── results/                   derived data, by campaign
│   ├── cp2k_master/ cp2k_feature_* cp2k_dhruv_*   Report-1 scaling CSVs
│   ├── maqao/                 MAQAO profile reports
│   ├── perf_cache/            perf stat outputs (per timestamp)
│   ├── figS4/                 MSD / viscosity outputs
│   ├── paper_fig2/            cnnp-paper RDF replication
│   ├── cheby_benchmark/       Chebyshev-kernel benchmark tables
│   ├── lammps_madrid/         Report-2 force-field pooled analysis (.npz/.csv)
│   ├── nacl_mp2_anchor/       MP2 two-cell diffusion anchor
│   ├── mp2_transport/ rpa_transport/   C-NNP segment analysis outputs
│   ├── transport_comparison/  cross-level series summaries (CSV)
│   └── dft_cost/              revPBE-D3 timing ladder
│
├── plots/                     rendered figures ONLY (no code, no data)
│   ├── thesis_figures/        thesis figure PDFs (output)
│   ├── scaling_csd3/          rendered scaling-table PDFs
│   └── ...                    cheby_benchmark_figs/, maqao_plots/, etc.
│                              (figure scripts live in scripts/benchmark_figures/,
│                               scripts/maqao_scripts/; derived tables in
│                               results/scaling_csd3/, results/cheby_benchmark/,
│                               results/maqao/)
│
└── logs/                      SLURM .out logs
```

The CP2K source trees themselves (`cp2k_master/`, `cp2k_optimized/`)
live outside this repo under `~/cp2k_master` and `~/cp2k_optimized`;
only the build glue lives here.

---

## Reproducing Report 1 (Part I: the CP2K NNP optimisation)

### 1. Build the binaries (one-time)

**On CSD3:**

```bash
source scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

# Full rebuild + cache of all three branches into /rds/.../cp2k_binaries/csd3/
bash scripts/CSD3_benchmark_scripts/rebuild_all_branches.sh
```

This drives `cp2k_master/CSD3_build_scripts/cp2k_CSD3_master_build.sh`
and `cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_build.sh`, then
checks out `feature/nnp-native-spline` and `feature/nnp-native-spline-omp`
in turn and runs an incremental rebuild against each.  Resulting
binaries land at
`/rds/user/$USER/hpc-work/cp2k_binaries/csd3/<branch>/cp2k.psmp`.

**On Cerberus:**

```bash
bash cp2k_master/cerberus_build_scripts/cp2k_cerberus_master_build.sh
bash cp2k_optimized/cerberus_build_scripts/cp2k_cerberus_opt_build.sh
```

### 2. Run the benchmarks

The four scaling workers in
`scripts/CSD3_benchmark_scripts/scaling/` are the canonical entry
points; each takes a branch name and writes a timestamped CSV under
`results/<branch>/NNP/`:

| Worker | Sweep | Feeds |
|---|---|---|
| `run_nnp_size_scaling_slurm.sh` | N = 64 → 4096 H₂O at full node | fig 1, fig 4 |
| `run_nnp_core_scaling_slurm.sh` | 1 → 76 cores at N = 1024 | fig 2 |
| `run_nnp_omp_thread_scaling_slurm.sh` | OMP = 1 → 76 at N = 64 | fig 3 (a) |
| `run_nnp_omp_size_scaling_slurm.sh` | OMP × N grid | fig 3 (b), fig 6 |

`benchmark_slurm.sh` chains the size + core sweeps across all three
branches for a from-scratch run.  On Cerberus the equivalent driver is
`scripts/cerberus_benchmark_scripts/benchmark_slurm.sh`.

For the hardware-counter cache-miss data feeding fig 7, run the
**login-node** perf driver:

```bash
bash scripts/CSD3_benchmark_scripts/perf_cache/run_perf_cache_login.sh
```

For the MAQAO CQA static-analysis data:

```bash
bash scripts/maqao_scripts/run_maqao_profile_slurm.sh
```

### 3. Generate the thesis figures

```bash          # pandas + scipy
python3 scripts/benchmark_figures/thesis_figures.py              # all 8 figures
python3 scripts/benchmark_figures/thesis_figures.py --only 7     # just one
```

Output PDFs land in `plots/thesis_figures/`:

| File | Content |
|---|---|
| `fig1_algorithmic_complexity.pdf` | per-step time vs N (linear + log + per-atom + speedup) |
| `fig2_strong_scaling.pdf` | strong-scaling sweep at N = 1024 |
| `fig3_openmp_threading.pdf` | OMP thread sweep + OMP×N grid |
| `fig4_statistical_significance.pdf` | per-replicate distribution at N = 1024 |
| `fig5_maqao_microarchitectural.pdf` | MAQAO hotspot + array-access efficiency |
| `fig6_omp_size_heatmap.pdf` | OMP×N heatmap preview |
| `fig7_perf_summary.pdf` | wall, instructions, IPC, DRAM traffic |
| `fig8_compiler_flags.pdf` | appendix compiler-flag table |

The MAQAO supplementary plots and the `maqao_*_table.csv` inputs used
by fig 5 are produced by:

```bash
python3 scripts/maqao_scripts/maqao_data_analysis.py
```

---

## Reproducing Report 2 (Part II: NaCl(aq) transport)

### Figures and numbers from the bundled data (any machine)

Three scripts regenerate every Report 2 figure from `results/`:

```bash
python3 scripts/cerberus_nacl_diffusion/lammps/make_figures.py      # -> plots/lammps_madrid/
python3 scripts/csd3_nacl_mp2_anchor/plot_nacl_diffusion.py         # -> plots/nacl_mp2_anchor/
python3 scripts/csd3_rpa_transport/make_series_figures.py           # -> plots/transport_comparison/
python3 scripts/benchmark_figures/plot_cheby_full_suite.py          # -> plots/cheby_benchmark_figs/
python3 scripts/benchmark_figures/plot_cheby_hybrid_mem.py          # -> plots/cheby_benchmark_figs/
```

Mapping to the figure numbers in Report 2:

| Report figure | File | Script |
|---|---|---|
| 8.1 | MD snapshots (`nacl_*molkg_md.png`) | rendered by hand in OVITO (only unscripted figure) |
| 8.2 | `fig1_msd_vs_vacf.pdf` | `make_figures.py` |
| 8.3 | `fig2_yeh_hummer.pdf` | `make_figures.py` |
| 8.4 | `fig35_kappa_tna.pdf` | `make_figures.py` |
| 8.5 | `nacl_diffusion_yh.pdf` | `plot_nacl_diffusion.py` |
| 9.1 | `fig4_ne_deviation.pdf` | `make_figures.py` |
| 9.2 | `fig67_decomposition.pdf` | `make_figures.py` |
| 9.3 | `fig8_three_way_1m.pdf` | `make_series_figures.py` |
| 9.4–9.8 | `fig9`–`fig13` series PDFs | `make_series_figures.py` |
| 9.9 | `fig15_carbon.pdf` | `make_series_figures.py` |
| A.1 | `fig5_size_scaling.pdf` | `plot_cheby_full_suite.py` |
| A.2 | `fig2_speed_vs_memory_N1024.pdf` | `plot_cheby_hybrid_mem.py` |
| A.3 | `fig3_cheby_omp_scaling.pdf` | `plot_cheby_hybrid_mem.py` |

`make_series_figures.py` also rewrites the pooled series tables
(`results/transport_comparison/transport_series_summary.csv` and
`pairing_summary_cube3.csv`) that back the numbers quoted in
Chapters 8–9, and the cost/carbon numbers of §9.4 come from
`scripts/csd3_rpa_transport/energy_cost_ledger.py` and
`green_algorithms_footprint.py`.

### Full re-simulation (HPC)

Each pipeline has a README written for exactly this purpose:

- **Madrid-2019 force field (LAMMPS)** —
  `scripts/cerberus_nacl_diffusion/lammps/README.md`. Three
  campaigns: (1) self-diffusion + Yeh–Hummer boxes, (2) the 50-run
  conductivity/Onsager series, (3) NVE velocity dumps for the
  Green–Kubo cross-checks; each is one sbatch driver plus one
  analysis script.
- **MP2 anchor (CP2K C-NNP)** — `scripts/csd3_nacl_mp2_anchor/`:
  `run_equil.sh` / `run_production.sh` submit the two cells
  (24.8 / 37.3 Å) on CSD3, `analyze_diffusion.py` produces
  `results/nacl_mp2_anchor/`.
- **RPA concentration series (CP2K C-NNP)** —
  `scripts/csd3_rpa_transport/README.md`: the 1/2/4 mol kg⁻¹
  segments and the `analyze_{onsager,viscosity,ion_pairing}_cp2k.py`
  chain feeding `results/{mp2,rpa}_transport/`.
- **On-the-fly DFT cost ladder (revPBE-D3)** —
  `scripts/csd3_dft_cost/README.md`: the four-cell timing ladder
  behind Table 9.4, analysed by `analyze_dft_ladder.py` into
  `results/dft_cost/`.

The committee models themselves are the published O'Neill et al.
parameters, redistributed unchanged under CC-BY-SA-4.0 in
`potentials/oneill-nacl-water-2024/` (see its README and LICENSE).

---

## Configuration

### CSD3 environment

Every CSD3 script sources `scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh`,
which loads the Intel oneAPI 2022 stack (`mpiifort` / `mpiicc` /
`mpiicpc`), MKL, Intel MPI, GCC 11 (for the C++ link), and Python.
The MKL module is mandatory — `cp2k.psmp` dynamically loads
`libmkl_intel_thread.so.2` at startup.

### Compiler flags

Both build scripts pass an **identical** flag set:

```
-O2 -g -xCORE-AVX512 -qopenmp                          # C, C++
-O2 -g -xCORE-AVX512 -qopenmp -funroll-loops -ftree-vectorize   # Fortran
```

See fig 8 (`fig8_compiler_flags.pdf`) for the full table of effective
flags and the explanation of the trailing `-O2` that lowers the
optimised branch's Fortran level to `-O2` (versus master at `-O3`).

### Cluster targets

- **CSD3 Peta4-IceLake**: Intel Xeon Platinum 8368Q, 2 × 38 cores (76
  total), 256 GB DDR4.  Used for every scaling, perf, and MAQAO run in
  the thesis.
- **Cerberus**: lab workstation; only used early in the project for
  initial development sweeps.  The Cerberus scripts are retained but
  not part of the thesis dataset.

---

## Provenance and reproducibility

- The H₂O-64 input deck (`H2O-64_NNP_MD.inp`) references the cnnp
  potential at `NNP/bulkH2O-jcp2020-cnnp/nnp-1` (the first committee
  member from Schran et al., JCP 2020), staged into run directories at
  runtime. The Morawietz reference NNPs live at
  `potentials/morawietz-water-2016/`.
- The compiler-flags table (fig 8) captures everything the build
  scripts choose; the trailing flags appended by CMake's `Release`
  preset and by CP2K's own `CMakeLists.txt` are build-system plumbing
  and are not enumerated.
- Every scaling CSV is paired with a `*_raw.csv` per-replicate dump so
  fig 4's box-and-jitter can run from the same source.

---

## Author

Claudio Malvino (crm98@cam.ac.uk).  MPhil, University of Cambridge.
