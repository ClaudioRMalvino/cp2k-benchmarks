# cp2k-benchmarks

Build, benchmark, and analysis harness for the MPhil thesis work on
optimising the CP2K Neural Network Potential (NNP) module.

The repository holds:

- per-cluster build / install / rebuild scripts for the **master** and
  **optimised** CP2K trees on CSD3 (Cambridge HPC) and Cerberus (lab
  workstation),
- SLURM driver scripts for the size, strong-scaling, and OMP scaling
  sweeps,
- MAQAO and `perf` profiling drivers,
- the H₂O-64 NNP input deck plus the bulk-water neural-network
  potentials (Schran et al., JCP 2020) used by every benchmark,
- a single Python entry point that regenerates every figure and
  appendix table in the thesis from the raw CSVs.

---

## Layout

```
cp2k-benchmarks/
├── H2O-64_NNP_MD.inp          headline benchmark input
├── NNPs/, NNPs.tar            bulk-H2O committee NN potentials
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
│   ├── CSD3_benchmark_scripts/
│   │   ├── cp2k_CSD3_env.sh           module + toolchain setup
│   │   ├── rebuild_all_branches.sh    full 3-branch rebuild + cache
│   │   ├── scaling/                   size/core/OMP scaling SLURM drivers
│   │   ├── perf_cache/                perf stat hardware-counter runs
│   │   ├── figS4/                     fig S4 (diffusion / viscosity)
│   │   └── nnp_acs_pub/               RDF reproduction of the cnnp paper
│   ├── cerberus_benchmark_scripts/    size + core scaling on Cerberus
│   └── maqao_scripts/                 MAQAO profile + CQA configs
│
├── results/                   raw CSVs / per-step traces, by tree
│   ├── cp2k_master/
│   ├── cp2k_optimized/
│   ├── cp2k_feature_native_spline/
│   ├── cp2k_feature_native_spline_omp/
│   ├── cp2k_feature_verlet_cells/
│   ├── maqao/                 MAQAO profile reports
│   ├── perf_cache/            perf stat outputs (per timestamp)
│   ├── figS4/                 MSD / viscosity outputs
│   └── paper_fig2/            cnnp-paper RDF replication
│
├── plots/                     analysis + figure generation
│   ├── thesis_figures.py      single entry point — produces all 8 figures
│   ├── maqao_data_analysis.py MAQAO CQA → CSVs + supplementary plots
│   ├── csd3_tables.py         scaling tables
│   ├── compiler_flags_table.py stand-alone copy of fig 8 (now ported)
│   ├── thesis_figures/        rendered PDFs (output)
│   └── scaling_csd3/          scaling-table outputs
│
└── logs/                      SLURM .out logs
```

The CP2K source trees themselves (`cp2k_master/`, `cp2k_optimized/`)
live outside this repo under `~/cp2k_master` and `~/cp2k_optimized`;
only the build glue lives here.

---

## Quick start

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
python3 plots/thesis_figures.py              # all 8 figures
python3 plots/thesis_figures.py --only 7     # just one
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
python3 plots/maqao_data_analysis.py
```

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
  potential at `NNPs/RPBE-vdW/nnp-1` (the first committee member from
  Schran et al., JCP 2020).
- The compiler-flags table (fig 8) captures everything the build
  scripts choose; the trailing flags appended by CMake's `Release`
  preset and by CP2K's own `CMakeLists.txt` are build-system plumbing
  and are not enumerated.
- Every scaling CSV is paired with a `*_raw.csv` per-replicate dump so
  fig 4's box-and-jitter can run from the same source.

---

## Author

Claudio Malvino (crm98@cam.ac.uk).  MPhil, University of Cambridge.
