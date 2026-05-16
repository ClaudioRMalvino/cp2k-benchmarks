# Methods, configuration, and provenance

This document records every parameter required for an independent scientist
to interpret (and as far as practical, reproduce) the CP2K-NNP performance
benchmarks reported in this thesis.  It is written following the experimental-
design rules of Hoefler & Belli, *"Scientific Benchmarking of Parallel
Computing Systems"*, SC'15 (DOI: 10.1145/2807591.2807644).

---

## 1. Hardware

All performance measurements were obtained on the **Peta4-IceLake**
partition of the University of Cambridge **CSD3** cluster.  Each compute
node has the following specification (verified via `lscpu` on `cpu-q-154` /
`login-q-1`, identical generation):

| Component | Specification |
|---|---|
| CPU model | 2 × Intel Xeon Platinum **8368Q**, "Ice Lake-SP" |
| Sockets per node | 2 |
| Cores per socket | 38 |
| Total physical cores per node | **76** |
| Hardware threads per core (SMT) | **1** (Hyper-threading disabled) |
| NUMA domains per node | 2 (one per socket) |
| Base clock | 2.60 GHz |
| Maximum turbo clock | 3.70 GHz |
| L1d cache (per core) | 48 KiB |
| L1i cache (per core) | 32 KiB |
| L2 cache (per core) | 1.25 MiB |
| L3 cache (shared, per socket) | ~57 MiB (LLC) |
| Vector ISA | AVX, AVX2, **AVX-512** (F/CD/DQ/BW/VL) |
| RAM per node | 256 GiB DDR4-3200 |
| Inter-node interconnect | Mellanox HDR-100 InfiniBand |

All benchmark jobs used `--nodes=1` (a single Ice Lake node), so the
inter-node interconnect does not affect the reported timings.

## 2. Software environment

| Component | Version |
|---|---|
| Operating system | Rocky Linux 8 |
| Kernel `perf_event_paranoid` | **2** (blocks unprivileged PMU access) |
| Kernel `yama/ptrace_scope` | **2** (since 2026-05-15 — Red Hat / Rocky Linux zero-day mitigation; blocks unprivileged `ptrace(PTRACE_ATTACH)`, the syscall MAQAO's LProf sampler uses) |
| Scheduler | SLURM 14.03+ |
| MPI | Intel MPI 2021.6.0 (Spack hash `guxuvcpm`) |
| BLAS / LAPACK / FFTW | Intel MKL 2022.1.0 (Spack hash `qqwrrcxw`) |
| Fortran compiler | Intel ifort **2021.6.0** (Build 20220226) — Intel oneAPI Classic |
| C / C++ compiler | Intel icc / icpc 2021.6.0 |
| Auxiliary C++ runtime | GCC 11.2.0 (`-cxxlib` link with explicit `-L${GCC11_LIB}`) |
| Python (for plotting / RDF) | 3.11.0 (CSD3 module `python/3.11.0-icl`), scipy 1.17.1, pandas 3.0.3, matplotlib 3.10.9, numpy 2.4.4 |
| Valgrind (for cachegrind backup path) | 3.22.0 (system `/usr/bin/valgrind`) |

## 3. CP2K branches under test

| Branch label | Repository | Git commit (HEAD at measurement time) | Description |
|---|---|---|---|
| upstream master | `~/cp2k_master` | `757bb76a80` ("Cut release version 2026.1") | CP2K 2026.1 reference release |
| `feature/nnp-native-spline` | `~/cp2k_optimized` | `7cbc8b3008` ("nnp: port cache-conscious optimisations from feature/nnp-native-spline-omp") | Pure-MPI optimised fork: cell-list neighbour search, persistent scratch buffers, Horner-form Hermite splines (via CP2K-native `spline_data_type`), precomputed ACSF scaling factors, sorted `active_atoms` via in-place heap sort |
| `feature/nnp-native-spline-omp` | `~/cp2k_optimized` | `57b54b31b0` | Hybrid MPI + OpenMP fork: all of the above + thread-private local-reduction OpenMP threading in `nnp_calc_acsf` and `nnp_scale_acsf` |

## 4. Build configuration

All three branches are built from clean (`rm -rf build`) via CMake.  Build
scripts: `cp2k-benchmarks/cp2k_master/CSD3_build_scripts/cp2k_CSD3_master_build.sh`
and `cp2k-benchmarks/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_build.sh`.

### Compiler flags actually applied at compile time

```
CMAKE_BUILD_TYPE      = Release   (CMake appends -O3 -DNDEBUG for C/CXX, -O3 for Fortran)
CMAKE_Fortran_FLAGS   = -O2 -g -xCORE-AVX512 -qopenmp -funroll-loops -ftree-vectorize
CMAKE_C_FLAGS         = -O2 -g -xCORE-AVX512 -qopenmp
CMAKE_CXX_FLAGS       = -O2 -g -xCORE-AVX512 -qopenmp
```

Final effective Fortran flag chain visible in the compile commands:
`-O2 -g -xCORE-AVX512 -qopenmp -funroll-loops -ftree-vectorize -O3 ...`

For the optimized branches *only*, the cp2k_optimized CMakeLists injects an
additional `-xHOST` per-target — this overrides our `-xCORE-AVX512` (the
compiler reports `warning #10121: overriding '-xCORE-AVX512' with '-xHOST'`).
Because the build hosts (Ice Lake login + compute) are themselves Ice Lake,
both flags resolve to the same instruction set and the resulting binaries are
functionally equivalent for AVX-512 vectorisation purposes.

### Linker

```
mpiifort <FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>
   -nofor-main -cxxlib -L${GCC11_LIB} -lstdc++
```

`${GCC11_LIB}` resolves to `/usr/local/software/master/gcc/11/lib64`
(verified by `dirname $(gcc -print-file-name=libstdc++.so)` after
`module load gcc/11`).  This explicit linkage is required to satisfy
`std::__throw_bad_array_new_length()` (`GLIBCXX_3.4.29`) calls in the Plumed
2.9.3 and COSMA 2.7.0 toolchain libraries.

### Binary verification

| Branch | `cp2k.psmp` size | md5sum | `.debug_info` present? | Branch-specific symbols verified |
|---|---:|---|---|---|
| master | 902 994 952 B | `17e0035070eadd4f91319c2d518158a1` | yes | none of `hermite_spline_value`, `nnp_sort_active_atoms`, `nnp_neighbor_ensure_persistent` |
| native-spline | 903 447 856 B | `26271f3c0820d30ea11d20bbe4c967b6` | yes | all three present |
| native-spline-omp | 903 492 744 B | `524e175955a9601d75b3469d9e6dd6e4` | yes | all three present |

Three distinct md5sums confirm that each measurement targets a genuinely
different binary; symbol presence confirms branch identity.

## 5. Runtime configuration

Every cp2k.psmp invocation is launched as:

```
srun --ntasks=<N>  --cpus-per-task=<M>  --hint=nomultithread  \
     <cp2k.psmp> -i run.inp
```

- `--hint=nomultithread` forces one MPI rank per physical core (Hyper-threading
  disabled, see Section 1)
- For pure-MPI branches: `N = total_cores`, `M = 1`
- For native-spline-omp: `N = total_cores / 2`, `M = 2`, `OMP_NUM_THREADS=2`

SLURM submission account: `MPHIL-NIKIFORAKIS-CRM98-SL2-CPU`.
Module environment (loaded by `cp2k_CSD3_env.sh` before every job):

```
module load rhel8/default-icl
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw
module load gcc/11
module load python/3.11.0-icl
source ~/cp2k_master/tools/toolchain/install/setup
```

Per-branch `LD_LIBRARY_PATH` is prepended with each branch's cached
`lib/` directory before that branch's `srun`.

## 6. Workload specification

All MD inputs are derived from `H2O-64_NNP_MD.inp` (METHOD NNP, bulk water).
Per-step trajectory / energy / force I/O is switched OFF at the input level
(by setting `&TRAJECTORY OFF`, `&RESTART OFF`, `&RESTART_HISTORY OFF`,
`&ENERGIES OFF`, `&FORCES OFF`) so the measured time reflects compute, not
disk traffic.  System size is varied via `MULTIPLE_UNIT_CELL` (each unit
cell = 64 H₂O):

| `N_MOLECULES` | `MULTIPLE_UNIT_CELL` | Atoms |
|---:|:---|---:|
| 64 | 1 1 1 | 192 |
| 256 | 2 2 1 | 768 |
| 512 | 2 2 2 | 1 536 |
| 1024 | 4 2 2 | 3 072 |
| 2048 | 4 4 2 | 6 144 |
| 4096 | 4 4 4 | 12 288 |

### NNP model used

The neural network potential used in every MD run is the bulk-water committee
NNP of **Schran, Brezina & Marsalek**, *J. Chem. Phys.* 153 (2020) 104105
(DOI: 10.1063/5.0016004), generation 4 (8-member committee, files at
`data/NNP/bulkH2O-jcp2020-cnnp/nnp-{1..8}/`).  The underlying reference DFT
was **revPBE-D3** (revised PBE [Zhang & Yang, PRL 80 (1998) 890] + Grimme D3
dispersion [Grimme et al., JCP 132 (2010) 154104]).  Dispersion is therefore
implicit in the NN parameters; no explicit `&VDW_POTENTIAL` block is included
at runtime.

For the figure-S4 viscosity-replication runs (separate workstream), the NNP
is the **Morawietz et al.** 2016 *RPBE-D3* model (PNAS 113, 8368, 2016),
at `cp2k-benchmarks/potentials/RPBE-vdW-2016/`.

## 7. Measurement methodology

### Timing metric

Time per MD step is extracted from CP2K's internal timing table as

```
t_per_step  =  qs_mol_dyn_low_AVG  /  STEPS
```

where `qs_mol_dyn_low_AVG` is the average wall time spent in the MD-loop
subroutine across all MPI ranks, as printed in the CP2K `--print-level low`
output.  Using `qs_mol_dyn_low` rather than the total job wall time excludes
startup costs (MPI_Init, NNP weight loading, MD-section setup) so that the
reported time reflects only the steady-state MD loop and is independent of
the chosen `STEPS` value.

### Repetitions and warm-up

For every (branch, configuration) data point the same `srun` is launched
**N_REPS + 1 = 6** times.  The first run is discarded as a warm-up (cold
filesystem cache, lazy MPI initialisation, etc.).  The remaining
`N_REPS = 5` runs constitute the timed sample.

### STEPS per sweep

| Sweep | STEPS | Rationale |
|---|---:|---|
| Size scaling (full-node, N = 64 → 4096) | 100 | Modest cost at full-node parallelism (largest size is ~5 min/rep on master) |
| Strong scaling at N = 1024 (1 → 76 cores) | 50 | Halved to keep the 1-core master point tractable: O(N²) in master makes single-core ~30 s/step at this size |
| Strong scaling at N = 256 | 100 | Smaller workload, full STEPS fits easily |
| Strong scaling at N = 2048 | 50 | Same reasoning as N = 1024 |
| OMP thread sweep (MPI = 1, N = 64) | 100 | Light workload |
| OMP × size 2D sweep (MPI = 1) | 20 | Reduced so the OMP=1 / large-N points fit in 1 h |

### Two CSV outputs per sweep

Each per-branch run produces:

- `results_<sweep>_<label>_<timestamp>_raw.csv` — every individual rep, columns
  `<config>,rep,time_per_step_s,walltime_s`
- `results_<sweep>_<label>_<timestamp>.csv` — summary, columns
  `<config>,n_reps,time_per_step_mean_s,time_per_step_std_s,time_per_step_min_s,
  walltime_mean_s,walltime_std_s,walltime_min_s[,speedup,parallel_efficiency]`

The plotting + statistical analysis scripts read both: summary CSVs for
visualisation, raw CSVs for the distribution and significance-testing sections
(see Section 8).

## 8. Statistical analysis

All quantitative claims about performance comparisons are made via the
following procedure, in accordance with Hoefler & Belli SC'15:

### Error bars (Rule 5)

Plotted error bars throughout this report are **95% confidence intervals**
on the mean, computed from the per-rep sample standard deviation `s` and
sample size `n = 5` via Student's t-distribution:

`CI₉₅ = t(n − 1, 0.025) · s / √n`

with `t(4, 0.025) ≈ 2.776`.  Raw sample standard deviations are
also reported in the summary tables for transparency, but the CIs are the
quantity the figure error bars depict.

### Distribution (Rule 8 / 12)

For the headline workload (N = 1024 H₂O on a full Ice Lake node) the raw
5-rep distribution is shown directly as a combined box-plus-violin plot
([`headline_distribution_N1024.pdf`](headline_distribution_N1024.pdf)), with
each individual rep overlaid as a scatter point.  This shows distribution
*shape* rather than just summarised mean ± CI.

### Pairwise significance testing (Rule 7)

For every pair of branches, the null hypothesis "the two branches have
equal mean time-per-MD-step" is tested using:

- **Welch's two-sample t-test** (parametric; allows unequal variances) — primary
- **Mann-Whitney U test** (non-parametric, distribution-free) — backup, robust
  if the per-rep samples are non-normal
- **Shapiro-Wilk normality test** on each branch's reps — informs which of the
  above to trust.  With `n = 5` the Shapiro-Wilk test has low power; the
  p-values are reported for transparency, not as definitive evidence of
  normality.

Effect sizes are reported as Cohen's d (with the pooled standard deviation):
`d = (mean_a − mean_b) / s_pooled`, classified per the conventional bins
(small d = 0.2, medium d = 0.5, large d = 0.8, very large d ≥ 1.2).

Full pairwise table at the headline workload:
[`statistical_comparison_N1024.txt`](statistical_comparison_N1024.txt).

## 9. Known caveats and blockers

| Issue | Effect on reported results | Mitigation |
|---|---|---|
| `kernel.perf_event_paranoid = 2` | MAQAO's hardware-counter metrics (`speedup_if_L1`, `nb_loops_80_if_L1`, DRAM throughput) appear as "Not Available" in the global metrics CSV | The static `array_access_efficiency` metric used in the report is derived by MAQAO's CQA module from assembly-level analysis, not PMU counters, and is therefore unaffected |
| `kernel.yama.ptrace_scope = 2` (since 2026-05-15) | Blocks `ptrace(PTRACE_ATTACH)`, which MAQAO's LProf sampler uses to attach to MPI child processes. All post-2026-05-15 MAQAO re-runs fail with `Critical: ptrace cannot attach to application process: Operation not permitted` | Report uses the May-14 MAQAO measurements (last successful run before the mitigation). Native-spline-omp is used as the cache-conscious representative because it already contained the cache opts at the time of that measurement (post-port pure-MPI native-spline has identical cache-relevant algorithmic content; see [`array_access_efficiency_disclaimer.md`](../maqao_plots/array_access_efficiency_disclaimer.md)) |
| Valgrind 3.22 cannot decode all AVX-512 subsets emitted by `-xCORE-AVX512` binaries | Cachegrind exits at startup on the new binaries with *"Please verify that both the operating system and the processor support Intel(R) AVX512F, AVX512DQ, ADX, AVX512CD, AVX512BW and AVX512VL instructions"* | Cachegrind path not pursued — MAQAO's `array_access_efficiency` (52.9 % → 78.2 %) carries the cache-locality argument |
| Optimized branches' `-xCORE-AVX512` overridden by `-xHOST` (CMake injection) | None — both flags resolve to identical AVX-512 ISA on the build host | Documented in Section 4 above |
| Strong-scaling N = 1024 1-core run on master takes ~3 h | First submission timed out at 2 h; required walltime of ≥ 6 h | `strong_scaling_only_slurm.sh` walltime bumped to 6 h; 1-core master point successfully captured at the second attempt (29430434) |

## 10. Reproducibility — script + data provenance

| Artifact | Path |
|---|---|
| Benchmark driver scripts | `cp2k-benchmarks/scripts/CSD3_benchmark_scripts/scaling/` |
| Build scripts | `cp2k-benchmarks/cp2k_*/CSD3_build_scripts/` |
| Cached binaries | `/rds/user/$USER/hpc-work/cp2k_binaries/csd3/<branch>/cp2k.psmp` |
| Raw + summary CSVs | `/rds/user/$USER/hpc-work/cp2k-benchmarks/results/cp2k_*/NNP/NNP_*_scaling_*/` (mirrored to `cp2k-benchmarks/results/`) |
| MAQAO experiment dirs | `/rds/user/$USER/hpc-work/cp2k-benchmarks/maqao/xp_*` |
| Plotting and statistics | `cp2k-benchmarks/plots/plot_csd3_scaling.py`, `maqao_data_analysis.py` |
| This document | `cp2k-benchmarks/plots/methods.md` |

This document generated 2026-05-16 from the CSD3 environment described above.
