# Report 2 — Plan: NaCl(aq) diffusion + Dhruv-vs-mine comparison

**Status: PLANNING ONLY. Nothing in here has been built, queued, or run.**
Drafted 2026-06-23. Supervisor follow-up pending (he was rushed; questions in Part C).

Context: supervisor has redirected away from further NNP optimization. Production
simulations for report 2 will use **Dhruv's code** (`pr/nnp-cpu-cell-list`, which is
going into a CP2K PR). My own optimization branch (`feature/nnp-mace`) is retired from
the production path but gets (a) a comparison write-up and (b) a head-to-head benchmark
to see whether his code is actually faster than mine.

---

## Part A — Code comparison: `feature/nnp-mace` (mine) vs `pr/nnp-cpu-cell-list` (Dhruv)

**Headline:** we independently built the *same* neighbour-search primitive — a linked-cell
+ Verlet-skin list with pre-expanded periodic images (shared `nnp_prepare_cell_list_cache`
name; my `feature/nnp-verlet-cells` lineage). The algorithmic core is shared, not opposed.
We diverge in everything built on top.

| Axis | Dhruv (`pr/nnp-cpu-cell-list`) | Mine (`feature/nnp-mace`) |
|---|---|---|
| **Cell list** | Factored into 2 new modules (`nnp_cell_list.F`, `nnp_neighbor_interface.F`); `VERLET_SKIN` keyword | Inlined into `nnp_acsf.F`; **+ halo-filtered image table** (~9× less memory, ~3N vs 27N images); **+ PBC correctness fix** for axes longer than 2·cutoff (>0.1 Ha discontinuities in the N=128 box); **+ walk-bound buffer sizing** |
| **Radial ACSF** | **Spline-tabulated** radials (`RAD_SPLINE_N=8192` knots) + `!$OMP SIMD` — memory for branch-free vectorized lookups | **Opposite choice: table-free minimax/Chebyshev** (Clenshaw recurrence, ≤1.3e-13 accuracy) — removes ~480 KiB tables/rank, more FLOPs, relies on inlining/`-ipo` |
| **Parallelism** | Centre-level `OMP DO` + SIMD on tables; MPI unchanged | Centre-level `OMP DO` + **shrinks per-step MPI all-reduce** (786 KiB → 64 B/step) |
| **Force/backprop** | Sparse per-neighbour force scatter | **Reverse-mode adjoint** backprop (DGEMV vs DGEMM) + scatter + build element-sorted index maps once |
| **Input keywords** | `RAD_SPLINE_N`, `VERLET_SKIN` | None for NNP; new surface is the `&MACE` backend |
| **Tests/docs** | **Bit-exact MD regtest + user docs + dedicated code-review commit** | `mace_unittest.F` + benchmark tree; no NNP regtest/docs |
| **Scope** | Tight, upstream-PR-shaped; no MACE | Broader/experimental + full MACE manybody backend; less PR-ready |

**Bottom line:** Dhruv built a *shippable* cell-list NNP (modular, tabulated, keyworded,
regtested). I built a *faster-and-broader* one (table-free, leaner MPI, memory-frugal, a
PBC fix he lacks, plus MACE) but inlined with weaker test/doc scaffolding.

**Key consequence for Part B:** because the neighbour core is shared, the benchmark is *not*
"his cell list vs mine." It mainly measures **his spline-table lookups vs my table-free
minimax**, plus my **MPI-collective shrink** (which should favour me at high rank counts /
across sockets). Clean, well-posed performance question.

---

## Part B — Performance benchmark: Dhruv vs `nnp-mace`

Two-series comparison only — **Dhruv's `pr/nnp-cpu-cell-list` vs my `feature/nnp-mace`**.
No master, no native-spline-omp baseline (`nnp-mace` ≡ `nnp-chebyshev` on the NNP path; the
MACE backend is idle on a pure-NNP run, so NNP performance is identical to chebyshev).

### B.0 — Build flags / comparability (the gating issue)

My existing `nnp-chebyshev` benchmark CSVs **cannot be reused as-is**, because Dhruv likely
built with different compilation flags. Mixing flag regimes would confound the
code comparison with a build comparison. So `nnp-mace` must be **rebuilt with a flag set
matched to Dhruv's** before any timing is comparable.

**Finding — his flags ARE documented, in his thesis Appendix A.4 (Compile-flag audit) + A.2.**
- **Toolchain (A.2):** Intel oneAPI **ifort 2022.1.0**, **MKL 2022.1.0**, **Intel MPI 2021.6.0**.
  CP2K master commit `4daa2013c2`; his optimised tree `b2e3918a2b`; upstreamed as **CP2K PR #5295**.
  External ref: LAMMPS `457ea8c1` (30Mar2026) + n2p2 `29b9c9f1`.
- **Platform (A.1):** CSD3 Icelake SL2, 2× Intel Xeon Platinum 8358 (76 cores), 256 GiB DDR4, HDR200.
- **Benchmark profile (A.4):** he built 4 bit-exact `ifort` profiles and benchmarked with the
  fastest — **`-O3 -xCORE-AVX2 -fp-model=precise`** (130 ms/step on H₂O-1024). The audit table:

  | Profile | flags | t (ms) | bit-exact |
  |---|---|---|---|
  | -O2 | `-O2 -g -DNDEBUG` | 169 | PASS |
  | **-O3 / AVX-2 / precise (USED)** | `-O3 -xCORE-AVX2 -fp-model=precise` | **130** | PASS |
  | -O3 / AVX-512 / zmm-high | `-O3 -xCORE-AVX512 -qopt-zmm-usage=high` | 218 | PASS |
  | -O3 / AVX-512 / zmm-high / precise | `-O3 -xCORE-AVX512 -qopt-zmm-usage=high -fp-model=precise` | 167 | PASS |

- **He deliberately used AVX-2, not AVX-512:** "AVX-512 underperformed on H₂O-1024, consistent
  with AVX-512 frequency behaviour on this Intel server class and a bandwidth-bound roofline."

**My existing build used** (`cp2k_optimized/build/CMakeCache.txt`) — **NOT comparable**:
- Compiler: Intel `mpiifort` **2021.6** (vs his 2022.1.0), `CMAKE_BUILD_TYPE=Release`
- `CMAKE_Fortran_FLAGS = -O2 -g -xCORE-AVX512 -qopenmp -funroll-loops -ftree-vectorize` (+ `-O3`)
- Three mismatches: compiler version, **ISA (my AVX-512 vs his AVX-2)**, opt flags. Worse, my
  AVX-512 choice is the *slower* one on this exact hardware per his audit — a naive head-to-head
  would hand him a win that is partly just his better ISA choice.

> ACTION: rebuild `nnp-mace` with **Dhruv's exact benchmark profile** — ifort 2022.1.0,
> `-O3 -xCORE-AVX2 -fp-model=precise`, MKL 2022.1.0, Intel MPI 2021.6.0 — before any timing.
> (Optionally also build an AVX-512 variant of mine to show whether my MPI-collective shrink
> changes the AVX-512-vs-AVX-2 trade-off, but the headline figure must be AVX-2 vs AVX-2.)

### B.1 — Benchmark protocol (reuse report-1 rig verbatim)
- System: H₂O-N NNP MD, same decks as report 1.
- Rigor: **100 steps, 5 timed reps + 1 discarded warm-up**, n=5 Student-t CI95.
- Core ladder (NUMA-aware): {1, 2, 4, 8, 16, 19, 32, 38, 76}.
- md5 binary logging; both binaries on the same nodes.

### B.2 — Axes (both engines support all three)
1. **Size scaling** — 76×1, N = 64 → 4096.
2. **Pure-MPI strong scaling** — N = 1024, the core ladder above.
3. **Hybrid MPI×OMP** — the fairest fight, since *both* have OpenMP (the old
   master/native-spline lacked it, which made those comparisons awkward).

### B.3 — Accuracy guard
Confirm both binaries still reproduce the n2p2 reference forces (Dhruv's bit-exact regtest;
my minimax branch diverges only at machine precision) so any speed delta is read against an
equal-accuracy baseline.

### B.4 — Output
Into existing `results/` + `plots/` convention; Cambridge palette; **two-series** figures
(size, strong-scaling, hybrid). Same aesthetic as `thesis_figures.py`.

---

## Part C — NaCl(aq) science campaign (needs supervisor input)

Supervisor's stated wants: **(i)** diffusion as a function of **concentration**, **(ii)**
test the **thermodynamics of an ionized system**. My task #1: a **Morawietz-Fig-S4-style
finite-size diffusion correction**, but for NaCl(aq). These compose into one 2-D grid
(concentration × box-size).

### C.-1 — What Dhruv ALREADY did (thesis §4.4) and why I'm not redundant
His thesis is a *performance* thesis; NaCl(aq) is the demo that his speed makes the
finite-size sweep affordable — NOT an electrolyte-science study. Concretely he did:
- **One state point only:** 1 mol/kg, 300 K. Four boxes N={1008,1685,3108,6218}, cell sides
  21.6–39.6 Å, 6/10/18/37 Na⁺/Cl⁻ pairs — purely to populate the 1/L axis at fixed conc.
- **Water-oxygen self-diffusion only** (D∞=1.68e-5 cm²/s, R²=0.97, ⟨D/D₀⟩=0.886 vs O'Neill 0.88).
  **He never computed Na⁺ or Cl⁻ diffusion** — the actual ions.
- **η was *fitted* from the 1/L slope** (η_fit=0.94 mPa·s), not computed — explicitly "in the
  absence of a separate Green-Kubo calculation," as an internal consistency check.
- **RPA committee only**; ~200 ps trajectories; 3 replicas/box except N=6218 (1 replica). He
  flags it as "a consistency check rather than a final transport benchmark."
- **No concentration dependence, no ion pairing, no thermodynamics.**

| Christoph's ask | Dhruv did | My extension |
|---|---|---|
| **Diffusion vs concentration** | Single 1 mol/kg point | Sweep molalities (e.g. 0.5→5 m) → D(c) curves |
| Ion transport | Water-O only | **Na⁺/Cl⁻ self-diffusion** + collective transport (Onsager L_ij, conductivity — Fong's domain; she coauthored the model) |
| **Thermodynamics of ionized system** | None | Ion pairing / PMF, RDFs & coordination, activity |
| Rigor | η fitted from slope | **Independent Green-Kubo η** (`compute_viscosity.py` already exists) |

**Pitch to Christoph:** "Dhruv proved the optimized path makes the finite-size diffusion sweep
affordable, but only as a single-concentration, water-diffusion consistency check with a fitted
viscosity. I'd use his fast engine to do the actual electrolyte science — ion diffusion and
transport vs concentration, with an independently-computed viscosity, plus ion-pairing
thermodynamics." Additive, not duplicative.

### C.0 — Engine & inputs (fixed regardless of his answers)
- Run on **Dhruv's binary** (his explicit instruction; PR-bound).
- O'Neill C-NNP deck: `METHOD MIXED` = short-range committee-NNP (8 members) + FIST/Ewald-
  SPME Coulomb baseline; H mass = 2 amu (deuterium; O'Neill/Morawietz protocol).
- Models: O'Neill et al. (niamh-on/nacl-water), 4 elements `O H Na Cl`, 320 SFs, cutoff_type 2.
  **Coauthor Kara Fong** — the exact author the supervisor flagged.
- Analysis: adapt the water-only figS4 suite (`compute_msd.py` etc.) to **per-species MSD**
  (Na⁺, Cl⁻, water-O), each with its own Yeh-Hummer correction
  `D₀ = D_PBC + k_BTξ/(6πηL)` (ξ=2.837297); reuse `compute_viscosity.py` (Green-Kubo η)
  unchanged.

### C.1 — Campaign shape (2-D grid)
1. Build NaCl(aq) boxes at each target concentration (insert Na⁺/Cl⁻ into equilibrated
   water, neutralize, set density).
2. For each concentration: a **box-size ladder** (the figS4 axis) → NVT equilibration →
   NVE production segments.
3. Per (concentration, size): per-species `D_PBC` + η from stress-ACF.
4. Yeh-Hummer extrapolate `D_PBC` vs 1/L → `D₀` per species → figS4 analogue, then collapse
   to **D₀(concentration)** curves for Na⁺ / Cl⁻ / water.
5. Thermodynamics block (definition TBD — see questions).

### C.2 — Open questions for supervisor (tomorrow)
1. **Concentrations** — which molalities, what range (dilute ~0.5 m → near-saturation ~5–6 m),
   how many points?
2. **Diffusion target** — self-diffusion `D_i` only, or also **collective transport**
   (Onsager `L_ij` / Maxwell-Stefan, ionic conductivity)? This is the Kara-Fong angle and
   changes the analysis substantially.
3. **Finite-size correction** — standard per-species Yeh-Hummer, or the **multicomponent /
   mixture** correction (Jamali–Vlugt) for electrolytes? Run the figS4 finite-size study at
   *every* concentration, or just one reference concentration?
4. **"Thermodynamics of an ionized system"** — which observables: ion-pairing free
   energy / PMF, RDFs & coordination numbers, activity coefficients, Kirkwood-Buff integrals?
5. **Functional** — revPBE0-D3 (O'Neill production hybrid) default, or another reference
   (MP2 / RPA / r2SCAN / optB88-vdW / revPBE-D3)?
6. **State point / ensemble** — NVT at experimental electrolyte density per concentration,
   or NPT to predict density?
7. **Sizes & trajectory lengths** — figS4 water ladder was N=64→1024; what minimum box and
   how many ps/segment for converged ion statistics (ions are dilute and slower → need
   longer/more sampling)?
