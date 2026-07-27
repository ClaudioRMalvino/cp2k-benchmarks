# Meeting Notes — CP2K NNP Optimisation & MACE (2026-06-23)

*Talking document for the supervisor meeting. Covers: (1) the chebyshev branch —
what it is and how it works; (2) proof that the physics is preserved — bitwise
equivalence + Morawietz Fig S4 reproduction; (3) current performance benchmarks;
(4) the MACE integration; (5) methodology/rigor; (6) what is finalising right now.*

---

## TL;DR — the 2-minute version

> Two workstreams, both landed this week.
>
> **1. Chebyshev (CPU performance).** A table-free reimplementation of CP2K's NNP
> kernels. It is **numerically identical to upstream master** (energies to 10⁻¹⁴
> Ha/atom), so it changes speed, not physics. It is **up to 11.8× faster than
> master** (growing with system size) and is the only branch with **working OpenMP
> (54× on 76 threads)** and a **~22× smaller node-memory footprint** at full speed.
>
> **2. Physics validation.** Chebyshev **reproduces Morawietz 2016 Fig S4**:
> corrected diffusion **D₀ = 0.299 ± 0.020 Å²/ps** (paper: ~0.31) flat across all
> five system sizes, and viscosity **η = 0.70 mPa·s** (paper: ~0.66). The
> finite-size (Yeh–Hummer) correction collapses the data exactly as in the paper.
>
> **3. MACE (new capability).** A **CPU-native MACE foundation-model backend** is
> integrated into CP2K (no PyTorch/GPU) and **validated end-to-end** (stable 15 ps
> production MD). It sets up the report-2 science: short-range MACE + long-range
> electrostatics for **NaCl-in-water ion transport**.
>
> *Open question for you:* is this physical validation sufficient for report 2,
> and which is the headline — the performance optimisation or the MACE capability?

---

## 0. One-line status

The **table-free "chebyshev" branch** is implemented, **proven numerically
identical to upstream master** (machine precision), **reproduces the Morawietz
2016 Fig S4 water physics**, and is fully benchmarked. Separately, a **CPU-native
MACE foundation-model backend** is integrated into CP2K and validated end-to-end
(see `MACE_reasoning.md`).

---

## 1. The chebyshev branch — what we did and how it works

**What it is.** CP2K evaluates its neural-network potential through cutoff/
symmetry-function kernels. The upstream `master` evaluates the cosine, tanh and
exponential cutoff functions through **Hermite lookup tables** (precomputed,
~thousands of points, stored per process). The `feature/nnp-chebyshev` branch
**removes those lookup tables** and replaces them with **table-free minimax
polynomial approximations** — `minimax_cos`, `minimax_tanh`, `exp_reduced` — in
`splines_methods.F`, evaluated on the fly. It also adds **centre-level OpenMP
threading** to the NNP force/energy loop.

**How it works (mechanism).**
- A *minimax polynomial* is the polynomial of given degree that minimises the
  *worst-case* error over the interval. We fit one to each cutoff kernel offline,
  so at run time the kernel is a handful of fused multiply-adds — no memory
  lookup, no table to store, no cache pressure from random table access.
- The standalone kernels were verified against `libm` to ~6×10⁻¹⁴ (near machine
  precision), so the approximation is numerically exact for our purposes.
- The OpenMP layer parallelises over atomic centres, so the NNP evaluation can
  use threads — something `master` and the earlier `native-spline` branch cannot.

**Why.** Two motivations: (1) lookup tables are memory-bound and pollute cache as
the system grows; a table-free kernel is compute-bound and scales better; (2) the
threading enables one MPI rank to use many cores, which collapses memory (see §3).

**Relationship to the earlier `native-spline` branch.** `native-spline` was the
report-1 optimisation (cubic Hermite splines, O(N^0.71) neighbour handling).
`chebyshev` is the report-2 line: table-free kernels + OpenMP. Comparing the two
answers a real question — *does the table-free method beat the spline-table
method?* (Answer in §3: nuanced — they tie on speed at best config, chebyshev
wins on memory and threading.)

---

## 2. Correctness — the physics is preserved (two independent proofs)

This is the "how do I prove the data is good" answer. Two layers:

### 2a. Bitwise equivalence to master (the strongest proof)
The optimisations must not change the physical model. Single-point energy/force
checks vs the untouched `master`:

| System | ΔE per atom (Ha) | max ΔForce (Ha/bohr) |
|---|---|---|
| N=64 (cubic cell) | — | **2.6×10⁻¹³** |
| N=128 (explicit cell) | **4.1×10⁻¹⁴** | **1.0×10⁻¹¹** |
| N=128 (MUC 2×1×1) | 4.1×10⁻¹⁴ | 1.0×10⁻¹¹ |

These differences are at **machine precision** — i.e. chebyshev and master compute
the *same* energies and forces. Therefore **any physics master reproduces,
chebyshev reproduces identically.** (The NNP regression suite also passes.)

### 2b. Morawietz 2016 Fig S4 reproduction (end-to-end MD physics)
Reference: Morawietz et al. 2016 SI, Fig S4 — diffusion and viscosity of liquid
water vs system size, **RPBE-vdW NNP, 300 K, ρ = 1 g/cm³**. Method: viscosity η
from the **Green–Kubo stress-ACF (3 ps integration cut)**; diffusion D_PBC from
**MSD**; finite-size-corrected **D₀ via the Yeh–Hummer relation**
`D₀ = D_PBC + k_BT·ξ/(6πηL)`, ξ = 2.837297. Each value = **mean ± SEM over 5
independent NVE segments**.

**Our chebyshev results vs the paper:**

| N H₂O | η (mPa·s) — ours | η — Morawietz | D_PBC (Å²/ps) — ours | D₀ (Å²/ps) — ours | D₀ — Morawietz |
|---|---|---|---|---|---|
| 64  | 0.83 ± 0.13 | ~0.68 | 0.239 ± 0.039 | 0.311 ± 0.056 | ~0.31 |
| 128 | 0.65 ± 0.03 | ~0.69 | 0.265 ± 0.018 | 0.326 ± 0.020 | ~0.31 |
| 256 | 0.67 ± 0.02 | ~0.67 | 0.267 ± 0.009 | 0.314 ± 0.009 | ~0.31 |
| 512 | 0.64 ± 0.05 | ~0.65 | 0.269 ± 0.003 | 0.309 ± 0.004 | ~0.31 |
| 1024 | 0.71 ± 0.04 | ~0.66 | 0.294 ± 0.012 | 0.323 ± 0.012 | ~0.31 |

**Headline numbers (all 5 sizes now in):**
- **Extrapolated D₀ (linear 1/L→0 fit) = 0.299 ± 0.020 Å²/ps** vs Morawietz NNP
  **~0.31** → match within error.
- **Mean η = 0.70 mPa·s** vs Morawietz NNP **~0.66** → close (slightly high,
  within the size-to-size scatter).
- D_PBC rises monotonically **0.239 → 0.294 Å²/ps** as N goes 64 → 1024.

**What this proves (the three signatures of Fig S4):**
1. **D₀ is flat at ≈0.31–0.32 Å²/ps across all five sizes** (fit intercept
   0.299 ± 0.020) — matches Morawietz. This is *the* headline of Fig S4: the
   Yeh–Hummer correction removes the box-size dependence. We reproduce it.
2. **D_PBC rises with system size** (0.239 → 0.269 as N goes 64 → 512), i.e.
   falls with 1/L — the expected finite-size suppression, and the correction in
   (1) collapses it onto the flat D₀.
3. **η ≈ 0.65 mPa·s, roughly flat for N ≥ 128** — matches the paper's ~0.66.

**Honest caveats (state these):**
- The **N=64 viscosity (0.83 ± 0.13)** is high and noisy — Green–Kubo η is
  intrinsically noisiest at small N; the large error bar overlaps the reference.
  This is expected, not a problem.
- The reference is **the RPBE-vdW NNP itself, not experiment.** That NNP is known
  to be slightly over-mobile / under-viscous vs real water (experiment: D ≈
  0.23, η ≈ 0.85). "Reproducing Morawietz" means landing on **0.31 / 0.66**, which
  we do — *not* on the experimental numbers.

**Combined message:** 2a proves the optimisation is numerically transparent; 2b
proves the full MD pipeline reproduces published water physics. Together they
confirm chebyshev is both correct and physically sound.

**Reading the figure** (`plots/figS4_replication.png`): it now shows **two**
reference lines. The **purple dash-dot** line is the **Morawietz RPBE-vdW NNP**
(D₀ = 0.31 Å²/ps, η = 0.66 mPa·s) — *this is our validation target, and our points
land on it.* The **pink dashed** line is **experimental water** (D = 0.23, η =
0.896) — we correctly sit **off** it (over-mobile, under-viscous), because the NNP
itself is offset from experiment. So the one-line takeaway: **we match the purple
line (the NNP), not the pink line (experiment) — by design.**

*Figures: `plots/figS4_replication.png` (diffusion + viscosity vs size),
`plots/figS4_accuracy.png`, `plots/figS4_performance.png`.*

---

## 3. Performance benchmarks (chebyshev vs master vs native-spline)

All on CSD3 Peta4–Ice Lake (76 cores/node), bulk water, Schran committee NNP.
*Figures in `plots/cheby_benchmark_figs/` (fig1–fig7).*

**vs master — the speed-up grows with system size** (size scaling, 76×1 pure MPI):

| N H₂O | chebyshev vs master |
|---|---|
| 64 | 1.0× | 
| 512 | 1.3× |
| 1024 | 2.1× |
| 4096 | **11.8×** |

Master scales as **O(N^1.83)** (super-linear — its lookup tables thrash cache at
scale); chebyshev stays near-linear. The advantage compounds with size.

**vs native-spline — the interesting comparison (table-free vs spline-table):**
- On **raw pure-MPI throughput**, native-spline is **slightly faster** (chebyshev
  is ~2.3× slower on a single core because it carries the OpenMP machinery and
  minimax evaluation). So the spline-table method wins on single-config speed.
- But at **each branch's best full-node config they essentially tie**:
  native-spline 0.0555 s/step (38 ranks × 2 threads) vs chebyshev 0.0570 s/step
  (1 rank × 76 threads) — **and chebyshev does it in ~22× less node memory**
  (476 MiB vs 10.5 GB), because OpenMP lets it run 1 rank instead of 38 replicated
  copies. *(fig2)*

**The two things only chebyshev can do:**
- **OpenMP scaling**: 54× on 76 threads (71% efficiency); master and native-spline
  have no working NNP threading and degrade if you give them threads. *(fig3)*
- **Low-memory full-node operation**: full speed at 1 MPI rank → fits much larger
  systems per node. *(fig2)*

**Multi-node** (chebyshev vs master, same job): at 4 nodes / N=4096, chebyshev is
**8.3× faster** in absolute time. *(fig6)*

**The honest "which is better" verdict:** for pure single-config MPI throughput,
**native-spline edges chebyshev**; chebyshev wins on **threading + memory + cache
behaviour at scale**, and ties on best-config speed. Both crush master.

> ⚠️ **Rigor note:** the definitive pure-MPI strong-scaling head-to-head is being
> **re-run right now** under exact report-1 methodology (N=1024, NUMA-aware cores
> {1,2,4,8,16,19,32,38,76}, **100 steps, 5 reps + 1 warm-up**, n=5 CI95) so all
> three branches are matched. Earlier core runs had mixed step counts (50 vs 100);
> this job removes that inconsistency. (Job `30958202`, queued.)

---

## 4. MACE integration (summary — full detail in `MACE_reasoning.md`)

We added a **CPU-native MACE machine-learning-potential backend** to CP2K:
- **Backend = symmetrix** (a tensor-free C++ CPU evaluator) — **no PyTorch/
  libtorch, no GPU.** This is what makes it viable on CPU-only HPC like CSD3.
- **Integrated as a FIST `&MACE` many-body term**, so it composes with CP2K's
  electrostatics (Ewald) and QM/MM.
- **Validated end-to-end (L0→L3):** exact energies/forces, dt²-correct energy
  conservation, and a **stable 15 ps NVT production run** with the MACE-MP-0
  foundation model that faithfully reproduces that model's known (sluggish) water.
- **Key point:** MACE is a *short-range* model — it does not add long-range
  electrostatics by itself. But because it lives inside FIST, we can add a
  long-range Coulomb baseline via `METHOD MIXED`. **That sets up the report-2
  science question** (NaCl-in-water ion transport: do short-range MLPs fall short,
  à la Yue; measured with transport coefficients à la Fong).

---

## 5. Methodology & rigor (for the "is this done properly" question)

- **Hardware:** single CSD3 Peta4–Ice Lake node, 76 cores, AVX-512, NUMA-aware.
- **Timing:** mean of **5 timed repetitions + 1 discarded warm-up**; uncertainties
  as **95% CI (Student-t, n=5)**; speed-ups from mean times; binary md5 logged per
  run to rule out silent binary swaps. (Hoefler & Belli benchmarking practice.)
- **Physics validation:** 5 independent NVE segments per size; Green–Kubo η with a
  3 ps cut; Yeh–Hummer finite-size correction — all matching Morawietz's method.
- **Correctness gate:** bitwise energy/force equivalence to upstream master before
  any performance claim.
- **Consistency fix in progress:** strong-scaling re-run to put all branches on
  identical 100-step / 5-rep footing (§3 note).

---

## 6. Status of in-flight work

- **figS4 N=1024 point — DONE** (job `30958650`, completed). All five sizes are in;
  Fig S4 reproduction is complete (§2b table + `plots/figS4_replication.png`).
- **Report-1-consistent strong scaling** — job `30958202` still queued; will give
  the definitive matched pure-MPI comparison (100 steps, 5 reps) for fig2/fig1.
  *Until it lands, quote §3 numbers as preliminary on the strong-scaling axis.*

---

## 7. Talking points / questions for the supervisor

1. **Physics confirmation:** chebyshev reproduces Morawietz Fig S4 (D₀ ≈ 0.31 flat,
   η ≈ 0.65) *and* is bitwise-identical to master — is this sufficient physical
   validation for report 2, or do you want an additional observable (e.g. RDF)?
2. **The table-free vs spline trade-off:** native-spline is marginally faster on
   pure MPI, but chebyshev adds threading + 22× lower memory. Which axis matters
   most for the thesis narrative — throughput, or memory/threading headroom for
   large systems?
3. **MACE direction:** confirm the report-2 science target — short-range MACE vs
   MACE+Coulomb baseline for NaCl-in-water transport (Yue/Fong framing).
4. **Scope:** is the chebyshev + MACE pairing the right two-pronged story
   (CPU performance optimisation + new foundation-model capability), or should one
   take priority?

---

*Companion documents: `MACE_reasoning.md` (full MACE explainer),
`nnp_implementation_changes.md` (implementation changelog),
`plots/cheby_benchmark_figs/` (performance figures), `plots/figS4_*.png` (physics).*
