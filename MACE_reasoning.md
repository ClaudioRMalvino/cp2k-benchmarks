# MACE in CP2K — What We Did, How, Why, and the Benefits

*A self-contained explainer of the `&MACE` machine-learning-potential backend
added to CP2K on the `feature/nnp-mace` branch. Written to be read end-to-end
and used to explain the work to a supervisor.*

---

## 0. TL;DR (the one-paragraph version)

We added a **new machine-learning-potential (MLP) backend to CP2K** that runs a
**MACE** model as a force engine for molecular dynamics. We did it *without*
PyTorch/libtorch by using **symmetrix**, a lightweight C++ CPU evaluator, and we
wired MACE in as a **FIST many-body term** (`&MACE`) so it plugs into CP2K's
existing force-field machinery. We validated it end-to-end — from a low-level
binding test up to a stable 15-picosecond molecular-dynamics production run with
a real foundation model — and it works. The point: this brings **foundation
machine-learning potentials into CP2K on CPU-only HPC**, and — because it lives
inside CP2K's force-field engine — it opens the door to **electrolyte / ion-in-
water science** (short-range MACE + long-range electrostatics), which is the
next research target.

---

## 1. What MACE actually is (and what it is *not*)

**MACE is a machine-learning interatomic potential.** It is a neural network
that predicts the energy and forces of a set of atoms directly from their
positions and chemical species, trained to reproduce expensive quantum-mechanical
(DFT/ab-initio) reference calculations. The payoff is *ab-initio-like accuracy at
a tiny fraction of the cost* — the whole reason MLPs exist.

Architecturally, MACE is an **equivariant message-passing graph neural network**
(higher body-order, rotation-equivariant features). This makes it more accurate
and more data-efficient than the older **Behler–Parrinello high-dimensional
neural network potential (HDNNP)** that CP2K already supports natively. **MACE
and HDNNP are different model families** — MACE is not "a better HDNNP," it is a
different, more expressive architecture.

**Crucial clarification — MACE is a SHORT-RANGE model.** Like HDNNP, MACE has a
finite cutoff radius (~5–6 Å). Its message passing extends the *effective* range
a little (≈ layers × cutoff), but it does **not** describe true long-range
electrostatics (the 1/r Coulomb tail, periodic Ewald sums, long-range polarisation).
This matters and is addressed in §8 — do not claim MACE "adds long-range
interactions" on its own. It does not.

**Foundation models.** The exciting part of MACE today is **MACE-MP-0**: a
*foundation model* pre-trained on a huge dataset spanning most of the periodic
table. It works "out of the box" on many elements (including Na, Cl, etc.)
**without training your own potential** — which is precisely what makes it useful
for new chemistry like salt water.

---

## 2. The problem — why we care

CP2K supports MLPs natively, but only through the **older HDNNP framework**, and
only for the chemistry someone has trained an HDNNP for (in our case: bulk water,
H and O). Two limitations follow:

1. **No access to modern foundation models.** State-of-the-art MLPs like MACE are
   far more accurate and transferable, and the foundation models cover the whole
   periodic table. CP2K had no way to run them.
2. **The standard way to run MACE needs PyTorch (libtorch) and is built for GPUs**
   — which is a poor fit for CPU-only HPC clusters (like CSD3's Ice Lake nodes)
   where building/matching libtorch + CUDA is painful or impossible.

Our project is a **CPU-focused performance-optimisation effort** for CP2K's NNP
machinery. Adding a **CPU-native MACE backend** is a natural extension: it gives
CP2K users modern foundation-model accuracy on the hardware they actually have,
and it positions us for the electrolyte science questions (salt in water) that
short-range MLPs are famously challenged by.

---

## 3. What we built — the high-level design

We added MACE as a **FIST `&NONBONDED` many-body term** in CP2K, activated by a
new `&MACE` input section and gated behind a `__MACE` compile flag. Two
deliberate design choices:

- **Backend = symmetrix, not PyTorch.** [symmetrix](https://github.com/) is a
  header-only, **tensor-free C++ evaluator** that runs a MACE model on CPU over a
  **CSR (compressed-sparse-row) neighbour graph**, using OpenMP for parallelism
  and the `sphericart` library for spherical harmonics. **No libtorch, no CUDA.**
  It consumes a MACE model exported to a **JSON** file (weights + architecture).

- **Integration point = FIST force field.** Rather than bolt MACE on as a
  standalone calculator, we made it a many-body potential *inside* CP2K's
  classical force-field engine. This is the key strategic choice (see §6 and §8):
  it means MACE composes with everything FIST already provides — Coulomb/Ewald
  electrostatics, QM/MM, `METHOD MIXED`, thermostats, enhanced sampling.

We patterned the whole thing on CP2K's existing NequIP/Allegro MLP interfaces, so
it looks like idiomatic CP2K code (built "as if for an upstream pull request").

---

## 4. How it works — the mechanics

The work is in three layers, committed in order on `feature/nnp-mace`:

### 4.1 Binding layer (commit `198408e1bd`)
A thin `ISO_C_BINDING` bridge between Fortran (CP2K) and C++ (symmetrix):
- `src/mace_api.F` + `src/mace_c_api.cpp` wrap symmetrix behind an opaque handle
  (`mace_model_type`). Entry points: load a model, query its cutoff, number of
  elements, atomic numbers, **compute energy/forces**, release.
- `cmake/modules/FindSymmetrix.cmake` finds the library; build with
  `-DCP2K_USE_MACE=ON`. A toolchain installer (`install_mace.sh`) fetches
  symmetrix + sphericart.
- symmetrix/sphericart are linked **statically** (`.a`), so the final `cp2k.psmp`
  needs no special runtime library path.

### 4.2 The force evaluator (commit `a83ef8eb15`)
`src/manybody_mace.F` is the heart of it. On every force evaluation it:
1. **Builds an edge list** from CP2K's neighbour list (centre atom *a*, neighbour
   *b*, displacement `r_b − r_a + cell_vector`, filtered by the cutoff).
2. **Packs it into symmetrix's CSR graph** (a counting-sort that groups edges by
   centre atom): node types, neighbour counts, neighbour indices, displacement
   vectors, distances.
3. Calls `mace_model_compute`.
4. **Scatters the result back**: per-edge forces onto the two atoms
   (`f[a] -= fe`, `f[b] += fe`), accumulates the virial (for pressure/stress),
   and sums per-atom energies into the total.

Element identity is resolved from the **model's own list of atomic numbers**, so
the order you list atoms in the input file does not matter — correctness is
automatic.

### 4.3 Units-consistency guard (commit `c3e371043d`)
A subtle, dangerous unit bug (see §7) led us to add a check: `read_mace_data`
**warns if `UNIT_FORCES ≠ UNIT_ENERGY / UNIT_LENGTH`**, because an inconsistent
force unit silently breaks energy conservation in MD while leaving single-point
energies looking correct.

### 4.4 Parallelism — an honest note
In the current design the **neighbour graph is replicated on every MPI rank**
(gathered via `allgatherv`), and the result is divided by the number of ranks.
Consequence: **MPI does not speed MACE up** — every rank does the whole system.
The only parallel lever today is **OpenMP inside symmetrix**. This is a known
limitation and the obvious target for future work (a particle/domain
decomposition would make MACE scale across ranks).

---

## 5. Validation — how we know it actually works

We validated in four escalating levels. **All passed.**

| Level | Test | Result |
|---|---|---|
| **L0/L1** | Binding + standalone unit test (symmetrix called directly) | Loads and evaluates correctly |
| **L2 — energy/forces** | 3-atom O–H–H single point, tiny model | Energy `0.031956406755371 Ha` — **exact** match to the symmetrix reference (15 digits). Forces **sum to zero** (momentum conserved). Out-of-plane forces vanish by symmetry. |
| **L2 — NVE drift** | 64-water energy-conserving MD, timestep ladder | Energy drift scales as **dt²** (ratios 3.98×, 4.00× per halving) — the signature that **forces are the exact analytic gradient of the energy** (F = −dE/dx). |
| **L3 — production** | **MACE-MP-0b3-medium** foundation model, 64-water, **15 ps NVT at 300 K** (30,000 steps), real run | **Stable.** ⟨T⟩ = 300.5 K with fluctuations matching equipartition theory (sd 18.1 K vs 17.7 K expected). Conserved-quantity drift **3.9×10⁻⁶ Ha/atom over 15 ps** — essentially flat. |

**The L3 physics result is worth understanding.** From the trajectory we measured
the water self-diffusion coefficient: **D ≈ 0.15×10⁻⁵ cm²/s**, about **15× lower
than experiment** (2.3×10⁻⁵). This is **not a bug** — it is the *well-documented*
behaviour of the MACE-MP-0 foundation model, which produces over-structured,
sluggish water. Our integration **faithfully reproduces the model's own known
(poor) dynamics**, which is exactly what you want from a correct implementation:
the pipeline is trustworthy, and the "wrong" number is the *right* wrong number.
(Caveat for honesty: D is indicative only — single run, no replicas, no
finite-size correction, NVT thermostat perturbs dynamics.)

**Bottom line: the `&MACE` path is validated end-to-end (L0 → L3).** It sustains
production-scale MD and yields real physical observables. The remaining work is
*scientific*, not *mechanical*.

---

## 6. How this differs from other MACE deployments (LAMMPS, ASE)

This is where the work is genuinely distinctive. The standard ways to run MACE —
**LAMMPS `pair_style mace`, the ASE calculator, OpenMM** — all run the actual
**PyTorch model through libtorch**, and are designed to shine on **GPUs**.

| | **Our CP2K + symmetrix** | **LAMMPS-MACE / ASE (libtorch)** |
|---|---|---|
| Engine | From-scratch C++ re-implementation of MACE inference | The real PyTorch model via **libtorch** |
| Dependency | symmetrix + sphericart, **no libtorch** | Full PyTorch C++ runtime + matching CUDA |
| Hardware | **CPU-native** (OpenMP + SIMD) | **GPU-first** (libtorch CUDA) |
| Model format | symmetrix **JSON** export | TorchScript `.pt` |
| Parallelism | Graph **replicated per rank** → OpenMP only (today) | **Spatial domain decomposition** → real multi-rank/GPU scaling |
| Ecosystem | Native **CP2K FIST**: QM/MM, `METHOD MIXED`, Ewald, samplers | LAMMPS ecosystem |

**The honest trade-off:**
- *What we give up vs LAMMPS:* GPU acceleration and genuine cross-node MPI
  scaling (their domain decomposition vs our replicated graph). For huge systems
  on GPU clusters, LAMMPS wins today.
- *What we gain:* (1) **no libtorch dependency** — huge on CPU-only HPC where
  building torch+CUDA is a nightmare; (2) a **lightweight, statically-linked
  binary**; and (3) the real prize — **direct coupling to CP2K's electrostatics
  and QM/MM machinery**, which is what makes the electrolyte science below
  possible *in one code*.

One line for your supervisor: **LAMMPS-MACE is the GPU/PyTorch route built for
scale; our symmetrix route is the CPU/torch-free route built for clean
integration into CP2K's physics ecosystem.**

---

## 7. The one trap we hit (and fixed)

MACE models can be exported in different unit conventions. During validation, the
single-point energy matched perfectly but an MD run failed to conserve energy.
Root cause: **`UNIT_FORCES` did not equal `UNIT_ENERGY / UNIT_LENGTH`** (we paired
Å lengths with bohr⁻¹ forces, off by a factor of ~1.89). The energy single-point
only uses `UNIT_ENERGY`, so the bug was **silent there** but broke F = −dE/dx in
dynamics. The tell-tale was a **dt-independent** energy drift (a correct
integrator drifts as dt²). Fix: use the self-consistent triple — for a standard
MACE model, `eV / Å / eV·Å⁻¹` — and we added the warning guard of §4.3 so nobody
else loses an afternoon to it.

---

## 8. Why this matters — the benefits and the science it unlocks

1. **Foundation-model accuracy in CP2K, on CPUs.** Users can now run MACE-MP-0
   for energies, forces, geometry optimisation, and MD — no PyTorch, no GPU
   required. That is new capability for CP2K's CPU user base.

2. **Multi-element coverage for free.** CP2K's existing water NNPs only know H and
   O. A MACE foundation model already knows Na, Cl, and most of the periodic
   table — so simulating *salt water* needs **no new training**.

3. **It enables the electrolyte science question — and this is the headline.**
   Because MACE is short-range (§1), it cannot capture the long-range
   electrostatics that dominate **ion transport in water**. But our integration
   lives inside CP2K's FIST engine, so we can run **`METHOD MIXED`: short-range
   MACE + an explicit long-range Coulomb/Ewald baseline** — exactly the setup the
   literature uses for electrolytes. This lets us pose a real, current research
   question:

   > *Does a short-range foundation MLP (MACE) reproduce NaCl-in-water ion
   > transport on its own, or only once long-range electrostatics are added back?*

   This connects directly to the literature your supervisor pointed at — **Shuwen
   Yue** (when do short-range MLPs fall short?) and **Kara Fong** (rigorous ion-
   transport coefficients). Our L3 water result already previews the answer:
   short-range MACE gets water dynamics wrong in the known way, so testing the
   long-range correction is the natural next experiment.

4. **It fits the project's CPU-optimisation theme.** A torch-free, OpenMP-parallel
   backend is consistent with the rest of the thesis (performance optimisation on
   CPU HPC), rather than a GPU detour.

---

## 9. Honest limitations (state these proactively)

- **MACE is short-range.** It does not add long-range electrostatics by itself;
  that must come from the FIST Coulomb/Ewald baseline (`METHOD MIXED`).
- **No MPI speed-up yet.** The neighbour graph is replicated per rank; only
  OpenMP parallelises the MACE evaluation today. Domain decomposition is the
  clear next optimisation.
- **No GPU.** symmetrix is CPU-only by design. For very large systems where GPUs
  dominate, LAMMPS-MACE is faster.
- **Validated for H/O water so far.** Na/Cl and other elements are *expected* to
  work (element mapping is automatic) but are **not yet validated** — that is the
  next scientific step.

---

## 10. One-paragraph script for the supervisor meeting

> "I added a MACE backend to CP2K. Instead of the usual PyTorch/GPU route, I used
> symmetrix — a CPU C++ evaluator — and integrated MACE as a FIST many-body term,
> so it composes with CP2K's electrostatics and QM/MM. I validated it end to end:
> exact energies and forces, dt²-correct energy conservation, and a stable 15 ps
> production run with the MACE-MP-0 foundation model that reproduces the model's
> known sluggish water. MACE is short-range, so on its own it won't capture ion
> transport in water — but because it's inside FIST, I can add a long-range
> Coulomb baseline via METHOD MIXED. That sets up the report-2 question: does a
> short-range foundation MLP reproduce NaCl-in-water transport, and what does
> adding long-range electrostatics fix — measured with proper transport
> coefficients à la Fong, testing the short-range-MLP limits à la Yue."

---

*Branch:* `feature/nnp-mace` · *commits:* `198408e1bd` (binding), `a83ef8eb15`
(FIST integration), `c3e371043d` (units guard). *Validated 2026-06-22.*
