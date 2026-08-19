# Implementation notes — Report 2 (NaCl(aq) transport)

What every piece of code in the Part-2 pipeline actually does: the simulation
setup, the algorithms behind each observable, how the statistics are formed,
and where the approximations are. Written so that any number in Report 2 can be
traced from the raw trajectory to the printed value.

Paths are relative to the repository root. Line numbers refer to the files as
committed.

---

## 0. Pipeline at a glance

```
cube_n3.xyz  ─┐
              ├─► run_equil.sh ──► 30 ps CSVR NVT + 5 snapshots 15 ps apart
model/RPA/   ─┘        │
                       ▼
              run_production.sh ──► 5 × 160 ps NVE segments (independent)
                       │              ├─ *-pos-1.xyz    positions, every 10 steps
                       │              ├─ *-1.stress     stress tensor, every step
                       │              └─ *-1.ener       energies, every 100 steps
                       ▼
              sbatch_analysis_160.sh
                       ├─ phase A  reduce_segment()          → reduced_traj_e20.npz
                       ├─ phase B  analyze_onsager_cp2k.py   → D, κ_NE, κ_Ons, t_Na
                       ├─ phase C  analyze_viscosity_cp2k.py → η
                       └─ phase D  analyze_ion_pairing_cp2k.py → g(r), w(r), K_A, n_hyd
                                        │
                                        ▼
                            make_series_figures.py → figures 8–15, summary CSV
```

Three concentrations (1 / 2 / 4 mol kg⁻¹), five segments each = 15 production
trajectories, all in the same 37.26 Å cubic cell.

**Engine.** Every production run used the `dhruv-cell-list` binary — CP2K with
PR #5295, the NNP cell-list optimisation that Part 1 of this work benchmarked
(`env_csd3.sh:14`). At the measured 0.2554 s/step against master's 4.0473, this
campaign is what the Part-1 optimisation *bought*: 2.4 M steps of
first-principles-quality dynamics that would otherwise have been unaffordable.
That is the link between the two halves of the report.

| directory  | NaCl pairs | waters | atoms | molality |
|------------|-----------:|-------:|------:|---------:|
| `cubic_1M` | 30 | 1668 | 5064 | 0.998 |
| `cubic_2m` | 58 | 1612 | 4952 | 1.997 |
| `cubic_4m` | 109 | 1510 | 4748 | 4.007 |

---

# Part A — the simulations

## A.1 The potential: a committee NNP on a Coulomb baseline

`scripts/csd3_nacl_mp2_anchor/nacl_diffusion_template.inp`

The force field is a **MIXED** CP2K calculation summing two force evaluators
(template lines 71–104):

```
&MULTIPLE_FORCE_EVALS
  FORCE_EVAL_ORDER 2 3        ! evaluator 2 = NNP, evaluator 3 = Fist
&END
&MIXED
  MIXING_TYPE GENMIX
  &GENERIC
    MIXING_FUNCTION a+b       ! E = E_NNP + E_Coulomb
  &END GENERIC
&END MIXED
```

- **Evaluator 3 (`METHOD Fist`)** is a pure electrostatic baseline: fixed
  charges (O −0.834, H +0.417, Na +1, Cl −1) evaluated with **SPME**
  (`EWALD_TYPE spme`, α = 0.35 Å⁻¹, `GMAX 120`, 6th-order splines). No
  Lennard-Jones, no bonded terms — the charges alone.
- **Evaluator 2 (`METHOD NNP`)** is an 8-member Behler–Parrinello committee
  (`train_001` … `train_008`), 320 atom-centred symmetry functions, 2 × 25
  hidden nodes, tanh cutoff.

This is the O'Neill *et al.* construction: the network is trained on the
**difference** between the reference electronic-structure energy and the
Coulomb baseline, so the long-range 1/r physics is handled analytically by SPME
and the network only has to learn the short-range remainder. That is why the
two evaluators are simply added — they are not two competing models, they are
two halves of one potential.

**Committee.** All eight members are evaluated every step. The mean is the
force used for dynamics; the spread between members is the model's own
uncertainty estimate. `&ENERGIES SILENT / &EACH MD 1` (lines 136–142) writes
the per-member energies every step so the committee disagreement can be
monitored as an extrapolation check — `nnp_check_extrapolation` inside CP2K
flags any step that leaves the training manifold.

**Why RPA.** MP2 and RPA committees share an identical architecture, so they
cost the same per step. MP2 dynamics is too sluggish to converge transport in
the available budget; RPA is closer to real water. The 1 m MP2 run is kept as a
free cross-reference through the identical analysis pipeline.

### Parallel decomposition

```
GROUP_PARTITION __NNPRANKS__ __FISTRANKS__      ! 68 8  at 76 ranks
```

CP2K's MIXED driver runs the two subsystems **concurrently on disjoint rank
sets**, so the step cost is `max(t_NNP, t_Fist)` rather than the sum. The 68/8
split was tuned by a sweep to balance the two. This is a performance choice
only — it cannot change the physics, since both evaluators see the same
configuration and their energies are summed afterwards.

## A.2 Protocol

`run_equil.sh` (stage 1) → `run_production.sh` (stage 2). Both substitute
placeholders into the one shared template with `sed`, so equilibration and
production provably use the same deck apart from the ensemble.

**Stage 1 — equilibration.** 30 ps NVT at 300 K with a **CSVR** thermostat
(canonical sampling through velocity rescaling, `TIMECON 200` fs), then five
restart snapshots written 15 ps apart (`run_equil.sh:28–31`; 210 000 steps
total). CSVR is a stochastic-velocity-rescaling thermostat that samples the
canonical ensemble correctly while perturbing the dynamics far less than
Nosé–Hoover chains.

**Stage 2 — production.** Each snapshot seeds an independent **NVE** segment.
The ensemble switch is done by deleting the thermostat block from the template
(`run_production.sh:78–81`):

```bash
-e "s|ENSEMBLE NVT|ENSEMBLE NVE|" ... | sed -e '/&THERMOSTAT/,/&END THERMOSTAT/d'
```

Transport coefficients are equilibrium properties and must be measured in a
*conservative* ensemble — a thermostat injects and removes momentum, which
directly corrupts both the stress autocorrelation (η) and the collective charge
displacement (κ). Only positions and velocities carry over from equilibration:

```
&EXT_RESTART
  RESTART_DEFAULT .FALSE.
  RESTART_POS .TRUE.
  RESTART_VEL .TRUE.
&END EXT_RESTART
```

Because the five snapshots are 15 ps apart and each segment then evolves under
its own conservative dynamics, the five segments are **statistically
independent** — this is what licenses treating them as independent samples in
§B.7.

**Numerics.** 0.5 fs timestep with natural hydrogen mass. (O'Neill *et al.*
used 1 fs with deuterated H; that variant is left commented in the template at
lines 98–102 so the paper protocol can be reproduced exactly if wanted.)
160 ps = 320 000 steps per segment.

**Sampling rates** (template lines 45–68) are set by what each analyzer needs:

| output | interval | why |
|---|---|---|
| stress tensor | **every step** (0.5 fs) | Green–Kubo needs the stress ACF resolved at the correlation time of the fluid |
| positions | every 10 steps (5 fs) | MSD only needs diffusive timescales |
| committee energies | every step | extrapolation monitoring |
| restart | every 10 000 steps (5 ps) | walltime survival |

`STRESS_TENSOR ANALYTICAL` on every evaluator: the virial is computed
analytically, not by finite differences.

## A.3 Checkpoint / resume machinery

A 160 ps segment exceeds a 12 h SL3 window, so `run_production.sh:99–122`
implements resume-from-checkpoint. Three details matter for correctness:

1. `RESTART_COUNTERS` is flipped to `T` on resume so step numbering stays
   continuous across the join (initial launches keep it off).
2. CP2K's `STEPS` is **per invocation** even with `STEP_START_VAL`, so the
   remaining budget is recomputed and the deck rewritten; the script then
   re-reads it and aborts if the cap did not take (`exit 99`). Without this a
   resumed segment would silently run the *full* length again.
3. Completion is detected by the `PROGRAM ENDED` banner **count increasing**,
   not its mere presence — CP2K appends to `prod.out` on resume, so a stale
   banner from the killed run would otherwise read as success.

The resume path duplicates up to one checkpoint interval of stress rows; the
viscosity analyzer removes them (§B.4).

## A.4 GPU track

`env_csd3_gpu.sh` swaps the module environment and binary; `run_production.sh`
injects `USE_GPU` into the `&NNP` block (lines 83–86). Same deck, same
protocol, separate tree (`runs_gpu/`) so the two never collide. CPU and GPU
agree to ~1e-9, not bit-exactly — different architectures reassociate
floating-point sums differently.

---

# Part B — the analysis

All analyzers are `numpy`-only, single-threaded, and driven by
`slurm/sbatch_analysis_160.sh`, which runs them as four phases (A serial,
then B/C/D concurrently across concentrations).

## B.1 Trajectory reduction

`analyze_onsager_cp2k.py:37–126`

Each position file is ~10.8 GB (15 files = 162 GB). Loading them into Python
is impossible, so frames are streamed through an **awk** one-liner
(`AWK`, lines 37–50) that:

- tracks position within each frame using `r = (NR-1) % (NAT+2)` — the XYZ
  block structure (count line, comment line, NAT atom lines);
- parses the step number from the comment line's `i = N` field;
- keeps every K-th frame (K = 20 → 0.1 ps sampling) and, within kept frames,
  emits only Na, Cl and O as `step elcode x y z`. Hydrogens are discarded —
  nothing downstream uses them.

The reduced array `(T, N, 3)` is cached as `reduced_traj_e20.npz`, so the
162 GB read happens **once** and phases B and D both reuse it. Frames are
checked to be strictly monotonic in step (line 121); a truncated tail frame
from a killed run is dropped (lines 117–119).

**Coordinates are unwrapped.** CP2K's XYZ trajectory is continuous (no
wrapping into the box), which is required for Einstein MSDs — a wrapped
coordinate would produce spurious jumps of one box length.

## B.2 The correlation engine (used by every observable)

`analyze_onsager_cp2k.py:53–99`

Every quantity in Part 2 is a time-correlation function, and all of them go
through one FFT routine. The naive double loop over time origins is O(T²); the
**Fast Correlation Algorithm** (Kneller; Calandrini *et al.*) does it in
O(T log T).

For the mean-square displacement, expand the square:

```
MSD(τ) = ⟨|r(t+τ) − r(t)|²⟩_t
       = ⟨|r(t+τ)|²⟩ + ⟨|r(t)|²⟩ − 2⟨r(t+τ)·r(t)⟩
       = S₁(τ) − S₂(τ)
```

- **S₁** is a pure sum of squares. It is computed from a *prefix-sum* array
  (`_prefix`, line 53) so each lag costs O(1):
  `S1 = (P[T−τ] + P[T] − P[τ]) / (T−τ)`.
- **S₂** is an autocorrelation, computed by FFT (`_crosscorr_fft`, line 59):
  transform zero-padded to 2T, multiply by the conjugate, inverse-transform.
  Zero-padding to 2T is essential — without it the FFT computes a *circular*
  correlation that wraps the end of the trajectory onto the beginning.

Both are normalised by `(T − τ)`, the number of available time origins at lag
τ. This is **multi-origin averaging**: every frame is used as a time origin,
which is what makes 160 ps of trajectory yield a usable MSD out to tens of ps.
The estimator degrades at large τ (few origins survive), which is why fits
never extend near the end of the trajectory.

`cross_msd_fft` (line 67) generalises this to *two different* signals A and B:

```
C_AB(τ) = ⟨[A(t+τ) − A(t)] · [B(t+τ) − B(t)]⟩_t
```

with `S₂ = crosscorr(a,b) + crosscorr(b,a)` — both orderings, because the cross
term is not symmetric in time. With A = B it reduces to the ordinary MSD. This
single function computes *both* the self-diffusion MSDs and the collective
Onsager correlations.

**Verification.** `--selftest` (line 91) builds correlated random walks and
compares the FFT result against the O(T²) brute-force implementation
(`cross_msd_brute`, line 81), asserting agreement to < 1e-8 for both the auto
and cross cases. This is the guard that the fast path is not subtly wrong.

## B.3 Diffusion and conductivity

`analyze_onsager_cp2k.py:129–236`

### Self-diffusion (Einstein)

For each species the per-atom MSDs are averaged (line 141), then

```
D = (1/6) · d MSD/dt
```

from a straight-line fit over a fixed window (`slope`, line 157;
`np.polyfit` degree 1, default 10–40 ps). The window matters: below ~10 ps the
motion is still ballistic/cage-rattling, and above ~40 ps the multi-origin
statistics thin out. **Window sensitivity is reported, not assumed** — the
whole Onsager analysis is run at 10–30, 10–40 and 10–60 ps
(`sbatch_analysis_160.sh`, phase B) and `window_sensitivity()` in
`make_series_figures.py:407` prints the shifts so the choice can be defended.

### Nernst–Einstein conductivity

```
κ_NE = (e²/(V k_B T)) · Σ_i z_i² N_i D_i
```

(line 201, z = ±1). This is the conductivity the solution *would* have if every
ion moved independently — it ignores all correlation between ions.

### Onsager (true) conductivity

The physically correct expression is the Einstein relation for the **collective
charge displacement**:

```
κ = (e²/(6 V k_B T)) · d/dt ⟨|Σ_i z_i Δr_i(τ)|²⟩
```

Expanding the square for a 1:1 salt gives three collective correlation
functions of the summed ionic positions `R_Na = Σ_Na r_i`, `R_Cl = Σ_Cl r_i`
(lines 145–149):

```
⟨|ΔR_Na − ΔR_Cl|²⟩ = C_NaNa + C_ClCl − 2 C_NaCl
```

so `s_sum = sNaNa + sClCl − 2·sNaCl` and `κ_Ons = pref · s_sum/6` (line 202).
The **Nernst–Einstein deviation**

```
Δ_NE = 1 − κ_Ons/κ_NE
```

measures how much ionic correlation (mostly anti-correlated Na–Cl motion, i.e.
pairing) suppresses conduction below the independent-ion estimate.

### Transference number

```
t_Na = (s_NaNa − s_NaCl) / s_sum          (line 203)
```

This is the fraction of current carried by Na⁺ in the **solvent-fixed
(barycentric)** frame. `make_series_figures.py:157` converts it to the
**Hittorf** frame — the one experiments actually measure — with the mass-flux
correction

```
t^H = t_Na + m·[M_Na(s_NaNa − s_NaCl) − M_Cl(s_ClCl − s_NaCl)] / s_sum
```

so the comparison against Smits–Duyvis (1966) Hittorf data is like-for-like.

### Channel decomposition

Lines 228–236 split κ into Na–Na, Cl–Cl and Na–Cl channels and subtract the
self (Nernst–Einstein) part of each to isolate the **distinct** correlations —
this is what shows *which* correlation is responsible for Δ_NE.

## B.4 Shear viscosity (Green–Kubo)

`analyze_viscosity_cp2k.py`

The user's question — "it's an integral, how was it solved" — is the crux of
this analyzer, so in order:

**1. Read and de-duplicate.** Stress rows are sliced to strictly increasing
step number (lines 44–50), removing the rows a resumed run rewrites.

**2. Symmetric traceless projection** (lines 51–55). The nine tensor components
are columns 2–10 of the stress file (`xx xy xz yx yy yz zx zy zz`, in bar,
converted to Pa). The tensor is symmetrised, `P_s = ½(P + Pᵀ)` — CP2K already
writes a symmetric tensor for this deck, so this is defensive rather than
load-bearing — and its isotropic part removed:

```
P°_ab = ½(P_ab + P_ba) − δ_ab·tr(P)/3
```

The trace is the hydrostatic pressure, which relates to the *bulk* viscosity,
not the shear viscosity — leaving it in would contaminate the result. This is
the **Daivis–Evans** estimator, which uses all five independent components of
the symmetric traceless tensor rather than only the three off-diagonal ones,
improving statistics by roughly √(5/3).

**3. Autocorrelation by FFT** (`acf_fft`, line 30) — same zero-padded,
`(n−τ)`-normalised machinery as §B.2, summed over all nine `ab` components
(lines 57–60).

**4. The integral.** Green–Kubo gives

```
η = V/(10 k_B T) · ∫₀^∞ Σ_ab ⟨P°_ab(0) P°_ab(t)⟩ dt
```

The integral is evaluated as a **cumulative trapezoidal rule** (lines 63–64):

```python
run = V/(10 kB T) * concatenate(([0.0], cumsum(0.5*(C[1:] + C[:-1]) * dt_s)))
```

giving the *running integral* η(t) — the value the estimate would take if the
upper limit were t. Prepending the explicit 0.0 keeps `run` aligned with the
time axis so `run[i]` is the integral up to `t[i]`.

**5. Plateau extraction.** The ∞ upper limit cannot be taken literally: the
stress ACF decays to zero within a few ps, but its *noise* does not, so the
running integral first plateaus and then random-walks away. η is taken as the
**mean of the running integral over a plateau window** (lines 85–86), 10–18 ps
for the 160 ps data. Averaging over the window rather than reading a single
point suppresses the noise without biasing the result, provided the window
really is flat.

This is the single most method-dependent number in the report, so it is
computed **twice**: `*_final` (ACF to 40 ps, plateau 10–18 ps — what 160 ps of
data buys) and `*_p0210` (ACF to 20 ps, plateau 2–10 ps — matching the earlier
80 ps round for a like-for-like comparison). The running-integral curves
themselves are stored in the `.npz` and plotted in fig 15, so a reader can see
the plateau rather than take it on trust.

## B.5 Structure: RDF, PMF, association

`analyze_ion_pairing_cp2k.py`

Reads the phase-A caches (so it costs no extra I/O) and pools all five
segments.

**Radial distribution function** (`pair_hist`, line 24). All A–B pair vectors
per frame via broadcasting, minimum-image convention applied as

```python
d -= box_l * np.round(d / box_l)
```

then histogrammed at 0.05 Å resolution out to L/2. The frame loop is chunked
(`chunk = 2e7 / (nA·nB)`, line 29) to bound peak memory. Normalisation is by
the ideal-gas expectation for each shell:

```
g(r) = h(r) / (N_frames · N_A · N_B / V · (4/3)π(r_{i+1}³ − r_i³))
```

**Potential of mean force.** Directly from the reversible-work theorem:

```
w(r) = −k_B T ln g(r)                     (line 126)
```

with `k_B T` in kcal/mol (`kt = 0.0019872041 · T`, line 95; 298.15 K to match
the MP2 anchor). Two implementation points:

- Bins with `g ≤ 0.02` are masked to NaN (line 125) — the log of a
  poorly-sampled near-empty bin is numerical noise, not a barrier.
- The curve is **tail-referenced**: the mean of w over 9 Å < r < L/2 is
  subtracted (lines 127–128), setting w(∞) = 0. Without this the absolute
  offset is arbitrary and PMFs from different concentrations cannot be
  overlaid.

The CIP minimum, the CIP/SSIP barrier and the SSIP minimum are located by
extremum search in fixed windows (lines 130–135, and again in
`make_series_figures.py:213` for the report's table conventions: CIP 2.4–3.2 Å,
barrier 3.2–4.2 Å, SSIP 4.2–5.6 Å).

**Coordination / contact-pair numbers** (`coord_number`, line 37):

```
n(r_c) = 4π ρ ∫₀^{r_c} g(r) r² dr
```

trapezoidal on the binned g. Applied with `r_c = r_b` (the pooled barrier
position) it gives n_CIP, the number of Cl⁻ in contact with each Na⁺; applied
to the ion–O RDFs out to their first minimum it gives hydration numbers
(lines 153–161). Note the deliberate split: the **cutoff comes from the pooled
curve**, but the integral is then evaluated **per segment**, so the resulting
SEM reflects genuine segment-to-segment scatter rather than jitter in the
cutoff.

**Association constant** (lines 144–151). Bjerrum-style, integrating the
Boltzmann factor of the PMF over the bound region:

```
K_A = 4π N_A ∫₀^{r_b} exp(−w(r)/k_B T) r² dr
```

(converted Å³ → L/mol by `N_A × 1e-27`). The paired fraction then follows from
the ideal mass-action equilibrium `A⁺ + B⁻ ⇌ AB`, `K_A = p/(c−p)²`, i.e. the
quadratic

```
K_A p² − (2 K_A c + 1) p + K_A c² = 0
```

solved with `np.roots`, taking the smaller (physical) root. **This is an ideal
treatment** — activity coefficients are set to 1 — so the paired fraction is a
model-internal diagnostic for comparing levels at fixed concentration, not a
prediction of the experimental association constant.

**Verification.** `--selftest` (line 53) runs the whole RDF and coordination
machinery on an ideal gas of uniformly random points, asserting g(r) = 1 to 1 %
and n(r_c) = (4/3)πρr_c³ to 2 %. That checks the shell normalisation, the
minimum-image convention and the integral together.

## B.6 Cross-level pooling and derived observables

`make_series_figures.py`

Loads the per-segment `.npz` products and forms the report's figures. Derived
quantities (`load_kappa`, line 143):

- **Molar conductivity** Λ = κ_Ons/c (line 163).
- **Ionic mobilities** in the Hittorf frame, `u_Na = t^H κ/(F c)`,
  `u_Cl = (1 − t^H) κ/(F c)` (lines 165–166) — the "Fong-style" observables,
  so the RPA series can be put on the same axes as Fong *et al.*
- **Experimental reference** (`expt_fong`, line 471) is *derived*, not measured
  as such: κ from Chambers–Stokes (1956) and t^H interpolated in log-molality
  from Smits–Duyvis (1966). The docstring says so explicitly, because a
  derived reference should never be presented as a direct measurement.

## B.7 How the expectation values and error bars are formed

Three distinct averaging layers, deliberately kept separate:

| layer | what is averaged | where |
|---|---|---|
| 1. time origins | every frame used as t = 0 in each correlation function | `_prefix` / FFT, `(T−τ)` normalisation |
| 2. particles / pooling | per-atom MSDs averaged; RDF histograms summed over all frames and segments | `segment_onsager:141`, `g_of:115` |
| 3. **segments** | the five independent NVE trajectories | every `sem()` in the pipeline |

**All quoted uncertainties are the standard error over the five segments:**

```python
sem = np.std(a, ddof=1) / np.sqrt(n)          # n = 5
```

`ddof=1` is the Bessel-corrected sample standard deviation — with n = 5 the
difference from `ddof=0` is 12 %, so it is not cosmetic.

The reason the error bar is taken at the *segment* level and not from the
scatter within a trajectory is that layers 1 and 2 average **correlated**
samples: successive time origins in one trajectory are not independent, so
their scatter badly underestimates the true uncertainty. The five segments,
started from snapshots 15 ps apart and evolved under independent conservative
dynamics, are the only genuinely independent replicates in the design — so
they are the only honest basis for an error bar. This is why the protocol
spends its budget on 5 × 160 ps rather than 1 × 800 ps.

Quantities that are non-linear functions of the primary observables (Δ_NE,
t_Na, ratios) are formed **per segment first and then averaged**, not computed
from the averaged inputs — for a non-linear function those differ, and only the
former propagates the uncertainty correctly.

## B.8 Finite-size effects — what is and is not corrected

Worth being precise about, since it is an obvious question.

The Yeh–Hummer correction for the periodic-image drag on self-diffusion,

```
D₀ = D_PBC + k_B T ξ / (6 π η L),    ξ = 2.837297
```

**is implemented and applied in the Madrid-2019 LAMMPS series**
(`scripts/cerberus_nacl_diffusion/lammps/analyze_replicas.py:18`), where
multiple box sizes were run specifically to fit the 1/L scaling.

It is **not applied to the CP2K RPA or MP2 results.** Those are single-cell
(cube3, 37.26 Å) — with one box size there is no 1/L fit to perform, and
applying the analytic correction would require assuming the η computed in the
same finite box. The cube2 partner cell that would have enabled a two-point
finite-size check was descoped from the campaign.

Consequences to state when quoting numbers:
- Absolute first-principles **D** values carry an uncorrected finite-size
  suppression (order 10–20 % for a box this size).
- **Ratios** D_Na/D_w, D_Cl/D_w largely cancel it, which is why the report
  leads with ratios when comparing levels.
- **κ and η** are collective properties with much weaker size dependence than
  D, and the Madrid comparison at matched box size is like-for-like.

## B.9 Verification built into the pipeline

Not incidental — these are the checks that make the numbers defensible:

| check | where | what it catches |
|---|---|---|
| FFT vs brute-force correlation, < 1e-8 | `analyze_onsager_cp2k.py --selftest` | a wrong fast correlation algorithm |
| ideal-gas g(r) = 1, n(r_c) exact | `analyze_ion_pairing_cp2k.py --selftest` | shell normalisation, minimum image |
| all 15 segments at step 320000, pos + stress present | `sbatch_analysis_160.sh` preflight | analysing a truncated campaign |
| 15/15 caches rebuilt, else abort | phase A gate | partial reduction silently thinning statistics |
| monotonic step check after reduction | `reduce_segment:121` | duplicated frames from resumes |
| fit window run at 10–30 / 10–40 / 10–60 ps | phase B | window-dependent conclusions |
| η at two plateau windows | phase C | plateau-choice dependence |
| committee energies every step | CP2K `&ENERGIES` | NNP extrapolation outside training data |

---

# Part C — cost and carbon accounting

**`energy_cost_ledger.py`** reads **`sacct`**, not the campaign's own
`timings.csv`, because `log_timing.sh` appends only at successful job end —
wall-killed jobs never log a row, and the CPU campaign is undercounted by ~10 %
if you trust the CSV. Energy is a **TDP-based estimate** (icelake node 700 W,
A100 550 W including host share): CSD3 has `AcctGatherEnergyType=(null)`, so no
measured per-job energy exists cluster-wide.

**`green_algorithms_footprint.py`** implements the Green Algorithms
methodology (Lannelongue *et al.*, 2021) directly rather than through their
sacct-parsing tool, because only one row of the comparison ever actually ran —
the master and 200 ns rows are *projections* from measured s/step, for which no
job records exist. Per-scenario:

```
E [kWh] = t[h] × (n_cores·TDP_core·usage + n_GPU·TDP_GPU + mem_GB·0.3725 W/GB)/1000 × PUE × PSF
carbon  = E × CI
```

Carbon intensity is fetched live from the UK National Grid API and averaged
over the actual campaign dates (122 gCO₂e/kWh over 93 daily means) rather than
using the tool's 2021-vintage default of 231 — UK grid intensity has roughly
halved since, so the default would overstate the footprint about twofold. PUE
(1.67) is an assumption with no CSD3-published source, and both it and CI are
pure multipliers, so the script prints an explicit sensitivity block instead of
presenting one number as exact.

The `ledger()` function in `make_series_figures.py:417` aggregates the campaign
into the cost table, pricing the identical step count at the measured master
s/step to get the speedup ratio. Rows with `wall_s ≤ 10` are requeue
bookkeeping and are dropped; smoke tests and the master reference row are
excluded from the campaign totals.

---

# Appendix — file map

| file | role |
|---|---|
| `nacl_diffusion_template.inp` | the one CP2K deck; MIXED NNP + SPME, all placeholders |
| `run_equil.sh` | stage 1: CSVR NVT, 5 snapshots |
| `run_production.sh` | stage 2: NVE segments, resume machinery, GPU injection |
| `slurm/sbatch_production_conc.sh` | array driver, one segment per task |
| `slurm/sbatch_analysis_160.sh` | 4-phase analysis orchestration + preflight |
| `analyze_onsager_cp2k.py` | trajectory reduction, FFT correlations, D / κ / t_Na |
| `analyze_viscosity_cp2k.py` | Green–Kubo η |
| `analyze_ion_pairing_cp2k.py` | g(r), PMF, K_A, coordination numbers |
| `make_series_figures.py` | pooling, derived observables, figures 8–15, summary |
| `energy_cost_ledger.py` | sacct-based cost/energy ledger |
| `green_algorithms_footprint.py` | carbon footprint, as-run and projected |

**Units convention throughout:** trajectories in Å and ps; conversions to SI
are applied at the point of use (`A2PS_TO_M2S = 1e-8`, `A3_TO_M3 = 1e-30`,
`BAR = 1e5`), never earlier, so intermediate arrays stay in the units the
files are written in.
