# RPA NaCl(aq) transport campaign — results and assessment

**Status: COMPLETE (2026-08-10).** All fifteen production segments finished;
every observable has been re-derived from the full 160 ps trajectories.

This document is the handoff for writing the report. It records what was
run, how the numbers were produced, what they say, and — importantly — which
claims the data will *not* support. Read the "Claims this data does not
support" section before writing anything quantitative.

Companion documents: `scripts/csd3_rpa_transport/README.md` (campaign brief
and submission sequence), `results/nacl_mp2_anchor/ANCHOR_SUMMARY.md` (the
1 m MP2 anchor), `scripts/csd3_dft_cost/README.md` (the AIMD cost ladder).

---

## 1. What was simulated

| | |
|---|---|
| System | NaCl(aq), cube3 geometry, 37.26 Å cubic box, all three molalities |
| Composition | 1 m: 30 NaCl + 1668 H₂O (5064 atoms) · 2 m: 58 + 1612 (4952) · 4 m: 109 + 1510 (4748) |
| Model | RPA 8-member committee NNP (`models/RPA/final_model/train_001..008`, O'Neill release), 320 symmetry functions, 2×25 nodes, tanh cutoff |
| Long range | MIXED: explicit SPME Coulomb alongside the short-range MLP |
| Engine | `dhruv-cell-list` (CP2K PR #5295) |
| Protocol | CSVR NVT 30 ps → 5 snapshots 15 ps apart → **5 × 160 ps NVE**, 0.5 fs, natural H |
| Output cadence | stress every step, positions every 10 steps, committee energies every step |
| Aggregate | 3 × 5 × 320,000 steps = **2.4 M MD steps = 1.2 ns** |

The 160 ps length is a **doubling of the original 5 × 80 ps design**, decided
once the CPU track had budget headroom. Section 6 covers what the doubling
did and did not buy.

### Run health (all 15 segments)

| | |
|---|---|
| Final step | 320,000 / 320,000 on every segment |
| Conserved-quantity drift | −0.0129 … +0.0109 meV/atom over 160 ps |
| Mean temperature | 296.7 … 304.9 K (target 300) |

Two segments (`cubic_1M/seg1`, `cubic_4m/seg1`) carry a few duplicate frames
where a restart overlapped. The awk reducer keeps only strictly increasing
step numbers, so they never reach the analysis; both reduce to exactly 1601
frames like everything else.

---

## 2. How the numbers were produced

One SLURM job did everything: `scripts/csd3_rpa_transport/slurm/sbatch_analysis_160.sh`
(job 33374194, `icelake-himem`, 16 cores, ran 14:45→14:53 on 2026-08-10 —
8 minutes wall).

**Everything was recomputed from the raw trajectories.** The extension
appended to the same `*-pos-1.xyz` files, so there is no "new part" spliced
onto an old result: phase A rebuilds all fifteen `reduced_traj_e20.npz`
caches with `--refresh` (the round-1 caches stop at 80 ps and would have
silently truncated everything downstream), and phases B/C/D run off those.

| phase | script | what it produces |
|---|---|---|
| A | `analyze_onsager_cp2k.reduce_segment` | reduced caches: Na/Cl/O positions every 0.1 ps, 1601 frames |
| B | `analyze_onsager_cp2k.py` | D_Na, D_Cl, D_w, κ_NE, κ_Ons, Δ_NE, t_Na — three fit windows |
| C | `analyze_viscosity_cp2k.py` | Green–Kubo η — two plateau windows |
| D | `analyze_ion_pairing_cp2k.py` | g(r), PMF, n_CIP, K_A, hydration numbers |

Sampling: positions every 10 steps × `--every 20` = 0.1 ps per analysis
frame; 1601 frames = exactly 160.0 ps per segment. Pairing reports
`frames 8005` = 5 × 1601, confirming all five replicas at full length.

Uncertainties are **SEM across the 5 independent segments** throughout.
There is no block averaging within a segment.

### Reproducing

```bash
cd ~/cp2k-benchmarks/scripts/csd3_rpa_transport/slurm
sbatch sbatch_analysis_160.sh          # ~8 min, self-contained, preflight-gated
```

The job refuses to start unless all 15 segments are at step 320,000 with pos
and stress files present.

---

## 3. Headline results (fit window 10–40 ps, η plateau 10–18 ps)

D in 10⁻⁹ m²/s, κ in S/m, η in mPa·s. ± is SEM over 5 segments.

| | D_w | D_Na | D_Cl | κ_NE | κ_Ons | Δ_NE | t_Na | η | n_CIP |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1.543±0.058 | 0.743±0.033 | 1.118±0.034 | 6.690±0.190 | 5.335±0.710 | 0.209±0.086 | 0.503±0.079 | 1.462±0.071 | 0.122±0.016 |
| 2 m | 1.404±0.047 | 0.701±0.021 | 0.939±0.025 | 11.398±0.303 | 7.879±2.615 | 0.302±0.234 | 0.442±0.110 | 1.378±0.212 | 0.231±0.019 |
| 4 m | 1.119±0.029 | 0.502±0.028 | 0.706±0.033 | 15.776±0.742 | 12.647±3.130 | 0.208±0.177 | 0.346±0.085 | 2.311±0.216 | 0.509±0.033 |

### Association / PMF (from the pooled g_NaCl(r))

| | r_peak Å | w_min kcal/mol | r_b Å | w_b kcal/mol | ΔW kcal/mol | n_CIP | K_A L/mol | paired frac | n_hyd(Na) |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 2.775 | −0.964 | 3.525 | +0.637 | 1.600 | 0.122±0.016 | 0.137 | 0.105 | 5.954 |
| 2 m | 2.775 | −0.936 | 3.575 | +0.744 | 1.680 | 0.231±0.019 | 0.130 | 0.167 | 5.839 |
| 4 m | 2.775 | −1.047 | 3.575 | +0.757 | 1.805 | 0.509±0.033 | 0.150 | 0.276 | 5.567 |

---

## 4. Assessment — what the data says

### 4.1 The model is uniformly too sluggish, but hydrodynamically self-consistent

Against the pinned experimental anchors (κ Chambers–Stokes 1956, t_Na
Smits–Duyvis 1966, η Kestin–Khalifa–Correia 1981, D_w Müller–Hertz 1996 via
Blázquez 2023 Table V):

| m | D_w RPA/exp | η RPA/exp | κ_Ons/exp | κ_NE/exp | **D_w·η RPA/exp** |
|---|---|---|---|---|---|
| 1 | 0.71× | 1.50× | 0.63× | 0.79× | **1.07×** |
| 2 | 0.70× | 1.28× | 0.54× | 0.79× | **0.89×** |
| 4 | 0.65× | 1.71× | 0.57× | 0.72× | **1.12×** |

This is the most interesting single result in the campaign. Water diffuses
~30% too slowly and the liquid is 30–70% too viscous — but **the product
D_w·η is right to within ~±12% at every concentration.** The model gets the
Stokes–Einstein relation between mobility and viscosity essentially right
while getting the absolute timescale wrong by a roughly constant factor.

That points at a global dynamical-timescale error (revPBE-D3 water is a
known-sluggish liquid, and these are classical nuclei — no NQE) rather than
at a structural failure or a broken long-range treatment. It is a much more
defensible framing than "the model underestimates diffusion", because it
identifies *what kind* of error it is.

### 4.2 Conductivity is underestimated, and correlation is the larger part of it

κ_Ons (the physically correct comparator — it includes ion–ion correlation)
lands at 54–63% of experiment. κ_NE (the uncorrelated Nernst–Einstein
estimate) lands at 72–79%. So roughly half the shortfall is inherited from
the too-slow single-ion diffusion, and the rest comes from over-correlated
ion motion, i.e. Δ_NE is too large.

Δ_NE itself (0.21, 0.30, 0.21) has error bars of 40–110% of its value and
does **not** vary monotonically with concentration in this data. Do not
build an argument on the concentration dependence of Δ_NE.

### 4.3 Ion pairing grows with concentration; K_A is roughly invariant

n_CIP rises 0.122 → 0.231 → 0.509 across 1/2/4 m, and the paired fraction
0.105 → 0.167 → 0.276. The PMF contact minimum sits at a constant 2.775 Å at
all three concentrations, with the CIP–SSIP barrier deepening modestly
(ΔW 1.60 → 1.68 → 1.81 kcal/mol).

K_A is 0.137 / 0.130 / 0.150 L/mol — essentially constant across a 4× change
in concentration. That is what an ideal-association model predicts and is a
genuine internal consistency check, worth stating: the association constant
extracted from three independent state points agrees to ~8%.

Na⁺ hydration number falls 5.95 → 5.84 → 5.57 as concentration rises, which
is the expected water-competition effect.

### 4.4 RPA was the right call over MP2

At 1 m, same cell, same pipeline, both 5 segments (MP2 anchor is 5 × 80 ps):

| | MP2 anchor | RPA 160 ps | experiment |
|---|---|---|---|
| D_w | 0.817±0.038 | 1.543±0.058 | 2.17 |
| D_Na | 0.471±0.024 | 0.743±0.033 | — |
| D_Cl | 0.605±0.056 | 1.118±0.034 | — |
| κ_Ons | 2.940±0.417 | 5.335±0.710 | 8.42 |
| η (2–10 ps plateau) | 2.129±0.113 | 1.365±0.055 | 0.9723 |
| n_CIP | 0.086 | 0.122 | — |

RPA moves every observable toward experiment — water diffusion nearly
doubles, conductivity nearly doubles, viscosity drops by a third. This
retrospectively justifies the 2026-07-30 supervisor decision to run
transport at RPA rather than MP2 ("MP2 dynamics is too sluggish for
transport"). It is a clean, quantitative vindication and worth a paragraph.

Note the η comparison uses the 2–10 ps plateau for both, because that is
what the MP2 anchor npz stores. Do not compare MP2's 2–10 ps η against RPA's
10–18 ps η.

---

## 5. Sensitivity — the windows

Fit-window sensitivity on the Onsager quantities (all three windows are in
the repo):

| m | window | D_w | κ_NE | κ_Ons | t_Na |
|---|---|---|---|---|---|
| 1 | 10–30 | 1.548±0.056 | 6.735±0.214 | 5.204±0.793 | 0.475±0.072 |
| 1 | 10–40 | 1.543±0.058 | 6.690±0.190 | 5.335±0.710 | 0.503±0.079 |
| 1 | 10–60 | 1.550±0.063 | 6.725±0.204 | 5.598±0.404 | 0.523±0.097 |
| 2 | 10–30 | 1.409±0.050 | 11.436±0.355 | 8.170±2.392 | 0.421±0.070 |
| 2 | 10–40 | 1.404±0.047 | 11.398±0.303 | 7.879±2.615 | 0.442±0.110 |
| 2 | 10–60 | 1.392±0.046 | 11.403±0.295 | 7.310±2.410 | 0.561±0.218 |
| 4 | 10–30 | 1.120±0.028 | 15.839±0.630 | 11.553±2.737 | 0.375±0.084 |
| 4 | 10–40 | 1.119±0.029 | 15.776±0.742 | 12.647±3.130 | 0.346±0.085 |
| 4 | 10–60 | 1.118±0.030 | 15.720±0.919 | 13.208±4.995 | 0.290±0.123 |

**D and κ_NE are window-insensitive** — they move by well under their SEM
across a 2× change in fit window. Quote them with confidence.

**κ_Ons and t_Na are not.** t_Na at 4 m runs 0.375 → 0.346 → 0.290 as the
window widens, a 23% drift, and its SEM grows with the window. κ_Ons at 1 m
runs 5.20 → 5.33 → 5.60. These are collective quantities built from
differences of large collective MSDs, and 5 × 160 ps is still not enough to
converge them. Always quote the window alongside them.

Viscosity plateau sensitivity (10–18 ps headline vs 2–10 ps round-1 default):

| m | η(10–18) | η(2–10) | per-segment spread at 10–18 |
|---|---|---|---|
| 1 | 1.462±0.071 | 1.365±0.055 | 1.234 … 1.599 |
| 2 | 1.378±0.212 | 1.311±0.123 | 0.900 … 1.989 |
| 4 | 2.311±0.216 | 2.056±0.113 | 1.665 … 2.813 |

The 10–18 ps window is the headline (the reasoning is recorded in
`make_series_figures.py`: the 2–10 ps running integral is still rising at
RPA 1 m, biased low by +0.31 mPa·s at 7.9σ). But note the per-segment spread
at 2 m — 0.90 to 1.99 mPa·s across five replicas is a factor of 2.2.

---

## 6. What doubling 80 → 160 ps actually bought

Like-for-like at the 10–30 ps window and the 2–10 ps η plateau (the only
settings round 1 produced):

| m | quantity | 80 ps | 160 ps | SEM change |
|---|---|---|---|---|
| 1 | κ_Ons | 7.755±0.859 | 5.204±0.793 | −8% |
| 1 | κ_NE | 7.234±0.455 | 6.735±0.214 | −53% |
| 1 | D_Na | 0.815±0.062 | 0.745±0.032 | −49% |
| 1 | η | 1.487±0.138 | 1.365±0.055 | −60% |
| 2 | κ_Ons | 9.415±2.482 | 8.170±2.392 | −4% |
| 2 | κ_NE | 11.272±0.368 | 11.436±0.355 | −3% |
| 2 | D_Na | 0.702±0.019 | 0.701±0.021 | +11% |
| 2 | η | 1.265±0.130 | 1.311±0.123 | −5% |
| 4 | κ_Ons | 14.284±5.318 | 11.553±2.737 | −49% |
| 4 | κ_NE | 16.253±0.283 | 15.839±0.630 | +123% |
| 4 | D_Na | 0.522±0.004 | 0.504±0.025 | +472% |
| 4 | η | 1.906±0.198 | 2.056±0.113 | −43% |

**The honest reading: the means are stable, the error bars are a coin flip.**

Every mean shifts by less than or about its own error bar. *That* is the
convergence evidence, and it is the claim to make.

The SEM column is not evidence of anything much. A SEM estimated from five
samples carries roughly 35% relative uncertainty of its own
(SEM/√(2(n−1))), so changes below about a factor of 2 are noise. The 4 m
D_Na entry makes this concrete: the 80 ps SEM of ±0.004 was five replicas
agreeing by luck, and ±0.025 at 160 ps is the more believable number — the
error bar got *more honest*, not worse.

---

## 7. Claims this data does not support

Please do not write any of these:

1. **"Doubling the trajectory halved the error bars."** It did not, and the
   table above is in the repo for anyone to check. Say the means are stable.
2. **"κ_Ons is converged."** Its SEM is 13%, 33% and 25% of its value at
   1/2/4 m, and it drifts with the fit window. It is the least converged
   quantity in the campaign.
3. **"η increases monotonically with concentration."** The measured order is
   1.462 (1 m), 1.378 (2 m), 2.311 (4 m) — 2 m sits *below* 1 m. The
   difference is within error bars, so it is not a real inversion, but the
   data does not demonstrate monotonicity either. Experiment *is* monotonic
   (0.9723, 1.0745, 1.3504), so state the comparison, not a trend.
4. **"t_Na reproduces experiment."** At 1 m the model gives 0.503±0.079
   against 0.374 experimental; at 2 m 0.442±0.110 vs 0.370; at 4 m
   0.346±0.085 vs 0.366. Only 4 m agrees. And t_Na is the most
   window-sensitive quantity measured.
5. **"Δ_NE decreases/increases with concentration."** 0.209, 0.302, 0.208 —
   non-monotonic with error bars of 40–110%.
6. **Any GPU transport result.** The GPU round-2 arrays were cancelled while
   still PENDING at 0:00 elapsed on 2026-08-09; that track produced no
   160 ps data. Everything here is the CPU track. Note this explicitly if
   the report's methods section still describes a dual-track design.

---

## 8. Timing, cost and energy

Everything needed to build the cost tables. All figures are measured on
CSD3 unless explicitly marked as an extrapolation or an estimate.

Regenerate with:

```bash
source ~/.fortran_env/bin/activate
python scripts/csd3_dft_cost/analyze_dft_ladder.py        # DFT ladder + fits
python scripts/csd3_rpa_transport/energy_cost_ledger.py   # consumption + kWh
```

### 8.1 Engine speedup — master vs optimised, measured

Job 32415445 (`mp2t_master_ref`) against job 32415444 (the optimised smoke).
Both on **the same physical node, cpu-q-102, all 76 cores, not overlapping**
(smoke 19:38:33–19:39:47, master 19:54:27–21:02:16, and `sacct` confirms no
other job on that node in the window).

| | optimised (PR #5295) | master 2026.1 |
|---|---|---|
| binary | `dhruv-cell-list` | pristine `757bb76a80` |
| steps | 200 | 1000 |
| s/step (ledger) | 0.2554 | **4.0473** |
| s/step (steady state, recomputed) | 0.2438 ± 0.0114 | 4.0196 ± 0.0233 |

**Speedup 15.85×** (16.49× on the steady-state recomputation).

The physics cross-check is exact — step-0 potential energy on identical
coordinates is **−30777.866694782 Ha** and conserved quantity
**−30770.651581044 Ha** in *both* binaries, to every printed digit. That is
what licenses treating them as the same calculation at different speeds.

Two caveats to state:

- the master deck is NVT/CSVR while production is NVE (otherwise identical —
  same `cube_n3.xyz`, same committee, same cell). A CSVR thermostat is a
  global velocity rescale, so the cost impact is negligible, but this is not
  a byte-identical deck and should not be claimed as one;
- the 68/8 MIXED rank split was tuned on the PR binary and never re-tuned for
  master. That biases *against* master, so 15.85× is if anything conservative.

### 8.2 On-the-fly DFT cost ladder — COMPLETE, 4 rungs

NaCl(aq), revPBE-D3, TZV2P-GTH, 1200 Ry, 1 node / 76 ranks × 1 thread,
pristine master binary. `results/dft_cost/dft_ladder_timings.csv`.

| rung | atoms | s/step | sd | partition | job |
|---|---|---|---|---|---|
| 1 | 188 | 35.502 | 0.490 | icelake | 33197085 |
| 2 | 376 | 75.648 | 1.011 | icelake | 33197086 |
| 3 | 752 | 250.244 | 10.719 | icelake | 33197087 |
| 3′ | 752 | 242.432 | 10.493 | icelake-himem | 33374344 |
| 4 | 1504 | 785.849 | 9.042 | icelake-himem | 33374342 |

Rung 3 was measured twice because rung 4 needed a different partition (see
below). The two agree to **+3.1%, t = 1.65 — not significant**, so the
partitions are timing-equivalent and rung 4 fits alongside rungs 1–3 without
an asterisk. The analyzer uses the last row per atom count, i.e. the himem
value, which is the consistent choice since rung 4 is also himem.

**The local exponent settles rather than steepening:**

| range | exponent |
|---|---|
| 188 → 376 | 1.091 |
| 376 → 752 | 1.680 |
| 752 → 1504 | 1.697 |

So rung 1 sits in a shallow small-N regime and 376–1504 atoms is a clean
power law. Two defensible fits:

| fit | exponent | R² | s/step at 5064 | h/step |
|---|---|---|---|---|
| rungs 1–4 (analyzer default) | 1.509 | 0.9906 | 4,506 | 1.25 |
| **rungs 2–4 (asymptotic regime)** | **1.688** | **0.99999** | **6,092** | **1.69** |

Prefer the rungs 2–4 fit if you want the physically meaningful exponent — R²
= 0.99999 across a 4× range in N, against a rungs 1–4 fit that averages two
regimes. Quote rungs 1–4 if you prefer the conservative number. Either way
the earlier 3-rung value (1.409, 3,414 s/step) is superseded and was an
**underestimate**, exactly as flagged.

**Rung 4 needed `icelake-himem`.** Its first attempt (33197088) converged the
step-0 wavefunction cleanly (E = −9065.617286 Ha) then OOM-killed on the
first force evaluation: 3.73 GB/rank × 76 ranks ≈ 283 GB against a 250 GB
icelake node, short by 13%. himem is the same Ice Lake silicon with 502 GB,
so the 76×1 layout the ladder depends on was preserved and only the DRAM
ceiling moved.

**Rung 5 (5064 atoms) is not obtainable.** It stalls before the first SCF
iteration in the realspace→planewave halo exchange. Four independent attempts
across three layouts (76×1, 19×4, replicated realspace grid) reproduce it
identically — a real CP2K limit at this cell size, not a misconfigured run.
The extrapolation from rungs 1–4 carries the argument.

### 8.3 Campaign cost comparison (2.4 M MD steps = 1.2 ns aggregate)

| method | node-h | core-h |
|---|---|---|
| On-the-fly RI-RPA (analytic) | infeasible | — |
| On-the-fly revPBE-D3 AIMD (rungs 1–4 fit) | 3,003,762 | 228,285,932 |
| NNP, master CP2K `757bb76a80` | 2,698.2 | 205,063 |
| NNP, optimised CPU (PR #5295) | 170.3 | 12,940 |
| NNP, GPU (PR #5295 + `nnp_gpu`) | — | 97.8 GPU-h |

RI-RPA is **memory**-infeasible, not merely slow: 1.63 PB for the (ia\|P)
integrals alone at 5064 atoms, 11× the total RAM of the entire 552-node
icelake partition (148 TB). That is a stronger statement than a time estimate
and should be made in those terms.

Two framings, both measured:

- the entire 84,000 core-h SL3 allocation buys **883 on-the-fly DFT steps =
  0.44 ps**; the campaign produced 1,200 ps, **2,718× more trajectory**;
- on master, one 12 h SL3 window holds 10,673 steps, so a single 160 ps
  segment needs **30 chained windows** and the campaign needs **450 job
  windows**. That is an operational impossibility rather than a budget
  overrun, and it is harder to argue with than a core-hour total.

### 8.4 What the campaign actually consumed

From `sacct` — the accounting database SLURM bills from — not from
`timings.csv`. The ledger written by `log_timing.sh` appends its row at the
*end* of a job, so a job killed by its wall limit never logs one and the
5-element extend arrays logged a single row between them; aggregating it
undercounts the CPU campaign by ~10%. Rows with elapsed ≤ 10 s are excluded
as requeue bookkeeping.

| group | jobs | node-h | core-h | GPU-h | kWh (est.) |
|---|---|---|---|---|---|
| CPU campaign | 43 | 389.7 | 29,614.0 | — | 272.8 |
| GPU campaign | 26 | 173.2 | 5,542.1 | 173.2 | 95.3 |
| DFT cost ladder | 13 | 19.3 | 1,441.4 | — | 13.5 |
| analysis | 1 | 0.1 | 2.1 | — | 0.1 |
| **total** | | **582.3** | **36,599.6** | **173.2** | **381.6** |

GPU jobs took `gres/gpu=1` (verified across all 43), so GPU-hours equal
elapsed hours.

**Do not mix these two quantities — they differ by 2.3×:**

| | value | correct use |
|---|---|---|
| idealised production-only | 170.3 node-h | the DFT comparison, because the DFT figure is idealised the same way |
| measured consumption | 389.7 node-h | any £ or kWh claim, because it is what the machine actually spent |

The gap is equilibration, insurance arrays, and re-runs after wall-clock
timeouts. Both are honest; using the wrong one in the wrong sentence is not.

### 8.5 Energy — TDP estimate, NOT measured

**Measured per-job energy is unavailable on CSD3.** `scontrol show config`
reports `AcctGatherEnergyType = (null)`, and `ConsumedEnergyRaw` is 0 for
every job — not only ours: a sweep of *all* cluster jobs over two days found
none with a nonzero value. This is a site configuration, not a permissions or
query problem, and the report should state it as such rather than leaving the
absence unexplained.

The kWh column above is therefore an estimate on these assumptions, which
must be quoted alongside any energy number:

| | |
|---|---|
| icelake node | **700 W** — 2 × Xeon Platinum 8368Q (38 cores, 270 W TDP each) = 540 W, plus DRAM, board, NIC, fans; midpoint of the 650–750 W band |
| A100 (ampere) | **550 W** per GPU including a quarter-node host share |
| PUE | **not applied** — the figures are IT load only |

Hardware pinned from the machine itself via `scontrol show node` / `lscpu` on
2026-08-10. To upgrade this to an authoritative figure, ask
support@hpc.cam.ac.uk for the standard per-node power draw and the
data-centre PUE — RCS hold both for sustainability reporting. Then re-run
with `--pue <value>`.

### 8.6 Adding £ — the arithmetic

The CSD3 rate card is behind Raven, so rates are inputs rather than repo
constants. Take the SL3 icelake £/core-hour and the ampere £/GPU-hour from
the internal service-charges page on hpc.cam.ac.uk, note the access date,
then:

```bash
python scripts/csd3_rpa_transport/energy_cost_ledger.py \
  --gbp-core-hour <rate> --gbp-gpu-hour <rate> \
  --accessed 2026-08-XX --pue <from RCS>
```

which multiplies the correct column per group (core-h for CPU groups, GPU-h
for the ampere group) and prints the £ column and the citation line. Cite as:
*University of Cambridge Research Computing Services, internal service
charges page, accessed 2026-08-XX.* Assessors are Cambridge members and can
verify it.

The three multiplications, if doing it by hand:

| | |
|---|---|
| CPU campaign | 29,614.0 core-h × £/core-h |
| GPU campaign | 173.2 GPU-h × £/GPU-h |
| master projection | 205,063 core-h × £/core-h |

### 8.7 Sampling context — we are short, and should say so

O'Neill's reference simulations ran **200 ns** (per user, 2026-08-10 —
*needs a citation before it goes in the report*). Ours is 160 ps per replica,
1.2 ns aggregate. Matching 200 ns at our cell size and engine would cost:

| | trajectory | steps | node-h at our cell/engine |
|---|---|---|---|
| one segment | 0.16 ns | 320,000 | 22.7 |
| our campaign | 1.20 ns | 2,400,000 | 170.3 |
| O'Neill-equivalent | 200 ns | 400,000,000 | **28,378** |

So their sampling is **167× our production cost** at our system size — which
is presumably why their cells are smaller. This is worth stating plainly
rather than hiding: it explains why κ_Ons and t_Na are unconverged here
(Section 5) while structural quantities like the PMF and K_A are solid, and
it frames our contribution as a large-cell, short-trajectory study rather
than pretending to compete on sampling.

### 8.8 Carbon footprint — Green Algorithms methodology

Produced by `scripts/csd3_rpa_transport/green_algorithms_footprint.py`, which
implements Lannelongue, Grealey & Inouye, *Adv. Sci.* **8** (2021) 2100707 —
the methodology the CSD3 docs point to at
`cambridge-ceu.github.io/csd3/systems/GreenAlgorithms.html`. Their configured
copy of the `GreenAlgorithms4HPC` tool sits in another project's RDS space and
is not readable; our own copy is installed and configured for CSD3 at
`/rds/user/crm98/hpc-work/GreenAlgorithms4HPC` if a cross-check is wanted, but
it derives node-hours by parsing `sacct` and **only the optimised-CPU row ever
ran** — `master` and both 200 ns rows are projections with no job records — so
the numbers below are computed from the Table 9.4 hour canon directly.

    E [kWh] = t[h] x (n_cores x TDP_core x usage + n_GPU x TDP_GPU
                      + mem_GB x 0.3725 W/GB) / 1000 x PUE x PSF

**Inputs, all taken from the machines rather than spec sheets.** icelake:
`AllocTRES` of a production job is `cpu=76,mem=256120M,node=1`; Xeon 8368Q is
270 W per 38-core socket, so 7.105 W/core → **633.2 W/node**. ampere:
`cpu=32,gres/gpu=1,mem=250G` — a *quarter node* per GPU; EPYC 7763 is
4.375 W/core and an A100-SXM4-80GB is 400 W (not the 250–300 W PCIe figure)
→ **633.1 W/GPU**.

> **Correction to Table 9.4.** The GPU row is costed there at "3 SPME cores".
> That is right for the £ column, since CSD3 bills GPU-hours, but wrong for
> energy: the jobs *reserved* 32 cores and 250 GB, which no one else could use
> while they held them. Carbon accounting charges reserved, not used.

**Carbon intensity is measured, not assumed.** The UK National Grid Carbon
Intensity API (`api.carbonintensity.org.uk`, run by National Grid ESO with EDF
Europe, Oxford and WWF), endpoint `/intensity/stats/{from}/{to}/24`, gives
**122.0 gCO2e/kWh** as the mean of 93 daily means over the campaign window
2026-05-10 → 2026-08-11 (daily range 57–209). The Green Algorithms UK default
of 231.12 is a 2021 vintage and would nearly double the answer — grid intensity
has roughly halved since. Cite the API and the window.

> This is the **national** figure. The same API serves regional data and
> Cambridge (CB3) is region 10, "East England"/UKPN East, whose generation mix
> differs. Recomputing regionally is a two-minute job if the report wants the
> tighter number.

**PUE = 1.67 — this is the weakest input in the whole calculation and must be
labelled as an assumption, not a value.** CSD3 publishes no PUE; RCS have not
replied. 1.67 is the number in the Green Algorithms *configuration template*,
where it appears as `PUE: <1.67>` — angle brackets, i.e. a placeholder to be
replaced with the real facility figure, not a recommended default. Web searches
suggest other work on CSD3 falls back to ~1.67 for the same reason, but **no
such publication has been verified**, so do not cite one. Modern facilities
target ~1.2. Since PUE is a pure multiplier, the honest presentation is to
state it in the caption and give the rescaling factor.

| | hours | kWh | kgCO2e | tree-months | car-km |
|---|---:|---:|---:|---:|---:|
| **as run, 2.4 ns** | | | | | |
| master (projected) | 5,396 node-h | 5,706 | 695.8 | 759 | 3,976 |
| optimised CPU (ran) | 340 node-h | 360 | 43.8 | 48 | 251 |
| optimised GPU (projected) | 227 GPU-h | 240 | 29.3 | 32 | 167 |
| **projected, 200 ns** | | | | | |
| master | 450,000 node-h | 475,826 | 58,030 | 63,283 | 331,602 |
| optimised CPU | 28,000 node-h | 29,607 | 3,611 | 3,938 | 20,633 |
| optimised GPU | 19,000 GPU-h | 20,089 | 2,450 | 2,672 | 14,000 |

**What the optimisation avoided.** At the campaign as run, 5,346 kWh /
0.65 tCO2e (CPU) and 5,466 kWh / 0.67 tCO2e (GPU) against master — about 13
Paris–London flights each. At O'Neill-scale sampling it is 54.4 and
55.6 tCO2e, roughly 1,100 flights.

**The finding that is easy to miss:** an icelake node and a quarter ampere node
draw within 0.1 W of each other (633.2 vs 633.1). The GPU track's advantage is
therefore *entirely* that it needs fewer hours — 227 vs 340, a factor 1.5 — and
not the ~2× that comparing 550 W/A100 against 700 W/node would suggest. Any
claim that the GPU port is dramatically greener per unit time is unsupported.

**Sensitivity.** Both PUE and CI are pure multipliers on a fixed,
job-dependent core, so state them and let the reader rescale: PUE 1.67 → 1.2
multiplies everything by 0.72; CI 122.0 → 231.12 multiplies carbon by 1.90.
Equilibration (30 ps × 3) adds ~4% to the as-run rows for an all-in number.

**Caveat that must survive into the report.** These are *modelled* figures. CSD3
has `AcctGatherEnergyType = (null)`, so `ConsumedEnergyRaw` is zero for every
job cluster-wide and no measured per-job energy exists on this machine for
anyone. Label the column "estimated".

---

## 9. Files to use

All under `results/rpa_transport/`. `{C}` is `` (1 m), `_2m`, `_4m`.

| file | contents |
|---|---|
| `rpa_cpu_cube3{C}_5x160_w1040_kappa.{csv,npz}` | **headline** D/κ/t_Na, 10–40 ps |
| `rpa_cpu_cube3{C}_5x160_w1030_kappa.{csv,npz}` | 10–30 ps, and the 80 ps comparison point |
| `rpa_cpu_cube3{C}_5x160_w1060_kappa.{csv,npz}` | 10–60 ps, new at 160 ps |
| `rpa_cpu_cube3{C}_5x160_final_eta.npz` | **headline** η, ACF to 40 ps, plateau 10–18 ps |
| `rpa_cpu_cube3{C}_5x160_p0210_eta.npz` | η on round-1's 2–10 ps plateau |
| `rpa_cpu_cube3{P}_5x160_pairing.csv` | one summary row per concentration |
| `rpa_cpu_cube3{P}_5x160_rdf_ions.csv` | r, g_NaCl, w_NaCl, g_NaO, g_ClO |
| `MANIFEST_5x160.txt` | generated inventory with job id and provenance |

Note the pairing labels use `_1M` where the kappa/eta labels use nothing —
`rpa_cpu_cube3_1M_5x160_pairing.csv`. That asymmetry is inherited from round
1 and is preserved deliberately so both rounds share a naming scheme.

Round-1 80 ps products remain untouched under `*_5x80_*` and `*_final_*`
labels. They are the input to the Section 6 comparison; they are **not**
inputs to anything in Sections 3–5.

### Figures

`scripts/csd3_rpa_transport/make_series_figures.py` has been repointed at
these labels (`SERIES`, `WINDOWS`, `PAIRING`) — it previously referenced the
abandoned GPU 80 ps files. Its `ETA_WINDOW` is already (10, 18) ps, matching
the headline η. It writes fig9–fig14 plus
`results/transport_comparison/transport_series_summary.csv`, and overlays the
Madrid-2019 classical curves and the MP2 anchor.

```bash
source ~/.fortran_env/bin/activate
python scripts/csd3_rpa_transport/make_series_figures.py
```

---

## 10. Still open

- **£ rates** — Section 8.6; needs the Raven-gated rate card and an access date.
- **Authoritative per-node power and PUE** — Section 8.5; one email to
  support@hpc.cam.ac.uk upgrades the kWh column from a TDP estimate.
- **Citation for O'Neill's 200 ns** — Section 8.7 currently attributes it to
  the user, which is not citable.
- **Yeh–Hummer finite-size corrections** on D — not applied to any number in
  this document. The reported D values are raw, single-box, 37.26 Å.
- **Two-cell finite-size check** on κ and η from the cube2 anchor segments —
  data is on disk, analysis not run.
- **Kirkwood–Buff activity derivatives** — flagged in the brief as
  order-of-magnitude only at this box size; not attempted.
- **VACF spot-check runs** — one 20 ps velocity run per concentration exists;
  not folded into this summary.
