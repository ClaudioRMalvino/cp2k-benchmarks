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

## 8. Cost argument (status: 3 of 5 rungs, one in flight)

From `scripts/csd3_dft_cost/analyze_dft_ladder.py` over
`results/dft_cost/dft_ladder_timings.csv`:

| rung | atoms | s/step | status |
|---|---|---|---|
| 1 | 188 | 35.502 ± 0.490 | measured |
| 2 | 376 | 75.648 ± 1.011 | measured |
| 3 | 752 | 250.244 ± 10.719 | measured |
| 4 | 1504 | — | **in flight** (job 33374342, `icelake-himem`) |
| 5 | 5064 | — | **not obtainable** — see below |

Current fit on rungs 1–3: `s/step = 2.065e-2 · N^1.409`, log-log R² = 0.983,
extrapolating to **3,414 s/step = 0.95 h per MD step** at 5064 atoms.

**Do not quote that exponent as well determined.** The local exponent between
rungs 1→2 is 1.09 and between 2→3 is 1.73 — the cost is *steepening*, so a
single power law through three points understates the extrapolation. Rung 4
is precisely the lever arm that would pin this down, which is why it is
being re-run. **Update this section once rung 4 lands.**

Rung 4 originally failed on 2026-08-09 (job 33197088): it converged the
step-0 wavefunction cleanly (E = −9065.617286 Ha) and was then OOM-killed on
the first force evaluation — 3.73 GB/rank × 76 ranks ≈ 283 GB against a
250 GB icelake node. It has been resubmitted to `icelake-himem` (same Ice
Lake silicon, same 76×1 layout, 502 GB). Rung 3 is being re-run on himem
alongside it as a cross-partition timing control.

Rung 5 (the 5064-atom production cell) stalls before the first SCF iteration
in the realspace→planewave halo exchange. Four independent attempts across
three layouts (76×1, 19×4, and a replicated realspace grid) reproduce it
identically. It is a real CP2K limit at this cell size, not a misconfigured
run, and the extrapolation from rungs 1–4 is what carries the argument.

Campaign cost comparison (2.4 M MD steps):

| | cost |
|---|---|
| On-the-fly RI-RPA (analytic) | memory-infeasible — 1.63 PB for the (ia\|P) integrals alone, 11× the entire 552-node icelake partition's RAM |
| On-the-fly revPBE-D3 AIMD | 2,276,222 node-h ≈ 173 M core-h |
| NNP, master CP2K (757bb76a80) | 2,698 node-h ≈ 205 k core-h |
| NNP, optimised CPU (PR #5295) | 170 node-h ≈ 12.9 k core-h |

Measured engine speedup on the exact campaign workload: **15.85×**
(0.2554 vs 4.0473 s/step, same node, same deck), with step-0 potential
energy identical to every printed digit between the two binaries
(−30777.866694782 Ha).

The framing that lands: the entire 84,000 core-h SL3 allocation would buy
**1,165 on-the-fly DFT steps = 0.58 ps**. The campaign produced 1,200 ps.

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

- **Rung 4 of the DFT ladder** — in flight; Section 8 needs updating with the
  measured value and the refitted exponent.
- **Yeh–Hummer finite-size corrections** on D — not applied to any number in
  this document. The reported D values are raw, single-box, 37.26 Å.
- **Two-cell finite-size check** on κ and η from the cube2 anchor segments —
  data is on disk, analysis not run.
- **Kirkwood–Buff activity derivatives** — flagged in the brief as
  order-of-magnitude only at this box size; not attempted.
- **VACF spot-check runs** — one 20 ps velocity run per concentration exists;
  not folded into this summary.
