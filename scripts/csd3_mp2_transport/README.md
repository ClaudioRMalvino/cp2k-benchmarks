# CSD3 MP2 transport campaign (2026-07-30 brief)

Full MP2-level transport + thermodynamics series for NaCl(aq) on the
optimised branch, plus the measured optimised-vs-master cost proof.
**Scope: 1 / 2 / 4 mol/kg only** (0.25/0.5 cut 2026-07-30; Madrid's
fractional points are simply dropped from the joint analysis).

All cells are cube3 geometry (37.26 Å cubic):

| conc dir  | pairs | waters | atoms | mol/kg |
|-----------|------:|-------:|------:|-------:|
| cubic_1M  |    30 |   1668 |  5064 |  0.998 |
| cubic_2m  |    58 |   1612 |  4952 |  1.997 |
| cubic_4m  |   109 |   1510 |  4748 |  4.007 |

Engines (staged at `/rds/user/crm98/hpc-work/cp2k_binaries/csd3/`):
- **campaign engine** = `dhruv-cell-list` (PR #5295, the branch this campaign
  showcases — anchor-validated on this exact deck: smoke 31514146, NNP
  regtest 10/10, cube3 0.251 s/step at the 68/8 split)
- **master** = pristine 2026.1 (757bb76a80) — reference timing only
- Claudio's own optimisation builds (chebyshev etc.) are NOT used here

Model (supervisor decision 2026-07-30): **RPA 8-member committee**
(`models/RPA/final_model/train_001..008`, staged from the O'Neill release).
Identical architecture to MP2 (320 symfunctions, 2x25 nodes, tanh cutoff)
=> identical s/step; chosen because MP2 dynamics is too sluggish for
transport while RPA is closer to real water. All runs live under
`runs/RPA/`. The MP2 anchor (runs/MP2/cubic_1M) stays as a free 1 m
cross-reference: MP2-vs-RPA model comparison + two-cell finite-size check
on kappa/eta from data already on disk.

Protocol per concentration: CSVR NVT 30 ps + 5 snapshots 15 ps apart →
**5 × 80 ps NVE segments** (0.5 fs, natural H), stress every step, positions
every 10 steps, committee energies every step; one 20 ps velocity run
(VACF spot-check). 80 ps (160k steps ≈ 11.4 h at the smoke-measured
0.2554 s/step) fits a single 12 h SL3 window — the resume arrays become
pure insurance. Statistical basis: the MP2 anchor validated κ at 14% SEM
with 5×160 ps at MP2's sluggish dynamics; RPA decorrelates ~2× faster, so
5×80 ps is expected to match — verified by the pilot before fan-out, and
segments are checkpoint-extendable if it doesn't hold. Every concentration (including 1 m) is equilibrated
fresh under `runs/RPA/` — the RPA model has no prior runs, so no seed
tricks or subdir separation are needed (the SEED_OVERRIDE / PROD_SUBDIR
machinery remains available in the run scripts).

Cost bookkeeping: every job appends to
`~/cp2k-benchmarks/results/mp2_transport/timings.csv` via `log_timing.sh`.

## Dual-track strategy (user decision 2026-07-30)

**GPU is the primary track for the report's findings** (ampere queue is much
faster than SL3-CPU); the CPU track keeps running and, if it completes in
time, the cost table becomes the three-way ladder master / optimised-CPU /
GPU. Both tracks run the identical deck, protocol (5x80 ps) and analysis;
the GPU track lives in `runs_gpu/` (own equil, own snapshots) so the trees
never collide. GPU-vs-CPU expected agreement ~1e-9 (not bit-exact across
architectures) — the pilot D/kappa cross-check between tracks doubles as the
hardware validation figure.

GPU plumbing: `env_csd3_gpu.sh` (ampere modules + GPU-binary libs, selected
via ENV_SCRIPT), `USE_GPU_OVERRIDE=ON` injected into the deck by
run_equil/run_production, 1 A100 + 3 CPU SPME ranks (GROUP_PARTITION 1 3),
binary from build 32427234 (CUDA-link + nnp_gpu-symbol gates PASSED).

GPU submission sequence (AFTER smoke 32428134 passes — check its s/step and
virial ON-vs-OFF diff first; trim --time to fit the measured speed):

```bash
cd ~/cp2k-benchmarks/scripts/csd3_mp2_transport/slurm
GPU=NIKIFORAKIS-CSC-FUNDS-SL3-GPU     # lowercase it if sbatch rejects
e=$(sbatch --parsable -A $GPU sbatch_gpu_equil_conc.sh cubic_1M)
er=$(sbatch --parsable -A $GPU --dependency=afterany:$e sbatch_gpu_equil_conc.sh cubic_1M)
p=$(sbatch --parsable -A $GPU --array=1-5 --dependency=afterok:$er sbatch_gpu_production_conc.sh cubic_1M)
sbatch -A $GPU --array=1-5 --dependency=afterok:$er,afterany:$p sbatch_gpu_production_conc.sh cubic_1M
# 2m/4m equil may start immediately; their production waits for the pilot
# kappa verdict (from whichever track's pilot lands first).
```

The VACF velocity run stays CPU-only (one cheap job per concentration,
already queued for 1 m) — no GPU duplication needed.

## 1 m strategy (updated for RPA, 2026-07-30)

The RPA pilot at 1 m is a fresh full run (equil + 5 segments + vacf) and
carries the κ/η error-bar verdict that gates the 2m/4m fan-out. The MP2
anchor is no longer a regression reference (different model, different
observables expected) — it becomes the 1 m MP2-vs-RPA comparison and the
two-cell (cube2/cube3) finite-size check on collective properties, both
free from data already on disk. The supervisor's final figure is FFMD
(Madrid) vs MLIP-MD (RPA) over concentration for the easiest-to-converge
properties (D's, PMF/K_A/ΔG_A, η), with κ/Δ_NE/t_Na as the zero-extra-cost
deeper analysis, plus the cost-savings breakdown.

## Science additions vs the brief (analysis-side only, 2026-07-30)

Mapped against the O'Neill outlook paragraph:
- **ΔG_A(c) = −kT ln K_A(c)** quoted explicitly next to K_A (their "standard
  ion pair association energies" first step — answered literally).
- **Finite-size check on collective properties, free**: κ and η extracted
  from the anchor *cube2* (24.84 Å) segments as well — two-cell comparison
  at 1 m addresses their "long-range/finite-size" open question for κ, not
  just for D (Yeh–Hummer) and η (YH-slope cross-check).
- **Kirkwood–Buff stretch item (optional, caveated)**: activity-coefficient
  derivatives from g_ij(r) at 3 concentrations; KB integrals converge poorly
  at 37 Å with 30–115 pairs, so order-of-magnitude only — never a headline.
- Model-side long-range answer stated in the report: MIXED carries explicit
  SPME Coulomb, the remedy for the short-range-MLP limitation they raise.

## Submission sequence

```bash
cd ~/cp2k-benchmarks/scripts/csd3_mp2_transport/slurm

# 0. smoke (SL2, ~30 min) - gates everything on the optimised engine
smoke=$(sbatch --parsable sbatch_smoke_cube3.sh)

# 2. master reference (SL3, independent of smoke)
sbatch sbatch_master_ref.sh

# 3. equil, all three concentrations (+resume insurance; resumes also
#    require afterok on the smoke so a failed gate cancels the whole tree)
e1=$(sbatch --parsable --dependency=afterok:$smoke sbatch_equil_conc.sh cubic_1M)
e1r=$(sbatch --parsable --dependency=afterok:$smoke,afterany:$e1 sbatch_equil_conc.sh cubic_1M)
e2=$(sbatch --parsable --dependency=afterok:$smoke sbatch_equil_conc.sh cubic_2m)
sbatch --dependency=afterok:$smoke,afterany:$e2 sbatch_equil_conc.sh cubic_2m
e4=$(sbatch --parsable --dependency=afterok:$smoke sbatch_equil_conc.sh cubic_4m)
sbatch --dependency=afterok:$smoke,afterany:$e4 sbatch_equil_conc.sh cubic_4m

# 1. pilot production 1 m (validates kappa error bars before fan-out)
p=$(sbatch --parsable --array=1-5 --dependency=afterok:$e1r sbatch_production_conc.sh cubic_1M)
sbatch --array=1-5 --dependency=afterok:$e1r,afterany:$p sbatch_production_conc.sh cubic_1M
sbatch --dependency=afterok:$e1r sbatch_vacf.sh cubic_1M

# 4. AFTER pilot analysis confirms the error bars: fan-out (same pattern)
# p2=$(sbatch --parsable --array=1-5 sbatch_production_conc.sh cubic_2m)
# sbatch --array=1-5 --dependency=afterany:$p2 sbatch_production_conc.sh cubic_2m
# ... and cubic_4m; then sbatch_vacf.sh cubic_2m / cubic_4m
```

QC per segment before it enters analysis:

```bash
source ~/.fortran_env/bin/activate
python ~/cp2k-benchmarks/scripts/csd3_nacl_mp2_anchor/check_run.py \
  /rds/user/crm98/hpc-work/nacl_mp2_anchor/runs/MP2/cubic_*/production*/cube3/seg*
```

Cross-binary check: the smoke and master_ref jobs both print their step-0
potential energy on the same coordinates — PR #5295's cell-list NNP vs
pristine master should agree to ~1e-10 Ha.

Walltime at the smoke-measured 0.2554 s/step (cube3, 68/8 split): equil
≈ 14.9 h (2 windows — resume pre-chained), production segment 80 ps ≈
11.4 h (single window; resume array = insurance only), vacf ≈ 2.9 h.
Campaign total ≈ 230 node-h ≈ 17k core-h (down from ~30k at 160 ps).

## Cost proof — MEASURED (2026-07-30, jobs 32415444 / 32415445)

Same 1 m RPA cube3 deck, same node (cpu-q-102), 76 ranks, 68/8 split:

| engine | s/step | 80 ps segment | campaign (~3.15M steps) |
|---|---|---|---|
| PR #5295 (dhruv-cell-list) | **0.2554** | 11.4 h — 1 window | ≈224 node-h ≈ **17k core-h** (20% of SL3) |
| master 2026.1 (757bb76a80) | **4.0473** | ≈180 h — 15 chained windows | ≈3,540 node-h ≈ **269k core-h** (3.2× the entire allocation) |

**Speedup 15.85× on the exact campaign workload** — the campaign is
infeasible on master with the available budget. Physics check: step-0
potential energy identical to every printed digit (−30777.866694782 Ha)
between the two binaries. Caveat: the 68/8 MIXED split was tuned on the
PR binary and not re-tuned for master (minor — the NNP term dominates
master's cost). Raw rows in results/mp2_transport/timings.csv.
