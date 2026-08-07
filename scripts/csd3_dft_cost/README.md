# On-the-fly (AIMD) cost ladder

Measures what this campaign would have cost with **on-the-fly** electronic
structure — i.e. Born–Oppenheimer AIMD, solving the SCF at every MD step —
instead of running MD on a neural network potential fitted once, offline.

This supplies the top rungs of the report's cost ladder. The three NNP rungs
below them are already measured (master CP2K 4.0473 s/step, optimised CPU
0.2554, GPU 0.1467, all on the same 5064-atom cell).

## Why a ladder at all, when rung 5 is the answer

The reference electronic structure was only ever run on **160–192 atom** cells
(that is what the training set is made of). The production cells are
**5064 / 4952 / 4748 atoms** — about 30× larger. So the cost gap has two
multipliers: cost per step at fixed size, *and* the scaling penalty on the way
up to production size. The NNP is O(N) and pays only the first; DFT and RPA
pay both.

Rung 5 measures the production cell directly, so the headline number needs no
extrapolation. Rungs 1–4 still earn their place: they measure the exponent
`s/step ∝ N^α`, which (a) corroborates rung 5 if the fit extrapolates onto it,
(b) shows *why* the gap widens with system size rather than just asserting it,
and (c) gives the cost at the other concentrations without rerunning them.

## What is run

`dft_ladder.inp.template` is the O'Neill revPBE-D3 training deck
(`potentials/oneill-nacl-water-2024/data-sets/revPBE-D3/revpbe-d3_nacl.inp`)
reproduced verbatim — TZV2P-GTH / GTH-PBE, 1200 Ry cutoff, NGRIDS 5, OT-DIIS
with FULL_ALL preconditioner, EPS_SCF 5e-7, revPBE + D3 (long-range corrected),
`EXTRAPOLATION USE_GUESS`, deuterium H mass — with a cell, an NVE MD block and
`MULTIPLE_UNIT_CELL` added.

System size is varied **only** through `MULTIPLE_UNIT_CELL`, replicating one
equilibrated 64-water box (192 atoms, 12.42 Å cubic, ambient density, frame 1
of their training trajectory):

The system is **NaCl(aq) throughout** — the campaign's own cell lineage, the
same one used for the MP2 anchor, the Madrid comparison and the RPA
production runs. There is no pure-water proxy anywhere in the ladder.

| rung | cell source | replication | cell (Å) | atoms | composition | mol/kg | MD steps |
|---|---|---|---|---|---|---|---|
| 1 | `cube_n1` | `1 1 1` | 12.42 cubic | 188 | 1 NaCl + 62 H₂O | 0.895 | 11 |
| 2 | `cube_n1` | `2 1 1` | 24.84×12.42×12.42 | 376 | 2 NaCl + 124 H₂O | 0.895 | 11 |
| 3 | `cube_n1` | `2 2 1` | 24.84×24.84×12.42 | 752 | 4 NaCl + 248 H₂O | 0.895 | 11 |
| 4 | `cube_n1` | `2 2 2` | 24.84 cubic | 1504 | 8 NaCl + 496 H₂O | 0.895 | 11 |
| **5** | **`cube_n3` (production)** | `1 1 1` | **37.26 cubic** | **5064** | **30 NaCl + 1668 H₂O** | **1.000** | 4 |

**Rung 5 is the headline: it is the production cell itself.** Not a
size-matched proxy — the very `cube_n3.xyz` coordinates the RPA campaign
started from, at 37.26 Å and 1 mol/kg. So the reported on-the-fly DFT cost is
*measured on the system under study*, exactly as the NNP rungs are, with no
extrapolation and no system substitution.

Rungs 1–4 exist only to supply a scaling exponent as an independent
cross-check on rung 5. They replicate `cube_n1` (built by
`scripts/cerberus_nacl_diffusion/make_base_cell.py`, the same builder that
produced `cube_n2`/`cube_n3`), which holds molality exactly constant at
0.895 mol/kg across the fit — replication cannot change molality — so the
fitted slope is pure size-scaling with no composition drift. The 0.895 vs
1.000 mol/kg difference between the fitted rungs and rung 5 shifts the
electron count by ~0.3% and is irrelevant to timing.

Step 1 pays for the ATOMIC-guess SCF and is discarded; the reported
`s_per_step` is the mean (± sd) over the remaining steady-state steps, taken
from the `.ener` UsedTime column — the same timing methodology used for every
NNP number in this campaign. Rung 5 gets only 4 steps because each one costs
hours. The `.ener` file is written incrementally, so a rung that hits its wall
limit still leaves usable per-step timings on disk; a TIMEOUT loses the CSV
row, not the measurement.

Binary: the **pristine upstream master build** (757bb76a80), the same one that
produced the master NNP reference, on the same icelake 1-node / 76-rank
layout. DFT and NNP numbers are therefore directly comparable.

## Running it

```bash
sbatch --time=01:00:00 sbatch_dft_ladder.sh 1     # 192 atoms
sbatch --time=02:00:00 sbatch_dft_ladder.sh 2     # 384
sbatch --time=04:00:00 sbatch_dft_ladder.sh 3     # 768
sbatch --time=12:00:00 sbatch_dft_ladder.sh 4     # 1536
sbatch --time=12:00:00 sbatch_dft_ladder.sh 5     # 5184 = production cell

source ~/.fortran_env/bin/activate
python analyze_dft_ladder.py
```

Rung 5 is the memory-risky one: `FULL_ALL` builds a dense matrix over ~69,000
basis functions. If it OOMs on one node, rerun it across 2–4 nodes and report
node-hours with a parallel-efficiency caveat — but prefer the single-node
result, since that is what makes it directly comparable to the NNP rungs.

Results land in `results/dft_cost/dft_ladder_timings.csv` (raw rungs) and
`dft_cost_summary.csv` (fitted extrapolation + full cost table).

## RPA is analytic, not measured

On-the-fly RPA is not benchmarked because it cannot be run at production size
at all. O'Neill's RPA deck is RI-RPA on PBE, cc-TZ + RI auxiliary basis, 20
minimax quadrature points, truncated-Coulomb exact exchange, on 168-atom cells.

Going 168 → 5064 atoms is a factor 30.1 in size:

- **time** scales O(N⁴) → 30.1⁴ ≈ **8.3×10⁵×** per force evaluation
- **memory** scales O(N³) → 30.1³ ≈ **2.7×10⁴×** for the `(ia|P)` three-centre
  integrals, which are ~60 GB at the 168-atom reference
  (N_occ 224 × N_virt 3024 × N_aux ~11000 doubles)

That is **~1.6 PB of integrals**, against 134.8 TB of RAM in the entire
552-node icelake partition — roughly **12× the whole machine**. On-the-fly RPA
at this system size is memory-infeasible, not merely expensive, which is the
cleanest possible statement of why the MLIP route is not just cheaper but
*enabling*: it delivers RPA-quality forces at a system size and timescale
where RPA itself cannot be evaluated even once.
