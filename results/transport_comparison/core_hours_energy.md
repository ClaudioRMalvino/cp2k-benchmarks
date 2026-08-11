# Core-hours for the RPA transport series — energy/CO2 calculator inputs

Written 2026-08-11. All rates are the measured values of the report
(Table 9.4): `master` 4.05 s/step, optimised CPU 0.26 s/step, optimised
GPU 0.17 s/step, all on the 5064-atom, 37.3 Å series cell. Only the
optimised CPU ran the full series; the `master` and GPU totals project
those measured rates over the same steps. Node-hour canon matches the
report after the 2026-08-11 factor-2 correction (340 / 5,396 / 227).

## Hardware (calculator inputs)

| | CPU track | GPU track |
|---|---|---|
| Node | CSD3 Ice Lake, 2× Intel Xeon Platinum 8368Q, 76 cores | CSD3 Wilkes3, NVIDIA A100-80GB |
| Unit billed | core-hour (76 per node-hour) | GPU-hour (+3 CPU cores for SPME) |
| TDP basis used in our estimates | 700 W per node (2× 270 W CPUs + board/memory) | 550 W per A100 incl. host share |
| Location | Cambridge, UK | Cambridge, UK |
| PUE | unknown (CSD3 publishes none) | unknown |

## 1. The series as run: 3 concentrations × 5 × 160 ps = 2.4 ns (4.8 M MD steps)

| Code | s/step | Wall-time on one unit | Node-hours | Core-hours | GPU-hours | kWh (TDP est.) |
|---|---|---|---|---|---|---|
| `master` (non-optimised) | 4.05 | 225 days (1 node) | 5,396 | 410,000 | — | ~3,780 |
| Optimised CPU | 0.26 | 14 days (1 node) | 340 | 25,840 | — | ~240 |
| Optimised GPU | 0.17 | 9.5 days (1 A100) | — | 680 (3 SPME cores) | 227 | ~125 |

Equilibration (30 ps × 3 ≈ 0.18 M steps) adds ~4% to each row if you
want the all-in number; the report's money table uses that basis.

## 2. Extended to 200 ns of trajectory (400 M MD steps, the O'Neill-scale sampling)

| Code | Serial time on one unit (derived) | Node-hours | Core-hours | GPU-hours | kWh (TDP est.) |
|---|---|---|---|---|---|
| `master` (non-optimised) | 51 years (1 node) | 450,000 | 34,200,000 | — | ~315,000 |
| Optimised CPU | 3.2 years (1 node) | 28,000 | 2,128,000 | — | ~19,600 |
| Optimised GPU | 2.2 years (1 A100) | — | 57,000 (SPME) | 19,000 | ~10,450 |

The node-hours, GPU-hours and £ equivalents (£342,000 / £21,600 /
£10,400 at CSD3 internal rates) are the report's §9.5 numbers. The
"serial time" column is NOT in the report; it is node-hours divided by
8,760 (hours per year), the time one unit would need running
non-stop. Real campaigns run independent segments in parallel, so no
one would ever wait this long; the column only illustrates scale.

## 3. If the whole series ran to 200 ns per concentration (600 ns, 1.2 G steps)

| Code | Node-hours | Core-hours | GPU-hours | kWh (TDP est.) |
|---|---|---|---|---|
| `master` (non-optimised) | 1,350,000 | 102,600,000 | — | ~945,000 |
| Optimised CPU | 84,000 | 6,384,000 | — | ~58,800 |
| Optimised GPU | — | 171,000 (SPME) | 57,000 | ~31,350 |

## Caveats

- The kWh column is a TDP estimate, not a measurement: CSD3's SLURM
  energy-accounting plugin is disabled site-wide, and the report
  deliberately makes no energy claim for this reason. Use it only to
  cross-check the external tool.
- `master` and GPU rows are projections of measured per-step rates;
  only the optimised CPU actually ran the 4.8 M production steps
  (actual all-in consumption 582 node-hours including equilibration,
  support runs and re-runs).
- If the tool asks for a carbon intensity, the UK grid is the right
  region; CSD3 publishes no PUE or emissions figure (ARCHER2, for
  contrast, publishes 0.014 kgCO2e per node-hour, embodied only, on a
  renewable supply).
- If a number from the tool goes into the report, cite the tool (e.g.
  Green Algorithms, Lannelongue et al., Adv. Sci. 2021,
  doi:10.1002/advs.202100707) and label the result an estimate.
