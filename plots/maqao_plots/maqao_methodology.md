# MAQAO CQA Methodology
This document records how the per-branch MAQAO CQA measurements that feed `maqao_data_analysis.py` were produced, and which Hoefler & Belli SC'15 rule each downstream table / plot satisfies.

## Environment
- **Cluster**: CSD3 Peta4-IceLake (Intel Xeon Platinum 8368Q, 2×38 = 76 cores/node, AVX-512)
- **Execution host**: CSD3 login node (same IceLake silicon as compute nodes; CQA is static disassembly so no compute-node ptrace is required).
- **Kernel mitigation**: `kernel.yama.ptrace_scope = 2` (LProf disabled cluster-wide; CQA unaffected).
- **MAQAO version**: `maqao 2026.0.0-b - d53714498d38428ad6a75949c25b07c813f07f11::20260206-105209`
- **CQA invocation**: `maqao cqa <binary> --fct-loops="nnp_" --output-format=csv --output-path=<csv>`

## Binaries analyzed (md5)
| Branch | Binary path | md5 |
|---|---|---|
| master | `/rds/user/crm98/hpc-work/cp2k_binaries/csd3/master/cp2k.psmp` | `17e0035070eadd4f91319c2d518158a1` |
| native-spline (post-port) | `/home/crm98/cp2k_optimized/install/bin/cp2k.psmp` | `57b36a22a08d154a6e2c4f67508680db` |
| native-spline-omp | `/rds/user/crm98/hpc-work/cp2k_binaries/csd3/feature-nnp-native-spline-omp/cp2k.psmp` | `54a578c1469174a4a761c208841d4140` |

## Fresh CSV inputs
| Branch | CSV | Last written |
|---|---|---|
| master | `/home/crm98/cp2k-benchmarks/maqao_login_cqa/master_nnp.csv` | 2026-05-26 19:49:44 |
| native-spline (post-port) | `/home/crm98/cp2k-benchmarks/maqao_login_cqa/native_spline_postport_nnp.csv` | 2026-05-26 19:34:34 |
| native-spline-omp | `/home/crm98/cp2k-benchmarks/maqao_login_cqa/omp_nnp.csv` | 2026-05-26 19:51:17 |

## Hoefler & Belli rule compliance

| Rule | Where it applies | How this script satisfies it |
|---|---|---|
| **R1** Speedup needs base-case absolute | `fresh_per_loop_table.csv`, `fresh_speedup_with_baseline.pdf` | `headroom (base/fullvec)` column always paired with the underlying `Nb cycles: min` it is computed from. The headroom plot shows absolute cycles, not the ratio. |
| **R2** Reason for subsets | `fresh_per_loop_table.csv` (top N = 20) | Sorted by `cycles/iter` and truncated to top 20 per branch. The unfiltered full per-loop CSV is also written so the truncation is auditable. |
| **R3** Arithmetic for costs, harmonic for rates | `fresh_branch_aggregate.csv` | Cycles/iter, FLOP, bytes use arithmetic mean. IPC and FLOP/cycle use harmonic mean. |
| **R4** Avoid summarizing ratios; geomean if must | `fresh_branch_aggregate.csv` | Vec ratio / vec efficiency report geometric mean across the positive tail, with the zero-tail count broken out separately so the distribution is not collapsed into one misleading number. |
| **R5** Deterministic vs nondeterministic | All fresh-CSV outputs | CQA is pure static disassembly of the compiled binary; same binary → identical numbers. No CIs are computed because the measurement has zero variance. |
| **R8** Investigate central-tendency choice | `fresh_branch_aggregate.csv`, `fresh_vec_ratio_distribution.pdf` | Cycles/iter is reported as mean + median + p25 + p75. Vec ratio is reported with full distribution (boxplot) plus zero-tail count. |
| **R9** Document factors and environment | this file | Cluster, kernel mitigation, MAQAO version, exact invocation, binary md5, CSV mtimes. |
| **R11** Upper performance bounds | `fresh_branch_aggregate.csv`, `fresh_speedup_with_baseline.pdf` | CQA's `Nb cycles if fully vectorized / FP arith vectorized / clean` columns are propagated as explicit upper-bound models alongside the achieved values. |
| **R12** Plot only valid trends | `fresh_*` plots | Bars and box plots throughout; no spurious line interpolation between non-comparable points. |

## What about LProf (coverage %) ?
LProf is the dynamic sampling component of MAQAO ONE View and is currently blocked on CSD3 by the `ptrace_scope = 2` kernel mitigation. The legacy `HOTSPOTS` / `CQA_LOOPS` dicts in this script (from the May-14 ONE View run, before the mitigation) retain coverage% values; the fresh CSV-driven outputs do **not** include coverage% and therefore cannot weight per-loop aggregates by execution frequency. Aggregates across loops are equal-weighted, with the loop-by-loop distribution always available in the per-loop CSV so the reader can re-weight by hand if desired.
