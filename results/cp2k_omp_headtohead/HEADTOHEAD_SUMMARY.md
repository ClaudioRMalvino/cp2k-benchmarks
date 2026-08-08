# OMP graft head-to-head: PR #5295 head vs head + bit-exact atom-level OpenMP

**Job 31837144** (CSD3 icelake, 1 node / 76 cores, 2026-07-23, 6 h 24 m).
Both binaries built from the **same configured CMake tree** (Intel 2022.1.0,
`-O3 -xCORE-AVX2 -fp-model=precise -qopenmp`, Release, same toolchain deps),
run in the **same job on the same node** — the only variable is the graft.

| binary | source | md5 |
|---|---|---|
| `dhruv-head-8be1dfa` | pristine PR head `8be1dfa50f` | `3a4307dbebdc47a9af6ce6ff318c29ff` |
| `dhruv-omp-bitexact-v2` | same head + uncommitted graft (5 files, 251+/120−) | `98a473a2f3f0aa9ffea416063d21a68c` |

## Leg 1 — accuracy (N = 1024 H₂O, 3072 atoms, cross-binary)

Forces of the graft binary vs the pristine binary, same input, both MD frames:

| comparison | max abs force diff | energy diff |
|---|---:|---:|
| pristine (OMP=1) vs graft (OMP=1) | **0.000e+00 — bit-identical** | 0.0 Ha |
| pristine (OMP=1) vs graft (OMP=76) | **0.000e+00 — bit-identical** | 0.0 Ha |

Backs the earlier regtest evidence (40/40 CORRECT at OMP = 1/2/4/8; Tier-2
force diff 0.0 between thread counts). The graft changes *no* physics at any
thread count — the deterministic ascending-thread-order fold makes the
reduction order identical to serial.

## Legs 2+3 — performance (N = 1024, 100 MD steps, 5 reps + warm-up, t/step from `qs_mol_dyn_low`)

| cores | head OMP t/step [s] | graft OMP t/step [s] | OMP ratio | head MPI t/step [s] | graft MPI t/step [s] |
|---:|---:|---:|---:|---:|---:|
| 1 | 2.1879 ± 0.0111 | 2.1367 ± 0.0126 | 1.02× | 2.1897 | 2.1410 |
| 2 | 2.4849 ± 0.0291 | 1.1713 ± 0.0080 | 2.12× | 1.2216 | 1.1892 |
| 4 | 2.4102 ± 0.0332 | 0.6441 ± 0.0026 | 3.74× | 0.6698 | 0.6443 |
| 8 | 2.4131 ± 0.0195 | 0.3275 ± 0.0007 | 7.37× | 0.3403 | 0.3280 |
| 16 | 2.4499 ± 0.0263 | 0.1700 ± 0.0005 | 14.41× | 0.1755 | 0.1710 |
| 19 | 2.4703 ± 0.0192 | 0.1464 ± 0.0004 | 16.88× | 0.1489 | 0.1445 |
| 32 | 2.5995 ± 0.0273 | 0.0964 ± 0.0001 | 26.96× | 0.0959 | 0.0926 |
| 38 | 2.7203 ± 0.0304 | 0.0890 ± 0.0014 | 30.56× | 0.0831 | 0.0804 |
| 76 | 2.7614 ± 0.0154 | **0.0543 ± 0.0021** | **50.88×** | 0.0835 | 0.0807 |

### Headline numbers

* **50.9× faster at 76 OMP threads** — the PR head's threading is a net
  *slowdown* (0.79× at 76 threads vs its own serial run); the graft reaches
  39.4× speedup (51.8 % efficiency) on the same node.
* **Graft OMP-76 beats the baseline's best possible full-node configuration**
  (76 MPI ranks, 0.0835 s/step) by **1.54×** — threading now outperforms pure
  MPI at full node, where MPI saturates beyond 38 ranks.
* **No serial cost**: graft OMP=1 is 2.3 % *faster* than the head (2.137 vs
  2.188 s/step).
* **No MPI regression**: graft pure-MPI matches or slightly beats the head at
  every rank count (e.g. 76 ranks: 0.0807 vs 0.0835 s/step).

Figure: `plots/omp_headtohead/omp_headtohead_scaling.[png|pdf]`.
Raw data: `results/cp2k_omp_headtohead/NNP/` (4 sweeps × 9 counts × 5 reps).
Job log (leg-1 verdicts): `logs/NNP_omp_h2h_31837144.out`.
