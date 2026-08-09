# Draft PR body — cp2k/cp2k, branch `ClaudioRMalvino/cp2k:nnp-omp-atom-loop`

Open as a **DRAFT** PR (Dhruv's request: "rebase it to the current cp2k master
and open a draft pr to cp2k. I can check over it there and comment in the open
pr if there are any issues"). Base `cp2k/cp2k:master`, head
`ClaudioRMalvino/cp2k:nnp-omp-atom-loop`.

Compare URL:
`https://github.com/cp2k/cp2k/compare/master...ClaudioRMalvino:cp2k:nnp-omp-atom-loop`

PR title (GitHub prefills this from the single commit's subject, so it should
need no typing): `NNP: Move OpenMP threading to the atom loop`

Attach `plots/omp_headtohead/omp_headtohead_scaling.png` where the text points
at it — GitHub will not render a repo-relative path from a fork.

Everything below the line is the PR body.

---

Follow-up to #5295. That PR made the NNP kernel scale under MPI; this one makes
it scale under OpenMP as well, by threading over atomic centres instead of over
symmetry-function groups, and it does so without changing a single computed
number at any thread count.

### What changes

The parallel region inside `nnp_calc_acsf` is replaced by threading over the
loop of atomic centres in `nnp_force.F`. Each thread accumulates into private
arrays — including a private copy of the committee networks — and the partial
results are folded together in ascending thread order once the region ends.
That fixed fold order is the whole point: the floating-point summation sequence
becomes identical to the serial loop, so results are bit-identical to serial for
any number of threads, not merely equal to within a tolerance.

This supersedes the angular-group OpenMP path, so `nnp_acsf_angular_loop_omp`
and its call site are removed. That path parallelised over symmetry-function
groups within a single centre and cannot coexist with threading over centres.
Its inline serial equivalent is retained and is what every thread now executes,
so the arithmetic performed per centre is unchanged.

Supporting changes: the persistent workspace in `nnp_environment_types` gains a
thread dimension so threads no longer share scratch space;
`nnp_neighbor_interface` default-initialises the widened fields for
`FTN_NO_DEFAULT_INIT` builds; and the helium call site in
`motion/helium_interactions.F` passes the shared networks explicitly — that path
stays serial and is numerically untouched.

Five files, 485 insertions, 559 deletions.

### Why

Measured on CSD3 icelake, one node, 76 cores, N = 1024 H₂O (3072 atoms), 100 MD
steps, 5 repetitions plus warm-up, t/step from `qs_mol_dyn_low`. Both binaries
were built from the **same configured CMake tree** (Intel 2022.1.0, Release,
`-O3 -xCORE-AVX2 -fp-model=precise -qopenmp`) and run **in the same job on the
same node**, so the only variable is this patch. The baseline is #5295's head
`8be1dfa50f`, whose NNP code is what merged as `6908c592bf` apart from cosmetic
review changes.

| cores | baseline OMP [s/step] | this PR OMP [s/step] | ratio | baseline MPI [s/step] | this PR MPI [s/step] |
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

- **50.9× at 76 threads.** The current threading does not speed the kernel up on
  this system — at 76 threads it is 0.79× its own serial run — so the ratio is
  large. Against serial, this PR is a 39.4× speedup, 51.8 % parallel efficiency
  on the node.
- **Threading now beats pure MPI at full node.** 0.0543 s/step at 76 threads vs
  0.0835 s/step for the best full-node MPI configuration, a factor 1.54. MPI
  saturates beyond about 38 ranks on this system.
- **No serial cost.** At one core this PR is 2.3 % faster (2.137 vs 2.188).
- **No MPI regression.** Pure-MPI timings match or slightly beat the baseline at
  every rank count (76 ranks: 0.0807 vs 0.0835).

![OMP and MPI scaling](omp_headtohead_scaling.png)

### Correctness

Bit-exactness is the design constraint, not a happy accident, so it is tested
directly rather than through tolerances.

**Cross-binary, N = 1024 (3072 atoms).** Forces from this PR against the
unmodified baseline binary, same input, both MD frames:

| comparison | max abs force difference | energy difference |
|---|---:|---:|
| baseline OMP=1 vs this PR OMP=1 | 0.000e+00 | 0.0 Ha |
| baseline OMP=1 vs this PR OMP=76 | 0.000e+00 | 0.0 Ha |

**On this branch, rebased onto current master:**

- `NNP/regtest-1` via `do_regtest.py` at OMP = 1, 2, 4, 8 → **10/10 CORRECT at
  every thread count (40/40)**, including the `H2O-256_C-NNP_MD-bitexact`
  fixture added by #5295.
- H₂O-64 step-0 per-atom forces, OMP=1 vs OMP=4, 192 atoms →
  **max_abs = 0.000e+00, max_rel = 0.000e+00**. Bit-identical, not merely within
  1e-14.

### Scope and limitations

Parallelism moves from symmetry-function groups to atomic centres, so a case
with very few centres per MPI rank but several angular groups loses parallelism
this PR does not replace. For any production-sized system centres per rank
greatly exceeds group count, which is the regime the measurements above cover.
I have not tried to make the two levels coexist; a fixed reduction order over
centres is what buys bit-exactness, and nesting a second level inside it would
put that at risk.

Reproduction material — raw sweeps, provenance and the plotting script — is
available on request; I did not want to add benchmark artefacts to the tree.
