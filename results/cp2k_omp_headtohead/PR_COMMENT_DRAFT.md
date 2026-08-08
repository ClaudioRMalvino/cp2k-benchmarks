# ARCHIVE — PR #5295 outreach, resolved 2026-07-27

**Outcome:** Claudio opened the mini-PR against Dhruv's branch, then CLOSED
it — #5295 was pending 1/2 approvals and Dhruv is unlikely to take new
commits pre-merge. He left a comment saying he'll open the patch as a
direct PR to cp2k/cp2k once #5295 merges.

**2026-07-27 comment cleanup (amended, force-pushed):** commit is now
`f74941bfb7` — comment-only changes (verified: no code lines touched, so
the validated binaries still correspond). Removed change-narration and
stale inner-OMP comments, added the missing doxygen \param lines
(arc/rad_y/ang_y/tid), added `\author Claudio Malvino (claudiormal@gmail.com)`
to the nnp_force.F/nnp_acsf.F module headers and to nnp_calc_energy_force /
nnp_calc_acsf. OPEN DECISION for the direct PR: the dead
nnp_acsf_angular_loop_omp helper behind IF (.FALSE.) — recommend deleting
it (and the IF branch) at rebase time, when regtests re-run anyway;
upstream reviewers will not accept dead code.

**NEXT STEP (when #5295 merges):**
1. `git remote add upstream https://github.com/cp2k/cp2k.git` (worktree) +
   `git fetch upstream master`.
2. `git rebase --onto upstream/master 8be1dfa50f omp-atom-bitexact-v2` —
   replays just the graft commit onto the new master.
3. Re-run the validation gates on the rebased build (NNP regtests at
   OMP 1/2/4/8; Tier-2 cross-binary force diff) before going public.
4. Push to the fork; Claudio opens a standalone PR against `cp2k/cp2k`
   master. Commit message stays as-is; the benchmark table, figure and
   rationale (drafted below) go in the PR body — per the data-in-PR rule.
   This route gives the cleanest credit: Claudio-authored commit merges
   upstream in his own PR, no squash-collapse inside #5295.

The comment below was drafted for #5295 and is kept as source material for
that future PR body (table, framing, bit-exactness argument). Figure:
`plots/omp_headtohead/omp_headtohead_scaling.png`, regenerate with
`python plots/plot_omp_headtohead.py`.

---

Hi @DhruvSkyy — thanks for this PR. I've been benchmarking it on our cluster
(CSD3, Icelake, 76 cores/node) as part of my MPhil project at Cambridge, and
the cell-list neighbour search is a really nice improvement over the current
quadratic search.

While testing hybrid MPI/OpenMP setups I noticed that the threaded path
doesn't currently get anything out of OpenMP: on a 3072-atom water box
(N = 1024 H₂O), the PR head at 76 threads actually runs a little slower than
its own serial time (2.76 vs 2.19 s/step). The parallel region inside
`nnp_calc_acsf` sits too deep — the work per atomic centre is small, so
fork/join overhead and the CRITICAL sections eat the gains.

I've put together a patch on top of your branch that moves the threading up
one level, to the loop over atomic centres in `nnp_force.F`. The catch at that
level is reproducibility, so each thread accumulates into its own private
arrays (including a private copy of the committee networks) and the partial
results are summed in thread order once the region ends. That keeps the
floating-point summation order identical to the serial loop, and the forces
come out **bit-identical to serial for any number of threads** — max abs
force diff 0.0 against your unpatched head at both 1 and 76 threads, and the
NNP regtests pass at 1, 2, 4 and 8 threads.

Timings below. Both binaries were built from the same configured CMake tree
(ifort 2022.1.0, `-O3 -xCORE-AVX2 -fp-model=precise -qopenmp`) and run in the
same job on the same node, so the only variable is the patch:

| cores | head, OMP [s/step] | patched, OMP [s/step] | ratio | head, pure MPI | patched, pure MPI |
|---:|---:|---:|---:|---:|---:|
| 1  | 2.188 ± 0.011 | 2.137 ± 0.013 | 1.02× | 2.190 | 2.141 |
| 2  | 2.485 ± 0.029 | 1.171 ± 0.008 | 2.12× | 1.222 | 1.189 |
| 4  | 2.410 ± 0.033 | 0.644 ± 0.003 | 3.74× | 0.670 | 0.644 |
| 8  | 2.413 ± 0.020 | 0.328 ± 0.001 | 7.37× | 0.340 | 0.328 |
| 16 | 2.450 ± 0.026 | 0.170 ± 0.001 | 14.4× | 0.176 | 0.171 |
| 32 | 2.600 ± 0.027 | 0.096 ± 0.000 | 27.0× | 0.096 | 0.093 |
| 76 | 2.761 ± 0.015 | **0.054 ± 0.002** | **50.9×** | 0.084 | 0.081 |

A few things worth pulling out of that:

* serial and pure-MPI timings are unchanged (if anything the patched build is
  marginally faster at OMP=1) — the patch only touches how the threaded path
  is organised;
* at full node, 76 threads now beats 76 ranks (0.054 vs 0.084 s/step), where
  pure MPI saturates beyond ~38 ranks — useful for memory-heavy committee
  models, since threads share one copy of the network;
* the helium path in `motion/helium_interactions.F` just gets the networks
  passed explicitly and stays serial — numerically untouched.

The patch is a single commit (5 files) on top of your current head, on my
fork here:
https://github.com/ClaudioRMalvino/cp2k/commit/f74941bfb73a6ef7b7f77ca8d033a654e246a7f3
To be clear, this isn't meant to hold anything up — the PR is a clear
improvement as it stands and this works perfectly well as a follow-up. I'm
happy to open it as a small separate PR once this merges, or as a PR against
your branch now if you'd rather fold it in — whatever suits how you want to
take this forward. And if you'd like the raw benchmark data or the run
inputs, just say.

---

# Posting checklist (for us, not part of the comment)

1. DONE — branch pushed to github.com/ClaudioRMalvino/cp2k (2026-07-27).
2. DONE — handle (@DhruvSkyy) filled in; his branch confirmed as
   `pr/nnp-cpu-cell-list`, still at 8be1dfa50f (our commit applies clean).
3. Claudio opens the mini-PR (step 1 above), verifying the base repo header.
4. Claudio posts the comment (step 2), filling in `<mini-PR-link>` and
   attaching the figure.
Direct commit link, if ever needed:
https://github.com/ClaudioRMalvino/cp2k/commit/f74941bfb73a6ef7b7f77ca8d033a654e246a7f3
