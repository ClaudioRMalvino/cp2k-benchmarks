**Disclaimer (CSD3 ptrace mitigation).** The MAQAO measurements presented in
this report were collected on 14 May 2026, prior to a kernel-level ptrace
mitigation applied across CSD3 on 15 May 2026 in response to a Red Hat / Rocky
Linux zero-day security advisory (HPC support communication, S. Rankin, 15 May
2026). The mitigation prevents `ptrace(PTRACE_ATTACH, ...)` -- the syscall
MAQAO's LProf sampler relies on -- from being used by unprivileged processes,
and a permanent fix awaits a Red Hat / Rocky Linux update at the time of
writing. As a consequence, MAQAO could not be re-run against the post-port
build of `feature/nnp-native-spline`. The cache-locality and hotspot
measurements reported below for "cache-conscious" code are therefore taken
from `feature/nnp-native-spline-omp`, which already contained the cache-conscious
optimisations (precomputed `scale_factor`/`inv_sigma_factor`, sorted
`active_atoms` via in-place heap sort) at the time of the May 14 measurement,
and is algorithmically identical to the post-port pure-MPI branch on the
cache-relevant code paths (only the OpenMP threading layer differs).

The static `array_access_efficiency` metric used here is derived by MAQAO's
CQA module from inspection of load/store stride patterns at the assembly level
and is independent of the blocked hardware performance counters; the
absolute number is therefore unaffected by the mitigation and is directly
comparable across the two branches shown.
