# NNP Implementation Changes — Cross-Branch Reference

This document records the implementation differences between the NNP-related
branches of CP2K under active benchmarking in this project (upstream master,
feature/nnp-native-spline, feature/nnp-native-spline-omp, and
feature/nnp-chebyshev), why each branch exists, and what the benchmark suite
is actually measuring.

> **Revision note (June 2026)** — Added Branch 4, `feature/nnp-chebyshev`,
> built on top of `nnp-native-spline-omp`.  It replaces the Hermite spline
> lookup tables with table-free Chebyshev/minimax polynomials and adds a
> post-thesis optimisation pass (MPI-collective shrink, one-time element
> sort, reverse-mode backpropagation, walk-bound neighbour-buffer sizing,
> halo-filtered image table).  This is now the active development branch.
>
> **Revision note (May 2026)** — This document was rewritten after commit
> `7cbc8b3008` ("nnp: port cache-conscious optimisations from
> feature/nnp-native-spline-omp") merged the precomputed-scaling-factor,
> sorted-active-atoms, and Horner-form Hermite-spline changes from the OMP
> fork back into the pure-MPI native-spline branch.  Earlier revisions of this
> document described those three changes as OMP-exclusive.  An earlier
> `feature/nnp-verlet-cells` branch is out of scope for the current
> benchmarks and is not described here; see git history if you need to
> reconstruct its state.

```
Branches under benchmark (all rooted at upstream master):

   upstream/master                       ← baseline (no NNP-specific optimisations)
       │
       ├── feature/nnp-native-spline     ← cell-list / Verlet-skin neighbour walk +
       │                                    persistent module-level buffers +
       │                                    native CP2K splines (splines_methods.F)
       │                                    in Horner form +
       │                                    cache-miss reductions in nnp_scale_acsf
       │                                    (precomputed scaling factors, sorted
       │                                     active_atoms via in-place heap sort) +
       │                                    inner-loop micro-optimisations (fc/dfc/dr
       │                                    hoisting, sqrt skip via cutoff_sq,
       │                                    ASSOCIATE aliases, local cut_type)
       │
       └── feature/nnp-native-spline-omp ← native-spline base + OpenMP threading
           │                                (THREADPRIVATE scratch, OMP PARALLEL DO
           │                                 with local-reduction + CRITICAL pattern,
           │                                 COLLAPSE(2) on rectangular angular j-k,
           │                                 SCHEDULE(dynamic, chunk) tuning)
           │
           └── feature/nnp-chebyshev     ← table-free minimax polynomials replacing
                                            the Hermite spline tables + post-thesis
                                            optimisation pass: shrunk MPI collectives,
                                            one-time element sort, reverse-mode
                                            backprop in nnp_gradients, FORCEINLINE on
                                            the minimax kernels, walk-bound neighbour
                                            buffer sizing, halo-filtered image table
```

---

## Why two feature branches, both benchmarked against master

The benchmark exposes two distinct optimisation strategies layered cleanly on
top of one another, with the comparison structure designed so each layer's
contribution is independently measurable:

| Concern | Affected branches | What the comparison answers |
|---|---|---|
| Combined effect of all serial-friendly optimisations | `master` vs `nnp-native-spline` | What does the optimisation pass as a whole buy?  Includes neighbour walk, spline evaluation, inner-loop hoisting, cache-miss reductions.  Largest expected jump. |
| Pure threading overhead vs threading benefit | `nnp-native-spline` vs `nnp-native-spline-omp` | Layering OMP on the already-optimised serial baseline keeps every non-OMP variable fixed across the comparison.  The measured difference is therefore *only* the OMP scaffolding (THREADPRIVATE, PARALLEL regions, CRITICAL reductions, schedule choices) net of the work-sharing speedup.  This is the right pairing for "does hybrid MPI+OMP recover scaling that pure MPI loses past ~16 ranks?" |
| Maintenance burden of new public modules | `nnp-native-spline` / `nnp-native-spline-omp` add public symbols to `splines_methods.F` | This change touches CP2K core and benefits the wider codebase; we want to see the real-world cost in the NNP context |
| Hyperthread contention | `nnp-native-spline-omp` | Icelake (76 cores per node) is sensitive to physical-vs-logical core binding; pure-MPI baselines on `nnp-native-spline` give us the same-implementation reference for "what should be possible" |

The benchmark harness builds all three binaries (master + two feature
branches) into separate per-branch directories so `LD_LIBRARY_PATH` cannot
accidentally pick up the wrong `libcp2k.so`.  Each branch is then run through
the same size-scaling and core-scaling sweeps.

---

## Branch 1 — `upstream/master` (baseline)

`src/nnp_acsf.F` = 1 414 lines.  Reference implementation as merged into CP2K
upstream.

Performance characteristics relevant to our work:

- **Neighbour search**: per centre atom calls `nnp_compute_pbc_copies` then
  loops over all `(j, c1, c2, c3)` periodic images, calling `pbc()` on every
  displacement.  Cost: O(N · (2p+1)³) per centre, plus ~9 ALLOCATE/DEALLOCATE
  calls per centre for buffer setup.
- **Spline evaluation**: the cutoff function uses direct `COS(π·r/rc)` and
  `TANH(1−r/rc)` libm calls.  The radial Gaussian uses `EXP(−η(r−rs)²)`.
  Each is ~30–80 cycles per evaluation and blocks vectorisation of the
  enclosing loop.
- **Force accumulation**: outer atom loop in `nnp_force.F` iterates over all
  `nnp%num_atoms` per centre even though only ~k actually contribute.
- **MPI atom assignment**: `allgather(mecalc, allcalc)` collective + prefix
  sum to compute each rank's `istart`.
- **Repeated derived-type field walks**: every read of
  `nnp%ang(ind)%symfgrp(s)%cutoff_sq` is a multi-level pointer chase pulling
  multiple cache lines per access.
- **Repeated divisions**: `drdx = rvect/r` is three separate divides in the
  inner loop; cutoff-argument divisions (`r/rc`, `π·r/rc`) are recomputed on
  every iteration.

The feature branches treat master as the truth source for correctness and
the "slow but always-right" reference for benchmarks.

---

## Branch 2 — `feature/nnp-native-spline`

`src/nnp_acsf.F` = 2 309 lines.  Combines the foundational neighbour-walk and
inner-loop optimisations with the cache-miss reductions and Horner-form
spline backported from the OMP fork.  This is the pure-MPI optimised baseline
against which OMP threading is measured.

### 2.1 New module-level types

#### `nnp_cell_list_cache_type`

A Verlet-skin / cell-list hybrid cache.  Each `(atom j, integer PBC shift)`
pair becomes one *image entry* stored explicitly, so the inner neighbour-query
loop contains no `pbc()` calls.  The coarse grid has `bin_width ≈ list_cutoff`
(one bin-span walk covers 27 bins).

| Field | Purpose |
|---|---|
| `coord_primary(3, N)` | PBC-wrapped primary atom positions, refreshed every step |
| `ref_coord_primary(3, N)` | Snapshot at last rebuild (Verlet reference) |
| `image_atom(n_img)` | Source atom index for each image entry |
| `image_translation(3, n_img)` | Cartesian translation vector for each image |
| `head(n_cells)` / `next(n_img)` | Linked-list cell grid (LIFO insertion) |
| `verlet_skin` | Skin in Å; `list_cutoff = max_cut + skin` |
| `exact_cutoff`, `hmat`, `perd` | Change-detection fields for rebuild triggers |

A `SAVE` instance `cell_list_cache` lives at module level and persists across
force evaluations.

### 2.2 New module-level SAVE variables

```fortran
TYPE(nnp_cell_list_cache_type), SAVE, PRIVATE :: cell_list_cache

TYPE(nnp_neighbor_type), SAVE, PRIVATE, TARGET :: persistent_neighbor
INTEGER, SAVE, PRIVATE :: persistent_n_capacity = -1
INTEGER, SAVE, PRIVATE :: persistent_n_rad_grp  = -1
INTEGER, SAVE, PRIVATE :: persistent_n_ang_grp  = -1

REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_sym(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_forcetmp(:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_force3tmp(:,:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE         :: scratch_fc1(:), scratch_dfc1(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE         :: scratch_fc2(:), scratch_dfc2(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE         :: scratch_dr1dx(:, :)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE         :: scratch_dr2dx(:, :)
```

`persistent_neighbor` replaces the per-centre-atom local `nnp_neighbor_type`.
In upstream, each call to `nnp_calc_acsf` allocated and deallocated 9 arrays
sized `~num_atoms × (2p+1)³`; for a 6 000-atom system those allocations cross
glibc's `mmap` threshold and become syscalls, giving O(N²) `malloc` traffic
per MD step.

### 2.3 `nnp_calc_acsf` — signature extended

```fortran
! upstream
SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz, stress)

! native-spline
SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz, stress, active_atoms, active_flag, n_active)
```

The three new optionals build a deduplicated active-atoms set (centre `i` plus
every distinct neighbour touched by the rad/ang lists), letting the caller
walk O(n_active) atoms instead of O(N) when accumulating forces.

### 2.4 `POINTER, CONTIGUOUS` scratch aliases

Module-level SAVE arrays accessed via per-call pointer aliases:

```fortran
REAL(KIND=dp), DIMENSION(:),     POINTER, CONTIGUOUS :: symtmp
REAL(KIND=dp), DIMENSION(:,:),   POINTER, CONTIGUOUS :: forcetmp
REAL(KIND=dp), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: force3tmp
! per symfgrp:
symtmp    => scratch_sym(1:n_symf)
force3tmp => scratch_force3tmp(:,:,1:n_symf)
```

`CONTIGUOUS` matters: without it the compiler treats the pointer as
potentially non-contiguous and potentially aliasing other `TARGET`/`POINTER`
variables in scope, blocking vectorisation of the inner `sf`/`l` loops.

### 2.5 Cell-list neighbour walk

Two-phase architecture:

- **`nnp_prepare_cell_list_cache`** — called once per force evaluation from
  `nnp_calc_energy_force`.  Refreshes `coord_primary`, checks six rebuild
  triggers (`!initialized` / `num_atoms` change / cutoff change / skin change /
  cell change / any-atom-moved-more-than-half-skin).  Skin defaults to
  `MIN(1.0 Å, 0.1·max_cut)`.
- **`nnp_compute_neighbors`** — called once per centre atom; walks only the
  27 bins around the centre's bin in the pre-built image table.  No
  allocations, no `pbc()` calls, lazy `SQRT` (sentinel `norm = -1.0`).

Cheap inner-loop filter chain:
1. `r² ≥ max_cut²` — scalar reject, no SQRT
2. Per-group `r² < cutoff_sq` — radial groups
3. Per-group `r² < cutoff_sq` — angular ele1 / ele2 separately

### 2.6 Bug fix — wrong `kk` index in hetero-angular path

In upstream, the hetero-angular loop (element 1 ≠ element 2) had a dead write
followed by a redundant inner-loop assignment:

```fortran
! upstream (buggy but numerically correct)
jj = neighbor%ind_ang1(j, s)
kk = neighbor%ind_ang1(k, s)   ! DEAD — immediately overwritten inside sf loop
CALL nnp_calc_ang(...)
DO sf = 1, n_symf
   m = off + ...
   jj = neighbor%ind_ang1(j, s)   ! redundant
   kk = neighbor%ind_ang2(k, s)   ! correct, but re-read every sf iteration
END DO
```

The outer `kk` was the wrong index (should be `ind_ang2`) but immediately
shadowed by the correct inner assignment, so the result was numerically
correct.  native-spline fixes the index and removes the redundant inner reads.

### 2.7 Native CP2K spline replacement of `EXP`, `COS`, `TANH`

The expensive intrinsic calls in `nnp_calc_sym_rad` and `nnp_calc_ang` are
replaced by Hermite-cubic splines living in `splines_methods.F` (the CP2K
core, not a module-local cache):

| Intrinsic | Spline range | Rationale |
|---|---|---|
| `EXP(x)` | `[−50, 0]` | Argument is always ≤ 0 (Gaussian decay) |
| `COS(π·r/rc)` | `[0, π]` | Cutoff arg ∈ `[0, 1)` so `tmp ∈ [0, π]` only — doubles grid resolution |
| `TANH(1−r/rc)` | `[0, 1]` | Argument lives in `[0, 1)` — 20× finer resolution |

5 000 grid points each.  For `EXP` specifically, the spline value (~`EXP`) is
used directly in the chain rule for the force rather than evaluating the
formal spline derivative — same answer, slightly more accurate (O(h⁴) value
vs O(h³) derivative).

The two new public symbols in `splines_methods.F`:

```fortran
PUBLIC :: init_hermite_spline   ! init Hermite cubic spline (Horner-form coefficients)
PUBLIC :: hermite_spline_value  ! evaluate Hermite cubic spline (+ optional derivative)
```

29 call sites in NNP code.  The initialiser populates
`spline_data%coef(4, n-1)` with per-interval Horner coefficients computed
once at init time:

```fortran
! For each interval [x_i, x_{i+1}] with t = (x - x_i)*invh ∈ [0, 1]:
!     p(t) = a_0 + a_1·t + a_2·t² + a_3·t³
!   where  a_0 = y_i,  a_1 = h·y'_i,
!          a_2 = -3·y_i + 3·y_{i+1} - 2·h·y'_i -   h·y'_{i+1},
!          a_3 =  2·y_i - 2·y_{i+1} +   h·y'_i +   h·y'_{i+1}.
```

The evaluator then uses pure Horner form:

```fortran
val = c1 + t*(c2 + t*(c3 + t*c4))               ! ~3 MUL + 3 ADD
y1  = invh*(c2 + t*(2.0_dp*c3 + t*3.0_dp*c4))   ! ~3 MUL + 3 ADD (if derivative needed)
```

vs ~22 MUL + 12 ADD per call for the basis-function form that was previously
used.  Mathematically identical (verified by symbolic expansion); the
per-evaluation cost is roughly 4× lower because the polynomial coefficients
are pre-collected once at init.

### 2.8 Division elimination

```fortran
! upstream
drdx(:) = rvect(:)/r

! native-spline
rinv = 1.0_dp/r
drdx(1) = rvect(1)*rinv; drdx(2) = rvect(2)*rinv; drdx(3) = rvect(3)*rinv
```

Cutoff-argument divisions move to `nnp_init_acsf_groups` as precomputed
reciprocals stored in the `symfgrp` struct: `pi_over_cut = π/funccut`,
`inv_cut = 1/funccut`, `cutoff_sq = funccut²`.

### 2.9 `nnp_calc_ang` — hoisted loop invariants

```fortran
r_sum_sq = rsqr1 + rsqr2 + rsqr3   ! EXP argument; only η varies per sf
g_sq_inv = 1.0_dp/(g*g)             ! denominator in dangular/dr
```

In upstream these were recomputed on every `sf` iteration.

### 2.10 Pre-computed `zeta_int` / `zeta_is_int`

Upstream did `i = NINT(zeta); IF (1.0_dp*i == zeta) ...` inside the angular
SF loop on every call.  native-spline precomputes `zeta_int(:)` and
`zeta_is_int(:)` once at init, eliminating the runtime NINT and FP-equality
test from every SF evaluation.

### 2.11 Hoist `f_c(r1)` / `f_c(r2)` out of the `k`-loop

The angular routine reads the cutoff function for `r1` and `r2` independently
of `k`, but every iteration of the inner `(j, k)` loop in upstream recomputed
them.  native-spline precomputes them per neighbour:

```fortran
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: scratch_fc1(:), scratch_dfc1(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: scratch_fc2(:), scratch_dfc2(:)

! Per s-iteration, before entering the (j,k) loop:
DO j = 1, neighbor%n_ang1(s)
   fc_tmp = neighbor%dist_ang1(4, j, s) * sg%pi_over_cut
   scratch_fc1(j)  = 0.5_dp*(hermite_spline_value(spl_cos, fc_tmp, y1=dfc_tmp) + 1.0_dp)
   scratch_dfc1(j) = 0.5_dp*dfc_tmp*sg%pi_over_cut
END DO
```

`nnp_calc_ang` accepts these via new optional args `fc1_in`, `fc2_in`,
`dfc1dr_in`, `dfc2dr_in`.  For typical N_k ~ 30 inner iterations this cuts
the cutoff-function spline-evaluation count by ~30×.

### 2.12 Hoist `dr1dx` / `dr2dx` out of the `k`-loop

Same idea, applied to the radial unit vectors.  Two more SAVE scratch arrays:

```fortran
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: scratch_dr1dx(:, :)  ! (3, n)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: scratch_dr2dx(:, :)
```

`nnp_calc_ang` accepts them via new `dr1dx_in` / `dr2dx_in` optionals.  Per
`(j,k)` pair this saves two divisions and six multiplies (the second `1/r2`
unit-vector computation).

### 2.13 Skip the SQRT in the `r3` cutoff check

```fortran
rvect3 = rvect2 - rvect1
r3sq = rvect3(1)**2 + rvect3(2)**2 + rvect3(3)**2
IF (r3sq < sg%cutoff_sq) THEN
   r3 = SQRT(r3sq)                             ! only on accepted pairs
   ...
```

For typical NNP cutoffs ~50–75 % of `(j,k)` pairs miss the cutoff, so the
SQRT is skipped for most of them.

### 2.14 `ASSOCIATE` aliases for hot derived-type fields

Both `nnp_calc_rad` and `nnp_calc_ang`, plus the `(j,k)` loops in
`nnp_calc_acsf`, wrap their bodies in:

```fortran
ASSOCIATE (sg => nnp%ang(ind)%symfgrp(s))
   ! ... sg%cutoff, sg%cutoff_sq, sg%pi_over_cut, sg%inv_cut, sg%n_symf, sg%symf(...)
END ASSOCIATE
```

Inside `nnp_calc_ang`'s `sf` loop a second alias `a => nnp%ang(ind)` collapses
the per-symmetry-function `lam(m)`, `zeta(m)`, `prefzeta(m)`, `eta(m)`,
`zeta_int(m)`, `zeta_is_int(m)` reads from a multi-step pointer chase down to
single-array indexing.

### 2.15 Local `cut_type`

```fortran
cut_type = nnp%cut_type   ! read once at top of nnp_calc_acsf (and rad/ang)
SELECT CASE (cut_type)    ! register-local instead of derived-type field walk
```

Trivial change individually; meaningful in aggregate because the SELECT CASE
appears in the fc precomputation, in `nnp_calc_rad`, and in `nnp_calc_ang` —
all hot.

### 2.16 `nnp_force.F` — accompanying changes

- **Allgather elimination**: `mecalc` and `istart` derived analytically from
  `(num_atoms, num_pe, mepos)`; no MPI collective needed.
- **SAVE scratch arrays**: `dsymdxyz`, `stress`, `denergydsym`, `active_atoms`,
  `active_flag` allocated once at first force eval, reused across steps.
- **Hoisted allocation**: `dsymdxyz` sized to `max_nsym` (max over all
  elements) at force-eval entry rather than re-allocating per centre atom.
- **Selective zeroing**: `dsymdxyz(1:n_nodes(ind), :, k) = 0` only for atoms
  in the active list, not the whole `(:,:,:)`.
- **O(n_active) force accumulation loop**: only walks the touched atoms.

### 2.17 `nnp_environment_types.F` — new fields

Per `symfgrp` (radial and angular, identical):
```fortran
REAL(KIND=dp) :: cutoff_sq   = -1.0_dp   ! eliminates SQRT in neighbour filter
REAL(KIND=dp) :: inv_cut     = -1.0_dp   ! eliminates / in cutoff arg
REAL(KIND=dp) :: pi_over_cut = -1.0_dp   ! eliminates / in cutoff arg
```

Per `nnp_rad_type` and `nnp_ang_type` — added in the May 2026 backport
(commit `7cbc8b3008`):
```fortran
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: scale_factor      ! (scmax-scmin) / (loc_max-loc_min)
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: inv_sigma_factor  ! (scmax-scmin) / sigma
```

Per `nnp_ang_type`:
```fortran
INTEGER, ALLOCATABLE :: zeta_int(:)        ! NINT(zeta), precomputed
LOGICAL, ALLOCATABLE :: zeta_is_int(:)
```

Per `nnp_rad_type` and `nnp_ang_type` (infrastructure for O(1) element→group
dispatch):
```fortran
INTEGER, ALLOCATABLE :: ele_grp_count(:)
INTEGER, ALLOCATABLE :: ele_grp_list(:,:)
```

### 2.18 Cache-miss reductions in `nnp_scale_acsf` (May 2026 backport)

Cachegrind on an earlier revision flagged `nnp_scale_acsf` as the largest
contributor to L1 / L2 misses.  Two fixes ship on this branch as of commit
`7cbc8b3008`:

#### Pre-computed scaling factors

The original inner body recomputed `(scmax − scmin)/(loc_max(j) − loc_min(j))`
for every `(centre, sym, neighbour)` triple — pulling `loc_max(j)`,
`loc_min(j)`, and either `sigma(j)` or `loc_av(j)` through the cache
alongside the hot `dsymdxyz` traffic.  All four values are static for the
whole run.

The two new per-element arrays (§2.17) are filled once in `nnp_init_model`
immediately after the scaling-data file is parsed.  The hot path becomes:

```fortran
dsymdxyz(j, :, k) = dsymdxyz(j, :, k) * nnp%rad(ind)%scale_factor(j)
```

— one multiply, one load.  No divide, no extra `loc_max`/`loc_min` reads.

#### Sorted `active_atoms`

The active-atoms list was previously inserted in BFS order (the order
neighbours appear in the rad/ang lists), so the `k` walk in `nnp_scale_acsf`
jumped randomly through `dsymdxyz(:, :, k)` — defeating prefetch and
thrashing TLB.  After the BFS construction, `active_atoms(1:n_active)` is
now sorted ascending via in-place heap sort:

```fortran
IF (n_active > 1) CALL nnp_sort_active_atoms(active_atoms, n_active)
```

`nnp_sort_active_atoms` is a `PURE` private helper — O(n log n) worst case,
zero `ALLOCATE`/`DEALLOCATE`, fully in place.  The heap-sort choice over a
generic merge sort (e.g. `cp_1d_i4_sort`) is deliberate: it avoids per-call
malloc traffic that would compound at large N, while keeping the optimal
asymptotic.  `n_active` itself stays geometrically bounded (cutoff sphere ×
density) regardless of system size, so total scaling remains O(N).

### Performance — what it measures vs master

This branch is the largest-jump comparison.  Eliminates the per-pair `pbc()`
call (master's dominant cost), the malloc traffic, the libm intrinsic calls,
the O(N) force-accumulation loop, and the worst cache-miss pattern in
`nnp_scale_acsf`.  Expected speedup 5–15× over master at moderate system
size.

---

## Branch 3 — `feature/nnp-native-spline-omp`

`src/nnp_acsf.F` = 2 491 lines.  Built directly on top of `nnp-native-spline`:
inherits the full optimisation stack (sections 2.1–2.18) and adds OpenMP
threading.  Net diff from `nnp-native-spline`: +433 / −251 lines in
`nnp_acsf.F`, zero changes to `nnp_force.F`, zero changes to environment or
spline source files.  Effectively *only* the OMP scaffolding differs.

The benchmark question this branch answers: *can we trade some MPI ranks for
OMP threads to recover scaling efficiency that pure-MPI loses past ~16 ranks?*
Layering on top of an otherwise-identical serial-optimised baseline keeps
the spline strategy, the cache-miss reductions, and every other optimisation
fixed across the MPI / MPI+OMP comparison, so the measured difference is
purely the threading overhead vs threading benefit.

### 3.1 THREADPRIVATE persistent buffers

```fortran
INTEGER, SAVE, PRIVATE :: scratch_max_sym = -1
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_sym(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_forcetmp(:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_force3tmp(:,:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: local_y(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: local_fi(:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: local_stress(:,:,:)
!$OMP THREADPRIVATE(scratch_max_sym, scratch_sym, scratch_forcetmp, scratch_force3tmp, &
!$OMP&              local_y, local_fi, local_stress)
```

Each thread keeps its own copies of:
- The same scratch buffers `nnp-native-spline` uses (`scratch_sym/forcetmp/force3tmp`)
- New per-thread reduction accumulators (`local_y`, `local_fi`, `local_stress`)
  used to defer race-prone shared-array writes until the end of the parallel loop.

These persist across calls and grow lazily (grow-to-fit) — no per-call
allocation.

The fc/dr precompute scratch (`scratch_fc1/dfc1/fc2/dfc2/dr1dx/dr2dx`)
intentionally stays **shared, not threadprivate**: one thread fills it, the
rest read it.  See §3.3.

### 3.2 Parallel-region structure

For each of the four hot paths (radial-force, angular-force, radial-no-force,
angular-no-force) the structure is:

```fortran
CALL timeset('nnp_acsf_radial', handle_sf)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP&   SHARED(neighbor, nnp, ind, i, dsymdxyz, has_stress, stress_p) &
!$OMP&   PRIVATE(j, jj, rvect1, r1, sf, m, l, s, n_symf_s, symtmp, forcetmp)
DO s = 1, nnp%rad(ind)%n_symfgrp                  ! sequential s-loop in each thread
   ! per-thread scratch grow if needed (THREADPRIVATE)
   ! reset per-thread reduction accumulators
   !$OMP DO SCHEDULE(static)
   DO j = 1, neighbor%n_rad(s)
      ! body: spline eval + accumulate into local_*  (no atomic in hot path)
      ! ATOMIC kept only for dsymdxyz(m, l, jj) writes (different jj per j)
   END DO
   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ! reduce local_y / local_fi / local_stress into shared dsymdxyz / nnp%y / stress_p
   !$OMP END CRITICAL
   !$OMP BARRIER                                  ! sync before next s-iteration
END DO
NULLIFY(symtmp); NULLIFY(forcetmp)
!$OMP END PARALLEL
CALL timestop(handle_sf)
```

One PARALLEL region wraps the entire s-loop (not one per `s`), to amortise
fork/join overhead across all symfgrps for a given centre.

### 3.3 Shared fc/dr precompute under `!$OMP SINGLE` + `!$OMP DO`

The `scratch_fc1/dfc1/fc2/dfc2/dr1dx/dr2dx` arrays (§2.11–2.12) are read by
every thread inside the `(j,k)` loop and so must be **shared**.  The
precompute pattern is:

```fortran
!$OMP SINGLE
   ! one thread does the realloc check + grow-to-fit on the SHARED arrays
!$OMP END SINGLE
! implicit barrier — alloc visible to every thread

SELECT CASE (cut_type)
CASE (nnp_cut_cos)
   !$OMP DO SCHEDULE(static)
   DO j = 1, neighbor%n_ang1(s)
      ! all threads share the precompute work — spline calls are the bulk of it
   END DO
   !$OMP END DO        ! implicit barrier publishes scratch_fc*/dr*
   IF (sg%ele(1) /= sg%ele(2)) THEN
      !$OMP DO SCHEDULE(static)
      DO k = 1, neighbor%n_ang2(s)
         ! same for n_ang2
      END DO
      !$OMP END DO
   END IF
CASE (nnp_cut_tanh)
   ...
END SELECT
```

Earlier iterations of this branch did the entire precompute under
`!$OMP SINGLE`, with one thread executing all the spline evaluations and the
rest spinning on the implicit barrier.  Splitting the realloc check (which
must be serial — only one allocator) from the precompute loops (which are
embarrassingly parallel and dominated by spline calls) gives the team useful
work during what was previously dead time.

### 3.4 Why the local-reduction + CRITICAL pattern

The naïve approach — `!$OMP ATOMIC UPDATE` on every shared write — produced
**high contention** on `dsymdxyz(m, l, i)` and `nnp%rad/ang(ind)%y(m)`,
because every `j` iteration writes to the same `i`-indexed slot and the same
`m`-indexed slot for a given `sf`.  The hardware bus lock serialised every
update, making the parallel loop slower than serial.

The fix is the local-reduction pattern:
- Accumulate per-thread partial sums in `local_fi`, `local_y`, `local_stress`
  (THREADPRIVATE — no contention).
- One `!$OMP CRITICAL` block per thread per s-iteration reduces the locals
  into the shared arrays.  CRITICAL serialises ~`n_threads` × ~10 simple
  scalar operations — negligible compared to the work saved.
- `!$OMP ATOMIC UPDATE` is kept only for `dsymdxyz(m, l, jj)` and `(:, :, kk)`
  writes, where different `j`/`k` iterations index different atoms (low
  contention) but PBC images of the same atom can occasionally collide.

### 3.5 Schedule + COLLAPSE tuning on the angular j-k loops

Two refinements on top of the §3.2 skeleton, applied to the angular
(triple) loops in both forces and no-forces variants:

| Branch of the angular loop | Original (initial port) | Tuned |
|---|---|---|
| Same-element (triangular: `k = j+1..n_ang1`) | `!$OMP DO SCHEDULE(dynamic)` | `!$OMP DO SCHEDULE(dynamic, 4)` |
| Mixed-element (rectangular: `k = 1..n_ang2`) | `!$OMP DO SCHEDULE(dynamic)` | `!$OMP DO COLLAPSE(2) SCHEDULE(dynamic, 16)` |

Triangular case can't COLLAPSE; the chunk size of 4 reduces scheduling
overhead vs the chunk-1 default of plain `dynamic`.  Rectangular case adds
COLLAPSE(2) so the `(j, k)` iteration space is fully distributed — keeps
threads busy when `n_ang1` alone is too small for the team — at the cost of
moving the j-only `rvect1`/`r1` reads inside the inner loop to satisfy
COLLAPSE's perfectly-nested-loop requirement.

### 3.6 Defensive grow-check

```fortran
IF (.NOT. ALLOCATED(scratch_sym) .OR. .NOT. ALLOCATED(local_y) .OR. &
    scratch_max_sym < n_symf_s) THEN
   ! grow all THREADPRIVATE arrays together
END IF
```

The `.NOT. ALLOCATED(local_y)` guard is non-obvious: a separate code path
(`nnp_neighbor_ensure_persistent` in serial code) pre-allocates `scratch_sym`
for the master thread only.  Without the `local_y` guard, the master thread
would skip the grow branch on first parallel entry — leaving `local_y`
unallocated for the master and segfaulting on `local_y(1:n_symf_s) = 0.0_dp`.

### 3.7 Parallelised `nnp_scale_acsf`

The hot `dsymdxyz` scaling loops (now using the precomputed scale_factor
from §2.18) are parallelised with
`!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(static)` over the `(idx, j)` pair —
each pair writes a unique `dsymdxyz(j, :, k)` slot, so no atomics are
needed.  COLLAPSE(2) handles the case where `n_active` or `n_rad/n_ang`
is small relative to the team.  Both the active-atoms branch and the
no-active-atoms (full `1..num_atoms`) fallback are parallelised.

### Performance — what it measures vs native-spline

This branch is the MPI+OMP comparison point against `nnp-native-spline`
(pure-MPI).  Since every non-OMP optimisation is now identical between the
two branches (per the May 2026 backport), the benchmark question reduces to
whether the OMP overhead in sections 3.1–3.6 is recovered by:

1. Reducing the MPI rank count (less collective overhead).
2. Spreading per-centre work across more cores when MPI ranks are limited
   by per-rank memory or by the number of NUMA domains on a node.

At small per-rank atom counts the synchronisation overhead is expected to
dominate; the question is at what system / rank ratio hybrid wins.

---

## Branch 4 — `feature/nnp-chebyshev` (active development)

Built directly on `nnp-native-spline-omp`: inherits the full optimisation
stack of Branches 2–3 (cell list, persistent buffers, OMP threading, etc.)
and replaces / extends it as below.  Motivated by the thesis profile of the
omp branch (MAQAO Fig 4.7): MPI collectives ~33%, OMP runtime 28%,
`hermite_spline_value` 6–11%, with the remaining ceiling identified as
L1→L2 cache transfer of the 480 KiB spline tables.

### 4.1 Table-free minimax polynomials (`3dbef3efb8`)

The three Hermite spline tables are removed entirely; the cutoff and
Gaussian evaluations go through Chebyshev-derived minimax polynomials
evaluated with the Clenshaw recurrence (new in `splines_methods.F`):

| Function | Replacement | Domain | Degree |
|---|---|---|---|
| `COS` (cutoff type 1) | `minimax_cos` | `[0, pi]` | 13 |
| `TANH` (cutoff type 2) | `minimax_tanh` | `[0, 1]` | 14 |
| `EXP` (radial Gaussian) | `exp_reduced` | range-reduced | 9 |

`exp_reduced` uses range reduction `exp(x) = 2^k * exp(r)` with
`r ∈ [-ln2/2, ln2/2]`, `k = NINT(x/ln2)`, and the exact `SCALE` intrinsic
for `2^k`, giving full RELATIVE precision down the Gaussian tail (a single
global fit would need degree ~35 and still be wrong in relative terms at
`exp(-50)`).  Accuracy vs libm: ≤1.3e-13 absolute, ≤1.8e-14 relative —
better than the splines replaced.  Derivatives are carried through the same
recurrence (`clenshaw_d`), so the returned slope is the exact analytic
derivative of the returned value (energy/force consistency).  Removes
~480 KiB of cache-resident tables per rank and the per-call coefficient
gathers — this attacks the L1→L2 ceiling directly.  `init_hermite_spline` /
`hermite_spline_value` remain in `splines_methods.F` for other callers but
the NNP module no longer uses them.

**Coverage note:** every NNP available locally (Morawietz, bulkH2O, all six
O'Neill NaCl committee models) uses `cutoff_type 2` (tanh), so MD validation
exercises `minimax_tanh` + `exp_reduced` only; `minimax_cos` is validated by
a standalone harness against libm.  Flagged for supervisor.

### 4.2 Shrunk per-step MPI collectives (`a9a17bc3e3`)

- **Energy:** only per-committee sums are consumed downstream, so the
  all-reduce carries `n_committee` scalars instead of the
  `(num_atoms, n_committee)` `atomic_energy` array (786 KiB → 64 B per step
  at 12k atoms × committee-of-8).  `atomic_energy` is now rank-local;
  `PRINT%ENERGIES` reads the reduced `committee_energy`.
- **Forces:** per-member forces are only consumed when BIAS is active or
  FORCES / FORCES_SIGMA printing fires.  Otherwise the committee average is
  formed locally and a single `(3, N)` array is reduced instead of
  `(3, N, n_committee)` — n_committee-fold less traffic on the collective
  the thesis identified as the scaling wall.
- `nnp_bias_sigma` no longer overwrites `committee_energy` with aligned
  values (local copy instead); the bias-energy print applies the shift
  itself.  Printed output is unchanged.

### 4.3 One-time element-sorted index maps (`286b675bc0`)

`sort / sort_inv / ele_ind / atoms / nuc_atoms` were rebuilt every MD step
with `n_ele*N` string comparisons plus one `get_ptable_info` periodic-table
lookup per atom.  They depend only on the particle kinds, which are fixed
during MD: built once on first force evaluation (per-instance
`sort_initialized` flag in `nnp_type`, safe for MIXED setups with several
NNP force_evals); only the coordinate array is refreshed per step.

### 4.4 Reverse-mode backpropagation in `nnp_gradients` (`c06fde41cd`)

The reference routine propagated full `(n_sym × n_nodes)` Jacobian blocks
(`tmp_der`) forward through every layer with per-layer DGEMMs, then kept
only the **first column** of the final block.  For a scalar-output network
the adjoint sweep gives the same dE/dG with vector-matrix products:

```
delta_L = 1
delta_{k-1} = W_k * (f'_k . delta_k)     (DGEMV per layer)
dE/dG = delta_1
```

Cost drops from O(n_sym · Σ n_{k-1}n_k) to O(Σ n_{k-1}n_k): ~12× fewer
flops on the [30,20,20,1] bulk-water network, in a kernel called
`num_atoms × n_committee` times per step.  The `tmp_der` buffers are
removed from `nnp_arc_layer_type`.  **Results are mathematically identical
but not bitwise** (different summation order) — this is the change that
moves validation from bitwise regtest equality to ε-tolerance +
branch-consistency + NVE drift.

### 4.5 FORCEINLINE on the minimax kernels (`2a05529bf4`)

`!DIR$ ATTRIBUTES FORCEINLINE` on `minimax_*` / `exp_reduced` /
`clenshaw(_d)`.  Inert without interprocedural optimisation — takes effect
under `-ipo` / `-flto`, which would let the enclosing symmetry-function
loops vectorise across the polynomial.  (PURE is not possible on the
wrappers: the optional INTENT(OUT) derivative argument is disallowed in
pure functions.)

### 4.6 Walk-bound persistent-buffer sizing (`73adb59162`)

`persistent_neighbor` was dimensioned `num_atoms*(2p+1)^3` rows per
symmetry-function group — ~330k rows ≈ ~100 MB/rank at 12,288 atoms — when
a centre's cell-list walk can only visit (walk bins × longest chain)
images, a few hundred in dense water.  The rebuild records
`max_bin_occupancy`; `nnp_prepare_cell_list_cache` derives the exact bound
(chains are static between rebuilds) and `nnp_neighbor_ensure_persistent`
now takes the capacity as an argument.  Defensive CPABORT guards on the
three append paths.  ~500× memory cut on the largest allocation.

### 4.7 Halo-filtered image table (`13cee34d90`)

The rebuild stored all `(2p+1)^3 * N` images.  Centres only live in the
primary cell, so images farther than `list_cutoff` from the primary
bounding box (per-axis conservative test) can never be paired; Verlet
drift is covered since `list_cutoff = max_cut + skin` and both partners
move ≤ skin/2 between rebuilds.  Count+fill two-pass keeps ~3N of 27N
images for the 1024-water box: ~9× less image-table memory, shorter bin
chains, cheaper rebuilds.

### 4.8 BUG FIX — missing periodic images along long axes (`50a179ea8a`)

Root cause of the May 2026 Fig. S4 "explosions" (N≥128 replicated boxes,
T runaway to 700–900 K, single-step PES jumps >0.1 Ha; N=64 cubic always
stable).  `nnp_compute_pbc_copies` answers the upstream question — image
shells needed by a per-pair `pbc()`-folding search, which is zero until
the cutoff exceeds the half-box.  The cell list enumerates images
explicitly with **no folding in the walk**, where a face-straddling pair
needs the ±1 image of its partner on every periodic axis regardless of
cutoff.  With the upstream count, every cross-face pair along any axis
longer than 2·cutoff was silently missing; when an atom wrapped across
such a face its pair set changed discontinuously → PES jump.  Fix: shift
bound for explicit enumeration is |T| ≤ 1 + cutoff/height.  Cubic boxes
with half-box < cutoff already got copies=1, which is why all N=64
validations (regtests, physics_check) were stable while every replicated
size failed.  Verified by `figS4/diag_n128_chebyshev.sh` (500-step NVT on
the exact failing configuration, master vs branch).

### 4.9 Centre-level OMP rework + same-TU minimax kernels (June 12)

Motivated by the Report 2 MAQAO ONE View compare (2000-step N64, 16 ranks,
I/O off, cpu-q-206): master 36.0 s vs chebyshev 51.1 s wall (1.42x SLOWER),
with `__kmpc_atomic_float8_add` alone at 18.4% of chebyshev samples and
OMP-runtime modules at 24.5% — in a run with OMP_NUM_THREADS=1.  An
`!$OMP ATOMIC` compiles to a locked RMW through libiomp5 whether there is
one thread or sixteen; the inner-OMP scaffolding (2 PARALLEL regions per
centre, per-group CRITICAL/BARRIER, ATOMIC on every dsymdxyz accumulation)
taxed every configuration, including all pure-MPI production runs.  The
same profile showed `exp_reduced` (17.5%) + `minimax_tanh` (5.3%) costing
more than master's libm calls (~11%): FORCEINLINE is inert without
-ipo/-flto, so each evaluation was a scalar CALL into splines_methods.

Three coordinated changes (one commit; all signatures backward-compatible,
`motion/helium_interactions.F` untouched and verified to still compile):

1. **Inner-OMP scaffolding deleted** from `nnp_calc_acsf` /
   `nnp_scale_acsf`: no PARALLEL/DO/SINGLE/CRITICAL/BARRIER/ATOMIC, no
   `local_y`/`local_fi`/`local_stress` fold buffers — accumulation into
   `dsymdxyz`/`stress`/`y` is direct, in the same summation order the
   1-thread path produced before (bitwise-equal modulo kernel inlining /
   FMA contraction).
2. **Centre-level threading** in `nnp_calc_energy_force`: one
   `!$OMP PARALLEL DO SCHEDULE(DYNAMIC)` over the rank's atom slice.
   Per-thread (THREADPRIVATE) state: `persistent_neighbor` + scratch +
   fc-precompute arrays + new `ws_y_rad`/`ws_y_ang` accumulators (replace
   the shared `nnp%rad/ang(ind)%y`, now unused) in nnp_acsf;
   `dsymdxyz`/`stress`/active set/`denergydsym`, a deep-cloned network
   state `my_arc` (~100 kB/thread; `layer` is a POINTER component so the
   clone is copied depth-2 explicitly), and `my_force`/`my_cstress`
   accumulators in nnp_force, folded under one named CRITICAL per force
   evaluation.  Sizing protocol: `nnp_neighbor_ensure_persistent` (serial,
   from prepare) now only PUBLISHES req_* sizes; each thread grows its own
   buffers in `nnp_acsf_ensure_thread_buffers` (integer compares on the
   fast path).  `nnp_calc_acsf` gained OPTIONAL `arc` and `expol` args:
   with them it is reentrant; without them (helium path) it falls back to
   the shared `nnp%arc(ind)` / `nnp%output_expol`.
   Memory: per-thread dsymdxyz is the dominant clone (max_nsym x 3 x N x 8B
   ~ 9 MB/thread at N=12288); fine for <=16 threads/rank on 256 GB nodes.
   Multi-threaded forces are bitwise-reproducible only for a fixed thread
   count (thread fold order); 1-thread runs stay deterministic.
3. **Minimax kernels moved into nnp_acsf.F** (deleted from
   splines_methods.F, which reverts to upstream + the Hermite routines).
   Same translation unit -> FORCEINLINE is honoured at -O2 without IPO.
   The OPTIONAL `y1` argument is gone: value-only PURE functions
   (`minimax_cos/minimax_tanh/exp_reduced`) plus a branch-free
   `minimax_cos_d` subroutine, so call sites carry no PRESENT() checks and
   the fc-precompute loops are vectorisable straight-line code.

Expected effect: removes the ~25% OMP-runtime tax at every size in every
pure-MPI run (the entire small-N deficit vs master at N64–N256 traced to
it), plus whatever the kernel inlining recovers of the ~25% kernel share;
OMP thread scaling should improve qualitatively (no per-pair atomics, no
per-centre fork/join — one region per force evaluation).
**Must re-run after rebuild: physics_check (force diff vs master + NVE
drift), NNP regtests, then the chebyshev-only scaling suite and the MAQAO
profile to quantify.**  Fig. S4 production data stays valid if the force
diff is still at ε (observables are functions of the forces); only quoted
t/step numbers become stale.

### Planned (scoped, not yet implemented)

- O(1) element→group dispatch in `nnp_compute_neighbors` via the
  `ele_grp_list` / `ele_grp_count` tables (built at init since the May
  backport but not yet consumed; matters most for the 4-element NaCl
  models).
- Committee-loop force assembly currently re-reads `dsymdxyz`
  `n_committee` times (hoist `denergydsym` for all members, single pass or
  small DGEMM).
- BLAS-3 batching of the per-atom forward pass (thesis §5.2 future work).
- Batched (array-at-a-time) evaluation of `exp_reduced` across the radial
  sf loop if MAQAO still shows it hot after the same-TU inlining.

### Validation status

NNP/regtest-1 passes on the minimax commit.  Physics gate for the whole
pass: `~/cp2k-benchmarks/validation/physics_check/run_check.sh` (sbatch)
— single-point master-vs-branch force diff + 2 ps NVE conserved-quantity
drift on H2O-64, both binaries.  4.4 changes summation order, so expect
agreement to ε, not bitwise.

---

## Cross-branch summary

### Spline strategy

| Branch | Where | Per-eval cost | Notes |
|---|---|---|---|
| `master` | libm (`COS`, `TANH`, `EXP`) | ~30–80 cycles | Direct intrinsic call; blocks vectorisation |
| `nnp-native-spline` | `splines_methods.F` Horner form | ~3 MUL + 3 ADD value, ~3 MUL + 3 ADD derivative | Public, available to any CP2K caller; coefficients pre-collected at init |
| `nnp-native-spline-omp` | `splines_methods.F` Horner form | same as native-spline | Identical spline strategy — only OMP layered on top |
| `nnp-chebyshev` | `splines_methods.F` minimax + Clenshaw | deg 9–14 polynomial, no table loads | Table-free: no index search, no coefficient gather, no 480 KiB cache footprint; exp via 2^k range reduction |

### Threading

| Branch | Granularity | Sync points |
|---|---|---|
| `master` | Single-threaded; MPI atom decomposition | None within a call |
| `nnp-native-spline` | Single-threaded; MPI atom decomposition | None within a call |
| `nnp-native-spline-omp` | OMP DO over j-loop per s-iteration in `nnp_calc_acsf`; OMP DO COLLAPSE(2) over (idx, j) in `nnp_scale_acsf` | `!$OMP DO`, `!$OMP SINGLE`, `!$OMP CRITICAL`, `!$OMP BARRIER`, `!$OMP ATOMIC UPDATE` |
| `nnp-chebyshev` | same as native-spline-omp (centre-level threading scoped, not yet implemented) | same as native-spline-omp; smaller MPI collectives (§4.2) |

### Key files touched (line counts, current state)

| File | master | native-spline | native-spline-omp | nnp-chebyshev |
|---|---:|---:|---:|---:|
| `nnp_acsf.F` | 1 414 | 2 309 | 2 491 | **2 522** |
| `nnp_force.F` | 674 | 710 | 710 | 761 |
| `nnp_environment_types.F` | 592 | 635 | 635 | 630 |
| `nnp_environment.F` | 798 | 832 | 832 | 830 |
| `splines_methods.F` | 244 | 339 | 339 | 519 |
| `splines_types.F` | (unchanged) | +8 lines | (same as native-spline) | (same) |
| `nnp_model.F` | (unchanged) | (unchanged) | (unchanged) | 257 (reverse-mode) |
| `common/splines.F` | (unchanged) | (unchanged) | (unchanged) | (unchanged) |

`nnp_acsf.F` is the only file that meaningfully differs between native-spline
and native-spline-omp — the +182 net lines on OMP are the threading
scaffolding (PARALLEL regions, CRITICAL reductions, COLLAPSE/SCHEDULE
directives, defensive grow-checks).  `nnp_force.F` has a 2-line cosmetic
difference between the two; everything else is identical.  `nnp-chebyshev`
additionally touches `nnp_model.F` (first branch to do so) and grows
`splines_methods.F` by the minimax/Clenshaw kernels and coefficient tables.

### Optimisation matrix (verified against the actual source state)

| Optimisation | master | native-spline | native-spline-omp |
|---|:-:|:-:|:-:|
| Cell-list / Verlet cache | — | ✓ | ✓ |
| Persistent neighbour buffer | — | ✓ | ✓ |
| Module SAVE scratch + `POINTER, CONTIGUOUS` aliases | — | ✓ | ✓ (THREADPRIVATE) |
| `kk` index correction (hetero-angular) | — | ✓ | ✓ |
| Native CP2K Hermite splines | — | ✓ | ✓ |
| Horner-form spline evaluator | — | ✓ | ✓ |
| Division → reciprocal-multiply | — | ✓ | ✓ |
| Precomputed `cutoff_sq` / `inv_cut` / `pi_over_cut` | — | ✓ | ✓ |
| `zeta_int` / `zeta_is_int` precompute | — | ✓ | ✓ |
| `r_sum_sq`, `g_sq_inv` hoists | — | ✓ | ✓ |
| MPI allgather → analytical `istart` | — | ✓ | ✓ |
| O(n_active) force-accumulation loop | — | ✓ | ✓ |
| Selective zeroing of `dsymdxyz` | — | ✓ | ✓ |
| Hoist `f_c(r1)` / `f_c(r2)` from k-loop | — | ✓ | ✓ |
| Hoist `dr1dx` / `dr2dx` from k-loop | — | ✓ | ✓ |
| Skip SQRT via `cutoff_sq` in angular | — | ✓ | ✓ |
| `ASSOCIATE` aliases for hot fields | — | ✓ | ✓ |
| Local `cut_type` instead of `nnp%cut_type` | — | ✓ | ✓ |
| Pre-computed `scale_factor` / `inv_sigma_factor` | — | ✓ | ✓ |
| Sorted `active_atoms` (in-place heap sort) | — | ✓ | ✓ |
| OpenMP threading in `nnp_calc_acsf` (4 hot loops) | — | — | ✓ |
| `!$OMP SINGLE` realloc + `!$OMP DO` fc/dr precompute | — | — | ✓ |
| COLLAPSE(2) on mixed-element angular j-k loop | — | — | ✓ |
| `SCHEDULE(dynamic, chunk)` tuning on angular | — | — | ✓ |
| OpenMP threading in `nnp_scale_acsf` (COLLAPSE(2)) | — | — | ✓ |
| THREADPRIVATE persistent scratch + reduction accumulators | — | — | ✓ |
| `has_stress` / `stress_p` hoist (skip `PRESENT` checks in hot loop) | — | — | ✓ |

`feature/nnp-chebyshev` inherits every ✓ of `nnp-native-spline-omp` except
the Hermite-spline rows (replaced by the table-free minimax kernels of
§4.1) and adds the §4.2–4.7 optimisations on top: shrunk MPI collectives,
one-time element sort, reverse-mode `nnp_gradients`, FORCEINLINE
directives, walk-bound neighbour-buffer sizing, and the halo-filtered
image table.  Commits `3dbef3efb8` → `13cee34d90` (June 2026).

### Verification of optimisation presence

The "✓" marks above were confirmed by grep counts on each branch's source
tree (as of May 2026):

| Symbol | master | native-spline | native-spline-omp |
|---|:-:|:-:|:-:|
| `nnp_sort_active_atoms` | 0 | 3 | 3 |
| `scale_factor` field references | 0 | 22 | 22 |
| `inv_sigma_factor` field references | 0 | 20 | 20 |
| `hermite_spline_value` call sites | 0 | 29 | 28 |
| `spl%coef` (Horner coefficient table) | 0 | 14 | 14 |
| `!$OMP THREADPRIVATE` | 0 | 0 | 4 |
| `!$OMP PARALLEL` | 0 | 0 | 12 |
| `!$OMP DO` | 0 | 0 | 17 |
| `!$OMP CRITICAL` | 0 | 0 | 5 |
| `COLLAPSE(2)` | 0 | 0 | 12 |

The single-call-site difference in `hermite_spline_value` (29 vs 28) between
native-spline and native-spline-omp is from a code-path consolidation in the
OMP rework of `nnp_calc_ang`; mathematically the same evaluations happen.

---

## How to read the benchmark output

The harness produces one CSV per branch per scaling type (size or core),
under:

```
/rds/user/crm98/hpc-work/cp2k-benchmarks/results/
  ├── cp2k_master/NNP/...
  ├── cp2k_feature_native_spline/NNP/...
  └── cp2k_feature_native_spline_omp/NNP/...
```

When comparing:

- **`master` vs `native-spline`** measures the *combined* effect of all the
  inner-loop optimisations, neighbour-walk optimisations, native-spline
  replacement of libm, and cache-miss reductions.  This is the largest jump
  in the benchmark suite.  Expected speedup 5–15× at moderate system size,
  larger at larger N (where malloc traffic dominated master).
- **`native-spline` vs `native-spline-omp`** is the pure-MPI vs hybrid-MPI+OMP
  comparison with every non-OMP optimisation held fixed.  The measured
  difference is the net of OMP synchronisation overhead (sections 3.1–3.6)
  against threaded speedup of the j-loops in `nnp_calc_acsf` and the
  `(idx, j)` loops in `nnp_scale_acsf`.  At small per-rank atom counts the
  synchronisation overhead is expected to dominate; the question is at what
  system / rank ratio hybrid wins.

Treating each branch as a separable benchmark makes it possible to point at
which optimisation strategy is responsible for which performance change, and
avoids the trap of "is this faster because of A, B, or both?" that a single
mega-branch would create.
