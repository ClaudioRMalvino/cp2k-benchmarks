# NNP Implementation Changes — Cross-Branch Reference

This document records the implementation differences between the four NNP-related
branches of CP2K maintained in this project, why each branch exists, and what the
benchmark suite is actually measuring.

```
Branches under test (all rooted at upstream master):

   upstream/master                       ← baseline
       │
       ├── feature/nnp-verlet-cells      ← cell-list neighbour walk +
       │                                    persistent buffers +
       │                                    local Hermite Horner spline cache +
       │                                    division/SQRT eliminations
       │
       ├── feature/nnp-native-spline     ← same neighbour walk +
       │                                    native CP2K splines (splines_methods.F) +
       │                                    additional micro-optimisations
       │                                    (fc/dfc/dr hoisting, sqrt skip,
       │                                     ASSOCIATE, cut_type local)
       │
       └── feature/nnp-native-spline-omp ← native-spline base + OpenMP threading
                                            (THREADPRIVATE scratch, OMP PARALLEL DO
                                             with local-reduction + CRITICAL) +
                                            cache-miss reductions in nnp_scale_acsf
                                            (precomputed scaling factors, sorted
                                             active_atoms via in-place heap sort)
```

---

## Why three feature branches, all benchmarked together

Each feature branch represents a **different optimisation strategy** layered on
top of the same neighbour-walk improvements.  Benchmarking them in isolation
matters because the strategies interact:

| Concern | Affected branches | Reason for separate measurement |
|---|---|---|
| Spline evaluation cost (Horner vs basis form) | `nnp-verlet-cells` (Horner) vs `nnp-native-spline` / `nnp-native-spline-omp` (basis) | Native CP2K splines use basis-function form (~22 MULs); local Horner cache uses ~5 MULs.  We need to know the absolute cost of going "native" before recommending it as the merge candidate |
| Synchronisation overhead vs work granularity | `nnp-native-spline-omp` | Adding OpenMP on top of `nnp-native-spline` rather than `nnp-verlet-cells` keeps the spline strategy constant across the MPI / MPI+OMP comparison — the measured difference is *only* the threading overhead, not the spline cost.  This is the right pairing for "does hybrid MPI+OMP recover scaling that pure MPI loses past ~16 ranks?" |
| Maintenance burden of new modules | `nnp-native-spline` / `nnp-native-spline-omp` add public symbols to `splines_methods.F` | This change touches CP2K core and benefits the wider codebase; we want to see the real-world cost in the NNP context |
| Hyperthread contention | `nnp-native-spline-omp` | Cerberus has 32 physical cores × 2 hyperthreads.  OMP behaviour is sensitive to physical-vs-logical core binding; pure-MPI baselines on `nnp-native-spline` give us the same-spline-strategy reference for "what should be possible" |

The benchmark harness builds **all four** binaries (master + three feature
branches) into separate per-branch directories so `LD_LIBRARY_PATH` cannot
accidentally pick up the wrong `libcp2k.so`.  Each branch is then run through
the same size-scaling and core-scaling sweeps.

---

## Branch 1 — `upstream/master` (baseline)

`src/nnp_acsf.F` ≈ 1 414 lines.  Reference implementation as merged into CP2K
upstream.

Performance characteristics relevant to our work:

- **Neighbour search**: per centre atom calls `nnp_compute_pbc_copies` then
  loops over all `(j, c1, c2, c3)` periodic images, calling `pbc()` on every
  displacement.  Cost: O(N · (2p+1)³) per centre, plus ~9 ALLOCATE/DEALLOCATE
  calls per centre for buffer setup.
- **Spline evaluation**: the cutoff function uses direct `COS(π·r/rc)` and
  `TANH(1−r/rc)` libm calls.  The radial Gaussian uses `EXP(−η(r−rs)²)`.
  Each is ~30–80 cycles per evaluation.
- **Force accumulation**: outer atom loop in `nnp_force.F` iterates over all
  `nnp%num_atoms` per centre even though only ~k actually contribute.
- **MPI atom assignment**: `allgather(mecalc, allcalc)` collective + prefix
  sum to compute each rank's `istart`.

The remaining branches treat this as the truth-source for correctness and the
"slow but always-right" reference for benchmarks.

---

## Branch 2 — `feature/nnp-verlet-cells`

`src/nnp_acsf.F` ≈ 2 193 lines.  This is the foundational performance branch.
Its changes are mostly about avoiding work the upstream couldn't avoid; almost
everything that follows in `nnp-native-spline` and `nnp-native-spline-omp`
inherits from this branch.

### 2.1 New module-level types

#### `nnp_cell_list_cache_type`
A Verlet-skin / cell-list hybrid cache.  Each `(atom j, integer PBC shift)`
pair becomes one *image entry* stored explicitly, so the inner neighbour-query
loop contains no `pbc()` calls.  The coarse grid has `bin_width ≈ list_cutoff`
(one bin-span walk covers 27 bins).

| Field | Purpose |
|---|---|
| `coord_primary(3,N)` | PBC-wrapped primary atom positions, refreshed every step |
| `ref_coord_primary(3,N)` | Snapshot at last rebuild (Verlet reference) |
| `image_atom(n_img)` | Source atom index for each image entry |
| `image_translation(3,n_img)` | Cartesian translation vector for each image |
| `head(n_cells)` / `next(n_img)` | Linked-list cell grid (LIFO insertion) |
| `verlet_skin` | Skin in Å; `list_cutoff = max_cut + skin` |
| `exact_cutoff`, `hmat`, `perd` | Change-detection fields for rebuild triggers |

A `SAVE` instance `cell_list_cache` lives at module level and persists across
force evaluations.

#### `nnp_spline_cache_type`
Stores Hermite cubic spline coefficients `(a, b, c, d)` for **Horner-method**
evaluation of `EXP`, `COS`, and `TANH`.  Coefficients are packed as
`coef(4, n_points)` so each evaluation reads four consecutive doubles.  A
`SAVE` instance `spline_cache` is initialised once in `nnp_init_acsf_groups`.

This local cache is the key reason this branch is so fast: per-evaluation cost
is **~5 MULs + 4 ADDs**, vs ~30–80 cycles for the libm intrinsic call.

### 2.2 New module-level SAVE variables

```fortran
TYPE(nnp_cell_list_cache_type), SAVE, PRIVATE :: cell_list_cache
TYPE(nnp_spline_cache_type),    SAVE, PRIVATE :: spline_cache

TYPE(nnp_neighbor_type), SAVE, PRIVATE, TARGET :: persistent_neighbor
INTEGER, SAVE, PRIVATE :: persistent_n_capacity = -1
INTEGER, SAVE, PRIVATE :: persistent_n_rad_grp  = -1
INTEGER, SAVE, PRIVATE :: persistent_n_ang_grp  = -1

REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_sym(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_forcetmp(:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_force3tmp(:,:,:)
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

! verlet-cells
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

`CONTIGUOUS` matters: without it the compiler treats the pointer as potentially
non-contiguous and potentially aliasing other `TARGET`/`POINTER` variables in
scope, blocking vectorisation of the inner `sf`/`l` loops.

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
! upstream (buggy)
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
correct.  Verlet-cells fixes the index and removes the redundant inner reads.

### 2.7 Spline replacement of `EXP`, `COS`, `TANH`

The expensive intrinsic calls in `nnp_calc_sym_rad` and `nnp_calc_ang` are
replaced by the local Horner-form Hermite cubic cache:

| Intrinsic | Spline range | Rationale |
|---|---|---|
| `EXP(x)` | `[−50, 0]` | Argument is always ≤ 0 (Gaussian decay) |
| `COS(π·r/rc)` | `[0, π]` | Cutoff arg ∈ `[0, 1)` so `tmp ∈ [0, π]` only — doubles grid resolution |
| `TANH(1−r/rc)` | `[0, 1]` | Argument lives in `[0, 1)` — 20× finer resolution |

5 000 grid points each.  For `EXP` specifically, the spline value (~`EXP`) is
used directly in the chain rule for the force rather than evaluating the
formal spline derivative — same answer, slightly more accurate (O(h⁴) value
vs O(h³) derivative).

### 2.8 Division elimination

```fortran
! upstream
drdx(:) = rvect(:)/r

! verlet-cells
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
SF loop on every call.  Verlet-cells precomputes `zeta_int(:)` and
`zeta_is_int(:)` once at init, eliminating the runtime NINT and FP-equality
test from every SF evaluation.

### 2.11 `nnp_force.F` — accompanying changes

- **Allgather elimination**: `mecalc` and `istart` derived analytically from
  `(num_atoms, num_pe, mepos)`; no MPI collective needed.
- **SAVE scratch arrays**: `dsymdxyz`, `stress`, `denergydsym`, `active_atoms`,
  `active_flag` allocated once at first force eval, reused across steps.
- **Hoisted allocation**: `dsymdxyz` sized to `max_nsym` (max over all
  elements) at force-eval entry rather than re-allocating per centre atom.
- **Selective zeroing**: `dsymdxyz(:, :, k) = 0` only for atoms in the active
  list, not the whole `(:,:,:)`.
- **O(n_active) force accumulation loop**: only walks the touched atoms.

### 2.12 `nnp_environment_types.F` — new fields

Per `symfgrp` (radial and angular, identical):
```fortran
REAL(KIND=dp) :: cutoff_sq   = -1.0_dp   ! eliminates SQRT in neighbour filter
REAL(KIND=dp) :: inv_cut     = -1.0_dp   ! eliminates / in cutoff arg
REAL(KIND=dp) :: pi_over_cut = -1.0_dp   ! eliminates / in cutoff arg
```

Per `nnp_rad_type` and `nnp_ang_type`:
```fortran
INTEGER, ALLOCATABLE :: ele_grp_count(:)    ! O(1) element→group dispatch (infrastructure)
INTEGER, ALLOCATABLE :: ele_grp_list(:,:)
```

Per `nnp_ang_type`:
```fortran
INTEGER, ALLOCATABLE :: zeta_int(:)        ! NINT(zeta), precomputed
LOGICAL, ALLOCATABLE :: zeta_is_int(:)
```

### Performance — observed

Verlet-cells is the fastest single-thread implementation we have.  In our
36-rank pure-MPI runs on 64–4 096 H₂O it produces the timing baseline that
the other two feature branches are measured against.

---

## Branch 3 — `feature/nnp-native-spline`

`src/nnp_acsf.F` ≈ 2 279 lines.  Same neighbour walk as verlet-cells but with
two distinguishing changes:

1. **Splines moved to CP2K-native `splines_methods.F`** (no module-local
   `nnp_spline_cache_type`).
2. **Additional inner-loop micro-optimisations** that are absent or
   incomplete on verlet-cells.

### 3.1 New public symbols in `splines_methods.F` (+69 lines)

```fortran
PUBLIC :: init_hermite_spline   ! init Hermite cubic spline (stores h*dy/dx in y2)
PUBLIC :: hermite_spline_value  ! evaluate Hermite cubic spline (+ optional derivative)
```

The new initialiser populates the standard `spline_data_type` struct, storing
`h·dy/dx` in `y2` so the evaluator needs no extra division per call.

The evaluator uses **Hermite basis-function form**:

```fortran
val = (2t³ − 3t² + 1)·ylo +
      (t³ − 2t² + t)·m0  +
      (−2t³ + 3t²)·yhi   +
      (t³ − t²)·m1
```

This is mathematically equivalent to the Horner form used by verlet-cells but
requires **~22 MULs + 12 ADDs per evaluation**, vs ~5 MULs + 4 ADDs in Horner
form.  This is the deliberate trade-off being measured: native splines mean
the optimisation lives where any future CP2K caller can use it, but each call
costs ~4× more than the local cache.

### 3.2 NNP usage of native splines

```fortran
USE splines_methods, ONLY: hermite_spline_value, init_hermite_spline, init_splinexy
USE splines_types,   ONLY: spline_data_type

TYPE(spline_data_type), SAVE, PRIVATE :: spl_exp, spl_cos, spl_tanh
LOGICAL,                SAVE, PRIVATE :: spline_initialized = .FALSE.
```

Call sites in `nnp_calc_rad`/`nnp_calc_ang` (and the fc precomputation block)
use `hermite_spline_value(spl_cos, x)` etc.  Three splines, same `[−50, 0]` /
`[0, π]` / `[0, 1]` ranges as verlet-cells.

### 3.3 Hoist `f_c(r1)` / `f_c(r2)` out of the `k`-loop

The angular routine reads the cutoff function for `r1` and `r2` independently
of `k`, but verlet-cells still computed them inside the `(j,k)` body.
Native-spline precomputes them per neighbour:

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

### 3.4 Hoist `dr1dx` / `dr2dx` out of the `k`-loop

Same idea, applied to the radial unit vectors.  Two more SAVE scratch arrays:

```fortran
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: scratch_dr1dx(:, :)  ! (3, n)
REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE :: scratch_dr2dx(:, :)
```

`nnp_calc_ang` accepts them via new `dr1dx_in` / `dr2dx_in` optionals.  Per
`(j,k)` pair this saves two divisions and six multiplies (the second `1/r2`
unit-vector computation).

### 3.5 Skip the SQRT in the `r3` cutoff check

In verlet-cells the inner angular loop body was:
```fortran
rvect3 = rvect2 - rvect1
r3 = NORM2(rvect3(:))                          ! always pays a sqrt
IF (r3 < sg%cutoff) THEN
   ...
```
Native-spline replaces it with:
```fortran
rvect3 = rvect2 - rvect1
r3sq = rvect3(1)**2 + rvect3(2)**2 + rvect3(3)**2
IF (r3sq < sg%cutoff_sq) THEN
   r3 = SQRT(r3sq)                             ! only on accepted pairs
   ...
```
For typical NNP cutoffs ~50–75 % of `(j,k)` pairs miss the cutoff, so the
sqrt is skipped for most of them.

### 3.6 `ASSOCIATE` aliases for hot derived-type fields

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

### 3.7 Local `cut_type`

```fortran
cut_type = nnp%cut_type   ! read once at top of nnp_calc_acsf (and rad/ang)
SELECT CASE (cut_type)    ! register-local instead of derived-type field walk
```

Trivial change individually; meaningful in aggregate because the SELECT CASE
appears in the fc precomputation, in `nnp_calc_rad`, and in `nnp_calc_ang` —
all hot.

### Performance — expected vs observed

Native-spline pays ~4× more per spline evaluation than verlet-cells (basis
form vs Horner) but avoids many spline calls outright via §3.3, §3.4, §3.5.
Whether the trade is profitable is exactly the question the benchmark answers.

---

## Branch 4 — `feature/nnp-native-spline-omp`

`src/nnp_acsf.F` ≈ 2 491 lines.  Built directly on top of `nnp-native-spline`:
inherits the full native-spline optimisation stack (sections 3.1–3.7) and
adds two distinguishable changes:

1. **OpenMP threading** of the four hot loops in `nnp_calc_acsf` and the
   `dsymdxyz` scaling loops in `nnp_scale_acsf`, using the THREADPRIVATE
   scratch + local-reduction + CRITICAL pattern.
2. **Cache-miss reductions in `nnp_scale_acsf`** — precomputed per-symfunc
   scaling factors (kills the inner divide), and a sort of `active_atoms`
   per centre so subsequent walks of `dsymdxyz(:, :, k)` hit memory in
   ascending `k` order.

The original goal was to answer: *can we trade some MPI ranks for OMP threads
to recover scaling efficiency that pure-MPI loses past ~16 ranks?*  Layering
on top of `nnp-native-spline` (rather than `nnp-verlet-cells`) keeps the
spline strategy fixed in the MPI / MPI+OMP comparison, so the measured
difference is purely the threading overhead.

### 4.1 THREADPRIVATE persistent buffers

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

These persist across calls and grow lazily (grow-to-fit) — no per-call allocation.

The fc/dr precompute scratch (`scratch_fc1/dfc1/fc2/dfc2/dr1dx/dr2dx` from
section 3.3–3.4) intentionally stays **shared, not threadprivate**: one thread
fills it, the rest read it.  See §4.3.

### 4.2 Parallel-region structure

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

### 4.3 Shared fc/dr precompute under `!$OMP SINGLE` + `!$OMP DO`

The `scratch_fc1/dfc1/fc2/dfc2/dr1dx/dr2dx` arrays added in `nnp-native-spline`
(§3.3, §3.4) are read by every thread inside the `(j,k)` loop and so must be
**shared**.  The precompute pattern is:

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

### 4.4 Why the local-reduction + CRITICAL pattern

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

### 4.5 Schedule + COLLAPSE tuning on the angular j-k loops

Two refinements on top of the §4.2 skeleton, applied to the angular
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

### 4.6 Defensive grow-check

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

### 4.7 Cache-miss reductions in `nnp_scale_acsf`

Cachegrind on the pure-MPI native-spline branch flagged `nnp_scale_acsf` as
the largest contributor to L1 / L2 misses.  Two fixes ship on this branch:

#### Pre-computed scaling factors

The original inner body recomputed `(scmax − scmin)/(loc_max(j) − loc_min(j))`
for every `(centre, sym, neighbour)` triple — pulling `loc_max(j)`,
`loc_min(j)`, and either `sigma(j)` or `loc_av(j)` through the cache
alongside the hot `dsymdxyz` traffic.  All four values are static for the
whole run.  Two new per-element arrays in `nnp_environment_types.F`
(`nnp_acsf_rad_type` / `nnp_acsf_ang_type`):

```fortran
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: scale_factor      ! (scmax-scmin) / (loc_max-loc_min)
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: inv_sigma_factor  ! (scmax-scmin) / sigma
```

are filled once in `nnp_init_model` immediately after the scaling-data file
is parsed.  The hot path becomes:

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

### 4.8 Parallelised `nnp_scale_acsf`

The hot `dsymdxyz` scaling loops are parallelised with
`!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(static)` over the `(idx, j)` pair —
each pair writes a unique `dsymdxyz(j, :, k)` slot, so no atomics are
needed.  COLLAPSE(2) handles the case where `n_active` or `n_rad/n_ang`
is small relative to the team.  Both the active-atoms branch and the
no-active-atoms (full `1..num_atoms`) fallback are parallelised.

### Performance — what it measures

This branch is the MPI+OMP comparison point against `nnp-native-spline`
(pure-MPI).  The benchmark question is whether the OMP overhead in
sections 4.2–4.6 is recovered by:

1. Reducing the MPI rank count (less collective overhead).
2. Spreading per-centre work across more cores when MPI ranks are limited
   by per-rank memory or by the number of NUMA domains on a node.
3. The cache-miss savings of §4.7, which are also visible in the pure-MPI
   path but matter more in the OMP path because the threadprivate
   reduction adds memory pressure.

Earlier OMP experiments at this granularity (the deprecated
`nnp-openmp-implementation`, layered on `nnp-verlet-cells`) were a net loss
for small systems — the j-loop body was too short to amortise the
synchronisation overhead.  The native-spline-omp branch retests that
hypothesis with a more expensive per-iteration body (basis-form spline
evaluation) and the cache-miss fixes; whether the threading wins is exactly
what the size-scaling and core-scaling sweeps answer.

---

## Cross-branch summary

### Spline strategy

| Branch | Where | Per-eval cost | Notes |
|---|---|---|---|
| `master` | libm | ~30–80 cycles | Direct intrinsic call |
| `nnp-verlet-cells` | local `nnp_spline_cache_type` Horner form | ~5 MULs + 4 ADDs | Module-local; doesn't extend the public framework |
| `nnp-native-spline` | `splines_methods.F` Hermite basis form | ~22 MULs + 12 ADDs | Public, available to any CP2K caller; partly compensated by §3.3–3.5 hoisting |
| `nnp-native-spline-omp` | `splines_methods.F` Hermite basis form | ~22 MULs + 12 ADDs | Same as `nnp-native-spline` (OMP layered on top — keeps spline strategy fixed across the MPI/OMP comparison) |

### Threading

| Branch | Granularity | Sync points |
|---|---|---|
| `master` | Single-threaded; MPI atom decomposition | None within a call |
| `nnp-verlet-cells` | Single-threaded; MPI atom decomposition | None within a call |
| `nnp-native-spline` | Single-threaded; MPI atom decomposition | None within a call |
| `nnp-native-spline-omp` | OMP DO over j-loop per s-iteration in `nnp_calc_acsf`; OMP DO COLLAPSE(2) over (idx, j) in `nnp_scale_acsf` | `!$OMP DO`, `!$OMP SINGLE`, `!$OMP CRITICAL`, `!$OMP BARRIER`, `!$OMP ATOMIC UPDATE` |

### Key files touched

| File | master | verlet-cells | native-spline | native-spline-omp |
|---|---:|---:|---:|---:|
| `nnp_acsf.F` | 1 414 | 2 193 | 2 279 | **2 491** |
| `nnp_force.F` | 674 | 710 | 710 | 710 |
| `nnp_environment_types.F` | 592 | 618 | 618 | **635** |
| `nnp_environment.F` | 798 | 800 | 800 | **832** |
| `splines_methods.F` | 244 | 244 | **313** | 313 |

`nnp-native-spline-omp` inherits the `splines_methods.F` extension from
`nnp-native-spline` (additive — two new public functions).  The bumps to
`nnp_environment_types.F` and `nnp_environment.F` are the new
`scale_factor` / `inv_sigma_factor` per-element fields and the once-at-init
fill in `nnp_init_model` (§4.7).

### Optimisation matrix

| Optimisation | master | verlet-cells | native-spline | native-spline-omp |
|---|:-:|:-:|:-:|:-:|
| Cell-list / Verlet cache | — | ✓ | ✓ | ✓ |
| Persistent neighbour buffer | — | ✓ | ✓ | ✓ |
| Module SAVE scratch + `POINTER, CONTIGUOUS` aliases | — | ✓ | ✓ | ✓ THREADPRIVATE |
| `kk` index correction (hetero-angular) | — | ✓ | ✓ | ✓ |
| Hermite spline for EXP/COS/TANH | — | ✓ Horner | ✓ basis | ✓ basis |
| Native CP2K splines (`splines_methods.F`) | — | — | ✓ | ✓ |
| Division → reciprocal-multiply | — | ✓ | ✓ | ✓ |
| Precomputed `cutoff_sq` / `inv_cut` / `pi_over_cut` | — | ✓ | ✓ | ✓ |
| `zeta_int` / `zeta_is_int` precompute | — | ✓ | ✓ | ✓ |
| `r_sum_sq`, `g_sq_inv` hoists | — | ✓ | ✓ | ✓ |
| MPI allgather → analytical `istart` | — | ✓ | ✓ | ✓ |
| O(n_active) force-accumulation loop | — | ✓ | ✓ | ✓ |
| Selective zeroing of `dsymdxyz` | — | ✓ | ✓ | ✓ |
| Hoist `f_c(r1)` / `f_c(r2)` from k-loop | — | partial | ✓ | ✓ |
| Hoist `dr1dx` / `dr2dx` from k-loop | — | — | ✓ | ✓ |
| Skip SQRT via `cutoff_sq` in angular | — | — | ✓ | ✓ |
| `ASSOCIATE` aliases for hot fields | — | — | ✓ | ✓ |
| Local `cut_type` instead of `nnp%cut_type` | — | — | ✓ | ✓ |
| OpenMP threading in `nnp_calc_acsf` (4 hot loops) | — | — | — | ✓ |
| `!$OMP SINGLE` realloc + `!$OMP DO` fc/dr precompute | — | — | — | ✓ |
| COLLAPSE(2) on mixed-element angular j-k loop | — | — | — | ✓ |
| `SCHEDULE(dynamic, chunk)` tuning on angular | — | — | — | ✓ |
| OpenMP threading in `nnp_scale_acsf` (COLLAPSE(2)) | — | — | — | ✓ |
| Pre-computed `scale_factor` / `inv_sigma_factor` | — | — | — | ✓ |
| Sorted `active_atoms` (in-place heap sort) | — | — | — | ✓ |
| `has_stress` / `stress_p` hoist (skip `PRESENT` checks in hot loop) | — | — | — | ✓ |

---

## How to read the benchmark output

The harness produces one CSV per branch per scaling type (size or core),
under:

```
/local/data/public/crm98/cp2k-benchmarks/results/
  ├── cp2k_master/NNP/...
  ├── cp2k_feature_verlet_cells/NNP/...
  ├── cp2k_feature_native_spline/NNP/...
  └── cp2k_feature_native_spline_omp/NNP/...
```

When comparing:

- **`master` vs `verlet-cells`** measures the *combined* effect of all the
  inner-loop and neighbour-walk optimisations.  This is the largest jump.
- **`verlet-cells` vs `native-spline`** measures the cost of moving the
  spline cache into the public CP2K framework, partially offset by the
  additional inner-loop micro-optimisations on this branch.  Expected:
  comparable to verlet-cells, sometimes slightly slower due to basis-form
  evaluator.
- **`native-spline` vs `native-spline-omp`** is the pure-MPI vs hybrid-MPI+OMP
  comparison with spline strategy held fixed.  The measured difference is the
  net of OMP synchronisation overhead (sections 4.2–4.6) against threaded
  speedup of the j-loops in `nnp_calc_acsf` and the `(idx, j)` loops in
  `nnp_scale_acsf`, plus the cache-miss savings from §4.7 (which appear in
  the OMP path more strongly because the threadprivate reduction pattern
  raises memory pressure on the shared arrays).  At small per-rank atom
  counts the synchronisation overhead is expected to dominate; the question
  is at what system / rank ratio hybrid wins.

Treating each branch as a separable benchmark makes it possible to point at
which optimisation strategy is responsible for which performance change,
and avoids the trap of "is this faster because of A, B, or both?" that a
single mega-branch would create.
