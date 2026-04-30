# NNP Implementation Changes: `feature/nnp-verlet-cells` vs `upstream/master`

This document records all key implementation differences between the
`feature/nnp-verlet-cells` branch and `upstream/master` across the four
modified NNP source files.  The intent is a reference for code review and
future maintenance — not an API guide.

---

## 1. `nnp_acsf.F`

This file received the largest changes.  The upstream file is ~1 414 lines;
the branch version is ~2 192 lines.

### 1.1 New module-level types

Two private derived types were added at module scope.

#### `nnp_cell_list_cache_type`
A Verlet-skin / cell-list hybrid cache.  Each `(atom j, integer PBC shift)`
pair becomes one *image entry* stored explicitly, so the inner neighbour-query
loop contains no `pbc()` calls.  The coarse grid has `bin_width ≈ list_cutoff`
(one bin-span walk covers 27 bins).

Key fields:

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
Stores Hermite cubic spline coefficients `(a, b, c, d)` for Horner's-method
evaluation of `EXP`, `COS`, and `TANH`.  Coefficients are packed as
`coef(4, n_points)` so each evaluation reads four consecutive doubles.  A
`SAVE` instance `spline_cache` is initialised once in `nnp_init_acsf_groups`.

---

### 1.2 New module-level SAVE variables

```fortran
TYPE(nnp_cell_list_cache_type), SAVE, PRIVATE :: cell_list_cache
TYPE(nnp_spline_cache_type),    SAVE, PRIVATE :: spline_cache

TYPE(nnp_neighbor_type), SAVE, PRIVATE, TARGET :: persistent_neighbor
INTEGER, SAVE, PRIVATE :: persistent_n_capacity = -1
INTEGER, SAVE, PRIVATE :: persistent_n_rad_grp  = -1
INTEGER, SAVE, PRIVATE :: persistent_n_ang_grp  = -1

INTEGER,       ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_max_sym
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_sym(:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_forcetmp(:,:)
REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_force3tmp(:,:,:)
```

`persistent_neighbor` replaces the per-centre-atom local `nnp_neighbor_type`.
In the upstream, each call to `nnp_calc_acsf` allocated and deallocated 9
arrays sized `~num_atoms × (2p+1)³`; for a 6 000-atom system those allocations
cross glibc's `mmap` threshold and become system calls, giving O(N²) `malloc`
traffic per MD step.

`scratch_sym/forcetmp/force3tmp` replace per-symfgrp-group
ALLOCATE/DEALLOCATE in the hot radial and angular loops.

---

### 1.3 `nnp_calc_acsf` — signature extended

```fortran
! upstream
SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz, stress)

! branch
SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz, stress, active_atoms, active_flag, n_active)
```

`active_atoms(:)`, `active_flag(:)`, and `n_active` are all `OPTIONAL`.  When
present, the routine builds an active-atoms set (centre atom `i` plus every
distinct neighbour touched by the radial and angular lists) as a deduplication
pass after `nnp_compute_neighbors`.  This is passed back to `nnp_force.F` so
the force-accumulation loop is restricted to O(n_active) atoms instead of O(N).

---

### 1.4 Scratch buffer management — `POINTER, CONTIGUOUS` local aliases

In the upstream, `symtmp`, `forcetmp`, and `force3tmp` were local
`ALLOCATABLE` arrays declared fresh in each call:

```fortran
! upstream
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)         :: symtmp
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:)       :: forcetmp
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: force3tmp
! ... per symfgrp:
ALLOCATE(symtmp(n_symf)); ALLOCATE(force3tmp(3,3,n_symf))
! ... use ...
DEALLOCATE(symtmp); DEALLOCATE(force3tmp)
```

In the branch they are `POINTER, CONTIGUOUS` aliases into the module-level
SAVE arrays:

```fortran
! branch
REAL(KIND=dp), DIMENSION(:),     POINTER, CONTIGUOUS :: symtmp
REAL(KIND=dp), DIMENSION(:,:),   POINTER, CONTIGUOUS :: forcetmp
REAL(KIND=dp), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: force3tmp
! ... per symfgrp:
symtmp    => scratch_sym(1:n_symf)
force3tmp => scratch_force3tmp(:,:,1:n_symf)
! ... use ...
NULLIFY(symtmp); NULLIFY(force3tmp)
```

**Why `CONTIGUOUS` matters:** without this attribute the compiler must treat
the pointer as potentially non-contiguous and potentially aliasing other
`TARGET` or `POINTER` variables in scope.  Both restrictions block
vectorisation of the inner `sf`/`l` loops in the angular double loop.
`CONTIGUOUS` tells the compiler the target is a unit-stride contiguous block,
restoring the same optimisation level as the original `ALLOCATABLE` version.

---

### 1.5 Neighbour building: O(N) cell-list replaces O(N·(2p+1)³) direct loop

#### Upstream — `nnp_neighbor_create` + `nnp_neighbor_release`

Each call to `nnp_calc_acsf` triggered:
1. `nnp_compute_pbc_copies` to find the periodic image count `p`
2. 9× ALLOCATE (arrays sized `num_atoms × (2p+1)³ × n_symfgrp`)
3. A triple PBC loop (`c1`, `c2`, `c3`) over all `(j, shift)` pairs, calling
   `pbc(dr, cell, nl)` for every displacement
4. `nnp_neighbor_release` — 9× DEALLOCATE

#### Branch — `nnp_compute_neighbors` (cell-list walk)

`nnp_prepare_cell_list_cache` is called **once per force evaluation** from
`nnp_calc_energy_force` before the per-centre loop.  `nnp_compute_neighbors`
is then called once per centre atom and walks only the 27 bins adjacent to the
centre's bin in the pre-built image table.  No allocations occur in the
per-atom hot path; `pbc()` is not called in the inner loop.

The inner loop filters candidates by:
1. Cheap `r² ≥ max_cut²` scalar reject (avoids `SQRT`)
2. Per-group `r² < cutoff_sq` for radial groups
3. Per-group `r² < cutoff_sq` for angular groups (ele1 and ele2 separately)

`SQRT` is computed **lazily** (only when at least one group accepts the
candidate), hoisted using a sentinel `norm = -1.0`.

---

### 1.6 `nnp_prepare_cell_list_cache` — Verlet skin logic

Called once per force evaluation.  The cheap path (O(N)) refreshes
`coord_primary` with PBC-wrapped coordinates and checks six rebuild triggers:

| Trigger | Why |
|---|---|
| `!initialized` | First call |
| `num_atoms` changed | System resized |
| `exact_cutoff` changed | Cutoff parameters changed |
| `verlet_skin` changed | Skin changed |
| `perd` or `hmat` changed | Cell geometry changed |
| Any atom moved > 0.5·skin | Verlet displacement criterion |

The displacement check uses raw `coord_primary − ref_coord_primary` without
`pbc()`.  An atom that crosses a box face in the wrapped coordinates will show
an apparent large displacement (~box size), causing a conservative rebuild.
This is safe: false-positive rebuilds are cheap; false-negative rebuilds would
produce wrong neighbour lists.

The expensive path (`nnp_build_cell_list_cache`, O(N·(2p+1)³)) is only entered
when a trigger fires.

Skin value: `skin = MIN(1.0_dp Å, 0.1·max_cut)`.

---

### 1.7 Bug fix — wrong `kk` index in hetero-angular path

In the upstream, the hetero-angular loop (element 1 ≠ element 2) contained
a dead write followed by a redundant inner-loop assignment:

```fortran
! upstream (buggy)
jj = neighbor%ind_ang1(j, s)
kk = neighbor%ind_ang1(k, s)   ! DEAD — immediately overwritten inside sf loop
CALL nnp_calc_ang(...)
DO sf = 1, n_symf
   m = off + ...
   jj = neighbor%ind_ang1(j, s)   ! redundant: same value as outer scope
   kk = neighbor%ind_ang2(k, s)   ! correct, but re-read every sf iteration
END DO
```

The outer `kk = ind_ang1(k,s)` was the wrong index (should be `ind_ang2`), but
it was immediately shadowed by the correct assignment inside the loop, so the
computation was numerically correct.  However, re-reading `ind_ang1` and
`ind_ang2` on every `sf` iteration adds unnecessary memory traffic.

Branch fix:

```fortran
! branch (fixed)
jj = neighbor%ind_ang1(j, s)
kk = neighbor%ind_ang2(k, s)   ! correct index, set once before CALL
CALL nnp_calc_ang(...)
DO sf = 1, n_symf
   m = off + ...
   ! jj and kk are already correct — no inner-loop re-reads
END DO
```

---

### 1.8 Spline replacement of `EXP`, `COS`, and `TANH`

In `nnp_calc_sym_rad` and `nnp_calc_ang`, the expensive intrinsic calls are
replaced by Hermite cubic spline lookups initialised in `nnp_init_acsf_groups`:

| Intrinsic | Spline range | Rationale |
|---|---|---|
| `EXP(x)` | `[−50, 0]` | Argument is always ≤ 0 (Gaussian decay) |
| `COS(π·r/rc)` | `[0, π]` | Cutoff argument lives in `[0, 1)` so `tmp ∈ [0, π]` only — doubles grid resolution vs `[−π, π]` |
| `TANH(1−r/rc)` | `[0, 1]` | Argument lives in `[0, 1)` — 20× finer resolution vs `[−10, 10]` |

Each spline uses 5 000 grid points.  Coefficients are stored as `coef(4, n)`
(column-major, Horner form) so each evaluation reads four contiguous doubles:
`a + t*(b + t*(c + t*d))`.

For the EXP term, the spline value `symtmp` (≈ EXP) is used directly in the
chain rule for the force (`dsymdr = symtmp * (−2η(r−rs))`), rather than
evaluating the spline derivative separately.  This is both faster and slightly
more accurate (O(h⁴) value vs O(h³) formal derivative).

---

### 1.9 Division elimination

Scalar reciprocal-multiply replaces division in three places:

```fortran
! upstream
drdx(:) = rvect(:)/r

! branch
rinv = 1.0_dp/r
drdx(1) = rvect(1)*rinv; drdx(2) = rvect(2)*rinv; drdx(3) = rvect(3)*rinv
```

The same pattern is applied to `r1`, `r2`, `r3` in the angular routine.

Division is typically 6× slower than multiplication on x86; three multiplies
share one division.

Cutoff argument divisions are moved to `nnp_init_acsf_groups` as precomputed
reciprocals stored in the symfgrp struct:

```fortran
! upstream
tmp = pi*r/funccut          ! division in hot loop
tmp = 1.0_dp - r/funccut   ! division in hot loop

! branch — uses precomputed fields set at init
tmp = r * symfgrp%pi_over_cut
tmp = 1.0_dp - r * symfgrp%inv_cut
```

`cutoff_sq = funccut²` is also precomputed to eliminate a `SQRT` in the
neighbour filter (branch uses `r² < cutoff_sq` instead of `norm < cutoff`).

---

### 1.10 `nnp_calc_ang` — hoisted loop invariants

Two sub-expressions are hoisted out of the per-symmetry-function loop:

```fortran
! branch — computed once per (j,k) pair, outside the sf loop
r_sum_sq = rsqr1 + rsqr2 + rsqr3   ! EXP argument; only eta varies per sf
g_sq_inv = 1.0_dp/(g*g)             ! denominator in dangular/dr
```

In the upstream these were recomputed (or reloaded) on every `sf` iteration.

---

### 1.11 Pre-computed `zeta_int` / `zeta_is_int`

The upstream checked integer-zeta at runtime on every angular SF evaluation:

```fortran
! upstream — inside sf loop
i = NINT(zeta)
IF (1.0_dp*i == zeta) THEN
   tmpzeta = tmp**(i - 1)   ! fast integer power
```

The branch precomputes these once at init and stores them in the `nnp_ang_type`
struct:

```fortran
! branch — computed in nnp_init_acsf_groups, used in hot loop
IF (nnp%ang(ind)%zeta_is_int(m)) THEN
   tmpzeta = tmp**(nnp%ang(ind)%zeta_int(m) - 1)
```

This removes one `NINT` call and one FP equality test from every angular SF
evaluation.

---

### 1.12 Element→group lookup tables (infrastructure, not yet used in hot path)

`nnp_init_acsf_groups` now builds reverse-index maps:

| Array | Meaning |
|---|---|
| `rad(i)%ele_grp_count(e)` | Number of radial symfgrps for centre element `i` that accept neighbour element `e` |
| `rad(i)%ele_grp_list(g, e)` | Group indices in that set |
| `ang(i)%ele1_grp_count(e)` / `ele2_grp_count(e)` | Same for the two angular partner roles |

These enable an O(1) dispatch "given neighbour element `e`, which groups does
it enter?" rather than the current O(n_symfgrp) scan in `nnp_compute_neighbors`.
The tables are allocated and populated but not yet wired into the neighbour
walk.

---

### 1.13 `nnp_scale_acsf` — active-atoms fast path

```fortran
! upstream — always O(num_atoms)
DO k = 1, nnp%num_atoms
   DO j = 1, nnp%n_rad(ind)
      dsymdxyz(:, j, k) = dsymdxyz(:, j, k) / (loc_max - loc_min) * (scmax - scmin)

! branch — O(n_active) when active_atoms present
IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
   DO idx = 1, n_active
      k = active_atoms(idx)
      DO j = 1, nnp%n_rad(ind) ...
```

The same pattern is applied to both the radial and angular halves of both the
`symtype_loc` and `symtype_std` scaling branches.

---

## 2. `nnp_force.F`

### 2.1 MPI parallelisation — allgather eliminated

#### Upstream
```fortran
mecalc = nnp%num_atoms/logger%para_env%num_pe + &
         MIN(MOD(nnp%num_atoms, logger%para_env%num_pe)/(logger%para_env%mepos + 1), 1)
ALLOCATE(allcalc(logger%para_env%num_pe))
allcalc(:) = 0
CALL logger%para_env%allgather(mecalc, allcalc)
istart = 1
DO i = 2, logger%para_env%mepos + 1
   istart = istart + allcalc(i - 1)
END DO
! ... end of routine:
DEALLOCATE(allcalc)
```

One MPI collective (`allgather`) per force evaluation, plus allocate/deallocate
of a `num_pe`-sized integer array.

#### Branch — analytical derivation
```fortran
n_base      = nnp%num_atoms/logger%para_env%num_pe
n_remainder = MOD(nnp%num_atoms, logger%para_env%num_pe)
IF (logger%para_env%mepos < n_remainder) THEN
   mecalc = n_base + 1
   istart = 1 + logger%para_env%mepos*(n_base + 1)
ELSE
   mecalc = n_base
   istart = 1 + n_remainder*(n_base + 1) + (logger%para_env%mepos - n_remainder)*n_base
END IF
```

The formula is mathematically identical to the upstream (same block
distribution), verified by expanding both expressions symbolically.  The
allgather is eliminated because `istart` is a deterministic closed-form
function of `(num_atoms, num_pe, mepos)`.

---

### 2.2 SAVE scratch arrays

```fortran
! upstream — re-allocated every force eval
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)         :: denergydsym
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: dsymdxyz, stress

! branch — SAVE; allocated once, reused across steps
INTEGER,       ALLOCATABLE, DIMENSION(:), SAVE :: active_atoms, active_flag
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), SAVE :: denergydsym
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: dsymdxyz, stress
```

---

### 2.3 Hoisted allocation guard

#### Upstream
`dsymdxyz` was allocated inside the per-centre loop with size fitted to the
current atom's element type:
```fortran
DO i = istart, istart + mecalc - 1
   ...
   IF (calc_forces) THEN
      ALLOCATE(dsymdxyz(3, nnp%arc(ind)%n_nodes(1), nnp%num_atoms))
      dsymdxyz(:,:,:) = 0.0_dp
      IF (calc_stress) THEN
         ALLOCATE(stress(3, 3, nnp%arc(ind)%n_nodes(1)))
         stress(:,:,:) = 0.0_dp
      END IF
```

This means: O(mecalc) allocations of an array of size `3 × n_nodes(1) ×
num_atoms` per force eval.

#### Branch
Allocated once at force-eval entry, sized to `max_nsym` (maximum `n_nodes(1)`
across all element types), guarded by a SAVE check:

```fortran
IF (calc_forces) THEN
   max_nsym = MAXVAL([nnp%arc(ie)%n_nodes(1), ie=1..n_ele])
   IF (.NOT. ALLOCATED(dsymdxyz) .OR. SIZE(dsymdxyz, 3) /= nnp%num_atoms) THEN
      IF (ALLOCATED(dsymdxyz)) DEALLOCATE(dsymdxyz, active_atoms, active_flag, denergydsym)
      IF (ALLOCATED(stress))   DEALLOCATE(stress)
      ALLOCATE(dsymdxyz(3, max_nsym, nnp%num_atoms), SOURCE=0.0_dp)
      ALLOCATE(active_atoms(nnp%num_atoms))
      ALLOCATE(active_flag(nnp%num_atoms), SOURCE=0)
      ALLOCATE(denergydsym(max_nsym))
   END IF
   IF (calc_stress .AND. .NOT. ALLOCATED(stress)) ALLOCATE(stress(3, 3, max_nsym))
END IF
```

`SOURCE=0` (Fortran 2003) handles the initial zero-fill at allocation time.

---

### 2.4 Selective zeroing — O(n_active) instead of O(N)

#### Upstream
After each centre atom, the entire `dsymdxyz` buffer was cleared:
```fortran
dsymdxyz(:,:,:) = 0.0_dp   ! O(3 × max_nsym × N) every centre atom
```

#### Branch
Only entries touched by the current centre's neighbour set are zeroed:
```fortran
DO idx = 1, n_active
   k = active_atoms(idx)
   dsymdxyz(:, 1:nnp%arc(ind)%n_nodes(1), k) = 0.0_dp
   active_flag(k) = 0
END DO
```

Cost drops from O(N) to O(n_active) ≈ O(k_neighbours) per centre atom.

---

### 2.5 Force accumulation — O(n_active) loop

#### Upstream — O(N) inner loop
```fortran
DO j = 1, nnp%arc(ind)%n_nodes(1)
   DO k = 1, nnp%num_atoms
      DO m = 1, 3
         nnp%myforce(m, k, i_com) = ... - denergydsym(j)*dsymdxyz(m, j, k)
      END DO
   END DO
   IF (calc_stress) ...
END DO
```

This iterates over all N atoms for each centre, most of which have
`dsymdxyz(:,:,k) = 0`.

#### Branch — O(n_active) inner loop, stress separated
```fortran
! Stress: separate j-loop, no per-atom inner loop
IF (calc_stress) THEN
   DO j = 1, nnp%arc(ind)%n_nodes(1)
      nnp%committee_stress(:,:,i_com) = nnp%committee_stress(:,:,i_com) &
                                        - denergydsym(j)*stress(:,:,j)
   END DO
END IF
! Force: walk only the active neighbour set
DO idx = 1, n_active
   k = active_atoms(idx)
   DO j = 1, nnp%arc(ind)%n_nodes(1)
      DO m = 1, 3
         nnp%myforce(m, k, i_com) = ... - denergydsym(j)*dsymdxyz(m, j, k)
      END DO
   END DO
END DO
```

Loop order is also reordered: `k` outermost, `j` middle, `m` innermost.  In
the upstream the j-outermost order accessed `dsymdxyz(m, j, k)` with varying
`k` in the inner loop — Fortran column-major layout means this scattered across
the third dimension.  The reordering keeps `k` fixed in the inner `j, m` loops,
accessing a contiguous `dsymdxyz(:, :, k)` block.

---

### 2.6 Cell-list integration

`nnp_prepare_cell_list_cache(nnp)` is called once per force evaluation before
the per-centre atom loop.  `nnp_calc_acsf` is called with the optional
`active_atoms`, `active_flag`, `n_active` arguments.

---

## 3. `nnp_environment_types.F`

### 3.1 New fields in `nnp_symfgrp_rad_type`

```fortran
REAL(KIND=dp) :: cutoff_sq   = -1.0_dp   ! funccut², avoids SQRT in neighbour filter
REAL(KIND=dp) :: inv_cut     = -1.0_dp   ! 1/funccut, avoids division in cutoff arg
REAL(KIND=dp) :: pi_over_cut = -1.0_dp   ! π/funccut, avoids division in cutoff arg
```

### 3.2 New fields in `nnp_rad_type`

```fortran
INTEGER, ALLOCATABLE, DIMENSION(:)    :: ele_grp_count   ! (n_ele)
INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: ele_grp_list    ! (n_symfgrp, n_ele)
```

Element→radial-group reverse index (see §1.12).

### 3.3 New fields in `nnp_symfgrp_ang_type`

Same three precomputed cutoff fields as the radial type.

### 3.4 New fields in `nnp_ang_type`

```fortran
INTEGER, ALLOCATABLE, DIMENSION(:) :: zeta_int(:)      ! NINT(zeta), precomputed
LOGICAL, ALLOCATABLE, DIMENSION(:) :: zeta_is_int(:)   ! .TRUE. when zeta is integer-valued

INTEGER, ALLOCATABLE, DIMENSION(:)   :: ele1_grp_count, ele2_grp_count
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ele1_grp_list,  ele2_grp_list
```

`zeta_int` / `zeta_is_int` eliminate the `NINT` + FP-equality check from every
angular SF evaluation (see §1.11).  The `ele1/ele2` arrays are the angular
equivalent of the radial element→group maps (see §1.12).

---

## 4. `nnp_environment.F`

### 4.1 Additional allocations at init

```fortran
ALLOCATE(nnp_env%ang(i)%zeta_int(nnp_env%n_ang(i)))
ALLOCATE(nnp_env%ang(i)%zeta_is_int(nnp_env%n_ang(i)))
```

These are populated in `nnp_init_acsf_groups` (in `nnp_acsf.F`) after the
zeta arrays are finalised by the sort.

### 4.2 Deallocate additions

Companion deallocations for the new arrays are wired into the env-release path
via the existing `nnp_env_release` cleanup in `nnp_environment_types.F`.

---

## Summary table

| Change | File | Upstream | Branch | Why |
|---|---|---|---|---|
| Cell-list / Verlet cache | `nnp_acsf.F` | None — O(N·(2p+1)³) per centre | `nnp_prepare_cell_list_cache` + `nnp_build_cell_list_cache` | Amortises rebuild; O(k) per centre atom |
| Persistent neighbour buffer | `nnp_acsf.F` | 9× alloc/free per centre atom | Single `SAVE` buffer, reused | Eliminates O(N²) `malloc` traffic for large N |
| Scratch sym/force buffers | `nnp_acsf.F` | ALLOCATE/DEALLOCATE per symfgrp | Module SAVE arrays + `POINTER, CONTIGUOUS` alias | Eliminates per-group alloc; CONTIGUOUS restores vectorisation |
| `kk` index bug | `nnp_acsf.F` | Dead outer write + redundant inner reads | Correct single write before CALL | Correctness fix; also removes memory traffic |
| `EXP` / `COS` / `TANH` splines | `nnp_acsf.F` | Direct intrinsic calls | Hermite cubic spline lookup | Replaces expensive libm calls in the hot SF loop |
| Division → reciprocal-multiply | `nnp_acsf.F` | `rvect/r`, `pi*r/funccut`, etc. | Precomputed reciprocals | `×` is ~6× faster than `/` on x86 |
| `zeta_int` / `zeta_is_int` | `nnp_acsf.F` + types | NINT + FP-equality per SF | Precomputed boolean + int | Removes runtime NINT from angular hot loop |
| `r_sum_sq`, `g_sq_inv` hoists | `nnp_acsf.F` | Recomputed per sf | Computed once per (j,k) pair | Loop-invariant code motion |
| MPI allgather removal | `nnp_force.F` | `allgather(mecalc, allcalc)` + prefix sum | Closed-form analytical formula | Eliminates 1 collective + alloc/free per force eval |
| SAVE scratch arrays | `nnp_force.F` | Re-alloc every force eval | SAVE ALLOCATABLE | Eliminates O(step) allocation traffic |
| Hoisted allocation guard | `nnp_force.F` | Alloc per centre inside loop | Alloc once at force-eval entry, guarded | Reduces alloc calls from O(N) to O(1) per step |
| Selective zeroing | `nnp_force.F` | `dsymdxyz(:,:,:)=0` — O(N) | Zero only active-atoms slice — O(n_active) | Avoids full N-atom clear every centre |
| O(n_active) force loop | `nnp_force.F` | Loop over all N atoms | Loop over n_active touched atoms | Skips all-zero rows; O(k) at fixed density |
| Stress loop placement | `nnp_force.F` | Inside per-atom inner loop | Separate j-loop before per-atom loop | Removes branch/per-atom overhead |
| Element→group maps | `nnp_acsf.F` + types | O(n_symfgrp) scan per neighbour | Reverse-index tables (infrastructure) | Future O(1) group dispatch in neighbour walk |
| Precomputed `cutoff_sq`, `inv_cut`, `pi_over_cut` | types | Computed at every call | Stored in symfgrp struct at init | Eliminates repeated multiply/divide at call time |
