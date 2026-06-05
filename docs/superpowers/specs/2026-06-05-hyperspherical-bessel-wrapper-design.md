# Design: Expose hyperspherical Bessel functions in the classy wrapper

**Date:** 2026-06-05
**Resolves:** issue #147 ("Hyperspherical Bessel functions in wrapper"); enables closing stale PR #153.

## Background

The hyperspherical Bessel functions Φ_l^β(x) are computed in `tools/hyperspherical.cpp`
and used internally by the transfer module, but they have never been exposed through
the Python (`classy`) wrapper. PR #153 (open since Feb 2023) was the original vehicle
for this issue, but:

- Its only correctness content — the `hit_the_ceiling` / `PhiL_plus_one` bugfix in the
  chunked recurrence path of `hyperspherical_HIS_create` — is **already in master**
  (`tools/hyperspherical.cpp:153–221`).
- It never actually touched any wrapper file.
- It is based on 2023 master (`22367274`) and has diverged so far it cannot be merged
  without reverting years of refactoring.

So PR #153 will be **closed** once this work lands. Issue #147's actual request (wrapper
exposure) is implemented here. The PR's intended OpenMP cleanup, never carried to master,
is included here as a small separate change.

## Design principle: "a special function that just works"

The numerical subtleties of Φ_l^β(x) are already solved inside the C++ code, but the
logic is **scattered** across `HIS_create`, `ClosedModY`, `get_CF1`, and the recurrence
routines. The wrapper must not re-implement or expose these foot-guns. Instead it presents
three robust special functions, with every numerical consideration kept in C++:

- **Closed-case symmetry (K=1).** `ClosedModY(l, β, &y, &phisign, &dphisign)` folds y mod
  2π and reflects about π and π/2 into the fundamental domain `[0, π/2]`, flipping
  `phisign` (parity of l, and of β−l) and `dphisign`. The recurrence and interpolation
  grid only ever run on `[0, π/2]`; all other x are obtained by symmetry. Therefore the
  user may pass **any** x > 0 in the closed case and get the correct value.
- **Turning-point stability.** Below the classical turning point `xfwd`, forward
  recurrence is unstable; backward recurrence (seeded by the `get_CF1` continued fraction
  / `CF1_from_Gegenbauer`, with overflow rescaling by `_HYPER_OVERFLOW_`) is used instead.
  Above `xfwd`, forward recurrence is used. `xfwd` must be computed from the **maximum**
  requested l (the most restrictive turning point), since one recurrence call fills the
  whole l-ladder. `xfwd` = `sqrt(lmax(lmax+1))/β` (flat), `asin(...)` (closed),
  `asinh(...)` (open).
- **Maximum l in the closed case.** `sqrtK[l] = sqrt(β²−l²)` requires **l < β**; l=β−1 is
  the maximum and is the `hit_the_ceiling` special case. Requests with l ≥ β are rejected.
- **Small-x cutoff.** Below `chi_at_phimin` the function is ≈0; the interpolation path
  returns exactly 0 there, and `direct` returns the (tiny) recurrence value.
- **Curvature support of WKB.** `hyperspherical_WKB` handles **only K = ±1** (returns
  `_FAILURE_` for K=0) and needs **l ≥ 1** (`e = 1/sqrt(l(l+1))`). It folds closed-case
  symmetry internally.

## Goals

1. Expose Φ_l^β(x) from Python through three robust evaluation paths: interpolation,
   direct recurrence, WKB approximation — each handling the regimes above transparently.
2. Remove the dead OpenMP scaffolding from `hyperspherical_HIS_create`.
3. Add a pytest validating the three paths against each other and a known ground truth,
   including the symmetry / turning-point / curvature edge cases.
4. Close PR #153 and issue #147.

## Non-goals

- No new physics, no change to the transfer-module code path or its numerical output.
- Not exposing the l≤9 closed-form `HypersphericalExplicit` (recurrence supersedes it).
- Not parallelising the module with `std::thread` (out of scope; functions are fast
  enough for interactive use, and the transfer module parallelises at the q-loop level).

## Public API

Three **module-level** functions in `classy.pyx` (pure functions of (K, β, l, x) — no
cosmology object, no `.compute()`):

```python
from classy import (
    hyperspherical_bessel_interpolate,
    hyperspherical_bessel_direct,
    hyperspherical_bessel_wkb,
)
```

### Common arguments

| Arg     | Type                       | Constraints |
|---------|----------------------------|-------------|
| `K`     | int                        | Curvature sign ∈ {-1, 0, 1}. `wkb` requires K = ±1. |
| `beta`  | float                      | Wavenumber ν > 0. For K=1, must be a positive integer. |
| `l`     | int or array-like of int   | Multipole(s), ≥ 0. For K=1, every l < β. `wkb` requires l ≥ 1. |
| `x`     | array-like of float        | Evaluation points > 0. For K=1, any x > 0 is accepted (folded by symmetry). |

### Return value

A NumPy `float64` array of Φ:

- scalar `l` → shape `(n_x,)`
- array `l` → shape `(n_l, n_x)`, rows ordered as the input `l`.

### Per-function signatures

```python
hyperspherical_bessel_interpolate(K, beta, l, x, sampling=None, derivatives=False)
hyperspherical_bessel_direct(K, beta, l, x)
hyperspherical_bessel_wkb(K, beta, l, x)
```

- **`interpolate`** — `HIS_create` over the range of `x`, then
  `hyperspherical_Hermite_interpolation_vector` at the requested points. The production
  (bugfixed) path used by the transfer module; already handles symmetry (`ClosedModY`),
  turning point, overflow, max-l, and small-x→0.
  - `sampling`: points per wavelength of the internal grid. `None` → CLASS defaults
    (`8.0` flat; `7.0` for ν<1000, `3.0` for ν≥1000). Overridable for convergence studies.
  - `derivatives=True` → returns `(Φ, Φ', Φ'')` (same shapes as Φ).
- **`direct`** — exact recurrence evaluated at each requested x via a **new thin C++
  orchestration helper** (see Implementation). Reuses `ClosedModY`, the fwd/bwd recurrence,
  and the `sqrtK` setup so it is as robust as the grid path but free of interpolation
  error. Works for all valid l. Returns Φ only.
- **`wkb`** — `hyperspherical_WKB` per (l, x). K = ±1 and l ≥ 1 only. Returns Φ only.

`derivatives` is accepted only by `interpolate`.

## Implementation

### New C++ helper for `direct` (`tools/hyperspherical.{h,cpp}`)

Add one function that encapsulates the direct-recurrence orchestration at an arbitrary
vector of x — keeping all numerical considerations in C++:

```c
int hyperspherical_bessel_direct_vector(int K, double beta,
                                        int* lvec, int nl,
                                        double* xvec, int nx,
                                        double* Phi,        // nl*nx, row-major by l
                                        ErrorMsg error_message);
```

Behaviour:
1. Validate K ∈ {-1,0,1}, β > 0; for K=1 require integer β and max(lvec) < β. On violation
   set `error_message` and return `_FAILURE_`.
2. Build `sqrtK[0..lmax+2]`, `one_over_sqrtK[...]` from (K, β) exactly as `HIS_create`.
3. Compute `xfwd` from `lmax = max(lvec)`.
4. Allocate `PhiL[lmax+2]`. For each x (one recurrence call per x, even in the closed
   case):
   - **Fold once.** For K=1, fold `y = x` into `[0, π/2]` (mod 2π, reflect about π, then
     about π/2) and record the two booleans "reflected about π" / "reflected about π/2".
     The folded y is **l-independent**; only the sign is l-dependent. For K≠1, y = x.
   - Compute `sinK`,`cotK` of the folded y.
   - Choose `hyperspherical_backwards_recurrence` if folded-y < `xfwd` else
     `hyperspherical_forwards_recurrence`; fill `PhiL[0..lmax]` in this single call.
   - For each requested l, derive `phisign` from the two reflection booleans (flip on
     reflect-about-π if l odd; flip on reflect-about-π/2 if (β−l) even) and write
     `Phi[il*nx + ix] = phisign * PhiL[l]`.
5. Return `_SUCCESS_`.

This reproduces `ClosedModY`'s fold + sign rules (`tools/hyperspherical.cpp`, the
`ClosedModY` definition; used per-`lnum` by the interpolation routine at line 326) but,
because the fold itself does not depend on l, folds once per x and applies the cheap
per-l sign afterward — one recurrence call per x serves the whole l-ladder in every
curvature case.

### Cython surface (`classy.pyx`)

Add an inline `cdef extern from "hyperspherical.h":` block (matching the existing
`cdef extern from "cosmology.h":` at classy.pyx:103). Declare only what is needed:
`ErrorMsg`, opaque `HyperInterpStruct`, `hyperspherical_HIS_create`,
`hyperspherical_HIS_free`, `hyperspherical_Hermite_interpolation_vector`,
`hyperspherical_bessel_direct_vector` (new), and `hyperspherical_WKB`.

**`cclassy.pxd` is NOT touched** — it is auto-generated by `generate_wrapper.py` (which
scans only `class.h`/`common.h`/`cosmology.h`).

### Shared Cython helpers (module-level)

- `_normalize_bessel_inputs(K, beta, l, x)` → validated C-contiguous `float64` x array,
  `int32` l array, scalar-l flag, lmax. Python-side validation produces clear errors
  before entering C: K ∈ {-1,0,1}; β > 0; K=1 ⇒ integer β and max(l) < β; l ≥ 0; x > 0.
  Method-specific extra checks: `wkb` ⇒ K ≠ 0 and min(l) ≥ 1.
- Error handling: an `ErrorMsg` buffer; any C call returning `_FAILURE_` raises
  `CosmoSevereError(<buffer>)`.

### `interpolate` path

1. Normalize; `lvec = sorted(unique(l))` int32.
2. `xmin = max(min(x), 1e-5)`, `xmax = max(x)` — for K=1 clamp to `[ε, π/2−ε]` and rely on
   `ClosedModY` to map requested x in. `phiminabs = 1e-10`, `l_WKB = lmax+1`, sampling
   resolved.
3. `cdef HyperInterpStruct his`; `hyperspherical_HIS_create(...)`.
4. `try`: for each requested l (its index in `lvec`), call
   `hyperspherical_Hermite_interpolation_vector(&his, n_x, l_index, x_ptr, Φ_row,
   dΦ_row_or_NULL, d2Φ_row_or_NULL)`. `finally`: `hyperspherical_HIS_free(&his, ...)`.
5. Assemble; collapse to 1-D if scalar l.

### `direct` path

Normalize, then a single call to `hyperspherical_bessel_direct_vector(...)` with `lvec`
= the requested l (order preserved); assemble; collapse if scalar l.

### `wkb` path

Normalize (enforce K=±1, l≥1). For each (l, x) call `hyperspherical_WKB(K, l, beta, x,
&Φ)`; assemble; collapse if scalar l.

### OpenMP cleanup (separate commit)

In `hyperspherical_HIS_create` (~lines 134–227) remove the three inert `#pragma omp`
directives, the `abort` variable + its init + the `if (abort == _TRUE_) return _FAILURE_;`
block, and the now-redundant braces. Pure dead-code removal (build uses no `-fopenmp`); the
`hit_the_ceiling` bugfix and chunked loops are preserved verbatim.

## Testing — `python/test_hyperspherical.py` (pytest)

1. **K=0 ground truth.** Φ_l^β(x) = j_l(βx); compare `interpolate` and `direct` against
   `scipy.special.spherical_jn(l, β·x)` over a grid of l, x (≈1e-6 relative away from x→0).
2. **Cross-method.** `interpolate` vs `direct` over the full l range agree to interpolation
   tolerance away from `xmin`.
3. **Closed-case symmetry (K=1).** Evaluate `direct`/`interpolate` at x ∈ (π/2, 2π) and
   verify they equal the symmetry-folded values in `[0, π/2]` (same machinery both should
   already encode — this checks the wrapper exposes it correctly, including sign flips for
   odd/even l).
4. **Turning-point continuity.** Across `xfwd` (bwd→fwd switch) `direct` is continuous and
   matches `interpolate` — guards against a stability/normalization seam.
5. **Curved WKB.** For K=±1, large l, `wkb` ≈ `direct` within WKB tolerance in the
   oscillatory region.
6. **Shape/dtype contract.** scalar l → 1-D; array l → 2-D with correct row order;
   `derivatives=True` → 3-tuple of matching shapes.
7. **Error paths.** K=1 non-integer β raises; K=1 with l ≥ β raises; `wkb` with K=0 raises;
   `wkb` with l=0 raises; invalid K raises; `derivatives=True` on `direct`/`wkb` raises.

Also verify the OpenMP cleanup did not change transfer output: run an existing
transfer/Cl scenario before and after; expect bit-identical results.

## Closeout

After merge to master: close PR #153 (note bugfix already in master, exposure delivered
here) and close issue #147 as completed.

## Risks / open questions

- **Closed-case fold in `direct_vector`.** The fold of y into `[0, π/2]` is l-independent;
  only the sign flips depend on l. So the helper folds once per x, runs one recurrence call
  per x, and applies the per-l sign afterward — efficient and correct in all curvature
  cases. (The existing interpolation routine calls `ClosedModY` per `lnum` only because it
  processes one l at a time.)
- **Interpolation grid bounds.** Uniform grid on `[xmin, xmax]`; very small `min(x)` is
  clamped to `1e-5`; out-of-support x return ≈0. Acceptable for an inspection/validation
  API.
- **WKB accuracy** is intentionally approximate near turning points; test tolerances and
  the docstring reflect this.
