# Hyperspherical Bessel Wrapper Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose the hyperspherical Bessel functions Φ_l^β(x) through the `classy` Python wrapper via three robust special functions (interpolation, direct recurrence, WKB), and remove dead OpenMP scaffolding from `hyperspherical_HIS_create`.

**Architecture:** Three module-level functions in `classy.pyx`. The interpolation and WKB paths call existing C++ routines directly; the direct-recurrence path calls a new thin C++ helper that keeps all numerical considerations (closed-case symmetry, turning-point stability, max-l) in C++. The Cython layer only validates inputs and marshals NumPy buffers.

**Tech Stack:** C++ (`tools/hyperspherical.cpp`), Cython (`classy.pyx`), NumPy, pytest, scipy (test ground truth).

**Spec:** `docs/superpowers/specs/2026-06-05-hyperspherical-bessel-wrapper-design.md`

**Build/test note for every task:** After editing `classy.pyx` or any `.cpp`/`.h`, rebuild with `pip install .` from the repo root (this regenerates `cclassy.pxd`, cythonizes `classy.pyx`, and recompiles — takes a few minutes). Run tests from the `python/` directory: `cd python && python -m pytest test_hyperspherical.py -v`. Never hand-edit `cclassy.pxd` (auto-generated).

---

## Task 1: Remove dead OpenMP scaffolding from `hyperspherical_HIS_create`

Pure dead-code removal. The build passes no `-fopenmp`, so the `#pragma omp` lines are inert and the `abort` machinery is unreachable. The `hit_the_ceiling` bugfix and chunked recurrence loops must be preserved verbatim. Must be bit-identical to before.

**Files:**
- Modify: `tools/hyperspherical.cpp` (function `hyperspherical_HIS_create`, ~lines 30–227)

- [ ] **Step 1: Capture a baseline before any change**

Build current master and record Cls for a curved model that exercises the hyperspherical path.

Run from repo root:
```bash
pip install . >/tmp/build_baseline.log 2>&1 && python - <<'PY'
import numpy as np, classy
def cls(omega_k):
    c = classy.Class()
    c.set({'output':'tCl,pCl','l_max_scalars':800,'Omega_k':omega_k})
    c.compute()
    cl = c.raw_cl(800)
    c.struct_cleanup(); c.empty()
    return cl['tt'].copy()
base = {k: cls(k) for k in (-0.05, 0.05)}
np.savez('/tmp/hyper_baseline.npz', **{str(k):v for k,v in base.items()})
print('baseline saved', {k: v.shape for k,v in base.items()})
PY
```
Expected: prints `baseline saved {...}` with two arrays.

- [ ] **Step 2: Read the current function body**

Run: `sed -n '30,228p' tools/hyperspherical.cpp` (read-only orientation; confirm the lines below match before editing).

- [ ] **Step 3: Remove the `abort` declaration**

In `hyperspherical_HIS_create`, delete the line:
```cpp
  int abort;
```
(It sits among the locals near the top: `int current_chunk, index_x;` stays.)

- [ ] **Step 4: Remove the `#pragma omp parallel` block opener and its scope brace**

Replace this block:
```cpp
  abort = _FALSE_;

#pragma omp parallel shared(nx,                                                         \
                                pHIS,                                                   \
                                xfwd,                                                   \
                                K,                                                      \
                                l_recurrence_max,                                       \
                                beta,                                                   \
                                sqrtK,                                                  \
                                one_over_sqrtK,                                         \
                                lvec,                                                   \
                                nl,                                                     \
                                xfwdidx,                                                \
                                abort,                                                  \
                                error_message) private(j, k, l, current_chunk, index_x) \
    firstprivate(lmax)
  {
    std::vector<double> PhiL((lmax + 2) * _HYPER_CHUNK_);
```
with:
```cpp
  std::vector<double> PhiL((lmax + 2) * _HYPER_CHUNK_);
```

- [ ] **Step 5: Remove the first `#pragma omp for` directive**

Delete the line (just before the first `for (j = 0; j < MIN(nx, xfwdidx); ...`):
```cpp
#pragma omp for schedule(dynamic)
```
(There is a blank line after it; leaving one blank line is fine.)

- [ ] **Step 6: Remove the second `#pragma omp for` directive**

Delete the second occurrence (just before `for (j = xfwdidx; j < nx; ...`):
```cpp
#pragma omp for schedule(dynamic)
```

- [ ] **Step 7: Remove the parallel-region closing brace and the `abort` check**

The backwards/forwards loops were wrapped in the `{ ... }` opened in Step 4. After the forwards loop's closing `}`, replace:
```cpp
    }
  }
  if (abort == _TRUE_)
    return _FAILURE_;

  for (k = 0; k < nl; k++) {
```
with:
```cpp
    }
  }

  for (k = 0; k < nl; k++) {
```
(Note: the forwards loop body keeps its own closing `}` and the outer `for` loop keeps its `}`. Only the extra brace that closed the `#pragma omp parallel` region is removed. After editing, verify brace balance by compiling in Step 8.)

- [ ] **Step 8: Rebuild and verify bit-identical output**

Run from repo root:
```bash
pip install . >/tmp/build_after.log 2>&1 && python - <<'PY'
import numpy as np, classy
base = np.load('/tmp/hyper_baseline.npz')
def cls(omega_k):
    c = classy.Class()
    c.set({'output':'tCl,pCl','l_max_scalars':800,'Omega_k':omega_k})
    c.compute()
    cl = c.raw_cl(800)
    c.struct_cleanup(); c.empty()
    return cl['tt'].copy()
for k in (-0.05, 0.05):
    after = cls(k)
    assert np.array_equal(base[str(k)], after), f"MISMATCH at Omega_k={k}"
print("bit-identical OK")
PY
```
Expected: prints `bit-identical OK`. If it fails to compile, re-check brace balance from Step 7.

- [ ] **Step 9: Commit**

```bash
git add tools/hyperspherical.cpp
git commit -m "Remove dead OpenMP scaffolding from hyperspherical_HIS_create

Inert without -fopenmp; pure dead-code removal. Bugfix and chunked
recurrence preserved; transfer output verified bit-identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Direct-recurrence path (C++ helper + Cython function)

Add the C++ orchestration helper and expose it as `classy.hyperspherical_bessel_direct`. TDD against `scipy.special.spherical_jn` for the flat case (Φ_l^β(x) = j_l(βx)).

**Files:**
- Modify: `tools/hyperspherical.h` (add declaration inside the `extern "C"` block)
- Modify: `tools/hyperspherical.cpp` (add helper definition)
- Modify: `classy.pyx` (extern block, `_normalize_bessel_inputs`, `hyperspherical_bessel_direct`)
- Create: `python/test_hyperspherical.py`

- [ ] **Step 1: Write the failing test**

Create `python/test_hyperspherical.py`:
```python
import numpy as np
import pytest
from scipy.special import spherical_jn
from classy import hyperspherical_bessel_direct


def test_direct_flat_matches_spherical_jn():
    """Flat case: Phi_l^beta(x) = j_l(beta * x)."""
    beta = 30.0
    x = np.linspace(0.05, 2.0, 50)
    for l in (0, 1, 2, 5, 10, 20):
        got = hyperspherical_bessel_direct(0, beta, l, x)
        ref = spherical_jn(l, beta * x)
        assert got.shape == (x.size,)
        assert np.allclose(got, ref, rtol=1e-6, atol=1e-9), f"l={l}"


def test_direct_array_l_shape():
    beta = 30.0
    x = np.linspace(0.05, 2.0, 40)
    l = np.array([0, 3, 7])
    got = hyperspherical_bessel_direct(0, beta, l, x)
    assert got.shape == (3, 40)
    for i, ll in enumerate(l):
        assert np.allclose(got[i], spherical_jn(ll, beta * x), rtol=1e-6, atol=1e-9)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd python && python -m pytest test_hyperspherical.py -v`
Expected: FAIL — `ImportError: cannot import name 'hyperspherical_bessel_direct'`.

- [ ] **Step 3: Declare the helper in the C++ header**

In `tools/hyperspherical.h`, inside the `extern "C" { ... }` block (alongside the other declarations, e.g. just after the `hyperspherical_HIS_free` declaration), add:
```c
int hyperspherical_bessel_direct_vector(int K,
                                        double beta,
                                        int* lvec,
                                        int nl,
                                        double* xvec,
                                        int nx,
                                        double* Phi,
                                        ErrorMsg error_message);
```

- [ ] **Step 4: Implement the helper in the C++ source**

In `tools/hyperspherical.cpp`, add this function (place it after `hyperspherical_HIS_free`, before `hyperspherical_Hermite_interpolation_vector`):
```cpp
int hyperspherical_bessel_direct_vector(int K,
                                        double beta,
                                        int* lvec,
                                        int nl,
                                        double* xvec,
                                        int nx,
                                        double* Phi,
                                        ErrorMsg error_message) {
  /** Evaluate Phi_l^beta(x) by direct forwards/backwards recurrence at each
      requested x, for every requested l. Keeps the numerical considerations
      that are otherwise spread across hyperspherical_HIS_create:
        - closed-case (K=1) symmetry fold into [0, pi/2] via ClosedModY,
        - turning-point stability (backwards below xfwd, forwards above),
        - maximum l in the closed case (l < beta).
      Phi is row-major: Phi[il*nx + ix]. */
  int il, ix, l, lmax, intbeta = 0;
  double beta2 = beta * beta, xfwd, folded_y, sinK, cotK;
  int phisign, dphisign;

  class_test((K != 0) && (K != 1) && (K != -1),
             error_message, "K must be -1, 0, or 1 (got %d)", K);
  class_test(beta <= 0.0, error_message, "beta must be positive (got %g)", beta);

  lmax = 0;
  for (il = 0; il < nl; il++) {
    class_test(lvec[il] < 0, error_message, "l must be non-negative (got %d)", lvec[il]);
    if (lvec[il] > lmax) lmax = lvec[il];
  }

  if (K == 1) {
    intbeta = (int) (beta + 0.2);
    class_test(fabs(beta - intbeta) > 1e-6,
               error_message, "closed case (K=1) requires integer beta (got %g)", beta);
    class_test(lmax >= intbeta,
               error_message, "closed case requires l < beta; got max l=%d, beta=%d", lmax, intbeta);
  }

  std::vector<double> sqrtK(lmax + 3), one_over_sqrtK(lmax + 3);
  for (l = 0; l <= lmax + 2; l++) {
    if (K == 0)      sqrtK[l] = beta;
    else if (K == 1) sqrtK[l] = sqrt(beta2 - (double) l * l);
    else             sqrtK[l] = sqrt(beta2 + (double) l * l);
    one_over_sqrtK[l] = 1.0 / sqrtK[l];
  }

  if (K == 0)      xfwd = sqrt(lmax * (lmax + 1.0)) / beta;
  else if (K == 1) xfwd = asin(sqrt(lmax * (lmax + 1.0)) / beta);
  else             xfwd = asinh(sqrt(lmax * (lmax + 1.0)) / beta);

  std::vector<double> PhiL(lmax + 2);

  for (ix = 0; ix < nx; ix++) {
    folded_y = xvec[ix];
    if (K == 1) {
      /* Fold y into [0, pi/2]; the fold is l-independent (signs handled below). */
      ClosedModY(lvec[0], intbeta, &folded_y, &phisign, &dphisign);
    }

    if (K == 0)      { sinK = folded_y;       cotK = 1.0 / folded_y; }
    else if (K == 1) { sinK = sin(folded_y);  cotK = 1.0 / tan(folded_y); }
    else             { sinK = sinh(folded_y); cotK = 1.0 / tanh(folded_y); }

    if (folded_y < xfwd)
      hyperspherical_backwards_recurrence(K, lmax, beta, folded_y, sinK, cotK,
                                          sqrtK.data(), one_over_sqrtK.data(), PhiL.data());
    else
      hyperspherical_forwards_recurrence(K, lmax, beta, folded_y, sinK, cotK,
                                         sqrtK.data(), one_over_sqrtK.data(), PhiL.data());

    for (il = 0; il < nl; il++) {
      l = lvec[il];
      double sign = 1.0;
      if (K == 1) {
        double tmp = xvec[ix];
        phisign = 1; dphisign = 1;
        ClosedModY(l, intbeta, &tmp, &phisign, &dphisign);
        sign = phisign;
      }
      Phi[il * nx + ix] = sign * PhiL[l];
    }
  }
  return _SUCCESS_;
}
```

- [ ] **Step 5: Add the Cython extern block in `classy.pyx`**

After the existing `cdef extern from "cosmology.h":` block (it ends before `cdef class Class(PyCosmology):` near line 127), add:
```cython
cdef extern from "hyperspherical.h":
    cdef cppclass HyperInterpStruct:
        pass
    int hyperspherical_HIS_create(int K, double beta, int nl, int* lvec,
                                  double xmin, double xmax, double sampling,
                                  int l_WKB, double phiminabs,
                                  HyperInterpStruct* pHIS, char* error_message)
    int hyperspherical_HIS_free(HyperInterpStruct* pHIS, char* error_message)
    int hyperspherical_Hermite_interpolation_vector(HyperInterpStruct* pHIS, int nxi,
                                                    int lnum, double* xinterp,
                                                    double* Phi, double* dPhi, double* d2Phi)
    int hyperspherical_bessel_direct_vector(int K, double beta, int* lvec, int nl,
                                            double* xvec, int nx, double* Phi,
                                            char* error_message)
    int hyperspherical_WKB(int K, int l, double beta, double y, double* Phi)
```

- [ ] **Step 6: Add the input-normalization helper in `classy.pyx`**

At module level (e.g. just after the `cdef extern from "hyperspherical.h":` block), add:
```cython
def _normalize_bessel_inputs(K, beta, l, x, method):
    """Validate and coerce (K, beta, l, x). Returns (lvec_int32, xvec_float64,
    scalar_l_bool). Raises ValueError on invalid input. method in
    {'interpolate','direct','wkb'} for method-specific checks."""
    if K not in (-1, 0, 1):
        raise ValueError("K must be -1, 0, or 1 (got %r)" % (K,))
    beta = float(beta)
    if beta <= 0.0:
        raise ValueError("beta must be positive (got %r)" % (beta,))
    scalar_l = np.isscalar(l)
    lvec = np.ascontiguousarray(np.atleast_1d(l), dtype=np.intc)
    if np.any(lvec < 0):
        raise ValueError("all l must be non-negative")
    xvec = np.ascontiguousarray(np.atleast_1d(x), dtype=np.double)
    if np.any(xvec <= 0.0):
        raise ValueError("all x must be positive")
    if K == 1:
        if abs(beta - round(beta)) > 1e-6:
            raise ValueError("closed case (K=1) requires integer beta (got %r)" % (beta,))
        if int(lvec.max()) >= int(round(beta)):
            raise ValueError("closed case requires l < beta")
    if method == 'wkb':
        if K == 0:
            raise ValueError("wkb is only defined for curved space (K = +/-1)")
        if int(lvec.min()) < 1:
            raise ValueError("wkb requires l >= 1")
    return lvec, xvec, bool(scalar_l)
```

- [ ] **Step 7: Add the `hyperspherical_bessel_direct` function in `classy.pyx`**

At module level, add:
```cython
def hyperspherical_bessel_direct(K, beta, l, x):
    """Hyperspherical Bessel function Phi_l^beta(x) by direct recurrence.

    Parameters
    ----------
    K : int          curvature sign, -1 (open), 0 (flat), or 1 (closed)
    beta : float     wavenumber nu (> 0; integer for K=1)
    l : int or array of int   multipole(s) (>= 0; < beta for K=1)
    x : array of float        evaluation points (> 0)

    Returns a float64 array: shape (n_x,) for scalar l, else (n_l, n_x).
    Closed-case symmetry and turning-point stability are handled internally.
    """
    cdef int Kc = K
    cdef double betac = float(beta)
    cdef int[::1] lvec
    cdef double[::1] xvec
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'direct')
    lvec = lvec_arr
    xvec = xvec_arr
    cdef int nl = lvec_arr.shape[0]
    cdef int nx = xvec_arr.shape[0]
    cdef np.ndarray[double, ndim=2] Phi = np.empty((nl, nx), dtype=np.double)
    cdef char errmsg[2048]
    if hyperspherical_bessel_direct_vector(Kc, betac, &lvec[0], nl,
                                           &xvec[0], nx, &Phi[0, 0], errmsg) != 0:
        raise CosmoSevereError(errmsg)
    if scalar_l:
        return np.ascontiguousarray(Phi[0])
    return Phi
```

- [ ] **Step 8: Rebuild and run the test**

Run: `pip install . >/tmp/build_t2.log 2>&1 && cd python && python -m pytest test_hyperspherical.py -v`
Expected: PASS for `test_direct_flat_matches_spherical_jn` and `test_direct_array_l_shape`.

- [ ] **Step 9: Commit**

```bash
git add tools/hyperspherical.h tools/hyperspherical.cpp classy.pyx python/test_hyperspherical.py
git commit -m "Expose hyperspherical_bessel_direct in classy

New C++ direct_vector helper (reuses ClosedModY + recurrence; handles
closed-case symmetry, turning point, max-l) plus thin Cython wrapper.
Validated against scipy.special.spherical_jn for the flat case.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Interpolation path

Expose `classy.hyperspherical_bessel_interpolate` using the production `HIS_create` +
Hermite interpolation pipeline, with optional derivatives.

**Files:**
- Modify: `classy.pyx` (add `hyperspherical_bessel_interpolate`)
- Modify: `python/test_hyperspherical.py` (add tests)

- [ ] **Step 1: Write the failing tests**

Append to `python/test_hyperspherical.py`:
```python
from classy import hyperspherical_bessel_interpolate


def test_interpolate_flat_matches_spherical_jn():
    beta = 50.0
    x = np.linspace(0.1, 2.0, 60)
    for l in (0, 2, 8, 25):
        got = hyperspherical_bessel_interpolate(0, beta, l, x)
        ref = spherical_jn(l, beta * x)
        assert np.allclose(got, ref, rtol=1e-4, atol=1e-7), f"l={l}"


def test_interpolate_matches_direct_flat():
    beta = 50.0
    x = np.linspace(0.1, 2.0, 60)
    l = np.array([1, 4, 12])
    interp = hyperspherical_bessel_interpolate(0, beta, l, x)
    direct = hyperspherical_bessel_direct(0, beta, l, x)
    assert np.allclose(interp, direct, rtol=1e-4, atol=1e-7)


def test_interpolate_derivatives_shapes():
    beta = 50.0
    x = np.linspace(0.1, 2.0, 30)
    l = np.array([2, 5])
    out = hyperspherical_bessel_interpolate(0, beta, l, x, derivatives=True)
    assert isinstance(out, tuple) and len(out) == 3
    for arr in out:
        assert arr.shape == (2, 30)
    # First derivative of j_l(beta x) is beta * j_l'(beta x).
    Phi, dPhi, d2Phi = out
    ref0 = beta * spherical_jn(2, beta * x, derivative=True)
    assert np.allclose(dPhi[0], ref0, rtol=1e-3, atol=1e-6)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd python && python -m pytest test_hyperspherical.py -k interpolate -v`
Expected: FAIL — `ImportError: cannot import name 'hyperspherical_bessel_interpolate'`.

- [ ] **Step 3: Implement `hyperspherical_bessel_interpolate` in `classy.pyx`**

At module level, add:
```cython
def hyperspherical_bessel_interpolate(K, beta, l, x, sampling=None, derivatives=False):
    """Hyperspherical Bessel function Phi_l^beta(x) by grid + Hermite interpolation
    (the path used internally by the transfer module).

    Parameters as in hyperspherical_bessel_direct, plus:
    sampling : float or None    grid points per wavelength (None -> CLASS default:
                                8.0 flat; 7.0 for nu<1000, 3.0 otherwise).
    derivatives : bool          if True return (Phi, dPhi, d2Phi).

    Returns a float64 array (or a 3-tuple of them if derivatives=True): shape
    (n_x,) for scalar l, else (n_l, n_x).
    """
    cdef int Kc = K
    cdef double betac = float(beta)
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'interpolate')
    cdef int nx = xvec_arr.shape[0]

    # Sorted unique l for the interpolation structure; remember each request's index.
    uniq = np.unique(lvec_arr)
    cdef int[::1] luniq = np.ascontiguousarray(uniq, dtype=np.intc)
    cdef int nl = luniq.shape[0]
    cdef double[::1] xv = xvec_arr

    if sampling is None:
        samp = 8.0 if K == 0 else (7.0 if betac < 1000.0 else 3.0)
    else:
        samp = float(sampling)

    cdef double pi = np.pi
    cdef double xmin, xmax
    if K == 1:
        # Closed case: requested x are folded by ClosedModY into the fundamental
        # domain, so the grid must always span all of [eps, pi/2 - eps]
        # regardless of the requested x range.
        xmin = 1.0e-5
        xmax = 0.5 * pi - 1.0e-5
    else:
        xmin = max(float(xvec_arr.min()), 1.0e-5)
        xmax = float(xvec_arr.max())
    cdef int l_WKB = int(luniq[nl - 1]) + 1
    cdef double phiminabs = 1.0e-10

    cdef HyperInterpStruct his
    cdef char errmsg[2048]
    if hyperspherical_HIS_create(Kc, betac, nl, &luniq[0], xmin, xmax, samp,
                                 l_WKB, phiminabs, &his, errmsg) != 0:
        raise CosmoSevereError(errmsg)

    cdef np.ndarray[double, ndim=2] Phi = np.empty((nl, nx), dtype=np.double)
    cdef np.ndarray[double, ndim=2] dPhi
    cdef np.ndarray[double, ndim=2] d2Phi
    cdef int i
    try:
        if derivatives:
            dPhi = np.empty((nl, nx), dtype=np.double)
            d2Phi = np.empty((nl, nx), dtype=np.double)
            for i in range(nl):
                hyperspherical_Hermite_interpolation_vector(&his, nx, i, &xv[0],
                                                            &Phi[i, 0], &dPhi[i, 0], &d2Phi[i, 0])
        else:
            for i in range(nl):
                hyperspherical_Hermite_interpolation_vector(&his, nx, i, &xv[0],
                                                            &Phi[i, 0], NULL, NULL)
    finally:
        hyperspherical_HIS_free(&his, errmsg)

    # Map sorted-unique rows back to the requested l order.
    idx = np.searchsorted(uniq, lvec_arr)

    def _shape(M):
        out = M[idx]
        return np.ascontiguousarray(out[0]) if scalar_l else out

    if derivatives:
        return _shape(Phi), _shape(dPhi), _shape(d2Phi)
    return _shape(Phi)
```

- [ ] **Step 4: Rebuild and run the tests**

Run: `pip install . >/tmp/build_t3.log 2>&1 && cd python && python -m pytest test_hyperspherical.py -k interpolate -v`
Expected: PASS for the three interpolate tests.

- [ ] **Step 5: Commit**

```bash
git add classy.pyx python/test_hyperspherical.py
git commit -m "Expose hyperspherical_bessel_interpolate in classy

Wraps HIS_create + Hermite interpolation (production path), with optional
derivatives. Validated against scipy and the direct recurrence.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: WKB path

Expose `classy.hyperspherical_bessel_wkb` (curved-only, l >= 1).

**Files:**
- Modify: `classy.pyx` (add `hyperspherical_bessel_wkb`)
- Modify: `python/test_hyperspherical.py` (add tests)

- [ ] **Step 1: Write the failing tests**

Append to `python/test_hyperspherical.py`:
```python
from classy import hyperspherical_bessel_wkb


def test_wkb_matches_direct_open_large_l():
    """Open case: WKB approximates the recurrence well at large l, away from
    the turning point."""
    K, beta = -1, 200.0
    l = 80
    # Turning point asinh(sqrt(l(l+1))/beta); sample in the oscillatory region.
    x = np.linspace(1.0, 3.0, 40)
    wkb = hyperspherical_bessel_wkb(K, beta, l, x)
    direct = hyperspherical_bessel_direct(K, beta, l, x)
    # WKB is asymptotic; compare on a robust normalized error.
    scale = np.maximum(np.abs(direct), 1e-3 * np.abs(direct).max())
    assert np.median(np.abs(wkb - direct) / scale) < 0.05


def test_wkb_rejects_flat():
    with pytest.raises(ValueError):
        hyperspherical_bessel_wkb(0, 100.0, 5, np.array([0.5, 1.0]))


def test_wkb_rejects_l_zero():
    with pytest.raises(ValueError):
        hyperspherical_bessel_wkb(-1, 100.0, 0, np.array([0.5, 1.0]))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd python && python -m pytest test_hyperspherical.py -k wkb -v`
Expected: FAIL — `ImportError: cannot import name 'hyperspherical_bessel_wkb'`.

- [ ] **Step 3: Implement `hyperspherical_bessel_wkb` in `classy.pyx`**

At module level, add:
```cython
def hyperspherical_bessel_wkb(K, beta, l, x):
    """Hyperspherical Bessel function Phi_l^beta(x) by WKB (Airy) approximation.

    Curved space only (K = +/-1) and l >= 1. Parameters as in
    hyperspherical_bessel_direct. Returns a float64 array: shape (n_x,) for
    scalar l, else (n_l, n_x). Closed-case symmetry is handled internally.
    """
    cdef int Kc = K
    cdef double betac = float(beta)
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'wkb')
    cdef int[::1] lvec = lvec_arr
    cdef double[::1] xvec = xvec_arr
    cdef int nl = lvec_arr.shape[0]
    cdef int nx = xvec_arr.shape[0]
    cdef np.ndarray[double, ndim=2] Phi = np.empty((nl, nx), dtype=np.double)
    cdef int il, ix
    cdef double val
    cdef char errmsg[2048]
    for il in range(nl):
        for ix in range(nx):
            if hyperspherical_WKB(Kc, lvec[il], betac, xvec[ix], &val) != 0:
                raise CosmoSevereError(b"hyperspherical_WKB failed")
            Phi[il, ix] = val
    if scalar_l:
        return np.ascontiguousarray(Phi[0])
    return Phi
```

- [ ] **Step 4: Rebuild and run the tests**

Run: `pip install . >/tmp/build_t4.log 2>&1 && cd python && python -m pytest test_hyperspherical.py -k wkb -v`
Expected: PASS for the three wkb tests.

- [ ] **Step 5: Commit**

```bash
git add classy.pyx python/test_hyperspherical.py
git commit -m "Expose hyperspherical_bessel_wkb in classy

WKB (Airy) approximation, curved-only (K=+/-1, l>=1). Validated against
the direct recurrence in the open case.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Edge-case, symmetry, and contract tests

Lock in the regime handling that makes these "just work": closed-case symmetry,
turning-point continuity, the shape/dtype contract, and the remaining error paths.

**Files:**
- Modify: `python/test_hyperspherical.py` (add tests)

- [ ] **Step 1: Write the tests**

Append to `python/test_hyperspherical.py`:
```python
def test_closed_case_symmetry():
    """Closed case: direct evaluation beyond pi/2 equals the symmetry-folded
    value. Phi_l(pi - y) = (-1)^(beta-l-1) Phi_l(y)."""
    K = 1
    beta = 40  # integer
    l = 5
    y = np.linspace(0.05, 0.5 * np.pi - 0.05, 30)
    base = hyperspherical_bessel_direct(K, beta, l, y)
    reflected = hyperspherical_bessel_direct(K, beta, l, np.pi - y)
    sign = (-1.0) ** (beta - l - 1)
    assert np.allclose(reflected, sign * base, rtol=1e-6, atol=1e-9)


def test_closed_case_interpolate_matches_direct_beyond_halfpi():
    K, beta, l = 1, 40, 4
    x = np.linspace(0.6 * np.pi, 1.3 * np.pi, 25)  # outside [0, pi/2]
    interp = hyperspherical_bessel_interpolate(K, beta, l, x)
    direct = hyperspherical_bessel_direct(K, beta, l, x)
    assert np.allclose(interp, direct, rtol=1e-4, atol=1e-7)


def test_direct_continuous_across_turning_point():
    """Open case: direct is continuous across xfwd (the bwd->fwd switch)."""
    K, beta, l = -1, 60.0, 30
    xfwd = np.arcsinh(np.sqrt(l * (l + 1.0)) / beta)
    x = np.linspace(xfwd - 0.02, xfwd + 0.02, 41)
    got = hyperspherical_bessel_direct(K, beta, l, x)
    # No jump: max adjacent difference is small relative to local scale.
    jumps = np.abs(np.diff(got))
    assert jumps.max() < 5e-3, jumps.max()
    # And it matches interpolation through the seam.
    interp = hyperspherical_bessel_interpolate(K, beta, l, x)
    assert np.allclose(got, interp, rtol=1e-3, atol=1e-6)


def test_scalar_l_returns_1d():
    out = hyperspherical_bessel_direct(0, 30.0, 3, np.linspace(0.1, 1.0, 10))
    assert out.ndim == 1 and out.shape == (10,)


def test_closed_case_rejects_noninteger_beta():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(1, 40.5, 3, np.array([0.5, 1.0]))


def test_closed_case_rejects_l_ge_beta():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(1, 10, 10, np.array([0.5, 1.0]))


def test_invalid_K_rejected():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(2, 30.0, 3, np.array([0.5]))


def test_derivatives_only_on_interpolate():
    with pytest.raises(TypeError):
        hyperspherical_bessel_direct(0, 30.0, 3, np.array([0.5]), derivatives=True)
```

- [ ] **Step 2: Run the full test suite**

Run: `cd python && python -m pytest test_hyperspherical.py -v`
Expected: PASS for all tests (Tasks 2–5).

Note: `test_derivatives_only_on_interpolate` passes because `hyperspherical_bessel_direct`
has no `derivatives` parameter, so Python raises `TypeError` for the unexpected keyword —
no code change needed.

- [ ] **Step 3: Commit**

```bash
git add python/test_hyperspherical.py
git commit -m "Add hyperspherical Bessel edge-case tests

Closed-case symmetry fold, turning-point continuity, shape/dtype
contract, and error paths.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Closeout — close PR #153 and issue #147

After the work is merged to master, close the stale PR and the issue.

**Files:** none (GitHub operations only).

- [ ] **Step 1: Confirm work is merged to master**

Run: `git log --oneline -6`
Expected: the five commits above are present on the integration branch / master.

- [ ] **Step 2: Close PR #153 with an explanatory comment**

Run:
```bash
gh pr comment 153 --body "Closing: the only correctness content here (the hit_the_ceiling / PhiL_plus_one bugfix) has long been in master, and the branch is based on 2023 master so it cannot be merged without reverting years of refactoring. The actual goal of issue #147 — exposing the hyperspherical Bessel functions through the classy wrapper — is delivered separately (hyperspherical_bessel_interpolate / _direct / _wkb), along with the OpenMP cleanup this PR intended."
gh pr close 153
```
Expected: PR #153 shows state CLOSED.

- [ ] **Step 3: Close issue #147**

Run:
```bash
gh issue comment 147 --body "Done: hyperspherical Bessel functions are now exposed in classy as module-level functions hyperspherical_bessel_interpolate, hyperspherical_bessel_direct, and hyperspherical_bessel_wkb. The dead OpenMP scaffolding in hyperspherical_HIS_create has also been removed. (std::thread parallelisation was considered out of scope.)"
gh issue close 147
```
Expected: issue #147 shows state CLOSED.

---

## Notes for the executor

- **Build cost:** every task that changes `classy.pyx` or C++ requires a full `pip install .`
  (a few minutes). Batch the build at the test step as shown.
- **`cclassy.pxd`** is regenerated on each build; never edit it by hand.
- **Brace balance** in Task 1 Step 7 is the one easy-to-break spot — if the build fails to
  compile right after Task 1, re-check that exactly one brace (the parallel-region closer)
  was removed.
- The `direct` and `wkb` functions deliberately have no `derivatives` parameter; only
  `interpolate` supports it (Φ only for the other two, per the spec).
