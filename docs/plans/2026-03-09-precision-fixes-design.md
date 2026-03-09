# Precision Fixes Port from class_public to CLASSpp

**Date:** 2026-03-09
**Status:** Approved

## Context

CLASSpp was forked from CLASS v2.9 in early 2020. Since then, class_public has introduced several numerical precision improvements through v3.3.4. This design covers porting the three most impactful precision fixes to CLASSpp.

## Prerequisite

The module rename PR (spectra -> harmonic, nonlinear -> fourier, matching class_public v3.x naming) must be merged before this work begins. This design assumes the rename is complete.

## Approach

Two independent PRs, one per fix, in order of complexity:

1. `fix/tensor-tca-initial-conditions`
2. `fix/full-limber-scheme`

Each PR includes its own validation against class_public v3.3.4 reference data.

**Note:** The cubic P(k) interpolation fix (class_public issue #142) was investigated and found to be unnecessary. CLASSpp already uses cubic spline interpolation of ln(P(k)) in ln(k) space at the C++ level (`nonlinear_module.cpp:810-817`), and the Python wrapper delegates directly to C++ without scipy interpolation. The class_public fix was Python-wrapper-only.

---

## PR1: Tensor TCA Initial Conditions Fix

**Origin:** class_public commit `e827ad32` (v3.3.0), by Cyril Pitrou and Thomas Tram.
**Impact:** Eliminates ~0.25% systematic error on tensor C_l.
**File:** `source/perturbations_module.cpp`

### Changes

#### 1. TCA switch-off initial conditions (~line 4656)

Correct coefficients when transitioning out of tight-coupling for tensor modes:

- `delta_g`: `-4/3 * gwdot/kappa'` -> `+sqrt(6)*4/3 * gwdot/kappa'`
- `pol0_g`: `+1/3 * gwdot/kappa'` -> `-sqrt(6)/3 * gwdot/kappa'`

#### 2. Tensor initial conditions with (k*tau)^2 corrections (~line 5500)

Add new code block (~40 lines):
- Compute free-streaming energy density `rho_fs` from ur + ncdm species
- Compute correction: `h_corr_2 = -h0 * (k^2 + 2K) / (6 + 8/5 * rho_fs/rho_r) * tau^2`
- Apply to `gw`, `gwdot`, `delta_ur`, and ncdm `psi0`

#### 3. Tensor pressure source (~line 7254)

Fix P coefficient in tight-coupling:
- `2/5*sqrt(6) * gwdot/kappa'` -> `-1/3 * gwdot/kappa'`

#### 4. Print variables diagnostics (~line 7724)

Match corrected coefficients in debug output.

---

## PR2: Full Limber Scheme

**Origin:** class_public commit `a7304a5e` (v3.2.2).
**Impact:** Improved accuracy for CMB lensing C_l^phiphi at high l.

### New data structures (`transfer.h`)

- `short do_lcmb_full_limber` -- flag
- `size_t q_size_limber` -- separate q-grid size for Limber path
- `double *q_limber` -- q-values for Limber path
- `double **k_limber` -- k-values per mode for Limber path
- `double **transfer_limber` -- transfer functions on the Limber grid

### New precision parameters (`precisions.h`)

- `q_logstep_limber = 1.025`
- `k_max_limber_over_l_max_scalars = 0.001`
- `perturbations_sampling_boost_above_age_fraction = 0.9`
- Modify default `k_max_tau0_over_l_max`: 2.4 -> 1.8

### New input parameter (`input_module.cpp`)

- `want_lcmb_full_limber` (default: true)

### Transfer module changes (`transfer_module.cpp`)

- New method: `TransferModule::transfer_get_q_limber_list()` to generate the coarser Limber q-grid
- Modified main loop: iterate over `max(q_size, q_size_limber)`, computing normal transfer functions for all types on the standard grid, and lcmb-only on the Limber grid
- Modified `transfer_compute_for_each_q()` and `transfer_compute_for_each_l()`: accept `use_full_limber` parameter to select which q/k arrays to use
- Memory allocation and deallocation for the new arrays

### Spectra/harmonic module changes (`spectra_module.cpp`)

- Allocate `cl_integrand_limber` alongside `cl_integrand`
- For l > l_switch_limber and phiphi spectrum: integrate using Limber transfer functions and the Limber q-grid with full spline integration
- For all other spectra/l-values: use the standard path unchanged

### Perturbations module changes (`perturbations_module.cpp`)

- Extend k_max when full Limber is enabled
- Boost time sampling by 2x when tau > 0.9 * conformal_age

### Design decisions

- New methods on `TransferModule` and `SpectraModule` classes (following CLASSpp conventions)
- New members on data structs
- Match existing CLASSpp allocation patterns

---

## Validation Script

**File:** `scripts/validate_precision_fixes.py`

Created in PR1, extended in PR2 and PR3. Uses matplotlib (matching repo notebook conventions).

### Interface

```
python validate_precision_fixes.py [--tensor-tca] [--cubic-pk] [--limber] [--all]
    [--generate-reference]  # run class_public to produce reference data
    [--reference-dir PATH]  # where reference .npz files live
```

### Reference data

Generated from class_public v3.3.4 with fixed LCDM parameters (r=0.01 for tensors). Saved as `.npz` files in `scripts/reference_data/` and committed to the repo.

### Acceptance criteria

| Test | Comparison | Criterion |
|------|-----------|-----------|
| `--tensor-tca` | Tensor B-mode C_l, l=2..500 | Relative diff to class_public < 0.05% |
| `--limber` | C_l^phiphi, l=2..3000 | Relative diff to class_public < 0.1% at l>500 |

### Output

- Prints pass/fail per test case
- Saves matplotlib comparison plots to `scripts/validation_plots/`
- Three panels per test: class_public reference, CLASSpp output, relative difference
