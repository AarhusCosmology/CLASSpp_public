# Precision Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Port two precision fixes from class_public (v3.2.2–v3.3.0) to CLASSpp: tensor TCA initial conditions and full Limber scheme.

**Architecture:** Two independent PRs targeting `source/perturbations_module.cpp` (PR1) and `source/transfer_module.{h,cpp}` + `source/spectra_module.cpp` + `include/precisions.h` (PR2). A shared validation script compares CLASSpp output against class_public v3.3.4 reference data.

**Tech Stack:** C++17, Cython, Python/matplotlib for validation. GNU Make build system.

**Prerequisite:** The module rename PR (spectra→harmonic, nonlinear→fourier) is assumed to be merged. If not yet merged, the file names in this plan use CLASSpp's current names (`spectra_module`, `nonlinear_module`) and should be adjusted after the rename.

**Note on cubic P(k) interpolation:** Investigation confirmed CLASSpp already uses cubic spline interpolation of ln(P(k)) in ln(k) space at the C++ level (`nonlinear_module.cpp:810-817`), and the Python wrapper delegates directly to C++ without scipy interpolation. No fix needed — CLASSpp is already correct.

---

## Task 1: Generate class_public Reference Data

**Files:**
- Create: `scripts/generate_reference_data.py`
- Create: `scripts/reference_data/tensor_bb_cl.npz`
- Create: `scripts/reference_data/lensing_phiphi_cl.npz`

This task generates reference output from class_public v3.3.4 that all subsequent validation compares against. Run it once; the `.npz` files are committed to the repo.

**Step 1: Write the reference data generation script**

```python
"""Generate reference data from class_public for validating precision fixes.

Requires class_public's classy wrapper to be importable (pip install from
/Users/ebh/Documents/study/cosmology/class_public).

Usage:
    python scripts/generate_reference_data.py
"""
import numpy as np
import sys

try:
    sys.path.insert(0, '/Users/ebh/Documents/study/cosmology/class_public/python')
    from classy import Class
except ImportError:
    print("ERROR: class_public's classy not importable. Install it first:")
    print("  cd /path/to/class_public && pip install .")
    sys.exit(1)

import os
outdir = os.path.join(os.path.dirname(__file__), 'reference_data')
os.makedirs(outdir, exist_ok=True)

# Shared LCDM baseline
baseline = {
    'output': 'tCl,pCl,lCl',
    'lensing': 'yes',
    'H0': 67.36,
    'omega_b': 0.02237,
    'omega_cdm': 0.12,
    'A_s': 2.1e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544,
}

# --- Tensor B-mode reference ---
print("Computing tensor B-mode reference...")
cosmo = Class()
tensor_params = {**baseline, 'modes': 's,t', 'r': 0.01, 'l_max_tensors': 500}
cosmo.set(tensor_params)
cosmo.compute()
cl = cosmo.raw_cl(500)
np.savez(os.path.join(outdir, 'tensor_bb_cl.npz'),
         ell=cl['ell'], bb=cl['bb'], params=tensor_params)
cosmo.struct_cleanup()
cosmo.empty()
print(f"  Saved tensor_bb_cl.npz (l=0..500)")

# --- Lensing potential reference ---
print("Computing lensing potential reference...")
cosmo = Class()
lensing_params = {**baseline, 'l_max_scalars': 3000}
cosmo.set(lensing_params)
cosmo.compute()
cl = cosmo.raw_cl(3000)
cl_lens = cosmo.lensed_cl(3000)
np.savez(os.path.join(outdir, 'lensing_phiphi_cl.npz'),
         ell=cl['ell'], pp=cl['pp'], tt_lensed=cl_lens['tt'],
         params=lensing_params)
cosmo.struct_cleanup()
cosmo.empty()
print(f"  Saved lensing_phiphi_cl.npz (l=0..3000)")

print("Done. Reference data saved to", outdir)
```

**Step 2: Run the script to generate reference data**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && python scripts/generate_reference_data.py`
Expected: Two `.npz` files created in `scripts/reference_data/`

**Step 3: Commit**

```bash
git add scripts/generate_reference_data.py scripts/reference_data/
git commit -m "Add reference data generation script for precision fix validation"
```

---

## Task 2: Write Validation Script Skeleton

**Files:**
- Create: `scripts/validate_precision_fixes.py`

This is the shared validation script. It starts with tensor TCA support and will be extended for Limber in Task 8.

**Step 1: Write the validation script**

```python
"""Validate CLASSpp precision fixes against class_public reference data.

Usage:
    python scripts/validate_precision_fixes.py --tensor-tca
    python scripts/validate_precision_fixes.py --limber
    python scripts/validate_precision_fixes.py --all
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(SCRIPT_DIR, 'reference_data')
PLOT_DIR = os.path.join(SCRIPT_DIR, 'validation_plots')

def load_reference(filename):
    path = os.path.join(REF_DIR, filename)
    if not os.path.exists(path):
        print(f"ERROR: Reference file not found: {path}")
        print("Run generate_reference_data.py first.")
        sys.exit(1)
    return np.load(path, allow_pickle=True)

def validate_tensor_tca():
    """Compare tensor B-mode Cl against class_public reference."""
    from classy import Class

    ref = load_reference('tensor_bb_cl.npz')
    ref_ell = ref['ell']
    ref_bb = ref['bb']
    params = ref['params'].item()

    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    cl = cosmo.raw_cl(500)
    cpp_ell = cl['ell']
    cpp_bb = cl['bb']
    cosmo.struct_cleanup()
    cosmo.empty()

    # Compare for l >= 2 where Cl is nonzero
    mask = ref_ell >= 2
    ell = ref_ell[mask]
    rel_diff = np.where(ref_bb[mask] != 0,
                        (cpp_bb[mask] - ref_bb[mask]) / ref_bb[mask],
                        0.0)

    max_diff = np.max(np.abs(rel_diff))
    passed = max_diff < 5e-4  # 0.05%

    # Plot
    os.makedirs(PLOT_DIR, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    factor = ell * (ell + 1) / (2 * np.pi)

    axes[0].semilogy(ell, factor * ref_bb[mask], label='class_public v3.3.4')
    axes[0].semilogy(ell, factor * cpp_bb[mask], '--', label='CLASSpp')
    axes[0].set_ylabel(r'$\ell(\ell+1)C_\ell^{BB}/(2\pi)$')
    axes[0].legend()
    axes[0].set_title('Tensor B-mode power spectrum')

    axes[1].plot(ell, rel_diff * 100)
    axes[1].axhline(0, color='gray', ls=':')
    axes[1].axhline(0.05, color='r', ls='--', alpha=0.5, label='0.05% threshold')
    axes[1].axhline(-0.05, color='r', ls='--', alpha=0.5)
    axes[1].set_xlabel(r'$\ell$')
    axes[1].set_ylabel('Relative difference (%)')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'tensor_tca_validation.png'), dpi=150)
    plt.close()

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Tensor TCA: max relative difference = {max_diff*100:.4f}% (threshold: 0.05%)")
    return passed

def validate_limber():
    """Compare lensing potential Cl^phiphi against class_public reference."""
    from classy import Class

    ref = load_reference('lensing_phiphi_cl.npz')
    ref_ell = ref['ell']
    ref_pp = ref['pp']
    params = ref['params'].item()

    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    cl = cosmo.raw_cl(3000)
    cpp_ell = cl['ell']
    cpp_pp = cl['pp']
    cosmo.struct_cleanup()
    cosmo.empty()

    # Compare for l >= 2
    mask = ref_ell >= 2
    ell = ref_ell[mask]
    rel_diff = np.where(ref_pp[mask] != 0,
                        (cpp_pp[mask] - ref_pp[mask]) / ref_pp[mask],
                        0.0)

    # Check high-l accuracy (l > 500)
    high_l_mask = ell > 500
    max_diff_high_l = np.max(np.abs(rel_diff[high_l_mask]))
    passed = max_diff_high_l < 1e-3  # 0.1%

    # Plot
    os.makedirs(PLOT_DIR, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    factor = (ell * (ell + 1))**2 / (2 * np.pi)

    axes[0].loglog(ell, factor * ref_pp[mask], label='class_public v3.3.4')
    axes[0].loglog(ell, factor * cpp_pp[mask], '--', label='CLASSpp')
    axes[0].set_ylabel(r'$[\ell(\ell+1)]^2 C_\ell^{\phi\phi}/(2\pi)$')
    axes[0].legend()
    axes[0].set_title('CMB lensing potential power spectrum')

    axes[1].semilogx(ell, rel_diff * 100)
    axes[1].axhline(0, color='gray', ls=':')
    axes[1].axhline(0.1, color='r', ls='--', alpha=0.5, label='0.1% threshold')
    axes[1].axhline(-0.1, color='r', ls='--', alpha=0.5)
    axes[1].set_xlabel(r'$\ell$')
    axes[1].set_ylabel('Relative difference (%)')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'limber_validation.png'), dpi=150)
    plt.close()

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Full Limber: max relative difference at l>500 = {max_diff_high_l*100:.4f}% (threshold: 0.1%)")
    return passed

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Validate CLASSpp precision fixes')
    parser.add_argument('--tensor-tca', action='store_true')
    parser.add_argument('--limber', action='store_true')
    parser.add_argument('--all', action='store_true')
    args = parser.parse_args()

    if not any([args.tensor_tca, args.limber, args.all]):
        parser.print_help()
        sys.exit(1)

    results = []
    if args.tensor_tca or args.all:
        results.append(('tensor_tca', validate_tensor_tca()))
    if args.limber or args.all:
        results.append(('limber', validate_limber()))

    all_passed = all(r[1] for r in results)
    sys.exit(0 if all_passed else 1)
```

**Step 2: Commit**

```bash
git add scripts/validate_precision_fixes.py
git commit -m "Add validation script for precision fix comparison against class_public"
```

---

## Task 3: Run Baseline Validation (Before Fix)

Before making any code changes, run the tensor TCA validation to capture the current (~0.25%) discrepancy.

**Step 1: Build CLASSpp and install the wrapper**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && pip install --no-build-isolation .`
Expected: Successful build

**Step 2: Run tensor TCA validation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && python scripts/validate_precision_fixes.py --tensor-tca`
Expected: `[FAIL] Tensor TCA: max relative difference = ~0.25%` — this confirms the known discrepancy exists.

**Step 3: Save the baseline plot**

The plot at `scripts/validation_plots/tensor_tca_validation.png` shows the pre-fix discrepancy. Keep it for before/after comparison.

---

## Task 4: Fix TCA Switch-Off Initial Conditions

**Files:**
- Modify: `source/perturbations_module.cpp:4566-4570`

Two locations need fixing: the zeroed-out TCA switch-off at lines 4566-4570, and the old-formula TCA switch-off at lines 4656-4658.

**Step 1: Fix the first TCA switch-off block (tensor, lines 4566-4570)**

The current code zeros out `delta_g` and `pol0_g` with a `//TBC` comment. Replace with the corrected Pitrou formulas.

Find (lines 4566-4570):
```cpp
        ppv->y[ppv->index_pt_delta_g] = 0.0; //TBC
        //-4./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];

        ppv->y[ppv->index_pt_pol0_g] = 0.0; //TBC
        //1./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];
```

Replace with:
```cpp
        /* Tight-coupling relation F_0^(2) = sqrt(6)*4/3 * H'/kappa' (C. Pitrou, class_public v3.3.0) */
        ppv->y[ppv->index_pt_delta_g] = _SQRT6_*4./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];

        /* Tight-coupling relation G_0^(2) = -sqrt(2/3) * H'/kappa' (C. Pitrou, class_public v3.3.0) */
        ppv->y[ppv->index_pt_pol0_g] = -_SQRT6_/3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];
```

**Step 2: Fix the second TCA switch-off block (scalar+tensor shared, lines 4656-4658)**

Find (lines 4656-4658):
```cpp
        ppv->y[ppv->index_pt_delta_g] = -4./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];

        ppv->y[ppv->index_pt_pol0_g] = 1./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];
```

Replace with:
```cpp
        /* Tight-coupling relation F_0^(2) = sqrt(6)*4/3 * H'/kappa' (C. Pitrou, class_public v3.3.0) */
        ppv->y[ppv->index_pt_delta_g] = _SQRT6_*4./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];

        /* Tight-coupling relation G_0^(2) = -sqrt(2/3) * H'/kappa' (C. Pitrou, class_public v3.3.0) */
        ppv->y[ppv->index_pt_pol0_g] = -_SQRT6_/3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];
```

**Important context:** Lines 4556-4570 are inside a tensor-mode block (`_tensors_` check). Lines 4647-4659 are inside a second tensor-mode TCA switch-off block. Both must be updated. Verify both are inside `if ((pa_old[ppw->index_ap_tca] == (int)tca_on) && (ppw->approx[ppw->index_ap_tca] == (int)tca_off))` conditionals for tensor modes.

**Step 3: Verify the file compiles**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`
Expected: Successful compilation

**Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "Fix tensor TCA switch-off initial conditions (Pitrou, class_public v3.3.0)

Correct coefficients for delta_g and pol0_g when switching off tight-coupling
approximation for tensor modes. Old formulas lacked sqrt(6) factor and had
wrong signs.

Reference: class_public commit e827ad32, C. Pitrou & T. Tram."
```

---

## Task 5: Add (k*tau)^2 Corrections to Tensor Initial Conditions

**Files:**
- Modify: `source/perturbations_module.cpp:4731-4739` (variable declarations)
- Modify: `source/perturbations_module.cpp:5402-5404` (after tensor gw initialization)

**Step 1: Add variable declarations**

Find (lines 4737-4739):
```cpp
  double delta_tot;
  double velocity_tot;
  double s2_squared;
```

Replace with:
```cpp
  double delta_tot;
  double velocity_tot;
  double s2_squared;
  double h_corr_2, rho_fs; /* Corrections to tensor mode initial conditions (C. Pitrou, class_public v3.3.0) */
```

**Step 2: Add the (k*tau)^2 correction block after tensor gw initialization**

Find (lines 5395-5404):
```cpp
    if (pba->sgnK == -1) {
      if (k*k+3*pba->K >= 0.) {
        ppw->pv->y[ppw->pv->index_pt_gw] *= sqrt(tanh(_PI_/2.*sqrt(k2+3*pba->K)/sqrt(-pba->K)));
      }
      else {
        ppw->pv->y[ppw->pv->index_pt_gw] = 0.;
      }
    }

  }
```

Replace with:
```cpp
    if (pba->sgnK == -1) {
      if (k*k+3*pba->K >= 0.) {
        ppw->pv->y[ppw->pv->index_pt_gw] *= sqrt(tanh(_PI_/2.*sqrt(k2+3*pba->K)/sqrt(-pba->K)));
      }
      else {
        ppw->pv->y[ppw->pv->index_pt_gw] = 0.;
      }
    }

    /**
     * Corrections of order (k*tau)^2 for h, order (k*tau) for h'.
     * Not including them results in 0.2-0.3% errors on final tensor Cl.
     * (Credits C. Pitrou, class_public v3.3.0)
     *
     * h = h0 + const*(k*tau)^2
     * const = -(1+2K/k^2)/(6 + 8/5 * sum_fs(3*P_fs)/rho_r) * h0
     *
     * where the sum is over free-streaming species (ur and ncdm).
     */

    /* Build energy density of free-streaming particles (3*P_fs) */
    rho_fs = 0.;

    if (pba->has_ur == _TRUE_)
      rho_fs += ppw->pvecback[background_module_->index_bg_rho_ur_];

    if (pba->has_ncdm == _TRUE_) {
      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++) {
        rho_fs += 3.*ppw->pvecback[background_module_->index_bg_p_ncdm1_ + n_ncdm];
      }
    }

    /* Correct h and h' */
    h_corr_2 = -ppw->pv->y[ppw->pv->index_pt_gw]*(k2 + 2.*pba->K)/(6. + 8./5.*rho_fs/rho_r)*tau*tau;
    ppw->pv->y[ppw->pv->index_pt_gw] += h_corr_2;
    ppw->pv->y[ppw->pv->index_pt_gwdot] = 2.*h_corr_2/tau;

    /* Set ur quadrupole F_0^(2) to order tau^2 value */
    if (evolve_tensor_ur_ == _TRUE_)
      ppw->pv->y[ppw->pv->index_pt_delta_ur] = _SQRT6_*h_corr_2;

    /* Set ncdm psi0 to order tau^2 value */
    if (evolve_tensor_ncdm_ == _TRUE_) {
      idx = ppw->pv->index_pt_psi0_ncdm1;
      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++) {
        for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++) {
          ppw->pv->y[idx] = _SQRT6_*h_corr_2*(-0.25*ncdm_->dlnf0_dlnq_ncdm_[n_ncdm][index_q]);
          idx += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
        }
      }
    }

  }
```

**Note on CLASSpp naming:** class_public uses `pba->index_bg_p_ncdm1` while CLASSpp uses `background_module_->index_bg_p_ncdm1_` (trailing underscore, module pointer). Similarly `pba->dlnf0_dlnq_ncdm` becomes `ncdm_->dlnf0_dlnq_ncdm_`. The `ncdm_` pointer is accessed via the NonColdDarkMatter helper. Verify the exact field name by checking other ncdm code in the same file (e.g., line 5239: `ncdm_->dlnf0_dlnq_ncdm_[n_ncdm][index_q]`). Also verify `index_bg_p_ncdm1_` exists in BackgroundModule — if it's named differently, grep for `p_ncdm` in `background_module.h`.

**Step 3: Verify compilation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`
Expected: Successful compilation

**Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "Add (k*tau)^2 corrections to tensor initial conditions (Pitrou, class_public v3.3.0)

Correct gravitational wave amplitude and derivative at order (k*tau)^2,
and set ur/ncdm quadrupoles to their order-tau^2 values. Without this,
tensor Cl has 0.2-0.3% systematic errors.

Reference: class_public commit e827ad32, C. Pitrou & T. Tram."
```

---

## Task 6: Fix Tensor Pressure Source and Print Diagnostics

**Files:**
- Modify: `source/perturbations_module.cpp:7254` (tensor pressure P in TCA)
- Modify: `source/perturbations_module.cpp:7724-7727` (print diagnostics)

**Step 1: Fix the tensor pressure source P coefficient**

Find (line 7254):
```cpp
        P = 2./5.*_SQRT6_*y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_]; //TBC
```

Replace with:
```cpp
        /* TCA solution P^(2) = -1/3 * H'/kappa', valid in both hierarchies (C. Pitrou, class_public v3.3.0) */
        P = -1./3.*y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[thermodynamics_module_->index_th_dkappa_];
```

**Step 2: Fix the print diagnostics**

Find (lines 7724-7727):
```cpp
        delta_g = -4./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/pvecthermo[thermodynamics_module_->index_th_dkappa_]; //TBC
        shear_g = 0.;
        l4_g = 0.;
        pol0_g = 1./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/pvecthermo[thermodynamics_module_->index_th_dkappa_]; //TBC
```

Replace with:
```cpp
        /* Corrected TCA relations (C. Pitrou, class_public v3.3.0) */
        delta_g = 4./3.*_SQRT6_*ppw->pv->y[ppw->pv->index_pt_gwdot]/pvecthermo[thermodynamics_module_->index_th_dkappa_];
        shear_g = 0.;
        l4_g = 0.;
        pol0_g = -_SQRT6_/3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/pvecthermo[thermodynamics_module_->index_th_dkappa_];
```

**Step 3: Verify compilation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`
Expected: Successful compilation

**Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "Fix tensor pressure source and diagnostics for corrected TCA (Pitrou, class_public v3.3.0)

Update P coefficient from 2/5*sqrt(6) to -1/3 in tight-coupling, and
match delta_g/pol0_g diagnostics to the corrected formulas."
```

---

## Task 7: Validate and Finalize PR1

**Step 1: Rebuild CLASSpp with all tensor TCA fixes**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && pip install --no-build-isolation .`
Expected: Successful build

**Step 2: Run the existing CLASSpp test suite**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public/python && TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py`
Expected: All tests pass. TEST_LEVEL=2 includes tensor modes.

**Step 3: Run tensor TCA validation against class_public reference**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && python scripts/validate_precision_fixes.py --tensor-tca`
Expected: `[PASS] Tensor TCA: max relative difference < 0.05%`

**Step 4: If validation fails, debug**

Check the validation plot at `scripts/validation_plots/tensor_tca_validation.png`. Common issues:
- Wrong sign on a coefficient: check the formulas against class_public `perturbations.c:5196-5200`
- Missing `_SQRT6_` factor: verify `_SQRT6_` is defined in `include/common.h`
- `rho_fs` computation wrong: verify `index_bg_p_ncdm1_` field name

**Step 5: Commit validation results**

```bash
git add scripts/validation_plots/tensor_tca_validation.png
git commit -m "Add tensor TCA validation plot showing <0.05% agreement with class_public"
```

---

## Task 8: Add Precision Parameters for Full Limber Scheme

**Files:**
- Modify: `include/precisions.h`

**Step 1: Add the new precision parameters**

Find the existing Limber-related precision parameters (around line 395):
```cpp
  double l_switch_limber;
```

Add the new parameters nearby (after the existing Limber parameters):
```cpp
  double q_logstep_limber; /**< logarithmic step for q-grid in full Limber scheme */
  double k_max_limber_over_l_max_scalars; /**< ratio k_max/l_max_scalars for full Limber scheme */
  double perturbations_sampling_boost_above_age_fraction; /**< time fraction above which source sampling is boosted for lensing accuracy */
```

**Step 2: Set default values in InputModule**

Find the input module file where precision defaults are set. Grep for where `l_switch_limber` is initialized (in `source/input_module.cpp`).

Add after the existing Limber defaults:
```cpp
  ppr->q_logstep_limber = 1.025;
  ppr->k_max_limber_over_l_max_scalars = 0.001;
  ppr->perturbations_sampling_boost_above_age_fraction = 0.9;
```

**Step 3: Modify the existing default for `k_max_tau0_over_l_max`**

Find where `k_max_tau0_over_l_max` is set to `2.4` and change it to `1.8`. This is safe because the full Limber scheme compensates for the reduced k_max in the standard path.

**Step 4: Verify compilation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`
Expected: Successful compilation

**Step 5: Commit**

```bash
git add include/precisions.h source/input_module.cpp
git commit -m "Add precision parameters for full Limber scheme

New parameters: q_logstep_limber, k_max_limber_over_l_max_scalars,
perturbations_sampling_boost_above_age_fraction. Also reduce default
k_max_tau0_over_l_max from 2.4 to 1.8 (compensated by Limber scheme)."
```

---

## Task 9: Add Full Limber Data Structures to Transfer Module

**Files:**
- Modify: `source/transfer_module.h` (add new public members and private method)
- Modify: `source/transfer.h` (add new struct fields, if the transfer struct is defined here)

**Step 1: Add new public data members to TransferModule class**

In `source/transfer_module.h`, add after `index_q_flat_approximation_` (around line 51):

```cpp
  short do_lcmb_full_limber_; /**< flag: use full Limber scheme for CMB lensing? */
  int q_size_limber_;         /**< number of q values for full Limber scheme */
  double * q_limber_;         /**< q-grid for full Limber scheme */
  double ** k_limber_;        /**< k-grid for full Limber scheme, per mode */
  double ** transfer_limber_; /**< transfer functions on the Limber grid */
```

**Step 2: Add new private method declaration**

In the private section of `TransferModule` (around line 69), add:

```cpp
  int transfer_get_q_limber_list(double K, int sgnK);
```

**Step 3: Update method signatures**

Update `transfer_compute_for_each_q` signature (line 76) to add `short use_full_limber`:

```cpp
  int transfer_compute_for_each_q(int ** tp_of_tt, int index_q, int tau_size_max, double tau_rec, double *** sources, double *** sources_spline, double * window, struct transfer_workspace * ptw, short use_full_limber);
```

Update `transfer_compute_for_each_l` signature (line 87) to add `short use_full_limber`:

```cpp
  int transfer_compute_for_each_l(struct transfer_workspace * ptw, int index_q, int index_md, int index_ic, int index_tt, int index_l, double l, double q_max_bessel, radial_function_type radial_type, short use_full_limber);
```

**Step 4: Verify compilation (expect failures — implementations not yet updated)**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class 2>&1 | head -20`
Expected: Compilation errors from signature mismatches. This is expected — we fix them in the next tasks.

**Step 5: Commit**

```bash
git add source/transfer_module.h
git commit -m "Add full Limber data structures and method signatures to TransferModule

New members: do_lcmb_full_limber_, q_size_limber_, q_limber_, k_limber_,
transfer_limber_. New method: transfer_get_q_limber_list(). Updated
signatures for transfer_compute_for_each_q/l with use_full_limber flag."
```

---

## Task 10: Add Input Parameter and Initialization Flag

**Files:**
- Modify: `include/perturbations.h` (add `want_lcmb_full_limber` to perturbations struct)
- Modify: `source/input_module.cpp` (read the parameter and set default)

**Step 1: Add flag to perturbations struct**

In `include/perturbations.h`, find the struct definition and add:

```cpp
  short want_lcmb_full_limber; /**< flag to enable full Limber scheme for CMB lensing (default: _TRUE_) */
```

**Step 2: Set default in input module**

In `source/input_module.cpp`, where perturbation defaults are set, add:

```cpp
  ppt->want_lcmb_full_limber = _TRUE_;
```

**Step 3: Add parameter reading**

In `source/input_module.cpp`, near where other Limber/lensing parameters are read, add:

```cpp
  class_read_flag("want_lcmb_full_limber", ppt->want_lcmb_full_limber);
```

**Step 4: Commit**

```bash
git add include/perturbations.h source/input_module.cpp
git commit -m "Add want_lcmb_full_limber input parameter (default: true)"
```

---

## Task 11: Implement transfer_get_q_limber_list()

**Files:**
- Modify: `source/transfer_module.cpp` (add new method + call from init)

**Step 1: Implement the q_limber grid generation**

Add the new method implementation. Reference: class_public `transfer.c:1296-1367`.

```cpp
int TransferModule::transfer_get_q_limber_list(double K, int sgnK) {

  double q, q_min, q_max;
  int index_q;

  /* q_max is determined by k_max for the full limber scheme */
  q_max = perturbations_module_->k_[index_md_scalars_]
          [perturbations_module_->k_size_[index_md_scalars_] - 1];
  if (sgnK != 0)
    q_max = sqrt(q_max*q_max + K); /* q^2 = k^2 + K for scalars */

  /* q_min same as standard grid */
  q_min = q_[0];

  /* Number of q values in logarithmic spacing */
  q_size_limber_ = (int)(log(q_max/q_min)/log(ppr->q_logstep_limber)) + 1;

  class_alloc(q_limber_, q_size_limber_*sizeof(double), error_message_);

  for (index_q = 0; index_q < q_size_limber_; index_q++) {
    q_limber_[index_q] = q_min*pow(ppr->q_logstep_limber, index_q);
  }

  /* Ensure last q value doesn't exceed q_max */
  if (q_limber_[q_size_limber_ - 1] > q_max)
    q_limber_[q_size_limber_ - 1] = q_max;

  return _SUCCESS_;
}
```

**Step 2: Call from transfer_init and set the flag**

In `transfer_init()` (near the beginning, after q_size_ and k_ are available), add:

```cpp
  /* Check whether full Limber scheme is needed */
  if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->want_lcmb_full_limber == _TRUE_)) {
    do_lcmb_full_limber_ = _TRUE_;
  } else {
    do_lcmb_full_limber_ = _FALSE_;
  }
```

And after `transfer_get_k_list()` is called, add:

```cpp
  if (do_lcmb_full_limber_ == _TRUE_) {
    class_call(transfer_get_q_limber_list(K, sgnK), error_message_, error_message_);

    /* Allocate k_limber arrays */
    class_alloc(k_limber_, md_size_*sizeof(double*), error_message_);
    for (int index_md = 0; index_md < md_size_; index_md++) {
      class_alloc(k_limber_[index_md], q_size_limber_*sizeof(double), error_message_);
      for (int iq = 0; iq < q_size_limber_; iq++) {
        if (sgnK == 0)
          k_limber_[index_md][iq] = q_limber_[iq];
        else
          k_limber_[index_md][iq] = sqrt(q_limber_[iq]*q_limber_[iq] - K);
      }
    }

    /* Allocate transfer_limber array */
    class_alloc(transfer_limber_, md_size_*sizeof(double*), error_message_);
    for (int index_md = 0; index_md < md_size_; index_md++) {
      class_alloc(transfer_limber_[index_md],
                  perturbations_module_->ic_size_[index_md]
                  * tt_size_[index_md]
                  * l_size_[index_md]
                  * q_size_limber_
                  * sizeof(double),
                  error_message_);
    }
  }
```

**Step 3: Add deallocation in transfer_free()**

In `transfer_free()`, add before the existing free calls:

```cpp
  if (do_lcmb_full_limber_ == _TRUE_) {
    for (int index_md = 0; index_md < md_size_; index_md++) {
      free(transfer_limber_[index_md]);
      free(k_limber_[index_md]);
    }
    free(q_limber_);
    free(k_limber_);
    free(transfer_limber_);
  }
```

**Step 4: Verify compilation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`
Expected: May still fail due to signature mismatches in compute functions. Fix in next task.

**Step 5: Commit**

```bash
git add source/transfer_module.cpp
git commit -m "Implement transfer_get_q_limber_list and Limber array allocation/deallocation"
```

---

## Task 12: Update Transfer Computation for Dual-Path Execution

**Files:**
- Modify: `source/transfer_module.cpp` (main q-loop, compute_for_each_q, compute_for_each_l)

This is the core of the Limber scheme. The main loop iterates over `max(q_size_, q_size_limber_)` and dispatches to either the normal path or the Limber path.

**Step 1: Update the main q-loop**

Find the main loop in `transfer_init()` (around line 291):
```cpp
for (index_q = 0; index_q < q_size_; index_q++)
```

Replace with a dual-path loop. The exact structure depends on the existing code — study the surrounding context (thread pool dispatch, HIS update, etc.) carefully. The key logic is:

```cpp
int q_loop_max = (do_lcmb_full_limber_ == _TRUE_) ?
    MAX(q_size_, q_size_limber_) : q_size_;

for (index_q = 0; index_q < q_loop_max; index_q++) {

    /* Normal path */
    if (index_q < q_size_) {
        class_call(transfer_update_HIS(ptw, index_q, tau0), error_message_, error_message_);
        class_call(transfer_compute_for_each_q(..., index_q, ..., ptw, _FALSE_),
                   error_message_, error_message_);
    }

    /* Full Limber path */
    if (do_lcmb_full_limber_ == _TRUE_ && index_q < q_size_limber_) {
        class_call(transfer_compute_for_each_q(..., index_q, ..., ptw, _TRUE_),
                   error_message_, error_message_);
    }
}
```

**Step 2: Update `transfer_compute_for_each_q` implementation**

Add `short use_full_limber` parameter. At the top, select q and k based on the flag:

```cpp
double q, k;
if (use_full_limber == _FALSE_) {
    q = q_[index_q];
    k = k_[index_md][index_q];
} else {
    q = q_limber_[index_q];
    k = k_limber_[index_md][index_q];
}
```

In the inner loop over transfer types, only compute lcmb type in Limber mode:

```cpp
if (use_full_limber == _TRUE_ && index_tt != index_tt_lcmb_)
    continue;
```

When storing the result, use the appropriate array:

```cpp
if (use_full_limber == _FALSE_) {
    transfer_[index_md][((index_ic * tt_size_[index_md] + index_tt)
                         * l_size_[index_md] + index_l)
                        * q_size_ + index_q] = transfer_function;
} else {
    transfer_limber_[index_md][((index_ic * tt_size_[index_md] + index_tt)
                                * l_size_[index_md] + index_l)
                               * q_size_limber_ + index_q] = transfer_function;
}
```

**Step 3: Update `transfer_compute_for_each_l` implementation**

Add `short use_full_limber` parameter. Select q/k and force Limber:

```cpp
double q, k;
short use_limber;
if (use_full_limber == _FALSE_) {
    q = q_[index_q];
    k = k_[index_md][index_q];
    class_call(transfer_use_limber(..., &use_limber), error_message_, error_message_);
} else {
    q = q_limber_[index_q];
    k = k_limber_[index_md][index_q];
    use_limber = _TRUE_;
}
```

**Step 4: Update all call sites**

Every existing call to `transfer_compute_for_each_q` and `transfer_compute_for_each_l` must pass `_FALSE_` as the `use_full_limber` argument (preserving existing behavior). Search for all call sites and update them.

**Step 5: Verify compilation and run tests**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`
Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && pip install --no-build-isolation . && cd python && python -m pytest -v -m test_scenario test_class.py`
Expected: Compiles and existing tests still pass (no behavioral change yet for standard runs).

**Step 6: Commit**

```bash
git add source/transfer_module.cpp
git commit -m "Implement dual-path q-loop for full Limber scheme in transfer module

Normal path unchanged. Full Limber path computes only lcmb transfer type
on a separate coarser q-grid, storing results in transfer_limber_ array."
```

---

## Task 13: Add Perturbations Module Changes for Limber

**Files:**
- Modify: `source/perturbations_module.cpp` (k_max extension, time sampling boost)

**Step 1: Extend k_max when full Limber is enabled**

Find where `k_max` for scalars is determined (grep for `k_max` assignments or `k_size_cl` in the perturbations module). Add:

```cpp
if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->want_lcmb_full_limber == _TRUE_)) {
    k_max = MAX(k_max, ppr->k_max_limber_over_l_max_scalars * ppt->l_scalar_max);
}
```

**Step 2: Add time sampling boost for lensing accuracy**

Find the section where `timescale_source` is computed for the adaptive time stepping (in the time integration loop). Add:

```cpp
/* Boost time sampling at late times for lensing accuracy (class_public v3.2.2) */
if (tau > pba->conformal_age * ppr->perturbations_sampling_boost_above_age_fraction) {
    timescale_source /= 2.;
}
```

**Note:** The exact variable names and location depend on CLASSpp's integration loop structure. Grep for `timescale_source` or the equivalent time step control variable.

**Step 3: Verify compilation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`

**Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "Extend k_max and boost time sampling for full Limber scheme

Increase k_max for scalars when CMB lensing + full Limber is enabled.
Refine time sampling by factor 2 at tau > 0.9 * conformal_age for
improved lensing source accuracy."
```

---

## Task 14: Add Full Limber Integration to Spectra Module

**Files:**
- Modify: `source/spectra_module.cpp` (Cl integration with Limber path)

This is the final computational piece: using `transfer_limber_` for the phiphi spectrum at l > l_switch_limber.

**Step 1: Allocate cl_integrand_limber**

In `spectra_cls()` (or equivalent), where `cl_integrand` is allocated, add:

```cpp
double * cl_integrand_limber = NULL;
if (transfer_module_->do_lcmb_full_limber_ == _TRUE_) {
    class_alloc(cl_integrand_limber,
                transfer_module_->q_size_limber_ * cl_integrand_num_columns * sizeof(double),
                error_message_);
}
```

**Step 2: Fill cl_integrand_limber for phiphi at high l**

In `spectra_compute_cl()` (or equivalent), after the normal cl_integrand is filled, add a block:

```cpp
if (transfer_module_->do_lcmb_full_limber_ == _TRUE_
    && _scalars_
    && has_pp_ == _TRUE_
    && index_ct == index_ct_pp_
    && l > ppr->l_switch_limber) {

    for (int index_q = 0; index_q < transfer_module_->q_size_limber_; index_q++) {
        double k = transfer_module_->k_limber_[index_md][index_q];
        cl_integrand_limber[index_q * cl_integrand_num_columns + 0] = k;

        /* Get primordial spectrum at k */
        class_call(primordial_spectrum_at_k(..., k, primordial_pk), ...);

        /* Get transfer functions from Limber array */
        double trsf_ic1 = transfer_module_->transfer_limber_[index_md]
            [((index_ic1 * transfer_module_->tt_size_[index_md] + transfer_module_->index_tt_lcmb_)
              * transfer_module_->l_size_[index_md] + index_l)
             * transfer_module_->q_size_limber_ + index_q];
        double trsf_ic2 = transfer_module_->transfer_limber_[index_md]
            [((index_ic2 * transfer_module_->tt_size_[index_md] + transfer_module_->index_tt_lcmb_)
              * transfer_module_->l_size_[index_md] + index_l)
             * transfer_module_->q_size_limber_ + index_q];

        double factor = 4. * _PI_ / k;
        cl_integrand_limber[index_q * cl_integrand_num_columns + 1 + index_ct_pp_] =
            primordial_pk[index_ic1_ic2] * trsf_ic1 * trsf_ic2 * factor;
    }
}
```

**Step 3: Use Limber integrand for phiphi at high l**

In the integration section, select the appropriate integrand:

```cpp
double * integrand;
int num_k;
if (transfer_module_->do_lcmb_full_limber_ == _TRUE_
    && _scalars_
    && has_pp_ == _TRUE_
    && index_ct == index_ct_pp_
    && l > ppr->l_switch_limber) {
    integrand = cl_integrand_limber;
    num_k = transfer_module_->q_size_limber_;
} else {
    integrand = cl_integrand;
    num_k = transfer_module_->q_size_;
}

/* Then use integrand and num_k for the spline + integration call */
```

**Step 4: Free cl_integrand_limber**

After the computation, add:

```cpp
if (cl_integrand_limber != NULL)
    free(cl_integrand_limber);
```

**Step 5: Verify compilation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && make -j class`

**Step 6: Commit**

```bash
git add source/spectra_module.cpp
git commit -m "Add full Limber integration path for phiphi spectrum in spectra module

For l > l_switch_limber and CMB lensing potential, use transfer_limber_
array with the coarser Limber q-grid for improved accuracy at high l."
```

---

## Task 15: Validate and Finalize PR2

**Step 1: Rebuild and install**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && pip install --no-build-isolation .`

**Step 2: Run existing test suite**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public/python && TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py`
Expected: All tests pass.

**Step 3: Run Limber validation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && python scripts/validate_precision_fixes.py --limber`
Expected: `[PASS] Full Limber: max relative difference at l>500 < 0.1%`

**Step 4: Run full validation**

Run: `cd /Users/ebh/Documents/study/cosmology/CLASSpp_public && python scripts/validate_precision_fixes.py --all`
Expected: Both tensor TCA and Limber pass.

**Step 5: If Limber validation fails, debug**

Common issues:
- Array indexing: verify the 4D index formula matches class_public's convention
- q_limber grid: verify it covers enough range (print `q_limber_[0]` and `q_limber_[q_size_limber_-1]`)
- transfer_limber_ not being filled: add debug prints in `transfer_compute_for_each_q` Limber path
- Spectra integration: verify the spline + integration call uses the right number of points

**Step 6: Commit validation results**

```bash
git add scripts/validation_plots/
git commit -m "Add Limber validation plot showing <0.1% agreement with class_public at l>500"
```
