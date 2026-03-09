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
    """Compare tensor B-mode Cl against class_public reference.

    CLASSpp was forked from class v2.9; class_public v3.3.4 has many other
    improvements beyond the TCA fix. We validate at low-l (l<=200) where the
    accumulated v2.9->v3.3.4 differences are small, so the TCA fix dominance
    is measurable. The pre-fix baseline reference captures the pre-fix state.
    """
    from classy import Class

    ref = load_reference('tensor_bb_cl.npz')
    ref_ell = ref['ell']
    ref_bb = ref['bb']
    params = ref['params'].item()

    # Also load pre-fix baseline if available
    baseline_path = os.path.join(REF_DIR, 'tensor_bb_cl_prefx.npz')
    has_baseline = os.path.exists(baseline_path)
    if has_baseline:
        baseline = np.load(baseline_path, allow_pickle=True)
        baseline_bb = baseline['bb']

    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    cl = cosmo.raw_cl(500)
    cpp_bb = cl['bb']
    cosmo.struct_cleanup()
    cosmo.empty()

    # Compare for l >= 2 where Cl is nonzero
    mask = ref_ell >= 2
    ell = ref_ell[mask]
    rel_diff = np.where(ref_bb[mask] != 0,
                        (cpp_bb[mask] - ref_bb[mask]) / ref_bb[mask],
                        0.0)

    # At low-l (2 <= l <= 200), the fix should bring agreement to < 1%
    # (residual is from other v2.9->v3.3.4 changes, not from TCA)
    low_l_mask = ell <= 200
    max_diff_low_l = np.max(np.abs(rel_diff[low_l_mask]))
    passed = max_diff_low_l < 1e-2  # 1%

    # If baseline available, verify improvement
    if has_baseline:
        baseline_diff = np.where(ref_bb[mask] != 0,
                                 (baseline_bb[mask] - ref_bb[mask]) / ref_bb[mask],
                                 0.0)
        max_baseline_low_l = np.max(np.abs(baseline_diff[low_l_mask]))
        improved = max_diff_low_l < max_baseline_low_l
        print(f"  Pre-fix low-l max diff: {max_baseline_low_l*100:.4f}%, post-fix: {max_diff_low_l*100:.4f}%")
        passed = passed and improved

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
    axes[1].axhline(1.0, color='r', ls='--', alpha=0.5, label='1% threshold (l<=200)')
    axes[1].axhline(-1.0, color='r', ls='--', alpha=0.5)
    if has_baseline:
        axes[1].plot(ell, baseline_diff * 100, ':', alpha=0.5, label='pre-fix')
    axes[1].set_xlabel(r'$\ell$')
    axes[1].set_ylabel('Relative difference (%)')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'tensor_tca_validation.png'), dpi=150)
    plt.close()

    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Tensor TCA: max relative difference at l<=200 = {max_diff_low_l*100:.4f}% (threshold: 1%)")
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
