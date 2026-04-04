#!/usr/bin/env python3
"""Compare current profile output against cached reference spectra.
Usage: python compare_profiles.py
Exits with non-zero status if any l*(l+1)*Cl/2pi differs by more than tol.
"""
import numpy as np
import sys

TOL = 1e-4  # relative tolerance on each Cl

scenarios = [
    ("output/profile_lcdm_ncdm_cl_lensed.dat",
     "output/reference/profile_lcdm_ncdm_cl_lensed.dat",
     "LCDM+NCDM (flat)"),
    ("output/profile_lcdm_ncdm_curved_cl_lensed.dat",
     "output/reference/profile_lcdm_ncdm_curved_cl_lensed.dat",
     "LCDM+NCDM (curved)"),
]

all_ok = True
for current_path, ref_path, label in scenarios:
    cur = np.loadtxt(current_path)
    ref = np.loadtxt(ref_path)
    if cur.shape != ref.shape:
        print(f"FAIL [{label}]: shape mismatch {cur.shape} vs {ref.shape}")
        all_ok = False
        continue
    # columns: l, TT, EE, BB, TE, ...
    ls = ref[:, 0]
    max_reldiff = 0.0
    worst_col = 0
    worst_l = 0
    for col in range(1, ref.shape[1]):
        denom = np.abs(ref[:, col])
        mask = denom > 1e-30
        if not mask.any():
            continue
        reldiff = np.abs(cur[mask, col] - ref[mask, col]) / denom[mask]
        idx = np.argmax(reldiff)
        if reldiff[idx] > max_reldiff:
            max_reldiff = reldiff[idx]
            worst_col = col
            worst_l = int(ls[mask][idx])
    if max_reldiff > TOL:
        print(f"FAIL [{label}]: max relative diff = {max_reldiff:.2e} at l={worst_l}, col={worst_col}")
        all_ok = False
    else:
        print(f"OK   [{label}]: max relative diff = {max_reldiff:.2e}")

sys.exit(0 if all_ok else 1)
