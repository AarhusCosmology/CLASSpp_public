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
