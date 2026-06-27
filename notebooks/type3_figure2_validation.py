"""
Figure-2 validation for the Type-3 (arXiv:1604.04222) CDM-scalar-field
momentum-transfer coupling implemented in CLASSpp.

Reproduces qualitatively:
  - cl_ratio_beta.pdf : Cl^TT(beta) / Cl^TT(0)  vs  ell
  - pk_ratio_beta.pdf : P(k, beta)   / P(k, 0)   vs  k

Reference baseline choice:
  The Type-3 composite only exists when scf_veta is set; setting scf_veta=0
  drops the composite entirely and places the run on a different code path
  (plain scalar field, different ICs).  To keep a controlled comparison we
  use scf_veta=-1e-4 as the "zero-coupling" reference — this is the smallest
  value that still triggers the composite while keeping the momentum-transfer
  negligible.  The same strategy is used by test_type3_suppresses_growth.

IC scheme: frozen-IC injectable 1EXP potential, scf_parameters='1, 1.22'
  (V0_placeholder=1, lambda=1.22; tuning index 0 shoots V0 to hit Omega_scf).
  attractor_ic_scf='no' avoids the singular-matrix failure at negative beta.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import classy
import os

# ---------------------------------------------------------------------------
# Common parameters
# ---------------------------------------------------------------------------
COMMON = {
    'output': 'tCl mPk',
    'gauge': 'synchronous',
    'Omega_fld': 0,
    'Omega_scf': 0.7,
    'attractor_ic_scf': 'no',
    'scf_parameters': '1, 1.22',
    # CMB / P(k) output settings (string-keyed entries)
    'l_max_scalars': 2500,
    'P_k_max_1/Mpc': 1.0,
    'z_pk': 0,
    '100*theta_s': 1.0420,
    # Cosmological parameters
    'omega_b': 0.0223,
    'omega_cdm': 0.119,
    'n_s': 0.967,
    'ln10^{10}A_s': 3.06,
    'tau_reio': 0.07,
}

BETAS = [-1, -1e1, -1e2, -1e3, -1e4]

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def run(beta):
    """Set up and compute a Type-3 run; returns the Class instance."""
    c = classy.Class()
    c.set(dict(COMMON, scf_veta=beta))
    c.compute()
    return c


def savefig_tight(fig, path):
    fig.tight_layout()
    fig.savefig(path)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Reference run (effectively beta=0 on the Type-3 code path)
# ---------------------------------------------------------------------------
print("Computing reference (beta=-1e-4) …")
ref = run(-1e-4)

ell = np.arange(2, 2501)
cl0 = ref.raw_cl(2500)['tt'][2:]         # shape (2499,), l=2..2500
ks = np.logspace(-3, 0, 200)             # 1e-3 .. 1  [1/Mpc]
pk0 = np.array([ref.pk(k, 0) for k in ks])

# ---------------------------------------------------------------------------
# Loop over beta values, collect ratios, generate plots
# ---------------------------------------------------------------------------
fig_cl, ax_cl = plt.subplots(figsize=(7, 4))
fig_pk, ax_pk = plt.subplots(figsize=(7, 4))

failed = []
for beta in BETAS:
    print(f"Computing beta={beta:g} …")
    try:
        c = run(beta)
    except Exception as e:
        print(f"  FAILED: {e}")
        failed.append(beta)
        continue

    cl = c.raw_cl(2500)['tt'][2:]
    pk = np.array([c.pk(k, 0) for k in ks])

    cl_ratio = cl / cl0
    pk_ratio = pk / pk0

    # Diagnostic: report key ratio values
    low_l_mean = float(np.mean(cl_ratio[:8]))    # l=2..9  (ISW region)
    high_l_mean = float(np.mean(cl_ratio[100:300]))  # l~102..302
    pk_min = float(np.min(pk_ratio))
    pk_at_k01 = float(c.pk(0.1, 0) / ref.pk(0.1, 0))
    print(f"  Cl_ratio low-l(2-9)={low_l_mean:.4f}  high-l(102-302)={high_l_mean:.6f}")
    print(f"  pk_ratio min={pk_min:.4f}  pk(0.1)={pk_at_k01:.4f}")

    label = rf'$\beta={beta:g}$'
    ax_cl.semilogx(ell, cl_ratio, label=label)
    ax_pk.semilogx(ks, pk_ratio, label=label)

    c.struct_cleanup()
    c.empty()

# Reference cleanup
ref.struct_cleanup()
ref.empty()

if failed:
    print(f"\nWARNING: the following beta values failed to compute: {failed}")

# ---------------------------------------------------------------------------
# Finalise Cl ratio plot
# ---------------------------------------------------------------------------
ax_cl.axhline(1.0, color='k', lw=0.7, ls='--')
ax_cl.set_xlabel(r'$\ell$')
ax_cl.set_ylabel(r'$C_\ell^{TT}(\beta)\,/\,C_\ell^{TT}(0)$')
ax_cl.set_title('Type-3: CMB TT ratio (paper Fig. 2 left)')
ax_cl.legend(fontsize=8)
ax_cl.set_ylim(0.95, 1.20)

script_dir = os.path.dirname(os.path.abspath(__file__))
savefig_tight(fig_cl, os.path.join(script_dir, 'type3_cl_ratio.pdf'))
plt.close(fig_cl)

# ---------------------------------------------------------------------------
# Finalise P(k) ratio plot
# ---------------------------------------------------------------------------
ax_pk.axhline(1.0, color='k', lw=0.7, ls='--')
ax_pk.set_xlabel(r'$k\;[\mathrm{Mpc}^{-1}]$')
ax_pk.set_ylabel(r'$P(k,\beta)\,/\,P(k,0)$')
ax_pk.set_title('Type-3: matter power ratio (paper Fig. 2 right)')
ax_pk.legend(fontsize=8)

savefig_tight(fig_pk, os.path.join(script_dir, 'type3_pk_ratio.pdf'))
plt.close(fig_pk)

print("\nDone.  Check notebooks/type3_cl_ratio.pdf and notebooks/type3_pk_ratio.pdf")
