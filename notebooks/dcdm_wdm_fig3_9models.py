"""Reproduce Fig. 3 of arXiv:2606.14849 (Bencke, Lee & Kamionkowski) with the
CLASSpp `dcdm_wdm` species — all 9 published models.

The paper's Fig. 3 legend (read off arXiv-2606.14849v1/DESIDDM/cmb_3panel.png):

    z_d = 133, v_kick/c = 0.029      z_d = 9.5, v_kick/c = 0.015
    z_d =  80, v_kick/c = 0.024      z_d = 5.3, v_kick/c = 0.017
    z_d =  48, v_kick/c = 0.020      z_d = 2.7, v_kick/c = 0.024
    z_d =  28, v_kick/c = 0.017      z_d = 1.1, v_kick/c = 0.048
    z_d =  17, v_kick/c = 0.015

with f = 0.1 fixed and the Planck 2018 best-fit fiducial. The decay rate follows
from the paper's definition Gamma = H(z_d); on this background the 9 z_d values
map to an (almost exactly 10^(1/3)-spaced) lifetime grid 0.016 ... 7.6 Gyr.

Free-streaming markers: the paper's Eq. (1) integral and its quoted
matter-domination limit lambda_fs ~ 22,000 (v/c)/sqrt(1+z_d) h^-1 Mpc differ by
a factor of 2; pixel-measuring the published figure shows the dashed lines
follow the *latter* normalization, with k_fs = 2pi/lambda_fs in the P(k) panel
and l_fs = pi chi_*/lambda_fs in BOTH C_l panels (measured to ~10%, consistent
with the legend's 2-significant-figure rounding). We therefore use
lambda_fs = 2 x [Eq. (1) integral].

Spectra are cached in dcdm_wdm_fig3_cache.npz next to this script; delete the
file to force a recompute.
"""

import os
import time

import numpy as np
import matplotlib

if __name__ == "__main__":  # headless as a script; keep the notebook's backend otherwise
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from classy import Class

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "dcdm_wdm_fig3_cache.npz")

GYR = 977.792  # 1 Gyr^-1 in km/s/Mpc
C_KMS = 299792.458

# Planck 2018 best fit (TT,TE,EE+lowE+lensing), one 0.06 eV neutrino.
# The neutrino is declared dot-style: legacy N_ncdm collides with the ddm
# dot-instance parsing on this branch.
H_FID = 0.6736
OMEGA_DM = 0.1200 / H_FID**2
BASE = {
    "h": H_FID, "omega_b": 0.02237,
    "ln10^{10}A_s": 3.044, "n_s": 0.9649, "tau_reio": 0.0544,
    "N_ur": 2.0328, "nu.type": "ncdm_standard", "nu.m": 0.06,
    "output": "tCl,pCl,lCl,mPk", "lensing": "yes",
    "l_max_scalars": 3200, "P_k_max_1/Mpc": 5.0,
    "gauge": "synchronous", "evolver": "2",
}
F_DDM = 0.1

# (z_d, v_kick/c) — the published legend, dark-viridis to yellow.
MODELS = [(133, 0.029), (80, 0.024), (48, 0.020), (28, 0.017), (17, 0.015),
          (9.5, 0.015), (5.3, 0.017), (2.7, 0.024), (1.1, 0.048)]

KK = np.logspace(np.log10(7.5e-6), np.log10(3.5), 500)  # 1/Mpc
ELL = np.arange(2, 3001)
ELL_PP = np.arange(2, 1001)


def run(extra):
    c = Class()
    c.set({**BASE, **extra})
    t0 = time.time()
    c.compute()
    out = {
        "tt": c.lensed_cl(3000)["tt"][2:3001],
        "pp": c.raw_cl(1000)["pp"][2:1001],
        "pk": np.array([c.pk(k, 0.0) for k in KK]),
    }
    bg = c.get_background()
    out["bg_z"] = bg["z"][::-1]                # ascending z
    out["bg_H"] = bg["H [1/Mpc]"][::-1]        # 1/Mpc
    out["chi_star"] = c.get_current_derived_parameters(["ra_rec"])["ra_rec"]
    out["t"] = time.time() - t0
    c.struct_cleanup()
    c.empty()
    return out


def freestream(bg_z, bg_H, chi_star, z_d, vk):
    """(k_fs [1/Mpc], l_fs) of the paper's dashed lines: lambda = 2 x Eq. (1)."""
    zz = np.linspace(0.0, z_d, 4001)
    Hz = np.interp(zz, bg_z, bg_H)  # 1/Mpc
    lam = 2.0 * vk / (1.0 + z_d) * np.trapezoid((1.0 + zz) / Hz, zz)  # Mpc
    return 2.0 * np.pi / lam, np.pi * chi_star / lam


# ── compute (or load cached) spectra ─────────────────────────────────────────
if os.path.exists(CACHE):
    data = dict(np.load(CACHE))
    print(f"loaded cached spectra from {CACHE}")
else:
    data = {}
    print("LCDM reference ...")
    ref = run({"Omega_cdm": OMEGA_DM})
    print(f"  {ref['t']:.1f} s")
    for key in ("tt", "pp", "pk", "bg_z", "bg_H"):
        data["ref_" + key] = ref[key]
    data["chi_star"] = ref["chi_star"]

    for i, (z_d, vk) in enumerate(MODELS):
        Gamma = np.interp(z_d, ref["bg_z"], ref["bg_H"]) * C_KMS  # km/s/Mpc
        r = run({
            "Omega_cdm": OMEGA_DM * (1.0 - F_DDM),
            "ddm.type": "dcdm_wdm",
            "ddm.Gamma": Gamma,
            "ddm.vkick": vk,
            "ddm.Omega_ini": OMEGA_DM * F_DDM,
            "ddm.momenta_bins": 96,
        })
        for key in ("tt", "pp", "pk"):
            data[f"m{i}_{key}"] = r[key]
        data[f"m{i}_Gamma"] = np.array(Gamma)
        print(f"  z_d = {z_d:>5}, v = {vk:.3f}: Gamma^-1 = {GYR / Gamma:7.4f} Gyr, "
              f"{r['t']:.1f} s")
    np.savez(CACHE, **data)
    print(f"cached spectra in {CACHE}")

# ── the figure ───────────────────────────────────────────────────────────────
colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(MODELS)))
fig, (ax_pk, ax_tt, ax_pp) = plt.subplots(3, 1, figsize=(7.0, 10.6))

for i, ((z_d, vk), col) in enumerate(zip(MODELS, colors)):
    k_fs, l_fs = freestream(data["ref_bg_z"], data["ref_bg_H"],
                            float(data["chi_star"]), z_d, vk)
    lab = f"$z_d = {z_d:g}$, $v_\\mathrm{{kick}}/c = {vk:.3f}$"
    ax_pk.semilogx(KK, data[f"m{i}_pk"] / data["ref_pk"] - 1.0, color=col,
                   lw=1.8, label=lab)
    ax_tt.semilogx(ELL, data[f"m{i}_tt"] / data["ref_tt"] - 1.0, color=col, lw=1.8)
    ax_pp.semilogx(ELL_PP, data[f"m{i}_pp"] / data["ref_pp"] - 1.0, color=col, lw=1.8)
    ax_pk.axvline(k_fs, color=col, ls="--", lw=1.2, alpha=0.6)
    ax_tt.axvline(l_fs, color=col, ls="--", lw=1.2, alpha=0.6)
    ax_pp.axvline(l_fs, color=col, ls="--", lw=1.2, alpha=0.6)
    print(f"z_d = {z_d:>5}: Gamma^-1 = {GYR / float(data[f'm{i}_Gamma']):7.4f} Gyr, "
          f"k_fs = {k_fs:.4f} 1/Mpc, l_fs = {l_fs:.0f}")

for ax in (ax_pk, ax_tt, ax_pp):
    ax.axhline(0.0, color="k", ls="--", lw=1.0)

ax_pk.set_title("Matter Power Spectrum $P(k)$", fontsize=14)
ax_pk.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=14)
ax_pk.set_ylabel("$P(k)^\\mathrm{DDM}/P(k)^{\\Lambda\\mathrm{CDM}} - 1$", fontsize=14)
ax_pk.legend(loc="center left", fontsize=9.5, frameon=False)

ax_tt.set_title("Temperature Power Spectrum $C_\\ell^{TT}$", fontsize=14)
ax_tt.set_xlabel("$\\ell$", fontsize=14)
ax_tt.set_ylabel("$C_\\ell^{TT,\\mathrm{DDM}}/C_\\ell^{TT,\\Lambda\\mathrm{CDM}} - 1$",
                 fontsize=14)
ax_tt.set_xlim(2, 3000)

ax_pp.set_title("Lensing Power Spectrum $C_L^{\\phi\\phi}$", fontsize=14)
ax_pp.set_xlabel("$L$", fontsize=14)
ax_pp.set_ylabel("$C_L^{\\phi\\phi,\\mathrm{DDM}}/C_L^{\\phi\\phi,\\Lambda\\mathrm{CDM}} - 1$",
                 fontsize=14)
ax_pp.set_xlim(10, 800)

fig.tight_layout()
for ext in ("png", "pdf"):
    out = os.path.join(HERE, f"dcdm_wdm_fig3_9models.{ext}")
    fig.savefig(out, dpi=110)
    print("wrote", out)
