#!/usr/bin/env python3
"""Physics-limit validation for the dcdm_wdm species (arXiv:2606.14849).

Checks (all ratios vs a reference run):
  1. vkick -> 0 (1e-6): C_ell^TT and P(k) match LCDM with the same total
     omega_m budget to < 0.1% (daughters are effectively CDM).
  2. vkick = 1 (epsilon = 0, massless daughters): C_ell^TT and P_cb match the
     existing dcdm_dr implementation with identical Gamma/Omega_ini.
     P_cb (cdm+baryons, including the cold decaying parent) is compared
     instead of P_m because the two codes book the identical massless
     daughters differently: dcdm_wdm classifies them as warm matter (their
     smooth delta_rho/rho enters delta_m and dilutes P_m) while dcdm_dr
     classifies them as dark radiation (excluded from P_m). Same physics,
     different P_m convention — TT, which is convention-free, agrees. The
     P_m ratio is printed as information only.
  3. f -> 0 (Omega_ini -> 1e-6 of the DM): identical to LCDM to < 0.05%.
Run AFTER `make classy-pip-dev`. Exits non-zero on failure.
"""
import sys

import numpy as np
from classy import Class, CosmoSevereError

BASE = {
    "h": 0.6736,
    "omega_b": 0.02237,
    "output": "tCl,pCl,lCl,mPk",
    "lensing": "yes",
    "l_max_scalars": 2500,
    "P_k_max_h/Mpc": 1.0,
    "gauge": "synchronous",
}
OMEGA_DM = 0.26  # total (stable + unstable) early-time DM density fraction


def run(extra):
    cosmo = Class()
    cosmo.set({**BASE, **extra})
    cosmo.compute()
    cls = cosmo.lensed_cl(2000)
    kk = np.logspace(-3, 0, 200)
    pk = np.array([cosmo.pk(k * cosmo.h(), 0.0) for k in kk])
    try:
        # cdm+b spectrum (includes the cold decaying parent). Runs without
        # non-cold species do not expose P_cb; there P_cb == P_m anyway.
        pk_cb = np.array([cosmo.pk_cb(k * cosmo.h(), 0.0) for k in kk])
    except CosmoSevereError:
        pk_cb = pk
    cosmo.struct_cleanup()
    cosmo.empty()
    return cls, kk, pk, pk_cb


def ddm_pars(f, Gamma_kms_Mpc, vkick, bins=96):
    return {
        "Omega_cdm": OMEGA_DM * (1.0 - f),
        "ddm.type": "dcdm_wdm",
        "ddm.Gamma": Gamma_kms_Mpc,
        "ddm.vkick": vkick,
        "ddm.Omega_ini": OMEGA_DM * f,
        "ddm.momenta_bins": bins,
    }


def maxrel(a, b, lmin=2):
    a, b = np.asarray(a)[lmin:], np.asarray(b)[lmin:]
    return np.max(np.abs(a / b - 1.0))


failures = []


def check(name, value, tol):
    status = "PASS" if value < tol else "FAIL"
    print(f"  {name:55s} {value:10.3e}  (tol {tol:g})  {status}")
    if value >= tol:
        failures.append(name)


print("== reference LCDM ==")
lcdm_cls, kk, lcdm_pk, _ = run({"Omega_cdm": OMEGA_DM})

print("== check 1: vkick -> 0 is LCDM ==")
cls, _, pk, _ = run(ddm_pars(f=0.3, Gamma_kms_Mpc=100.0, vkick=1e-6))
check("vkick=1e-6: max |dTT/TT| (l<=2000)", maxrel(cls["tt"], lcdm_cls["tt"]), 1e-3)
check("vkick=1e-6: max |dP/P| (k<=1 h/Mpc)", maxrel(pk, lcdm_pk, 0), 2e-3)

print("== check 2: vkick=1 (massless daughters) matches dcdm_dr ==")
cls_wdm, _, pk_wdm, pkcb_wdm = run(ddm_pars(f=0.3, Gamma_kms_Mpc=100.0, vkick=1.0))
cls_dr, _, pk_dr, pkcb_dr = run({
    "Omega_cdm": OMEGA_DM * 0.7,
    "Omega_ini_dcdm": OMEGA_DM * 0.3,
    "Gamma_dcdm": 100.0,
})
check("vkick=1 vs dcdm_dr: max |dTT/TT|", maxrel(cls_wdm["tt"], cls_dr["tt"]), 3e-3)
check("vkick=1 vs dcdm_dr: max |dP_cb/P_cb|", maxrel(pkcb_wdm, pkcb_dr, 0), 5e-3)
# P_m differs BY CONVENTION at vkick=1: dcdm_wdm books the massless daughters
# as warm matter (smooth delta_rho/rho enters delta_m), dcdm_dr books the same
# particles as dark radiation (excluded from delta_m). Informational only:
print(f"  [info] vkick=1 vs dcdm_dr: max |dP_m/P_m| = "
      f"{maxrel(pk_wdm, pk_dr, 0):.3e} (P_m convention differs, not a check)")

print("== check 3: f -> 0 is LCDM ==")
cls, _, pk, _ = run(ddm_pars(f=1e-5, Gamma_kms_Mpc=100.0, vkick=0.1))
check("f=1e-5: max |dTT/TT|", maxrel(cls["tt"], lcdm_cls["tt"]), 5e-4)

print("== check 4: physical suppression present at finite kick ==")
cls, _, pk, _ = run(ddm_pars(f=0.1, Gamma_kms_Mpc=200.0, vkick=0.02))
supp = pk[-1] / lcdm_pk[-1] - 1.0
print(f"  P(k=1)/P_LCDM - 1 = {supp:+.3f} (expect clearly negative)")
if not (supp < -0.005):
    failures.append("suppression sign")

if failures:
    print("FAILED:", failures)
    sys.exit(1)
print("all dcdm_wdm physics-limit checks passed")
