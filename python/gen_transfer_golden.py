"""One-shot generator for transfer-function golden tables (issue #309).

Run ONCE on the pre-refactor build:  python gen_transfer_golden.py
Do NOT re-run after the refactor starts -- the goldens are the oracle.
"""
import os
import numpy as np
from classy import Class

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "transfer_golden")

COMMON = {
    "output": "mTk,vTk",
    "gauge": "synchronous",
    "P_k_max_1/Mpc": 1.0,
    "z_pk": 0.0,
}

CASES = {
    "cdm_ur_ncdm": {
        "N_ur": 1.0, "N_ncdm": 2, "m_ncdm": "0.04,0.06", "T_ncdm": "0.71611,0.71611",
    },
    "fluid": {
        "Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.0, "N_ur": 3.044,
    },
    "fluid_nonppf": {
        "Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.1, "use_ppf": "no",
        "N_ur": 3.044,
    },
    "scf": {
        "Omega_Lambda": 0.0, "Omega_fld": 0.0,
        "Omega_scf": -1, "attractor_ic_scf": "yes",
        "scf_parameters": "10.0, 0.0, 0.0, 0.0", "N_ur": 3.044,
    },
    "dcdm_dr": {
        "Omega_ini_dcdm": 0.01, "Gamma_dcdm": 100.0, "N_ur": 3.044,
    },
    "dncdm_dr": {
        "N_ur": 2.0,
        "YHe": 0.25,
        "dncdm1.type": "ncdm_decay_dr",
        "dncdm1.m": 1.0,
        "dncdm1.T": 0.71611,
        "dncdm1.Gamma": 10.0,
        "dncdm1.Omega_ini": 0.001,
    },
    "idm_dr_idr": {
        "Omega_idm_dr": 0.12, "a_idm_dr": 1000.0, "nindex_idm_dr": 4,
        "xi_idr": 0.3, "N_ur": 2.0,
    },
    "idm_drmd_idr_drmd": {
        "z_stop": 1.0e4, "G_over_aH_drmd_ini": 1.0,
        "f_idm_drmd": 0.1, "delta_Neff_drmd": 0.5, "N_ur": 2.0,
    },
}


def flatten(tk):
    out = {}
    sample = next(iter(tk.values())) if tk else None
    # ic-nested branch: only reached when output includes isocurvature modes
    # (n_iso > 0). None of the CASES request them, but keep it so a future
    # iso case flattens correctly instead of crashing.
    if isinstance(sample, dict):
        for ic, d in tk.items():
            for name, arr in d.items():
                out[f"{ic}::{name}"] = np.asarray(arr, dtype=np.float64)
    else:
        for name, arr in tk.items():
            out[name] = np.asarray(arr, dtype=np.float64)
    return out


def transfer_for(extra):
    cosmo = Class()
    cosmo.set({**COMMON, **extra})
    try:
        cosmo.compute()
        return flatten(cosmo.get_transfer())
    finally:
        # Free the C++ structs even if compute() raises, so a crash in one
        # refactor commit doesn't leak state into later cases (repo convention).
        cosmo.struct_cleanup()
        cosmo.empty()


if __name__ == "__main__":
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    for name, extra in CASES.items():
        flat = transfer_for(extra)
        np.savez(
            os.path.join(GOLDEN_DIR, name + ".npz"),
            __order__=np.array(list(flat.keys())),
            **flat,
        )
        print(f"wrote {name}.npz  ({len(flat)} columns)")
