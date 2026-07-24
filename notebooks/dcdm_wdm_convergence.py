"""Momentum-grid convergence sweep for dcdm_wdm: run once per classy module.

Usage:
    python dcdm_wdm_convergence.py <module> <tag> [out.npz]

<module> is the python module to drive ("classy" = this branch's
injection-adapted grid, "classyref" = master's uniform-ln-q grid). Each
(model x bins) run executes in a SUBPROCESS with a hard timeout, and the
output .npz is rewritten after every run — a hang or crash costs one point,
never the sweep. Identical inputs on both modules (same evolver=2, same
k / l grids): an apples-to-apples accuracy-vs-bins test.
"""
import os
import subprocess
import sys
import time

import numpy as np

MODULE = sys.argv[1] if len(sys.argv) > 1 else "classy"
TAG = sys.argv[2] if len(sys.argv) > 2 else ("new" if MODULE == "classy" else "old")
OUT = sys.argv[3] if len(sys.argv) > 3 else f"conv_{TAG}.npz"
HERE = os.path.dirname(os.path.abspath(__file__))
TIMEOUT = 2400  # s per run; the 192-bin reference is the slowest

KK = np.logspace(-3, 0, 300)  # h/Mpc
ELL = np.arange(2, 2001)
MODELS = {
    "g0p1_v0p03": (0.1, 0.03),   # Gamma^-1 = 0.1 Gyr, v = 0.03  (early, hard)
    "g1_v0p011":  (1.0, 0.011),  # Gamma^-1 = 1.0 Gyr, v = 0.011 (notebook default)
}
BINS = [16, 24, 32, 48, 96, 192]  # 192 = reference ("truth")

RUN_SRC = r'''
import importlib, sys, time
import numpy as np
module, ig, vk, bins, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
cls = importlib.import_module(module).Class
GYR = 977.792
c = cls()
c.set({"h": 0.6736, "omega_b": 0.02237, "output": "tCl,pCl,lCl,mPk", "lensing": "yes",
       "l_max_scalars": 2500, "P_k_max_h/Mpc": 1.0, "gauge": "synchronous", "evolver": "2",
       "Omega_cdm": 0.26 * 0.9, "ddm.type": "dcdm_wdm", "ddm.Gamma": GYR / ig,
       "ddm.vkick": vk, "ddm.Omega_ini": 0.26 * 0.1, "ddm.momenta_bins": bins})
t0 = time.time()
c.compute()
t = time.time() - t0
h = c.h()
kk = np.logspace(-3, 0, 300)
pk = np.array([c.pk(k * h, 0.0) for k in kk])
cl = np.array(c.lensed_cl(2000)["tt"][2:2001])
np.savez(out, pk=pk, cl=cl, t=t)
print(f"OK {t:.1f}")
'''


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] [{TAG}] {msg}", flush=True)


data = {"KK": KK, "ELL": ELL, "BINS": np.array(BINS), "TAG": np.array(TAG)}
out_path = os.path.join(HERE, OUT)
for name, (ig, vk) in MODELS.items():
    for b in BINS:
        shard = os.path.join(HERE, f"_shard_{TAG}_{name}_{b}.npz")
        t0 = time.time()
        r = None  # cleared every iteration: never let a stale prior-run's
        # stdout/stderr get logged under this run's FAILED message.
        try:
            r = subprocess.run([sys.executable, "-c", RUN_SRC, MODULE,
                                str(ig), str(vk), str(b), shard],
                               capture_output=True, text=True, timeout=TIMEOUT)
            ok = r.returncode == 0 and os.path.exists(shard)
        except subprocess.TimeoutExpired:
            ok = False
            log(f"{name} bins={b}: TIMEOUT after {TIMEOUT}s")
        if ok:
            s = np.load(shard)
            data[f"{name}_b{b}_pk"] = s["pk"]
            data[f"{name}_b{b}_cl"] = s["cl"]
            data[f"{name}_b{b}_t"] = float(s["t"])
            os.remove(shard)
            log(f"{name} bins={b:4d}: {data[f'{name}_b{b}_t']:6.1f} s compute, "
                f"{time.time() - t0:6.1f} s wall")
        else:
            data[f"{name}_b{b}_pk"] = np.full_like(KK, np.nan)
            data[f"{name}_b{b}_cl"] = np.full(len(ELL), np.nan)
            data[f"{name}_b{b}_t"] = np.nan
            if r is not None:
                log(f"{name} bins={b}: FAILED -- {(r.stdout + r.stderr)[-300:]}")
            if os.path.exists(shard):  # partial write (e.g. killed mid-savez)
                os.remove(shard)
        np.savez(out_path, **data)  # incremental: saved after EVERY run
log(f"saved -> {out_path}")
