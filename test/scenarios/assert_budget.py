#!/usr/bin/env python3
"""Run ./class on a scenario with background_verbose=2 and assert the printed
budget equation. Usage:
    assert_budget.py <ini> [--tol 0.02] [--positive-line "Neutrino Species"] [--expect-error SUBSTR]
- default: assert the "TOTAL ... Omega = X" line has |X-1| <= tol.
- --positive-line NAME: also assert a "-> NAME ... Omega = Y" line exists with Y > 0.
- --expect-error SUBSTR: instead expect ./class to FAIL with SUBSTR in its output.
Exits 0 on success, 1 on failure (prints reason + captured output tail).
"""
import argparse, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLASS = os.path.join(ROOT, "class")

def run_class(ini_text):
    with tempfile.TemporaryDirectory() as d:
        ini = os.path.join(d, "run.ini")
        with open(ini, "w") as f:
            f.write(ini_text)
        p = subprocess.run([CLASS, ini], capture_output=True, text=True)
        return p.returncode, p.stdout + p.stderr

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ini")
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--positive-line", default=None)
    ap.add_argument("--expect-error", default=None)
    a = ap.parse_args()

    with open(a.ini) as f:
        base = [ln for ln in f.read().splitlines()
                if not re.match(r"\s*(background_verbose|root)\s*=", ln)]
    with tempfile.TemporaryDirectory() as d:
        ini_text = "\n".join(base) + f"\nbackground_verbose = 2\nroot = {d}/out_\n"
        rc, out = run_class(ini_text)

    if a.expect_error is not None:
        if rc == 0 or a.expect_error not in out:
            print(f"FAIL: expected error containing {a.expect_error!r}; rc={rc}\n{out[-1500:]}")
            return 1
        print(f"OK: rejected with {a.expect_error!r}")
        return 0

    if rc != 0:
        print(f"FAIL: ./class exited {rc}\n{out[-1500:]}")
        return 1
    # Shooting prints one budget block per iteration; assert on the LAST (converged) block.
    blocks = out.split("Budget equation")
    budget = blocks[-1] if len(blocks) > 1 else out
    m = re.search(r"TOTAL\s+Omega = ([-0-9.eE+]+)", budget)
    if not m:
        print(f"FAIL: no TOTAL line in budget output\n{out[-1500:]}")
        return 1
    total = float(m.group(1))
    if abs(total - 1.0) > a.tol:
        print(f"FAIL: TOTAL Omega = {total}, |1-total| > tol {a.tol}\n{out[-1500:]}")
        return 1
    if a.positive_line:
        pm = re.search(rf"->\s+{re.escape(a.positive_line)}\s+Omega = ([-0-9.eE+]+)", budget)
        if not pm or float(pm.group(1)) <= 0.0:
            print(f"FAIL: line {a.positive_line!r} missing or <= 0\n{out[-1500:]}")
            return 1
    print(f"OK: TOTAL Omega = {total}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
