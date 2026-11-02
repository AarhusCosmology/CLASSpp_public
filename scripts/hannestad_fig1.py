#!/usr/bin/env python3
"""Hannestad & Madsen 1997 (astro-ph/9702125), fig 1: two-level decay H -> F + phi.

Runs the standalone decay harness (build target `decay_harness`, tools/decay_harness.h)
once per parent lifetime, overlays the final light-neutrino spectrum q^2 f_F on the
paper's digitized curves (scripts/data/hannestad_madsen_1997/fig1.csv), and prints the
rms and max |Delta(q^2 f)| per lifetime. Issue #389, step 2.

    python3 scripts/hannestad_fig1.py [--m_H 2] [--r_l 0.1] [--r_phi 0.0025]
                                      [--caption | --tau ... [--paper 0.2 1 20]]
                                      [--out hannestad_fig1]
                                      [--harness build/cmake/decay_harness]

CONVENTIONS (issue #389, the maintainer's comment of 2026-11-01; each was checked
against the digitized data):
  * The paper's x-axis is labelled p/T_gamma but is q = p/T_nu, the harness's comoving
    momentum: int q^2 f_F dq = 3.55 for all three curves, against 3.553 on a p/T_nu
    axis and 1.29 on a post-annihilation p/T_gamma axis. y = q^2 f_F, f_F the bare
    per-dof occupation. The harness prints p/T_gamma, so x is converted with
    q = (p/T_gamma) T_gamma / T_nu, from the block's own header.
  * Every neutrino starts Fermi-Dirac in ENERGY at T_start = 5 MeV, F included; the
    harness does exactly that with T_start=5.
  * Three flavours, H and F two of them, so one spectator flavour: n_eff_spectator=1.
  * m_H = 2, m_F = 0.2, m_phi = 0.005 MeV, i.e. r_l = 0.1, r_phi = 0.0025.

Writes <out>_tau<tau>.dat (the harness's own table) for every lifetime and, when
matplotlib is installed, <out>.png. Further harness settings pass through unchanged
with --extra key=value, e.g. --extra daughter_n=192 for a resolution study.

THE LIFETIMES. The caption says tau_H = 0.2, 1, 20 s, but the curves are the paper's
own lifetime parameter alpha = (m_H/3 MeV)^2 (tau/1 s) = 0.1, 1, 10, i.e. tau = 0.225,
2.25, 22.5 s at m_H = 2 MeV: those reproduce all three curves to <= 0.004, while the
caption's "1 s" misses the dotted curve by 0.11 (PR #492 has the measurements and the
controls). So by default the script runs alpha = 0.1, 1, 10 against the curves the
caption labels 0.2, 1, 20 s; --caption runs the caption's lifetimes instead. --tau runs
any lifetimes, and --paper names the digitized curve each is compared with (default:
the same lifetime).

The comparison is made at the digitized q in [0.1, 13]. Our curve is interpolated there
linearly in ln f against q (f is exponential in q on the tail, where a linear
interpolation of q^2 f between log-spaced nodes would itself be ~1e-2 off).
"""

import argparse
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "hannestad_madsen_1997", "fig1.csv")
COLUMN = {0.2: "F_tau0.2s", 1.0: "F_tau1s", 20.0: "F_tau20s"}
STYLE = {0.2: "-", 1.0: ":", 20.0: "--"}


def read_final_block(text):
    """The last block of the harness output: (header dict, daughter rows)."""
    blocks = [b for b in text.strip().split("\n\n\n") if b.strip()]
    header = {}
    rows = []
    for line in blocks[-1].splitlines():
        if line.startswith("#"):
            if "=" in line:
                key, value = line[1:].split("=", 1)
                header[key.strip()] = float(value.split()[0])
        elif not line.startswith("parent"):
            rows.append([float(x) for x in line.split()])
    return header, rows


def read_paper(path):
    """{column: [(q, y)]} from the digitized table; empty cells are skipped."""
    if not os.path.exists(path):
        sys.exit(f"digitized data not found: {path}")
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(line for line in fh if not line.startswith("#"))
        for row in reader:
            q = float(row["q"])
            for key, value in row.items():
                if key != "q" and value.strip():
                    out.setdefault(key, []).append((q, float(value)))
    if not out:
        sys.exit(f"no curves in {path}")
    return out


def interp_q2f(q_nodes, f_nodes, q):
    """q^2 f at q, interpolating ln f linearly in q between the bracketing nodes."""
    if q <= q_nodes[0] or q >= q_nodes[-1]:
        raise ValueError(f"q = {q} outside the harness grid")
    lo, hi = 0, len(q_nodes) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if q_nodes[mid] <= q:
            lo = mid
        else:
            hi = mid
    f0, f1 = f_nodes[lo], f_nodes[hi]
    t = (q - q_nodes[lo]) / (q_nodes[hi] - q_nodes[lo])
    if f0 > 0.0 and f1 > 0.0:
        f = math.exp((1.0 - t) * math.log(f0) + t * math.log(f1))
    else:
        f = (1.0 - t) * f0 + t * f1
    return q * q * f


def compare(q_nodes, f_nodes, paper):
    """rms and max |Delta(q^2 f)| at the digitized points with q in [0.1, 13]."""
    diffs = []
    for q, y in paper:
        if 0.1 <= q <= 13.0:
            diffs.append(interp_q2f(q_nodes, f_nodes, q) - y)
    rms = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    return rms, max(abs(d) for d in diffs), len(diffs)


def peak(xs, ys):
    """Height and position of the maximum, refined by a parabola through 3 points."""
    i = max(range(len(ys)), key=lambda k: ys[k])
    if 0 < i < len(ys) - 1:
        x0, x1, x2 = xs[i - 1], xs[i], xs[i + 1]
        y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
        den = (x0 - x1) * (x0 - x2) * (x1 - x2)
        a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
        b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den
        if a < 0.0:
            xv = -b / (2.0 * a)
            c = y0 - a * x0 * x0 - b * x0
            return a * xv * xv + b * xv + c, xv
    return ys[i], xs[i]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--m_H", type=float, default=2.0, help="parent mass, MeV (fig 1: 2)")
    ap.add_argument("--r_l", type=float, default=0.1, help="m_F / m_H (fig 1: 0.1)")
    ap.add_argument(
        "--r_phi", type=float, default=0.0025, help="m_phi / m_H (fig 1: 0.0025)"
    )
    ap.add_argument(
        "--tau",
        type=float,
        nargs="+",
        default=None,
        help="lifetimes, s (default: alpha = 0.1, 1, 10, i.e. (3 MeV/m_H)^2 x those)",
    )
    ap.add_argument(
        "--caption",
        action="store_true",
        help="run the caption's lifetimes, 0.2, 1, 20 s, instead of alpha = 0.1, 1, 10",
    )
    ap.add_argument(
        "--paper",
        type=float,
        nargs="+",
        default=None,
        help="the paper's curve (0.2, 1 or 20) to compare each --tau run against; "
        "default: the same lifetime. E.g. --tau 0.225 2.25 22.5 --paper 0.2 1 20 "
        "tests alpha = (m_H/3 MeV)^2 (tau/1 s) = 0.1, 1, 10",
    )
    ap.add_argument("--out", default="hannestad_fig1")
    ap.add_argument("--harness", default="build/cmake/decay_harness")
    ap.add_argument("--data", default=DATA, help="digitized fig 1 table")
    ap.add_argument(
        "--extra", nargs="*", default=[], help="further key=value harness arguments"
    )
    args = ap.parse_args()

    paper = read_paper(args.data)
    if args.tau is None:
        if args.caption:
            args.tau = [0.2, 1.0, 20.0]
        else:
            args.tau = [alpha * (3.0 / args.m_H) ** 2 for alpha in (0.1, 1.0, 10.0)]
            if args.paper is None:
                args.paper = [0.2, 1.0, 20.0]
    elif args.caption:
        sys.exit("--caption and --tau are exclusive")
    against = args.paper if args.paper is not None else args.tau
    if len(against) != len(args.tau):
        sys.exit("--paper needs one entry per --tau")
    curves = []
    for tau in args.tau:
        cmd = [
            args.harness,
            f"m_H={args.m_H}",
            f"tau={tau}",
            f"r_l={args.r_l}",
            f"r_phi={args.r_phi}",
            "T_start=5",
            "n_eff_spectator=1",
        ] + args.extra
        run = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if run.returncode != 0:
            sys.exit(f"{' '.join(cmd)} failed:\n{run.stderr}")
        sys.stderr.write(f"tau = {tau:g} s: {run.stderr.strip()}\n")
        with open(f"{args.out}_tau{tau:g}.dat", "w") as fh:
            fh.write(run.stdout)
        header, rows = read_final_block(run.stdout)
        to_q = header["T_gamma"] / header["T_nu"]
        q = [r[0] * to_q for r in rows]
        f = [r[1] for r in rows]
        curves.append((tau, q, f, against[len(curves)]))

    print("tau [s]   rms|dy|   max|dy|   n   peak ours (q)      peak paper (q)")
    for tau, q, f, ref in curves:
        y = [qi * qi * fi for qi, fi in zip(q, f)]
        h, x = peak(q, y)
        col = COLUMN.get(ref)
        if col is None or col not in paper:
            print(f"{tau:7g}   (no digitized curve)                 {h:.3f} ({x:.2f})")
            continue
        rms, mx, n = compare(q, f, paper[col])
        ph, px = peak([p[0] for p in paper[col]], [p[1] for p in paper[col]])
        print(
            f"{tau:7g}   {rms:.4f}    {mx:.4f}   {n}  {h:.3f} ({x:.2f})       "
            f"{ph:.3f} ({px:.2f})"
        )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr.write("matplotlib not installed: wrote the .dat files only\n")
        return

    fig, ax = plt.subplots(figsize=(6, 4.2))
    for tau, q, f, ref in curves:
        style = STYLE.get(ref, "-")
        (line,) = ax.plot(
            q, [qi * qi * fi for qi, fi in zip(q, f)], style, label=rf"$\tau_H = {tau:g}$ s"
        )
        col = COLUMN.get(ref)
        if col in paper:
            ax.plot(
                [p[0] for p in paper[col]],
                [p[1] for p in paper[col]],
                "o",
                ms=2.5,
                color=line.get_color(),
                label=f"HM97 {ref:g} s",
            )
    ax.set_xlim(0.0, 13.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$q = p/T_\nu$  (labelled $p/T_\gamma$ in the paper)")
    ax.set_ylabel(r"$q^2 f_F$")
    ax.set_title(
        rf"$m_H = {args.m_H:g}$ MeV, $m_F/m_H = {args.r_l:g}$, "
        rf"$m_\phi/m_H = {args.r_phi:g}$"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=150)
    sys.stderr.write(f"wrote {args.out}.png\n")


if __name__ == "__main__":
    main()
