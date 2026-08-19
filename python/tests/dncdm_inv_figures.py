#!/usr/bin/env python3
"""dncdm_inv_figures.py -- regenerate the arXiv:2011.01502 BACKGROUND figures.

NOT a `make test` / CTest target (same status as dncdm_inv_ab.py): this is a
documented driver you run by hand. It exists so the branch's inverse-decay sector
can be checked against a *published* result rather than only against itself, and
so that check is repeatable by someone who was not present when it was first done.

Reproduces, from `inverse_decays=yes` background runs:

  --fig2   Figure 2: background energy densities of nu_H / nu_l / phi and their
           sum vs. scale factor, for g = 3e-13, 1e-13, 1e-14 at m_nuH = 0.3 eV
           (left column), plus the final (a=1) phase-space distributions of nu_l
           and phi (right column).
  --fig3   Figure 3: background phase-space distributions q^3 f(q) of nu_H, phi,
           nu_l at a = 0.0005 and a = 0.001, for g = 3e-13, m_nuH = 0.3 eV. This
           is the quasi-equilibrium figure; the paper overlays thermal fits with
           hand-tuned (T, mu) which are reproduced here from the caption.
  --fig11  Figure 11(a): the exact transition-network identity
           n_nuH + 1/2 (n_phi + n_nul) = const, which is what the conservative
           formulation exists to deliver.

PARAMETER DICTIONARY (verified against the reference implementation's own printed
`decay_rate_eq` to six significant figures):

    Gamma_input[km/s/Mpc] = g^2 * m[eV] * 3.730568e33

from the paper's eq. (2.5) Majorana rest-frame rate Gamma_0 = g^2 m / (4 pi),
combined with our input convention `Gamma_ = Gamma_raw * (1e3/_c_)`. So the three
published couplings map to Gamma = 1.007e8, 1.119e7, 1.119e5.

INITIAL CONDITION. The published runs use the reference's `ncdm_psd_parameters`
amplitude vector A = [phi, nu_H, nu_l] = [0, 1, 1]: phi starts EMPTY while nu_l
starts at a full Fermi-Dirac. Here that is `dncdm1.dr_f_ini_phi=0` and
`dncdm1.dr_f_ini_l=1`. (Before those per-daughter keys existed both daughters
shared one seeding flag and this configuration was unreachable, so figures 2/3/11
could not be reproduced at all.)

GRID. `dncdm1.dr_q_sampling=linear`, `dr_q_max=22`, `dr_N_q=240` reproduces the
published grid EXACTLY -- the reference notebook builds it as
`np.linspace(0, q_max, N+1)[1:]`, which is quadrature.c's `qm_trapz` qmin==0
branch. The parent is put on the same grid via `quadrature_strategy=3`,
`momenta_bins=240`, `max_q=22`.

KNOWN GAP -- the "dec" (decay-only) curves. The paper's figures 2 and 11 show a
three-way ladder: dec / dec+inv / dec+inv+qs. The first is decay-only *within the
momentum-resolved framework*, i.e. the reference's `has_decay_sector=1` with
`has_inverse_bg=0`. This branch cannot express that: `inverse_decays=no` does not
switch the kernel's inverse term off, it declines to build the momentum-resolved
composite at all (species/dncdm_dr_species.cpp's CreateAll dispatch) and falls
back to the integrated-DR `DNCDM_DR_Species` path -- a different discretisation,
not the same calculation with one term removed. The kernel itself already carries
a `Config::inverse_decays` flag, so the missing piece is only a way to build the
composite with it false; that is an interface change and is deliberately NOT made
here. Consequences:
  * dec+inv and dec+inv+qs are produced from this branch.
  * dec is produced ONLY when --reference-bin points at a build of the older
    Tram-Oldengott code (ThomasTram/class @ decaying_neutrinos), which has the
    four-way background switch ladder (has_decay_sector / has_collision /
    has_inverse_bg / has_qs_bg) and can therefore do it.
Do not read a missing dotted line as agreement; it means the curve was not run.

OPEN QUESTION the figures surface. In the g=1e-14 panel the nu_H comoving number
does not go to zero: it freezes at ~2e-5 of its initial value, so a^4 rho_nuH turns
around and climbs again. This is the high-momentum tail, whose rest-frame rate is
time-dilated by (m/E) and which therefore never decays; at a=1 the survivors have
median q = 17.5 with 90% of the population above q = 10, so the residual sits in
the last few bins of a q_max = 22 grid and its size is grid-sensitive. The published
panel shows nu_H leaving the frame instead. Whether that is a genuine disagreement
or only a plotting range has not been settled -- run with --reference-bin and
compare before treating either as physics.

REFERENCE OVERLAY. `--reference-bin PATH` additionally runs the older code at the
matched coupling and overlays it dashed. That code is a peer research code, not an
oracle: at low Gamma (late, non-relativistic decay) its daughters are produced at
comoving q ~ a m /(2 T_0) far beyond its `qmax=22` heuristic and its own
conservation identity then drifts by O(1). At the couplings used here decay is
early and relativistic, the daughters stay well inside the grid, and it is sound.

USAGE

    python python/tests/dncdm_inv_figures.py --all
    python python/tests/dncdm_inv_figures.py --fig2 --outdir /tmp/figs
    python python/tests/dncdm_inv_figures.py --all \
        --reference-bin ~/src/class_decaying_neutrinos/class

Requires numpy + matplotlib. Runs are cached in --workdir by config, so re-running
only replots (the g=3e-13 point costs ~45 s the first time).
"""

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from dncdm_inv_ab import parse_dat, repo_root  # noqa: E402

# ── the dictionary ──────────────────────────────────────────────────────────
GAMMA_PER_G2_M = 3.730568e33  # Gamma_input[km/s/Mpc] = g^2 * m[eV] * this
M_NU_H = 0.3                  # eV, the published mass
COUPLINGS = [3e-13, 1e-13, 1e-14]

# Paper fig. 3 caption: hand-tuned thermal fits at the two epochs.
# (T/T_0, mu) per species, keyed by scale factor.
#
# CAVEAT on the nu_H overlay. The caption defines the fit with the species' own
# energy, q^3 f_Th = q^3/(exp[(eps_i(q) - mu_i)/T_i] +- 1), and that is what is
# implemented here. Taken literally it puts the nu_H dotted curve at a = 0.001 near
# 0.6 at its peak, whereas the published panel shows ~1.2; evaluating the same fit
# with eps = q instead of sqrt(q^2 + (aM)^2) gives ~1.09, i.e. the published curve.
# At a = 0.001, aM ~ 1.78 is not small next to a thermal q ~ 3, so the two differ by
# e^{0.49/0.8} ~ 1.85 -- the whole discrepancy. This affects only the ANALYTIC
# OVERLAY: the computed nu_H distribution itself matches the published solid curve
# (peak ~0.95 at q ~ 3 in both). Left caption-literal rather than tuned to match.
THERMAL_FITS = {
    0.0005: {"nu_H": (0.88, 0.15), "phi": (0.80, -0.33), "nu_l": (0.90, 0.51)},
    0.0010: {"nu_H": (0.80, 0.46), "phi": (0.84, -0.20), "nu_l": (0.85, 0.70)},
}

Q_MAX, N_Q = 22.0, 240

BASE = {
    "h": "0.67556",
    "T_cmb": "2.7255",
    "omega_b": "0.022032",
    "omega_cdm": "0.12038",
    "N_ur": "1",
    "YHe": "0.245",  # daughters push DeltaNeff past the BBN table; bypass it
    "dncdm1.type": "ncdm_decay_dr",
    "dncdm1.m": str(M_NU_H),
    "dncdm1.T": "0.71611",
    "dncdm1.deg": "1",
    "dncdm1.inverse_decays": "yes",
    "dncdm1.quadrature_strategy": "3",   # qm_trapz -- same family as the daughters
    "dncdm1.momenta_bins": str(N_Q),
    "dncdm1.max_q": str(Q_MAX),
    "dncdm1.dr_q_sampling": "linear",    # the published grid
    "dncdm1.dr_q_max": str(Q_MAX),
    "dncdm1.dr_N_q": str(N_Q),
    "dncdm1.dr_bg_refine": "1",
    "dncdm1.dr_f_ini_l": "1",            # A = [0, 1, 1]
    "dncdm1.dr_f_ini_phi": "0",
    "evolver": "2",                      # rkdp45; ndf15's Jacobian is dense here
    "write background": "yes",
    "background_verbose": "1",
}


def gamma_of(g, m=M_NU_H):
    """Published Yukawa coupling -> our Gamma input (km/s/Mpc)."""
    return g * g * m * GAMMA_PER_G2_M


def q_grid():
    """The published linear grid: q_i = (i+1) q_max/N, last cell half-weighted."""
    h = Q_MAX / N_Q
    q = np.array([(i + 1) * h for i in range(N_Q)])
    q[-1] = Q_MAX
    dq = np.full(N_Q, h)
    dq[-1] *= 0.5
    return q, dq


# ── running ─────────────────────────────────────────────────────────────────
def _run(binary, params, workdir, tag, timeout=1800):
    """Run `binary` on an ini built from `params`; cached on the param hash."""
    key = hashlib.sha1(
        (str(sorted(params.items())) + str(binary)).encode()
    ).hexdigest()[:12]
    root = workdir / f"{tag}_{key}_"
    out = Path(str(root) + "background.dat")
    if out.exists():
        return out
    ini = workdir / f"{tag}_{key}.ini"
    ini.write_text(
        "\n".join(f"{k} = {v}" for k, v in params.items()) + f"\nroot = {root}\n"
    )
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run([str(binary), str(ini)], capture_output=True, text=True,
                          env=env, timeout=timeout, cwd=str(workdir))
    if not out.exists():
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-6:]
        raise RuntimeError(f"{tag}: no background.dat.\n  " + "\n  ".join(tail))
    return out


def run_branch(binary, g, qs, workdir):
    """This branch, inverse decays ON, quantum statistics per `qs`."""
    p = dict(BASE)
    p["dncdm1.Gamma"] = f"{gamma_of(g):.6e}"
    p["dncdm1.quantum_statistics"] = "yes" if qs else "no"
    return _run(binary, p, workdir, f"br_g{g:.0e}_qs{int(qs)}")


REF_BASE = {
    "h": "0.67556", "T_cmb": "2.7255", "omega_b": "0.022032",
    "omega_cdm": "0.12038", "N_ur": "1", "N_ncdm": "3",
    "m_ncdm": f"1e-10, {M_NU_H}, 1e-10",
    "has_decay_sector": "1", "has_collision": "0",
    "ncdm_psd_parameters": "0.0, 1, 1",            # A = [0, 1, 1]
    "Quadrature strategy": "3, 3, 3",              # qm_trapz
    "Number of momentum bins": f"{N_Q},{N_Q},{N_Q}",
    "Maximum q": f"{int(Q_MAX)},{int(Q_MAX)},{int(Q_MAX)}",
    "background_method": "1", "evolver": "0", "evolver_bg": "0",
    "write background": "yes", "background_verbose": "1",
}


def run_reference(binary, g, inverse, qs, workdir):
    """Older Tram-Oldengott code; it *can* do the decay-only rung of the ladder."""
    p = dict(REF_BASE)
    p["coupling_constant"] = f"{g:.6e}"
    p["has_inverse_bg"] = "1" if inverse else "0"
    p["has_qs_bg"] = "1" if qs else "0"
    return _run(binary, p, workdir, f"ref_g{g:.0e}_i{int(inverse)}_qs{int(qs)}")


# ── reading ─────────────────────────────────────────────────────────────────
def cols(path):
    names, data = parse_dat(path)
    return {n: data[:, i] for i, n in enumerate(names)}


def branch_series(path):
    c = cols(path)
    a = 1.0 / (1.0 + c["z"])
    o = np.argsort(a)
    return {
        "a": a[o],
        "nu_H": c["(.)rho_dncdm1"][o],
        "nu_l": c["(.)rho_dr_dncdm1_l"][o],
        "phi": c["(.)rho_dr_dncdm1_phi"][o],
        "n_nu_H": c["(.)number_dncdm1"][o],
        "n_nu_l": c["(.)number_dr_dncdm1_l"][o],
        "n_phi": c["(.)number_dr_dncdm1_phi"][o],
        "f_nu_l": np.vstack([c[f"f_dr_dncdm1_l[{i}]"] for i in range(N_Q)]).T[o],
        "f_phi": np.vstack([c[f"f_dr_dncdm1_phi[{i}]"] for i in range(N_Q)]).T[o],
        # parent stores f suppressed by kappa = deg_H (2 pi)^3 / 2 (issue #385);
        # multiply back out to get the BARE per-dof occupation the paper plots.
        "f_nu_H": (
            (1.0 * 8.0 * np.pi ** 3 / 2.0)
            * np.exp(np.vstack([c[f"lnf_dncdm1[{i}]"] for i in range(N_Q)]).T)
        )[o],
    }


def reference_series(path):
    c = cols(path)
    a = 1.0 / (1.0 + c["z"])
    o = np.argsort(a)
    out = {"a": a[o],
           "nu_H": c["(.)rho_ncdm[1]"][o],
           "nu_l": c["(.)rho_ncdm[2]"][o],
           "phi": c["(.)rho_ncdm[0]"][o],
           "n_nu_H": c["(.)number_ncdm[1]"][o],
           "n_nu_l": c["(.)number_ncdm[2]"][o],
           "n_phi": c["(.)number_ncdm[0]"][o]}
    for key, tag in (("f_phi", "phi"), ("f_nu_H", "nu1"), ("f_nu_l", "nu2")):
        k0 = f"f_{tag}__0"
        if k0 in c:
            out[key] = np.vstack([c[f"f_{tag}__{i}"] for i in range(N_Q)]).T[o]
    return out


def dimensionless_mass():
    """M = m/(k_B T_ncdm T_cmb): the parent's mass on the dimensionless q axis."""
    k_B, eV = 1.3806504e-23, 1.602176487e-19
    return M_NU_H / (k_B * 0.71611 * 2.7255 / eV)


def thermal_f(q, T, mu, fermion, a=None):
    """f_Th = 1/(exp((eps - mu)/T) +- 1), paper fig. 3 caption.

    eps is the species' OWN energy: eps = q for the massless daughters, but
    eps = sqrt(q^2 + (a M)^2) for the massive parent -- at a = 5e-4, aM ~ 0.89 is
    comparable to the thermal q, so using eps = q there would misplace the curve.
    """
    eps = q if a is None else np.sqrt(q ** 2 + (a * dimensionless_mass()) ** 2)
    x = np.clip((eps - mu) / T, -600.0, 600.0)
    if fermion:
        return 1.0 / (np.exp(x) + 1.0)
    # Bose: regular only for x > 0; the fits have mu < 0 for phi so this holds.
    return 1.0 / np.expm1(np.maximum(x, 1e-12))


# ── figures ─────────────────────────────────────────────────────────────────
COLOR = {"nu_H": "tab:blue", "nu_l": "tab:orange", "phi": "tab:green",
         "sum": "tab:red"}
SYMBOL = {"nu_H": r"$\nu_H$", "nu_l": r"$\nu_l$", "phi": r"$\phi$",
          "sum": r"$\nu_H+\nu_l+\phi$"}


def figure2(branch_bin, ref_bin, workdir, outdir):
    import matplotlib.pyplot as plt
    q, _ = q_grid()
    fig, axes = plt.subplots(len(COUPLINGS), 2, figsize=(11, 3.1 * len(COUPLINGS)))
    for row, g in enumerate(COUPLINGS):
        axL, axR = axes[row]
        curves = []  # (series, linestyle, label)
        if ref_bin:
            curves.append((reference_series(run_reference(ref_bin, g, False, False,
                                                          workdir)), ":", "dec (ref)"))
        curves.append((branch_series(run_branch(branch_bin, g, False, workdir)),
                       "--", "dec+inv"))
        curves.append((branch_series(run_branch(branch_bin, g, True, workdir)),
                       "-", "dec+inv+qs"))

        plateau = 0.0
        for s, ls, lab in curves:
            a = s["a"]
            tot = np.zeros_like(a)
            # NOTE the nu_H curve does not vanish: after the plunge its comoving
            # number freezes at ~2e-5 of its initial value and a^4 rho then RISES
            # (frozen comoving number, gone non-relativistic => rho ~ a^-3). That is
            # NOT the kFFloor seed (which would sit ~23 decades lower) and it is not
            # masked here, because it is a real feature of the solution: the
            # survivors are the highest-momentum particles, whose decay rate
            # Gamma_0 (m/E) is the most time-dilated. At g=1e-14, a=1 the residual
            # population has median q = 17.5 with 90% of it above q = 10 -- i.e. it
            # lives in the last few bins and its size is therefore sensitive to
            # dr_q_max. Worth checking against --reference-bin before it is read as
            # physics.
            for sp in ("nu_H", "nu_l", "phi"):
                axL.loglog(a, a ** 4 * s[sp], ls, color=COLOR[sp], lw=1.3,
                           label=SYMBOL[sp] if ls == "-" else None)
                tot = tot + s[sp]
            axL.loglog(a, a ** 4 * tot, ls, color=COLOR["sum"], lw=1.3,
                       label=SYMBOL["sum"] if ls == "-" else None)
            plateau = max(plateau, np.max(a ** 4 * tot))
        # Published axes: a from 1e-5, and ~5 decades of a^4 rho below the plateau.
        # Both cuts are load-bearing rather than cosmetic. Below a ~ 1e-5 nothing has
        # happened yet and the daughters sit at the kFFloor positivity seed, which in
        # a^4 rho reaches ~1e-113 and would otherwise set the y-scale. At the top end
        # the parent underflows to denormals once it has fully decayed (its rho runs
        # off to ~1e-80 non-monotonically), so its curve legitimately leaves the frame
        # rather than being drawn as noise.
        axL.set_xlim(1e-5, 1.3)
        axL.set_ylim(plateau * 1e-5, plateau * 3)
        axL.set_xlabel("$a$")
        axL.set_ylabel(r"$a^4\bar\rho$  [arb.]")
        axL.set_title(rf"$\mathfrak{{g}} = {g:g}$,  $\Gamma = {gamma_of(g):.3g}$",
                      fontsize=9)
        axL.legend(fontsize=7, loc="lower right", ncol=2)
        axL.grid(alpha=0.25, which="both", lw=0.4)

        # right: final (a=1) phase-space distributions of the two daughters
        for s, ls, lab in curves:
            if "f_nu_l" not in s:
                continue
            k = int(np.argmin(np.abs(s["a"] - 1.0)))
            axR.plot(q, q ** 3 * s["f_nu_l"][k], ls, color=COLOR["nu_l"], lw=1.3,
                     label=rf"$\nu_l$ {lab}")
            axR.plot(q, q ** 3 * s["f_phi"][k], ls, color=COLOR["phi"], lw=1.3,
                     label=rf"$\phi$ {lab}")
        axR.set_xlabel("$q$  [$T_0$]")
        axR.set_ylabel(r"$q^3\bar f(q)$")
        axR.set_title("final phase-space distributions", fontsize=9)
        axR.set_xlim(0, Q_MAX)
        axR.set_ylim(bottom=0)
        axR.legend(fontsize=6, loc="upper right")
        axR.grid(alpha=0.25, lw=0.4)
    fig.suptitle(r"arXiv:2011.01502 figure 2  —  $m_{\nu_H}=0.3$ eV, massless $\nu_l$",
                 y=1.002)
    fig.tight_layout()
    p = outdir / "figure2.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def figure3(branch_bin, ref_bin, workdir, outdir):
    import matplotlib.pyplot as plt
    q, _ = q_grid()
    g = 3e-13  # the coupling the published figure uses
    s = branch_series(run_branch(branch_bin, g, True, workdir))
    epochs = sorted(THERMAL_FITS)
    # (name, label, is_fermion, is_massive). nu_H is a MASSIVE FERMION: both flags
    # matter -- Bose statistics would diverge as q->0 and eps=q would misplace the
    # curve at a=5e-4 where a*M ~ 0.89 is comparable to the thermal momentum.
    species = [("nu_H", r"$\nu_H$", True, True), ("phi", r"$\phi$", False, False),
               ("nu_l", r"$\nu_l$", True, False)]
    fig, axes = plt.subplots(len(epochs), 3, figsize=(12, 3.3 * len(epochs)))
    for row, a_t in enumerate(epochs):
        k = int(np.argmin(np.abs(s["a"] - a_t)))
        for col, (sp, lab, fermion, massive) in enumerate(species):
            ax = axes[row, col]
            ax.plot(q, q ** 3 * s[f"f_{sp}"][k], "-", color=COLOR[sp], lw=1.4,
                    label="CLASSpp")
            T, mu = THERMAL_FITS[a_t][sp]
            ax.plot(q, q ** 3 * thermal_f(q, T, mu, fermion,
                                          a=s["a"][k] if massive else None),
                    ":", color="k", lw=1.1, label=rf"thermal $T={T},\ \mu={mu}$")
            ax.set_title(rf"{lab},  $a={s['a'][k]:.3g}$", fontsize=9)
            ax.set_xlabel("$q$  [$T_0$]")
            ax.set_xlim(0, Q_MAX)
            ax.set_ylim(bottom=0)
            if col == 0:
                ax.set_ylabel(r"$q^3\bar f(q)$")
            ax.legend(fontsize=6)
            ax.grid(alpha=0.25, lw=0.4)
    fig.suptitle(r"arXiv:2011.01502 figure 3  —  quasi-equilibrium, "
                 rf"$\mathfrak{{g}}=3\times10^{{-13}}$, $m_{{\nu_H}}=0.3$ eV", y=1.002)
    fig.tight_layout()
    p = outdir / "figure3.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def figure11(branch_bin, ref_bin, workdir, outdir):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5))
    report = []
    for g in COUPLINGS:
        entries = []
        if ref_bin:
            entries.append((reference_series(run_reference(ref_bin, g, True, True,
                                                           workdir)), ":", "ref dec+inv+qs"))
        entries.append((branch_series(run_branch(branch_bin, g, True, workdir)),
                        "-", "dec+inv+qs"))
        for s, ls, lab in entries:
            a = s["a"]
            C = (s["n_nu_H"] + 0.5 * (s["n_phi"] + s["n_nu_l"])) * a ** 3
            m = C > 0
            rel = np.abs(C[m] / C[m][0] - 1.0)
            # The first few output rows carry a start-up transient of the ODE
            # solver that is identical across couplings (same ICs, same grid) and
            # sits at a ~ 1e-13, decades before any decay. Reporting the raw max
            # therefore measures the initial step, not the transition network.
            # Quote the END-TO-END drift as the headline and the early max
            # separately so neither is mistaken for the other.
            final = rel[-1]
            early = np.max(rel)
            ax.loglog(a[m], rel + 1e-18, ls, lw=1.3,
                      label=rf"$\mathfrak{{g}}={g:g}$ {lab} "
                            rf"(a=1: {final:.1e}; peak {early:.1e})")
            report.append((g, lab, final, early, a[m][int(np.argmax(rel))]))
    ax.set_xlabel("$a$")
    ax.set_ylabel(r"$|\,C(a)/C(a_{\rm ini}) - 1\,|$")
    ax.set_title(r"arXiv:2011.01502 figure 11(a):  "
                 r"$C = n_{\nu_H} + \frac{1}{2}(n_\phi + n_{\nu_l})$ comoving")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25, which="both", lw=0.4)
    fig.tight_layout()
    p = outdir / "figure11a.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("\n  conservation identity C = n_nuH + 1/2(n_phi + n_nul), comoving:")
    print(f"    {'coupling':<12s} {'variant':<17s} {'drift at a=1':<14s} "
          f"{'peak':<11s} {'peak at a':<11s}")
    for g, lab, final, early, a_peak in report:
        print(f"    g={g:<10g} {lab:<17s} {final:<14.3e} {early:<11.3e} {a_peak:<11.3e}")
    print("    (the peak sits at the solver's first steps, decades before decay;\n"
          "     the drift at a=1 is the transition network's actual conservation)")
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fig2", action="store_true")
    ap.add_argument("--fig3", action="store_true")
    ap.add_argument("--fig11", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--outdir", default=None, help="where to write the PDFs")
    ap.add_argument("--workdir", default=None, help="run cache (default: <outdir>/runs)")
    ap.add_argument("--class-bin", default=None, help="branch `class` (default: repo ./class)")
    ap.add_argument("--reference-bin", default=None,
                    help="older ThomasTram/class@decaying_neutrinos build; enables the "
                         "decay-only ('dec') rung this branch cannot express")
    args = ap.parse_args()
    if not (args.fig2 or args.fig3 or args.fig11 or args.all):
        args.all = True

    root = Path(repo_root())
    branch_bin = Path(args.class_bin) if args.class_bin else root / "class"
    if not branch_bin.exists():
        sys.exit(f"no CLASS binary at {branch_bin} -- run `make class` first")
    ref_bin = Path(args.reference_bin) if args.reference_bin else None
    if ref_bin and not ref_bin.exists():
        sys.exit(f"--reference-bin {ref_bin} does not exist")
    if not ref_bin:
        print("note: no --reference-bin, so the decay-only ('dec') curves of figures 2\n"
              "      and 11 are OMITTED -- this branch cannot express that rung (see the\n"
              "      module docstring). Absent curves are not agreement.")

    outdir = Path(args.outdir) if args.outdir else root / "generated" / "dncdm_inv_figures"
    outdir.mkdir(parents=True, exist_ok=True)
    workdir = Path(args.workdir) if args.workdir else outdir / "runs"
    workdir.mkdir(parents=True, exist_ok=True)

    made = []
    if args.all or args.fig2:
        made.append(figure2(branch_bin, ref_bin, workdir, outdir))
    if args.all or args.fig3:
        made.append(figure3(branch_bin, ref_bin, workdir, outdir))
    if args.all or args.fig11:
        made.append(figure11(branch_bin, ref_bin, workdir, outdir))
    print("\n  wrote:")
    for p in made:
        print(f"    {p}")


if __name__ == "__main__":
    main()
