# %% [markdown]
# # Oscillations in the CMB from axion monodromy inflation
#
# Reproduces the numerical/analytic comparison of
#
# > R. Flauger, L. McAllister, E. Pajer, A. Westphal, G. Xu,
# > *Oscillations in the CMB from Axion Monodromy Inflation*,
# > [arXiv:0907.2916](https://arxiv.org/abs/0907.2916), JCAP 1006:009 (2010)
#
# using CLASS++'s inflation module (`P_k_ini type = inflation_V`), which solves the
# inflaton background and the Mukhanov--Sasaki equation mode by mode.
#
# The model is a linear (monodromy) potential with a periodic modulation from
# non-perturbative effects, the paper's eq. (2.2):
#
# $$V(\phi) = \mu^3\phi + \Lambda^4\cos(\phi/f) = \mu^3\left[\phi + b\,f\cos(\phi/f)\right],
#   \qquad b \equiv \frac{\Lambda^4}{\mu^3 f}.$$
#
# To first order in $b$, to leading order in $M_p/\phi_*$ and in slow roll, the paper
# derives (eqs. 2.4, 2.5 and 2.9)
#
# $$\Delta^2_{\mathcal R}(k) = \Delta^2_{\mathcal R}(k_*)\left(\frac{k}{k_*}\right)^{n_s-1}
#   \left[1 + \delta n_s\cos\!\left(\frac{\phi_k}{f}\right)\right],\qquad
#   \phi_k = \sqrt{\phi_*^2 - 2\ln(k/k_*)},$$
#
# $$\delta n_s = \frac{12\,b}{\sqrt{1+(3f\phi_*)^2}}
#   \sqrt{\frac{\pi}{8}\coth\!\left(\frac{\pi}{2f\phi_*}\right)f\phi_*}.$$
#
# **Figure 1 of the paper** plots that $\delta n_s(f)$ against a numerical solve at
# $b = 0.08$. That is what the cells below rebuild, plus the supporting plots and the
# convergence work needed to trust it.
#
# ## Two conventions to get right
#
# 1. **Units.** The paper uses reduced Planck units ($M_p^2 = 1/8\pi G$, $H^2 = \rho/3$);
#    the CLASS inflation module uses the non-reduced Planck mass ($H^2 = 8\pi\rho/3$, see
#    `primordial_inflation_derivs_member`). Hence $\phi_\mathrm{CLASS} =
#    \phi_\mathrm{paper}/\sqrt{8\pi}$ and $V_\mathrm{CLASS} = V_\mathrm{paper}/(8\pi)^2$.
#    The potential is form-invariant under that rescaling, so only $\mu^3$ and $f$ move;
#    $b$ and the shape are unchanged.
# 2. **Direction of roll.** The module requires $\phi$ to *grow* with $\mathrm{d}V/\mathrm{d}\phi<0$,
#    while the paper's inflaton rolls down towards zero. The `monodromy` shape therefore
#    works in $x = V_4 - \phi$, and since `inflation_V` pins $\phi_\mathrm{pivot} = 0$,
#    $V_4$ *is* $\phi_*$ (in CLASS units).
#
# ## Three precision parameters this model needs
#
# CLASS's defaults are tuned for featureless potentials and are *not* adequate here.
# Cell "Convergence" measures each of these; the rules are
#
# | parameter | default | needed | why |
# |---|---|---|---|
# | `primordial_inflation_ratio_min` | 100 | $\gtrsim 1.5/(2f\phi_*)$ | each mode must be started *before* it resonates with the oscillating background, which happens at $-k\tau = 1/(2f\phi_*)$ |
# | `k_per_decade_primordial` | 10 | $\gtrsim 3.7/(f\phi_*)$ | ten samples per oscillation of $P(k)$ |
# | `primordial_inflation_tol_integration` | $10^{-3}$ | $10^{-5}$ | the default leaves ~4% mode-to-mode scatter once the background oscillates |

# %%
import math
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from classy import Class

SQRT_8PI = math.sqrt(8.0 * math.pi)  # phi_paper = SQRT_8PI * phi_CLASS
K_PIVOT = 0.05  # 1/Mpc, CLASS default
PHI_STAR = 11.0  # M_p (reduced), ~60 e-folds before the end of inflation
A_S_TARGET = 2.1e-9

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "axion_monodromy_flauger2009_cache.npz")

plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "figure.autolayout": True})

FIGDIR = os.path.dirname(os.path.abspath(__file__))


def show(fig, name):
    """Display in the interactive window, and drop a PDF/PNG next to this file."""
    for suffix in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"axion_monodromy_{name}.{suffix}"))
    plt.show()


# %% [markdown]
# ## The analytic results of the paper

# %%
def delta_ns_analytic(f, b, phi_star=PHI_STAR):
    """arXiv:0907.2916 eq. (2.9). f, phi_star in REDUCED Planck units."""
    y = np.asarray(f, dtype=float) * phi_star
    return (12.0 * b / np.sqrt(1.0 + (3.0 * y) ** 2)
            * np.sqrt(np.pi / 8.0 / np.tanh(np.pi / (2.0 * y)) * y))


def phi_of_k(k, phi_star=PHI_STAR):
    """arXiv:0907.2916 eq. (2.5): the field value when k left the horizon."""
    return np.sqrt(phi_star**2 - 2.0 * np.log(k / K_PIVOT))


def resonance_depth(f, phi_star=PHI_STAR):
    """-k tau at which a mode resonates with the oscillating background (sec. 2.3).

    A mode has to be *inside* the horizon by this factor when the module starts
    integrating it, i.e. this is the floor on primordial_inflation_ratio_min.
    """
    return 1.0 / (2.0 * f * phi_star)


# %% [markdown]
# ## Running the CLASS inflation module

# %%
def monodromy_parameters(f, b, p=1.0, phi_star=PHI_STAR):
    """Translate (f, b, phi_*) in reduced Planck units into CLASS inputs.

    V_0 is fixed by the first-order slow-roll normalisation A_s = 128 pi/3 V^3/V'^2
    evaluated at the pivot; it only sets the overall amplitude, and delta n_s -- a
    *ratio* -- does not depend on it.
    """
    x_star = phi_star / SQRT_8PI
    return {
        "P_k_ini type": "inflation_V",
        "potential": "monodromy",
        "V_0": A_S_TARGET * 3.0 / (128.0 * math.pi * x_star**3),
        "V_1": f / SQRT_8PI,   # axion decay constant
        "V_2": b,              # modulation amplitude
        "V_3": p,              # monodromy power (1 = the paper's linear case)
        "V_4": x_star,         # x at the pivot scale
    }


def primordial_settings(f, phi_star=PHI_STAR, k_per_oscillation=10.0,
                        ratio_min_margin=3.0, tol=1e-5):
    """Precision parameters sized from the model, not from CLASS's defaults."""
    y = f * phi_star
    return {
        "k_per_decade_primordial": max(10.0, k_per_oscillation * math.log(10.0)
                                       / (2.0 * math.pi * y)),
        "primordial_inflation_ratio_min": max(100.0,
                                              ratio_min_margin * resonance_depth(f, phi_star)),
        "primordial_inflation_tol_integration": tol,
    }


def primordial_spectrum(f, b, extra=None):
    """P_s(k) from the inflation module. Returns (k, P_scalar).

    `output = tCl` and `modes = s,t` are not a choice: the module refuses to run
    without both scalars and tensors (input_module.cpp), even though we stop the
    pipeline at the primordial module and never touch a transfer function.
    """
    params = {"output": "tCl", "modes": "s, t",
              "l_max_scalars": 2500, "l_max_tensors": 500}
    params.update(monodromy_parameters(f, b))
    params.update(primordial_settings(f))
    if extra:
        params.update(extra)
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute(level=["primordial"])
    table = cosmo.get_primordial()
    k = table["k [1/Mpc]"].copy()
    P = table["P_scalar(k)"].copy()
    cosmo.struct_cleanup()
    cosmo.empty()
    return k, P


def modulation(f, b, extra=None):
    """P_s(k)/P_s^{b=0}(k) - 1, which to first order in b is delta n_s cos(phi_k/f).

    Dividing by the b=0 run of the *same* model is what makes the smooth part drop
    out exactly -- no detrending of the tilt or its running is needed.
    """
    k_smooth, P_smooth = primordial_spectrum(f, 0.0, extra=extra)
    k, P = primordial_spectrum(f, b, extra=extra)
    assert np.allclose(k, k_smooth)
    return k, P / P_smooth - 1.0


# %% [markdown]
# ## Measuring $\delta n_s$ from a computed spectrum
#
# The template is $\delta n_s\cos(\phi_k/f + \beta)$ with $\phi_k$ from eq. (2.5). The
# phase $\beta$ is free (the `cos`/`sin` columns), and so is a quadratic in $\ln k$ that
# soaks up any residual smooth mismatch.
#
# $\phi_*$ is **scanned rather than fixed at 11**. Eq. (2.5) is a slow-roll relation
# between $k$ and the field, so it is only good to $O(\epsilon)\sim0.4\%$; at $f=10^{-3}$
# that is $0.04/f \approx 40$ radians of accumulated phase across the fitting window, and
# a fit at fixed $\phi_*=11$ silently returns almost zero. What the scan optimises over is
# the *chirp* of the template, not its overall phase, which the free $\beta$ already covers.
#
# Maximising an amplitude over a nuisance parameter biases it upwards, so every call
# also reports (a) the residual left after subtracting the best fit and (b) the amplitude
# the same estimator returns for a deliberately wrong frequency. Both have to be small
# compared with the answer, or the answer means nothing.

# %%
def fit_delta_ns(k, ratio, f, window=2.0,
                 phi_scan=(10.5, 12.0), n_scan=1501, frequency_scale=1.0):
    """Project the modulation onto cos(phi_k/f), scanning phi_* for the best chirp.

    Returns (amplitude, phi_star_fit, residual_rms).
    """
    ln_k = np.log(k / K_PIVOT)
    inside = np.abs(ln_k) < window
    ln_k, y = ln_k[inside], ratio[inside]
    best = (0.0, np.nan, np.nan)
    for phi_star_fit in np.linspace(*phi_scan, n_scan):
        argument = phi_star_fit**2 - 2.0 * ln_k
        if np.any(argument <= 0.0):
            continue
        theta = np.sqrt(argument) / (f * frequency_scale)
        design = np.column_stack([np.cos(theta), np.sin(theta),
                                  np.ones_like(ln_k), ln_k, ln_k**2])
        coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
        amplitude = math.hypot(coefficients[0], coefficients[1])
        if amplitude > best[0]:
            best = (amplitude, phi_star_fit, (y - design @ coefficients).std())
    return best


def measure(f, b=0.08, extra=None, window=2.0):
    """delta n_s for one (f, b), with its two controls."""
    k, ratio = modulation(f, b, extra=extra)
    amplitude, phi_fit, residual = fit_delta_ns(k, ratio, f, window=window)
    # Control: the same estimator at 1/3 of the true frequency must find nothing.
    null, _, _ = fit_delta_ns(k, ratio, f, window=window, frequency_scale=1.0 / 3.0)
    return dict(f=f, b=b, delta_ns=amplitude, phi_star_fit=phi_fit,
                residual=residual, null=null,
                analytic=float(delta_ns_analytic(f, b)))


# %% [markdown]
# ## The spectrum itself
#
# Before any fitting: does the computed $P_s(k)$ actually look like eq. (2.4)?

# %%
F_SHOW, B_SHOW = 0.006, 0.08
k_show, ratio_show = modulation(F_SHOW, B_SHOW)
fit_show = fit_delta_ns(k_show, ratio_show, F_SHOW)
k_smooth_show, P_smooth_show = primordial_spectrum(F_SHOW, 0.0)
_, P_show = primordial_spectrum(F_SHOW, B_SHOW)

fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.2), sharex=True,
                         gridspec_kw={"height_ratios": [1.0, 1.0]})
band = (k_show > 3e-4) & (k_show < 5e-2)

axes[0].plot(k_show[band], 1e9 * P_show[band], lw=0.8, color="C0",
             label=fr"CLASS++, $b={B_SHOW}$, $f={F_SHOW}\,M_p$")
axes[0].plot(k_show[band], 1e9 * P_smooth_show[band], lw=1.4, color="0.4",
             label="same model with $b=0$")
axes[0].set_xscale("log")
axes[0].set_ylabel(r"$10^{9}\,\Delta^2_{\mathcal{R}}(k)$")
axes[0].legend(frameon=False, loc="lower left")
axes[0].set_title("Axion monodromy: the modulation in the primordial spectrum")

theta_show = phi_of_k(k_show[band], fit_show[1]) / F_SHOW
axes[1].plot(k_show[band], ratio_show[band], lw=0.8, color="C0",
             label=r"$P_s^{b}/P_s^{b=0}-1$ (CLASS++)")
# The best-fit phase, recovered the same way fit_delta_ns does it.
ln_k_band = np.log(k_show[band] / K_PIVOT)
design_show = np.column_stack([np.cos(theta_show), np.sin(theta_show),
                               np.ones_like(ln_k_band), ln_k_band, ln_k_band**2])
coef_show, *_ = np.linalg.lstsq(design_show, ratio_show[band], rcond=None)
axes[1].plot(k_show[band], design_show @ coef_show, lw=1.0, ls="--", color="C3",
             label=r"$\delta n_s\cos(\phi_k/f+\beta)$, eq. (2.4)")
axes[1].axhline(0.0, lw=0.5, color="0.7")
for sign in (+1, -1):
    axes[1].axhline(sign * fit_show[0], lw=0.7, ls=":", color="C3")
axes[1].set_xscale("log")
axes[1].set_xlabel(r"$k\ [\mathrm{Mpc}^{-1}]$")
axes[1].set_ylabel(r"$\Delta P_s/P_s$")
axes[1].set_ylim(-0.27, 0.21)
axes[1].legend(frameon=False, loc="lower center", ncol=2, fontsize=8)
show(fig, "spectrum")

print(f"f = {F_SHOW}  b = {B_SHOW}")
print(f"  delta n_s  numerical = {fit_show[0]:.4f}")
print(f"             analytic  = {float(delta_ns_analytic(F_SHOW, B_SHOW)):.4f}")
print(f"  phi_* recovered = {fit_show[1]:.3f}  (input 11.0; eq. 2.5 is slow-roll)")
print(f"  residual after subtracting the template = {fit_show[2]:.5f}")

# %% [markdown]
# ## Figure 1: $\delta n_s$ against $f$ at $b = 0.08$
#
# The paper's figure runs over $0 < f < 0.1\,M_p$ on a linear axis, and its numerical
# points fall off the analytic curve and plunge towards zero below $f\approx 2\times10^{-3}$.
# The right-hand panel is the same data on a log axis, extended down to $f = 2\times10^{-4}$
# — a decade below where the paper's numerics gave up — to show that the breakdown was
# never in the physics.

# %%
F_GRID = np.array([0.1, 0.08, 0.06, 0.05, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015,
                   0.012, 0.01, 0.008, 0.006, 0.005, 0.004, 0.003, 0.002, 0.0015,
                   0.001, 0.0007, 0.0005, 0.0003, 0.0002])
B_FIG1 = 0.08


def figure1_scan(f_grid=F_GRID, b=B_FIG1, use_cache=True):
    if use_cache and os.path.exists(CACHE):
        with np.load(CACHE) as data:
            if data["f"].shape == f_grid.shape and np.allclose(data["f"], f_grid):
                print(f"loaded {CACHE}")
                return {key: data[key] for key in data.files}
    rows = []
    for f in f_grid:
        t0 = time.time()
        # Fewer than ~5 oscillations in the window and the phi_* scan has nothing to
        # bite on; there, fix phi_* and use the whole computed range instead.
        oscillations = 2.0 * 2.0 / (2.0 * math.pi * f * PHI_STAR)
        if oscillations > 5.0:
            row = measure(f, b)
        else:
            k, ratio = modulation(f, b)
            amplitude, phi_fit, residual = fit_delta_ns(
                k, ratio, f, window=1e9, phi_scan=(PHI_STAR, PHI_STAR), n_scan=1)
            null, _, _ = fit_delta_ns(k, ratio, f, window=1e9,
                                      phi_scan=(PHI_STAR, PHI_STAR), n_scan=1,
                                      frequency_scale=1.0 / 3.0)
            row = dict(f=f, b=b, delta_ns=amplitude, phi_star_fit=phi_fit,
                       residual=residual, null=null,
                       analytic=float(delta_ns_analytic(f, b)))
        rows.append(row)
        print(f"f={f:8.5f}  dns={row['delta_ns']:.4f}  analytic={row['analytic']:.4f}  "
              f"ratio={row['delta_ns']/row['analytic']:.3f}  "
              f"resid={row['residual']:.4f}  null={row['null']:.4f}  "
              f"({time.time()-t0:5.1f} s)")
    out = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    np.savez(CACHE, **out)
    return out


scan = figure1_scan()

# %%
f_dense = np.logspace(math.log10(1.5e-4), math.log10(0.11), 600)
curve = delta_ns_analytic(f_dense, B_FIG1)

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))

axes[0].plot(f_dense, curve, color="C3", lw=1.2, label="analytic, eq. (2.9)")
axes[0].plot(scan["f"], scan["delta_ns"], "o", ms=3.2, color="C0",
             label="CLASS++ inflation module")
axes[0].set_xlim(0.0, 0.1)
axes[0].set_ylim(0.0, 0.27)
axes[0].set_xlabel(r"$f\ [M_p]$")
axes[0].set_ylabel(r"$\delta n_s$")
axes[0].set_title(f"Flauger et al. fig. 1  ($b={B_FIG1}$)")
axes[0].legend(frameon=False, loc="lower right")

axes[1].plot(f_dense, curve, color="C3", lw=1.2)
axes[1].plot(scan["f"], scan["delta_ns"], "o", ms=3.2, color="C0")
axes[1].plot(scan["f"], scan["null"], "x", ms=3.0, color="0.6",
             label="control: wrong template frequency")
axes[1].axvspan(1.5e-4, 2e-3, color="C1", alpha=0.10)
axes[1].text(5.5e-4, 0.225, "the paper's numerics\nbreak down here",
             fontsize=7, color="C1", ha="center")
axes[1].set_xscale("log")
axes[1].set_xlabel(r"$f\ [M_p]$")
axes[1].set_ylabel(r"$\delta n_s$")
axes[1].set_title("same data, log axis")
axes[1].set_ylim(-0.01, 0.27)
axes[1].legend(frameon=False, loc="center left", fontsize=7)
show(fig, "fig1")

ratio = scan["delta_ns"] / scan["analytic"]
print(f"numerical/analytic over the whole range: "
      f"min {ratio.min():.3f}, median {np.median(ratio):.3f}, max {ratio.max():.3f}")
print(f"controls: residual/amplitude <= {(scan['residual']/scan['delta_ns']).max():.3f}, "
      f"wrong-frequency/amplitude <= {(scan['null']/scan['delta_ns']).max():.3f}")

# %% [markdown]
# ## Why the paper's numerics broke down at small $f$
#
# The resonance between a mode and the oscillating background happens at
# $-k\tau = 1/(2f\phi_*)$, i.e. that many horizon radii *inside* the horizon
# (paper, sec. 2.3). The inflation module starts each mode at
# $-k\tau = $ `primordial_inflation_ratio_min`, whose default is **100**. Once
# $f\phi_* < 1/200$, that is *after* the resonance — the code integrates a mode that has
# already missed the only event that imprints the oscillation, and returns a spectrum with
# no feature in it.
#
# This is a precision-parameter floor, not a limitation of the method: the plot below
# collapses onto one curve in `ratio_min / resonance_depth`, and everything above
# $\approx1.5$ is converged.

# %%
DEPTH_F = [0.002, 0.0005, 0.0002]
DEPTH_RATIO_MIN = [100.0, 200.0, 300.0, 1000.0, 3000.0]

fig, ax = plt.subplots(figsize=(5.4, 3.4))
for colour, f in zip(("C0", "C1", "C2"), DEPTH_F):
    xs, ys = [], []
    for ratio_min in DEPTH_RATIO_MIN:
        row = measure(f, B_FIG1,
                      extra={"primordial_inflation_ratio_min": ratio_min})
        xs.append(ratio_min / resonance_depth(f))
        ys.append(row["delta_ns"] / row["analytic"])
        print(f"f={f:7.4f}  ratio_min={ratio_min:6.0f}  "
              f"ratio_min/depth={xs[-1]:6.2f}  dns/analytic={ys[-1]:.3f}")
    ax.plot(xs, ys, "o-", ms=3.5, lw=1.0, color=colour,
            label=fr"$f={f}$, resonance at $-k\tau={resonance_depth(f):.0f}$")
ax.axhline(1.0, lw=0.7, color="0.6")
ax.axvline(1.0, lw=0.7, ls=":", color="0.4")
ax.set_xscale("log")
ax.set_xlabel(r"$\mathtt{primordial\_inflation\_ratio\_min} \times 2f\phi_*$")
ax.set_ylabel(r"$\delta n_s^{\rm num}/\delta n_s^{\rm analytic}$")
ax.set_title("A mode must be started before it resonates")
ax.legend(frameon=False, fontsize=7)
show(fig, "resonance_depth")

# %% [markdown]
# ## Convergence
#
# Every precision knob that could plausibly matter, at the hardest $f$ actually used.
# The estimator is insensitive to the integration tolerance because white mode-to-mode
# scatter does not project onto the template — but the *spectrum* is not, and that
# matters for anything downstream (a $C_\ell$, a likelihood).

# %%
F_CONV = 0.0005
CONVERGENCE_CASES = [
    ("reference", {}),
    ("k_per_decade_primordial x2", {"k_per_decade_primordial":
                                    2 * primordial_settings(F_CONV)["k_per_decade_primordial"]}),
    ("ratio_min x3", {"primordial_inflation_ratio_min":
                      3 * primordial_settings(F_CONV)["primordial_inflation_ratio_min"]}),
    ("tol_integration 1e-3 (default)", {"primordial_inflation_tol_integration": 1e-3}),
    ("tol_integration 1e-7", {"primordial_inflation_tol_integration": 1e-7}),
    ("pt_stepsize 0.002", {"primordial_inflation_pt_stepsize": 0.002}),
    ("bg_stepsize 0.001", {"primordial_inflation_bg_stepsize": 0.001}),
    ("ratio_max 0.005", {"primordial_inflation_ratio_max": 0.005}),
]

print(f"f = {F_CONV}, b = {B_FIG1};  analytic delta n_s = "
      f"{float(delta_ns_analytic(F_CONV, B_FIG1)):.5f}")
print(f"{'case':34s} {'delta n_s':>10s} {'/analytic':>10s} {'residual':>10s} {'time':>7s}")
for label, extra in CONVERGENCE_CASES:
    t0 = time.time()
    row = measure(F_CONV, B_FIG1, extra=extra)
    print(f"{label:34s} {row['delta_ns']:10.5f} "
          f"{row['delta_ns']/row['analytic']:10.4f} {row['residual']:10.5f} "
          f"{time.time()-t0:6.1f}s")

# %% [markdown]
# ## Where the last 1--2% comes from
#
# The numerics sit consistently $1$--$2\%$ **above** eq. (2.9), which is exactly what a
# derivation carried to first order in $b$ and to leading order in slow roll should do.
# Splitting the excess by varying $b$ at fixed $f$ separates the two:
#
# * a piece $\propto b$ — the $O(b^2)$ term the derivation drops;
# * a $b$-independent floor of $\approx 0.85\%$ — the slow-roll corrections dropped when
#   the Hankel index was set to $\nu_0 = 3/2$ ($\epsilon = 1/2\phi_*^2 = 0.41\%$ here).
#
# If the excess did not vanish as $b\to0$ there would be something wrong with either the
# module or this notebook.

# %%
F_BSCAN = 0.006
B_SCAN = [0.16, 0.08, 0.04, 0.02, 0.01]
print(f"f = {F_BSCAN}")
print(f"{'b':>7s} {'delta n_s':>11s} {'analytic':>11s} {'ratio':>8s}")
b_values, ratios = [], []
for b in B_SCAN:
    row = measure(F_BSCAN, b)
    b_values.append(b)
    ratios.append(row["delta_ns"] / row["analytic"])
    print(f"{b:7.3f} {row['delta_ns']:11.5f} {row['analytic']:11.5f} {ratios[-1]:8.4f}")

slope, intercept = np.polyfit(b_values, ratios, 1)
print(f"\nfit: numerical/analytic = {intercept:.4f} + {slope:.3f} * b")
print(f"  -> O(b) correction   {slope:.3f} b   ({100*slope*0.08:.1f}% at b=0.08)")
print(f"  -> b-independent floor {100*(intercept-1):.2f}%  "
      f"(slow roll: epsilon = {1/(2*PHI_STAR**2)*100:.2f}%)")

# %% [markdown]
# ## What it does to the CMB
#
# The paper's fig. 5 shows $C_\ell^{TT}$ for its best fit to WMAP5. Without redoing the
# fit, what is worth seeing is how the primordial modulation transfers through the full
# pipeline. It is periodic in $\ln k$ and $\ell\propto k$, so it reaches the $C_\ell$ as
# rapid wiggles at low $\ell$ and as slow, wide features at high $\ell$ -- and it survives
# projection and lensing at the ~10% level, which is why WMAP could constrain it at all.
#
# The printed check re-runs both models with the primordial $k$-sampling doubled and the
# source sampling halved: an oscillating $P(k)$ is exactly the case where the spline onto
# the perturbation $k$-grid could alias, and a $C_\ell$ nobody checked for that is not
# worth plotting.

# %%
F_CL, B_CL = 0.01, 0.08


def lensed_tt(f, b, extra=None):
    params = {"output": "tCl,pCl,lCl", "lensing": "yes", "modes": "s, t",
              "l_max_scalars": 2500, "l_max_tensors": 500}
    params.update(monodromy_parameters(f, b))
    params.update(primordial_settings(f))
    if extra:
        params.update(extra)
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    cls = cosmo.lensed_cl(2500)
    cosmo.struct_cleanup()
    cosmo.empty()
    return cls


cl_modulated = lensed_tt(F_CL, B_CL)
cl_smooth = lensed_tt(F_CL, 0.0)
ell = cl_modulated["ell"][2:]
delta_cl = cl_modulated["tt"][2:] / cl_smooth["tt"][2:] - 1.0

fig, ax = plt.subplots(figsize=(6.4, 3.2))
ax.plot(ell, 100 * delta_cl, lw=0.8, color="C0")
ax.axhline(0.0, lw=0.5, color="0.7")
ax.set_xlim(2, 2500)
ax.set_xlabel(r"$\ell$")
ax.set_ylabel(r"$\Delta C_\ell^{TT}/C_\ell^{TT}$  [%]")
ax.set_title(fr"Lensed TT, $f={F_CL}\,M_p$, $b={B_CL}$ "
             fr"($\delta n_s={float(delta_ns_analytic(F_CL, B_CL)):.2f}$)")
show(fig, "cl_tt")

print(f"max |Delta C_l/C_l| = {100*np.abs(delta_cl).max():.2f}% at "
      f"l = {int(ell[np.argmax(np.abs(delta_cl))])}")
print(f"rms over 2 <= l <= 2500 = {100*delta_cl.std():.2f}%")

finer = {"k_per_decade_primordial": 2*primordial_settings(F_CL)["k_per_decade_primordial"],
         "perturb_sampling_stepsize": 0.04}
delta_cl_finer = (lensed_tt(F_CL, B_CL, extra=finer)["tt"][2:]
                  / lensed_tt(F_CL, 0.0, extra=finer)["tt"][2:] - 1.0)
print(f"sampling check: max |change| in Delta C_l/C_l when the primordial k-grid is "
      f"doubled and perturb_sampling_stepsize halved = "
      f"{100*np.abs(delta_cl_finer - delta_cl).max():.3f}% (signal is "
      f"{100*np.abs(delta_cl).max():.2f}%)")
