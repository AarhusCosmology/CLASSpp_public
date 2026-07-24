"""Old (uniform ln-q, master) vs new (injection-adapted, this branch) momentum-grid
convergence, reframed on the FIXED erf deposit (commit 1586d250 completed the
dJ/dlnq measure term; the earlier 311337bc numbers here were stale, generated
before that fix).

Truth = new-192 (this branch's grid at 192 bins). This is no longer a
"does old agree with new" consistency gate -- it is independently validated
two other ways (see .superpowers/sdd/task-3-diagnosis.md and task-3-report.md
Sec 5):
  (1) background rho_wdm AND p_wdm match analytic monochromatic-kick
      predictions from the branch's own history (rho_dev, p_dev <~ 2e-3,
      most <1e-4);
  (2) P_m(k) matches a LambdaCDM control at k<=0.01 h/Mpc to <0.1%
      (task-3-report.md Sec.5: R-1 = -0.037%/-0.074% for the two fiducial
      models, -0.015% for the vkick->0 CDM limit).

Old-192 vs new-192 is now a MEASUREMENT of the old scheme's known bias, not a
gate: master's D-renormalized Gaussian deposit carries a flat +1.9% <beta^2>
velocity excess relative to its own analytic w(z) (task-3-diagnosis.md Sec.2,
beta^2_ratio pinned at 1.916-1.927 across z and models) -- the same q^3-tilt
placement-bias mechanism documented in the design spec's v2 history. That
excess velocity shifts the free-streaming cutoff in P(k); it is invisible to
the energy budget (rho_wdm agrees to <~1e-4) and to the P(k) itself at
k<=0.01 (below the cutoff), which is why the low-k sanity check below is
expected to land near zero even though the two schemes disagree substantially
at high k.
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

HERE = os.path.dirname(os.path.abspath(__file__))
new = np.load(os.path.join(HERE, "conv_new.npz"))
old = np.load(os.path.join(HERE, "conv_old.npz"))

KK = new["KK"]
BINS = list(new["BINS"])
REF = 192
SWEEP = [b for b in BINS if b != REF]
MODELS = {"g0p1_v0p03": r"$\Gamma^{-1}=0.1$ Gyr, $v=0.03$ (early)",
          "g1_v0p011":  r"$\Gamma^{-1}=1$ Gyr, $v=0.011$"}

# Low-k anchor: below the free-streaming cutoff for both fiducial models,
# where the old scheme's velocity bias should be invisible in P(k).
K_LOWK = 0.01
I_LOWK = int(np.argmin(np.abs(KK - K_LOWK)))
LOWK_SANITY_TOL = 0.01  # 1%: if the old-scheme bias at k<=0.01 exceeds this,
# the "bias is a high-k-only free-streaming effect" premise is contradicted --
# BLOCK rather than silently report a reframed story that doesn't hold up.

# --- dataviz: documented default palette, slots assigned in fixed order ---
# 2-series nominal categorical -> slot 1 (blue) then slot 2 (orange); both
# checks (CVD, normal-vision, contrast) validated for this pair (see
# scripts/validate_palette.js "#2a78d6,#eb6834" --mode light -> ALL PASS).
C_OLD, C_NEW = "#2a78d6", "#eb6834"          # slot 1 (old), slot 2 (new)
INK = "#0b0b0b"                              # primary text
INK_SECONDARY = "#52514e"                    # secondary text
INK_MUTED = "#898781"                        # axis / tick text
GRID = "#e1e0d9"                             # hairline gridline
TARGET_LINE = "#c3c2b7"                      # baseline/axis gray, for target refs
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "text.color": INK,
    "axes.edgecolor": TARGET_LINE,
    "axes.labelcolor": INK_SECONDARY,
    "xtick.color": INK_MUTED,
    "ytick.color": INK_MUTED,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def max_err(d, name, b, ref_pk, ref_cl):
    pk = d[f"{name}_b{b}_pk"]
    cl = d[f"{name}_b{b}_cl"]
    return (np.nanmax(np.abs(pk / ref_pk - 1.0)),
            np.nanmax(np.abs(cl / ref_cl - 1.0)))


def bins_for_target(bins, errs, target):
    """Smallest bin count whose error <= target (log-log interp between points).

    Assumes error is (roughly) monotonically decreasing in bin count -- true
    for the new grid's own self-convergence, but the old grid's error floors
    at its ~1.9% bias (see module docstring) and can be locally non-monotone
    there; the interpolation below is only meaningful above that floor.
    """
    bins = np.asarray(bins, float)
    errs = np.asarray(errs, float)
    ok = np.isfinite(errs)
    bins, errs = bins[ok], errs[ok]
    below = errs <= target
    if below.all():
        return bins[0]
    if not below.any():
        return np.inf
    # first index (ascending bins) that meets target
    i = np.argmax(below)
    if i == 0:
        return bins[0]
    x0, x1 = np.log(bins[i - 1]), np.log(bins[i])
    y0, y1 = np.log(errs[i - 1]), np.log(errs[i])
    xt = x0 + (np.log(target) - y0) * (x1 - x0) / (y1 - y0)
    return np.exp(xt)


fig, axes = plt.subplots(1, len(MODELS), figsize=(10.5, 4.4), squeeze=False)
print("=== old-scheme bias vs new truth (max|dP/P| over k, and value at k=0.01) ===")
print("# measures master's known +1.9% <beta^2> velocity excess (q^3-tilt placement")
print("# bias in its D-renormalized Gaussian deposit), NOT a build-consistency gate.")
print("# See .superpowers/sdd/task-3-diagnosis.md for the background-level proof.")
print(f"{'model':<16} {'max|dP/P|':>12} {'at k=0.01':>12}")
summary = {}
blocked = []
for ax, (name, label) in zip(axes[0], MODELS.items()):
    ref_pk, ref_cl = new[f"{name}_b192_pk"], new[f"{name}_b192_cl"]
    old_pk_192 = old[f"{name}_b192_pk"]
    bias_full = np.nanmax(np.abs(old_pk_192 / ref_pk - 1.0))
    bias_lowk = abs(old_pk_192[I_LOWK] / ref_pk[I_LOWK] - 1.0)
    print(f"{name:<16} {bias_full:>12.2e} {bias_lowk:>12.2e}")
    if bias_lowk > LOWK_SANITY_TOL:
        blocked.append((name, bias_lowk))

    e_old = np.array([max_err(old, name, b, ref_pk, ref_cl)[0] for b in SWEEP])
    e_new = np.array([max_err(new, name, b, ref_pk, ref_cl)[0] for b in SWEEP])
    summary[name] = (e_old, e_new)

    ax.loglog(SWEEP, 100 * e_old, "o-", color=C_OLD, lw=2, ms=7,
               mfc=C_OLD, mec=SURFACE, mew=1.2, label="uniform ln $q$ (old, biased)")
    ax.loglog(SWEEP, 100 * e_new, "s-", color=C_NEW, lw=2, ms=7,
               mfc=C_NEW, mec=SURFACE, mew=1.2, label="injection-adapted (new, truth)")
    for tgt in (1.0, 0.3, 0.1):
        ax.axhline(tgt, color=TARGET_LINE, lw=1, ls=":")
        ax.text(SWEEP[-1] * 1.03, tgt, f"{tgt:g}%", fontsize=7.5,
                 color=INK_MUTED, va="center", ha="left")

    ax.set_title(label, fontsize=10.5, color=INK, pad=8)
    ax.set_xlabel("momentum bins $N_q$", fontsize=9.5)
    ax.set_ylabel(r"max$_k\,|\Delta P_m/P_m|$ vs new-192  [%]", fontsize=9.5)
    ax.grid(True, which="major", color=GRID, lw=0.7)
    ax.grid(True, which="minor", color=GRID, lw=0.4, alpha=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xticks(SWEEP)
    ax.tick_params(labelsize=8.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(TARGET_LINE)
    leg = ax.legend(fontsize=8.5, frameon=False, loc="upper right")
    for text in leg.get_texts():
        text.set_color(INK_SECONDARY)

fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.suptitle("Momentum-grid convergence, measured against the validated truth (new-192)",
             fontsize=12.5, color=INK, y=1.04)
fig.text(0.5, 0.955,
          r"$P_m(k)$ relative error vs. momentum bins $N_q$; old's residual floor traces "
          r"its own $+1.9\%\ \langle\beta^2\rangle$ deposit bias, diluted by $P_m$ (not resolution)",
          fontsize=9, color=INK_SECONDARY, ha="center")
out = os.path.join(HERE, "dcdm_wdm_convergence.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print(f"\nsaved figure -> {out}")

if blocked:
    print("\nBLOCKED: old-scheme bias at k=0.01 exceeds the 1% sanity tolerance for:")
    for name, v in blocked:
        print(f"  {name}: {v:.2%} (expected <~1%, contradicts the low-k-anchor premise)")
    sys.exit(1)

print("\n=== bins needed to reach a P_m(k) precision target (new grid, vs new-192) ===")
print("# 'old bins' rows are the SAME old-scheme error curve, which floors at its own")
print("# ~1.9%-velocity-bias plateau (see header) rather than continuing to descend --")
print("# 'inf' / large values below the target than the floor are that bias, not a")
print("# resolution failure; annotated per-row.")
print(f"{'model':<16} {'target':>7} {'old bins':>9} {'new bins':>9} {'reduction':>10}  note")
for name in MODELS:
    e_old, e_new = summary[name]
    old_floor = float(np.nanmin(e_old))
    for tgt in (0.01, 0.003, 0.001):
        nb_old = bins_for_target(SWEEP, e_old, tgt)
        nb_new = bins_for_target(SWEEP, e_new, tgt)
        red = nb_old / nb_new if np.isfinite(nb_old) and np.isfinite(nb_new) and nb_new > 0 else np.nan
        note = "old below its bias floor -> unreachable" if tgt < old_floor else ""
        red_s = f"{red:>9.1f}x" if np.isfinite(red) else f"{'nan':>9}x"
        nb_old_s = f"{nb_old:>9.0f}" if np.isfinite(nb_old) else f"{'inf':>9}"
        print(f"{name:<16} {100*tgt:>6.1f}% {nb_old_s} {nb_new:>9.0f} {red_s}  {note}")

print("\n=== per-run compute time [s] (evolver=2) ===")
print(f"{'model':<16} {'bins':>5} {'old t':>8} {'new t':>8}")
for name in MODELS:
    for b in BINS:
        t_old = float(old[f"{name}_b{b}_t"])
        t_new = float(new[f"{name}_b{b}_t"])
        print(f"{name:<16} {b:>5} {t_old:>8.1f} {t_new:>8.1f}")
