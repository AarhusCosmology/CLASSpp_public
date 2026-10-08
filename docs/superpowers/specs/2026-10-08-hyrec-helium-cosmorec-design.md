# HyRec's helium, calibrated to CosmoRec — and two driver faults found on the way

**Status (2026-10-08).** Shipped on branch `hyrec-helium-cosmorec` (PR #449): the
recombination driver reports the ODE solution instead of a blend back to Saha, with
the Saha→ODE switches at x = 0.999 (§4); HyRec's helium relaxes to zero below
XHEII_MIN instead of freezing there (§5); and `hyrec_calibration = cosmorec` (§6).
**Later the same day** both CosmoRec calibrations became the compiled defaults —
`hyrec_calibration = cosmorec` and `recfast_calibration = cosmorec` — on the
maintainer's rule that accurate settings which cost nothing should be the default
(§8.1). That supersedes `2026-09-30-act-dr6-precision-design.md` §4.2's "`planck`
(the default)". The old values are in `pre/precision_pre_2026_08_21.pre`. **Not
done:** the open residuals of §8.2–8.5 are measured, not explained. The harness is on
the not-for-merge branch `benchmarks/hyrec-helium`.

Evidence tags: **[M]** measured here, **[C]** the code says so, **[L]** published,
quoted with its source, **[?]** no evidence, do not act on it.

## 1. The question

`docs/superpowers/specs/2026-09-30-act-dr6-precision-design.md` §3 recorded that
HYREC-2 does not track CosmoRec through He I recombination: 4.5e-3 in x_e, against
7.8e-5 for CAMB 2's refitted RECFAST, and Δχ² = 3.8 on an SO-like likelihood [L].
Can HyRec's helium be brought to CosmoRec — ideally by adding a few equations?

## 2. Sources

- **CAMB 2**: Lewis, arXiv:2607.14854 §V.B and fig. 1: the CosmoRec reference
  (accuracy 6, ten H shells, nS_2γ = 8, nS_Raman = 7, two diffusion iterations), the
  RECFAST refit and its He I rate correction F_He (eq. 23), and the HyRec-2 numbers.
- **CFS12**: Chluba, Fung & Switzer, arXiv:1110.0247, helium radiative transfer.
- **CosmoRec** v3.0 from `cmbant/CosmoRec`, branch `camb` (6f28c15) — the version
  Lewis used. Its batch API exposes the helium physics switches used in §3.
- **HYREC-2** at the vendored commit 09e8243d, built standalone.

## 3. What makes HyRec's helium differ

**The reference reproduces Lewis [M].** Standalone HYREC-2 against our CosmoRec
reference (massless ν, no reionization, ω_b = 0.0225, ω_c = 0.1193, h = 0.676,
Y_He = 0.2455): max |Δx_e/x_e| = 3.50e-4 over 800–1600 and 4.52e-3 over 1600–3000.
Lewis quotes 3.46e-4 and 4.51e-3.

**Attribution [M].** CosmoRec's `accuracy = -1` is documented as "closest
equivalent of HyRec". Starting there and switching one piece of helium physics on
at a time (He I window, max |Δx_e/x_e| to the reference):

| CosmoRec helium setting | He I window |
|---|---|
| HyRec-like: 2 He shells, feedback to n = 2, H I absorption on | 4.13e-3 |
| + H I absorption with the diffusion correction | 1.61e-3 |
| + H I absorption by full radiative transfer | 4.45e-4 |
| + 5 He shells (alone) | 3.88e-3 |
| + 5 He shells and full radiative transfer | 4.2e-6 |
| standalone HYREC-2 | 4.52e-3 |

About 90% of the gap is how H I continuum absorption of the 584 Å photons is
treated. HYREC-2 uses an analytic escape probability (`helium.c`: Sobolev escape
enhanced by a fit in the continuum optical depth, plus the 591 Å intercombination
line) [C]. CosmoRec solves the transfer: resonance scattering diffuses photons into
the line wings, where H I absorbs them, and helium recombines faster. HYREC-2's
x_e therefore runs high, by up to 4.7e-3 near z = 1850 [M].

**Why not "a few differential equations".** The missing physics is a PDE in
frequency and time (line overlap, electron scattering, delayed feedback). Its
escape-probability reduction is what HYREC-2 already has. Promoting the n = 2
helium levels from Saha to evolved populations would be a few equations, but the
table shows it is the *shells* row — a small part of the gap.

## 4. Driver fault 1: x_e blended back towards Saha

Found while checking that CLASS++'s HyRec is HYREC-2 [M]. Both CLASS++ recombination
models, but neither CosmoRec nor HYREC-2, showed x_e low by up to 1.3e-3 at
z ≈ 1500–1640 and 5e-4 at z ≈ 2630–2850.

**Cause [C].** After each Saha→ODE switch (at Saha x = 0.995), `thermodynamics_recfast_output`
reported a smoothstep blend of the ODE value with the Saha value until x had fallen
another 0.05. Recombination lags Saha, so the blend biased x_e low across exactly
those ranges. Moving the switch alone does not help; it moves the blend window.

**Fix.** Report y[0] + f_He y[1] from the first step (each phase starts from the
Saha value, so there is nothing to smooth); retire `recfast_x_{H,He}0_trigger_delta`
with an error; move the switches to 0.999, which removes the kink the ODE start left
in x_e(z) (max |Δ² ln x_e| 2e-5 → background curvature) at no cost [M].

**Effect [M].** CLASS++ HyRec's hydrogen window 1.31e-3 → 2.5e-4 (standalone
3.5e-4); RECFAST's 1.33e-3 → 2.2e-4. Lensed C_ℓ move by 3e-4 on ACT-lite and
0.036 on the SO-like distance, alike for RECFAST, RECFAST cosmorec and HyRec.

**A robustness fault, and its real cause [M].** With `tol_thermo_integration = 1e-9`,
master HyRec failed in ndf15 (singular LU, or a bogus z_rec < 500) on **36 of 100**
random ΛCDM cosmologies; with `reio_none` at the default tolerance, 2 of 200. RECFAST:
0. On master's wrapper, moving the hydrogen switch to 0.999 alone took the first to
0/100 (the helium switch and the blend changed nothing) — which first read as the fix.
It is not: with §5's relaxation in place, all three switches back at 0.995 also give
0/100, while the freezing wrapper with the hydrogen switch at 0.995 gives 37/100. The
cause is §5's zero-derivative trap; the earlier hydrogen start only steered ndf15's
path around it. `test_hyrec_at_tight_thermo_tolerance` pins two of the failing
cosmologies, checking z_rec against RECFAST's.

## 5. Driver fault 2: helium frozen at XHEII_MIN

**Cause [C].** Below XHEII_MIN = 1e-6 HeII per hydrogen, HYREC-2 takes helium to be
entirely neutral (`history.c`); CLASS++'s wrapper returned a zero derivative, freezing
it at 1e-6. That added 1e-6 electrons per hydrogen to the whole tail (CLASS++ HyRec
+8.4e-4 at z = 200 against CosmoRec where standalone HYREC-2 is −7.2e-4) [M], and made
the cut a trap: a stiff step landing below it, even below zero, was self-consistent
and stayed. The calibrated rate of §6 hit it at N_eff = 2.5: x_He = −7e-3, then a
reionization-sampling failure [M].

**Fix.** dx_HeII/dlna = −100 x_HeII below the cut, blended linearly with HYREC-2's
rate over [XHEII_MIN, 2 XHEII_MIN]. The tail now follows standalone HYREC-2 to 5e-5.

## 6. The calibration

**Form.** CAMB 2's: dx_HeII/dlna ← F dx_HeII/dlna with
F = 1 + (a0 + a1 u + a2 u²)/(1 + d2 u² + d4 u⁴), u = (z_scale − z0)/w, and z_scale
the redshift at which a 2.7255 K sky has the current photon temperature, in atomic
units so it follows He I under varying α and m_e. Applied in `source/hyrec_model.cpp`;
the vendored HYREC-2 is untouched. No z-window: with CAMB's hard 1500–3000 window
the fit's starting constants put a 10–20% step in the rate at z = 1500 [M], and one
fit evaluation then failed in ndf15 at 1e-9 — whether because of the step is [?],
since the trap of §5 was found only afterwards. F → 1 smoothly anyway (1.005 at
z = 3000), so the window buys nothing.

**Fit [M].** Least squares on x_e/x_e(CosmoRec) − 1 over 1550 ≤ z ≤ 3000 (Δz = 5) at
the fiducial cosmology of §3, CLASS++ at `tol_thermo_integration = 1e-9`, six
starts, F > 0.3 enforced everywhere. All starts land at max 1.5–2e-4, rms 8–10e-5:
the form is degenerate at this level. Kept: (a0, a1, a2, d2, d4, z0, w) =
(0.635977, 0.298024, 0.239789, 2.91828, 3.31521, 1751.94, 286.606). F peaks at 1.65
near z = 1780. Reproduced by the hard-coded constants at the default tolerance:
max 1.55e-4, rms 8.1e-5.

**What it leaves [M].** A +1.5e-4 plateau over z ≈ 2300–2700, where helium is still
near Saha and a rate factor barely moves x_e (F = 1.02 there moves it 2.0e-4 → 1.5e-4).
Standalone HYREC-2 has it too; CosmoRec's HyRec-like mode does not (−1e-4 at 2400).
So it is something in HYREC-2's quasi-equilibrium helium, not radiative transfer. [?]
on what.

**Transfer across cosmologies [M].** Each against its own direct CosmoRec run, He I
window max |Δx_e/x_e|:

| cosmology | HyRec | HyRec, calibrated | RECFAST cosmorec |
|---|---|---|---|
| fiducial | 4.70e-3 | 1.55e-4 | 3.38e-4 |
| ω_b 0.020 / 0.025 | 4.69 / 4.71e-3 | 4.79 / 4.73e-4 | 3.09 / 8.23e-4 |
| ω_c 0.10 / 0.14 | 4.67 / 4.73e-3 | 1.89 / 1.67e-4 | 5.41 / 1.46e-4 |
| Y_He 0.22 / 0.28 / 0.35 | 3.99 / 5.76 / 8.35e-3 | 2.26 / 4.25e-4 / 1.41e-3 | 1.76 / 7.91e-4 / 2.41e-3 |
| N_eff 2.5 / 3.8 | 4.68 / 4.71e-3 | 1.58 / 1.55e-4 | 4.02 / 2.53e-4 |
| h 0.60 / 0.76 | 4.70e-3 | 1.55 / 1.56e-4 | 3.38e-4 |
| α 0.99 / 1.01 | 4.80 / 4.59e-3 | 1.74 / 2.03e-4 | 3.26 / 8.30e-4 |
| m_e 0.95 / 1.05 | 4.57 / 4.82e-3 | 1.65 / 1.63e-4 | 8.95 / 5.26e-4 |
| corners lo / hi | 3.94 / 5.82e-3 | 1.54 / 1.91e-4 | 1.22 / 9.04e-4 |

10–30× everywhere; ω_b and Y_He carry the residual cosmology dependence, which a
temperature-only factor cannot follow.

**C_ℓ [M].** Lensed TT/TE/EE at `pre/act_dr6.pre` precision, one massive ν; the
calibrated HyRec stands in for CosmoRec (CLASS++ cannot be fed an external x_e):

| | ACT-lite | SO-like (Lewis' proxy, ℓ ≤ 4000) |
|---|---|---|
| HyRec (this branch) vs calibrated | 0.013 | 2.35 |
| HyRec (master) vs calibrated | 0.008 | 1.47 |
| RECFAST planck vs calibrated | 0.045 | 0.78 |
| RECFAST cosmorec vs calibrated | 0.003 | 0.062 |

Lewis' corresponding numbers against direct CosmoRec are 3.8, —, 0.56 and 0.03 [L].
The last row is the useful one: two independent routes to CosmoRec agree to 0.06 on
an SO-like likelihood. Master's HyRec looks closer than this branch's uncalibrated
HyRec only because its two driver faults partly cancel its helium error.

**Robustness [M].** 0 failures in 300 draws from a wide prior (ω_b 0.017–0.028, ω_c
0.08–0.16, h 0.55–0.8, Y_He 0.18–0.35, N_ur 2–4.5) for each of: calibrated and
uncalibrated with `reio_none` at the default tolerance, calibrated with reionization at
1e-9, and RECFAST cosmorec at 1e-9.

## 7. Alternatives rejected

- **The physics.** §3: a radiative-transfer PDE, which is CosmoRec.
- **CosmoRec inside CLASS++.** ~5 s per cosmology at the reference settings, GSL, global
  solver state, and a licence that does not plainly allow redistribution in our wheels.
  It is the reference here, run out of process (`cr_driver.cpp` in the harness).
- **Fitting to CAMB's `recfast_cosmorec`.** A fit to a fit; it matches CosmoRec only at
  Lewis' fiducial. Direct CosmoRec was 5 s a run.
- **Correcting the 2¹P escape probability inside `helium.c`.** More physical, possibly
  better on ω_b and Y_He, but it edits the vendored code (PROVENANCE.md keeps one
  portability patch) and the temperature-only factor already transfers 10–30×.
- **Lewis' constants.** They correct RECFAST's helium, which errs differently.

## 8. Open

1. **The default — decided.** `cosmorec` for both codes (Status). It moves HyRec
   C_ℓ by 0.013 (ACT-lite) and 2.35 (SO-like), RECFAST's by 0.027 (ACT-lite) and
   0.46 (SO-like; Lewis: 0.56 against direct CosmoRec). At default precision,
   master's lensed C_ℓ move by at most 6.5e-4 (RECFAST) and 1.7e-3 (HyRec, low-ℓ
   EE). class_public, which runs HyRec as distributed, now differs from CLASS++ by
   that much at default settings.
2. **The 2300–2700 plateau** (§6), 1.5e-4: unexplained.
3. **CLASS++'s RECFAST cosmorec** sits 3.4e-4 from CosmoRec in the He I window where
   CAMB's sits at 7.8e-5, and its hydrogen window reaches 1–3e-3 off the fiducial
   (7e-3 at m_e ± 5%) where HyRec stays near 2e-4. Not run down; a CLASS++-specific
   refit would be the obvious step.
4. **HyRec's tail** is still −7e-4 from CosmoRec at z = 200: HYREC-2's own hydrogen,
   shared with the standalone code.
5. **Varying constants** compare CLASS++'s instantaneous transition with CosmoRec's
   "rescale everything" mode; their conventions were not checked against each other.
