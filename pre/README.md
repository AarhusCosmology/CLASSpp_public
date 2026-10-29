# Precision files

CLASS++ ships one set of compiled defaults (`include/precision.h`) and a small
number of `.pre` files that override them. Pass one as the second argument:

```
./class my_model.ini pre/precision_precise.pre
```

Every number quoted in these files was measured, not estimated. The metric is

> **S** = the standard deviation of Δχ² against a converged reference, over 120
> draws from the Planck 2018 `base_mnu` posterior, through the Planck 2018
> `plik_lite` TTTEEE likelihood.

A constant offset in Δχ² cancels in the posterior and is not penalised; what
biases parameters is how Δχ² *varies* across the parameter space, and S measures
exactly that. S ≈ 2 × the induced parameter bias in σ, so S = 0.5 is roughly a
0.1σ shift. The reference is itself converged to S = 0.008.

| file | S | speed vs defaults | use it for |
|---|---|---|---|
| *(compiled defaults)* | 0.387 ³ ⁴ ⁵ | 1.00× | everything |
| `precision_precise.pre` | 0.009 | ~60× slower | reference-grade, for convergence tests |
| `precision_pre_2026_08_21.pre` | 1.725 | 2.2× slower | reproducing results from before 2026-08-21 |
| `precision_pre_2026_10_09.pre` | 0.387 | 2% faster on one thread | reproducing results from 2026-08-21 to 2026-10-09 |
| `precision_pre_2026_10_21.pre` | 0.387 ³ | 15% slower on one thread, 20% with halofit | reproducing results from 2026-10-09 to 2026-10-21 |
| `precision_pre_2026_10_27.pre` | 0.387 ³ | 8% slower on one thread, 61% on w0waCDM | reproducing results from 2026-10-21 to 2026-10-27 |
| `act_dr6.pre` | — ¹ | 3.3× slower ² | ACT DR6 and SPT-3G D1: temperature, polarisation and lensing |

¹ Not scored on `plik_lite`. Its own metric is the distance to a converged reference
through ACT-lite, SPT-3G D1 lite, ACT DR6 + Planck PR4 lensing and the φφ part of SPT-3G
MUSE, at a ΛCDM point and at the best fit of self-interacting neutrinos + DRMD to the H0
World Cup data, as the bound it gives on the parameter shift. For the four likelihoods
together: 0.12σ and 0.16σ, against 0.19σ and 0.22σ for the defaults. The version before
2026-10-29 was 0.14σ and 0.54σ: its neutrino lines, the defaults of 2026-08-21, put all
three self-interacting neutrinos on a short hierarchy. ACT-lite uses no multipole above
ℓ = 6325 although it asks for 9001; the file's header says how to cut the request. See
`docs/superpowers/specs/2026-10-29-act-dr6-small-scale-preset-design.md` (and, for the
first version, `docs/superpowers/specs/2026-09-30-act-dr6-precision-design.md`).
² ΛCDM at `l_max_scalars = 6327`, one thread, both sides: 4.8 s against 1.5 s (7.6 s before 2026-10-29).
³ Measured on the defaults of 2026-08-21. On 2026-10-09 three neutrino knobs went
back to their earlier values (see the table under **Tensors** below); S has not been
re-measured since.
⁴ On 2026-10-21 the cuts on the scalar CMB transfer functions became l-dependent: the
transfer stage is 25% faster and, on Planck PR4 CamSpec + ACT DR6 lensing + DESI DR2, the
cosmological bias of the defaults went from 0.175 to 0.153σ. Not re-measured on
`plik_lite`. Every preset above pins the constant cuts it was measured with. See
`docs/superpowers/specs/2026-10-21-transfer-neglect-per-l-design.md`. The same day the
linear sources were continued past the last perturbation mode: the CMB lensing source,
which had stopped there and left C_L^φφ 2–3% low at L ~ 3000, and the linear spectrum
halofit and HMcode integrate, which let the non-linear floor come down from 5 to 2/Mpc.
On the same likelihoods the bias went to 0.148σ. In runs without a non-linear method it
moves lensed TT and EE by at most 4e-5 and 7e-6, so the measured presets are not pinned
for it; the `precision_pre_*` files are. See
`docs/superpowers/specs/2026-10-21-lensing-source-extrapolation-design.md`.
⁵ On 2026-10-27 four things changed. The C_l are interpolated in ℓ by a spline built for the
acoustic oscillation, and `l_linstep` went from 40 to 50 (4% faster on one thread). The
radiation-streaming and visibility triggers went from 26 and 0.008 to 45 and 2e-3: lensed EE
had been 8.3e-3 off a converged reference, 2.3σ for a cosmic-variance-limited measurement, and
is 1.2e-3 off now. The flat-space Bessel functions are sampled at 12 points per period instead
of 8. Dark energy fluids run the explicit rkdp45 perturbations evolver, as accurate as ndf15
and 35% faster on w0waCDM.

Against a reference converged further than the one above, on the same likelihoods, the
cosmological bias of the defaults is now 0.138σ for ΛCDM (0.077 before) and 0.144σ for
w0waCDM (0.137). The spread of Δχ² is 0.42 (0.72) and 0.41 (0.84). ΛCDM runs at the same speed
on one thread. The lower bias before came from the cheap triggers cancelling the error of the
Bessel interpolation and of three other transfer-stage discretizations.

The same day the k range of the scalar C_l stopped following the lensing headroom, which went
from 500 to 1000 multipoles (see **Tensors** below). `act_dr6` and `precise` set their own
`l_linstep`, triggers and headroom. The
`precision_pre_*` files pin the old values, and their headers say how to restore the old k range
and the fluid's evolver. See `docs/superpowers/specs/2026-10-27-precision-sync-design.md`.

## Two things these files do not cover

**Low-ℓ EE.** `plik_lite` starts at ℓ = 30, so nothing below that entered the
objective these settings were tuned on. `precise` leaves `k_step_super` at a
cheap rung, which puts the ℓ = 2–29 EE residual at 1.14× the full-sky cosmic
variance of a multipole. τ is measured there. Set `k_step_super = 0.002` before
pairing it with a low-ℓ EE likelihood (Planck lowE, SRoll2); it costs about 10%. **The compiled defaults do not have
this problem** — they carry `k_step_super = 0.002` for exactly this reason.

**Nonlinear.** `non linear = none` throughout. The halofit and HMcode knobs were
never exercised and are at their own defaults. Linear P(k) *was* measured: the
defaults are accurate to ~0.2% rms over 10⁻³ < k < 10 h/Mpc. An rms hides the high-k
end, though: with massive neutrinos the defaults of 2026-08-21 were 3.8e-3 low at
k = 5/Mpc, which is why three neutrino knobs went back on 2026-10-09. The *nonlinear*
spectrum is a different matter. Until #399, `nonlinear_min_k_max = 5` truncated the
integrals halofit needs, costing 5–9% at k ≳ 0.5 /Mpc. Since 2026-10-21 the linear
spectrum is continued past the last mode in its high-k form and trusted up to
`k_max_extrapolated` = 70/Mpc, so the floor came down to 2/Mpc and P_NL stays within
6e-4 of modes solved to 70/Mpc up to z = 3. Only `act_dr6.pre` touches these knobs: it
solves the modes to 10/Mpc, because the continuation has the high-k form of cold dark matter,
which dark acoustic oscillations break.

**What the presets still carry.** The compiled defaults hold seven knobs back from
the campaign's configuration (table above). `precise` no longer overrides
`k_max_tau0_over_l_max`, because a reference-grade preset that is 12% wrong on
tensor BB is not a reference; that cost S 0.0083 → 0.0087. Its header names the
low-ℓ EE knob it still leaves cheap, and the line that fixes it.

`precision_fast.pre` and `precision_tuned.pre` were removed on 2026-10-27. They
were the campaign's own settings (S 0.554 at 2.4× the speed of the defaults of the
time, and 0.301 at 1.7× slower), so they reintroduced the low-ℓ EE, P(k)-truncation,
lensing-headroom and tensor-BB issues, and since then the defaults have become fast
enough that they no longer earned their place.

**Tensors.** The campaign scored scalars only, at ℓ ≥ 30, with no `mPk`. **Seven**
knobs are therefore held at their pre-2026-08-21 values in the compiled defaults,
because that objective could not constrain them:

| knob | what it would have broken |
|---|---|
| `k_step_super` | low-ℓ EE, hence τ — the residual reached 1.6× cosmic variance |
| `k_min_tau0` | P(k) below k ≈ 3×10⁻⁴ h/Mpc, to 2.5% |
| `delta_l_max` | lensing headroom — the lensed spectrum degrades to 7% at the top of the requested range (1000 since 2026-10-27) |
| `k_max_tau0_over_l_max` | tensors need a wider k range per ℓ than scalars; a tensor scenario went to 2.5% |

| `ncdm_fluid_trigger_tau_over_tau_k`, `l_max_ncdm`, `tol_ncdm_synchronous` | P(k) with massive neutrinos: 1.3e-3 low at k = 1/Mpc, 3.8e-3 at k = 5/Mpc (restored 2026-10-09) |

Each was caught by running the reference-comparison workflow and asking, for
every failing scenario, whether the new defaults are *closer to a converged
reference* than the old ones — not merely different from them.

Since 2026-10-27 the scalar k range follows the multipoles asked for, not the lensing headroom
above them: k_max τ0 = max(`k_max_tau0_over_l_max` l_max, l_max + `k_max_tau0_above_l_max`).
`act_dr6` and `precise` take the new rule (`precise` 30% faster, the same spectra to 0.1 CV).

A fifth, `tight_coupling_trigger_tau_c_over_tau_k`, was **split by mode** instead
of held back. It cost 14% on tensor BB over ℓ = 100–300 at the campaign's value,
because while the tight-coupling approximation is on the tensor equations set the
photon polarisation quadrupole — the entire B-mode source — identically to zero,
and unlike the scalar case there is no higher-order closure to fall back on
(every `tight_coupling_approximation` scheme gives bit-identical tensor output;
the selector is only ever consulted from the scalar branch). Scalars keep the
loose, measured-good 0.04; tensors get their own
`tight_coupling_trigger_tau_c_over_tau_k_ten`, still 0.01. Splitting rather than
reverting keeps 26% of the runtime, and follows the `_ten` convention CLASS
already uses for `l_max_g_ten` and `l_max_pol_g_ten`.

## Reproducing these numbers

The harness lives on the **`benchmarks/precision-campaign`** branch, under
`benchmarks/precision/` — it is deliberately kept off `master`, like
`benchmarks/class-public-comparison`. `s23_write_pre.py` regenerates this
directory; `paths.py` takes its locations from the environment. The campaign
report explains the method and the measurement behind every claim above:

<https://claude.ai/code/artifact/bc8ddc19-fd2d-4de7-a05c-58c974af45b9>
