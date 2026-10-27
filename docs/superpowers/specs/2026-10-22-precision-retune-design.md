# Precision retune of 2026-10-22, and the knobs shared with class_public

**Status (2026-10-27).** Not shipped as written: the pull request that carried it (#473) was closed
unmerged, and `2026-10-27-precision-sync-design.md` replaced it. From this note, three things shipped there:
the dark energy fluid's opt-in to rkdp45 (section "The fluid and the explicit evolver"; re-measured to the same
numbers), and radiation streaming 45 with the visibility cut at 2e-3 rather than 1e-3 (last section).

Two things overtook the rest. The acoustic spline (#477) does at `l_linstep` 50 what the cubic spline did at
30. And a reference converged further (codecmp `ref2`; the reference used below is 0.030 sigma from it)
scores this retune, rebased onto master, at 0.165 sigma against master's 0.077: the loosenings of 2026-10-26
cost more than this note says. Against `ref2`, the error that the cheap triggers cancel turned out to be four
small discretizations in the transfer stage, of which the cubic interpolation of the flat Bessel functions is
the cheapest to fix (`hyper_sampling_flat` 8 -> 12 there). The four free knobs from class_public were left
out. The class_public side was never posted. Everything below is a dated record.

Evidence tags: **[M]** measured for this note, **[C]** read from the code, **[P]** an earlier
campaign, not re-measured, **[?]** no evidence.

## Metric

The three-code benchmark (`benchmarks/codecmp/`, untracked): Planck PR4 CamSpec TTTEEE + low-l
TT + LoLLiPoP EE + ACT DR6/PR4 lensing + DESI DR2, 160 LambdaCDM and 112 w0waCDM posterior
draws, each scored against a converged CLASS++ reference **[P]**. b_cosmo is the largest shift,
in posterior sigma, of any combination of cosmological parameters; S the spread of delta chi^2,
which includes what the nine CamSpec nuisance parameters absorb. Bootstrap errors on b_cosmo are
0.004-0.006. Screening drivers: `screen.py` with `variants_screen{1..5}.json` **[M]**.

## What moves the defaults

37 single-knob variants from master (b_cosmo 0.148 LCDM) **[M]**. Only the C_l sampling in l
matters by more than the noise:

| l_linstep | 40 | 35 | 30 | 25 | 20 |
|---|---|---|---|---|---|
| LCDM b_cosmo | 0.148 | 0.109 | 0.088 | 0.084 | 0.085 |

The C_l are computed at sampled multipoles and splined in between; the spline error at a linear
step of 40 is oscillatory in l, with the acoustic peaks, which CamSpec's 9915 multipoles see as a
shift of cosmological parameters. It converges monotonically towards the reference's 15, so this
is accuracy bought, not borrowed. The next-largest single levers, `l_logstep` 1.06 (0.138),
`hyper_sampling_flat` 12 (0.141) and `perturb_sampling_stepsize` 0.05 (0.125, +43%), add little
on top of l_linstep = 30 (0.082-0.093) **[M]**.

## The approximation triggers: cancellation (2026-10-22; reversed on 2026-10-26, see the last section)

class_public's values of four knobs are worse in CLASS++ than CLASS++'s own (LCDM, from master)
**[M]**:

| knob | CLASS++ | class_public | b_cosmo with class_public's |
|---|---|---|---|
| `radiation_streaming_trigger_tau_over_tau_k` | 26 | 45 | 0.219 |
| `tight_coupling_trigger_tau_c_over_tau_k` | 0.04 | 0.01 | 0.160 |
| `ur_fluid_trigger_tau_over_tau_k` | 20 | 30 | 0.158 |
| `neglect_CMB_sources_below_visibility` | 0.008 | 1e-3 | 0.177 |

The reference sits on the accurate side of all of them (240, 0.01, 130, 1e-4). With
l_linstep = 30, moving the triggers all the way there gives b_cosmo 0.147 instead of 0.088, but
halves S (0.75 -> 0.52), and with `perturb_sampling_stepsize` 0.05 as well b_cosmo 0.112 and
S 0.40 **[M]**. So the cheap triggers carry errors that cancel another error in the cosmological
directions; it is not the source time sampling, and finding it is beyond a retune. The triggers
keep their values. The ACT-lite measurement that blamed `radiation_streaming_trigger = 26` for
most of the defaults' error at l > 3000 still stands **[P]**.

## Free knobs taken from class_public

Neutral on b_cosmo and S within the noise, and free on one thread (interleaved, 3 rounds, median
of per-draw minima over 16 draws) **[M]**: `start_small_k_at_tau_c_over_tau_h` 0.006 -> 0.0015
(low-l S 0.012 -> 0.010), `q_logstep_limber` 1.1 -> 1.025 (the full-Limber grid noise in
C_L^phiphi, 3.1e-3 at 1.1 **[P]**), `tol_tau_approx` 1e-6 -> 1e-10, `hyper_phi_min_abs`
1e-8 -> 1e-10 (+0.5%). Not taken: `k_max_limber_over_l_max_scalars` 0.001 (+6% on one thread,
nothing gained: since #470 the lensing source is continued past the last mode),
`q_logstep_spline` 170, `l_switch_limber` 30, the thermodynamics and background tolerances, whose
meaning differs between the two codes, and `reionization_sampling`, which is a relative variation
here and a fixed dz in class_public **[C]**.

## The fluid and the explicit evolver

The automatic choice of rkdp45 needs every species to opt in, and `FluidSpecies` never had, so
every model with a fluid ran ndf15 at the tolerance tuned for rkdp45 **[C]**. Measured against
ndf15 at 1e-5, at the default tolerance 1e-4 **[M]**:

| fluid | ndf15: TT / EE / TE / P(k) | rkdp45: TT / EE / TE / P(k) | wall, 11 threads (nd 1e-4 -> rk) |
|---|---|---|---|
| w0wa, PPF | 2.5e-6 / 6.9e-6 / 1.9e-5 / 5.5e-5 | 2.0e-6 / 6.4e-6 / 8.1e-6 / 1.7e-5 | 0.228 -> 0.166 s |
| w = -0.9, cs2 = 1 | 1.9e-6 / 5.7e-6 / 1.6e-5 / 2.9e-5 | 4.0e-6 / 6.7e-6 / 6.1e-6 / 1.8e-5 | 0.294 -> 0.193 s |
| w = -0.5, cs2 = 0.01 | 4.9e-6 / 1.6e-5 / 1.4e-5 / 2.6e-5 | 4.7e-6 / 6.0e-6 / 6.3e-6 / 3.0e-5 | 0.234 -> 0.155 s |
| w0wa, cs2 = 1e-6 | 3.0e-6 / 6.7e-6 / 2.1e-5 / 3.0e-5 | 5.0e-6 / 5.8e-6 / 7.6e-6 / 1.9e-5 | 0.228 -> 0.152 s |
| w0wa PPF, newtonian | 4.7e-6 / 7.2e-6 / 2.2e-5 / 6.2e-4 | 6.7e-6 / 6.5e-6 / 1.1e-5 / 7.7e-5 | 0.327 -> 0.180 s |
| PEDE | 1.4e-6 / 6.4e-6 / 1.8e-5 / 5.7e-5 | 2.1e-6 / 6.6e-6 / 8.4e-6 / 3.3e-5 | 0.301 -> 0.204 s |
| GPEDE, PPF | 9.6e-7 / 7.0e-6 / 1.8e-5 / 4.3e-5 | 2.2e-6 / 6.2e-6 / 7.7e-6 / 1.8e-5 | 0.235 -> 0.169 s |
| pheno-axion, n = 3 | 1.3e-6 / 9.5e-6 / 1.6e-5 / 9.9e-6 | 1.3e-6 / 4.6e-6 / 5.2e-6 / 1.3e-5 | 0.320 -> 0.271 s |
| pheno-axion, n = 2 | 1.6e-6 / 5.4e-6 / 1.9e-5 / 1.1e-5 | 3.3e-6 / 4.7e-6 / 6.2e-6 / 2.3e-5 | 0.326 -> 0.263 s |

(C_l as the largest difference over l = 30-2500 relative to the peak; P(k) as the largest
relative difference at z = 0 up to 2/Mpc.) On the w0wa posterior draws the forced explicit
evolver scores b_cosmo 0.147 against 0.142 for ndf15, inside the noise, and the tolerance
1e-5 changes neither. `FluidSpecies` (and with it the PPF and pheno-axion fluids) now opts in;
the EDE fluid equation of state stops in the background, as before, with either evolver.

## Results

| | LCDM b_cosmo | S | S_lens | S_lowl | w0wa b_cosmo | S_lens | 1 thread LCDM | 1 thread w0wa |
|---|---|---|---|---|---|---|---|---|
| master | 0.148 | 0.783 | 0.055 | 0.012 | 0.158 | 0.070 | 1.048 s | 1.699 s |
| this | 0.088 | 0.747 | 0.026 | 0.010 | 0.145 | 0.031 | 1.144 s | 1.300 s |

**[M]**, timing interleaved over 3 rounds on 16 (LCDM) and 12 (w0wa) draws, laptop on mains
power. LambdaCDM is 9% slower, w0waCDM 24% faster.

## Not validated

Tensors, low-l EE and P(k) are unaffected by construction: the new l step only adds sampled
multipoles, and the earlier start of small k is the accurate direction **[C]**. The lvl2
reference comparison will show the l_linstep change in every C_l scenario.

## 2026-10-26: accurate triggers, loosened elsewhere

**Why the triggers changed after all.** The posterior metric cannot see the price of the
cancellation above; a converged reference can. Over the 160 LambdaCDM draws (median of the
largest error over l = 30-2500, and the cosmic-variance significance sqrt(sum (2l+1)/2 (dC/C)^2)),
radiation streaming 26 and the visibility cut 0.008 put lensed EE 8.1e-3 off the reference, a
significance of 2.3, and TE 6.2e-3 (1.1); each trigger alone does it **[M]**. Radiation streaming
26 was also most of the defaults' error on ACT DR6 **[P]**. Back at 45 and 1e-3 the bias rises to
0.156 sigma (LCDM) and 0.164 (w0wa) for 3% on one thread **[M]**. Buying that bias back costs
+27% (2026-10-22 option C: `hyper_sampling_flat` 12, the transfer cuts x1.3,
`perturb_sampling_stepsize` 0.085, 0.089 sigma). The target here is the opposite trade: keep the
accurate spectra, accept a bias below 0.2 sigma, and be faster than before.

**Where the time goes**, one thread, 2026-10-22 defaults, LCDM (`stages.py`): thermodynamics
0.027 s, perturbations 0.413, non-linear 0.025, transfer 0.530, spectra 0.017, lensing 0.070, of
1.08 s **[M]**.

**Every loosening was checked twice**: on the posterior, and against a converged reference on
eight models (`adjudicate_cpp/run2.py`, `compare2.py`: LCDM; 3 x 0.06 eV, 1.5 x 0.03 eV and
1 x 0.3 eV massive neutrinos; Omega_k = +-0.01; a scalar field; w0wa; lensed and unlensed C_l to
l = 3000 in cosmic-variance units, linear and halofit P(k) to 5/Mpc at z = 0 and 1), plus tensor BB
and l = 2-29 (`adjudicate_cpp/sectors_cpp.py`) **[M]**. The reference is `precision_precise.pre`
with `k_step_super` 2e-3, accurate lensing and `delta_l_max` 1800; for the tensors and low l,
`k_step_super` 5e-4 and `k_max_tau0_over_l_max` 3.

| knob | one thread | what the reference says |
|---|---|---|
| `k_step_transition` 0.4 -> 0.6 | -6% | unchanged |
| `k_step_super` 0.002 -> 0.003 | -2.8% | unchanged at l = 2-29 (largest EE residual 0.007 of cosmic variance per multipole); scalar-field raw TT 1.20 -> 1.31 |
| `l_switch_limber` 40 -> 30 (class_public's) | -2% | unchanged |
| `ncdm_fluid_trigger_tau_over_tau_k` 31 -> 25 | -2% | closer in P(k) on all three neutrino models (3 x 0.06 eV, k > 1/Mpc: 5.3e-3 -> 0.9e-3); C_l^phiphi 0.74 -> 0.91 in cosmic-variance units (1.06 on master) |
| `perturb_sampling_stepsize` 0.1 -> 0.12 | -7.5% | bias 0.152 -> 0.170 sigma; EE and TT unchanged, TE 0.49 -> 0.56 |

Rejected, though most looked free on the posterior **[M]**: `tol_perturb_integration` 4e-4 and
1e-3 (the scalar field's P(k) 16-94% off, C_L^phiphi 30x worse); `l_max_ncdm` 14 (3 x 0.06 eV:
P(k) 2e-2); `delta_l_max` 400 and `k_max_tau0_over_l_max` 1.2 (TT at l = 2500-3000, significance
1.81 -> 2.94 and 2.83); `tight_coupling_trigger_tau_c_over_tau_k` 0.06 (bias 0.129, but worse at
high l: another cancellation); `hyper_sampling_flat` 6 (0.243), the transfer cuts x0.8 (0.258),
`start_sources_at_tau_c_over_tau_h` 0.012 (0.246), `q_linstep` 0.55 (1.63), `ur_fluid_trigger`
15 (0.196); `l_max_g` 10 and 8 (EE worse, no time saved); radiation streaming 35 and visibility
2e-3 (EE worse); `start_large_k_at_tau_h_over_tau_k` 0.1 (nothing). `l_linstep` 35 instead of
the sampling step costs about as much bias for less time (0.166, -5%); both together give 0.176.

**`k_max_tau0_over_l_max` 1.5, measured and not taken.** It was free in everything above, which
always asks for 500 multipoles more than it compares, and it is the largest saving (-9% on one
thread; on the final defaults 1.5 / 1.65 / 1.8 are -15.5% / -11.5% / -6.4% against master).
Outside that it costs **[M]**: tensor BB 1.2% -> 2.6% at l = 30-300; without a margin the top
~250 multipoles of lensed TT (`l_max_scalars` 2000: largest error 3.0e-3 against 1.6e-3,
significance 0.82 -> 1.43, master 0.97); and at small `l_max_scalars` the TT near l_max, which is
already far off at 1.8 (`l_max_scalars` 200: 10% low at l = 200, 3% at l = 60; 17% and 5% at 1.5).
That last error is master's too and the converged reference shows it as well (3.4% at l = 200
when asked for l_max 200). The likely cause [?] is that k_max ∝ l_max cuts the j_l tail of sources
Silk damping has not yet killed. CAMB takes k eta_0 up to 2.5 l_max. Left for the CAMB comparison
(`benchmarks/camb_anatomy/BRIEF.md`).

**Results**, same metric as above, timing interleaved (3 rounds x 16 draws LCDM, 2 x 16 w0wa):

| | LCDM bias | S | w0wa bias | S | lensed EE max / significance | TE significance | 1 thread LCDM | w0wa |
|---|---|---|---|---|---|---|---|---|
| master | 0.148 | 0.783 | 0.158 | 0.821 | 8.5e-3 / 2.30 | 1.16 | 1.041 s | 1.70 s |
| 2026-10-22 step | 0.088 | 0.747 | 0.145 | 0.861 | 8.1e-3 / 2.33 | 1.12 | 1.141 s | 1.29 s |
| option C (not shipped) | 0.089 | 0.554 | 0.097 | 0.502 | 1.5e-3 / 0.98 | 0.40 | 1.448 s | 1.60 s |
| **this** | 0.170 | 0.595 | 0.179 | 0.566 | 1.3e-3 / 0.89 | 0.56 | 0.975 s | 1.098 s |

**[M]**, master and this measured in one interleaved run on the compiled defaults (master = its
values passed as input; for w0wa master itself, with the fluid on ndf15; master's values with the
fluid on rkdp45 take 1.200 s). The bias is master's plus 0.02 sigma with no cancellation behind it.
On the eight models every spectrum is at least as close to the reference as master's (to 2%), with
two exceptions: lensed TT for Omega_k = -0.01 (significance 2.20 -> 2.42, from the triggers) and the
scalar field's halofit P(k) at z = 1 above 1/Mpc (6.5e-3 -> 1.0e-2, from the sampling step).
One thread is 6% faster than master on LCDM, 35% on w0wa.

## Open: CAMB

CAMB 2.0.5 reaches b_cosmo 0.083 sigma on the same metric in 0.60 s on one thread, measured on
2026-10-21, when CLASS++ master (before #469 and #470) took 1.28 s in the same harness **[P]**;
CAMB against these defaults has not been timed directly. Where that comes from (k-grid placement, steps and failed steps per mode, equations
per step, the absence of interpolation in CAMB's stepping, the transfer stage) is the next
investigation, before the same values go to class_public.
