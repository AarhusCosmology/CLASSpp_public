# Interpolating the CMB spectra in l with their acoustic oscillation (2026-10-27)

**Status 2026-10-27:** implemented (branch `spectra/acoustic-l-interpolation`). No precision value changes in this
branch; what it allows for `l_linstep` is measured below and left to a separate decision. Evidence tags: [M]
measured here, [C] stated by the code, [?] not measured.

**Status 2026-10-27, later:** merged (#477), and `l_linstep` 50 is the default since the precision change that
followed (`2026-10-27-precision-sync-design.md`: bias 0.077 -> 0.090 sigma against a converged reference, 4%
faster on one thread).

## Problem

CLASS computes C_l at sampled multipoles (log steps, then every `l_linstep`) and splines them in l. The cubic
spline treats the acoustic oscillation as interpolation error, which grows as (h omega)^4 with omega = 2 theta_s.
`l_linstep` 40 -> 30 was the single largest accuracy gain of the 2026-10 retune (cosmological bias
0.148 -> 0.088 sigma), and it costs time: step 50 instead of 30 is 11% faster [M].

CAMB samples every 50 and stays accurate because it splines the difference to a fixed fiducial LCDM template
[C] (`InterpolateClArrTemplated`; `benchmarks/camb_anatomy/FINDINGS.md`). A fixed template was not wanted: it is
a data file, and models whose peaks move could interpolate worse with it. Thomas's suggestion was a template made
from the spectrum's own samples: a power law times a cosine with an l-dependent phase.

## Approach

`AcousticSpline` (`include/acoustic_spline.h`) interpolates u_l = l(l+1) C_l / S(l) with a trigonometric spline:

* **S(l):** a smooth positive envelope, exp of a degree-6 Chebyshev polynomial in log l, fitted by least squares
  to the log of the running maximum of |l(l+1) C_l| over one period. It takes out the damping, so what is left
  oscillates with a nearly constant amplitude. The fit is the "power law" part of the suggestion. It only
  normalises: the result passes through the samples whatever S is.
* **The spline:** on each interval, a function in span{1, l, cos(omega l), sin(omega l)} (a spline under
  imaginary tension), C2 at the nodes, with omega = 2 theta_s. It is exact for a sinusoid of the peaks'
  frequency plus a linear trend, so the "cosine with a phase" needs no fitting at all. A slowly varying phase
  and amplitude are what is left for it to interpolate. It is built and solved like a cubic spline
  (tridiagonal), and becomes the cubic spline as omega -> 0.
* **End conditions:** u'' + omega^2 u collinear over the last three nodes, which every function of the space
  satisfies; on cos + linear the spline is then exact to 3e-14 over the whole range [M, bench].
* **Where it is used:** scalar TT, EE, TE in the spectra module, and lensed TT, EE, TE in the lensing module.
  Everything else (phiphi, BB, number counts, tensors) keeps the cubic spline. It falls back to the cubic
  spline for a zero spectrum, for fewer than 8 samples, and when omega h > 2 (the equations become singular as
  omega h -> pi).

## Alternatives tried (`benchmarks/lint_template/bench.py`, on CLASS++'s own samples against an l_linstep 5 truth)

* **A spline of D_l = l(l+1) C_l** instead of C_l: no gain.
* **Global least-squares fits** of the oscillation, as the suggestion put it, then a spline of the residual.
  With frequency 2 theta_s only, the gain was 1.2x at step 30 and the result was worse at step 50. Temperature
  also oscillates at theta_s (the odd/even peak modulation from baryon loading): with both frequencies and an
  envelope, ~1.7x. Scanning the frequency or adding a chirp changed nothing. Fitting from lower l was fragile:
  a polynomial fit that is poor somewhere makes the residual worse than the spectrum itself, and the errors
  depended on where the fit started. Rejected for robustness.
* **The trigonometric spline without the envelope:** 3x at step 50 (EE 1.81 -> 0.42). The envelope takes it to
  6.5x (EE -> 0.28).

## Verification

* `unittests/tools/acoustic_spline_test.cpp`: passes through every sample; on a damped, phase-drifting
  acoustic-like spectrum sampled every 50, the largest relative error is 1.4e-4 against the cubic spline's
  3.2e-3 [M]. It also refuses a zero spectrum, samples too sparse for omega, and too few samples.
* CLASS++ against the l_linstep 5 truth, CV distance l = 2..3000, LCDM [M] (EDE, scf, Omega_k = +0.01 alike):

  | | step | unlensed TT / EE / TE | lensed TT / EE / TE |
  |---|---|---|---|
  | master (cubic) | 30 | 0.104 / 0.265 / 0.139 | 0.540 / 0.578 / 0.189 |
  | acoustic | 30 | 0.081 / 0.189 / 0.054 | 0.587 / 0.536 / 0.114 |
  | acoustic | 40 | 0.080 / 0.171 / 0.073 | 0.583 / 0.485 / 0.128 |
  | acoustic | 50 | 0.120 / 0.295 / 0.123 | 0.492 / 0.464 / 0.155 |

  Lensed TT has a floor of ~0.5 in every configuration, set by something other than this interpolation.
* What is left is concentrated at l = 300–500 [M]: the end of the log-spaced sampling, around the first EE peaks,
  where the phase departs from 2 theta_s l. Finer sampling there alone would remove most of it.
* Cosmological bias against Planck PR4 (CamSpec, lollipop, ACT DR6 lensing) + DESI DR2
  (`benchmarks/codecmp`, LCDM, 160 draws, the precision of the 2026-10 retune with `l_max_scalars` 3000)
  [M]:

  | | step 30 | step 40 | step 50 | step 60 |
  |---|---|---|---|---|
  | cubic | 0.175 ± 0.003 | | 0.363 ± 0.003 | |
  | acoustic | 0.163 ± 0.003 | 0.172 ± 0.003 | 0.183 ± 0.003 | 0.207 ± 0.003 |

  The acoustic spline at step 50 is as good as the cubic at 30; the cubic at 50 doubles the bias.
* Time, 1 thread, the same 20 draws, fastest of 3 interleaved rounds [M]: cubic step 30 0.848 s, acoustic step 30
  0.845 s (the spline costs nothing measurable), acoustic step 40 0.782 s (-7.8%), step 50 0.750 s (-11.6%).
  Which step to use is left to the precision defaults.
* `test_scalar_tensor_matches_hi_class_brans_dicke` failed: EE at l = 810, 1.12e-3 against atol 1e-3. Its
  pins were hi_class at its default step 40 (cubic). Both codes converged in l (step 5) [M]: they agree to 2e-4
  except EE at 810, 1.1e-3, where the BD/LCDM ratio is steepest. That is a difference between the codes, at the
  level the scalar-tensor spec measured for EE (1.1e-3 rms). Master passed only because its cubic spline erred
  by -3.6e-4 there, towards hi_class. The pins are now hi_class converged in l, and the C_l tolerance 1.5e-3.
  The rebuilt hi_class setup reproduces the old pins to 7e-5, and to 1.8e-4 at TT l = 2000, the band edge [M].
* `TEST_LEVEL=2` against master's classyref: 1512/1512 [M]. The other Python tests fail the same 18 golden
  comparisons as master on macOS.
