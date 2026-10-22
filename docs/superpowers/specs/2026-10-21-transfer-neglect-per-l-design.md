# Transfer functions neglected by a cut that grows with l

**Status (2026-10-21).** Shipped: the four scalar CMB cuts are l-dependent and are the
defaults below. Old behaviour: `pre/precision_pre_2026_10_21.pre`, and every measured
preset pins the constant cuts it was scored with. Not done: the lensing-source
extrapolation at the end of this note, and any change to the k grids (rejected below).

Evidence tags: **[M]** measured for this note, **[C]** read from the code, **[P]** an
earlier campaign, not re-measured, **[?]** no evidence.

## Problem

The transfer stage is as expensive as the perturbations: 626 of 1350 ms for one
Planck-like cosmology on one thread (lensed C_l to l = 3000, l_max_scalars = 3500, halofit)
**[M]**. Almost all of it is the line-of-sight integral, and the hot spot is the Hermite
interpolation of the Bessel functions, about a quarter of the whole run **[M]**.

That interpolation already does what makes it cheap per point. It keeps the coefficients
of the current interval of the Bessel table and only evaluates the cubic while the argument
stays inside it **[C]**, at about 2.3 ns a point **[M]**. Order 4 is cubic Hermite, the
lowest order that uses the stored derivatives **[C]**. What is expensive is the count:
163 million evaluations per cosmology **[M]**, two thirds of them at l < 250.

The count comes from which (q, l) pairs are integrated at all. A scalar CMB type X is set
to zero for l < (q - delta_k_X) r_a(tau_rec) **[C]**, with one constant margin per type.
Since #397, E has delta_k = 1/Mpc, so it is never cut **[C]**. The Bessel functions are
interpolated once per (q, l) for all types of a time sampling, over the widest range any
of them needs **[C]**, so E keeps every q at every l for T0, T1 and T2 as well.

## What each type needs

From transfer tables with nothing cut (one Planck-like cosmology, l to 4000), the margin
q_c - l/r_a(tau_rec) beyond which cutting that type alone moves TT, TE and EE by less
than 1e-4 of their scale **[M]**:

| l | T0 | T1 | T2 | E |
|---|---|---|---|---|
| 2 | 0.032 | 0.007 | 0.001 | 0.002 |
| 23 | 0.098 | 0.021 | 0.051 | 0.031 |
| 143 | 0.128 | 0.020 | 0.085 | 0.070 |
| 387 | 0.167 | 0.039 | 0.126 | 0.100 |
| 1000 | 0.176 | 0.045 | 0.137 | 0.123 |
| 2500 | 0.179 | 0.040 | 0.135 | 0.148 |
| 4000 | 0.177 | 0.020 | 0.124 | 0.156 |

(1/Mpc.) At low l every type needs a small fraction of what a constant cut keeps. At high l
the constant T0 (0.15) and T2 (0.12) cuts were short, which is the TT residual of the
transfer integrals noted in `pre/act_dr6.pre` **[P]**.

## Approach

delta_k(l) = transfer_neglect_delta_k_S_X sqrt(l / (l + transfer_neglect_l_half_S_X)): one
asymptotic margin and one scale per type. l_half = 0 restores the constant cut exactly
(byte-identical with fast-math off) **[M]**. Defaults, chosen to cover the table above at
every l:

| type | delta_k (1/Mpc) | l_half |
|---|---|---|
| T0 | 0.19 | 50 |
| T1 | 0.05 | 20 |
| T2 | 0.15 | 100 |
| E | 0.165 | 400 |

## Results

Scored with the posterior metric of the three-code benchmark (`benchmarks/codecmp/`,
untracked): Planck PR4 CamSpec TTTEEE + low-l TT + LoLLiPoP EE + ACT DR6/PR4 lensing + DESI
DR2, against a converged CLASS++ reference. b_cosmo is the largest shift, in posterior sigma,
of any combination of cosmological parameters; S the spread of delta chi^2 **[M]**:

| | LCDM b_cosmo | S | S_lensing | w0wa b_cosmo | S |
|---|---|---|---|---|---|
| constant cuts | 0.175 | 0.844 | 0.106 | 0.190 | 0.944 |
| per-l cuts | 0.153 | 0.770 | 0.083 | 0.169 | 0.820 |

(160 LCDM and 112 w0wa posterior draws.) Against a run with nothing cut, unlensed to
l = 2500: TT 2.8e-4 at the median and 4.7e-4 at worst before, 3.5e-5 and 1.2e-4 now
**[M]**; `unittests/source/transfer_neglect_test.cpp` holds it there. Bessel evaluations
163M -> about 100M, transfer stage 650 -> 490 ms on one thread, total -13% **[M]**.

## Alternatives measured and rejected

- **The first fit (E l_half 400 with delta_k 0.16, T2 0.14 with l_half 200).** T2 fell short
  of the table from l ~ 400 up (TT 2.4e-4 at l = 403 against nothing cut). b_cosmo 0.159,
  15 ms cheaper **[M]**.
- **Every type at T0's cut.** The Bessel functions are shared, so this looked free; it is
  not, because T2 and E need the second derivative and widen where the costlier
  interpolation runs: +70 ms, b_cosmo 0.150 but S 0.798 **[M]**.
- **Smaller asymptotic margins** (x0.8): b_cosmo 0.239. **A slower ramp** (l_half x2):
  0.181 for no speed. The asymptotic margin is what sets the accuracy **[M]**.
- **Sparser k grids.** k_step_sub 0.08 -> 0.12 gives b_cosmo 0.308 for 7%; q_linstep
  0.45 -> 0.6 breaks (3.6) **[M]**. CAMB's source grid (`SetkValuesForSources`) is not
  sparser where it matters: 0.75/tau_rec ~ 0.0027/Mpc through the acoustic range against
  CLASS++'s 0.0035, and sparser only at low k, where modes are cheap **[C]**.

## Not validated

Curved models, isocurvature modes, and species with unusual late-time sources **[?]**. The
cut is a property of the scalar CMB integrals, not of a species, but its margins were read
off LCDM and checked on w0wa only. Vector and tensor cuts are unchanged.

## Next: the lensing source beyond the last mode

The lensing potential cannot be cut this way: it needs q up to the full grid at every l
**[M]**. Its cost at high L is the other way round. The full-Limber grid stops at the last
perturbation mode, so at L = 3000 it misses everything below z ~ 0.14, which is why
C_L^phiphi is 2.4% low there unless halofit's k range is raised to 70/Mpc (at 2.6x the
runtime) **[P]**. Continuing the lensing source past the last mode with its high-k shape,
times halofit's correction on the linear spectrum it already extrapolates (#399), is the
candidate.
