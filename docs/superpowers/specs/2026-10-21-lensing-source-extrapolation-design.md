# The CMB lensing source beyond the last perturbation mode

**Status (2026-10-21).** Implemented, stacked on #469 (per-l transfer cuts). The linear
sources are continued past the last perturbation mode up to `k_max_extrapolated` = 70/Mpc;
the full-Limber CMB lensing source uses the continuation, times the non-linear correction
there, and halofit and HMcode trust it, so `nonlinear_min_k_max` came down from 5 to
2/Mpc. `precision_pre_2026_10_21.pre` (and the two older history files) restore the
previous behaviour byte for byte. Follows `2026-10-21-transfer-neglect-per-l-design.md`.

Evidence tags: **[M]** measured for this note, **[C]** read from the code, **[P]** an
earlier campaign, not re-measured, **[?]** no evidence.

## Problem

C_L^phiphi is computed in the full-Limber scheme on a q grid that ended at the last
perturbation mode **[C]**. In the Limber integral, multipole L at wavenumber k comes from
comoving distance chi = (L + 1/2)/k, so stopping at k_max drops everything closer than
(L + 1/2)/k_max: at L = 3000 with k_max = 5/Mpc, everything below z ~ 0.14 **[C]**.
With halofit that is where the non-linear power is largest. C_L^phiphi was 0.5% low at
L ~ 900 and 2.4% low at L ~ 2850 **[P]**, and lensed BB and the ACT DR6 lensing likelihood
carry it.

A second truncation sat in halofit itself: it declares itself computable only at
redshifts where its smallest smoothing radius, set by the last computed mode, still has
sigma > 1. With modes to 5.8/Mpc that stopped it at z = 3.3 **[M]**, and the correction
dropped abruptly to 1 above, at all k.

Solving more modes fixes both at a high price: nonlinear_min_k_max = 70/Mpc costs 2.6x the
runtime **[P]**. The 5/Mpc floor itself was not cheap either: in a CAMB-matched run its
modes above 1.84/Mpc were a quarter of the perturbation stage **[M]**.

## Approach

**The continuation** (`AsymptoticExtrapolation`, `tools/asymptotic_extrapolation.cpp`).
Deep inside the horizon in matter domination, phi + psi = -3 a^2 rho_m delta_m / k^2
(Poisson), and delta_m grows as ln k with a 1/k correction. So s(k) k^p = A ln k + B + C/k,
p = 2 for the potential and 0 for delta_m. A and C are fitted by least squares to the modes
in [k_max/3, k_max), B so that the continuation passes through the last mode. Against
solved modes (LCDM + one 0.06 eV neutrino, z = 0, 1, 3, phi + psi) **[M]**:

| fitted on | k = 5 | 10 | 20 | 40 | 70 /Mpc |
|---|---|---|---|---|---|
| [0.33, 1] | -0.4% | -0.6% | -0.7% | -0.6..-0.7% | +0.5..-0.1% |
| [0.67, 2] | -0.06% | -0.13% | -0.15..-0.19% | +0.07..-0.08% | +0.6..+1.2% |
| [1.67, 5] | - | 0.00% | +0.04% | +0.2..+0.3% | +0.9..+1.5% |

The excess at 70/Mpc grows towards z = 0, as baryon pressure would **[?]**; at L = 3000,
k = 70 is chi = 43 Mpc, where the lensing kernel has almost no weight. Rejected **[M]**:
the two-parameter form A ln k + B (the old `extrap_max_scaled`), 1.0-1.5% low at
10-40/Mpc from k_max = 2; the Eisenstein-Hu transfer-function ratio, 1-3% high and growing.

**Where it is used** **[C]**:
- `extrapolation_method` defaults to the new `extrap_asymptotic` (was `extrap_max_scaled`)
  for the linear P(k) that halofit and HMcode integrate beyond the last mode. It is now
  read for every run, so that a precision file can set it.
- halofit's smallest radius and HMcode's k_NL test use `k_max_trusted_` =
  max(last mode, `k_max_extrapolated`) instead of the last mode, so halofit now reaches
  z = 7.4, as it does with modes solved to 70/Mpc **[M]**.
- the non-linear module computes P_NL, and stores sqrt(P_NL/P_L) in
  `nl_corr_density_extra_`, on its extended k list up to `k_max_extrapolated`.
- the transfer module fits the phi + psi source per initial condition and time once
  (`transfer_lensing_source_extrapolation`); the full-Limber q grid runs to
  `k_max_extrapolated` and evaluates the continuation times the correction there
  (`transfer_extrapolate_lensing_source`).
- the primordial table reaches that k too. An external table (`external_Pk`) need only
  cover the computed modes, as before; past its last point ln P is continued linearly in
  ln k, as the non-linear module does for its own tail. The bundled example stops at
  10/Mpc, and without this a lensing run with it stopped with k = 10.4 out of range **[M]**;
  a power law tabulated to 10/Mpc now gives the analytic C_l to 4e-15 **[M]**.
- `k_max_extrapolated = 0` turns all of it off.

## Results

**Against modes solved to 70/Mpc**, per model, max relative error over L = 2-3000 and, for
P_NL, over k < 0.67/Mpc at z = 0-3; "old" is the #469 code (floor 5, no continuation)
**[M]** (`benchmarks/lensing_extrap/robust.py`):

| model | C_L^phiphi new / old | lensed BB new / old | P_NL worst z new / old |
|---|---|---|---|
| LCDM, 0.06 eV | 8.9e-5 / 2.7e-2 | 2.1e-5 / 1.1e-2 | 5.2e-4 / 1.1e-2 (z = 3) |
| 3 x 0.2 eV neutrinos | 7.0e-5 / 4.5e-2 | 2.7e-5 / 1.7e-2 | 9.1e-4 / 4.0e-1 (z = 3) |
| w0wa (-0.8, -0.6) | 1.1e-4 / 2.6e-2 | 1.9e-5 / 1.1e-2 | 5.2e-4 / 1.2e-2 (z = 3) |
| Brans-Dicke, omega_BD = 50 | 9.4e-5 / 1.9e-2 | 1.0e-5 / 6.9e-3 | 1.6e-4 / 8.1e-7 (z = 3) |
| Omega_k = -0.05 | 1.2e-4 / 2.5e-2 | 1.4e-5 / 1.0e-2 | 3.7e-4 / 3.5e-3 (z = 3) |
| Omega_k = +0.05 | 3.4e-3 / 3.1e-2 | 3.2e-4 / 1.1e-2 | 6.4e-4 / 4.1e-2 (z = 3) |
| DCDM, Gamma = 100 km/s/Mpc | 2.4e-3 / 2.6e-2 | 1.8e-4 / 9.5e-3 | 3.5e-4 / 6.4e-4 (z = 3) |
| HMcode 2016 | 1.9e-4 / 3.2e-2 | 3.5e-5 / 1.0e-2 | 2.8e-4 / 3.9e-1 (z = 3) |
| 1 keV WDM, 1/4 of the DM | 5.0e-2 / 1.1e-1 | 1.7e-2 / 4.6e-2 | 5.0e-2 / 2.8e-2 (z = 3) |
| IDM-DRMD, f = 0.7 (DAO) | 7.8e-2 / 1.3e-1 | 4.8e-2 / 8.2e-2 | 5.2e-2 / 3.4e-2 (z = 1) |

- The open-universe 3.4e-3 is not the continuation: it is the same with the floor at 5 or
  10, and two runs that both solve modes to 70/Mpc differ by 3.1e-3 when their Limber q
  grids end at different q (q = sqrt(k^2 + K) moves every point). That is the documented
  ~3e-3 of q_logstep_limber = 1.1, which cancels in flat and closed A/Bs only because the
  grids coincide there **[M]**.
- **Warm dark matter** breaks the asymptotic form: the free-streaming cutoff of a keV
  relic sits at a few /Mpc, inside the fitting range. Lensing is still better than before,
  P_NL at z = 3 is not. Raising the floor past the cutoff fixes it: 1.2e-3 at 5/Mpc, 1.0e-6
  at 10/Mpc (C_L^phiphi 8.5e-3) **[M]**. The precision.h comment on nonlinear_min_k_max
  says so. The code cannot see it: refitting on [k_max/4.5, k_max/1.5] and predicting the
  last mode misses by 7.9e-4 for WDM and 4.8e-4 for LCDM **[M]**, so no self-test separates
  them.
- **IDM-DRMD** at f_idm_drmd = 0.7 has dark acoustic oscillations aliased by the k grid
  (#467); both codes are percent-level off there, and the lensing source improves.

**Unit tests.** `test-asymptotic-extrapolation` pins the fit (exact for the form, through
the last point, stride, one and two points). `test-lensing-extrapolation` holds, on the
defaults, C_L^phiphi to 3e-3 of modes solved to 70/Mpc (it is 5.0e-5) and halofit P_NL at
z <= 3 to 1e-3 (4.8e-4); without the trusted continuation the second is 0.46 **[M]**.

**Cosmology** (`benchmarks/codecmp/variant.py`, Planck PR4 CamSpec + low-l + ACT DR6
lensing + DESI DR2, against the converged CLASS++ reference) **[M]**:

| | b_cosmo | S | S_camspec | S_lens |
|---|---|---|---|---|
| LCDM (160 draws), #469 | 0.153 | 0.770 | 0.768 | 0.083 |
| LCDM, continuation | 0.148 | 0.783 | 0.781 | 0.055 |
| w0wa (112 draws), #469 | 0.169 | 0.820 | 0.809 | 0.107 |
| w0wa, continuation | 0.158 | 0.821 | 0.815 | 0.070 |

S_camspec rising is what solving the modes does too: before #469, nonlinear_min_k_max = 70
took b_cosmo 0.175 -> 0.169, S_camspec 0.830 -> 0.838, S_lens 0.106 -> 0.078 **[M]**. The
truncated lensing had cancelled part of another CamSpec error.

**Runtime** (`class_profiled`, best of 3, M-series laptop) **[M]**:

| | 1 thread | 11 threads |
|---|---|---|
| fiducial (l_max 3500, halofit), #469 | 1165 ms | 262 ms |
| fiducial, continuation | 1084 ms (-7%) | 254 ms (-3%) |
| CAMB-matched (l_max 2500, halofit + mPk), #469 | 905 ms | 203 ms |
| CAMB-matched, continuation | 821 ms (-9%) | 197 ms (-3%) |

Perturbations lose 90-100 ms on one thread (17 ms on 11); the non-linear stage gains
11-13 ms because halofit now runs to z = 7.4 instead of 3.3, and that stage is serial.
Parallelising it over time would recover that at 11 threads **[?]**.

**Against master's classyref** (`COMPARE_OUTPUT_REF=1 TEST_LEVEL=2`, what the PR CI runs)
**[M]**: 24 of 1512 scenarios fail, all of them `Omega_scf = 0.1` with halofit and `lCl`, all
on C_L^phiphi. There master is the one that is off: against modes solved to 70/Mpc it is
6.6e-3 low at L = 300, 1.1e-2 at 500 and 4.5e-2 at 2000; this branch is within 3.1e-5,
2.5e-4 and 2.5e-3 (with q_logstep_limber = 1.02 on both sides, at most 1.5e-3, so the rest is
the q grid). LCDM moves as much at high L but passes, because the comparison's absolute
tolerance of 1e-20 exceeds C_L^phiphi above L ~ 1000, so only low L is really compared.
Master with modes solved to 70/Mpc agrees with this branch's solved run to 7e-6.

**Presets.** The history files `precision_pre_2026_10_21.pre`, `_10_09` and `_08_21` set
nonlinear_min_k_max = 5, k_max_extrapolated = 0 and extrapolation_method = 3, and
reproduce the #469 code byte for byte with fast-math off (halofit, HMcode and linear) **[M]**.
The measured presets are not pinned: without a non-linear method the continuation moves
lensed TT and EE by at most 3.5e-5 and 6.1e-6 (fast, tuned, precise, act_dr6), C_L^phiphi
by 4e-3 and BB by 8e-4; act_dr6.pre with halofit moves by 1e-7, since it solves modes to
70/Mpc **[M]**.

## Open

- The non-linear stage is serial and now 26 ms; parallel over time it would shrink at
  11 threads **[?]**.
- A species with a small-scale cutoff (warm dark matter) could raise the floor itself.
  Not done: the free-streaming scale of a relic is not a quantity the species exposes
  today **[C]**.
