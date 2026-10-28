# Halofit and the dark energy: w(a) summed over the species, and a Lambda that was dropped (2026-10-27)

**Status 2026-10-27:** implemented on branch `nonlinear/halofit-w-de`. Evidence tags: [M] measured here, [C] stated
by the code, [P] from the literature.

## Problem

The three-code comparison after #484 showed the C_L^phiphi of the CLASS++ and CAMB references disagreeing by 3.5e-3
(L = 40-763) for w0waCDM, against 3.5e-4 for LCDM [M]. With halofit off the two agree to 3e-4 [M]. The halofit of
CLASS++ differed from class_public's and CAMB's in two details, both corrected in class_public by Lesgourgues in
946b9a15 (2020-06-10, v2.9.4, "corrected two details in halofit formulas") after CLASS++'s non-linear module had
split off (2020-04-08) [C]:

* **w.** Takahashi et al. (2012) fitted constant w, so w enters the formula through Omega_v (1 + w). CLASS++ used the
  fluid's w today at every redshift. class_public uses w(a) at each redshift, and so does CAMB (`halofit.f90`:
  `w_hf + wa_hf (1 - a)`) [C]. Mead et al. (2016, sec. 3) ran five w(a)CDM N-body simulations (w0 = -1,
  wa = -1 ... 0.7) and found w(a) gave the more accurate halofit [P].
* **Neutrinos.** CLASS++ multiplied the halo term by (1 + fnu (0.977 - 18.015 (Omega0_m - 0.3))), a correction of
  Bird et al. that CAMB dropped and class_public removed with the comment saying so [C].

At (w0, wa) = (-0.45, -1.6), near DESI DR2's best fit, the CLASS++ halofit boost P_nl/P_lin differed from class_public's
by 7.5% at z = 0.5 and 5.2% at z = 1; at (-0.9, -0.3) by 1.1% [M].

Looking at what halofit should read for w also turned up three more problems:

* **Omega_v** was 1 - Omega_m - Omega_r (class_public's choice), which includes the curvature, so a curved model was
  interpolated as flat (Takahashi's frac = Omega_v / (1 - Omega_m) is meant to interpolate open against flat) [C].
* **Omega_m of a scalar field.** The background split every EnergyType::Other species into radiation (3p) and matter
  (rho - 3p). That suits massive neutrinos, and class_public does the same for its scalar field ("not very
  meaningful for oscillatory models"). For a quintessence field near w = -1 it puts 4 rho into rho_m and -3 rho into
  rho_r: Omega_m(z = 0) was 0.486 instead of 0.312 for Omega_scf = 0.1 at w = -0.245 [M]. Halofit read that
  Omega_m. The axion EDE fluid (also Other) did the same before it thaws, around matter-radiation equality.
* **Pk_equal** interpolated the equivalent model's w and Omega_m, then overwrote Omega_m and Omega_v with the real
  model's [C]. class_public devel does the same (`external/Halofit/halofit.c` lines 207-208), against its own comment
  ("Omega_m(tau_i) have been pre-computed ... we just need to interpolate").

And the one checked in passing, a plain bug: **Lambda existed only as the closure species** (`LambdaSpecies::CreateAll`
returned nothing without a closure value). Omega_Lambda = 0.3 with a fluid gave a fluid of 0.69 and no Lambda, without
a warning; class_public gives Lambda 0.3 and a fluid of 0.39 [M].

## Approach

* **Dark energy is what has negative pressure** (Thomas's rule: simple, and no new species property). The background
  sums rho and p over every species of type DarkEnergy and over every other species while its pressure is negative,
  and stores Omega_de = sum rho / rho_crit and w_de = sum p / sum rho (-1 when there is none) in the long-format
  background vector, next to Omega_m. A species of type Other is split into radiation and matter only while p >= 0.
* **Halofit** reads Omega_m, Omega_v = Omega_de and w = w_de at tau, and applies (1 + 0.977 fnu).
* **Pk_equal** uses the equivalent model's w, Omega_m and Omega_v = 1 - Omega_m, as described in 0810.0190 and
  1601.07230.
* **Lambda** given in the input is built in ConstructSpecies's first pass like any fixed component; only an absent
  Omega_Lambda makes it the closure species.

## Alternatives considered

* **A species property (`IsDarkEnergy()`)**, with the scalar field's role set by its potential (quintessence potentials
  count, the axion potential does not). Rejected for simplicity: halofit is a best-effort fit in these models anyway.
  The pressure rule keeps an oscillating axion field that dominates the dark matter from being badly wrong, and it
  needs nothing new on the species.
* **Keeping class_public's Omega_v** (1 - Omega_m - Omega_r): identical for flat models, wrong for curved ones.
* **A continuous split for p < 0** (Lambda-like -p plus matter-like rho + p): continuous where a field's pressure changes
  sign, but it counts a w = -0.25 field as 75% matter and loses w.

## Measurements [M]

* Halofit boost against class_public devel (same formula now), largest |ratio - 1| over k = 0.01-5/Mpc at z = 0, 0.5,
  1, 2: LCDM + 0.06 eV 1.5e-3 / 1.3e-3 / 1.1e-3 / 7e-4 -> 6.0e-4 / 3.8e-4 / 1.8e-4 / 5.1e-4; w0wa (-0.45, -1.6)
  1.5e-3 / 7.5e-2 / 5.2e-2 / 9.3e-3 -> 5.7e-4 / 3.8e-4 / 2.0e-4 / 5.2e-4; constant w = -0.9 as LCDM. What is left is
  common to every model and is not w.
* C_L^phiphi on two codecmp w0wa draws at reference precision, against CAMB at AccuracyBoost 2: 3.5e-3 and 3.2e-3 ->
  3.6e-4 and 3.5e-4 (L = 40-763), the LCDM level. The halofit boost of C_L^phiphi at L = 700: CLASS++ 1.1255 -> 1.1218,
  CAMB 1.1219. class_public devel at its defaults is 1.1193, low by its halofit k range (5/Mpc), not by w.
* `test-dark-energy-budget`: Omega_de and w_de against the per-species background columns for LCDM, a CPL fluid,
  Lambda 0.3 with a fluid closure, and a quintessence field with Lambda (400 of 4605 rows with p < 0): exact to rounding,
  Omega_m + Omega_r + Omega_de = 1 to 4e-16, Omega_m(0) = 0.3099 = cdm + b. Fails on master (no Lambda, Omega_m 0.486).
* `test_halofit_dark_energy_follows_w_of_a` pins the boost at (-0.45, -1.6), z = 0-2, k = 0.2-2/Mpc, to class_public
  devel at 1.5e-3: this code is within 5.3e-4, master 4.2% off at z = 0.5.
* TEST_LEVEL=2 against master's classyref: 9 of 1512 scenarios differ, all with the scalar field and halofit. There
  halofit's Omega_m goes from 0.486 to 0.312 (z = 0) and 0.864 to 0.674 (z = 1), and the boost moves by up to 1.5%
  (z = 0), 1.0% (z = 1), 2e-4 (z = 2). Massive-neutrino, curved and fluid scenarios stay within the comparison's
  tolerances. The other Python tests pass, and the golden comparisons fail the same 18 as master on macOS.
* Equality (a_eq, H_eq) feeds only the non-linear k scaling at nonlinear_module.cpp:1810 and classy's `k_eq()`; for the
  axion EDE fluid it now comes from a budget without EDE's -3 rho in radiation [C].

## Not done

* HMcode 2020 integrates its own growth in a cosmology of matter, Lambda and the fluid, and refuses other smooth
  components (NonClusteringOmega0) [C]; its dark energy is still the fluid's.
* Closing the budget with a scalar field (Omega_scf < 0) returned an unclosed cosmology for the explanatory.ini
  potential, before and after this change (a tracker cannot close it; class_public fails there differently) [M].
