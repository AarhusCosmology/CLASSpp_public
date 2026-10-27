# The k range of the scalar CMB C_l's follows the multipoles asked for (2026-10-27)

**Status 2026-10-27:** the k range: merged (#480). The headroom: `delta_l_max` = 1000, Thomas's choice ("it
makes sense to be precise in the range asked for"; otherwise people ask for a higher l_max than they need,
which costs more), branch `perturbations/lensing-headroom-1000`. Both change the default spectra. Evidence
tags: [M] measured here,
[C] stated by the code.

## Problem

The scalar CMB transfer functions were computed up to k_max tau0 = 1.8 l_scalar_max [C], where
l_scalar_max = l_max_scalars + delta_l_max when lensed C_l's are asked for. Two defects follow
(`benchmarks/camb_anatomy/FINDINGS.md`, defects a and b):

* **The lensing headroom drags the k range up.** The multipoles above l_max_scalars only reach the lensed
  C_l's through the lensing convolution, which smooths them. Resolving them in k cost 7% of a lensed run at
  l_max_scalars = 3000 (k_max tau0 6300 against 5400) for no change in the lensed spectra [M]. Widening the
  headroom, the fix for the top of the lensed range (TT 1% off at l = 2750-3000), cost more time again.
* **Small l_max.** k_max tau0 = 1.8 l_max is 360 at l_max_scalars = 200: unlensed TT came out 11% low at
  l = 200 [M], and up to 12% and 4.9% low near l_max for l_max_scalars = 500 and 1000 [M, FINDINGS].

## Approach

k_max tau0 = max(`k_max_tau0_over_l_max` l, l + `k_max_tau0_above_l_max`), with l = l_scalar_max minus the
lensing headroom (`perturbs::l_scalar_max_lensing_margin`, set by the input module), and
`k_max_tau0_above_l_max` = 2500.

Why this form [M, `benchmarks/topband/lmax_scan.py`, against references at precise settings with
k_max tau0 = max(3 l, 6000), headroom 2000, l_linstep 10]: 1.8 l alone (no headroom) is enough at
l_max_scalars = 3000, but not at 2000. There the transfer functions still carry power above k tau0 = 1.8 l,
and the lensed TT at the top of the range had twice the error. Above l ~ 3000 the damping tail removes it.
l + 2000 is not enough at 3000; l + 2500 matches 1.8 l there. The second term also gives small l_max a floor
(2700 at l_max_scalars = 200).

Tensors and vectors keep k_max tau0 = `k_max_tau0_over_l_max` l_max.

## Measurements [M]

Lensed, master precision, 2 codecmp draws, CV distance to the reference (TT body l = 2..l_max-250 / top band
l_max-250..l_max):

| l_max | rule, headroom | TT body | TT top | EE top | time (1 thread, 12 draws x 3 rounds) |
|---|---|---|---|---|---|
| 2000 | old, 500 | 0.45 | 0.93 | 0.51 | 0.584 s |
| 2000 | new, 500 | 0.45 | 0.93 | 0.51 | -0.5% |
| 2000 | new, 1000 | 0.44 | 0.14 | 0.35 | +3.5% |
| 2500 | old, 500 | 0.72 | 2.58 | 1.33 | 0.666 s |
| 2500 | new, 500 | 0.72 | 2.60 | 1.33 | -2.9% |
| 2500 | new, 1000 | 0.57 | 0.62 | 0.50 | +0.6% |
| 3000 | old, 500 | 2.03 | 6.21 | 2.52 | 0.763 s |
| 3000 | new, 500 | 2.04 | 6.24 | 2.52 | -7.1% |
| 3000 | new, 750 | 1.23 | 3.34 | 1.19 | -5.1% |
| 3000 | new, 1000 | 0.89 | 1.91 | 0.72 | -3.8% |

* b_cosmo against Planck PR4 + DESI DR2 at l_max_scalars 3000 (160 draws, `benchmarks/codecmp`): old 0.148,
  new with headroom 500/750/1000: 0.148/0.147/0.147 (+-0.004). Per spectrum against the codecmp reference, the
  wider headroom brings TT l = 30-2500 from 0.886 to 0.750 and the top band from 6.42 to 1.92; phiphi does not
  move. The lensing likelihood's S goes 0.023 -> 0.055 with headroom 1000: a cancellation undone. With the
  reference's phiphi put in, it is 0.036 at headroom 500 and 0.006 at 1000 (b_cosmo 0.143). So at 500 the
  lensed CMB's top-band error, which the likelihood's normalisation reads, partly cancelled the phiphi error
  (-4.7e-4 on average at L = 40-763); at 1000 that error is gone and the phiphi error shows.
* Unlensed, small l_max (CV distance to the reference, body / top quarter): l_max 200: TT 6.32 / 8.86 -> 0.04 /
  0.02 (+65% time, 0.21 -> 0.35 s); 500: 7.91 / 25.0 -> 0.06 / 0.67 (+48%); 1000: 3.80 / 12.2 -> 0.37 / 0.23
  (+35%). Lensed at l_max 500: TT 0.30 / 0.86 -> 0.06 / 0.11, phiphi max error 3.1e-3 -> 1.3e-3 (+15%).
* Presets (lensed, l_max 2000-3000): `precise` takes the new rule, 22-32% faster (46 -> 32 s at 3000) with the
  same spectra (TT body 0.37 -> 0.39, top 1.42 -> 1.50). `act_dr6` keeps the k range it was tuned with
  through `k_max_tau0_above_l_max` = 10840 (exact at its l_max 9500); under the new default it would move
  TT/EE by up to 1.8e-4 and phiphi by 1.4e-3 at l <= 8500. `fast` and `tuned` are being removed (Thomas,
  2026-10-27; separate PR); for the record, most of their lensed TT error came from their k range: CV 13 at
  l_max 2500, 1.6 under the default, for +22% and +70% time.

## Verification

`test_cmb_tt_at_small_l_max` (python/test_class.py): unlensed TT asked for up to l = 200 and up to l = 1000
must agree at l <= 200 to 1e-3. It agrees to 5e-5 here and fails on master (11%).

Against master's classyref, as CI runs it [M]:

* `TEST_LEVEL=2`: 1508/1512. The four failures (1161, 1167-1169) are the closed universe with mPk, tCl, lCl
  and lensing, on P(k) at k ~ 1.8/Mpc (1.2e-2 against rtol 1e-2). Both codes are 0.6-1.4% from a converged run
  there, by where the 10-per-decade P(k) nodes land above the CMB range, which now ends elsewhere: master is
  closer for the closed universe (6.6e-3 against 9.7e-3), this branch for the flat (7.1e-3 against 9.8e-3) and
  open ones. At 20 per decade both give the same P(k) to the digit, 6.6-7.1e-3 from converged (the precision of
  the high-k modes).
* `test_dcdm_dr_matches_reference`, `test_theta_s_shooting_matches_reference`: l_max_scalars 200, where master
  is 10-11% low at l = 200 against its own l_max 2000 run and this branch 5e-5. The other Python tests fail the
  same 18 golden comparisons as master on macOS. ctest passes.

## The headroom

`delta_l_max` goes from 500 to 1000 (the "new, 1000" rows above): at l_max_scalars 2000/2500/3000, +3.5%/+0.6%/
-3.8% time against the old rule, the lensed TT of the top 250 multipoles 0.93/2.58/6.21 -> 0.14/0.62/1.91 CV.
The date presets (`precision_pre_*`) are not updated here; the precision sync that follows (#473) is.
* The codecmp reference config relies on the old rule for its k range (1.8 x (3000 + 1800)): set it
  explicitly when regenerating the reference.
