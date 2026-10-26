# Transfer: Bessel interpolation located once per wavenumber, node derivatives tabulated (2026-10-27)

**Status 2026-10-27:** implemented as described below (branch `transfer/los-hermite-kernel`, on top of #474).
Evidence tags: [M] measured here, [C] stated by the code, [?] not measured.

## Problem

The transfer stage was the largest gap to CAMB: 412 ms against 133 ms on one thread for LCDM lensed to
l = 3000, for the same number of (l, q) pairs and time points (`benchmarks/camb_anatomy/FINDINGS.md`).
That is 4.7 ns per point against 1.7 [M]. About 60% of it was the cubic Hermite interpolation of the flat
Bessel functions, `hyperspherical_Hermite_interpolation<4, ...>`, called once per (l, q, time sampling)
with the arguments x = q (tau0 - tau). For each argument it:

* found the interval: a division for every point where it jumped, and a division for z = (x - x_i)/dx at
  every point;
* evaluated Phi'' and Phi''' at the interval's nodes from the radial equation (divisions by sin^2) each
  time it entered a new interval;
* built the coefficients of up to three cubics (Phi, Phi', Phi'') from those node values.

Over 80.6M points, 36% stayed in the previous interval, 20% moved to the next one and 44% jumped [M]. The
location and the Hermite basis depend on x only, yet were redone for each of the ~155 multipoles at the
same x. The node derivatives depend on (l, node) only, yet were redone for every q.

`transfer_select_method` also scanned backwards from the last time point, for every (l, q, type), to the
Bessel cut, over trailing zero sources and to the late-source cut.

## Approach

All of it gives the same result up to rounding.

* **Node tables.** `HyperInterpStruct::OnFixedSpacing`, which builds the table kept by `BesselTableCache`
  (#474), also stores Phi'' and Phi''' at every node (`d2phi`, `d3phi`), with the expressions the
  interpolation evaluated on the fly. Built once per cached table.
* **Locate once.** `hyperspherical_Hermite_locate` finds, for every argument, the left node and the four
  cubic Hermite weights (value and derivative weights at both nodes). The transfer module calls it once
  per (q, time sampling), where it builds `chi_full_reverse`. Each multipole uses a suffix of those
  arguments, cut by its Bessel cut.
* **Evaluate.** `hyperspherical_Hermite_evaluate` gives each of Phi, Phi', Phi'' as four multiply-adds of
  node values with those weights. These are the same cubic polynomials as before; only the rounding
  differs. The weights are separate arrays rather than one struct per point, so they stream; that is worth
  a few percent.
* **The flat `rescale_function`** stays at the ones set in `transfer_workspace_init`. It was refilled for
  every (l, q).
* **`transfer_select_method`** finds the Bessel cut and the late-source cut by binary search on the
  decreasing `tau0_minus_tau`. It finds the last nonzero source from a table, built once per (q, type),
  of the last nonzero index at or before each time. The source can have interior zeros, so "the last
  nonzero overall" would not be equivalent; the table is. The indices are the same as the scans'.
* **Curved models** keep the old path. There the arguments depend on l through the rescaling, and the
  per-q tables are not built by `OnFixedSpacing`.

## Cost

The node tables double the cached table: about +20 MB at l_max_scalars = 3000 (155 multipoles x ~7–8 thousand
nodes x 2 arrays) [?, estimated from the sizes, not measured]. They are built once per classy object, and once per run for the
executable.

## Alternatives not taken

* **CAMB's tables** [C]: four cubic-spline coefficients per (l, interval). CAMB needs only j_l, because it
  integrates the temperature source by parts. CLASS needs Phi, Phi' and Phi'' (types t0, t1, t2, E), which
  would be 12 coefficients per (l, interval), three times this table.
* **Integrating by parts as CAMB does**, so that one j_l serves every type. That changes the sources (it
  needs their time derivatives), a physics change for another time.
* **Node derivatives on the fly**, without tables but with the location hoisted. That is about as many
  operations per point as before, because the node work is then done at every point.
* **Phi'' from the radial equation at the argument**, from the interpolated Phi and Phi'. No third
  derivative needed, but it changes the results at the interpolation-error level.
* **Fusing interpolation, radial function and integration** into one pass over the points, or evaluating
  several multipoles per point. Possible next steps, at the cost of restructuring
  `transfer_compute_for_each_q`.

## Verification

* `unittests/tools/hermite_evaluate_test.cpp` compares, for every multipole, Phi, Phi' and Phi'' with
  `hyperspherical_Hermite_interpolation<4, true, true, true>` on the same table, at arguments inside and
  outside the table, on nodes and at its last node, and through the suffix (`first`) path. The largest
  difference is 2.6e-12 of the largest value [M]. Near x_min the radial equation gives Phi''' as a
  difference of terms ~1e4 larger, and the two codes round it differently. Tolerance 1e-10. With the sign
  of one weight flipped, the test fails [M].
* Against master (92ab2bb9, which has #474) [M], max |dC_l/C_l|: lensed TT 7.6e-15, EE 1.2e-15, BB 1e-13
  for LCDM; Omega_k = +-0.01 bit-identical (old path); tensors and number counts <= 3e-15.
* Timing [M], 1 thread, LCDM lensed to 3000. `codecmp/stages.py` (median of 4 draws, best of 3, a fresh
  Class each time, so the table is rebuilt each time): transfer 411 -> 300 ms (-27%), whole run
  0.905 -> 0.797 s (-12%). One classy object over 16 draws, the table cached: 1 thread 890 -> 774 ms
  (-13%), 11 threads 172 -> 157 ms (-9.5%).
  - By step, transfer stage: locate + node tables 416 -> 317 ms; lookups in `transfer_select_method`
    -> 311 ms; separate weight arrays and the rescale fill -> 300 ms.
* ctest 60/60, TEST_LEVEL=1 84/84, TEST_LEVEL=2 as on master. The python non-scenario tests fail the same 18
  golden comparisons on master and on this branch, on this Mac.

Where the rest goes [M, timers in a scratch build, taken before the last step; the timers add ~10 ms]: one
thread spends 164 ms interpolating (about 150 now), 57 ms in radial functions and integrals, 20 ms in method
selection, 22 ms in per-q set-up and 7 ms in Limber.
