# Restructure of the transfer module

**Status (2026-10-08).** Implemented on branch `transfer-restructure`, as
designed below, except for these points.

* **The eight partial sums were not ported.** Under `-ffast-math` the compiler
  already vectorises the convolution: 0.440 s vs 0.444 s (CMB) and 5.05 s vs
  5.04 s (number counts), one thread. **[M]**
* **The vector `cotKgen` fix has no separate commit.** It came with the shared
  Bessel interpolation. Vector modes are not implemented in the CLASS++
  perturbation module, so it changes no output. **[C]**
* **The source reordering suggested by Nils is not needed.** Source
  interpolation is 2.1% of the CMB transfer time and below 1% for number counts
  (one thread, `sample`). **[M]**
* **Defaults of the #267 sampling: class_public's, kept.** With factor 1 and
  fraction 0.001, C_L^φφ is within 2.1e-4 of a converged reference at L ≤ 400.
  Factor 0.5 and fraction 1e-4 give 4e-5 at L ≤ 40 for +5% transfer time. That
  is not worth it: the full-Limber q grid (`q_logstep_limber = 1.1`) already
  limits φφ to ~3e-3 at L > 40. **[M]**
* **Commented-out debugging code was also removed,** on request.

Results **[M]** (13 harness cases plus references; details in the commit
messages):

* **Output-neutral steps:** CMB Δ_l(q) ≤ 7e-14, C_l ≤ 1e-12 against the
  previous step. One pre-existing fragility shows up instead: the g4/g5 windows
  of a bin touching z = 0, or of Dirac bins, count sources strictly behind the
  lens. A 1-ULP shift of a grid (here an FMA contraction of `mean + width·cut`
  under `-ffast-math`) moves the points where a lens time coincides with a
  source time by one step. That is 22 of 3480 points for a gaussian bin at
  0.5 ± 0.1, and C_l move by ≤ 6e-5. A smooth integrated window would remove
  this; not done here.
* **#267 port.**
  * C_L^φφ at L ≤ 40: from 2–7.5% off on master to ≤ 2.1e-4.
  * Lensed BB: 1.5e-3 → 1.7e-4. Lensed EE: 5.6e-4 → 2.0e-4.
  * `output = lCl` alone was 37–94% wrong on master (sources started at
    z_max_pk + 1); now ≤ 6.5e-4.
  * The two converged references (old scheme at step 0.0025, new scheme at step
    0.0025 with factor 0.25) agree to ≤ 3.8e-5.
* **Anti-aliasing.**
  * Reproduced the class_public low-z case exactly: Δ_l = 15.7 instead of ~0
    at q = 2.98/Mpc.
  * C_l: 1.7e-4 → 8.5e-6 against 5× denser references.
  * Other bins move by ≤ 6e-6. CMB bit-identical.
* **`selection = tabulated`:** reproduces analytic gaussian bins to a median of
  1.6e-6 (worst 2.5e-4). Per-bin biases scale C_l exactly.
* **Line count:** the transfer module went from 4615 to 4023 lines, new features
  included.

## Origin

In class_public, Julien Lesgourgues audited `transfer.c` (note of 17.09.2026,
"Audit of the transfer module", v3.4.0) and Nils Schöneberg added quick comments
by e-mail (29.09.2026). Two class_public PRs followed:

* **lesgourg/class#267** (Julien, merged into `devel` on 2026-10-04): a
  logarithmic late-time sampling of the CMB lensing source.
* **lesgourg/class#268** (stacked on #267): a restructure along the lines of the
  audit, an anti-aliasing fix for number counts, `selection = tabulated`, and two
  performance fixes in the hottest loop.

CLASS++'s `source/transfer_module.cpp` descends from the pre-#267 `transfer.c`. It
has the same functions, the same macros and the same per-wavenumber logic. The
C-to-C++ modernisation changed the memory management, the error handling and
the threading, but not the structure. **[C]** This spec ports both PRs, adapted
to CLASS++. Julien's audit is the reference; Nils's comments were quick remarks
and are weighed as such.

## Problems in CLASS++ master (532e52af)

1. **The properties of a transfer type are scattered.** What distinguishes a
   type X is its perturbation source Y(X), its redshift bin, its radial function,
   its window W(τ), its power of k, its Limber threshold and its neglect rules.
   These live in four macros (`_integrated_ncl_`, `_nonintegrated_ncl_`,
   `_get_bin_..._`) and in chains of `if`s in seven functions:
   `transfer_get_source_correspondence`, `transfer_source_tau_size`,
   `transfer_use_limber`, `transfer_can_be_neglected`,
   `transfer_late_source_can_be_neglected`, `transfer_select_radial_function`
   and `transfer_get_l_list`. Adding a type means editing all of them. **[C]**
2. **Work that does not depend on q is redone for every q.** For each (q,
   initial condition, type), `transfer_sources` rebuilds the time sampling:
   `transfer_source_tau_size` calls `transfer_selection_times`, which calls
   `background_tau_of_z` three times. It also recomputes the trapezoidal weights
   and resamples the source with one bisection per sampled time. The windows are
   precomputed once, but their sampling is rebuilt every time. **[C]**
3. **The Bessel functions are interpolated once per type and initial
   condition,** although all types of a bin share one time sampling. A bin with
   all number-count contributions has seven local types (density, RSD, d0, d1,
   g1, g2, g3) and four integrated ones. **[C]**
4. **CMB lensing is a special case** inside `transfer_sources`, with its kernel
   computed inline, instead of an ordinary type with a window of support
   [τ_rec, τ0], as the audit asks. **[C]**
5. **C_L^φφ is under-resolved at low L.** Below `l_switch_limber` (40 in CLASS++)
   the line-of-sight integral uses the perturbation time sampling. That sampling
   does not resolve the 1/χ lensing kernel near τ0. At default precision C_L^φφ
   was measured 3–8% high at L < 40 against CAMB 2 and against a CLASS++ run with
   `perturbations_sampling_stepsize = 0.01`, which agree to 1e-3. **[P]**
   (measured 2026-10-08 in an uncommitted harness; re-measured here.) #267 fixes
   this in class_public.
6. **Number counts can alias.** For evenly spaced samplings, the number of
   points resolves the Bessel oscillations up to k ≈ l_limber/(τ0 − τ_mean). The
   line-of-sight integral is nevertheless computed up to the largest q of the
   module. Near q ≈ 2π/Δτ the sampled Bessel function stops oscillating, and the
   trapezoidal sum gives large spurious values. #268 shows this at default
   precision in class_public `devel` for low-redshift thin bins (C_l off by
   6.9e-4). **[P]** The same code is in CLASS++; to be reproduced here.
7. **The vector-mode radial functions T2 and E index `cotKgen` forward while
   `Phi` is reversed.** class_public fixed this typo in #256; CLASS++ still has
   it. **[C]**
8. **Dead code:** `transfer_limber2` (never called), `transfer_get_q_list_v1`
   and `transfer_init_HIS_from_bessel` (declared, never defined), and
   commented-out debugging blocks. **[C]**
9. **Users cannot give the dN/dz of each bin.** This was requested repeatedly
   (Julien's audit, Nils, the class_public branch `class_window_read`).

## Design

The approach is data-driven, not polymorphic. Most of the variation between
types is per type, not per family. The eleven number-count and galaxy-lensing
types differ in source, radial function, power of k, Limber threshold and window
formula. A class hierarchy would only replace two small switches (window kind
and sampling kind), at the cost of a base class, four derived classes and a
factory. All the per-type formulas would remain. We estimated it at 150–300
lines longer, with no gain in clarity.

### Type properties, defined once

One function fills a record for each (mode, type):

* the index of the perturbation source and the redshift bin (−1 for CMB types);
* the radial function (Bessel function or a combination of its derivatives);
* the window kind: none, CMB lensing, local selection, integrated selection;
* the factor (k/k_pivot)^n: n = −2 for d0, −1 for d1, +1 for g5, `lcmb_tilt`
  for CMB lensing;
* the Limber threshold `l_limber` (∞ if never used);
* `delta_k_neglect`, for the rule that sets Δ_l(q) = 0 when
  l < (q − δk) r_a(τ_rec) (∞ if never applied);
* whether late sources may be cut for l > `transfer_neglect_late_source`;
* whether the type is also computed on the full-Limber q grid (only the scalar
  CMB lensing potential);
* the index of its time sampling, and its window W(τ) tabulated on it.

Every per-type decision reads this record. The macros and the `if` chains go.

### Time samplings, defined once

There are three kinds:

* **from perturbations:** the perturbation times, for CMB T, E and B;
* **uniform:** evenly spaced times, for number counts and galaxy lensing. The
  support is that of the bin for local types, and [τ_min, τ0] for integrated
  types. The number of points is unchanged from today's formula;
* **logarithmic:** for CMB lensing (#267). It keeps the perturbation times after
  τ_rec as long as their step is smaller than (τ0 − τ)·s. Then it continues with
  Δτ = (τ0 − τ)·s up to (1 − f)·τ0, and ends at τ0. The step is
  s = `cmb_lensing_sampling_factor` / min(`l_switch_limber`, `l_scalar_max`),
  and f = `cmb_lensing_sampling_fraction`. CLASS++'s `l_switch_limber` (40) is
  larger than class_public's (30), so the defaults must be checked here rather
  than copied.

Types whose samplings are identical share one. For each sampling we precompute,
before the loop over q:

* the times;
* the trapezoidal weights;
* for each time, the bracketing perturbation time and the linear-interpolation
  weight, so that resampling a source needs no search;
* the largest q at which it resolves the Bessel functions.

The logarithmic grid is built in one pass with `push_back`. This is the only
place where C counted the points first and then filled them, and the C++
container answers Nils's "growing buffer" remark there. Everywhere else the
sizes are known before the loop over q.

When φφ is requested, the perturbation source sampling starts at
(1 − 10ε)·τ_rec. This guarantees a perturbation time below τ_rec for the
logarithmic grid to interpolate from (#267's change in `perturbations.c`).

### Windows, defined once

W(τ) is tabulated for every type on its sampling:

* 1 for CMB types;
* the lensing kernel `sin_K(τ − τ_rec) / (sin_K(τ0 − τ) sin_K(τ0 − τ_rec))`
  times `lcmb_rescale` for CMB lensing;
* the selection function times a function of the background for local types;
* the selection function convolved with a kernel for integrated types.

The transfer source of every type is then S(k(q), τ) × W(τ) × (k/k_pivot)^n.

### Per wavenumber

1. Interpolate each needed perturbation source once, at k(q).
2. For each time sampling, build the transfer sources of all its types and
   initial conditions: resample with the precomputed weights, then multiply by
   W and the k factor.
3. For each l, choose a method for each (initial condition, type): zero, Limber,
   or a line-of-sight integral with its last contributing time. Interpolate only
   the needed Φ, Φ′ and Φ″, once, over the union of the needed ranges. Then
   build each radial function from them and integrate.

The full-Limber pass runs the same code and skips the types without the
full-Limber flag.

### Anti-aliasing

For each uniform sampling with more than one point,

    q_max_sampled = (N − 1) / (τ_max − τ_min) / selection_sampling_bessel_min.

This is the largest q at which the Bessel functions have at least
`selection_sampling_bessel_min` points per radian (new precision parameter,
default 1). Above it the Limber approximation is used, which there gives ≈ 0, as
the true transfer function does. #268 found identical results for thresholds
between Nyquist (1/π) and 1. **[P]**

Two samplings are exempt. The logarithmic one is built to resolve the Bessel
functions down to τ0. The perturbation one is not chosen for the line-of-sight
integral at all. In class_public, applying the cut to the logarithmic grid by
accident changed C_l^φφ. #268 caught this semantic conflict, which git did not
report, while rebasing.

### `selection = tabulated`

A file has z in the first column and dN/dz of each bin in the next columns.
Lines starting with `#` are ignored. The number of bins is the number of columns
minus one, and the bins are ordered by growing mean redshift.

* The input module reads the file once. It derives the number of bins, their
  mean z (used for the Limber thresholds and k_max) and their non-zero z range
  (used for `z_max_pk` and the time sampling).
* The transfer module reads the file again and splines it.
* `selection_bias` and `selection_magnification_bias` work as for other bins.
* A global `dNdz_selection` is refused with an explicit error, since it would
  multiply a dN/dz that is already complete.

Different dN/dz for different tracers is out of scope. nCl and sCl would need
separate sets of bins in the perturbation and harmonic modules.

### Not ported from #268

* **Eight independent partial sums in `array_trapezoidal_convolution`.** In C,
  the single sum formed a dependent chain of additions. CLASS++ builds with
  `-ffast-math`, which lets the compiler reassociate and vectorise the
  reduction. This is ported only if it is measured to pay.
* **The crash with `output = pCl,lCl` plus tensors, and the out-of-bounds read
  on `q_limber`.** CLASS++'s full-Limber pass already skips non-scalar modes and
  reads `q_limber_`/`k_limber_`. **[C]** To be confirmed by running the case.

## Alternatives rejected

* **Polymorphic transfer types:** see the start of Design.
* **Caching only** (keep the macros, add caches): smaller diff, but it leaves
  problems 1 and 4 in place.
* **A q-dependent time sampling** (the audit's question about
  `transfer_source_tau_size`). #268's sampling study found the fixed sampling
  accurate, provided the line-of-sight integral is not used above
  `q_max_sampled`. **[P]**
* **Reordering the perturbation sources to put k first** (Nils). This touches
  the perturbation and nonlinear modules and classy. In class_public, source
  interpolation was 4.6% of the CMB transfer time after the restructure. **[P]**
  To be measured here; changed only if it is large.

## Verification

* **Output-neutral steps** (properties, samplings, windows, shared Bessel
  functions, Limber bisection): A/B against master on raw Δ_l(q) tables and C_l
  files, built in one build directory.
  * About 13 cases: CMB with and without full Limber, tensors, vectors; ad+cdi
    with halofit; number counts with gaussian, top-hat, Dirac and tabulated
    bins, global dN/dz from file and analytic; open and closed; `pCl,lCl` with
    tensors.
  * Expected: rounding-level differences. Anything above ~1e-6 (scale-relative)
    is investigated.
* **Steps that change the output,** each in its own commit:
  * vector `cotKgen`: vector modes only;
  * #267 sampling: C_L^φφ against a converged reference (CLASS++ with a fine
    perturbation sampling, and CAMB 2);
  * anti-aliasing: number counts against densely sampled references, including
    the low-redshift case of #268.
* **AddressSanitizer** on all cases, and `make test`.
* **`python/test_class.py`** at TEST_LEVEL 1 and 2 with `COMPARE_OUTPUT_REF`
  against a fresh master `classyref`. The φφ, vector and aliasing differences are
  expected there and must be the only ones.
* **Timing of `transfer_init` alone,** against master, at 1 and 11 threads.
