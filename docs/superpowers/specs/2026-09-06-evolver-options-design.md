# An options struct for the shared evolver interface

**Status — 2026-09-06. SHIPPED as designed, with one contract tightened.**
`EvolverOptions` is in `include/evolver_options.h`; five evolvers, six call sites
and three test files converted; `evolver_erk_configure`/`evolver_erk_config`
deleted. `test-evolver-options` covers the validation. Evidence tags: **[M]**
measured here, **[C]** stated by the code, **[P]** earlier work, **[?]** none.

**Acceptance met.** Output is **bit-identical** before and after, across all five
evolvers on both the background and the perturbations (`cl.dat` and `pk.dat`,
`cmp` byte-equal), and on a lensed TT+PP+mPk reference run **[M]**. Both binaries
were built in the SAME build directory, because a different one costs ~1e-3 in
P(k) on identical source. Runtime on the reference run: 184.4 ms before,
178.2 ms after **[M]** — no regression; the 3% is run-to-run noise.

**Four contracts tightened, three of them from review.** `output` is required —
every evolver calls it without a null check. `x_sampling` needs at least two
points AND at least one at or beyond `x_ini`: each evolver opens with
`while (x_sampling[i] < x_ini) ++i` and no bound check, so a grid entirely below
the start runs off the end of the array — and requiring two points does not
prevent that, since both can be below, which is why reachability is checked
directly rather than just size. An evolver honouring `Timescale` *requires* it
rather than merely accepting it: the legacy `rk` takes its whole step size from
there and calls the function without a null check.

**And one section-3 claim was false.** `abstol` was documented as "honoured by
every evolver" while `evolver_erk_impl.h` and `evolver_etd.cpp` each hardwired
`1e-15` of their own — the exact silent ignore this struct exists to prevent,
inside the struct's own documentation. Both now read it from the options, and a
test asserts that raising it costs fewer steps in both evolver families.

**Contract tightened during implementation: `used_in_output` is REQUIRED and
non-null**, like `x_sampling`. Section 3 originally let it be null meaning "all
components". That is not implementable: `evolver_erk_impl.h` indexes it without
checking, so a null segfaults — found by the new test on its first run **[M]**.
Callers pass an all-ones array to mean all.

**An unplanned benefit.** Because validation rejects an option the chosen evolver
does not consume, each call site now builds its options *inside* the same
if/else that picks the evolver. That puts in one visible place what used to be
buried in five implementations: the Jacobian diagonal matters to `etd` alone, the
timescale pair to the legacy `rk` alone, the controller to `rkdp45`/`tsit5`. The
call sites previously passed all three to whichever evolver was selected.

Deleted along the way: `BbnTimescale` in `bbn_solver.cpp`, a no-op function that
existed only to fill a slot `ndf15` ignores.

## 1. The problem

Every evolver takes the same sixteen positional parameters. Several are ignored
by most of them, and which is which is documented only in prose, per header:

| parameter | ndf15 | rkck (`rk`) | rkdp45 / tsit5 | etd |
|---|---|---|---|---|
| `minimum_variation` | uses as `hmin` **[C]** | uses **[C]** | ignores | `(void)` cast **[C]** |
| `evaluate_timescale` + `timestep_over_timescale` | ignores | **uses** **[C]** | ignores | `(void)` cast **[C]** |
| `derivs_diagonal` | ignores | ignores | ignores | **uses** **[C]** |

Note the first row: `evolver_ndf15.h` says `minimum_variation` is "unused here",
while `evolver_ndf15.cpp:295` assigns it to `hmin` **[C]**. The prose is already
out of step with the code, which is the failure mode a table in the type system
would prevent.

Three consequences, all observed rather than hypothesised:

* **Setting an option that does nothing is silent.** Passing
  `timestep_over_timescale` to `ndf15` is accepted and ignored. This is the same
  shape as #395, where a knob that silently did nothing carved a hole in a
  posterior.
* **The interface has already failed to hold something.** `ErkControllerConfig`
  did not fit, so it escaped into process-wide mutable state. Three modules each
  call `evolver_erk_configure(ppr->erk_controller_config())` before evolving
  (`background_module.cpp:692`, `thermodynamics_module.cpp:3260`,
  `perturbations_module.cpp:658`), each with the same three-line comment warning
  that `Cosmology` is lazy and a second object parsed in between would otherwise
  decide this run's settings **[C]**. A fourth evolving module that forgets the
  call silently inherits another cosmology's controller.
* **Sixteen positional arguments invite errors that compile.** In this session:
  a null `x_sampling` passed to `ndf15`, which dereferences it unconditionally
  (`evolver_ndf15.cpp:182`) and segfaults **[M]**; and a no-op timescale function
  written purely to fill a slot `ndf15` ignores **[M]**.

Things the interface cannot currently express, which are wanted:

* **`abstol`.** Hardwired to `1e-15` in `evolver_ndf15.cpp:92` **[C]**. It is the
  error-weight floor, so it decides which variables are resolved relatively.
* **Step statistics.** `ndf15` collects six counters in a local `stepstat[6]`
  array and discards them **[C]**. ERK has a separate, process-wide
  enable/reset/get API — more global state.
* **Reduced maximum order.** `maxk = 5` is a local in `evolver_ndf15.cpp:93` **[C]**.

## 2. What is dropped, not moved

**`minimum_variation` leaves the interface entirely.** It is
`ppr->smallest_allowed_variation`, whose default is `DBL_EPSILON` and which is
**never read from the input file** — `input_module.cpp` validates it is
non-negative and nothing assigns it **[M]**. It is therefore a compile-time
constant threaded through every call site and every evolver. The two evolvers
that use it take `std::numeric_limits<double>::epsilon()` directly instead.

This is a behaviour change only for a caller that passed something other than
`DBL_EPSILON`, and no such caller exists **[M]**.

## 3. The struct

```cpp
struct EvolverOptions {
  // Honoured by every evolver.
  double rtol   = 1e-6;
  double abstol = 1e-15;   // error-weight floor; ndf15 previously hardwired it

  // Dense output, honoured by every evolver.
  const double* x_sampling = nullptr;   // REQUIRED and non-null; see below
  int x_sampling_size      = 0;
  const int* used_in_output = nullptr;
  void (*output)(...)          = nullptr;
  void (*print_variables)(...) = nullptr;

  // Evolver-specific. Each is "unset" at its default, and setting one an
  // evolver does not consume is an ERROR rather than an ignore.
  void (*derivs_diagonal)(...)      = nullptr;  // etd
  void (*evaluate_timescale)(...)   = nullptr;  // rkck
  double timestep_over_timescale    = 0.0;      // rkck
  int max_order                     = 0;        // ndf15; 0 = evolver's own
  std::optional<ErkControllerConfig> erk;       // rkdp45, tsit5
  EvolverStats* stats               = nullptr;  // where supported
};
```

Every evolver-specific field has an unset state that is distinguishable from any
value a caller would deliberately pass. That is what makes validation possible
without an `operator==` on each field, and it is why `max_order` is `0`-means-
default rather than `5`.

`x_sampling` stays **required and non-null**: `ndf15` dereferences it before its
first step, so a null is a segfault rather than "no dense output" **[M]**. The
struct makes that a checked precondition instead of a latent crash.

## 4. Validation

Each evolver declares what it consumes, and a shared helper rejects the rest:

```cpp
EvolverOptionsCheck(options, "ndf15", {EvolverFeature::MaxOrder,
                                       EvolverFeature::Stats});
```

`class_test_severe`, not `class_test`: these are run-configuration choices a
sampler never varies, and section 8 of the error conventions puts those in the
severe channel. A silently-ignored option is exactly the failure #395 punished.

## 5. Scope

**In:** the struct; the validation helper; `abstol`, `max_order` and `stats` made
reachable; `ErkControllerConfig` moved out of process-wide state into the struct;
`minimum_variation` dropped; five evolvers and six call sites updated.

**Out, deliberately:**

* ~~The ERK **stats** globals stay.~~ **REVERSED — they went too.** The
  reasoning that kept them was wrong; see section 7.
* No `std::variant` of per-evolver option structs. Making invalid combinations
  unrepresentable is stronger than rejecting them, but it forces every call site
  to name its evolver at compile time, which `generic_evolver`'s runtime function
  pointer exists precisely to avoid.
* No analytic-Jacobian hook. There is no consumer.

## 6. Acceptance

The refactor is **numerically inert**: every default reproduces today's
behaviour, so a run before and after must agree bit-for-bit. That is the
acceptance criterion, and it is checkable — `abstol` defaults to the hardwired
`1e-15`, `max_order = 0` means the existing `maxk = 5`, `erk` unset means the
existing defaults, and dropping `minimum_variation` restores `DBL_EPSILON`.

Cost must not regress measurably in `perturbations_module`, which is the hot
path. The struct is passed by `const&`, but that is an assumption until measured.

---

## 7. The stats globals went too (2026-09-06, second pass)

Section 5 deferred `evolver_erk_stats_enable/reset/get` and
`evolver_erk_histograms_enable/get` on the grounds that moving them meant
deciding how counters aggregate across the perturbation module's threads. That
turned out to be a smaller question than it looked, and leaving them created
exactly the asymmetry the struct exists to remove: `ndf15` collected six counters
into a local `int stepstat[6]` and discarded them, while the explicit evolvers
accumulated into process-wide atomics behind an enable/reset/get API.

**What the threading question actually was.** The perturbation workspace is
constructed per WAVENUMBER TASK, not per thread
(`perturbations_module.cpp:684`), so there is no per-thread slot to accumulate
into. But that makes the answer easier, not harder: each task accumulates its own
per-call counters into the module under one lock **per interval**, replacing a
relaxed atomic **per step**. That is thousands of times fewer synchronisations.

**It also removes a measurement hazard.** `class_profiled` carried
`CLASS_ERK_NOSTATS` with a comment explaining that the per-step atomics tilted
wall-time comparisons towards whichever integrator took fewer steps. Per-call
counters are plain increments needing no synchronisation, so there is nothing
left to switch off and the environment variable is gone. `CLASS_ERK_HIST` is gone
too: histograms are collected on the `evolver_histograms` input and printed when
there is something to print, rather than gated by a second hidden switch.

**Histograms kept, separately.** Counters are free; histograms cost a `log10` per
step to find a bin. So they are a distinct field (`EvolverOptions::histograms`)
with its own feature flag, and a null pointer is the gate. They remain
explicit-RK only — asking `ndf15` for them is an error, not an empty result that
would read as "no steps were taken".

**`stepstat[6]` is gone.** `ndf15` accumulates into a named `EvolverStats`, so
`stepstat[3]` is `counters.jacobians` and the index-to-meaning table that lived
in a comment nearly two hundred lines above its uses is unnecessary. Both evolver
families now report through one struct; each fills the fields its method actually
produces (`jacobians` and `lu_decompositions` for the implicit one,
`dense_points` and `exact_points` for those with a continuous extension) and
leaves the rest at zero. That residual asymmetry is a property of the methods,
not of the interface.

`tools/evolver_erk.cpp` went from ~110 lines of globals and atomics to 43 lines
of binning arithmetic with no state at all **[M]**.

`evolver_erk_config_test` is kept although the hazard it guards is now
structurally impossible: it still checks the surviving requirement -- that each
cosmology's controller reaches its own evolver calls -- and would catch a
reintroduction.

---

## 8. A regression the acceptance criterion should have caught (2026-09-06, review)

Section 6 makes the refactor's inertness the acceptance criterion, and section 7
reported it met. It was — **for the commit that introduced the struct.** The
later commit that moved the stats did not re-run the five-evolver sweep, and its
message still claimed inertness on the strength of an ndf15-only reference run.

What that hid: the perturbation module set `EvolverOptions::stats`
unconditionally, while the legacy `rk` and `etd` did not declare
`EvolverFeature::Stats`. Every run selecting either **aborted before taking a
step** **[M]**. Outputs cannot differ if the run does not start, so a bit-identity
check on a third evolver says nothing about the other two.

`tools/evolver_selection_test.cpp` is the standing version of that sweep: it runs
a small cosmology through every evolver the input accepts. It was confirmed to
FAIL on the reintroduced bug and pass without it **[M]**, because a regression
test that does not fail on its own bug is worthless.

Five more from the same review, all real:

* **`abstol` is feature-gated, not universal.** Section 3 called it "honoured by
  every evolver"; `evolver_erk_impl.h` and `evolver_etd.cpp` each hardwired their
  own 1e-15, and the legacy `rk` builds its error scale inside
  `generic_integrator` from `_TINY_` = 1e-30 — a different quantity that cannot
  be swapped without changing results. The first two now read the option; `rk`
  declines it through `EvolverFeature::AbsTol`, so setting it there is an error
  rather than a silent no-op.
* **`max_order` was unchecked.** It indexes five-element coefficient arrays, so 6
  was an out-of-bounds read and a negative value silently meant "default".
* **`timestep_over_timescale` was only checked non-zero.** A negative value sends
  `rk` backwards while its loop still asks for `x1 < x_end` — away from the
  endpoint — and NaN passed too. Spelled `!(x > 0)` now, which rejects both.
* **`etd` reports statistics** rather than declining them, so only the legacy
  `rk` abstains; its stepping lives in `generic_integrator` and has no counters.
* **The profiler's label lied.** It read "Explicit-RK perturbation steps" while
  now receiving ndf15's too.

The through-line is the one this document is about: every one of these is an
option that would have been accepted and quietly ignored, or accepted and used
out of range. The struct made them visible; it did not make them impossible.

## 8.1 The NaN check could not work where it was written

The `timestep_over_timescale > 0` guard added in section 8 was spelled
`!(x > 0.)` so that NaN would fail it too. It does under Apple clang. It does
**not** under the GCC on Linux CI, where `test-evolver-options` aborted on
exactly that case while the local suite stayed green **[M]**.

`EvolverOptionsCheck` is in `classpp`, which is built with `-ffast-math`, i.e.
`-ffinite-math-only`: the compiler may ASSUME no NaN occurs and fold any test for
one. `include/nonfinite.h` already records this, from measurement, and is blunt
about it -- *"No check written in a fast-math TU can work, however it is
spelled"*, the bit-pattern-through-memcpy form included -- which is why
`tools/nonfinite.cpp` is the single file CMake compiles without the flag.

So the root cause is not the spelling; it is that the check was written in the
wrong translation unit. `rtol`, `abstol`, `timestep_over_timescale` and `x_ini`
now all go through `IsNonFinite()`. Verified locally rather than by reasoning:
`nm -u` on the built object shows `IsNonFinite` as an **undefined** symbol, i.e.
a genuine cross-TU call the optimiser could not fold — there is no LTO in this
build, which is the property the mechanism rests on **[M]**.

The wider lesson matches section 8's: a green local suite said nothing about the
platform that mattered, and only the finite half of a check can be trusted to
mean what it says inside this library.