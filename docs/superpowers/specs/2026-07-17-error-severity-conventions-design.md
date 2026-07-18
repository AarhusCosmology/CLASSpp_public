# Error-severity conventions: `class_test_severe` and removal of the InputModule rewrap

**Date:** 2026-07-17
**Status:** Approved
**Scope:** One combined PR.

## Problem

CLASSpp signals errors through two exception types that the Python wrapper maps
to sampler-visible behaviour (`classy.pyx` `raise_my_py_error`):

- `std::invalid_argument` → `CosmoSevereError` — the sampler aborts the chain.
- `std::runtime_error` → `CosmoComputationError` — the sampler rejects the
  point (−inf likelihood) and the chain continues.

Today the codebase expresses this distinction inconsistently:

1. `class_test`/`class_stop` (`include/errors.h`) always throw `runtime_error`,
   so severe errors cannot be raised through the macros. Consequently
   `species/` uses ~107 raw `throw std::invalid_argument(...)` with verbose
   string concatenation (`std::to_string` formats doubles poorly) and no
   function/line context, while `source/` uses ~323 `class_test`.
2. `InputModule`'s constructor blanket-rewraps every `runtime_error` raised
   during input into `invalid_argument` (added in PR #204, before hook-based
   shooting and the species plugins existed). The input phase now runs real,
   point-dependent computations — NCDM Omega↔m Newton iteration, grey-body
   bracketing, and full trial cosmologies inside `ShootingResidual` — whose
   failures the rewrap wrongly promotes to chain aborts.
3. A third mechanism, `ThrowRuntimeErrorIf`/`ThrowRuntimeError`
   (`tools/exceptions.cpp`, 7 call sites), duplicates `class_test` minus the
   context.
4. `SpeciesCollection` throws `logic_error`/`out_of_range`, which match neither
   branch of the classy mapping and surface as `NotImplementedError`.
5. `readDoubleList` (`input_module.cpp`) catches the parser's correctly-typed
   `invalid_argument` and demotes it to `runtime_error` via `class_stop`,
   relying on the constructor rewrap to promote it back.

## The convention (the semantic core)

**Severe (`invalid_argument`, macros `class_test_severe`/`class_stop_severe`):**
a check that depends only on the *structure* of the input — which keys are
present, mutual exclusivity ("give either T or gstar_dec, not both"), missing
required fields, unparseable strings, unknown enum names ("incomprehensible
input '%s'"), and list-length/count consistency (e.g. `ncdm_psd_filenames` has
N entries but `N_ncdm` = M). Structure is identical at every point of an MCMC
chain, so aborting is always correct.

**Computation (`runtime_error`, macros `class_test`/`class_stop`):** anything
that depends on a parsed *numeric value* — including plain range checks such as
"T must be positive", "a_c must lie in (0,1)", "gstar_dec must exceed 43/11" —
and every numerical failure (convergence, bracketing, integration), wherever it
occurs, **including during the input phase and inside shooting trial builds**.

The governing rule: **a `class_test_severe` must never depend on
possibly-varying parameters.** A sampler can propose any numeric value; a bad
value at one proposal is a bad point, not a bad file. Aborting a long chain on
one unlucky proposal destroys work; rejecting is always the safe direction.
This rule deliberately includes discretization-style numeric parameters
(`momenta_bins`, tolerances): the simple "any numeric value → computation"
line is mechanically checkable and errs on the side that cannot destroy a
chain.

**Accepted trade-off:** a user who *fixes* a bad value in their `.ini` (e.g.
`gstar_dec = 2` held constant) gets a chain that rejects every point with the
same message rather than an immediate abort. CLASS cannot distinguish
fixed-from-varied from the inside; the repeated rejection message in the
sampler log is the accepted, observable cost.

**Programmer-error invariants** (mis-implemented hooks, `SpeciesCollection`
misuse) keep `std::logic_error`. The classy fallback branch changes from
`NotImplementedError` to `CosmoSevereError` so unknown exception types abort
loudly.

**API-argument validation** (e.g. requesting `lmax` beyond what was computed in
`SpectraModule`/`LensingModule`) is severe: such arguments are fixed
per-analysis configuration, cannot be proposed by a sampler, and classy already
raises `CosmoSevereError` for the same conditions in its own Python-level
checks.

The convention is documented in a comment in `include/errors.h` and a short
section in `STYLE.md`.

## Design

### 1. Macros and plumbing

- `tools/common.cpp`: refactor `ThrowFormatted` so the existing formatting body
  (identical `"<func>(L:<line>) :<prefix><message>"` format) is shared by two
  `[[noreturn]]` entry points:
  - `ThrowFormatted(...)` — throws `std::runtime_error` (unchanged signature
    and behaviour; zero churn at existing call sites);
  - `ThrowFormattedSevere(...)` — throws `std::invalid_argument`.
- `include/errors.h` adds, with bodies identical to their existing
  counterparts (per decision: uniform message format, including the
  `condition (...) is true;` prefix):
  - `class_test_severe(condition, args, ...)`
  - `class_stop_severe(args, ...)`
- Headers stay plain C++ (no C-compat guards). `cclassy.pxd` is
  auto-generated; do not hand-edit.

### 2. InputModule

- Delete the constructor's `catch (runtime_error) → throw invalid_argument`
  rewrap (`input_module.cpp:216–218`).
- Delete `readDoubleList`'s try/catch entirely; the parser's
  `invalid_argument` is already correctly typed and now propagates unchanged.
- Audit all ~90 `class_test` and 12 `class_stop` in `input_module.cpp` under
  the structure/value rule: structural checks (exclusivity, incomprehensible
  strings — all `class_stop`s remaining after the `readDoubleList` deletion
  qualify) become `_severe`; numeric-value checks stay plain and are otherwise
  untouched.
- The shooting hook-contract guard (`input_module.cpp:2738`, mis-implemented
  species hook) becomes `std::logic_error`.
- `tools/parser.cpp` already throws `invalid_argument` for parse/structure
  errors — unchanged, now correct end-to-end.

### 3. Species port

- Classify each of the ~107 raw `invalid_argument` throws in `species/`:
  - **Structural** → `class_test_severe`/`class_stop_severe` one-liners with
    printf formatting (`%g` for doubles instead of `std::to_string`).
  - **Value-range** → plain `class_test`. This flips the exception type for
    checks currently misclassified as severe; the clearest cases are sampled
    parameters: `fluid.cpp` pheno-axion `a_c`, `w_fld_i`, `nu_fld`,
    `Theta_initial_fld`; axion `T`/`gstar_dec` range checks;
    `wdm_decay_product` numeric ranges; grey-body `alpha`/`x`/`q0` and
    `M2`/`M3` positivity.
- Grey-body bracketing/moment-ratio failures
  (`greybody_ncdm_species.cpp:89,109`) stay `runtime_error` (become
  `class_test`/`class_stop`): M2/M3 are sampleable.
- The six existing `class_test` in `species/` (numerics) are already correct —
  unchanged.
- `species_collection.cpp` keeps `logic_error`/`out_of_range` (programmer
  invariants), which now map to `CosmoSevereError` via the classy fallback.
- `ncdm_base_species.cpp:491` `ThrowRuntimeError` (Newton mass) →
  `class_stop` (computation).

### 4. Retire the third mechanism

- Delete `ThrowRuntimeErrorIf`/`ThrowRuntimeError` from
  `tools/exceptions.{h,cpp}`; keep `get_my_py_error_message` (classy needs it).
- Convert the 7 call sites: `lensing_module.cpp`/`spectra_module.cpp` `lmax`
  guards → `class_test_severe` (API-argument validation, see convention);
  `ncdm_base_species.cpp:491` → `class_stop` as above.

### 5. classy.pyx

- `raise_my_py_error`'s fallback branch: `NotImplementedError` →
  `CosmoSevereError`.

## Behavioural consequences

- A structurally bad parameter file still aborts a chain — now because
  validation itself throws `invalid_argument`, not via blanket promotion.
- An extreme but well-formed point whose input-phase computation fails (NCDM
  mass Newton, thermo failure inside a shooting trial build) now rejects the
  point instead of aborting the chain.
- Newton overshoot during shooting into a range rejected by a value check now
  throws `runtime_error`: `fzero_Newton`'s Jacobian retry loop absorbs it, and
  if it escapes the main step it becomes a rejected point. (A damped/step-
  halving retry in `fzero_Newton`'s main loop remains a separate follow-up
  numerics improvement, no longer semantically required.)
- Value-range errors surfaced through classy change type:
  `CosmoSevereError` → `CosmoComputationError`. This is the intended behaviour
  change of the PR.
- Message texts change formatting (macro prefix + printf) where raw throws are
  converted; no numerical output changes anywhere.

## Tests and verification

- New `errors_test.cpp`: each of the four macros throws the correct type with
  the expected message shape (function name, line, condition text, formatted
  args). Registered in **both** `CMakeLists.txt` and the Makefile
  `TEST_TARGETS` list (CI only builds targets named there).
- Existing species factory tests that provoke *structural* errors and catch
  `invalid_argument` keep passing and become convention guards. Tests that
  provoke *value-range* errors are updated to expect `runtime_error`. Audit
  tests for message-substring matching during the port.
- End-to-end classy check: a structurally bad input raises `CosmoSevereError`;
  a bad value and a forced input-phase computation failure raise
  `CosmoComputationError`.
- No golden/`classyref` regeneration: error paths only.

## Out of scope (deliberate)

- Damped-Newton step-halving in `fzero_Newton`'s main loop (file follow-up
  issue).
- The ~130 string-name species branches (separate backlog).
- Any change to hyrec (C code, own error handling).
