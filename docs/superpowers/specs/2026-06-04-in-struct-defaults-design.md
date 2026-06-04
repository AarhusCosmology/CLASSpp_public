# Design: Move all struct defaults into the structs; delete `input_default_params` / `input_default_precision`

**Date:** 2026-06-04
**Status:** Approved (brainstorming)

## Motivation

Each `struct <foo>` (config/params known at input) is paired with a
`class <foo>_module` (quantities computed during construction). Default values
for the input structs should live **in the struct itself** rather than being
scattered through `InputModule::input_default_params()` and
`InputModule::input_default_precision()`.

Benefits:
1. Fewer lines of code (delete two functions).
2. All default settings — the "magic numbers" — are localized next to the
   fields they belong to.
3. Enables (later, not now) generating default-valued input structs on the fly
   from the Python wrapper, so all defaults can be inspected from Python.

## Current state

This migration is **already partly done**: `precision` fully uses in-struct
member initializers, and `background` partially does (`T_cmb`, `h`, `Omega0_k`,
`a_today`, `background_verbose`). What remains in the two default functions:

- `input_default_params()` — computed/derived scalar defaults and
  array-element defaults.
- `input_default_precision()` — one validation `class_test` and runtime
  `class_dir`-prefixed file-path concatenation.

## Key constraints (verified)

- **No constructors.** Every remaining default is expressed as an in-class
  member initializer. Where a default is derived from other fields, the
  initializer references earlier-declared members (C++ evaluates default member
  initializers in declaration order), so derived fields must be **declared after
  the fields they depend on** (reordering is permitted).
- **`generate_wrapper.py` already tolerates in-struct initializers.** Its parser
  keeps only `words[1] + words[2]` of a member line, so `double x = 5;`,
  `double H0 = h * 1.e5 / _c_;`, and `double a[N] = {1.};` all parse correctly
  (the `= …` tail is ignored). A free function defined in a header (e.g. in
  `common.h`) sits outside any `struct`, so the parser skips it entirely.
- The structs are owned by value as members of `InputModule`
  (`input_module.h:37-43`) and default-constructed in C++. The Cython layer
  never stack-allocates them (no `cdef background` in `classy.pyx`); it only
  reads fields through pointers.
- The derived "defaults" in `input_default_params` were only ever the
  *no-user-input* fallback: the readers in `input_read_parameters` recompute
  the same quantities when the user supplies `h`/`T_cmb`/`omega_b`/`n_t`/etc.
  Those readers are **unchanged**. In particular `Omega0_g` is recomputed
  *unconditionally* by the reader (`input_module.cpp:730-732`, the
  none-of-three branch), so its in-struct default never affects a real run — it
  exists only so a default-constructed struct holds the physical value (future
  goal 3).
- `FileName` path fields are **not** exported to the wrapper (`allowed_types`
  has `FileArg` but not `FileName`), so adding string-literal initializers to
  them is invisible to Cython.

## Changes

### Bucket 1 — Derived scalars → in-class member initializers (no constructor)

These derive purely from other in-struct defaults; reproduce the exact
expressions currently in `input_default_params()`. The depended-on field must be
declared first.

- `background` (reorder so `h`, `T_cmb` precede the derived fields):
  - `double H0 = h * 1.e5 / _c_;` (`_c_` is a `common.h` macro, in scope)
  - `double Omega0_b = 0.022032 / (h * h);` (avoids `pow`/`<cmath>`)
  - `double Omega0_g = Omega0gFromTcmb(T_cmb, h);` (see below)
- `primordial`: `r` (line 75) and `n_s` (line 72) already precede `n_t`/`alpha_t`
  (lines 76-77), so no reorder is needed:
  - `double n_t = -r/8. * (2. - r/8. - n_s);`
  - `double alpha_t = r/8. * (r/8. + n_s - 1.);`

**`Omega0_g` and the photon-density helper.** `Omega0gFromTcmb` is currently a
`static` method of `PhotonsSpecies`, and `species/photons.h` `#include`s
`background.h` — so calling it from a `background.h` initializer would be a
circular include. Resolution: **move both `Omega0gFromTcmb` and its inverse
`TcmbFromOmega0g` out of `PhotonsSpecies` into free functions in
`include/common.h`** (no new file; `common.h` already defines every constant
they combine — `_c_`, `_k_B_`, `_h_P_`, `_G_`, `_Mpc_over_m_`, `_PI_`).
  - Both become plain **`inline` free functions with their bodies copied
    verbatim** (keep the `pow()` calls). A runtime non-static member initializer
    does **not** require `constexpr`, so there is no reason to change the
    arithmetic.
  - **Do not** rewrite to `constexpr`. `constexpr` would force replacing `pow()`
    (not `constexpr` in C++17) with explicit multiplications, and `pow(x, 4.)` is
    not guaranteed bit-identical to `x*x*x*x`. Because the *reader*
    (`input_module.cpp:731,738`) calls the same `Omega0gFromTcmb` on real runs,
    that ULP shift would propagate into actual results and churn baselines, for
    no functional gain. Keeping `pow` makes the change bit-identical.
  - Update the call sites (`input_module.cpp:731,738,744,749`; the
    `input_default_params` line 2660 is deleted) plus
    `species/photons_formula_test.cpp` to drop the `PhotonsSpecies::` qualifier;
    remove the two declarations from `species/photons.h` and the bodies from
    `species/photons.cpp`.

Trade-off accepted: the photon-density physics moves from `PhotonsSpecies` to
`common.h`. This was chosen over (a) a `background` constructor and (b) inlining
the formula directly into the struct (formula duplication), in order to keep a
single source of truth with no constructor and no new file.

### Bucket 2 — Array & cross-struct scalars → plain in-struct initializers

No constructor needed (aggregate init is invisible to the wrapper parser):

- `transfers`: `double selection_bias[_SELECTION_NUM_MAX_] = {1.};`,
  `double selection_magnification_bias[_SELECTION_NUM_MAX_] = {0.};`
- `perturbs`: `double selection_mean[_SELECTION_NUM_MAX_] = {1.};`,
  `double selection_width[_SELECTION_NUM_MAX_] = {0.1};`
- `z_max_pk`: `double z_max_pk = 0.;` on **both** `perturbs` and `spectra`.
  (`z_pk[0]`'s default is already `0.`, so this matches
  `psp->z_max_pk = pop->z_pk[0]`. The reader still does the real
  `MAX(...)` derivation when P(k) output is requested.)

### Bucket 3 — Runtime string paths → relative default + `ResolveDataPaths()`

Chosen approach: relative-path defaults in the struct, plus one small member
method that prepends the runtime `class_dir` in place. Preserves exact current
ordering/semantics (default resolved, *then* `parse()` override replaces
verbatim) with **zero** changes to the `thermodynamics_module` use-sites.

- In `include/common.h`, give each path field a relative default:
  - `FileName sBBN_file = "/bbn/sBBN_2017.dat";`
  - `FileName hyrec_Alpha_inf_file = "/hyrec/Alpha_inf.dat";`
  - `FileName hyrec_R_inf_file = "/hyrec/R_inf.dat";`
  - `FileName hyrec_two_photon_tables_file = "/hyrec/two_photon_tables.dat";`
- Add `void precision::ResolveDataPaths();` (declared in `include/common.h`,
  defined in `input_module.cpp` alongside `precision::parse()`) — prepends
  `this->class_dir` to each of the four path fields in place (same
  `strncpy`/`strcat` logic now in `input_default_precision`).
- In `input_read_precisions`, replace the `input_default_precision()` call
  (currently line 608, after `class_dir` is determined at 604) with
  `precision_.ResolveDataPaths();`.

### Deletions / relocation

- Delete `InputModule::input_default_params()` and its declaration.
- Delete `InputModule::input_default_precision()` and its declaration; remove its
  call at line 176 (`input_default_params()`).
- Relocate the lone `class_test(smallest_allowed_variation < 0, ...)` validation
  to just after `ppr->parse(file_content_)` (line 618) — where it actually
  belongs (it validates a possibly-overridden value; today it ran before
  `parse()` and only ever saw the default).

## Out of scope

- Python-side default-struct generation (point 3) — enabled by this work but
  **not implemented now**.
- Any change to reader logic, recomputation, or numerical results.

## Verification

- Build the C++ `class` binary and run the test suite.
- Run `generate_wrapper.py`; confirm `cclassy.pxd` regenerates without the new
  initializers corrupting struct parsing, and the wrapper compiles.
- Build `species/photons_formula_test.cpp` (or the equivalent target) to confirm
  the relocated `Omega0gFromTcmb`/`TcmbFromOmega0g` round-trip still holds.
- Numerical equivalence: a default-config run should be **bit-identical**, since
  every initializer reproduces the exact same default expressions and the photon
  formula is moved verbatim (`pow` kept). Still spot-check a few scenarios within
  ~0.1% tolerance (handle Cl^TE zero-crossings; no blind max-rel-diff).
