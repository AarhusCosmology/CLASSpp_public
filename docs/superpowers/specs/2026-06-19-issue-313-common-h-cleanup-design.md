# Issue #313 — `common.h` split, typed thrower, macro removal, `enum class`

**Date:** 2026-06-19
**Issue:** #313 (`v4 prep: split common.h; replace sprintf-buffer error macros with a typed thrower; enum class`)
**Slicing decision:** single sweep (one PR), like #312/#324.

## Problem

`include/common.h` (1064 lines) is a grab-bag mixing unrelated concerns: physics
constants, function-like macros (`MIN`/`MAX`/`SIGN`/`NRSIGN`, `typeof`), the error
machinery (`sprintf` into a `char[2048] ErrorMsg`, then either `return _FAILURE_`
or `throw`), output-formatting macros, the ~600-line `precision` struct, the
conversion enums, and the `Omega0gFromTcmb` helpers. Unscoped `enum`s leak very
generic names (`newtonian`, `synchronous`, …) into the global namespace. It is the
first header anyone opens, and for a "designed-from-scratch C++" pitch it makes the
opposite impression.

## Scope-finding results (verified against the real tree, 2026-06-19)

These findings shrink the work and justify deleting a whole subsystem rather than
porting it:

- **`hyrec/*.c` is fully self-contained** — it includes no CLASS headers and uses
  none of our macros (`class_*`, `MIN`/`MAX`/`SIGN`/`NRSIGN`).
- **No `.c` translation unit includes `common.h`.** `arrays`/`hyperspherical` are
  `.cpp` now; the only `.c` files in the tree are `hyrec/*.c`. Therefore the entire
  C (`return _FAILURE_`) macro branch — the `#ifdef __cplusplus … #else` dual
  definition — serves translation units that no longer exist. **Dead.**
- **`ErrorMsg`** has 8 occurrences, **all inside `common.h` itself** (the typedef +
  the macro-body buffers). The struct members it once described were removed in
  PRs #38/#303. Once the macros stop using a buffer, the typedef is unused.
- **`class_protect_sprintf`** is only *defined* (`tools/common.cpp`); its only
  callers are the macro bodies. **`class_protect_fprintf`** and
  **`class_protect_memcpy`** have **zero** users.
- **`class_call_parallel` / `class_test_parallel`** have **zero** call sites (only
  the `#define`s).
- **`typeof`**: the only non-definition occurrence is a *comment* in
  `input_module.cpp`. No code uses the macro.
- **Macro→std sweep size:** `MAX(` ~116 uses, `MIN(` ~93, `SIGN(` ~4 real uses
  (nonlinear_module, hyperspherical), `NRSIGN(` 3 uses (all hyperspherical),
  across ~13 real source files plus the `hermite*_interpolation_csource.h`
  textual includes. `class_test` (~429), `class_stop` (~66), `class_open` (~16)
  call sites are **untouched** — only the macro *bodies* change.
- **`class_call`**: 5 sites, all `class_call(get_xmin_generic(…), "", "")` in
  `transfer_module.cpp`; the message args are always empty and the 3rd is ignored.

## Design

### 1. Header split

Split `common.h` into focused headers; keep `common.h` as a thin umbrella that
`#include`s them so every existing `#include "common.h"` keeps compiling.

| Header | Contents |
|---|---|
| `include/constants.h` | Physics constants & conversion factors (`_PI_`, `_PIHALF_`, `_TWOPI_`, `_SQRT2_`, `_SQRT6_`, `_SQRT_PI_`, `_E_`, `_Mpc_over_m_`, `_Gyr_over_Mpc_`, `_c_`, `_G_`, `_eV_`, `_k_B_`, `_h_P_`) + the `Omega0gFromTcmb` / `TcmbFromOmega0g` inline helpers (pure constant conversions). |
| `include/errors.h` | **C++-only, throw-only.** Status/flags (`_SUCCESS_`, `_FAILURE_`, `_TRUE_`, `_FALSE_`), the `ThrowFormatted` declaration, and the throwing `class_call` / `class_test` / `class_stop` / `class_open` macros. No `ErrorMsg`, no `#ifdef __cplusplus` C branch, no `*_parallel`, no `class_*_message` / `class_build_error_string`. |
| `include/precision.h` | The precision-adjacent enums (`evolver_type`, `pk_def`, `file_format`, `tca_method`, `rsa_method`, `rsa_idr_method`, `ufa_method`, `ncdmfa_method`) + `struct precision`. |
| `include/common.h` | Thin umbrella: system includes, forward decls + `…ModulePtr` typedefs, `_VERSION_`, generic non-physics defines (`_HUGE_`, `_EPSILON_`, `_MAX_IT_`, `_QUADRATURE_MAX_`, `_QUADRATURE_MAX_BG_`, `_TOLVAR_`, `_OUTPUTPRECISION_`, `_COLUMNWIDTH_`, `_DELIMITER_`, `__CLASSDIR__`), `index_symmetric_matrix`, `class_define_index`, the output-format macros (`class_fprintf_*`, `class_store_*`), `get_number_of_titles`, and `#include "constants.h" / "errors.h" / "precision.h"`. The `#define typeof` and `MIN`/`MAX`/`SIGN`/`NRSIGN` are **removed**. |

`errors.h` is self-contained (defines `_SUCCESS_`/`_FAILURE_`/`_TRUE_`/`_FALSE_`
itself) so it can be included directly. `precision.h` includes what it needs
(`constants.h`/`errors.h`, `<float.h>` for `DBL_EPSILON`).

### 2. Typed thrower

```cpp
// errors.h (declaration); defined in tools/common.cpp
[[noreturn]] void ThrowFormatted(const char* func, int line, const char* fmt, ...);
```

Implementation: `vsnprintf` into a dynamically sized `std::string` (no 2048 cap),
prefix `"<func>(L:<line>) :"`, then `throw std::runtime_error(...)`. The error
message format is preserved.

The macros stay (they still do `#condition` stringizing + `__func__`/`__LINE__`)
but each body collapses to a single `ThrowFormatted` call, using compile-time
string-literal concatenation:

- `class_test(cond, args, …)` → `ThrowFormatted(__func__, __LINE__, "condition (" #cond ") is true; " args, ##__VA_ARGS__)`
- `class_stop(args, …)` → `ThrowFormatted(__func__, __LINE__, "error; " args, ##__VA_ARGS__)`
- `class_open(ptr, filename, mode)` → on `fopen` failure, `ThrowFormatted(__func__, __LINE__, "could not open %s with name %s and mode %s", #ptr, filename, #mode)`
- `class_call(function)` → if `(function) == _FAILURE_`, `ThrowFormatted(__func__, __LINE__, "error in " #function)`. The 5 call sites drop their always-empty `"", ""` args.

### 3. Dead-code deletions (the cleanup payoff)

- `class_call_parallel`, `class_test_parallel` (0 uses).
- The entire C `return _FAILURE_` branch of `class_call`/`class_test`/`class_stop`/`class_open`, plus `class_call_message`, `class_test_message`, `class_build_error_string`.
- `ErrorMsg` typedef and `_ERRORMSGSIZE_`.
- `class_protect_sprintf`, `class_protect_fprintf`, `class_protect_memcpy` — declarations in `common.h` and definitions in `tools/common.cpp`.
- `#define typeof` (the lone comment reference in `input_module.cpp` is reworded).

### 4. Macro → standard-library sweep

- `MAX(a,b)` → `std::max(a, b)`, `MIN(a,b)` → `std::min(a, b)` (`<algorithm>`).
  **Type gotcha:** the macros are type-agnostic; `std::max`/`std::min` require
  matching argument types. Mixed-type sites (e.g. `MAX(0, <double>)`) get the
  literal typed (`std::max(0.0, x)`), audited per site. This is a compile concern,
  not a numerical one.
- `SIGN(a)` → `std::copysign(1.0, a)`.
- `NRSIGN(a, b)` → `std::copysign(a, b)`.
  - These are deliberately **not** contorted to preserve the old sign-of-zero
    semantics (`SIGN(0.) == -1.`, `NRSIGN(a, -0.) == +fabs(a)`); `copysign` is the
    clean, standard intent. The inputs at which they differ are measure-zero
    (exact `0.`/`-0.`), and the only consumers are sign-change predicates in
    hyperspherical bisection and a correlation-sign product in nonlinear.
- `index_symmetric_matrix` is **kept** (not named by the issue; genuinely useful).

### 5. `enum` → `enum class`

Convert: `possible_gauges` (`source/perturbations.h`) and the `precision.h` set
(`evolver_type`, `pk_def`, `file_format`, `tca_method`, `rsa_method`,
`rsa_idr_method`, `ufa_method`, `ncdmfa_method`). Ripple:

- All use sites become scoped (`pk_def::delta_m_squared`, `tca_method::compromise_CLASS`, …).
- `precision` fields that stored these as `int` (`tight_coupling_approximation`,
  `radiation_streaming_approximation`, `idr_streaming_approximation`,
  `ur_fluid_approximation`, `ncdm_fluid_approximation`, …) are retyped to the enum
  with `enum::value` defaults.
- `(int) first_order_CLASS`-style comparisons in `input_module.cpp` and the modules
  are rewritten to scoped names.
- The `.ini` input path keeps its **integer** public interface: read an `int` into a
  temporary and `static_cast` to the enum (e.g. in `ReadPrecision`). User-facing
  keys like `tight_coupling_approximation = 4` are unchanged.

**Out of scope:** `quadrature_method` (`include/quadrature.h`) — freshly designed,
not named by the issue. Leave it.

## Verification

Per the project's standing rule, verify within **~0.1% tolerance** (handle Cl^TE
zero-crossings, never a blind max-rel-diff), **not** bit-identical. Almost all of
this refactor is numerically inert (header moves, thrower, `enum class`); only the
`copysign` swaps can perturb output, and only at exact `0.`/`-0.` inputs.

Plan:
1. Build (CMake) clean.
2. Run the full `TEST_LEVEL=2` suite; expect agreement within tolerance.
3. A handful of representative scenario diffs (as in #312/#324).
4. If `copysign` shifts anything, accept it and regenerate goldens + `classyref`.

## Non-goals / possible follow-ups

- Converting `get_xmin_generic` (and similar `_FAILURE_`-returning leaf helpers) to
  throw directly, retiring `class_call` entirely. Out of scope here.
- Splitting the output-formatting macros into their own header. Left in `common.h`.
- `quadrature_method` → `enum class`.
