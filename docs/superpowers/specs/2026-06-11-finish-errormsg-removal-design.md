# Finish the `ErrorMsg` removal — everything throws

**Date:** 2026-06-11
**Follow-up to:** [#38 / PR #303](https://github.com/AarhusCosmology/CLASSpp/pull/303) — [`2026-06-08-remove-error-message-member-design.md`](2026-06-08-remove-error-message-member-design.md)
**Branch:** `finish-errormsg-removal`

## Problem

PR #303 removed the mutable `error_message_` *member* and made the C++ error
macros throw from a stack-local buffer. But it deliberately stopped at the
"less mechanical" boundary: the `ErrorMsg error_message` **parameter** is still
threaded through the entire `tools/` numerical layer and the function-pointer
contracts that reach it.

Because every `tools/*.cpp` file compiles as C++, the `class_test` / `class_stop`
macros *inside* those functions already throw — they never write the passed-in
buffer. So the remaining `ErrorMsg error_message` / `ErrorMsg errmsg` parameters
are **vestigial**: dead weight kept only to satisfy the C-style
`int (*derivs)(..., ErrorMsg)` function-pointer signatures, which in turn force
every callback in `background` / `perturbations` / `thermodynamics` /
`primordial` / `input` to carry the buffer too.

The result is ~155 `class_call_failure(leaf(...), buf)` call sites plus their
local `ErrorMsg buf;` declarations, ~40 `array_*` signatures, the evolver /
integrator typedefs, and `generic_integrator_workspace::error_message` — all
plumbing for a buffer that is never read.

## Goal

Delete the `ErrorMsg` **parameter** and its scratch buffers everywhere they are
still threaded, so that error reporting is exclusively via thrown
`std::runtime_error`. Collapse the now-redundant call-site wrappers and macros.

## Non-goals / explicitly kept

- **The `ErrorMsg` *type* stays.** The throwing macros in `common.h` still use a
  stack-local `ErrorMsg _class_err_` to format the message before
  `throw std::runtime_error(_class_err_)`. We remove `ErrorMsg` as a *parameter*,
  not as a type.
- **Genuine numerical control-flow `int` returns stay as `int`** (decision:
  *keep as int return codes*). These carry **no** error message and are not
  errors — they are numerical status used for fallback/iteration:
  - `get_CF1` non-convergence (hot loop; falls back to Gegenbauer)
  - singular-matrix `ludcmp` / `sp_ludcmp`
  - `quadrature` `qm_auto` dispatch sentinel
  - internal evolver status
  These functions already take no `ErrorMsg`, so they are untouched except where
  they happen to sit in a header alongside renamed signatures.
- **No `int` → `void` conversion.** Pure-error functions keep returning `int`
  (always `_SUCCESS_`); we drop only the `ErrorMsg` param. Converting ~40
  `array_*` functions and 150+ call sites to `void` is a larger, riskier diff
  with no `ErrorMsg` benefit — left as a possible future cleanup. The
  function-pointer typedefs also keep `int` return, so the evolver's
  `class_call((*derivs)(...))` continues to act as a safety net.
- **HyRec (`hyrec/*.c`) is out of scope.** It does not include `common.h` and
  uses none of the `class_*` macros — it does its own
  `fprintf(stderr, ...); exit(1)`. Untouched.
- **`class_alloc` / `class_calloc` / `class_realloc` — nothing to do.** Already
  removed in earlier refactors (growTable → `std::vector`, RAII workspaces); not
  defined in `common.h`, zero usages.

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Control-flow `_FAILURE_` returns (get_CF1, ludcmp, qm_auto) | Keep as plain `int` return codes |
| PR boundary | Single PR, all layers |
| Pure-error `int` returns → `void` | No — keep `int` this PR |
| `class_call_except` (and sibling `class_test_except`) | **Remove the macros entirely**; convert sites to bare calls / explicit `try`/`catch` |

## Scope, layer by layer

### 1. Tools leaf signatures + bodies
Drop `ErrorMsg error_message` / `ErrorMsg errmsg` from the signatures and bodies of:
`tools/arrays.cpp` + `include/arrays.h` (largest surface, ~40 fns),
`tools/quadrature.cpp` + `include/quadrature.h`,
`tools/parser.cpp` + `include/parser.h`,
`tools/hyperspherical.cpp` + `include/hyperspherical.h`,
`tools/trigonometric_integrals.cpp` + `include/trigonometric_integrals.h`,
`tools/dei_rkck.cpp` + `include/dei_rkck.h`,
`tools/evolver_ndf15.cpp` + `include/evolver_ndf15.h`,
`tools/evolver_rkck.cpp` + `include/evolver_rkck.h`,
`tools/evolver_rkdp45.cpp` + `include/evolver_rkdp45.h`.
Error paths already throw via `class_test` / `class_stop`. Functions keep their
`int` return type. Internal control-flow `return _FAILURE_` stays.

### 2. Function-pointer typedefs
Remove the trailing `ErrorMsg error_message` from the `derivs` / `output` /
`timescale_and_approximation` / `print_variables` / `fzero_Newton` callback
typedefs in the evolver and integrator headers, and from `numjac` / `new_linearisation`
/ `initialize_*` helper signatures.

### 3. `generic_integrator_workspace::error_message`
Delete the member from `include/dei_rkck.h` and any writes to it.

### 4. Callbacks in modules
Update the matching callback definitions and declarations to drop the param:
`background_module` (.cpp/.h), `perturbations_module` (.cpp/.h),
`thermodynamics_module` (.cpp/.h), `primordial_module` (.cpp/.h),
`input_module` (.cpp/.h, the `fzero_Newton` shooting callback),
`species/ncdm_base_species.cpp`, `species/dncdm_species.cpp`, `species/thermo_context.h`.

### 5. Call sites
- ~155 `class_call_failure(leaf(...), buf)` → bare `leaf(...);`, delete the local
  `ErrorMsg buf;` / `ErrorMsg errmsg;`.
- `class_call(leaf(..., errmsg), errmsg, errmsg)` (parser reads in `input_module`,
  etc.) → bare `leaf(...);`.
- Inside the evolvers, `class_call((*derivs)(..., error_message), ...)` → bare
  `(*derivs)(...);`.
- Then audit the remaining `class_call(` sites (170 today): bare-ify any that now
  wrap a throwing function with no buffer to pass. (Several already-bare C++ member
  wrappers exist from #303; this finishes the few stragglers, e.g.
  `background_module.cpp:1590` `background_functions(...)` inside `background_derivs`.)

### 6. Remove the `*_except` cleanup macros
Decision: remove `class_call_except` **and** its sibling `class_test_except`
(identical "run cleanup, then rethrow" pattern; one site shares the same `gi`
cleanup, so they must go together). Convert the sites:

`class_call_except` (6):
- `thermodynamics_module.cpp:450`, `:455` — empty cleanup → bare call.
- `primordial_module.cpp:274`, `:351`, `:371` — `primordial_free()` cleanup →
  `try { call(); } catch (...) { primordial_free(); throw; }`.
- `primordial_module.cpp:1799` — `cleanup_generic_integrator(&gi)` cleanup →
  `try`/`catch` with that cleanup.

`class_test_except` (4):
- `nonlinear_module.cpp:2328`, `primordial_module.cpp:312`, `:1802`,
  `thermodynamics_module.cpp:2343` — rewrite as `if (cond) { cleanup; class_stop(...); }`
  (which throws) or wrap in `try`/`catch`, matching the local cleanup.

Then delete both `class_call_except` and `class_test_except` definitions
(C and C++ blocks) from `include/common.h`.

### 7. Macro cleanup in `common.h`
- Delete `class_call_failure` (no users left).
- If `class_call` ends with zero users after the audit, remove it; otherwise leave
  its throwing C++ form intact.
- Keep `class_call_try` only if still used (verify); remove if dead.

### 8. Binding layer
- `classy.pyx` (hand-edited): update the extern declarations and call sites for
  the hyperspherical functions (~lines 132–139, `hyperspherical_HIS_create` /
  `_free` / interpolate) to drop `char* error_message`.
- `generate_wrapper.py`: the `ErrorMsg` ctypedef (line 123) and its
  `allowed_types` entry (line 260) may stay (harmless — `ErrorMsg` is still a real
  type); remove only if wrapper generation breaks. **`cclassy.pxd` is generated —
  never hand-edit** (see `feedback_cclassy_pxd_generated`).
- `get_my_py_error_message()` (the exception→Python bridge) is unrelated — keep.

## Risk & verification

This is pure error-path plumbing with **zero numerical change** — no math, no
control flow, no data is touched. Expectation: **bit-identical output** (same bar
PR #303 met).

Verification:
1. Build all relevant manifests (Makefile / `setup.py` / Xcode project).
2. Run the 96 `nose_tests`.
3. Diff a reference scenario's output for bit-identity.

If any drift appears (it should not), apply the standard rule: verify TT < 0.1%,
handle Cl^TE zero-crossings, never blind max-rel-diff
(`feedback_no_bit_identical_requirement`, `feedback_vectorization_reduction_drift`).

## Watch-outs

- **Don't convert control-flow `return _FAILURE_` to a throw.** `get_CF1`,
  `ludcmp`/`sp_ludcmp`, `qm_auto` are numerical status, not errors. They have no
  `ErrorMsg` param, so they should simply not be touched.
- **Function-pointer signature changes ripple.** A typedef edit must land together
  with the evolver-internal call, every callback definition/declaration, and every
  site that takes the function's address — or it won't compile. Do each callback
  family (derivs/output/timescale) as one atomic edit.
- **`numjac` / `fzero_Newton` pass the derivs pointer onward** — their internal
  `(*derivs)(...)` calls and their own `ErrorMsg` params must change in lockstep.
- **`parser_read_*` return value is real but already throwing on malformed input**;
  the `flag` out-param (found/not-found) is the control signal and is unaffected.
