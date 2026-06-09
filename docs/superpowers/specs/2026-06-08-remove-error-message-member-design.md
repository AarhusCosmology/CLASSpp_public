# Remove `error_message_` member; throw error strings from locals

**Date:** 2026-06-08
**Issue:** [#38 — Thread safety](https://github.com/AarhusCosmology/CLASSpp/issues/38)
**Branch:** `fix-issue-38-remove-error-message-member`

## Problem

Issue #38 states that all public CLASS module methods are `const` except for a
single mutable buffer, `error_message_`, which cannot easily be mutex-protected,
so "if two threads call a public function to the same module, and both fail,
there is a possibility of a crash."

Two findings from investigating the current code:

1. **The hazard is real but narrow, and mostly benign.** The concurrency model
   in this codebase is the `tools/thread_pool.h` task system (not OpenMP). Five
   module methods submit tasks that capture `this` and use the shared
   `mutable ErrorMsg error_message_` member as the format-and-throw buffer:
   - `source/perturbations_module.cpp:778` and `:815`
   - `source/transfer_module.cpp:266`
   - `source/primordial_module.cpp:1414`
   - `source/spectra_module.cpp:861`
   - `source/lensing_module.cpp:303`

   If two concurrent tasks fail at the same time, they `vsnprintf` into the same
   2048-byte buffer while one reads it to build the `std::runtime_error`
   `std::string`. This is a formal data race (UB). The realistic symptom is a
   garbled/interleaved/truncated error message, **not** a segfault — the buffer
   is fixed-size and `vsnprintf` always null-terminates within bounds, so the
   `std::string` ctor's `strlen`+copy stays in range. A true overrun needs a
   read before any null was ever written (the ctor only sets
   `error_message_[0] = '\n'`), a narrow window. And it only triggers when two
   concurrent operations *both* hit an error path, which is exceptional — hence
   "never verified."

2. **The member is unnecessary in the C++ throw model.** `error_message_` is a
   leftover from the C `return _FAILURE_` model, where a function wrote its
   message into a persistent struct field so callers up the chain could read it.
   In the C++ build the error macros (`include/common.h`, `#ifdef __cplusplus`
   block) format the message and immediately `throw std::runtime_error(...)`,
   which copies the text into the exception's `std::string`. The buffer is never
   needed again — it is pure scratch. The codebase already demonstrates the
   correct pattern in `tools/exceptions.cpp` (`ThrowRuntimeError` /
   `ThrowRuntimeErrorIf`), which use a **stack-local** `char error_message[2048]`
   and throw — inherently per-call and race-free.

This work is adjacent to
[`2026-06-04-eliminate-fixed-size-string-buffers-design.md`](2026-06-04-eliminate-fixed-size-string-buffers-design.md)
but distinct: that effort targeted fixed-size string buffers generally; this one
removes the specific module error-scratch buffers that create the #38 race.

## Goal

Remove the `error_message_` module member (and the analogous CLASSpp C-struct
`error_message` scratch fields), so error messages are formatted into stack
locals at the throw site. This eliminates the shared mutable state behind #38
while keeping the `const`-method design intact.

## Non-goals

- **Do not** convert the `tools/` numerical layer off its C return-code
  contract. Functions in `arrays`, `hyperspherical`, `quadrature`, `sparse`,
  `evolver_ndf15`, and the `generic_integrator` (dei_rkck) keep their
  `ErrorMsg`/`return _FAILURE_` API. Some of their `_FAILURE_` returns are
  **control flow, not errors** (e.g. `get_CF1` non-convergence falls back to
  `CF1_from_Gegenbauer`; the ndf15 evolver converts failures into step-size
  rejection). Converting these to `throw` is risky and orthogonal to #38.
- **Do not** touch HyRec `.c` files or the C (non-`__cplusplus`) macro
  definitions. HyRec includes no CLASSpp module headers, so the return-code
  macros stay byte-for-byte unchanged.
- **Do not** remove `generic_integrator_workspace::error_message` (dei_rkck) —
  it is the integrator's internal C-style scratch, part of the `tools/` layer we
  are leaving alone, and carries the same control-flow-`_FAILURE_` risk.

## Design

### Macro layer (`include/common.h`, `#ifdef __cplusplus` block only)

The throwing macros stop taking a caller-owned buffer and instead format into a
macro-local `ErrorMsg` before throwing:

- `class_test(condition, args, ...)` — drop `error_message_output`; build a
  local message and `throw std::runtime_error(...)`.
- `class_stop(args, ...)` — drop the leading `error_message_output`; same.
- `class_open(pointer, filename, mode)` — drop `error_output`; throw from a
  local on open failure.
- `class_call_except(function, list_of_commands)` — drop both buffer args; keep
  the `try { function } catch (...) { list_of_commands; throw; }` shape so
  cleanup-on-failure still runs (6 sites in thermodynamics/primordial).
- **New** `class_call_failure(function, msg_buffer)` — for the `tools/` leaf
  calls that genuinely `return _FAILURE_`:
  `{ if ((function) == _FAILURE_) { build chain naming #function and msg_buffer; throw std::runtime_error(...); } }`.
- `class_call(...)` on throwing members → replaced by **bare calls** at the call
  sites (members never `return _FAILURE_` in the all-C++ build; they throw, so
  the wrapper is a no-op). Remove the C++ `class_call` macro once unreferenced.

The C definitions (the non-`__cplusplus` branch) are left unchanged.

### Call-site transformation (~1483 module refs + 6 struct-field refs)

- **~1455 throwing-member sites:**
  `class_call(f(...), error_message_, error_message_);` → `f(...);`
  Done by script, with manual review of multi-line invocations. Cross-module
  forms collapse the same way, e.g.
  `class_call(background_module_->background_at_tau(...), background_module_->error_message_, error_message_);`
  → `background_module_->background_at_tau(...);` — the buffer references vanish.

- **~28 leaf-boundary sites** (where `error_message_` is also passed *into* a
  `tools/` leaf function): declare a local `ErrorMsg buf;` in scope and use the
  new macro, e.g.
  ```cpp
  ErrorMsg buf;
  class_call_failure(array_spline_table_lines(..., buf), buf);
  ```

- **`class_test` / `class_stop` sites:** drop the `error_message_` argument, e.g.
  `class_test(q_size_ < 2, error_message_, "buggy q-list");` →
  `class_test(q_size_ < 2, "buggy q-list");`

- **6 CLASSpp C-struct scratch sites** (`psp->error_message` ×4,
  `ppr->error_message` ×1, `pba->error_message` ×1): replace with a local
  `ErrorMsg`.

### Threaded lambdas

The five thread-pool lambdas (listed under Problem) currently capture and write
`this->error_message_`. After the transformation they make bare throwing calls;
the captured shared buffer is gone, closing the race at its source. An exception
thrown inside a task is stored in its `std::future` and re-thrown at
`future.get()` as today.

### Removals

- `mutable ErrorMsg error_message_;` and the `error_message_[0] = '\n';` ctor
  line in `source/base_module.h`.
- `ErrorMsg error_message_;` in `source/input_module.h`.
- `ErrorMsg error_message;` field from the `spectra` struct (`source/spectra.h:50`),
  the `precision` struct (`include/common.h:1080`), and the `background` struct.

### Left untouched (consistent C contract)

- All `tools/` leaf functions and their `ErrorMsg errmsg` out-params.
- `generic_integrator_workspace::error_message` (`include/dei_rkck.h:44`, 38 uses
  in/under thermodynamics recfast integration).
- HyRec `.c` files and the C (return-code) macro definitions.

## Boundaries between units

- **Macros** (`include/common.h`): the only place that formats and throws. After
  the change, the throw buffer is a macro-local; no caller state is involved.
- **Module code** (`source/`, `species/`): makes bare throwing calls; declares a
  local `ErrorMsg` only at the ~28 `tools/`-leaf boundaries and 6 struct-field
  sites.
- **`tools/` numerical layer**: unchanged C return-code API; the boundary
  contract is "leaf returns `_FAILURE_` and writes a caller-provided `ErrorMsg`;
  the module-side caller wraps it in `class_call_failure` with a local buffer."

## Error handling

Behavior is unchanged from the caller's perspective: failures throw
`std::runtime_error` with the same message content. The only difference is that
the message text lives on the stack of the throwing call instead of in a shared
member, so concurrent failures no longer race.

## Testing / verification

This is a pure error-path refactor: it does not touch any numerical code, so
**bit-identical** output is the expectation, to be verified rather than assumed.

1. Build all three manifests, which must stay in sync (Makefile, `setup.py`,
   `CLASS.xcodeproj`).
2. Rebuild `classy`. No public-header signatures change, so
   `generate_wrapper.py` output is unaffected; confirm the build regardless.
3. Run `nose_tests` (including `python/test_greybody.py`, now wired into CI).
4. Compare CMB and P(k) output against a pre-change reference at ~0.1% tolerance
   with TE zero-crossing handling; expect bit-identical here.
5. Sanity-check the error path itself: trigger a known failure (e.g. an invalid
   input) and confirm the thrown message is well-formed and complete.

## Risks

- **Multi-line `class_call` invocations.** Many sites span several lines; the
  scripted transform must handle these and be reviewed. Mitigation: stage the
  edit, compile, and diff; the compiler flags any malformed residue.
- **Hidden `class_call`-on-leaf sites beyond the ~28 estimate.** The count comes
  from a grep heuristic. Mitigation: after removing the member, any remaining
  reference to `error_message_` is a compile error that pinpoints a missed site.
- **`background` struct error field declaration.** The grep for the field
  declaration matched 3 sites; the single `pba->error_message` use implies a
  declaration that may be formatted differently. Confirm and remove it during
  implementation.
