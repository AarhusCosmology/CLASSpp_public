# Finish ErrorMsg Removal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Delete the vestigial `ErrorMsg` parameter and its scratch buffers from the `tools/` numerical layer and the function-pointer contracts that reach it, so all error reporting is via thrown `std::runtime_error`.

**Architecture:** `tools/*.cpp` already compile as C++, so `class_test`/`class_stop` inside them already throw — the threaded `ErrorMsg` buffers are dead. Each task changes one leaf family's signatures *together with all its call sites* so the build stays green after every task. Function-pointer families (integrator, evolvers) change atomically with their callbacks. Genuine control-flow `int` returns (`get_CF1`, `ludcmp`, `qm_auto`) and the `ErrorMsg` *type* are kept.

**Tech Stack:** C++17, CLASS Makefile (`make class`) + Cython wrapper (`pip install .`), pytest nose_tests in `python/`.

**Spec:** `docs/superpowers/specs/2026-06-11-finish-errormsg-removal-design.md`

---

## Conventions used in every task

**The two call-site transformations** (apply throughout):

```cpp
// BEFORE                                   // AFTER
ErrorMsg buf;                               (delete the local)
class_call_failure(leaf(args, buf), buf);   leaf(args);

// BEFORE                                   // AFTER
class_call(leaf(args, errmsg), errmsg, errmsg);   leaf(args);
```

`leaf(...)` throws on error (C++ build), so no wrapper is needed. Functions keep their `int` return; the ignored return value is fine (no `[[nodiscard]]`, no `-Werror` — see Makefile, only `-g -fPIC -std=c++17 -O3`).

**Remove the vestigial scope block + reformat (apply throughout):** the `#303` idiom wrapped each call in an explicit `{ ErrorMsg buf; class_call_failure(...); }` block to scope the local. Once `buf` is deleted and the block holds a single statement with no declarations, **remove the block's braces and dedent** so the bare call sits in its parent scope (`if (cond) { { foo(); } }` → `if (cond) { foo(); }`). Leave any block that still has a declaration or 2+ statements. Then run `clang-format -i` on every file touched (master is clang-format-conformant; the strip misaligns continuation args, so reformatting is required). Verify conformance with `clang-format --dry-run <file>` (zero warnings). `.clang-format` is NOT modified (`RemoveBracesLLVM` was evaluated and rejected — it doesn't remove standalone blocks and churns ~37 unrelated files).

**Unused-local note:** A scope may call several leaves sharing one `ErrorMsg buf;`/`ErrorMsg errmsg;`. After converting one leaf, the local may still be used by another not-yet-converted leaf — leave it until its last user is gone. Leftover unused locals are harmless warnings (no `-Werror`); a final sweep (Task 11) removes any stragglers.

**Per-task build gate:**
```bash
cd /Users/au192734/Projects/class_claude && make class -j4
```
Expected: links `class` with no errors. (Tasks touching `classy.pyx` additionally run `pip install .`.)

**Do NOT touch:** the `ErrorMsg` *type* in `common.h`; bare control-flow `return _FAILURE_` in `get_CF1`, `ludcmp`, `sp_ludcmp`, `qm_auto`, internal evolver status; `hyrec/*.c`; `cclassy.pxd` (generated).

---

### Task 0: Baseline — confirm green + capture bit-identical reference

**Files:** none modified.

- [ ] **Step 1: Confirm clean build on the branch**

```bash
cd /Users/au192734/Projects/class_claude && git checkout finish-errormsg-removal && make class -j4
```
Expected: builds `class` with no errors (branch is functionally identical to master).

- [ ] **Step 2: Capture reference outputs for the bit-identical check**

```bash
cd /Users/au192734/Projects/class_claude
./class explanatory.ini
mkdir -p python/baseline_ref && cp output/explanatory*.dat python/baseline_ref/
ls python/baseline_ref/
```
Expected: `explanatory*.dat` files copied. (`python/baseline*/` is gitignored — never committed.)

- [ ] **Step 3: Confirm the lvl1 nose_tests pass as a starting point**

```bash
cd /Users/au192734/Projects/class_claude && pip install . -q
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py
```
Expected: all pass (this is the green baseline to preserve).

- [ ] **Step 4: No commit** (nothing changed).

---

### Task 1: `arrays` layer

**Files:**
- Modify: `include/arrays.h` (all `ErrorMsg errmsg` params), `tools/arrays.cpp` (signatures + bodies)
- Modify call sites: `source/background_module.cpp`, `source/lensing_module.cpp`, `source/nonlinear_module.cpp`, and any other `array_*` callers.

- [ ] **Step 1: Enumerate the surface**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg" include/arrays.h
git grep -n "array_" -- 'source/*' 'species/*' 'tools/*' | grep -E "class_call_failure|class_call\(" 
```
This lists every signature and every call site to convert.

- [ ] **Step 2: Drop `ErrorMsg errmsg` from `include/arrays.h`**

Remove the trailing `ErrorMsg errmsg` parameter (and the preceding comma) from every `array_*` declaration. Example:

```cpp
// BEFORE
int array_search_bisect(int n_lines, double* array, double c, int* last_index, ErrorMsg errmsg);
// AFTER
int array_search_bisect(int n_lines, double* array, double c, int* last_index);
```

- [ ] **Step 3: Drop the param from every body in `tools/arrays.cpp`**

Remove `ErrorMsg errmsg` from each function definition's parameter list. The `class_test`/`class_stop` inside already throw (they ignore `errmsg`); the closing `return _SUCCESS_;` stays. If any internal call passes `errmsg` onward to another `array_*`, drop that argument too.

- [ ] **Step 4: Convert all call sites**

Apply the two transformations from Conventions to each `class_call_failure(array_*(...), buf)` and `class_call(array_*(...), errmsg, errmsg)`. Delete each `ErrorMsg buf;`/`errmsg;` local that no longer has any user in its scope (re-grep the enclosing function).

- [ ] **Step 5: Build**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4
```
Expected: no errors. If "too few arguments to function `array_*`" appears, a call site was missed — re-grep and fix.

- [ ] **Step 6: Fast functional check**

```bash
cd /Users/au192734/Projects/class_claude && ./class explanatory.ini && diff -q output/explanatory_background.dat python/baseline_ref/explanatory_background.dat
```
Expected: files identical (no diff output).

- [ ] **Step 7: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add include/arrays.h tools/arrays.cpp source/
git commit -m "Drop ErrorMsg param from arrays leaf layer (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: `parser` layer

**Files:**
- Modify: `include/parser.h`, `tools/parser.cpp`
- Modify call sites: `source/input_module.cpp` (the `class_call(parser_read_*(..., errmsg), errmsg, errmsg)` sites), any other `parser_*` callers.

- [ ] **Step 1: Enumerate**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg" include/parser.h
git grep -n "parser_" -- 'source/*' 'tools/*' | grep -E "class_call"
```

- [ ] **Step 2: Drop `ErrorMsg errmsg` from `include/parser.h` and `tools/parser.cpp`** (signatures + bodies), as in Task 1. The `flag` out-param (found/not-found) is unchanged — only the error buffer goes.

- [ ] **Step 3: Convert call sites** in `source/input_module.cpp` using the two transformations. Note several functions there open with `ErrorMsg errmsg;` shared across many `parser_read_*` calls — remove that local only once its last user in the function is converted.

- [ ] **Step 4: Build**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4
```
Expected: no errors.

- [ ] **Step 5: Functional check + commit**

```bash
cd /Users/au192734/Projects/class_claude
./class explanatory.ini && diff -q output/explanatory_background.dat python/baseline_ref/explanatory_background.dat
git add include/parser.h tools/parser.cpp source/input_module.cpp
git commit -m "Drop ErrorMsg param from parser leaf layer (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: `quadrature` layer

**Files:**
- Modify: `include/quadrature.h`, `tools/quadrature.cpp`
- Modify call sites: `source/lensing_module.cpp`, `species/ncdm_base_species.cpp`, others calling `quadrature_*`/`get_qsampling*`.

- [ ] **Step 1: Enumerate**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg" include/quadrature.h
git grep -n "quadrature_\|get_qsampling" -- 'source/*' 'species/*' 'tools/*' | grep -E "class_call"
```

- [ ] **Step 2: Drop `ErrorMsg errmsg`/`error_message` from signatures + bodies.** **Keep** `get_qsampling_manual` returning `int` — its `case (qm_auto): return _FAILURE_;` is a control-flow sentinel, not an error. Only the `errmsg` parameter is removed; the `return _FAILURE_` stays and its internal caller's return-value check (if any) is unchanged.

- [ ] **Step 3: Convert call sites** (two transformations).

- [ ] **Step 4: Build + functional check + commit**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4 && ./class explanatory.ini && diff -q output/explanatory_background.dat python/baseline_ref/explanatory_background.dat
git add include/quadrature.h tools/quadrature.cpp source/lensing_module.cpp species/ncdm_base_species.cpp
git commit -m "Drop ErrorMsg param from quadrature leaf layer (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: `trigonometric_integrals` layer

**Files:**
- Modify: `include/trigonometric_integrals.h`, `tools/trigonometric_integrals.cpp`
- Modify call sites: grep below.

- [ ] **Step 1: Enumerate + convert**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg" include/trigonometric_integrals.h
git grep -n "cosine_integral\|sine_integral\|trig" -- 'source/*' 'tools/*' | grep -E "class_call"
```
Drop the param from signatures + bodies; convert call sites (two transformations).

- [ ] **Step 2: Build + functional check + commit**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4 && ./class explanatory.ini && diff -q output/explanatory_cl.dat python/baseline_ref/explanatory_cl.dat
git add include/trigonometric_integrals.h tools/trigonometric_integrals.cpp source/
git commit -m "Drop ErrorMsg param from trigonometric_integrals (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: `hyperspherical` layer (incl. classy binding)

**Files:**
- Modify: `include/hyperspherical.h`, `tools/hyperspherical.cpp`
- Modify call sites: `source/transfer_module.cpp`
- Modify binding: `classy.pyx`

- [ ] **Step 1: Enumerate**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg error_message" include/hyperspherical.h
git grep -n "hyperspherical_" -- 'source/*' 'tools/*' | grep -E "class_call"
grep -n "hyperspherical_\|errmsg\|error_message" classy.pyx
```

- [ ] **Step 2: Drop `ErrorMsg error_message` from signatures + bodies.** **Keep** all 9 control-flow `return _FAILURE_` sites in `tools/hyperspherical.cpp` (e.g. `get_CF1`, range checks) — they carry no message and are not errors. Only the parameter and the `class_test`/`class_stop`-driven error paths change (those already throw).

- [ ] **Step 3: Convert C++ call sites** in `source/transfer_module.cpp` (two transformations).

- [ ] **Step 4: Update `classy.pyx` extern declarations** — drop `char* error_message` and add `except +` (the project's standard C++→Python exception bridge, matching the generated `cclassy.pxd` methods):

```cython
    int hyperspherical_HIS_create(int K, double beta, int nl, int* lvec,
                                  double xmin, double xmax, double sampling,
                                  int l_WKB, double phiminabs,
                                  HyperInterpStruct* pHIS) except +
    int hyperspherical_HIS_free(HyperInterpStruct* pHIS) except +
    int hyperspherical_bessel_direct_vector(int K, double beta, int* lvec, int nl,
                                            double* xvec, int nx, double* Phi) except +
```

- [ ] **Step 5: Update `classy.pyx` call sites** — remove the `cdef char errmsg[2048]` locals and the `!= 0 ... raise CosmoSevereError(errmsg)` checks; call the functions directly (a thrown C++ error now surfaces as a Python exception via `except +`):

```cython
    # direct (was lines ~199-202)
    hyperspherical_bessel_direct_vector(Kc, betac, &lvec[0], nl, &xvec[0], nx, &Phi[0, 0])

    # interpolate (was lines ~286-306): drop errmsg, drop the != 0 check
    hyperspherical_HIS_create(Kc, betac, nl, &luniq[0], xmin, xmax, samp, l_WKB, phiminabs, his)
    ...
    finally:
        hyperspherical_HIS_free(his)
        del his
```

- [ ] **Step 6: Build C++ and Cython**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4 && pip install . -q && python -c "import classy"
```
Expected: both succeed.

- [ ] **Step 7: Hyperspherical wrapper tests**

```bash
cd /Users/au192734/Projects/class_claude/python && python -m pytest -q test_hyperspherical.py
```
Expected: all pass (covers create/free/direct/interpolate + the error-raising paths).

- [ ] **Step 8: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add include/hyperspherical.h tools/hyperspherical.cpp source/transfer_module.cpp classy.pyx
git commit -m "Drop ErrorMsg param from hyperspherical layer + classy bindings (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: `generic_integrator` (dei_rkck) family

**Files:**
- Modify: `include/dei_rkck.h` (typedef + signatures + `generic_integrator_workspace::error_message` member), `tools/dei_rkck.cpp`
- Modify callbacks + call sites: `source/primordial_module.cpp` / `.h` (inflation derivs + `generic_integrator` calls), any other `generic_integrator` users.

- [ ] **Step 1: Enumerate**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg\|error_message" include/dei_rkck.h
git grep -n "generic_integrator\|initialize_generic_integrator\|cleanup_generic_integrator\|rkqs\|rkck" -- 'source/*' 'tools/*'
```

- [ ] **Step 2: Edit `include/dei_rkck.h` atomically:**
  - Drop trailing `ErrorMsg error_message` from the `derivs` typedef inside `generic_integrator` / `rkqs` / `rkck`.
  - Drop it from the function signatures themselves.
  - Delete the `ErrorMsg error_message;` member from `struct generic_integrator_workspace` (and its doc comment).

- [ ] **Step 3: Edit `tools/dei_rkck.cpp`:** drop the param from bodies; remove any `pgi->error_message` writes; change internal `class_call((*derivs)(..., pgi->error_message), ...)` to bare `(*derivs)(...);`.

- [ ] **Step 4: Update the inflation derivs callback** in `source/primordial_module.cpp`/`.h` (the function passed to `generic_integrator`) to drop its `ErrorMsg error_message` param, matching the typedef. Convert the `class_call(generic_integrator(...), ...)` / `class_call_failure(...)` call sites (two transformations).

- [ ] **Step 5: Build**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4
```
Expected: no errors. A signature/typedef mismatch shows as "cannot convert ... to int (*)(...)" — fix the callback or call site flagged.

- [ ] **Step 6: Functional check + commit**

```bash
cd /Users/au192734/Projects/class_claude && ./class explanatory.ini && diff -q output/explanatory_cl.dat python/baseline_ref/explanatory_cl.dat
git add include/dei_rkck.h tools/dei_rkck.cpp source/primordial_module.cpp source/primordial_module.h
git commit -m "Drop ErrorMsg from generic_integrator + workspace member (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 7: Evolver family (ndf15 / rkck / rkdp45 + callbacks)

The three evolvers share one signature (assigned to one `auto generic_evolver` pointer in `background_module.cpp:881` and `perturbations_module.cpp:2696`), so they and every callback they receive must change together.

**Files:**
- Modify typedefs + signatures: `include/evolver_ndf15.h`, `include/evolver_rkck.h`, `include/evolver_rkdp45.h`
- Modify impls: `tools/evolver_ndf15.cpp`, `tools/evolver_rkck.cpp`, `tools/evolver_rkdp45.cpp`
- Modify callbacks + decls: `source/background_module.cpp`/`.h`, `source/perturbations_module.cpp`/`.h`, `source/thermodynamics_module.cpp`/`.h`, `source/primordial_module.cpp`/`.h`, `source/input_module.cpp`/`.h` (the `fzero_Newton` shooting callback), `species/thermo_context.h`.

- [ ] **Step 1: Enumerate every callback and call site**

```bash
cd /Users/au192734/Projects/class_claude
grep -n "ErrorMsg" include/evolver_ndf15.h include/evolver_rkck.h include/evolver_rkdp45.h
git grep -n "ErrorMsg error_message)" -- 'source/*.h' 'source/*.cpp' 'species/*'
git grep -n "generic_evolver\|evolver_ndf15\|evolver_rk\|evolver_rkdp45\|fzero_Newton\|numjac" -- 'source/*' 'tools/*'
```

- [ ] **Step 2: Edit the three headers** — drop trailing `ErrorMsg error_message` from:
  - the `derivs`, `output`, `timescale_and_approximation`, `print_variables` callback typedefs inside each evolver signature,
  - the evolver function signatures themselves (`evolver_ndf15`/`evolver_rk`/`evolver_rkdp45`),
  - the helper signatures in `evolver_ndf15.h`: `numjac`, `new_linearisation`, `initialize_jacobian`, `initialize_numjac_workspace`, `fzero_Newton` (both its own param and its `func` callback typedef).

- [ ] **Step 3: Edit the three `tools/evolver_*.cpp`** — drop the param from bodies; change every internal `class_call((*derivs)(..., error_message), ...)` / `class_call(numjac(..., error_message), ...)` to bare calls. Keep `ludcmp`'s `return _FAILURE_` (singular matrix) and `evolver_ndf15.cpp:1097` control flow untouched.

- [ ] **Step 4: Update every callback definition + declaration** to drop `ErrorMsg error_message`, matching the typedefs, in:
  - `background_module` (derivs/output for `background_solve`),
  - `perturbations_module` (`perturb_derivs`, `perturb_sources`, `perturb_print_variables`, `perturb_timescale`, `perturb_tca_slip_and_shear`, etc. — all in the `.h` lines ~285–310),
  - `thermodynamics_module` (recombination derivs/output),
  - `primordial_module` (any evolver-fed callback),
  - `input_module` (`fzero_Newton` shooting target callback, `.h:71`),
  - `species/thermo_context.h` (the 1 `ErrorMsg`).

- [ ] **Step 5: Convert the `generic_evolver`/`fzero_Newton`/`numjac` call sites** in the modules (two transformations), e.g. `class_call_failure(generic_evolver(background_derivs_loga, ...), buf)` → `generic_evolver(background_derivs_loga, ...);`.

- [ ] **Step 6: Build**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4
```
Expected: no errors. Function-pointer mismatch errors point directly at any callback/typedef still out of sync.

- [ ] **Step 7: Functional check (Cl + Pk, since this is the perturbation hot path)**

```bash
cd /Users/au192734/Projects/class_claude && ./class explanatory.ini
for f in python/baseline_ref/*.dat; do diff -q "output/$(basename $f)" "$f"; done
```
Expected: no diffs (bit-identical — error plumbing only).

- [ ] **Step 8: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add include/evolver_*.h tools/evolver_*.cpp source/ species/thermo_context.h
git commit -m "Drop ErrorMsg from evolver function-pointer contracts + callbacks (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 8: Remove `class_call_except` and `class_test_except`

**Files:**
- Modify: `source/primordial_module.cpp`, `source/thermodynamics_module.cpp`, `source/nonlinear_module.cpp`
- Modify: `include/common.h` (delete both macros, C and C++ blocks)

- [ ] **Step 1: Convert the 6 `class_call_except` sites**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "class_call_except(" -- source/
```
- `thermodynamics_module.cpp:450`, `:455` (empty cleanup) → bare call:
```cpp
thermodynamics_recombination(preco, pvecback.data());
...
thermodynamics_reionization(preco, preio, pvecback.data());
```
- `primordial_module.cpp:274`, `:351`, `:371` (`primordial_free()` cleanup) → try/catch:
```cpp
try {
  primordial_analytic_spectrum_init();   // (resp. solve_inflation / external_spectrum_init)
} catch (...) {
  primordial_free();
  throw;
}
```
- `primordial_module.cpp:1799` (`cleanup_generic_integrator(&gi)` cleanup) → try/catch:
```cpp
try {
  primordial_inflation_get_epsilon(y[index_in_phi_], &epsilon);
} catch (...) {
  cleanup_generic_integrator(&gi);
  throw;
}
```

- [ ] **Step 2: Convert the 4 `class_test_except` sites**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "class_test_except(" -- source/
```
Rewrite each `class_test_except(cond, cleanup, "msg", args)` as a guarded throw that runs the cleanup first:
```cpp
if (cond) {
  cleanup;          // e.g. cleanup_generic_integrator(&gi);  (omit if empty)
  class_stop("msg", args);
}
```
Sites: `nonlinear_module.cpp:2328`, `primordial_module.cpp:312`, `:1802`, `thermodynamics_module.cpp:2343`. For `primordial_module.cpp:1802` the cleanup is `cleanup_generic_integrator(&gi)`, matching its sibling from Step 1.

- [ ] **Step 3: Delete both macros from `include/common.h`**

Remove the C definitions of `class_call_except` (lines ~157–167) and `class_test_except` (~208–215), and the C++ `#undef`+redefine blocks for both (~281–291 and ~317–326). Leave `class_call`, `class_test`, `class_stop`, `class_call_try`, `class_open` intact for now.

- [ ] **Step 4: Build**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4
```
Expected: no errors, and no remaining references:
```bash
git grep -n "class_call_except\|class_test_except" -- source/ species/ tools/ include/
```
Expected: only zero matches (macros gone, no callers).

- [ ] **Step 5: Functional check + commit**

```bash
cd /Users/au192734/Projects/class_claude && ./class explanatory.ini && for f in python/baseline_ref/*.dat; do diff -q "output/$(basename $f)" "$f"; done
git add source/primordial_module.cpp source/thermodynamics_module.cpp source/nonlinear_module.cpp include/common.h
git commit -m "Remove class_call_except / class_test_except; inline cleanup-rethrow (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 9: Macro cleanup — delete `class_call_failure`, audit `class_call`/`class_call_try`

**Files:** `include/common.h`, plus any straggler call sites the audit finds.

- [ ] **Step 1: Confirm `class_call_failure` has no users, then delete it**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "class_call_failure" -- source/ species/ tools/ main/
```
Expected: zero matches. Delete the `class_call_failure` definition from `include/common.h` (the `#define class_call_failure ...` block, ~303–315).

- [ ] **Step 2: Audit remaining `class_call(` sites**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "class_call(" -- source/ species/ tools/ | grep -vE "class_call_|class_calls"
```
For each remaining site: it now wraps a C++ function that throws on error. Replace `class_call(fn(...), msg, msg)` with bare `fn(...);` and drop any leftover buffer. (Expected residuals are member-call wrappers like `background_module.cpp:1590`.)

- [ ] **Step 3: Audit `class_call_try`**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "class_call_try" -- source/ species/ tools/ main/
```
If zero matches, delete its definition from `common.h`. If it has users, leave it.

- [ ] **Step 4: If `class_call` itself now has zero users, delete it too**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "class_call(" -- source/ species/ tools/ main/ | grep -vE "class_call_|class_calls"
```
If zero matches, remove the C and C++ `class_call` definitions from `common.h`. Otherwise leave the throwing C++ form.

- [ ] **Step 5: Build C++ + Cython**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4 && pip install . -q && python -c "import classy"
```
Expected: both succeed.

- [ ] **Step 6: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add include/common.h source/ species/ tools/
git commit -m "Remove class_call_failure and any now-dead class_call/class_call_try (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 10: Sweep residual `ErrorMsg` references

**Files:** any with leftover `ErrorMsg` params/locals; `generate_wrapper.py`.

- [ ] **Step 1: Find every remaining `ErrorMsg` reference**

```bash
cd /Users/au192734/Projects/class_claude && git grep -n "ErrorMsg" -- source/ species/ tools/ include/ main/ generate_wrapper.py
```
Expected remaining (legitimate): the `ErrorMsg` typedef in `common.h`; `ErrorMsg _class_err_`/`_class_opt_args_`/`FMsg`/`Optional_arguments` inside the throwing macros; `tools/exceptions.*`; `generate_wrapper.py`'s ctypedef/allowed-type. Anything else is a missed param or unused local — remove it.

- [ ] **Step 2: Remove now-unused local buffers** flagged by the grep (e.g. an `ErrorMsg buf;`/`errmsg;`/`errmsg = "";` with no remaining user — confirm per scope, then delete).

- [ ] **Step 3: `generate_wrapper.py`** — the `ErrorMsg` ctypedef (line 123) and `allowed_types` entry (line 260) stay (the type still exists). Only touch if regeneration breaks (next step).

- [ ] **Step 4: Build C++ + regenerate/build Cython**

```bash
cd /Users/au192734/Projects/class_claude && make class -j4 && pip install . -q && python -c "import classy"
```
Expected: both succeed (this re-runs `generate_wrapper.py` → `cclassy.pxd`).

- [ ] **Step 5: Commit (if anything changed)**

```bash
cd /Users/au192734/Projects/class_claude
git add -A
git commit -m "Sweep residual ErrorMsg locals/params (#38 follow-up)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 11: Final verification

**Files:** none modified (verification + cleanup of baseline dir).

- [ ] **Step 1: Clean rebuild, both manifests**

```bash
cd /Users/au192734/Projects/class_claude && make clean && make class -j4 && pip install . -q && python -c "import classy; print('classy OK')"
```
Expected: clean build, import OK.

- [ ] **Step 2: Full lvl1 nose_tests + greybody + hyperspherical**

```bash
cd /Users/au192734/Projects/class_claude/python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py
python -m pytest -v test_greybody.py
python -m pytest -v test_hyperspherical.py
```
Expected: all pass (matches Task 0 baseline — 96 scenario tests green).

- [ ] **Step 3: Bit-identical confirmation against the captured reference**

```bash
cd /Users/au192734/Projects/class_claude && ./class explanatory.ini
for f in python/baseline_ref/*.dat; do diff -q "output/$(basename $f)" "$f" || echo "DRIFT in $(basename $f)"; done
```
Expected: no "DRIFT" lines (pure plumbing change → bit-identical). If drift appears, apply the ~0.1% / Cl^TE-zero-crossing rule from the spec before proceeding.

- [ ] **Step 4: Verify the Xcode manifest still references all files** (no added/removed source files, so it should be untouched — confirm no `tools/` file was renamed):

```bash
cd /Users/au192734/Projects/class_claude && git diff --name-status master --diff-filter=AD
```
Expected: only the spec + plan added (`A`), no deletions/renames of source files.

- [ ] **Step 5: Final confirmation grep**

```bash
cd /Users/au192734/Projects/class_claude
git grep -n "class_call_failure\|class_call_except\|class_test_except" -- source/ species/ tools/ include/ main/
```
Expected: zero matches.

- [ ] **Step 6: Open the PR**

```bash
cd /Users/au192734/Projects/class_claude && git push -u origin finish-errormsg-removal
gh pr create --title "Finish ErrorMsg removal: error reporting is exclusively via exceptions" \
  --body "$(cat <<'EOF'
Follow-up to #303 (closed #38). Removes the vestigial \`ErrorMsg\` parameter
threaded through the \`tools/\` numerical layer and the function-pointer
contracts that reach it, so error reporting is exclusively via thrown
exceptions. Also removes \`class_call_failure\`, \`class_call_except\`, and
\`class_test_except\` (cleanup-rethrow inlined as try/catch).

Kept deliberately: the \`ErrorMsg\` type (used by throwing macros for stack
locals), genuine control-flow \`int\` returns (get_CF1, ludcmp/sp_ludcmp,
qm_auto), and \`int\` return types. HyRec untouched (self-contained C).

Pure error-path plumbing — verified bit-identical, lvl1 nose_tests +
greybody + hyperspherical all green.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Self-review

- **Spec coverage:** layer 1 (tools sigs) → Tasks 1–7; layer 2 (typedefs) → Tasks 6–7; layer 3 (workspace member) → Task 6; layer 4 (callbacks) → Tasks 6–7; layer 5 (call sites) → every task + Task 9 audit; layer 6 (`*_except` removal) → Task 8; layer 7 (macro cleanup) → Task 9; layer 8 (bindings) → Task 5 (classy.pyx) + Task 10 (generate_wrapper). Non-goals (keep type, control-flow returns, int returns, HyRec) called out in each relevant task. Verification → Tasks 0 + 11. **No gaps.**
- **Placeholder scan:** transformations and cleanup patterns are given as concrete code; high-fan-out mechanical sites are located by exact grep + a stated rule + representative example (appropriate for an identical-pattern sweep). No TBD/TODO.
- **Type consistency:** `class_call_failure`/`class_call`/`class_call_except`/`class_test_except` names used consistently; `generic_evolver`/`generic_integrator`/`fzero_Newton`/`numjac` match the headers inspected.
