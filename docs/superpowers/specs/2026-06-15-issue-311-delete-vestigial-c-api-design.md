# Issue #311 — Delete vestigial C-era API (+ complete HIS RAII)

**Date:** 2026-06-15
**Issue:** [#311](https://github.com/AarhusCosmology/CLASSpp/issues/311) — *v4 prep: delete vestigial C-era API (empty `*_free()`, lying `int` returns, `_TRUE_`/`_FALSE_` shorts, stale Doxygen)*
**Label:** refactoring
**Packaging:** one PR closing #311, four commits.

## Goal

Remove the leftovers from the C→C++ transition that a `class_public` reviewer notices in the
first ten minutes of judging the v4 proposal. This is pure mechanical polish: **no float math is
touched**, so every commit must produce **bit-identical** numerical output. Any drift is a bug.

## Scope

In scope (the four parts of #311, plus one RAII completion the user asked for):

1. **Complete HIS RAII + delete vestigial `*_free()`** — turn `HyperInterpStruct` into a proper
   RAII class (constructor replaces the `*_create` free-function; delete the no-op `*_free`), and
   delete the remaining empty module/workspace `*_free()` methods and their call sites.
2. **`int → void` full sweep** — convert every module method that only ever returns `_SUCCESS_`.
3. **`short`/`_TRUE_`/`_FALSE_` → `bool`/`true`/`false`** in the C++-only structs.
4. **Doxygen lifecycle sweep** — rewrite the module-header blocks describing the C init/free call
   sequence.

Out of scope:

- C-linkage leaf tools (`parser_*`, `get_xmin_generic`, the surviving 6 `class_call` sites, the
  `_TRUE_`/`_FALSE_` macros themselves). These keep their C contract.
- Methods bound to evolver/quadrature function-pointer typedefs (must keep `int`).
- Any behavioral/numerical change.

### Sweep finding (RAII across the codebase)

A full sweep of `include/` + `tools/` + `source/` confirms:

- `hyperspherical_HIS_create` / `hyperspherical_HIS_free` is the **only** remaining named
  `*_create`/`*_free` free-function pair.
- There are **zero** raw `malloc`/`calloc`/`realloc` and **zero** `class_alloc` uses in `source/`,
  `tools/`, `include/`, `species/`. The lone textual `free(pvecback); free(integrand_array)` in
  `nonlinear_module.cpp` is inside a commented-out old `class_test_except` block; both variables are
  already `std::vector`. No action.

So "RAII across the codebase" reduces, concretely, to the HIS conversion in commit 1.

## Commit 1 — Complete HIS RAII + delete vestigial `*_free()`

### 1a. `HyperInterpStruct` becomes an RAII class

`HyperInterpStruct` already owns its memory through `std::vector` members, so its destructor is
already correct and implicit. The remaining C-ism is **two-phase init**: `hyperspherical_HIS_create`
fills a default-constructed struct, and `hyperspherical_HIS_free` is a dead no-op.

- **`include/hyperspherical.h`:** add a constructor to the struct:
  `HyperInterpStruct(int K, double beta, int nl, const int* lvec, double xmin, double xmax,
  double sampling, int l_WKB, double phiminabs);`
  Remove the `hyperspherical_HIS_create` and `hyperspherical_HIS_free` declarations from the
  `extern "C"` block (they are only ever called from C++ / Cython-as-C++, so no C linkage is lost).
- **`tools/hyperspherical.cpp`:** move the body of `hyperspherical_HIS_create` into the constructor
  (member access via `this`/direct names instead of `pHIS->`). The single failure path
  (`class_test` on invalid `K`) becomes a throw from the constructor — exception-safe because the
  vectors clean themselves up. Delete `hyperspherical_HIS_free`.

### 1b. Transfer module uses the RAII class

- **`source/transfer.h`:** delete the `HIS_allocated` flag (it exists *only* to gate the now-deleted
  no-op free; it is never used as real logic). Keep `HyperInterpStruct HIS;` as a plain member.
- **`source/transfer_module.cpp`:**
  - The flat `BIS` local becomes a directly-constructed object
    `HyperInterpStruct BIS(0, 1., l_size_max_, l_.data(), …);` (its declaration moves down to where
    `xmax` is known). Delete the trailing `hyperspherical_HIS_free(&BIS)`.
  - In `transfer_update_HIS`, replace create+flag with move-assignment:
    `ptw->HIS = HyperInterpStruct(ptw->sgnK, nu, l_size_max, …);`. Delete the
    free-then-clear-flag block at the top and the `HIS_allocated = _TRUE_` at the bottom.
  - The three read sites (`ptw->HIS.l_size` at 1478, `ptw->HIS.chi_at_phimin[...]` at 2540,
    `pHIS = &(ptw->HIS)` at 3052) are **unchanged** — `HIS` stays a plain member.

  Rationale for plain member + move-assignment over `std::optional<HyperInterpStruct>`: it preserves
  current semantics exactly (the member is repopulated per curved `q`, never observed empty), keeps
  the read sites untouched, and avoids any empty-deref UB risk.

### 1c. classy.pyx

`HyperInterpStruct` is hand-declared as a `cppclass` in `classy.pyx` (not in the auto-generated
`cclassy.pxd`). Update the wrapper:

- Declare the constructor in the `cdef extern from "hyperspherical.h"` block.
- Construct with `new HyperInterpStruct(Kc, betac, nl, &luniq[0], xmin, xmax, samp, l_WKB,
  phiminabs)`; `del his` (calls `delete` → destructor) in `finally`.
- Drop the `hyperspherical_HIS_create` / `hyperspherical_HIS_free` extern declarations and calls.

### 1d. Delete the remaining empty `*_free()` methods

All eight module-level `*_free()` are empty (RAII no-ops) with no real call sites (only Doxygen /
comments). Delete declaration + definition for each:

- `background_free` (+ `background_free_noinput`, its only callee, also empty)
- `lensing_free`, `nonlinear_free`, `perturb_free`, `primordial_free`, `spectra_free`,
  `thermodynamics_free`, `transfer_free`

Two empty workspace helpers have a single call site each — delete method + call site:

- `perturb_workspace_free` (call site `perturbations_module.cpp:752`)
- `nonlinear_hmcode_workspace_free` (call site `nonlinear_module.cpp:1227`)

`transfer_workspace_free` is fully absorbed by 1b (its only real work was the HIS free); delete the
method and its call site (`transfer_module.cpp:283`).

## Commit 2 — `int → void` full sweep

Convert every module method that only ever returns `_SUCCESS_`: signature `int X(...)` → `void
X(...)` in header and definition; drop the early/trailing `return _SUCCESS_;` (a `void` function
falls off the end). Roughly ~180 of the ~199 `int`-returning declarations.

**Exclusions (keep `int`):**

- Methods bound to the evolver/quadrature function-pointer typedefs — they must match
  `int (*)(...)`: the `*_derivs`, `*_timescale*` / `evaluate_timescale`, `*_output`,
  `*_print_variables`, `*_sources` callbacks, and any static thunk passed to `generic_evolver`,
  `quadrature_*`, `transfer_get_lmax`'s `get_xmin_generic`, etc.
- Any method that returns a *real* value (not `_SUCCESS_`).
- C-linkage leaf tools (not module methods anyway).

**Why this is safe:** zero callers check `== _SUCCESS_` / `!= _SUCCESS_`, and the surviving 6
`class_call` sites all wrap C-leaf tools (`get_xmin_generic`, `parser_read_string`), not module
methods. The work is volume, not risk — but the callback exclusions must be classified per method,
not blanket-applied.

## Commit 3 — `short`/`_TRUE_`/`_FALSE_` → `bool`/`true`/`false`

In the C++-only structs (`background.h`, `perturbations.h`, `thermodynamics.h`, `transfer.h`,
`primordial.h`, `nonlinear.h`, `spectra.h`, `lensing.h` and the species headers — all use
`class`/`enum class`, confirming C++-only):

- `short has_*` (and sibling boolean flags) members → `bool`.
- Their assignments `= _TRUE_`/`= _FALSE_` → `= true`/`= false`.
- Their comparisons `flag == _TRUE_` → `flag`; `flag == _FALSE_` → `!flag`.

**Classify per use, not blanket-replace.** Leave `_TRUE_`/`_FALSE_` where they are *macro arguments*
that expect an `int` flag (notably `class_define_index(index, _TRUE_, …)`) and anywhere in the
C-leaf tools. Only the struct members that are genuinely booleans and their direct reads/writes
convert.

## Commit 4 — Doxygen lifecycle sweep

Rewrite the ~11 module-header Doxygen blocks that still describe the C call sequence ("call
`X_init()` … then `X_free()` at the end, when no more calls to … are needed", "`X_free()` has not
been called yet"). Replace with prose that matches the RAII/module-object lifecycle. No code change.

## Verification

Pure mechanical refactor → **bit-identical** output is the bar (this is achievable here and far
stronger than the usual ~0.1% tolerance, because no float operations are added, removed, or
reordered).

**Build:**
- `cmake -S . -B build/cmake && cmake --build build/cmake --parallel` (CLI binary; `make class` is a
  thin shim).
- For commit 1 (touches `classy.pyx`): also `pip install . --no-build-isolation` and confirm the
  `classy` module imports.
- `make test` (C++ unit tests via CTest). Note: some unit tests are NDEBUG-vacuous, so they are a
  weak signal — the bit-identical run below is the real check.

**Bit-identical reference runs (the load-bearing check):**
1. *Before any change*, on a clean HEAD build, run a representative set of `.ini` files and archive
   their output `.dat` files.
2. After each commit, rebuild and re-run the same `.ini` files; `cmp`/`diff` each output file
   against the archived baseline — expect **byte-identical**.
3. The set **must** include curvature-varied scenarios so the HIS path (`sgnK != 0`) is exercised:
   - a flat scenario (`explanatory.ini`),
   - an **open** universe (`Omega_k > 0`, K = −1),
   - a **closed** universe (`Omega_k < 0`, K = +1).
   A flat-only run would not catch a HIS regression.

**Python suite (final gate):** `pytest python/test_class.py` (the scenario sweep). Optionally with
`COMPARE_OUTPUT_REF=1` to compare against `classyref`; the within-tolerance gates
(`COMPARE_CL_RELATIVE_ERROR = 3e-3`, etc.) should pass unchanged.

## Risks

- **Commit 1 (HIS):** the only path that changes generated code is curvature-dependent — mitigated by
  the open/closed bit-identical runs above. The Cython binding must still compile and import.
- **Commit 2 (int→void):** misclassifying a callback method would break the build (typedef
  mismatch) — caught immediately by the compiler, not silently. Low residual risk.
- **Commit 3 (bool):** accidentally converting a `_TRUE_` that is a macro `int` argument would change
  behavior — mitigated by per-use classification and the bit-identical runs.
- **Commit 4 (Doxygen):** comment-only, no risk.

## Definition of done

- Four commits on a feature branch, PR opened that closes #311.
- Clean build (CLI + classy), `make test` green.
- All representative reference runs (flat + open + closed) byte-identical to the pre-change baseline
  at every commit.
- `pytest python/test_class.py` green.
