# Issue #311 — Delete Vestigial C-era API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove C→C++ transition leftovers (empty `*_free()`, lying `int` returns, `short`/`_TRUE_`/`_FALSE_` flags, stale Doxygen) and complete `HyperInterpStruct` RAII, in one PR closing #311.

**Architecture:** Four sequential commits on branch `issue-311-vestigial-c-api`. This is a pure mechanical refactor — no float math changes — so the load-bearing test for every commit is **byte-identical** CLASS output versus a pre-change baseline, captured across flat + open + closed curvature (the HIS path is curvature-gated). The compiler is the second safety net: misclassifying a callback breaks the build loudly.

**Tech Stack:** C++ (CLASS modules in `source/`, tools in `tools/`, headers in `include/`), Cython (`classy.pyx`), CMake build via `make` shims (`make class`, `make classy`, `make test`), pytest (`python/test_class.py`).

**Spec:** `docs/superpowers/specs/2026-06-15-issue-311-delete-vestigial-c-api-design.md`

---

## Task 0: Build the bit-identical verification harness

This harness is reused by every later task. It captures a baseline from the current (pre-change) tree, then each commit re-runs and compares byte-for-byte.

**Files:**
- Create: `/tmp/bitcheck/flat.ini`, `/tmp/bitcheck/open.ini`, `/tmp/bitcheck/closed.ini` (scratch, not committed)

- [ ] **Step 1: Confirm a clean tree at the branch tip**

Run: `git status --short && git rev-parse --abbrev-ref HEAD`
Expected: working tree clean (only untracked unrelated files), branch `issue-311-vestigial-c-api`.

- [ ] **Step 2: Create three scratch ini files exercising flat + curved HIS paths**

`/tmp/bitcheck/flat.ini`:
```
output = tCl,pCl,lCl,mPk
root = /tmp/bitcheck/out_flat_
l_max_scalars = 1000
P_k_max_h/Mpc = 1.
write warnings = no
```
`/tmp/bitcheck/open.ini` — identical but add `Omega_k = 0.01` and `root = /tmp/bitcheck/out_open_`.
`/tmp/bitcheck/closed.ini` — identical but add `Omega_k = -0.01` and `root = /tmp/bitcheck/out_closed_`.

(Open `Omega_k>0` ⇒ K=−1; closed `Omega_k<0` ⇒ K=+1 — both drive `transfer_update_HIS` through `sgnK != 0`.)

- [ ] **Step 3: Build the current binary and capture the baseline**

Run:
```bash
make class
mkdir -p /tmp/bitcheck/baseline
for s in flat open closed; do ./class /tmp/bitcheck/$s.ini; done
cp /tmp/bitcheck/out_*_*.dat /tmp/bitcheck/baseline/
ls /tmp/bitcheck/baseline/
```
Expected: a non-empty set of `out_{flat,open,closed}_*.dat` files (cl, pk, lensed cl).

- [ ] **Step 4: Define the compare helper (used after every commit)**

The check after each commit is:
```bash
for s in flat open closed; do ./class /tmp/bitcheck/$s.ini; done
for f in /tmp/bitcheck/baseline/*.dat; do
  cmp "$f" "/tmp/bitcheck/$(basename "$f")" || echo "DRIFT: $(basename "$f")";
done
echo "compare done"
```
Expected after a correct refactor commit: `compare done` with **no** `DRIFT:` lines.

No commit in this task (scratch files only).

---

## Task 1: Commit 1 — Complete HIS RAII + delete vestigial `*_free()`

**Files:**
- Modify: `include/hyperspherical.h` (struct + extern "C" block)
- Modify: `tools/hyperspherical.cpp:82` (create→constructor), delete `hyperspherical_HIS_free` at `:300`
- Modify: `source/transfer.h:137-143` (drop `HIS_allocated`)
- Modify: `source/transfer_module.cpp` (`BIS` local, `transfer_update_HIS`, `transfer_workspace_free`, call sites)
- Modify: `classy.pyx:129-133, 281-302`
- Modify: 8 module headers + cpp: delete empty `*_free()` declarations/definitions and call sites

### 1a — Add the HIS constructor, delete the no-op free

- [ ] **Step 1: Declare the constructor inside the struct**

In `include/hyperspherical.h`, inside `typedef struct HypersphericalInterpolationStructure { … }` (after the data members, before the closing `}`), add:
```cpp
  HyperInterpStruct() = default;
  HyperInterpStruct(int K, double beta, int nl, const int* lvec, double xmin,
                    double xmax, double sampling, int l_WKB, double phiminabs);
```
(The defaulted default-constructor is kept because `transfer_workspace` holds `HyperInterpStruct HIS;` as a default-constructed member and classy uses `new HyperInterpStruct()` for the empty case is removed in 1c — but the move-assignment in 1b needs a default-constructible member.)

- [ ] **Step 2: Remove the create/free declarations from the `extern "C"` block**

In `include/hyperspherical.h`, delete the two declarations:
```cpp
int hyperspherical_HIS_create( … HyperInterpStruct* pHIS);
int hyperspherical_HIS_free(HyperInterpStruct* pHIS);
```
Leave every other declaration in that block unchanged.

- [ ] **Step 3: Convert the definition to a constructor**

In `tools/hyperspherical.cpp`, change the definition header at line 82 from:
```cpp
int hyperspherical_HIS_create(int K, double beta, int nl, int* lvec, double xmin,
                              double xmax, double sampling, int l_WKB, double phiminabs,
                              HyperInterpStruct* pHIS) {
```
to:
```cpp
HyperInterpStruct::HyperInterpStruct(int K, double beta, int nl, const int* lvec, double xmin,
                                     double xmax, double sampling, int l_WKB, double phiminabs) {
```
Then in the body, replace every `pHIS->` with direct member access (delete the `pHIS->` prefix), and delete the trailing `return _SUCCESS_;` (a constructor returns nothing). The `class_test(1, "K must be -1, 0, or 1 …")` throw is unchanged and is the legitimate constructor failure path.

- [ ] **Step 4: Delete `hyperspherical_HIS_free`**

In `tools/hyperspherical.cpp`, delete the whole function at line 300 (`int hyperspherical_HIS_free(HyperInterpStruct* pHIS) { … return _SUCCESS_; }`).

### 1b — Transfer module uses RAII

- [ ] **Step 5: Drop the `HIS_allocated` flag**

In `source/transfer.h`, delete the member:
```cpp
  int HIS_allocated; /**< flag specifying whether the previous structure has been allocated */
```
Keep `HyperInterpStruct HIS;` and `HyperInterpStruct* pBIS;` as-is.

- [ ] **Step 6: Construct `BIS` directly**

In `source/transfer_module.cpp`, delete the early declaration `HyperInterpStruct BIS;` (line 137). At the create site (~line 205) replace:
```cpp
  hyperspherical_HIS_create(0, 1., l_size_max_, l_.data(), ppr->hyper_x_min, xmax,
                            ppr->hyper_sampling_flat, l_[l_size_max_ - 1] + 1,
                            ppr->hyper_phi_min_abs, &BIS);
```
with:
```cpp
  HyperInterpStruct BIS(0, 1., l_size_max_, l_.data(), ppr->hyper_x_min, xmax,
                        ppr->hyper_sampling_flat, l_[l_size_max_ - 1] + 1,
                        ppr->hyper_phi_min_abs);
```
(`BIS` is still captured by `&BIS` in the task lambda and passed as `pBIS` — unchanged.) Delete the trailing `hyperspherical_HIS_free(&BIS);` (~line 293).

- [ ] **Step 7: Move-assign in `transfer_update_HIS`**

In `transfer_update_HIS` (`source/transfer_module.cpp:~3505`), delete the top block:
```cpp
  if (ptw->HIS_allocated == _TRUE_) {
    hyperspherical_HIS_free(&(ptw->HIS));
    ptw->HIS_allocated = _FALSE_;
  }
```
At the create site (~line 3583) replace:
```cpp
    hyperspherical_HIS_create(ptw->sgnK, nu, l_size_max, l_.data(), xmin, xmax, sampling,
                              l_[l_size_max - 1] + 1, ppr->hyper_phi_min_abs, &(ptw->HIS));
    ptw->HIS_allocated = _TRUE_;
```
with:
```cpp
    ptw->HIS = HyperInterpStruct(ptw->sgnK, nu, l_size_max, l_.data(), xmin, xmax, sampling,
                                 l_[l_size_max - 1] + 1, ppr->hyper_phi_min_abs);
```
The three read sites (`ptw->HIS.l_size`, `ptw->HIS.chi_at_phimin[...]`, `pHIS = &(ptw->HIS)`) are unchanged.

- [ ] **Step 8: Remove the now-empty `HIS_allocated = _FALSE_` init**

In `transfer_workspace_init` (`source/transfer_module.cpp:~3449`), delete the line `ptw->HIS_allocated = _FALSE_;`.

### 1c — classy.pyx

- [ ] **Step 9: Update the Cython extern block and wrapper**

In `classy.pyx`, in the `cdef extern from "hyperspherical.h"` block, replace the `cppclass HyperInterpStruct: pass` + the two free-function declarations with a constructor declaration:
```cython
    cdef cppclass HyperInterpStruct:
        HyperInterpStruct(int K, double beta, int nl, int* lvec, double xmin,
                          double xmax, double sampling, int l_WKB, double phiminabs) except +
```
Delete the `int hyperspherical_HIS_create(...)` and `int hyperspherical_HIS_free(...)` extern lines (keep `hyperspherical_Hermite_interpolation_vector`).

In the wrapper (~line 282), replace `cdef HyperInterpStruct* his = new HyperInterpStruct()` + the `hyperspherical_HIS_create(...)` call with a single constructing `new`:
```cython
    cdef HyperInterpStruct* his = new HyperInterpStruct(Kc, betac, nl, &luniq[0], xmin, xmax,
                                                        samp, l_WKB, phiminabs)
```
In the `finally:` block, delete `hyperspherical_HIS_free(his)` and keep `del his` (which calls `delete` → destructor). Update the heap-allocation comment above to note the constructor now fills the struct directly.

### 1d — Delete the remaining empty `*_free()` methods

- [ ] **Step 10: Delete the 8 empty module `*_free()` + `background_free_noinput`**

Delete declaration (header) and definition (cpp) for each — none have real call sites:
`background_free` + `background_free_noinput` (background_module.{h,cpp}); `lensing_free`; `nonlinear_free`; `perturb_free`; `primordial_free`; `spectra_free`; `thermodynamics_free`; `transfer_free`.

- [ ] **Step 11: Delete the two empty workspace helpers + their call sites**

`perturb_workspace_free`: delete declaration (`perturbations_module.h:255`), definition (`:2355`), and call site `perturbations_module.cpp:752` (`perturb_workspace_free(index_md, &pw);`).
`nonlinear_hmcode_workspace_free`: delete declaration (`nonlinear_module.h:127`), definition (`:2908`), and call site `nonlinear_module.cpp:1227`.

- [ ] **Step 12: Delete `transfer_workspace_free` + its call site**

Its only real work (the HIS free) is gone via 1b. Delete declaration (`transfer_module.h:198`), definition (`transfer_module.cpp:3489`), and call site (`transfer_module.cpp:283`, `transfer_workspace_free(&ptw);`).

### Verify + commit

- [ ] **Step 13: Build CLI + classy + unit tests**

Run: `make class && make classy && make test`
Expected: all succeed; `python -c "import classy"` works.

- [ ] **Step 14: Bit-identical check (flat + open + closed)**

Run the Task 0 Step-4 compare helper.
Expected: `compare done`, **no** `DRIFT:` lines. (Open/closed are essential here — this commit changes the curvature-gated HIS path.)

- [ ] **Step 15: Commit**

```bash
git add include/hyperspherical.h tools/hyperspherical.cpp source/transfer.h \
        source/transfer_module.cpp classy.pyx source/*_module.h source/*_module.cpp
git commit -m "Complete HIS RAII; delete vestigial *_free() (#311)

Turn HyperInterpStruct into an RAII class (constructor replaces
hyperspherical_HIS_create; delete the no-op hyperspherical_HIS_free).
Transfer workspace move-assigns HIS and drops the HIS_allocated flag.
Delete the eight empty module *_free() and the empty workspace helpers
and their call sites. Output byte-identical (flat+open+closed).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 2: Commit 2 — `int → void` full sweep

Convert every module method that only ever returns `_SUCCESS_` to `void`. Procedure is per-module; one commit at the end.

**Keep `int` (DO NOT convert) — the callback machinery:**
- Static thunks passed as function pointers: `background_derivs_loga`, `background_add_line_to_bg_table`, `background_print_variables`, `perturb_timescale`, `perturb_sources`, `perturb_print_variables`, `perturb_derivs`, `primordial_inflation_derivs`, `thermodynamics_derivs_with_recfast`, `InputModule::ShootingResidual`.
- Their `*_member` partners (invoked via `return …->X_member(...)` inside a thunk): `background_derivs_loga_member`, `background_add_line_to_bg_table_member`, `background_derivs_member`, `perturb_timescale_member`, `perturb_sources_member`, `perturb_print_variables_member`, `perturb_derivs_member`, `primordial_inflation_derivs_member`, `thermodynamics_derivs_with_recfast_member`.
- Any method whose return value is *used* (returned, assigned, or compared). Note: `transfer_get_xmin_generic`-style C-leaf tools and `parser_*` are not module methods and are untouched.

**General rule per candidate method:** if its body's only `return` statements are `return _SUCCESS_;` AND its address is never taken (`grep -rn "&ClassName::method\|(method," source/`) AND it is not in the keep-int list, convert it. When unsure, leave it `int` — the compiler will flag any wrong conversion (typedef mismatch or `return void`), so err toward leaving it.

- [ ] **Step 1: Sweep one module at a time**

For each of `background, thermodynamics, perturbations, primordial, transfer, nonlinear, spectra, lensing` (and `input`, `output`) modules, for every eligible method:
1. In the header: change `int methodName(` → `void methodName(`.
2. In the cpp definition: change `int ClassName::methodName(` → `void ClassName::methodName(`.
3. In the body: delete the trailing `return _SUCCESS_;` and convert any early `return _SUCCESS_;` to bare `return;`.
4. Leave call sites unchanged — they are already bare `methodName();` (verified: zero callers compare `== _SUCCESS_`).

Work module-by-module and rebuild after each module to catch a misclassified callback early:
Run (after each module): `make class`
Expected: compiles. If a `return void` or function-pointer-type error appears, that method is a callback — revert it to `int` and move on.

- [ ] **Step 2: Full build (CLI + classy + tests)**

Run: `make class && make classy && make test`
Expected: all succeed.

- [ ] **Step 3: Bit-identical check**

Run the Task 0 Step-4 compare helper.
Expected: `compare done`, no `DRIFT:` (a signature/return-type change cannot alter output).

- [ ] **Step 4: Commit**

```bash
git add source/
git commit -m "int -> void for module methods that never fail (#311)

These methods always returned _SUCCESS_ while real errors throw, so the
return value was a lie. No caller inspected it. Function-pointer
callbacks (*_derivs/*_timescale/*_sources/*_print_variables and their
_member partners) keep int to match the evolver/quadrature typedefs.
Output byte-identical.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 3: Commit 3 — `short`/`_TRUE_`/`_FALSE_` → `bool`/`true`/`false`

Convert boolean flag members in the C++-only structs. Per-header; one commit.

**Headers (short-member counts):** `background.h` (4), `perturbations.h` (26), `thermodynamics.h` (5), `transfer.h` (6), `primordial.h` (1), `nonlinear.h` (2), `spectra.h` (1), `lensing.h` (2), `species/fluid.h` (3), `species/scalar_field.h` (3).

**PRESERVE `_TRUE_`/`_FALSE_` where they are macro int-arguments** — notably the 48 `class_define_index(index, _TRUE_, …)` calls and any other macro that takes a `short`/int flag. Only struct members that are genuine booleans and their direct reads/writes convert.

- [ ] **Step 1: Convert members header-by-header**

For each header above, for each `short` flag member (e.g. `short has_cls;`, `short initialise_HIS_cache = _FALSE_;`):
1. Change the type `short` → `bool`.
2. Change any in-struct default `= _TRUE_` → `= true`, `= _FALSE_` → `= false`.

- [ ] **Step 2: Convert assignments and comparisons in the cpp files**

For each converted member `X` (e.g. `ppt->has_cls`):
- Assignments: `X = _TRUE_;` → `X = true;`, `X = _FALSE_;` → `X = false;`.
- Comparisons: `X == _TRUE_` → `X`; `X == _FALSE_` → `!X`; `X != _TRUE_` → `!X`; `X != _FALSE_` → `X`.

Find them with: `grep -rn "\b<member>\b" source/ species/`. Do **not** touch `_TRUE_`/`_FALSE_` that are arguments to `class_define_index` or other macros, nor any C-leaf tool usage.

- [ ] **Step 3: Build after each header's members are done**

Run: `make class`
Expected: compiles. A type error pinpoints a missed comparison/assignment site.

- [ ] **Step 4: Full build + bit-identical check**

Run: `make class && make classy && make test`, then the Task 0 Step-4 compare helper.
Expected: all succeed; `compare done`, no `DRIFT:` (`bool` and `short==_TRUE_` evaluate identically).

- [ ] **Step 5: Commit**

```bash
git add source/ species/
git commit -m "short/_TRUE_/_FALSE_ -> bool/true/false in C++ structs (#311)

Boolean flag members in the C++-only structs are now bool. _TRUE_/_FALSE_
remain only as macro arguments (class_define_index) and in C-leaf tools.
Output byte-identical.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 4: Commit 4 — Doxygen lifecycle sweep

Rewrite the module-header Doxygen blocks that describe the C init/free call sequence. Comment-only.

**Exact locations:**
- `background_module.cpp:72,74` — drop the `background_free()` bullet; keep the `background_init()` bullet.
- `lensing_module.cpp:13,15,34-35`
- `perturbations_module.cpp:22,24`
- `primordial_module.cpp:12,14,80`
- `spectra_module.cpp:10,12,231`
- `thermodynamics_module.cpp:70,72`
- `transfer_module.cpp:17,22,62-63`

- [ ] **Step 1: Edit each block**

For each file: remove the `-# X_free() at the end …` bullet and rewrite any "provided that X_init() has been called before, and X_free() has not been called yet" prose to describe the module-object lifetime (the module is constructed and destroyed as an object; results are valid for the lifetime of the module instance). Keep the `X_init()`/dependency-ordering bullets that are still accurate.

- [ ] **Step 2: Build (sanity — comments only)**

Run: `make class`
Expected: compiles (no functional change).

- [ ] **Step 3: Commit**

```bash
git add source/
git commit -m "Sweep stale C init/free lifecycle text from module Doxygen (#311)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 5: Final verification + PR

- [ ] **Step 1: Full python test suite**

Run: `cd python && python -m pytest test_class.py -x -q`
Expected: pass (this exercises the scenario sweep including `Omega_k = ±0.01` curvature and within-tolerance comparisons).

- [ ] **Step 2: Final bit-identical confirmation on the merged result**

Run the Task 0 Step-4 compare helper one more time against the four-commit tip.
Expected: `compare done`, no `DRIFT:`.

- [ ] **Step 3: Push and open the PR**

```bash
git push -u origin issue-311-vestigial-c-api
gh pr create --title "v4 prep: delete vestigial C-era API; complete HIS RAII (closes #311)" \
  --body "$(cat <<'BODY'
Closes #311.

Four commits, each output byte-identical (flat + open + closed curvature):
1. Complete HIS RAII (constructor replaces hyperspherical_HIS_create;
   delete the no-op hyperspherical_HIS_free); delete the empty module/
   workspace *_free() methods and their call sites.
2. int -> void for module methods that only ever returned _SUCCESS_
   (callbacks bound to evolver/quadrature typedefs keep int).
3. short/_TRUE_/_FALSE_ -> bool/true/false in the C++-only structs.
4. Sweep stale C init/free lifecycle text from module Doxygen.

Sweep confirmed HIS was the only remaining C-era create/free pair; no
raw malloc/calloc/class_alloc remain in source/tools/include/species.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Self-review notes

- **Spec coverage:** Part 1 → Task 1 (1a–1d); Part 2 → Task 2; Part 3 → Task 3; Part 4 → Task 4; HIS RAII addition → Task 1 (1a–1c); verification (bit-identical + curvature + pytest) → Task 0 + per-task checks + Task 5. All spec sections map to a task.
- **No placeholders:** every code step shows the exact before/after; large sweeps (Tasks 2–3) give an exact rule + explicit keep-int / preserve-`_TRUE_` lists rather than enumerating ~180/~50 identical edits, with per-module rebuilds as the catch.
- **Type consistency:** the constructor signature `HyperInterpStruct(int, double, int, const int*, double, double, double, int, double)` is identical in `include/hyperspherical.h` (decl), `tools/hyperspherical.cpp` (def), and `classy.pyx` (extern); `HIS` stays a plain member so the three read sites are untouched.
