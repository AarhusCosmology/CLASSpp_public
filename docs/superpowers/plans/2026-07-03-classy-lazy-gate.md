# Lazy Dirty-Check Gate in classy.pyx — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Cython wrapper as lazy as the C++ core: any method call always reflects the current parameter dict, construction no longer runs shooting, and `set()`/`compute()` stay as fully compatible MontePython shims.

**Architecture:** A single `_cosmo()` accessor becomes the only path to the C++ `Cosmology` object; it rebuilds (cheap parse-only `reset()`) when the parameter dict is dirty. The ten cached input-struct pointers (`self.pr`, `self.ba`, …) are replaced by accessor methods that fetch through `GetInputModule()` on every use, because `DoShooting` *replaces* the `InputModule` instance and cached pointers go stale.

**Spec:** `docs/superpowers/specs/2026-07-03-classy-lazy-wrapper-design.md` (read it first).

**Tech Stack:** Cython (`classy.pyx`), CMake via `make classy` (pip install + MontePython layout), pytest + unittest (`python/test_class.py`).

## Global Constraints

- **No C++ changes.** All code changes live in `classy.pyx`; tests in `python/test_class.py`.
- **Never hand-edit `cclassy.pxd`** — it is generated at build time by `generate_wrapper.py`.
- **Never `git add -A` / `git add .`** in this repo (in-source build artifacts get swept in). Stage explicit paths only.
- **No bit-identity assertions** in tests — use relative tolerances.
- **Rebuild command:** `make classy` from the repo root (runs `pip install .` and recreates the `python/build` layout). Every `classy.pyx` change requires a rebuild before tests see it.
- **Test command:** run pytest from inside `python/` so `pytest.ini` is picked up: `cd python && python -m pytest test_class.py -k <filter> -v`.
- Branch: `classy-lazy-gate` (already created; spec committed as `7b1a3e2b`).
- Module ptr typedefs are `std::shared_ptr<const T>` (`include/common.h:34-43`), so Cython locals may copy them freely — the existing `input_module = deref(...)...; deref(input_module).member` idiom is the one to follow.

---

### Task 1: Lifecycle regression tests

**Files:**
- Modify: `python/test_class.py` (append a new test class immediately before the `if __name__ == '__main__':` block at the end of the file)

**Interfaces:**
- Consumes: `classy.Class`, `classy.CosmoSevereError` (already imported at the top of `test_class.py`), `numpy as np`, `unittest` (already imported).
- Produces: test class `TestLazyLifecycle` with six tests. Two of them are decorated `@unittest.expectedFailure` in this task; Task 2 removes those decorators when the implementation makes them pass. Task 2 relies on these exact test names:
  - `test_set_without_compute_returns_fresh_results` (expectedFailure until Task 2)
  - `test_unread_parameter_raises_at_first_access` (expectedFailure until Task 2)
  - `test_methods_work_without_compute`
  - `test_bad_parameter_raises_at_construction`
  - `test_unread_parameter_raises_at_compute_not_set`
  - `test_montepython_flow`

- [ ] **Step 1: Ensure the installed `classy` matches this branch**

Run from the repo root:
```bash
make classy
```
Expected: exits 0, ends with a successful pip install. (The branch is master + one docs commit, so this build represents current behavior.)

- [ ] **Step 2: Append the test class**

Insert the following into `python/test_class.py`, immediately **before** the final `if __name__ == '__main__':` block:

```python
class TestLazyLifecycle(unittest.TestCase):
    """Wrapper lifecycle: methods must always reflect the current parameter
    dict, with construction = validity and set()/compute() as legacy shims.
    See docs/superpowers/specs/2026-07-03-classy-lazy-wrapper-design.md.

    The two @unittest.expectedFailure markers document the pre-gate staleness
    bug; the gate implementation removes them.
    """

    def _fresh_angular_distance(self, h, z=1.0):
        cosmo = Class({'h': h})
        try:
            return cosmo.angular_distance(z)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    @unittest.expectedFailure
    def test_set_without_compute_returns_fresh_results(self):
        # The historic footgun: set() then a method call WITHOUT compute()
        # must return results for the new parameters, not stale ones.
        cosmo = Class({'h': 0.67})
        try:
            d_before = cosmo.angular_distance(1.0)
            cosmo.set({'h': 0.70})
            d_after = cosmo.angular_distance(1.0)  # no compute() in between
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
        d_expected = self._fresh_angular_distance(0.70)
        self.assertNotEqual(d_before, d_after)
        self.assertAlmostEqual(d_after / d_expected, 1.0, places=10)

    @unittest.expectedFailure
    def test_unread_parameter_raises_at_first_access(self):
        # After set() with a bogus parameter, the first access (not just
        # compute()) must surface the input error instead of silently
        # serving results for the previous parameters.
        cosmo = Class()
        try:
            cosmo.set({'this_parameter_does_not_exist': 1})
            with self.assertRaises(CosmoSevereError):
                cosmo.angular_distance(1.0)
        finally:
            cosmo.empty()
            cosmo.struct_cleanup()

    def test_methods_work_without_compute(self):
        # Construction = validity: no compute() call is ever needed.
        cosmo = Class({'output': 'tCl', 'l_max_scalars': 100})
        try:
            cl = cosmo.raw_cl(100)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
        self.assertTrue(np.all(np.isfinite(cl['tt'])))
        self.assertGreater(np.max(np.abs(cl['tt'])), 0.0)

    def test_bad_parameter_raises_at_construction(self):
        # Constructor-style use validates input at construction (parse-time).
        with self.assertRaises(CosmoSevereError):
            Class({'this_parameter_does_not_exist': 1})

    def test_unread_parameter_raises_at_compute_not_set(self):
        # MontePython compatibility: set() never raises; the error surfaces
        # at compute(), inside the caller's try block.
        cosmo = Class()
        try:
            cosmo.set({'this_parameter_does_not_exist': 1})  # must not raise
            with self.assertRaises(CosmoSevereError):
                cosmo.compute()
        finally:
            cosmo.empty()
            cosmo.struct_cleanup()

    def test_montepython_flow(self):
        # The canonical MontePython sequence: repeated set(), compute(),
        # likelihood-style method calls, then a second point.
        cosmo = Class()
        try:
            cosmo.set({'output': 'tCl', 'l_max_scalars': 100})
            cosmo.set({'h': 0.70})
            cosmo.compute()
            cl1 = cosmo.raw_cl(100)
            self.assertTrue(cosmo.state)
            self.assertGreater(np.max(np.abs(cl1['tt'])), 0.0)

            cosmo.struct_cleanup()
            cosmo.set({'h': 0.68})
            cosmo.compute()
            cl2 = cosmo.raw_cl(100)
            # atol=0: raw dimensionless Cl values are ~1e-10, far below
            # np.allclose's default atol=1e-8, which would compare any two
            # Cl arrays as "close".
            self.assertFalse(np.allclose(cl1['tt'][2:], cl2['tt'][2:],
                                         rtol=1e-5, atol=0.0))
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
```

Note the two error-timing tests use `cosmo.empty()` **before** `struct_cleanup()` in the `finally`: after the gate lands, any access on an object whose dict still contains the bogus parameter re-raises; `empty()` clears the dict first. `struct_cleanup()` is a no-op but kept for MontePython-style symmetry.

- [ ] **Step 3: Run the new tests against the current build**

```bash
cd python && python -m pytest test_class.py -k TestLazyLifecycle -v
```
Expected: **4 passed, 2 xfailed** — the two `expectedFailure` tests (`test_set_without_compute_returns_fresh_results`, `test_unread_parameter_raises_at_first_access`) fail against the current wrapper (stale results, no raise) and are reported as XFAIL; the other four pass, demonstrating that today's contract is preserved.

If anything reports as plain FAILED or XPASS, stop and investigate — the premise of the plan is wrong.

- [ ] **Step 4: Commit**

```bash
git add python/test_class.py
git commit -m "Add lifecycle regression tests for classy wrapper laziness

Two tests are marked expectedFailure: they document the set()-without-
compute() staleness bug that the lazy gate (next commit) fixes."
```

---

### Task 2: Implement the `_cosmo()` gate and remove the struct-pointer cache

**Files:**
- Modify: `classy.pyx:346-472` (lifecycle block) plus mechanical global replacements across the file (~170 call sites)
- Modify: `classy.pyx:1853` (CosmoHammer `__call__`)
- Modify: `python/test_class.py` (remove the two `@unittest.expectedFailure` decorators)
- Test: `python/test_class.py::TestLazyLifecycle`

**Interfaces:**
- Consumes: test names from Task 1; `Cosmology`, `GetInputModule()` and struct members `precision_, background_, thermodynamics_, perturbations_, primordial_, nonlinear_, transfers_, spectra_, lensing_, output_` already declared in the generated `cclassy.pxd`.
- Produces (wrapper-internal API used throughout `classy.pyx`):
  - `cdef Cosmology* _cosmo(self) except NULL` — the only sanctioned path to the C++ object
  - `cdef const precision* pr(self) except NULL` and analogous `ba, th, pt, pm, nl, tr, sp, le, op` — struct accessors replacing the former cached pointer members of the same names

**ORDERING WARNING:** Steps 1–5 (manual edits) must be completed **before** Step 6 (global sed). In particular, `reset()`'s `deref(self._thisptr)` line must be *deleted* in Step 2 before the sed rewrites `deref(self._thisptr)` → `deref(self._cosmo())` file-wide — otherwise `reset()` would call `_cosmo()`, which calls `reset()`: infinite recursion.

- [ ] **Step 1: Remove the ten cached pointer members**

In the `cdef class PyCosmology` attribute block (`classy.pyx:346-360`), replace:

```cython
    cdef unique_ptr[Cosmology] _thisptr
    cdef dict _pars
    cdef bool parameters_changed
    cdef FileContent _fc

    cdef const precision* pr
    cdef const background* ba
    cdef const thermo* th
    cdef const perturbs* pt
    cdef const primordial* pm
    cdef const nonlinear* nl
    cdef const transfers* tr
    cdef const spectra* sp
    cdef const lensing* le
    cdef const output* op
```

with:

```cython
    cdef unique_ptr[Cosmology] _thisptr
    cdef dict _pars
    cdef bool parameters_changed
    cdef FileContent _fc
```

- [ ] **Step 2: Slim down `reset()`**

In `reset()` (`classy.pyx:384-414`), replace the tail:

```cython
        input_module = deref(self._thisptr).GetInputModule()
        self.pr = &deref(input_module).precision_
        self.ba = &deref(input_module).background_
        self.th = &deref(input_module).thermodynamics_
        self.pt = &deref(input_module).perturbations_
        self.pm = &deref(input_module).primordial_
        self.nl = &deref(input_module).nonlinear_
        self.tr = &deref(input_module).transfers_
        self.sp = &deref(input_module).spectra_
        self.le = &deref(input_module).lensing_
        self.op = &deref(input_module).output_
        return self
```

with:

```cython
        return self
```

`reset()` now does parse-time work only (rebuild `_fc`, construct `Cosmology`, unread-parameter check) — no shooting. It must never call `_cosmo()` or the struct accessors.

- [ ] **Step 3: Add the gate and the struct accessors**

Insert immediately after `reset()` (before `cdef _update_fc_from_pars`):

```cython
    # Single gate to the C++ object: rebuilds the Cosmology if the parameter
    # dict changed since it was last built. Everything must reach the C++
    # object through _cosmo() (or the struct accessors below), never through
    # _thisptr directly, so results always reflect the current parameters.
    cdef Cosmology* _cosmo(self) except NULL:
        if self.parameters_changed or self._thisptr.get() == NULL:
            self.reset()
        return self._thisptr.get()

    # The input structs live inside the InputModule, and DoShooting replaces
    # the InputModule instance on first use; fetching through GetInputModule()
    # on every access is the only way a cached pointer cannot go stale.
    cdef const precision* pr(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).precision_

    cdef const background* ba(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).background_

    cdef const thermo* th(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).thermodynamics_

    cdef const perturbs* pt(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).perturbations_

    cdef const primordial* pm(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).primordial_

    cdef const nonlinear* nl(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).nonlinear_

    cdef const transfers* tr(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).transfers_

    cdef const spectra* sp(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).spectra_

    cdef const lensing* le(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).lensing_

    cdef const output* op(self) except NULL:
        input_module = deref(self._cosmo()).GetInputModule()
        return &deref(input_module).output_
```

(Name collisions were checked: no existing `def pr(/ba(/th(/pt(/pm(/nl(/tr(/sp(/le(/op(` in `classy.pyx`.)

- [ ] **Step 4: Make `compute()` a pure trigger**

In `compute()` (`classy.pyx:450-472`), replace:

```cython
    cpdef compute(self, level=None):
        if level is None:
            level = ['lensing']
        if self.parameters_changed:
            self.reset()
        final_level = level[0].lower()
```

with:

```cython
    cpdef compute(self, level=None):
        # Legacy shim (MontePython): modules are built lazily on first access,
        # so compute() only pre-triggers the requested level. Its remaining
        # value is that computation errors surface here, where callers expect
        # to catch them.
        if level is None:
            level = ['lensing']
        final_level = level[0].lower()
```

The `deref(self._thisptr).GetXxxModule()` lines inside `compute()` are rewritten by the global sed in Step 6 — leave them alone here.

- [ ] **Step 5: Fix the CosmoHammer `__call__`**

At `classy.pyx:1851-1854`, replace:

```cython
        self._pars = data.cosmo_arguments
        self.reset()
        self.compute()
```

with:

```cython
        self._pars = data.cosmo_arguments
        self.parameters_changed = True
        self.compute()
```

- [ ] **Step 6: Global mechanical replacements**

Run from the repo root (macOS sed):

```bash
sed -i '' 's/deref(self\._thisptr)/deref(self._cosmo())/g' classy.pyx
for s in pr ba th pt pm nl tr sp le op; do
  sed -i '' "s/self\.${s}\./self.${s}()./g" classy.pyx
  sed -i '' "s/deref(self\.${s})/deref(self.${s}())/g" classy.pyx
done
```

The first sed converts ~110 module-getter sites (including the `Omega_nu` property and all of `compute()`). The loop converts ~58 struct-field reads (`self.ba.h` → `self.ba().h`; Cython auto-derefs pointer member access) and the ten `get_input_*` bodies (`deref(self.ba)` → `deref(self.ba())`).

- [ ] **Step 7: Guard greps — verify the rewrite is complete and did not leak**

```bash
grep -n 'deref(self\._thisptr)' classy.pyx
```
Expected: no output.

```bash
grep -n 'self\._thisptr' classy.pyx
```
Expected: exactly 4 lines — the `cdef unique_ptr[Cosmology] _thisptr` declaration, the `self._thisptr.reset(...)` in `reset()`, and the two `self._thisptr.get()` in `_cosmo()`.

```bash
for s in pr ba th pt pm nl tr sp le op; do grep -n "self\.${s}[.)]" classy.pyx | grep -v "self\.${s}()" ; done
```
Expected: no output (every former pointer use now goes through the accessor).

Then review the full diff by eye for collateral damage (docstrings, comments):
```bash
git diff classy.pyx | head -400
```

- [ ] **Step 8: Remove the two `@unittest.expectedFailure` decorators**

In `python/test_class.py::TestLazyLifecycle`, delete the `@unittest.expectedFailure` lines above `test_set_without_compute_returns_fresh_results` and `test_unread_parameter_raises_at_first_access`, and delete the last sentence of the class docstring ("The two @unittest.expectedFailure markers ... removes them.").

- [ ] **Step 9: Rebuild**

```bash
make classy
```
Expected: exits 0. Cython/C++ compile errors here mean a botched edit in Steps 1–6 — fix before proceeding.

- [ ] **Step 10: Run the lifecycle tests**

```bash
cd python && python -m pytest test_class.py -k TestLazyLifecycle -v
```
Expected: **6 passed** (no xfail markers remain).

- [ ] **Step 11: Commit**

```bash
git add classy.pyx python/test_class.py
git commit -m "Route all classy.pyx access through a lazy dirty-check gate

Adds _cosmo(): rebuilds the C++ Cosmology iff the parameter dict changed,
making every method call reflect current parameters (fixes silent staleness
after set() without compute()). Replaces the ten cached input-struct
pointers with accessors fetched through GetInputModule() on each use --
DoShooting replaces the InputModule, so cached pointers forced reset() to
shoot eagerly at construction. reset() is now parse-only; compute() is a
legacy shim that pre-triggers a level so errors keep surfacing where
MontePython catches them."
```

---

### Task 3: Full-suite verification

**Files:** none created or modified (verification only; plan checkboxes updated).

**Interfaces:**
- Consumes: the built wrapper from Task 2.
- Produces: evidence that the legacy contract is intact.

- [ ] **Step 1: Run the wrapper unit-test files**

```bash
cd python && python -m pytest test_greybody.py test_hyperspherical.py -v
```
Expected: all pass (these exercise wrapper entry points unrelated to the lifecycle; a failure means the sed leaked).

- [ ] **Step 2: Run the full scenario matrix at the default test level**

```bash
cd python && python -m pytest test_class.py -v
```
Expected: same pass/fail/skip profile as master (TEST_LEVEL=0 default). This takes a while (parameterized scenario matrix). Known pre-existing failures documented in memory (stale committed goldens under ffast-math) live in `test_background_columns.py`/`test_transfer_columns.py`, which are *not* part of this run — any failure here is new and must be investigated.

- [ ] **Step 3: Reference comparison (only if `classyref` is installed)**

```bash
python -c "import classyref" 2>/dev/null && (cd python && COMPARE_OUTPUT_REF=1 python -m pytest test_class.py -k TestReviewRegressions -v) || echo "classyref not installed - skipping"
```
Expected: pass or explicit skip. This change is wrapper-only (zero numerical code touched), so any numerical drift vs the reference build is a bug.

- [ ] **Step 4: Confirm branch state**

```bash
git log --oneline master..classy-lazy-gate
git status --short
```
Expected: three commits (spec, tests, implementation); no unstaged changes to tracked files.

---

## Self-review notes

- Spec coverage: gate (§1) → Task 2 Steps 3/6; pointer-cache removal (§2) → Task 2 Steps 1/3/6; slim `reset()` (§3) → Task 2 Step 2; lifecycle methods (§4) → Task 2 Steps 4/5 (`__init__`, `set()`, `empty()`, `struct_cleanup()`, `state` intentionally unchanged); error semantics (§5) → Task 1 tests 2/4/5; testing (§Testing) → Tasks 1 and 3.
- The `except NULL` on `_cosmo()` and the accessors lets Python exceptions raised inside `reset()` propagate through pointer-returning cdef methods.
- `_cosmo()` must never be called from `reset()` / `_update_fc_from_pars()` — enforced by Step 2 ordering and the recursion warning.
