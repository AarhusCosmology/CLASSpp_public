# Species-registered Transfer-Source Indices (#309) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let each species register its own transfer-source (`index_tp`) slots via a new `RegisterTransferSourceIndices(int&, const SourceRequestContext&)` hook, dissolving the per-species `index_tp_*` / `has_source_*` members and both `SetSourceSlot` `dynamic_cast` loops from `PerturbationsModule`.

**Architecture:** This is the **background-index pattern** (`RegisterBackgroundIndices(int&)`) applied to a third index space. Transfer-source indices are computed once before the k-loop and are identical for every k, so they live as **plain members on the species** (like `index_bg_rho_`), written once at registration and read read-only thereafter (thread-safe in the parallel k-loop for the same reason background indices are). Source-writing **leaf** species own their `index_tp_delta_`/`index_tp_theta_` members; `CompositeSpecies` forwards registration to `children_`; the two interacting composites keep their channel members on the composite because that is where the `has_idr()`/`has_idm_dr()` gating lives.

**This is a behavior-preserving refactor.** No arithmetic changes — only where the source slot integer is stored and who allocates it. Output must stay **byte-identical**. The internal `index_tp` ordering changes (a species' δ and θ become adjacent instead of living in separate blocks), but this is invisible: Tk/vTk **column order** is produced by the species-iteration column writer in two passes (density section then velocity section), independent of the numeric `index_tp` values. The characterization test in Task 1 locks this invariant.

**Tech Stack:** C++17 (CLASS++ modules in `source/`, species plugins in `species/`); Python `classy` Cython wrapper; `pytest`; build via `pip install . --no-build-isolation`; reference comparison via `classyref` + the `TEST_LEVEL`/`COMPARE_OUTPUT_REF` env knobs in `python/test_class.py`.

**Spec:** `docs/superpowers/specs/2026-06-16-issue-309-species-registered-transfer-source-types-design.md`

---

## Conventions used by every migration task

These are referenced (not repeated) in tasks 3–13. Each task still lists its own exact symbols.

- **Member declarations** (added to the owning species' header, in its `private:` section, default `-1` like `index_bg_rho_`):
  ```cpp
  int index_tp_delta_ = -1;  // transfer-source slot; -1 = not registered
  int index_tp_theta_ = -1;
  ```
  (Composites that own multiple channels use channel-suffixed names, e.g. `index_tp_delta_idr_`; each task gives the exact names.)

- **`class_define_index` semantics** (`include/common.h:294`): writes the index and bumps the running index **only when the condition is true**; leaves the member at `-1` otherwise. So an unrequested source naturally stays `-1`.

- **Repoint rule:** in the species' `FillSources` (writes), `p_mod->index_tp_<x>_` becomes the species' own `index_tp_<x>_` member. In `WriteOutputColumns` (columns), `mod.index_tp_<x>_` becomes the same member. The `active` argument that was `mod.has_source_<x>_` (or a literal) becomes `index_tp_<x>_ >= 0`. Literal `_TRUE_`/`true` actives that encode a *gauge* condition (CDM `theta`) are preserved as-is because the member already encodes that gating.

- **Module deletion rule (symbol-based, drift-proof):** for each migrated source type `<x>`, remove from `source/perturbations_module.h` the `int index_tp_<x>_;` declaration and the `bool has_source_<x>_;` declaration; remove from `source/perturbations_module.cpp` the `has_source_<x>_ = true;` assignment(s) in `perturb_indices_of_perturbs()` and the matching `class_define_index(index_tp_<x>_, has_source_<x>_, index_type, …);` line. Do **not** touch module-owned types (`t0/t1/t2/p`, `delta_m/cb/tot`, `theta_m/cb/tot`, metric/`phi/psi/h/eta/H_T_Nb_prime/k2gamma_Nb`).

- **Per-task verification loop** (the "GREEN GATE"):
  1. `pip install . --no-build-isolation` — expect success, no new warnings on the touched files.
  2. `cd python && python -m pytest test_transfer_columns.py -v` — expect **all PASS** (byte-identical to the Task-1 golden).
  3. Commit.

- **NEVER regenerate the Task-1 golden files after Task 1.** They encode pre-refactor behavior and are the oracle.

---

## File Structure

**New files:**
- `python/test_transfer_columns.py` — characterization (golden) test, the per-task oracle.
- `python/gen_transfer_golden.py` — one-shot generator for the golden `.npz` files.
- `python/transfer_golden/*.npz` — committed golden transfer-function tables (one per cosmology case).
- `species/perturb_source_context.h` — gains `struct SourceRequestContext` (append to the existing file).

**Modified — module:**
- `source/perturbations_module.h` — delete per-species `index_tp_*` + `has_source_*` members (keep module-owned).
- `source/perturbations_module.cpp` — replace the per-species `has_source_*` gating block, the per-species `class_define_index` block, and both `SetSourceSlot` loops with one registration loop.

**Modified — species (one cohesive responsibility each):**
- `species/base_species.h` — new `RegisterTransferSourceIndices` virtual (default no-op).
- `species/composite_species.{h,cpp}` — generic child-forwarding override.
- Leaves: `species/{fluid,ultra_relativistic,scalar_field,cdm,photons,baryons}.{h,cpp}`.
- NCDM family: `species/ncdm_base_species.{h,cpp}` (members + register), `species/ncdm_species.cpp` (repoint + drop `source_slot_`), `species/dncdm_species.*` (registers, no fill).
- DR channel: `species/dark_radiation_species.{h,cpp}` (members + register), `species/dcdm_dr_species.cpp`, `species/dncdm_dr_species.cpp` (repoint to `dr().index_tp_delta_`).
- DCDM channel: `species/dcdm.{h,cpp}` (member + register + own column), `species/dcdm_dr_species.cpp` (repoint fill to `dcdm().index_tp_delta_`).
- IDM composites: `species/idm_dr_idr_species.{h,cpp}`, `species/idm_drmd_idr_drmd_species.{h,cpp}` (composite-owned channel members + register override).

---

## Task 1: Characterization golden test (the safety net)

**Files:**
- Create: `python/gen_transfer_golden.py`
- Create: `python/test_transfer_columns.py`
- Create: `python/transfer_golden/` (golden `.npz` files, generated)

**This task runs on the current branch (spec commit only), so it captures true pre-refactor behavior.**

- [ ] **Step 1: Build the current (pre-change) code**

Run: `pip install . --no-build-isolation`
Expected: builds and installs `classy` successfully.

- [ ] **Step 2: Write the shared case list + flatten helper**

Create `python/gen_transfer_golden.py`:

```python
"""One-shot generator for transfer-function golden tables (issue #309).

Run ONCE on the pre-refactor build:  python gen_transfer_golden.py
Do NOT re-run after the refactor starts — the goldens are the oracle.
"""
import os
import numpy as np
from classy import Class

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "transfer_golden")

# Each case turns on a set of source-writing species and asks for both the
# density (mTk) and velocity (vTk) transfer functions in synchronous gauge.
# Cases are split so mutually-exclusive species don't collide.
COMMON = {
    "output": "mTk,vTk",
    "gauge": "synchronous",
    "l_max_scalars": 200,
    "P_k_max_1/Mpc": 1.0,
    "z_pk": 0.0,
}

CASES = {
    # photons, baryons, cdm, ur, ncdm (×2 flavors)
    "cdm_ur_ncdm": {
        "N_ur": 1.0, "N_ncdm": 2, "m_ncdm": "0.04,0.06", "T_ncdm": "0.71611,0.71611",
    },
    # fluid dark energy
    "fluid": {
        "Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.0, "N_ur": 3.044,
    },
    # scalar field (synchronous only)
    "scf": {
        "Omega_Lambda": 0.0, "Omega_fld": 0.0,
        "Omega_scf": -1, "attractor_ic_scf": "yes",
        "scf_parameters": "10.0, 0.0, 0.0, 0.0, 0.0, 0.0", "N_ur": 3.044,
    },
    # decaying CDM -> dark radiation (DCDM_DR + DR child)
    "dcdm_dr": {
        "Omega_ini_dcdm": 0.01, "Gamma_dcdm": 100.0, "N_ur": 3.044,
    },
    # decaying NCDM -> dark radiation (DNCDM_DR: dncdm child + DR child)
    "dncdm_dr": {
        "N_ncdm": 1, "m_ncdm": 0.06, "Omega_ini_dcdm": 0.0,
        "N_dncdm": 1, "m_dncdm": 1.0, "Gamma_dncdm": 10.0, "N_ur": 2.0,
    },
    # interacting DM + dark radiation
    "idm_dr_idr": {
        "f_idm": 0.1, "a_idm_dr": 1000.0, "nindex_idm_dr": 4,
        "xi_idr": 0.3, "N_ur": 2.0,
    },
    # interacting DM + dark radiation (DRMD variant)
    "idm_drmd_idr_drmd": {
        "f_idm": 0.1, "Gamma_0_drmd": 1.0, "xi_idr": 0.3, "N_ur": 2.0,
    },
}


def flatten(tk):
    """get_transfer() may return {name: array} or {ic: {name: array}}.
    Flatten to {fullkey: 1d-array} with deterministic key order preserved."""
    out = {}
    sample = next(iter(tk.values())) if tk else None
    if isinstance(sample, dict):
        for ic, d in tk.items():
            for name, arr in d.items():
                out[f"{ic}::{name}"] = np.asarray(arr, dtype=np.float64)
    else:
        for name, arr in tk.items():
            out[name] = np.asarray(arr, dtype=np.float64)
    return out


def transfer_for(extra):
    cosmo = Class()
    cosmo.set({**COMMON, **extra})
    cosmo.compute()
    flat = flatten(cosmo.get_transfer())
    cosmo.struct_cleanup()
    cosmo.empty()
    return flat


if __name__ == "__main__":
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    for name, extra in CASES.items():
        flat = transfer_for(extra)
        # store the key order explicitly so the test can assert column order
        np.savez(
            os.path.join(GOLDEN_DIR, name + ".npz"),
            __order__=np.array(list(flat.keys())),
            **flat,
        )
        print(f"wrote {name}.npz  ({len(flat)} columns)")
```

- [ ] **Step 3: Generate the goldens from the pre-refactor build**

Run: `cd python && python gen_transfer_golden.py`
Expected: prints one line per case (e.g. `wrote dcdm_dr.npz  (NN columns)`) and creates `python/transfer_golden/*.npz`. If any case errors on input keys, adjust that case's parameters until it computes (the exact physics values are not important — only that the species are active and the run reproduces).

- [ ] **Step 4: Write the characterization test**

Create `python/test_transfer_columns.py`:

```python
"""Characterization test for issue #309 (species-registered transfer sources).

Locks the Tk/vTk transfer-function output (column names, column order, and
values) against goldens captured from the pre-refactor build. The #309 refactor
must keep every assertion green: it changes only where source-slot integers are
stored, never the physics.
"""
import os
import numpy as np
import pytest

from gen_transfer_golden import CASES, transfer_for  # reuse the case list + runner

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "transfer_golden")


@pytest.mark.parametrize("name", sorted(CASES))
def test_transfer_columns_match_golden(name):
    golden = np.load(os.path.join(GOLDEN_DIR, name + ".npz"), allow_pickle=False)
    expected_order = list(golden["__order__"])

    flat = transfer_for(CASES[name])

    # 1. column identity + order unchanged
    assert list(flat.keys()) == expected_order, (
        f"{name}: transfer column set/order changed.\n"
        f"  expected: {expected_order}\n  got:      {list(flat.keys())}"
    )

    # 2. values byte-identical (pure plumbing refactor => exact equality)
    for key in expected_order:
        np.testing.assert_array_equal(
            flat[key], golden[key], err_msg=f"{name}: column {key!r} changed"
        )
```

- [ ] **Step 5: Run the test — expect green on the pre-refactor build**

Run: `cd python && python -m pytest test_transfer_columns.py -v`
Expected: all parametrized cases PASS (the goldens were just generated from this same build).

- [ ] **Step 6: Commit**

```bash
git add python/gen_transfer_golden.py python/test_transfer_columns.py python/transfer_golden
git commit -m "Add transfer-column characterization test + goldens for #309"
```

---

## Task 2: Add the hook, the context, and the (no-op) registration loop

No behavior change: the new virtual defaults to no-op, the old per-species block stays, so the loop allocates nothing yet.

**Files:**
- Modify: `species/perturb_source_context.h` (append struct)
- Modify: `species/base_species.h` (new virtual)
- Modify: `species/composite_species.h`, `species/composite_species.cpp` (forwarding override)
- Modify: `source/perturbations_module.cpp` (`perturb_indices_of_perturbs`, scalar branch)

- [ ] **Step 1: Add `SourceRequestContext`**

Append to `species/perturb_source_context.h` (after the existing `PerturbScalarContext`):

```cpp
/**
 * Output requests a species needs in order to decide which transfer-source
 * slots to register. Plain const input, not per-k state (see issue #309).
 */
struct SourceRequestContext {
  bool wants_density;   // ppt->has_density_transfers
  bool wants_velocity;  // ppt->has_velocity_transfers
  int  gauge;           // possible_gauges; CDM registers theta only when != synchronous
};
```

- [ ] **Step 2: Add the virtual to `BaseSpecies`**

In `species/base_species.h`, next to `RegisterPerturbationIndices` (around line 257), add:

```cpp
/**
 * Register this species' transfer-source (index_tp) slots. Set once in
 * perturb_indices_of_perturbs() before the k-loop; mirrors
 * RegisterBackgroundIndices (k-independent, cached in a plain species member).
 * Takes no PerturbLayout — that absence marks it as non-per-k state.
 * Default: no-op (species that emit no transfer functions).
 */
virtual void RegisterTransferSourceIndices(int& /*index_tp*/,
                                           const SourceRequestContext& /*ctx*/) {}
```

Ensure `perturb_source_context.h` is included where `SourceRequestContext` is visible (it is already included for the existing context types; confirm with a build).

- [ ] **Step 3: Add the generic composite forwarding**

In `species/composite_species.h`, declare (near the other `Register*Indices` overrides):

```cpp
void RegisterTransferSourceIndices(int& index_tp,
                                   const SourceRequestContext& ctx) override;
```

In `species/composite_species.cpp`, add (mirroring `RegisterBackgroundIndices`):

```cpp
void CompositeSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                     const SourceRequestContext& ctx) {
  for (auto& child : children_)
    child->RegisterTransferSourceIndices(index_tp, ctx);
}
```

- [ ] **Step 4: Add the registration loop after the module-owned head**

In `source/perturbations_module.cpp`, inside `perturb_indices_of_perturbs()`, in the `if (_scalars_)` branch, immediately **after** the existing `class_define_index(index_tp_theta_tot_, …);` line (the last module-owned `theta` head type, currently near line 1091) and **before** `class_define_index(index_tp_theta_g_, …)`, insert:

```cpp
      // Per-species transfer sources register themselves (issue #309).
      // Migrated species allocate here; not-yet-migrated species are still
      // allocated by the legacy block below until their task lands.
      const SourceRequestContext src_ctx{ppt->has_density_transfers,
                                         ppt->has_velocity_transfers,
                                         ppt->gauge};
      for (auto& [name, sp] : all_species_)
        sp->RegisterTransferSourceIndices(index_type, src_ctx);
```

(Placement note: putting the loop here means a migrated species' δ and θ land adjacently at this point; unmigrated species keep their legacy slots. Both are internally consistent and output-invariant. The loop moves to its final position in Task 14.)

- [ ] **Step 5: GREEN GATE**

Run: `pip install . --no-build-isolation`
Run: `cd python && python -m pytest test_transfer_columns.py -v`
Expected: builds clean; all cases PASS (loop is a no-op, nothing migrated yet).

- [ ] **Step 6: Commit**

```bash
git add species/perturb_source_context.h species/base_species.h \
        species/composite_species.h species/composite_species.cpp \
        source/perturbations_module.cpp
git commit -m "#309: add RegisterTransferSourceIndices hook + no-op registration loop"
```

---

## Task 3: Migrate Fluid (template — both δ and θ unconditional)

**Files:**
- Modify: `species/fluid.h` (add members)
- Modify: `species/fluid.cpp` (`FillSources` repoint, `WriteOutputColumns` repoint, new hook)
- Modify: `source/perturbations_module.h`, `source/perturbations_module.cpp` (delete `fld` symbols)

- [ ] **Step 1: Add members to `species/fluid.h`**

In the `private:` section (after `int index_bg_rho_fld_ = -1;`, ~line 212):

```cpp
  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
```

- [ ] **Step 2: Add the hook to `species/fluid.cpp`**

Add the method (declare it in `fluid.h` next to the other perturbation overrides as
`void RegisterTransferSourceIndices(int&, const SourceRequestContext&) override;`):

```cpp
void FluidSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                 const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3: Repoint `FillSources` and `WriteOutputColumns`**

In `species/fluid.cpp`:
- line ~174: `p_mod->index_tp_delta_fld_` → `index_tp_delta_`
- line ~187: `p_mod->index_tp_theta_fld_` → `index_tp_theta_`
- line ~225: `w.Add("d_fld", mod.index_tp_delta_fld_, _TRUE_);` → `w.Add("d_fld", index_tp_delta_, index_tp_delta_ >= 0);`
- line ~227: `w.Add("t_fld", mod.index_tp_theta_fld_, _TRUE_);` → `w.Add("t_fld", index_tp_theta_, index_tp_theta_ >= 0);`

- [ ] **Step 4: Delete the `fld` symbols from the module**

- `source/perturbations_module.h`: remove `int index_tp_delta_fld_;`, `int index_tp_theta_fld_;`, `bool has_source_delta_fld_;`, `bool has_source_theta_fld_;`.
- `source/perturbations_module.cpp`: remove `if (all_species_.count("Fluid")) has_source_delta_fld_ = true;`, the matching `… has_source_theta_fld_ = true;`, and the two `class_define_index(index_tp_delta_fld_, …)` / `(index_tp_theta_fld_, …)` lines.

- [ ] **Step 5: Verify no stragglers**

Run: `grep -rn "index_tp_delta_fld_\|index_tp_theta_fld_\|has_source_delta_fld_\|has_source_theta_fld_" source/ species/`
Expected: **no output** (all references gone).

- [ ] **Step 6: GREEN GATE + commit**

Run: `pip install . --no-build-isolation`
Run: `cd python && python -m pytest test_transfer_columns.py -v`
Expected: all PASS (esp. the `fluid` case).

```bash
git add species/fluid.h species/fluid.cpp source/perturbations_module.h source/perturbations_module.cpp
git commit -m "#309: migrate Fluid to self-registered transfer sources"
```

---

## Task 4: Migrate UltraRelativistic (UR)

Same recipe as Task 3. Symbols: `index_tp_delta_ur_`, `index_tp_theta_ur_`, `has_source_delta_ur_`, `has_source_theta_ur_`.

**Files:** `species/ultra_relativistic.h`, `species/ultra_relativistic.cpp`, `source/perturbations_module.{h,cpp}`.

- [ ] **Step 1:** In `species/ultra_relativistic.h` `private:` (~line 143) add `int index_tp_delta_ = -1;` and `int index_tp_theta_ = -1;`, and declare `void RegisterTransferSourceIndices(int&, const SourceRequestContext&) override;`.

- [ ] **Step 2:** In `species/ultra_relativistic.cpp` add:

```cpp
void UltraRelativisticSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                             const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint in `species/ultra_relativistic.cpp`:
  - line ~276: `p_mod->index_tp_delta_ur_` → `index_tp_delta_`
  - line ~288: `p_mod->index_tp_theta_ur_` → `index_tp_theta_`
  - line ~340: `w.Add("d_ur", mod.index_tp_delta_ur_, _TRUE_);` → `w.Add("d_ur", index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~342: `w.Add("t_ur", mod.index_tp_theta_ur_, _TRUE_);` → `w.Add("t_ur", index_tp_theta_, index_tp_theta_ >= 0);`
  - line ~346: `w.Add("-T_ur/k2", mod.index_tp_delta_ur_, _TRUE_);` → `w.Add("-T_ur/k2", index_tp_delta_, index_tp_delta_ >= 0);`

- [ ] **Step 4:** Delete `ur` symbols from `source/perturbations_module.{h,cpp}` per the Module deletion rule.

- [ ] **Step 5:** `grep -rn "index_tp_delta_ur_\|index_tp_theta_ur_\|has_source_delta_ur_\|has_source_theta_ur_" source/ species/` → no output.

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate UltraRelativistic to self-registered transfer sources`).

---

## Task 5: Migrate ScalarField (scf)

Same recipe. Symbols: `index_tp_delta_scf_`, `index_tp_theta_scf_`, `has_source_delta_scf_`, `has_source_theta_scf_`.

**Files:** `species/scalar_field.h`, `species/scalar_field.cpp`, `source/perturbations_module.{h,cpp}`.

- [ ] **Step 1:** In `species/scalar_field.h` add the two members + the override declaration.

- [ ] **Step 2:** In `species/scalar_field.cpp` add:

```cpp
void ScalarFieldSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                       const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint in `species/scalar_field.cpp`:
  - line ~245: `p_mod->index_tp_delta_scf_` → `index_tp_delta_`
  - line ~257: `p_mod->index_tp_theta_scf_` → `index_tp_theta_`
  - line ~301: `w.Add("d_scf", mod.index_tp_delta_scf_, _TRUE_);` → `w.Add("d_scf", index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~303: `w.Add("t__scf", mod.index_tp_theta_scf_, _TRUE_);` → `w.Add("t__scf", index_tp_theta_, index_tp_theta_ >= 0);` (keep the existing `"t__scf"` spelling — changing it would move a column and break the golden)

- [ ] **Step 4:** Delete `scf` symbols from `source/perturbations_module.{h,cpp}`.

- [ ] **Step 5:** `grep -rn "index_tp_delta_scf_\|index_tp_theta_scf_\|has_source_delta_scf_\|has_source_theta_scf_" source/ species/` → no output.

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate ScalarField to self-registered transfer sources`).

---

## Task 6: Migrate CDM (θ is gauge-conditional)

Same recipe, with one twist: `theta_cdm` is registered only in Newtonian gauge.

**Files:** `species/cdm.h`, `species/cdm.cpp`, `source/perturbations_module.{h,cpp}`. Symbols: `index_tp_delta_cdm_`, `index_tp_theta_cdm_`, `has_source_delta_cdm_`, `has_source_theta_cdm_`.

- [ ] **Step 1:** In `species/cdm.h` add the two members + override declaration.

- [ ] **Step 2:** In `species/cdm.cpp` add (note the gauge gate, mirroring the legacy `has_source_theta_cdm_` condition at perturbations_module.cpp:992-993):

```cpp
void CDMSpecies::RegisterTransferSourceIndices(int& index_tp,
                                               const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_,
                     ctx.wants_velocity && ctx.gauge != synchronous, index_tp, 1);
}
```

(`synchronous` is the `possible_gauges` enumerator already visible in this translation unit via `cdm.cpp`'s existing `ppt->gauge != synchronous` use at line 182.)

- [ ] **Step 3:** Repoint in `species/cdm.cpp`:
  - line ~153: `p_mod->index_tp_delta_cdm_` → `index_tp_delta_`
  - line ~165: `p_mod->index_tp_theta_cdm_` → `index_tp_theta_`
  - line ~180: `w.Add("d_cdm", mod.index_tp_delta_cdm_, _TRUE_);` → `w.Add("d_cdm", index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~182: `w.Add("t_cdm", mod.index_tp_theta_cdm_, (ppt->gauge != synchronous));` → `w.Add("t_cdm", index_tp_theta_, index_tp_theta_ >= 0);` (the member is set iff Newtonian, so `>= 0` reproduces the old gauge active)
  - line ~186: `w.Add("-T_cdm/k2", mod.index_tp_delta_cdm_, _TRUE_);` → `w.Add("-T_cdm/k2", index_tp_delta_, index_tp_delta_ >= 0);`

- [ ] **Step 4:** Delete `cdm` (the `index_tp_delta_cdm_/theta_cdm_` + `has_source_delta_cdm_/theta_cdm_`) symbols from `source/perturbations_module.{h,cpp}`. **Do not** touch `index_tp_delta_cb_`/`index_tp_theta_cb_` — those are module-owned combined cdm+baryon types and stay.

- [ ] **Step 5:** `grep -rn "index_tp_delta_cdm_\|index_tp_theta_cdm_\|has_source_delta_cdm_\|has_source_theta_cdm_" source/ species/` → no output. (Note the `_cb_` variants must still be present — that is correct.)

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate CDM to self-registered transfer sources`).

---

## Task 7: Migrate Photons

Same recipe; both δ_g/θ_g gated on density/velocity. The legacy flags `has_source_delta_g_`/`has_source_theta_g_` were set unconditionally inside the density/velocity blocks (perturbations_module.cpp:955, 990), so density/velocity gating reproduces them.

**Files:** `species/photons.h`, `species/photons.cpp`, `source/perturbations_module.{h,cpp}`. Symbols: `index_tp_delta_g_`, `index_tp_theta_g_`, `has_source_delta_g_`, `has_source_theta_g_`.

- [ ] **Step 1:** In `species/photons.h` add the two members + override declaration.

- [ ] **Step 2:** In `species/photons.cpp` add:

```cpp
void PhotonsSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                   const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint in `species/photons.cpp`:
  - line ~508: `p_mod->index_tp_delta_g_` → `index_tp_delta_`
  - line ~517: `p_mod->index_tp_theta_g_` → `index_tp_theta_`
  - line ~576: `w.Add("d_g", mod.index_tp_delta_g_, true);` → `w.Add("d_g", index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~578: `w.Add("t_g", mod.index_tp_theta_g_, true);` → `w.Add("t_g", index_tp_theta_, index_tp_theta_ >= 0);`

- [ ] **Step 4:** Delete `g` symbols from `source/perturbations_module.{h,cpp}`.

- [ ] **Step 5:** `grep -rn "index_tp_delta_g_\|index_tp_theta_g_\|has_source_delta_g_\|has_source_theta_g_" source/ species/` → no output.

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate Photons to self-registered transfer sources`).

---

## Task 8: Migrate Baryons

Same recipe. Symbols: `index_tp_delta_b_`, `index_tp_theta_b_`, `has_source_delta_b_`, `has_source_theta_b_`.

**Files:** `species/baryons.h`, `species/baryons.cpp`, `source/perturbations_module.{h,cpp}`.

- [ ] **Step 1:** In `species/baryons.h` add the two members + override declaration.

- [ ] **Step 2:** In `species/baryons.cpp` add:

```cpp
void BaryonsSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                   const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint in `species/baryons.cpp`:
  - line ~159: `p_mod->index_tp_delta_b_` → `index_tp_delta_`
  - line ~170: `p_mod->index_tp_theta_b_` → `index_tp_theta_`
  - line ~184: `w.Add("d_b", mod.index_tp_delta_b_, true);` → `w.Add("d_b", index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~186: `w.Add("t_b", mod.index_tp_theta_b_, true);` → `w.Add("t_b", index_tp_theta_, index_tp_theta_ >= 0);`

- [ ] **Step 4:** Delete `b` symbols from `source/perturbations_module.{h,cpp}`. (Leave `delta_cb`/`theta_cb` — module-owned.)

- [ ] **Step 5:** `grep -rn "index_tp_delta_b_\|index_tp_theta_b_\|has_source_delta_b_\|has_source_theta_b_" source/ species/` → no output.

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate Baryons to self-registered transfer sources`).

---

## Task 9: Migrate IDM_DR_IDR (composite owns two channels)

The composite writes both channels' sources/columns itself and holds the `has_idm_dr()`/`has_idr()` gating, so the channel members live **on the composite** and it **overrides** the register hook (rather than using the generic child-forwarding).

**Files:** `species/idm_dr_idr_species.h`, `species/idm_dr_idr_species.cpp`, `source/perturbations_module.{h,cpp}`. Module symbols: `index_tp_delta_idm_dr_`, `index_tp_theta_idm_dr_`, `index_tp_delta_idr_`, `index_tp_theta_idr_` (+ the four `has_source_*` flags).

- [ ] **Step 1:** In `species/idm_dr_idr_species.h` `private:` add:

```cpp
  int index_tp_delta_idm_dr_ = -1;
  int index_tp_theta_idm_dr_ = -1;
  int index_tp_delta_idr_    = -1;
  int index_tp_theta_idr_    = -1;
```

and declare `void RegisterTransferSourceIndices(int&, const SourceRequestContext&) override;`.

- [ ] **Step 2:** In `species/idm_dr_idr_species.cpp` add (gating mirrors the legacy flags: idm_dr channel exists iff `has_idm_dr()`, idr channel iff `has_idr()`):

```cpp
void IDM_DR_IDR_Species::RegisterTransferSourceIndices(int& index_tp,
                                                       const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_idm_dr_, ctx.wants_density  && has_idm_dr(), index_tp, 1);
  class_define_index(index_tp_delta_idr_,    ctx.wants_density  && has_idr(),    index_tp, 1);
  class_define_index(index_tp_theta_idm_dr_, ctx.wants_velocity && has_idm_dr(), index_tp, 1);
  class_define_index(index_tp_theta_idr_,    ctx.wants_velocity && has_idr(),    index_tp, 1);
}
```

(Order δ_idm_dr, δ_idr, θ_idm_dr, θ_idr matches the legacy `class_define_index` order so any same-pass adjacency is preserved; column order is set separately by `WriteOutputColumns` and is unaffected regardless.)

- [ ] **Step 3:** Repoint in `species/idm_dr_idr_species.cpp` (replace `p_mod->index_tp_<x>_` / `mod.index_tp_<x>_` with the composite member of the same name; the existing `has_idm_dr()`/`has_idr()` actives stay):
  - line ~88: `p_mod->index_tp_delta_idm_dr_` → `index_tp_delta_idm_dr_`
  - line ~94: `p_mod->index_tp_theta_idm_dr_` → `index_tp_theta_idm_dr_`
  - lines ~101 and ~105: `p_mod->index_tp_delta_idr_` → `index_tp_delta_idr_`
  - lines ~112 and ~115: `p_mod->index_tp_theta_idr_` → `index_tp_theta_idr_`
  - line ~127: `mod.index_tp_delta_idm_dr_` → `index_tp_delta_idm_dr_`
  - line ~128: `mod.index_tp_delta_idr_` → `index_tp_delta_idr_`
  - line ~131: `mod.index_tp_theta_idm_dr_` → `index_tp_theta_idm_dr_`
  - line ~132: `mod.index_tp_theta_idr_` → `index_tp_theta_idr_`
  - line ~137: `mod.index_tp_delta_idm_dr_` → `index_tp_delta_idm_dr_`

  (Leave the `has_idm_dr() ? _TRUE_ : _FALSE_` / `has_idr()` active arguments unchanged — they already encode channel existence.)

- [ ] **Step 4:** Delete the four `idm_dr`/`idr` `index_tp_*` + `has_source_*` symbols from `source/perturbations_module.{h,cpp}`.

- [ ] **Step 5:** `grep -rn "index_tp_delta_idm_dr_\|index_tp_theta_idm_dr_\|index_tp_delta_idr_\|index_tp_theta_idr_\|has_source_delta_idm_dr_\|has_source_theta_idm_dr_\|has_source_delta_idr_\|has_source_theta_idr_" source/` → no output. (Hits remain only in `species/idm_dr_idr_species.*`, which is correct — they are now the composite's own members.)

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate IDM_DR_IDR to self-registered transfer sources`).

---

## Task 10: Migrate IDM_DRMD_IDR_DRMD (composite owns two channels)

Same shape as Task 9. Module symbols: `index_tp_delta_idm_drmd_`, `index_tp_theta_idm_drmd_`, `index_tp_delta_idr_drmd_`, `index_tp_theta_idr_drmd_` (+ four `has_source_*`). Gating predicates: `has_idm_drmd()`, `has_idr_drmd()` (the accessors used at idm_drmd_idr_drmd_species.cpp:119-124).

**Files:** `species/idm_drmd_idr_drmd_species.h`, `species/idm_drmd_idr_drmd_species.cpp`, `source/perturbations_module.{h,cpp}`.

- [ ] **Step 1:** In `species/idm_drmd_idr_drmd_species.h` add the four members + override declaration:

```cpp
  int index_tp_delta_idm_drmd_ = -1;
  int index_tp_theta_idm_drmd_ = -1;
  int index_tp_delta_idr_drmd_ = -1;
  int index_tp_theta_idr_drmd_ = -1;
```

- [ ] **Step 2:** In `species/idm_drmd_idr_drmd_species.cpp` add:

```cpp
void IDM_DRMD_IDR_DRMD_Species::RegisterTransferSourceIndices(int& index_tp,
                                                              const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_idm_drmd_, ctx.wants_density  && has_idm_drmd(), index_tp, 1);
  class_define_index(index_tp_delta_idr_drmd_, ctx.wants_density  && has_idr_drmd(), index_tp, 1);
  class_define_index(index_tp_theta_idm_drmd_, ctx.wants_velocity && has_idm_drmd(), index_tp, 1);
  class_define_index(index_tp_theta_idr_drmd_, ctx.wants_velocity && has_idr_drmd(), index_tp, 1);
}
```

- [ ] **Step 3:** Repoint in `species/idm_drmd_idr_drmd_species.cpp`:
  - line ~87: `p_mod->index_tp_delta_idm_drmd_` → `index_tp_delta_idm_drmd_`
  - line ~93: `p_mod->index_tp_theta_idm_drmd_` → `index_tp_theta_idm_drmd_`
  - line ~99: `p_mod->index_tp_delta_idr_drmd_` → `index_tp_delta_idr_drmd_`
  - line ~107: `p_mod->index_tp_theta_idr_drmd_` → `index_tp_theta_idr_drmd_`
  - line ~119: `mod.index_tp_delta_idm_drmd_` → `index_tp_delta_idm_drmd_`
  - line ~120: `mod.index_tp_delta_idr_drmd_` → `index_tp_delta_idr_drmd_`
  - line ~123: `mod.index_tp_theta_idm_drmd_` → `index_tp_theta_idm_drmd_`
  - line ~124: `mod.index_tp_theta_idr_drmd_` → `index_tp_theta_idr_drmd_`

- [ ] **Step 4:** Delete the four `idm_drmd`/`idr_drmd` `index_tp_*` + `has_source_*` symbols from `source/perturbations_module.{h,cpp}`.

- [ ] **Step 5:** `grep -rn "index_tp_.*_idm_drmd_\|index_tp_.*_idr_drmd_\|has_source_.*_idm_drmd_\|has_source_.*_idr_drmd_" source/` → no output (hits remain only in the species `.cpp`/`.h`).

- [ ] **Step 6:** GREEN GATE + commit (`#309: migrate IDM_DRMD_IDR_DRMD to self-registered transfer sources`).

---

## Task 11: Migrate the NCDM family (dissolve the ncdm1 shared block + its SetSourceSlot loop)

Each NCDM instance registers its own δ/θ; the member lives on `NCDMBaseSpecies` so the decaying child (`DNCDMSpecies`) also registers a slot it never writes (preserving the legacy "DNCDM consumes a slot silently" behavior). The `index_tp_delta_ncdm1_ + slot` base-block, `source_slot_`/`SetSourceSlot`, and the NCDM `SetSourceSlot` loop all go away.

**Files:** `species/ncdm_base_species.h`, `species/ncdm_base_species.cpp`, `species/ncdm_species.cpp`, `source/perturbations_module.{h,cpp}`. (Confirm whether `SetSourceSlot`/`source_slot_` are declared on `NCDMSpecies` (ncdm_species.h:118,150) or the base — remove wherever declared.)

- [ ] **Step 1:** In `species/ncdm_base_species.h` `private:`/`protected:` add:

```cpp
  int index_tp_delta_ = -1;  // #309 transfer-source slot (this instance)
  int index_tp_theta_ = -1;
```

and declare `void RegisterTransferSourceIndices(int&, const SourceRequestContext&) override;`.

- [ ] **Step 2:** In `species/ncdm_base_species.cpp` add:

```cpp
void NCDMBaseSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                    const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint NCDM source-writing in `species/ncdm_species.cpp` (drop the `+ n`/`source_slot_` offset — each instance owns its slot):
  - line ~269: delete `const int n = source_slot_;` (and any later `n` use) — read the instance member instead.
  - line ~279: `p_mod->index_tp_delta_ncdm1_ + n` → `index_tp_delta_`
  - line ~291: `p_mod->index_tp_theta_ncdm1_ + n` → `index_tp_theta_`
  - line ~437: `const int n = source_slot_;` (second occurrence in `WriteOutputColumns`) → remove.
  - line ~443: `w.Add("d_" + nm, mod.index_tp_delta_ncdm1_ + n, _TRUE_);` → `w.Add("d_" + nm, index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~445: `w.Add("t_" + nm, mod.index_tp_theta_ncdm1_ + n, _TRUE_);` → `w.Add("t_" + nm, index_tp_theta_, index_tp_theta_ >= 0);`
  - line ~450: `w.Add("-T_ncdm/k2", mod.index_tp_delta_ncdm1_, _TRUE_);` → `w.Add("-T_ncdm/k2", index_tp_delta_, index_tp_delta_ >= 0);` (read the file around 437-452 to confirm this is inside the per-instance `WriteOutputColumns`; if it was relying on the base index for instance 0, it now correctly uses this instance's member)

- [ ] **Step 4:** Remove `SetSourceSlot(int)` and `int source_slot_` and the `source_slot()` accessor from `species/ncdm_species.h` (lines ~116-122, ~150) — or from `NCDMBaseSpecies` if that is where they live. `grep -rn "SetSourceSlot\|source_slot_\|source_slot()" species/ncdm* ` to locate all.

- [ ] **Step 5:** In `source/perturbations_module.cpp` delete the **NCDM `SetSourceSlot` loop** (the block at ~639-649 iterating `all_species_` casting to `NCDMSpecies*`/`DNCDM_DR_Species*`). Delete the `index_tp_delta_ncdm1_`/`index_tp_theta_ncdm1_` `class_define_index` lines (the ones sized by `NcdmFamily(all_species_).size()`), the `has_source_delta_ncdm_`/`has_source_theta_ncdm_` gating assignments, and in `source/perturbations_module.h` the `index_tp_delta_ncdm1_`, `index_tp_theta_ncdm1_`, `has_source_delta_ncdm_`, `has_source_theta_ncdm_` declarations.

- [ ] **Step 6:** Verify:

Run: `grep -rn "index_tp_delta_ncdm1_\|index_tp_theta_ncdm1_\|has_source_delta_ncdm_\|has_source_theta_ncdm_\|SetSourceSlot\|source_slot_" source/ species/`
Expected: **no output** for the ncdm1/has_source/SetSourceSlot symbols. (The DR `SetSourceSlot` loop is removed in Task 13; if `SetSourceSlot` still appears for DR here, that is expected until Task 13 — narrow the grep to `ncdm` if needed.)

- [ ] **Step 7:** GREEN GATE (esp. `cdm_ur_ncdm` and `dncdm_dr` cases) + commit (`#309: migrate NCDM family to self-registered transfer sources; drop ncdm SetSourceSlot loop`).

---

## Task 12: Migrate the DCDM channel (DCDM_DR composite + DCDM child column)

The DCDM child owns the member (it writes its own `d_dcdm`/`t_dcdm` column); the `DCDM_DR` composite reads it via `dcdm().index_tp_delta_` for the source fill. `DCDM_DR` keeps the generic `CompositeSpecies` forwarding (added in Task 2), which now reaches `DCDMSpecies::RegisterTransferSourceIndices`.

**Files:** `species/dcdm.h`, `species/dcdm.cpp`, `species/dcdm_dr_species.cpp`, `source/perturbations_module.{h,cpp}`. Module symbols: `index_tp_delta_dcdm_`, `index_tp_theta_dcdm_`, `has_source_delta_dcdm_`, `has_source_theta_dcdm_`.

- [ ] **Step 1:** In `species/dcdm.h` add `int index_tp_delta_ = -1;`, `int index_tp_theta_ = -1;` (+ override declaration). Add public const accessors so the composite can read them (mirroring the `dcdm()` accessor pattern), or rely on `DCDM_DR_Species` friendship/existing `dcdm()` returning a non-const ref — simplest is two public getters:

```cpp
  int transfer_delta_index() const { return index_tp_delta_; }
  int transfer_theta_index() const { return index_tp_theta_; }
```

- [ ] **Step 2:** In `species/dcdm.cpp` add:

```cpp
void DCDMSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint the DCDM **column** in `species/dcdm.cpp`:
  - line ~163: `w.Add("d_dcdm", mod.index_tp_delta_dcdm_, _TRUE_);` → `w.Add("d_dcdm", index_tp_delta_, index_tp_delta_ >= 0);`
  - line ~165: `w.Add("t_dcdm", mod.index_tp_theta_dcdm_, _TRUE_);` → `w.Add("t_dcdm", index_tp_theta_, index_tp_theta_ >= 0);`

- [ ] **Step 4:** Repoint the DCDM **source fill** in `species/dcdm_dr_species.cpp`:
  - line ~207: `p_mod->index_tp_delta_dcdm_` → `dcdm().transfer_delta_index()`
  - line ~219: `p_mod->index_tp_theta_dcdm_` → `dcdm().transfer_theta_index()`

  (Leave lines ~236/250/266/268 — the DR channel — untouched; they are Task 13.)

- [ ] **Step 5:** Delete `dcdm` (`index_tp_delta_dcdm_/theta_dcdm_` + `has_source_*`) symbols from `source/perturbations_module.{h,cpp}`.

- [ ] **Step 6:** `grep -rn "index_tp_delta_dcdm_\|index_tp_theta_dcdm_\|has_source_delta_dcdm_\|has_source_theta_dcdm_" source/ species/` → no output.

- [ ] **Step 7:** GREEN GATE (`dcdm_dr` case) + commit (`#309: migrate DCDM channel to self-registered transfer sources`).

---

## Task 13: Migrate the DR channel (DarkRadiationSpecies + both DR-emitting composites; drop the DR SetSourceSlot loop)

`DarkRadiationSpecies` owns the member; `DCDM_DR` and `DNCDM_DR` read it via `dr().index_tp_delta_` (they already call `dr_sp_->source_slot()` today). The DR `SetSourceSlot` loop and `source_slot_`/`SetSourceSlot` on `DarkRadiationSpecies` are removed.

**Files:** `species/dark_radiation_species.h`, `species/dark_radiation_species.cpp`, `species/dcdm_dr_species.cpp`, `species/dncdm_dr_species.cpp`, `source/perturbations_module.{h,cpp}`. Module symbols: `index_tp_delta_dr_`, `index_tp_theta_dr_`, `has_source_delta_dr_`, `has_source_theta_dr_`.

- [ ] **Step 1:** In `species/dark_radiation_species.h` add `int index_tp_delta_ = -1;`, `int index_tp_theta_ = -1;`, public getters `transfer_delta_index()`/`transfer_theta_index()`, the override declaration, and remove `SetSourceSlot(int)` / `int source_slot_` / `source_slot()` (lines ~133-146).

- [ ] **Step 2:** In `species/dark_radiation_species.cpp` add:

```cpp
void DarkRadiationSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                         const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density,  index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}
```

- [ ] **Step 3:** Repoint in `species/dcdm_dr_species.cpp`:
  - line ~236: `p_mod->index_tp_delta_dr_ + dr_sp_->source_slot()` → `dr_sp_->transfer_delta_index()`
  - line ~250: `p_mod->index_tp_theta_dr_ + dr_sp_->source_slot()` → `dr_sp_->transfer_theta_index()`
  - line ~264: `const int slot = dr_sp_->source_slot();` → remove (and its uses below).
  - line ~266: `w.Add("d_" + dr_sp_->name(), mod.index_tp_delta_dr_ + slot, mod.has_source_delta_dr_);` → `w.Add("d_" + dr_sp_->name(), dr_sp_->transfer_delta_index(), dr_sp_->transfer_delta_index() >= 0);`
  - line ~268: `w.Add("t_" + dr_sp_->name(), mod.index_tp_theta_dr_ + slot, mod.has_source_theta_dr_);` → `w.Add("t_" + dr_sp_->name(), dr_sp_->transfer_theta_index(), dr_sp_->transfer_theta_index() >= 0);`

- [ ] **Step 4:** Repoint identically in `species/dncdm_dr_species.cpp`:
  - line ~347: `p_mod->index_tp_delta_dr_ + dr_sp_->source_slot()` → `dr_sp_->transfer_delta_index()`
  - line ~361: `p_mod->index_tp_theta_dr_ + dr_sp_->source_slot()` → `dr_sp_->transfer_theta_index()`
  - line ~377: `w.Add("d_" + dr_sp_->name(), mod.index_tp_delta_dr_ + slot, mod.has_source_delta_dr_);` → `w.Add("d_" + dr_sp_->name(), dr_sp_->transfer_delta_index(), dr_sp_->transfer_delta_index() >= 0);`
  - line ~379: `w.Add("t_" + dr_sp_->name(), mod.index_tp_theta_dr_ + slot, mod.has_source_theta_dr_);` → `w.Add("t_" + dr_sp_->name(), dr_sp_->transfer_theta_index(), dr_sp_->transfer_theta_index() >= 0);`
  - remove the `const int slot = dr_sp_->source_slot();` line in this file too if present.

- [ ] **Step 5:** In `source/perturbations_module.cpp` delete the **DR `SetSourceSlot` loop** (the block at ~654-664 casting to `DCDM_DR_Species*`/`DNCDM_DR_Species*` and calling `d->dr().SetSourceSlot`). Delete the `class_define_index(index_tp_delta_dr_, …, n_dr_species)` / `index_tp_theta_dr_` lines, the `has_source_delta_dr_`/`has_source_theta_dr_` gating, and the `h` declarations. If `int n_dr_species = DrSpeciesCount(all_species_);` is now unused in this function, remove it too.

- [ ] **Step 6:** Verify:

Run: `grep -rn "index_tp_delta_dr_\|index_tp_theta_dr_\|has_source_delta_dr_\|has_source_theta_dr_\|SetSourceSlot\|source_slot_\|source_slot()" source/ species/`
Expected: **no output** anywhere (both `SetSourceSlot` loops and all per-species `index_tp`/`has_source` source symbols are now gone).

- [ ] **Step 7:** GREEN GATE (`dcdm_dr`, `dncdm_dr` cases) + commit (`#309: migrate DR channel to self-registered transfer sources; drop dr SetSourceSlot loop`).

---

## Task 14: Final module cleanup + move the loop to its resting place

All per-species sources now self-register. Tidy the module: relocate the registration loop to its clean position and remove anything left dead.

**Files:** `source/perturbations_module.cpp`, `source/perturbations_module.h`.

- [ ] **Step 1:** Confirm the legacy per-species block is empty

Run: `grep -n "class_define_index(index_tp_" source/perturbations_module.cpp`
Expected: only **module-owned** types remain (`t0,t1,t2,p, delta_m,delta_cb,delta_tot, theta_m,theta_cb,theta_tot, phi,phi_prime,phi_plus_psi,psi,h,h_prime,eta,eta_prime,H_T_Nb_prime,k2gamma_Nb`, plus the perturbed-recombination types). No per-species δ/θ remain.

- [ ] **Step 2:** Position the loop cleanly

Move the `for (auto& [name, sp] : all_species_) sp->RegisterTransferSourceIndices(index_type, src_ctx);` loop (with its `src_ctx` construction) so it sits **after** the full module-owned head (all the `delta_*`/`theta_*` module-owned `class_define_index` lines) and **before** the metric/potential tail (`class_define_index(index_tp_phi_, …)` onward). The final layout: module-owned scalars head → species loop → metric tail → `tp_size_[index_md] = index_type;`.

- [ ] **Step 3:** Remove now-dead helpers if unused

Run: `grep -n "NcdmFamily\|DrSpeciesCount\|n_dr_species\|HasNcdm" source/perturbations_module.cpp`
For each hit inside `perturb_indices_of_perturbs()` that is now unused, remove it. **Keep** any uses elsewhere in the file (e.g. tensor/IC logic, `evolve_*` flags) — these belong to #308, not #309. If `DrSpeciesCount`/`NcdmFamily` become entirely unreferenced across the module, remove their free-function definitions too; otherwise leave them.

- [ ] **Step 4:** Sweep for any orphaned per-species source symbols

Run: `grep -rEn "index_tp_(delta|theta)_(g|b|cdm|dcdm|fld|scf|dr|ur|idr|idm_dr|idr_drmd|idm_drmd|ncdm1)_|has_source_(delta|theta)_(g|b|cdm|dcdm|fld|scf|dr|ur|idr|idm_dr|idr_drmd|idm_drmd|ncdm)_" source/`
Expected: **no output**.

- [ ] **Step 5:** GREEN GATE

Run: `pip install . --no-build-isolation`
Run: `cd python && python -m pytest test_transfer_columns.py -v`
Expected: builds clean (no unused-variable/function warnings on perturbations_module); all cases PASS.

- [ ] **Step 6:** Commit (`#309: relocate registration loop, drop dead source-sizing helpers`).

---

## Task 15: Full regression suite

**Files:** none (verification only).

- [ ] **Step 1:** Ensure `classyref` is importable (needed for `COMPARE_OUTPUT_REF`). If missing, install per the project's reference-build instructions; if it cannot be built, note it and rely on Steps 3–4.

- [ ] **Step 2:** Full scenario suite with reference comparison

Run: `cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest test_class.py -m test_scenario -q`
Expected: all scenarios pass (no `Reference comparison failed` assertions). Note: this compares `Cl`/`Pk` against `classyref`; the spec's expectation is no new failures vs. the master baseline.

- [ ] **Step 3:** Named multi-species guards

Run: `cd python && python -m pytest test_class.py -k "dcdm_dr or dncdm_dr or tensor_massive_ncdm or idm_dr_idr or idr_without_idm_dr or drmd_without_idr_drmd or scalar_field" -v`
Expected: all PASS.

- [ ] **Step 4:** Characterization test (final confirmation)

Run: `cd python && python -m pytest test_transfer_columns.py -v`
Expected: all PASS — byte-identical transfer output end to end.

- [ ] **Step 5:** Commit any notes; the branch is ready for PR.

```bash
git commit --allow-empty -m "#309: full regression green (scenarios + guards + characterization)"
```

---

## Self-Review (completed during planning)

- **Spec coverage:** the new hook + `SourceRequestContext` (Task 2); module-owned head/tail kept, per-species block + both `SetSourceSlot` loops deleted (Tasks 3–14); all 12 source-writing units migrated across the three structural shapes — leaves (3–8), interacting composites (9–10), multi-instance NCDM/DR (11–13); behavior-preservation locked by the characterization test (Task 1) and the full suite (Task 15); `background_indices()` cleanup and residual `HasNcdm`/`NcdmFamily` casts left to #308/sibling follow-up per spec.
- **Placeholder scan:** no TBD/"handle edge cases"; every code step shows code; repoints are exact symbol+location edits.
- **Type consistency:** the hook signature `RegisterTransferSourceIndices(int&, const SourceRequestContext&)` and the member names `index_tp_delta_`/`index_tp_theta_` (channel-suffixed on the IDM composites) are used identically across all tasks; the DCDM/DR child accessors `transfer_delta_index()`/`transfer_theta_index()` are defined where introduced (Tasks 12–13) and used consistently.
- **Known residual to confirm at execution:** Task 11 Step 3 line ~450 (`-T_ncdm/k2` using the base index) — read `ncdm_species.cpp:437-452` to confirm it sits in the per-instance `WriteOutputColumns`; the repoint to the instance member is correct in that case.
