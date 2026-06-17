# Issue #319 — clean-loop `background_indices()` + post-accumulation hooks — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the stringly-typed gates + `dynamic_cast` loops in `BackgroundModule::background_indices()` with a flat per-species registration loop, and dissolve the module-side residue (dead indices, scf duplicate indices, the scf/DRMD post-H exception blocks, the DRMD interaction physics) into proper per-species background lifecycle hooks.

**Architecture:** Every species already implements `RegisterBackgroundIndices(int&)`; the registration call site becomes `head → flat loop → tail`. Three new `BaseSpecies` lifecycle hooks — `PPrime()` (post-Friedmann pressure-derivative accumulator, replacing `DpDloga()`), `FinalizeBackground()` (post-Friedmann per-step owned-slot writes), `ProcessBackgroundTable()` (post-table-fill analysis) — absorb the work the module's hardcoded post-H blocks did. The DRMD interaction physics (`background_idm_drmd` + its IC constants + the decoupling scan) relocates into `IDM_DRMD_IDR_DRMD_Species`, mirroring the existing `FluidSpecies::ComputeWFld` precedent, with a thin module wrapper kept for cross-module callers.

**Tech Stack:** C++17 (`source/background_module.*`, species plugins in `species/`); Python `classy` Cython wrapper; `pytest`. Build: `pip install . --no-build-isolation` from repo root. Reference: the regenerated `classyref` + `COMPARE_OUTPUT_REF`/`TEST_LEVEL` knobs in `python/test_class.py`.

**Design doc:** `docs/superpowers/specs/2026-06-17-issue-319-background-indices-clean-loop-design.md`

---

## File Structure

| File | Responsibility | Tasks |
|------|----------------|-------|
| `python/gen_background_golden.py` (new) | One-shot golden generator: compute each cosmology, dump `get_background()` columns to `.npz` | 1 |
| `python/background_golden/*.npz` (new) | Captured pre-refactor background-table oracle | 1 |
| `python/test_background_columns.py` (new) | Assert candidate `get_background()` matches golden — exact except `(.)p_tot_prime` (tolerance) | 1 |
| `species/scalar_field.{h,cpp}` | scf owns its bg indices; `PPrime`/`FinalizeBackground`; drop `DpDloga`/`ComputePPrimeAndWrite` | 2, 4 |
| `source/background_module.{h,cpp}` | `background_indices()` head+loop+tail; delete dead/duplicate members; post-H accumulate loop; `ProcessBackgroundTable` loop; DRMD wrapper + IC delegation; banner accessors | 3, 4, 5, 6 |
| `species/base_species.h` | Remove `DpDloga()`; add `PPrime()`, `FinalizeBackground()`, `ProcessBackgroundTable()` | 4, 6 |
| All other `species/*.{h,cpp}` | `DpDloga` → `PPrime` (only where non-zero) | 4 |
| `species/idm_drmd_idr_drmd_species.{h,cpp}` | Own `index_bg_G_over_aH_drmd_`, `ComputeIdmDrmd`, IC constants, `FinalizeBackground`, `ProcessBackgroundTable` | 5, 6 |

**Build/test cadence (used in every task):**
- Build: `pip install . --no-build-isolation` (from `/Users/au192734/Projects/class_claude`) — expect success.
- Golden: `cd python && python -m pytest test_background_columns.py -v`.
- Scenarios (heavier; Tasks 4 & 7): `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -m test_scenario test_class.py -q`.

---

## Task 1: Characterization golden test (the safety net)

Capture the current-master background table for a species-spanning set of cosmologies, so every later task is guarded. Mirrors the #309 `test_transfer_columns.py` pattern.

**Files:**
- Create: `python/gen_background_golden.py`
- Create: `python/test_background_columns.py`
- Create: `python/background_golden/*.npz` (generated)

- [ ] **Step 1: Write the golden generator**

Create `python/gen_background_golden.py`:

```python
"""One-shot generator for background-table golden columns (issue #319).

Run ONCE on the pre-refactor build:  python gen_background_golden.py
Do NOT re-run after the refactor starts -- the goldens are the oracle.
"""
import os
import numpy as np
from classy import Class

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "background_golden")

COMMON = {"output": "mPk", "P_k_max_1/Mpc": 1.0}

CASES = {
    "lcdm": {},
    "curvature_open": {"Omega_k": 0.01},
    "curvature_closed": {"Omega_k": -0.01},
    "cdm_ur_ncdm": {
        "N_ur": 1.0, "N_ncdm": 2, "m_ncdm": "0.04,0.06", "T_ncdm": "0.71611,0.71611",
    },
    "fluid": {
        "Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.0, "N_ur": 3.044,
    },
    "scf": {
        "Omega_Lambda": 0.0, "Omega_fld": 0.0,
        "Omega_scf": -1, "attractor_ic_scf": "yes",
        "scf_parameters": "10.0, 0.0, 0.0, 0.0", "N_ur": 3.044,
    },
    "dcdm_dr": {
        "Omega_ini_dcdm": 0.01, "Gamma_dcdm": 100.0, "N_ur": 3.044,
    },
    "dncdm_dr": {
        "N_ur": 2.0, "YHe": 0.25,
        "dncdm1.type": "ncdm_decay_dr", "dncdm1.m": 1.0, "dncdm1.T": 0.71611,
        "dncdm1.Gamma": 10.0, "dncdm1.Omega_ini": 0.001,
    },
    "idm_dr_idr": {
        "Omega_idm_dr": 0.12, "a_idm_dr": 1000.0, "nindex_idm_dr": 4,
        "xi_idr": 0.3, "N_ur": 2.0,
    },
    "idm_drmd_idr_drmd": {
        "z_stop": 1.0e4, "G_over_aH_drmd_ini": 1.0,
        "f_idm_drmd": 0.1, "delta_Neff_drmd": 0.5, "N_ur": 2.0,
    },
}


def background_for(extra):
    cosmo = Class()
    cosmo.set({**COMMON, **extra})
    try:
        cosmo.compute()
        bg = cosmo.get_background()
        return {k: np.asarray(v, dtype=np.float64) for k, v in bg.items()}
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()


if __name__ == "__main__":
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    for name, extra in CASES.items():
        flat = background_for(extra)
        np.savez(
            os.path.join(GOLDEN_DIR, name + ".npz"),
            __order__=np.array(list(flat.keys())),
            **flat,
        )
        print(f"wrote {name}.npz  ({len(flat)} columns)")
```

- [ ] **Step 2: Write the characterization test**

Create `python/test_background_columns.py`:

```python
"""Characterization test for issue #319 (clean-loop background_indices()).

Locks the get_background() output (column names, order, values) against goldens
captured from the pre-refactor build. The refactor must keep every assertion
green: it relocates index registration and per-step writes, never the physics.

Exception: '(.)p_tot_prime' is checked at a tolerance, not exactly. The PPrime
refactor distributes the a*H factor across the per-species sum (Sum aH*d_i instead
of aH*Sum d_i) and folds the scalar field into the sum rather than adding it last;
both are ULP-level reorderings. Per project policy bit-identical is the wrong bar.
"""
import os
import numpy as np
import pytest

from gen_background_golden import CASES, background_for

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "background_golden")
TOL_COLUMNS = {"(.)p_tot_prime"}   # checked at rtol below, not exactly
PPRIME_RTOL = 1e-9


@pytest.mark.parametrize("name", sorted(CASES))
def test_background_columns_match_golden(name):
    golden = np.load(os.path.join(GOLDEN_DIR, name + ".npz"), allow_pickle=False)
    expected_order = list(golden["__order__"])

    bg = background_for(CASES[name])

    assert list(bg.keys()) == expected_order, (
        f"{name}: background column set/order changed.\n"
        f"  expected: {expected_order}\n  got:      {list(bg.keys())}"
    )

    for key in expected_order:
        if key in TOL_COLUMNS:
            np.testing.assert_allclose(
                bg[key], golden[key], rtol=PPRIME_RTOL, atol=0.0,
                err_msg=f"{name}: column {key!r} drifted beyond tolerance",
            )
        else:
            np.testing.assert_array_equal(
                bg[key], golden[key], err_msg=f"{name}: column {key!r} changed"
            )
```

- [ ] **Step 3: Build current master and generate the goldens**

Run: `pip install . --no-build-isolation`
Expected: success.
Run: `cd python && python gen_background_golden.py`
Expected: prints `wrote lcdm.npz (19 columns)` … one line per case, no exceptions.

- [ ] **Step 4: Run the characterization test — must PASS on master**

Run: `cd python && python -m pytest test_background_columns.py -v`
Expected: all 10 cases PASS (candidate == golden, both current master).

- [ ] **Step 5: Commit**

```bash
git add python/gen_background_golden.py python/test_background_columns.py python/background_golden
git commit -m "test(#319): characterization golden for get_background() columns"
```

---

## Task 2: Route scalar-field background-index reads to the species' own members

`ScalarFieldSpecies` already registers and owns `index_bg_phi_scf_`/`phi_prime`/`V`/`dV`/`ddV`/`p_prime_scf` (`scalar_field.h:48-55`). Stop reading the `BackgroundModule`'s duplicate copies. Pure read-redirect — byte-identical.

**Files:**
- Modify: `species/scalar_field.cpp` (reads via `bgm->index_bg_*_scf_` / `mod.GetBackgroundModule()->index_bg_*_scf_`)

- [ ] **Step 1: Redirect the module-qualified reads to own members**

In `species/scalar_field.cpp`, every read of the scf indices currently goes through the background module. All such reads are inside `ScalarFieldSpecies` methods, so the species' own members are in scope. Replace each `bgm->index_bg_X_scf_` and `mod.GetBackgroundModule()->index_bg_X_scf_` with the bare member `index_bg_X_scf_`, for `X` in {`phi`, `phi_prime`, `V`, `dV`, `ddV`}. Concretely the sites are:

- `:132` `pvecback[bgm_->index_bg_H_]` — **leave** (`index_bg_H_` is module-owned, not scf).
- `:194-195` `ppw->pvecback[bgm->index_bg_phi_prime_scf_]`, `…index_bg_ddV_scf_` → `index_bg_phi_prime_scf_`, `index_bg_ddV_scf_`.
- `:222-223` `pvecback[bgm->index_bg_phi_prime_scf_]`, `…index_bg_dV_scf_` → own members.
- `:289-290` `pvecback[bgm->index_bg_phi_prime_scf_]`, `…index_bg_phi_scf_` → own members.
- `:330-331` `pvecback[mod.GetBackgroundModule()->index_bg_phi_prime_scf_]`, `…index_bg_dV_scf_` → own members.

(The bare-member reads at `:107-114`, `:131,133`, `:155-159`, `:376-377`, `:403`, `:415-416` already use the species' members — leave them.)

- [ ] **Step 2: Build**

Run: `pip install . --no-build-isolation`
Expected: success, no new warnings.

- [ ] **Step 3: Golden test — exact PASS**

Run: `cd python && python -m pytest test_background_columns.py -v`
Expected: all PASS (byte-identical; the `scf` case in particular).

- [ ] **Step 4: Commit**

```bash
git add species/scalar_field.cpp
git commit -m "refactor(#319): scalar_field reads its own bg indices, not the module copies"
```

---

## Task 3: Flat registration loop + delete dead/duplicate module indices

Rewrite `background_indices()` as `head → flat loop → tail`, and delete the module members that are now dead (`number_ncdm1`, `pseudo_p_ncdm1`, `Gamma0_drmd` index) or duplicated (the 5 scf copies, now unused after Task 2). The live `index_bg_G_over_aH_drmd_` stays module-owned for now (moved to the species in Task 5), parked in the tail.

**Files:**
- Modify: `source/background_module.cpp:586-707` (the `background_indices()` body)
- Modify: `source/background_module.h:52-62` (member declarations)

- [ ] **Step 1: Rewrite the species-registration section of `background_indices()`**

In `source/background_module.cpp`, replace the whole block from the photons line (`:611`) through the UR block (`:682`) — i.e. everything between `bg_size_short_ = index_bg;` and the `// Module aggregate indices` comment — with:

```cpp
  // ── Every species registers its own background indices (flat loop) ──────────
  for (auto& [name, sp] : all_species_)
    sp->RegisterBackgroundIndices(index_bg);
```

This deletes: the photons/baryons explicit lines, the CDM/IDM_DRMD/Lambda/Fluid/UR `count()`-gates, the `index_bg_G_over_aH_drmd_`/`index_bg_Gamma0_drmd_` module-define block, the NCDM `dynamic_cast` loop (incl. the `index_bg_number_ncdm1_`/`pseudo_p_ncdm1_` bookkeeping), the DNCDM_DR `dynamic_cast` loop, the DCDM_DR `static_cast` block, and the scf module-copy block (`:660-670`).

- [ ] **Step 2: Re-add the still-module-owned indices to the tail**

The module aggregate block (currently `:684-699`) keeps `rho_tot`/`p_tot`/`p_tot_prime`/`Omega_r`. `index_bg_G_over_aH_drmd_` and `IDM_DR_IDR` were registered outside the flat loop; `IDM_DR_IDR` now registers *inside* the flat loop (delete its old `:697-699` block). Keep `index_bg_G_over_aH_drmd_` module-owned until Task 5 by adding it to the tail right after `Omega_r`:

```cpp
  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(index_bg_Omega_r_, _TRUE_, index_bg, 1);

  /* DRMD interaction conductance G/(aH) — module-owned until #319 Task 5 moves
     it onto IDM_DRMD_IDR_DRMD_Species. */
  class_define_index(index_bg_G_over_aH_drmd_,
                     all_species_.count("IDM_DRMD_IDR_DRMD"), index_bg, 1);
```

Delete the old `IDM_DR + IDR composite (optional)` block at `:697-699`.

- [ ] **Step 3: Delete the dead/duplicate member declarations**

In `source/background_module.h`, delete these lines:
- `:53` `int index_bg_Gamma0_drmd_;`
- `:55-59` the five `int index_bg_phi_scf_; … index_bg_ddV_scf_;`
- `:61-62` `int index_bg_number_ncdm1_;` and `int index_bg_pseudo_p_ncdm1_;`

Keep `:52` `index_bg_G_over_aH_drmd_` (used until Task 5).

- [ ] **Step 4: Build**

Run: `pip install . --no-build-isolation`
Expected: success. If the compiler flags a leftover reference to any deleted member, that reference is a bug to remove (there should be none — verified: those members had no readers).

- [ ] **Step 5: Golden test — exact PASS**

Run: `cd python && python -m pytest test_background_columns.py -v`
Expected: all PASS. Output is by-name, so the reordered internal indices and the removed dead slots are invisible to `get_background()`.

- [ ] **Step 6: Commit**

```bash
git add source/background_module.cpp source/background_module.h
git commit -m "refactor(#319): flat per-species background registration loop; drop dead module indices"
```

---

## Task 4: `PPrime()` accumulator + `FinalizeBackground()` hook; scf post-H rework

Replace `DpDloga()` with a post-Friedmann `PPrime()` accumulator, add the void `FinalizeBackground()` hook, and dissolve the scf post-H exception. This is the one task that changes `(.)p_tot_prime` at ULP level (golden test's tolerance branch).

**Files:**
- Modify: `species/base_species.h:189` (the virtual)
- Modify: `source/background_module.cpp:316-411` (`background_functions` accumulation)
- Modify: every `species/*` with a `DpDloga` override (see enumeration)
- Modify: `species/scalar_field.{h,cpp}` (`PPrime`, `FinalizeBackground`, delete `ComputePPrimeAndWrite`)

- [ ] **Step 1: Swap the base-class virtual**

In `species/base_species.h`, delete:
```cpp
  virtual double DpDloga(const double* pvecback) const = 0;
```
and add (near the other background hooks, after `ComputeBackground`):
```cpp
  /**
   * Post-Friedmann. This species' conformal-time pressure derivative p'
   * (a' / a = a*H). pvecback_B is the ODE state (IN); pvecback is fully
   * populated through H. Default 0 (matter / Lambda).
   */
  virtual double PPrime(double a, double H,
                        const double* pvecback_B, const double* pvecback) const { return 0.; }

  /**
   * Post-Friedmann. Write this species' H-dependent owned output slots.
   * H passed explicitly so the hook need not re-read it. Default no-op.
   */
  virtual void FinalizeBackground(double a, double H,
                                  const double* pvecback_B, double* pvecback) {}
```

- [ ] **Step 2: Delete the zero-returning `DpDloga` overrides**

These returned 0 (matter / Lambda); with the new default-0 `PPrime` they need no override. Delete the `DpDloga` override (declaration + any out-of-line definition) from:
`species/cdm.{h,cpp}`, `species/baryons.h`, `species/lambda.{h,cpp}`, `species/dcdm.{h,cpp}`, `species/idm_dr.h`, `species/idm_drmd.h`, `species/scalar_field.h` (scf's `DpDloga` returned 0 — its real p' moves to `PPrime` in Step 4).

- [ ] **Step 3: Convert the non-zero `DpDloga` overrides to `PPrime`**

For each file below, replace the `DpDloga(const double* pvecback)` override with a `PPrime` override returning `a*H*(old body)`. Exact replacements:

`species/ultra_relativistic.h` (decl) + `species/ultra_relativistic.cpp`:
```cpp
double UltraRelativisticSpecies::PPrime(double a, double H,
                                        const double* /*pvecback_B*/, const double* pvecback) const {
  return a * H * (-4. / 3. * pvecback[index_bg_rho_]);
}
```
`species/photons.h`:
```cpp
  double PPrime(double a, double H,
                const double* /*pvecback_B*/, const double* pvecback) const override {
    return a * H * (-4. / 3. * pvecback[index_bg_rho_]);
  }
```
Identical body for `species/dark_radiation_species.h`, `species/idr.h`, `species/idr_drmd.h` (all `-4/3 rho`).

`species/ncdm_species.h`:
```cpp
  double PPrime(double a, double H,
                const double* /*pvecback_B*/, const double* pvecback) const override {
    return a * H * (pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_]);
  }
```
Identical body for `species/dncdm_species.h`.

`species/fluid.h` (decl) + `species/fluid.cpp`:
```cpp
double FluidSpecies::PPrime(double a, double H,
                            const double* /*pvecback_B*/, const double* pvecback) const {
  const double w_fld          = pvecback[index_bg_w_fld_];
  const double dw_over_da_fld = pvecback[index_bg_dw_over_da_fld_];
  return a * H * (a * dw_over_da_fld - 3. * (1. + w_fld) * w_fld) * pvecback[index_bg_rho_fld_];
}
```
(Drops the `bgm_->index_bg_a_` read — `a` is now a parameter.)

`species/composite_species.{h,cpp}`:
```cpp
double CompositeSpecies::PPrime(double a, double H,
                                const double* pvecback_B, const double* pvecback) const {
  double pp = 0.;
  for (const auto& child : children_)
    pp += child->PPrime(a, H, pvecback_B, pvecback);
  return pp;
}
```

- [ ] **Step 4: scf — `PPrime` (folds in the old exception) + `FinalizeBackground` (writes its slot)**

In `species/scalar_field.h`, delete the `ComputePPrimeAndWrite` declaration (`:177`) and add:
```cpp
  double PPrime(double a, double H,
                const double* pvecback_B, const double* pvecback) const override;
  void FinalizeBackground(double a, double H,
                          const double* pvecback_B, double* pvecback) override;
```
In `species/scalar_field.cpp`, delete `ComputePPrimeAndWrite` (`:130-137`) and add:
```cpp
double ScalarFieldSpecies::PPrime(double a, double H,
                                  const double* pvecback_B, const double* /*pvecback*/) const {
  const double phi       = pvecback_B[index_bi_phi_scf_];
  const double phi_prime = pvecback_B[index_bi_phi_prime_scf_];
  return phi_prime * (-phi_prime * H / a - 2. / 3. * dV_scf(phi));
}

void ScalarFieldSpecies::FinalizeBackground(double a, double H,
                                            const double* pvecback_B, double* pvecback) {
  pvecback[index_bg_p_prime_scf_] = PPrime(a, H, pvecback_B, pvecback);
}
```
Also in `ComputeBackground` (`scalar_field.cpp:114`) leave `pvecback[index_bg_p_prime_scf_] = 0.;` as the pre-H placeholder (overwritten by `FinalizeBackground`).

- [ ] **Step 5: Rewrite the accumulation in `background_functions`**

In `source/background_module.cpp`:
- Delete the `dp_dloga` declaration (`:323`) and the `dp_dloga += sp.DpDloga(pvecback);` line in the `accumulate` lambda (`:341`).
- Replace the `p_tot_prime` write + scf exception block (`:393-400`) with:
```cpp
  /* Derivative of total pressure w.r.t. conformal time: accumulated post-H so
     each species applies a' / a = a*H itself (the scalar field's p' is not of
     the a*H*dp/dlna form). FinalizeBackground writes any H-dependent owned slot
     (e.g. scf p_prime_scf, DRMD G_over_aH). */
  double p_tot_prime = 0.;
  for (const auto& [name, sp] : all_species_) {
    p_tot_prime += sp->PPrime(a, H, pvecback_B, pvecback);
    sp->FinalizeBackground(a, H, pvecback_B, pvecback);
  }
  pvecback[index_bg_p_tot_prime_] = p_tot_prime;
```
Here `H = pvecback[index_bg_H_]` (computed at `:382`) — bind `const double H = pvecback[index_bg_H_];` just after the Friedmann line, and reuse it for the `H_prime` line too if convenient. `pvecback_B` is the function's first argument.

Leave the DRMD `if (all_species_.count("IDM_DRMD_IDR_DRMD"))` block (`:402-411`) **in place** for now (Task 5 removes it).

- [ ] **Step 6: Build**

Run: `pip install . --no-build-isolation`
Expected: success. Any remaining `DpDloga` reference is a compile error — grep `grep -rn DpDloga species source` must return nothing.

- [ ] **Step 7: Golden test — `(.)p_tot_prime` within tolerance, all else exact**

Run: `cd python && python -m pytest test_background_columns.py -v`
Expected: all PASS. `(.)p_tot_prime` passes via the `rtol=1e-9` branch; the `scf` case's `(.)p_prime_scf` is exact (same formula).

- [ ] **Step 8: Scenario suite — TT/P(k) unchanged at the 0.1% bar**

Run: `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -m test_scenario test_class.py -q`
Expected: all PASS (no newtonian stale-ref fails against the regenerated `classyref`; `p_tot_prime` ULP shifts stay far under the Cl/Pk tolerance).

- [ ] **Step 9: Commit**

```bash
git add species source/background_module.cpp source/background_module.h
git commit -m "refactor(#319): PPrime() post-H accumulator replaces DpDloga(); add FinalizeBackground hook"
```

---

## Task 5: Relocate the DRMD interaction physics + `G_over_aH` index into the species

Move `background_idm_drmd` (body), the IC constants (`Gamma0_drmd`, `f_idr_drmd`; `z_stop` is already a species input), and the `index_bg_G_over_aH_drmd_` slot onto `IDM_DRMD_IDR_DRMD_Species`. Keep a thin `BackgroundModule::background_idm_drmd` wrapper for the two `perturbations_module` callers. Mirrors `FluidSpecies::ComputeWFld`.

**Files:**
- Modify: `species/idm_drmd_idr_drmd_species.{h,cpp}`
- Modify: `source/background_module.{h,cpp}`

- [ ] **Step 1: Add the species' owned state, index, and accessors**

In `species/idm_drmd_idr_drmd_species.h`, in the `private:` block add:
```cpp
  int index_bg_G_over_aH_drmd_ = -1;   // bg-table slot (this species owns it)
  double Gamma0_drmd_ic_       = 0.;   // interaction rate today, computed at IC time
  double f_idr_drmd_           = 0.;   // idr/idr_tot fraction at IC (verbose diagnostic)
```
and in the `public:` block add:
```cpp
  int bg_G_over_aH_index() const { return index_bg_G_over_aH_drmd_; }
  double f_idr_drmd() const { return f_idr_drmd_; }

  void RegisterBackgroundIndices(int& index_bg) override;
  void ComputeIdmDrmd(double a, double rho_idm_over_rho_idr,
                      double* Rint, double* csp2, double* Gint) const;
  void InitializeDrmdBackground(double rho_tot, double H, double a, double a_today);
  void FinalizeBackground(double a, double H,
                          const double* pvecback_B, double* pvecback) override;
```
(`z_stop_` and `G_over_aH_drmd_` already exist as input-parameter members with accessors `z_stop()`/`G_over_aH_drmd()` — do not duplicate them.)

- [ ] **Step 2: Implement the species methods**

In `species/idm_drmd_idr_drmd_species.cpp` add:
```cpp
void IDM_DRMD_IDR_DRMD_Species::RegisterBackgroundIndices(int& index_bg) {
  CompositeSpecies::RegisterBackgroundIndices(index_bg);  // children: idm_drmd, idr_drmd
  class_define_index(index_bg_G_over_aH_drmd_, _TRUE_, index_bg, 1);
}

void IDM_DRMD_IDR_DRMD_Species::ComputeIdmDrmd(double a, double rho_idm_over_rho_idr,
                                               double* Rint, double* csp2, double* Gint) const {
  const double z         = 1.0 / a - 1.0;
  const double R_int_tmp = 3.0 / 4.0 * rho_idm_over_rho_idr;
  *Rint = R_int_tmp;
  *csp2 = 1.0 / 3.0 / (1.0 + R_int_tmp);
  if ((1.0 + z_stop_) / (1.0 + z) > 100)  // avoid exp() overflow
    *Gint = 0.;
  else
    *Gint = Gamma0_drmd_ic_ / R_int_tmp * exp(-(1.0 + z_stop_) / (1.0 + z));
}

void IDM_DRMD_IDR_DRMD_Species::InitializeDrmdBackground(double rho_tot, double H,
                                                         double a, double a_today) {
  const double rho_idr = idr_drmd_->Rho_at_a(a);   // see note below
  const double rho_idm = idm_drmd_->Rho_at_a(a);
  f_idr_drmd_     = rho_idr / rho_tot;
  Gamma0_drmd_ic_ = 0.;
  if (rho_idm > 0. && rho_idr > 0.)
    Gamma0_drmd_ic_ = 3. / 4. * G_over_aH_drmd_ * rho_idm / rho_idr * a / a_today * H;
}

void IDM_DRMD_IDR_DRMD_Species::FinalizeBackground(double a, double H,
                                                   const double* /*pvecback_B*/, double* pvecback) {
  double Rint, csp2, Gint;
  const double a_rel = a / pba_.a_today;
  ComputeIdmDrmd(a, idm_drmd_->Rho(pvecback) / idr_drmd_->Rho(pvecback), &Rint, &csp2, &Gint);
  pvecback[index_bg_G_over_aH_drmd_] = Gint / (H * a_rel);
}
```
**Note on `InitializeDrmdBackground`:** the original module code read the densities from `pvecback` (`drmd_ic.idr_drmd().Rho(pvecback)`). Keep that to stay byte-identical — pass `pvecback` instead of re-deriving:
```cpp
  void InitializeDrmdBackground(double rho_tot, double H, double a, double a_today,
                                const double* pvecback);
```
with body using `idr_drmd_->Rho(pvecback)` / `idm_drmd_->Rho(pvecback)` (drop the `Rho_at_a` calls — that accessor may not exist; the `Rho(pvecback)` path is what the module used). Adjust Step 3's call accordingly.

- [ ] **Step 3: Module — thin wrapper, delete hardcoded block, delegate IC, update banner**

In `source/background_module.cpp`:
- Replace `background_idm_drmd`'s body (`:462-473`) with a delegating wrapper:
```cpp
void BackgroundModule::background_idm_drmd(
    double a, double rho_idm_over_rho_idr, double* Rint, double* csp2, double* Gint) const {
  static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"))
      .ComputeIdmDrmd(a, rho_idm_over_rho_idr, Rint, csp2, Gint);
}
```
- Delete the per-step DRMD block in `background_functions` (`:402-411`) — `FinalizeBackground` now writes `G_over_aH`.
- Replace the IC block (`:1059-1072`) with:
```cpp
  if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
    static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"))
        .InitializeDrmdBackground(pvecback[index_bg_rho_tot_], pvecback[index_bg_H_],
                                  a, pba->a_today, pvecback);
  }
```
- In the verbose banner (`:949-950` region), replace `Gamma0_drmd_`/`f_idr_drmd_` reads with the species accessors (fetch `auto& drmd = static_cast<…>(*all_species_.at("IDM_DRMD_IDR_DRMD"));` and use a public `Gamma0_drmd_ic()` accessor — add one returning `Gamma0_drmd_ic_` — and `drmd.f_idr_drmd()`).
- Remove the tail `class_define_index(index_bg_G_over_aH_drmd_, …)` added in Task 3 Step 2 (now species-owned).

In `source/background_module.h`: delete `index_bg_G_over_aH_drmd_` (`:52`), `Gamma0_drmd_` (`:114`), `f_idr_drmd_`, and the module `z_stop_` if module-local (keep `z_dec_drmd_`/`G_over_aH_tmp_` until Task 6). Keep the `background_idm_drmd` wrapper declaration (`:22`).

- [ ] **Step 4: Point the species' own perturbation-side calls at `ComputeIdmDrmd`**

In `species/idm_drmd_idr_drmd_species.cpp`, the two self-calls through `bgm->background_idm_drmd(...)` (`:59`, `:337`) can call `ComputeIdmDrmd(...)` directly (same result via the wrapper; direct is cleaner). The `perturbations_module.cpp` callers (`:3363`, `:4256`) keep using `background_module_->background_idm_drmd(...)` (the wrapper) — leave them.

- [ ] **Step 5: Build**

Run: `pip install . --no-build-isolation`
Expected: success. `grep -rn "index_bg_G_over_aH_drmd_" source` must show only the wrapper-free module (zero hits); the symbol now lives in the species.

- [ ] **Step 6: Golden test — exact PASS**

Run: `cd python && python -m pytest test_background_columns.py -v`
Expected: all PASS, `idm_drmd_idr_drmd` included (G_over_aH computed by the identical formula, just relocated).

- [ ] **Step 7: Commit**

```bash
git add species/idm_drmd_idr_drmd_species.cpp species/idm_drmd_idr_drmd_species.h source/background_module.cpp source/background_module.h
git commit -m "refactor(#319): relocate DRMD interaction physics + G_over_aH index into the species"
```

---

## Task 6: `ProcessBackgroundTable()` hook + relocate the decoupling-redshift scan

Move the post-table-fill decoupling scan (`z_dec_drmd_`/`G_over_aH_tmp_`) out of the module's generic table loop into a species hook. `z_dec_drmd_` feeds only a verbose printf, so this is byte-identical.

**Files:**
- Modify: `species/base_species.h` (add the virtual)
- Modify: `species/idm_drmd_idr_drmd_species.{h,cpp}`
- Modify: `source/background_module.{h,cpp}`

- [ ] **Step 1: Add the base-class hook**

In `species/base_species.h`, after `WriteBackgroundData`, add:
```cpp
  /**
   * Called once after the full background table is built. Lets a species run
   * table-scope analysis over its own columns. Default no-op.
   */
  virtual void ProcessBackgroundTable(const double* background_table, int n_rows,
                                      int row_stride, const double* z_table) {}
```

- [ ] **Step 2: Species owns the scan + result**

In `species/idm_drmd_idr_drmd_species.h` `private:` add:
```cpp
  double z_dec_drmd_    = -1.0;   // decoupling redshift (G_over_aH closest to 1); -1 = none
  double G_over_aH_tmp_ = 1e20;   // running |G_over_aH - 1| scan tracker
```
`public:` add:
```cpp
  double z_dec_drmd() const { return z_dec_drmd_; }
  void ProcessBackgroundTable(const double* background_table, int n_rows,
                              int row_stride, const double* z_table) override;
```
In `species/idm_drmd_idr_drmd_species.cpp` add:
```cpp
void IDM_DRMD_IDR_DRMD_Species::ProcessBackgroundTable(const double* background_table, int n_rows,
                                                       int row_stride, const double* z_table) {
  for (int i = 0; i < n_rows; i++) {
    const double g = background_table[i * row_stride + index_bg_G_over_aH_drmd_];
    if (pow(g - 1.0, 2.0) < pow(G_over_aH_tmp_ - 1.0, 2.0)) {
      G_over_aH_tmp_ = g;
      z_dec_drmd_    = z_table[i];
    }
  }
}
```

- [ ] **Step 3: Module — call the hook, delete the in-loop scan + init, update banner**

In `source/background_module.cpp`:
- Delete the DRMD scan inside the table loop (`:888-894`).
- Delete the `G_over_aH_tmp_ = 1e20;` / `z_dec_drmd_ = -1.0;` init block (`:524-531` region).
- After the background table is fully built (just after the table-fill loop completes, before the verbose-output section in `background_solve`/`background_solve_evolver`), add:
```cpp
  for (auto& [name, sp] : all_species_)
    sp->ProcessBackgroundTable(background_table_.data(), bt_size_, bg_size_, z_table_.data());
```
- In the decoupling-redshift banner (`:956-957`), read via accessor:
```cpp
      auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
      if (drmd.z_dec_drmd() > 0)
        printf("     -> decoupling occurred at z=%f \n", drmd.z_dec_drmd());
```

In `source/background_module.h`: delete `z_dec_drmd_` (`:117`) and `G_over_aH_tmp_` (`:113`).

- [ ] **Step 4: Build**

Run: `pip install . --no-build-isolation`
Expected: success. `grep -rn "z_dec_drmd_\|G_over_aH_tmp_" source` must show zero hits.

- [ ] **Step 5: Golden test — exact PASS**

Run: `cd python && python -m pytest test_background_columns.py -v`
Expected: all PASS (scan result feeds only a printf; numerical output unchanged).

- [ ] **Step 6: Commit**

```bash
git add species/base_species.h species/idm_drmd_idr_drmd_species.cpp species/idm_drmd_idr_drmd_species.h source/background_module.cpp source/background_module.h
git commit -m "refactor(#319): ProcessBackgroundTable hook; DRMD decoupling scan moves to the species"
```

---

## Task 7: Full-suite verification

No code change unless a regression surfaces. Confirms the whole refactor against the regenerated `classyref`.

- [ ] **Step 1: Full characterization + core unit tests**

Run: `cd python && python -m pytest test_background_columns.py test_transfer_columns.py test_class.py -q`
Expected: all PASS (transfer golden from #309 untouched; background golden green; core class tests green).

- [ ] **Step 2: Full scenario suite, both gauges, against fresh classyref**

Run: `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -m test_scenario test_class.py -q`
Expected: all PASS. (If any newtonian fail appears, it is a real regression — not stale classyref, which was regenerated 2026-06-17.)

- [ ] **Step 3: Confirm no dead symbols remain**

Run:
```bash
grep -rn "DpDloga\|index_bg_number_ncdm1_\|index_bg_pseudo_p_ncdm1_\|index_bg_Gamma0_drmd_\|ComputePPrimeAndWrite" source species
```
Expected: zero hits.
Run:
```bash
grep -rn "index_bg_phi_scf_\|index_bg_G_over_aH_drmd_\|z_dec_drmd_" source
```
Expected: zero hits in `source/` (these symbols now live only in `species/`).

- [ ] **Step 4: Final commit (if any baseline/test adjustments were needed)**

```bash
git add -A
git commit -m "test(#319): full-suite verification green against regenerated classyref"
```

---

## Self-Review

**Spec coverage:**
- §A flat loop → Task 3. §B dead indices → Task 3. §C scf ownership → Task 2 (reads) + Task 3 (delete copies). §D1 PPrime → Task 4. §D2 FinalizeBackground (scf + DRMD) → Task 4 (scf) + Task 5 (DRMD). §E DRMD physics relocation → Task 5. §F1 ProcessBackgroundTable → Task 6. Verification → Tasks 1, 4, 7. ✔ All sections covered.
- **Spec correction folded in:** §E loosely said the DRMD IC constants compute in `SetBackgroundInitialConditions`; that hook runs (line 1038) *before* `pvecback` has `rho_tot`/`H` (populated at line 1042), so the constants instead compute in a dedicated `InitializeDrmdBackground(...)` the module calls at the post-population point (line 1059). Same intent (physics + state on the species), accurate mechanism.

**Placeholder scan:** No TBD/TODO; every code step shows the code. The one "see note" in Task 5 Step 2 is an explicit correction (use `Rho(pvecback)` to stay byte-identical, with the adjusted signature) — resolved inline, not deferred.

**Type/signature consistency:** `PPrime(double a, double H, const double* pvecback_B, const double* pvecback) const` is identical across base, all overrides, and the call site. `FinalizeBackground(double a, double H, const double* pvecback_B, double* pvecback)` consistent (base default no-op, scf + DRMD override). `ProcessBackgroundTable(const double*, int, int, const double*)` consistent. `InitializeDrmdBackground` signature reconciled to include `const double* pvecback` (Task 5 Step 2 note → Step 3 call passes it). `bg_G_over_aH_index()`, `f_idr_drmd()`, `z_dec_drmd()`, `Gamma0_drmd_ic()` accessors named consistently between definition and use.
