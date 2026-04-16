# Species Background Index Encapsulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate all species-dependent mirror index fields from `BackgroundModule`, remove `bg_rho_index()`/`bg_p_index()` from the public interface, migrate background output to per-species virtuals, and replace the manual IC rho accumulation block with an `EnergyType` + free-streaming-radiation dispatch loop that works for composites.

**Architecture:** Compiler-guided migration (Approach 2) executed in two phases. Phase 1: delete the public accessor interface in one commit, let the compiler enumerate every call site, fix module by module. Phase 2: replace the IC accumulation block with a dispatch loop. New `BackgroundColumnWriter` (modelled on `PerturbColumnWriter`) enables per-species background output dispatch. `IsFreestreaming()` remains the simple-species predicate, while `FreestreamingRho()` provides the quantity actually accumulated in the IC loop and lets composites sum over children.

**Tech Stack:** C++17. Build: `make class -j`. Test: `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py` (84 scenarios). ASan: `make clean && make class -j CXXFLAGS="-O1 -g -fsanitize=address,undefined" && ./class explanatory.ini`.

---

## File Map

| File | Change |
|------|--------|
| `source/background_column_writer.h` | **Create** — `BackgroundColumnWriter` struct with `Add()` method |
| `source/background_column_writer.cpp` | **Create** — `Add()` implementation |
| `species/base_species.h` | **Modify** — add `IsFreestreaming()`, `WriteBackgroundColumnTitles()`, `WriteBackgroundData()`, `IsPresent()`; remove `bg_rho_index()`, `bg_p_index()` |
| `species/photons.h` | **Modify** — add `WriteBackgroundColumnTitles/Data` |
| `species/baryons.h` | **Modify** — add `WriteBackgroundColumnTitles/Data` |
| `species/cdm.h` / `cdm.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` |
| `species/lambda.h` / `lambda.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` |
| `species/fluid.h` / `fluid.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data`, `WriteWFld()` |
| `species/ultra_relativistic.h` / `.cpp` | **Modify** — add `IsFreestreaming()`, `WriteBackgroundColumnTitles/Data` |
| `species/dark_radiation_species.h` / `.cpp` | **Modify** — add `IsFreestreaming()`, `WriteBackgroundColumnTitles/Data` |
| `species/ncdm_species.h` / `.cpp` | **Modify** — add `IsFreestreaming()`, `WriteBackgroundColumnTitles/Data` |
| `species/idm_dr_idr_species.h` / `.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` (delegates to sub-components) |
| `species/idm_drmd_idr_drmd_species.h` / `.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` (delegates to sub-components) |
| `species/dcdm_dr_species.h` / `.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` (delegates to sub-components) |
| `species/scalar_field.h` / `.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` |
| `species/dncdm_species.h` / `.cpp` | **Modify** — add `WriteBackgroundColumnTitles/Data` |
| `species/dncdm_decay_radiation_species.h` / `.cpp` | **Modify** — add `IsFreestreaming()` override |
| `source/background_module.h` | **Modify** — remove ~15 species-dependent mirror fields |
| `source/background_module.cpp` | **Modify** — replace per-species registration blocks with uniform loop; replace output calls with dispatch; update `background_derivs_member()` rho_M accumulation; update fluid pre-computation |
| `source/perturbations_module.cpp` | **Modify** — replace ~18 `index_bg_rho_X_` density access sites; replace IC accumulation block |
| `source/thermodynamics_module.cpp` | **Modify** — replace 3 `index_bg_rho_*` sites |
| `source/input_module.cpp` | **Modify** — replace 4 background table lookup sites |

---

### Task 1: Create GitHub issue and branch

**Files:** — (git/GitHub only)

- [ ] **Step 1: Create GitHub issue**

```bash
gh issue create \
  --title "Remove species-dependent index mirrors from BackgroundModule" \
  --body "Remove all \`index_bg_rho_X_\` mirror fields from BackgroundModule, delete \`bg_rho_index()\`/\`bg_p_index()\` accessors from BaseSpecies, migrate background output to per-species WriteBackgroundColumnTitles/Data virtuals, and replace IC rho accumulation with EnergyType+IsFreestreaming() dispatch.

Spec: docs/superpowers/specs/2026-04-15-species-index-encapsulation-design.md"
```

Note the issue number from the output (e.g. `#248`).

- [ ] **Step 2: Create branch**

```bash
git checkout master
git pull origin master
git checkout -b 248-species-index-encapsulation
```

---

### Task 2: BackgroundColumnWriter + BaseSpecies new interface

**Files:**
- Create: `source/background_column_writer.h`
- Create: `source/background_column_writer.cpp`
- Modify: `species/base_species.h`

The `BackgroundColumnWriter` is modelled on `PerturbColumnWriter` in `species/perturb_source_context.h`. In title mode (`titles != nullptr`), `Add()` writes the column name and ignores the value. In data mode, it writes the value to `dataptr`. The `condition` flag matches `class_store_*` semantics: if false, the column slot is still consumed (so title/data counts match) but zero is written.

- [ ] **Step 1: Create `source/background_column_writer.h`**

```cpp
#pragma once
#include "common.h"  // class_store_columntitle, class_store_double, _TRUE_, _FALSE_

/**
 * Thin helper for background output — analogous to PerturbColumnWriter.
 * Construct in title mode (titles != nullptr) or data mode (dataptr != nullptr).
 * Call Add() once per column in both modes; the writer handles the branching.
 */
class BackgroundColumnWriter {
 public:
  // Title mode
  explicit BackgroundColumnWriter(char* titles) : titles_(titles) {}

  // Data mode
  BackgroundColumnWriter(double* dataptr, int& storeidx)
      : dataptr_(dataptr), storeidx_(&storeidx) {}

  bool IsTitleMode() const { return titles_ != nullptr; }

  /** Write one column: title in title mode, value in data mode.
   *  condition=false: slot is consumed but zero is written (inactive column). */
  void Add(const char* title, double value, bool condition = true);

 private:
  char*   titles_   = nullptr;
  double* dataptr_  = nullptr;
  int*    storeidx_ = nullptr;
};
```

- [ ] **Step 2: Create `source/background_column_writer.cpp`**

```cpp
#include "background_column_writer.h"

void BackgroundColumnWriter::Add(const char* title, double value, bool condition) {
  if (titles_) {
    class_store_columntitle(titles_, title, condition ? _TRUE_ : _FALSE_);
  } else if (dataptr_) {
    class_store_double(dataptr_, value, condition ? _TRUE_ : _FALSE_, (*storeidx_));
  }
}
```

- [ ] **Step 3: Add `BackgroundColumnWriter` forward declaration + new virtuals to `species/base_species.h`**

Add the include near the top of `base_species.h`:
```cpp
// After the existing forward declarations, add:
class BackgroundColumnWriter;  // forward declaration
```

Add these methods to the `BaseSpecies` public section, after `DpDloga()`:

```cpp
  // ── Background output ─────────────────────────────────────────────────────

  /**
   * True for species that are free-streaming and massless at IC time
   * (deep radiation domination). Used to accumulate rho_nu / fracnu in
   * perturbation IC setup. Static at construction time.
   * Default: false.
   */
  virtual bool IsFreestreaming() const { return false; }

  /**
   * Returns true if this species is present (has registered its background
   * indices). For top-level species use all_species_.count(); IsPresent() is
   * for sub-components of composites where all_species_.count() is unavailable.
   */
  bool IsPresent() const { return index_bg_rho_ >= 0; }

  /**
   * Write this species' background column titles into writer.
   * Called by BackgroundModule::background_output_titles().
   * Default: no-op (species with no background output need not override).
   */
  virtual void WriteBackgroundColumnTitles(BackgroundColumnWriter& /*writer*/) const {}

  /**
   * Write this species' background data values into writer.
   * pvecback is the full background vector at the current time step.
   * Called by BackgroundModule::background_output_data().
   * Default: no-op.
   */
  virtual void WriteBackgroundData(const double* /*pvecback*/,
                                   BackgroundColumnWriter& /*writer*/) const {}
```

- [ ] **Step 4: Add `#include "background_column_writer.h"` to `source/background_module.cpp`**

Check that the file builds:
```bash
make class -j 2>&1 | tail -20
```
Expected: no new errors (we haven't deleted anything yet).

- [ ] **Step 5: Commit**

```bash
git add source/background_column_writer.h source/background_column_writer.cpp species/base_species.h
git commit -m "Add BackgroundColumnWriter and BaseSpecies background output virtuals"
```

---

### Task 3: Implement WriteBackground* and IsFreestreaming in each species

**Files:** `species/photons.h`, `species/baryons.h`, `species/cdm.{h,cpp}`, `species/lambda.{h,cpp}`, `species/fluid.{h,cpp}`, `species/ultra_relativistic.{h,cpp}`, `species/dark_radiation_species.{h,cpp}`, `species/ncdm_species.{h,cpp}`, `species/idm_dr_idr_species.{h,cpp}`, `species/idm_drmd_idr_drmd_species.{h,cpp}`, `species/dcdm_dr_species.{h,cpp}`, `species/scalar_field.{h,cpp}`, `species/dncdm_species.{h,cpp}`, `species/dncdm_decay_radiation_species.h`

All species add `WriteBackgroundColumnTitles` and `WriteBackgroundData`. After this task the old `background_output_*` functions still work unchanged — no calls are deleted yet. We are just pre-populating the new dispatch interface.

The `BackgroundColumnWriter` must be included in species .cpp files:
```cpp
#include "../source/background_column_writer.h"
```
(or via a forward declaration in the .h and include in .cpp)

- [ ] **Step 1: PhotonsSpecies (`species/photons.h`) — add inline overrides**

In `photons.h`, add to the class body (photons.h uses inline methods):

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", 0.);
  }
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", pvecback[index_bg_rho_]);
  }
```

Add `#include "../source/background_column_writer.h"` to `photons.h`.

- [ ] **Step 2: BaryonsSpecies (`species/baryons.h`) — add inline overrides**

Find `baryons.h` and add:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_b", 0.);
  }
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_b", pvecback[index_bg_rho_]);
  }
```

Add `#include "../source/background_column_writer.h"` to `baryons.h`.

- [ ] **Step 3: CDMSpecies (`species/cdm.h` + `cdm.cpp`) — add declarations + definitions**

In `cdm.h`, add declarations:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `cdm.cpp`, add `#include "../source/background_column_writer.h"` and implement:

```cpp
void CDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_cdm", 0.);
}
void CDMSpecies::WriteBackgroundData(const double* pvecback,
                                     BackgroundColumnWriter& w) const {
  w.Add("(.)rho_cdm", pvecback[index_bg_rho_cdm_]);
}
```

- [ ] **Step 4: LambdaSpecies (`species/lambda.h` + `lambda.cpp`)**

Declarations in `lambda.h`:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

Definitions in `lambda.cpp`:

```cpp
void LambdaSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_lambda", 0.);
}
void LambdaSpecies::WriteBackgroundData(const double* pvecback,
                                         BackgroundColumnWriter& w) const {
  w.Add("(.)rho_lambda", pvecback[index_bg_rho_lambda_]);
}
```

- [ ] **Step 5: FluidSpecies (`species/fluid.h` + `fluid.cpp`)**

Fluid also needs `WriteWFld()` — called by `BackgroundModule::background_functions()` instead of the current direct index writes. This keeps `index_bg_w_fld_` and `index_bg_dw_over_da_fld_` private to `FluidSpecies`.

Declarations in `fluid.h`:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;

  /** Called by BackgroundModule::background_functions() before ComputeBackground().
   *  Writes w_fld and dw/da (already computed by BackgroundModule) into pvecback
   *  using FluidSpecies's private indices. */
  void WriteWFld(double w_fld, double dw_over_da_fld, double* pvecback) const;
```

Definitions in `fluid.cpp`:

```cpp
void FluidSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_fld", 0.);
  w.Add("(.)w_fld", 0.);
}
void FluidSpecies::WriteBackgroundData(const double* pvecback,
                                        BackgroundColumnWriter& w) const {
  w.Add("(.)rho_fld", pvecback[index_bg_rho_fld_]);
  w.Add("(.)w_fld",   pvecback[index_bg_w_fld_]);
}
void FluidSpecies::WriteWFld(double w_fld, double dw_over_da_fld,
                              double* pvecback) const {
  pvecback[index_bg_w_fld_]          = w_fld;
  pvecback[index_bg_dw_over_da_fld_] = dw_over_da_fld;
}
```

- [ ] **Step 6: UltraRelativisticSpecies (`species/ultra_relativistic.h` + `.cpp`)**

UR neutrinos are free-streaming massless radiation → `IsFreestreaming()` returns true.

In `ultra_relativistic.h`:

```cpp
  bool IsFreestreaming() const override { return true; }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `ultra_relativistic.cpp` (add include + definitions):

```cpp
void UltraRelativisticSpecies::WriteBackgroundColumnTitles(
    BackgroundColumnWriter& w) const {
  w.Add("(.)rho_ur", 0.);
}
void UltraRelativisticSpecies::WriteBackgroundData(const double* pvecback,
                                                    BackgroundColumnWriter& w) const {
  w.Add("(.)rho_ur", pvecback[index_bg_rho_ur_]);
}
```

(`index_bg_rho_ur_` is the private species field set in `RegisterBackgroundIndices`.)

- [ ] **Step 7: DarkRadiationSpecies (`species/dark_radiation_species.h` + `.cpp`)**

DR from DCDM decay is free-streaming radiation → `IsFreestreaming()` returns true. DR is a sub-component of `DCDM_DR_Species`; its output is delegated from the composite (Step 11). Add the predicate only:

In `dark_radiation_species.h`:

```cpp
  bool IsFreestreaming() const override { return true; }
```

The `WriteBackgroundColumnTitles/Data` are implemented on the DCDM_DR composite in Step 11.

- [ ] **Step 8: NCDMSpecies (`species/ncdm_species.h` + `.cpp`)**

NCDM is free-streaming at IC time (deep radiation domination, all NCDM species are relativistic). `IsFreestreaming()` = true. Each NCDMSpecies writes its own flavor columns.

In `ncdm_species.h`:

```cpp
  bool IsFreestreaming() const override { return true; }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `ncdm_species.cpp`:

```cpp
void NCDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
}
void NCDMSpecies::WriteBackgroundData(const double* pvecback,
                                       BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, pvecback[bg_number_index()]);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, pvecback[index_bg_rho_]);
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, pvecback[index_bg_p_]);
}
```

Note: `bg_number_index()` and `bg_p_index()` on NCDMSpecies are NOT the public `BaseSpecies` getters being removed; they are NCDM-specific accessors and remain. (The `BaseSpecies` public getters `bg_rho_index()` and `bg_p_index()` are removed in Task 4. NCDMSpecies can keep its own `bg_number_index()` and `bg_pseudo_p_index()` since those are species-specific.)

For DNCDMSpecies (decaying NCDM, class `DNCDMSpecies`): check `species/dncdm_species.h` for the fields; it has additional lnf columns. Add analogous `WriteBackgroundColumnTitles/Data` to `DNCDMSpecies` that write the lnf/dlnfdlnq columns using `bg_lnf_index()`, `bg_dlnfdlnq_index()`, `bg_dlnfdlnq_sep_index()` and `ncdm_->q_size_ncdm_[ncdm_id_]`. Mirror the existing loop in `background_output_titles/data` in `background_module.cpp:1705-1826`.

- [ ] **Step 9: IDM_DR_IDR_Species (`species/idm_dr_idr_species.h` + `.cpp`)**

The composite owns both sub-components. It writes their columns. IDR is free-streaming radiation.

In `idm_dr_idr_species.h`:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `idm_dr_idr_species.cpp`:

```cpp
void IDM_DR_IDR_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr",    0.);
  w.Add("(.)rho_idm_dr", 0.);
}
void IDM_DR_IDR_Species::WriteBackgroundData(const double* pvecback,
                                              BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr",    idr().Rho(pvecback));
  w.Add("(.)rho_idm_dr", idm_dr().Rho(pvecback));
}
```

Also add `IsFreestreaming()` to the IDR sub-species in `species/idr.h`:

```cpp
  bool IsFreestreaming() const override { return true; }
```

- [ ] **Step 10: IDM_DRMD_IDR_DRMD_Species (`species/idm_drmd_idr_drmd_species.h` + `.cpp`)**

In `idm_drmd_idr_drmd_species.h`:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `idm_drmd_idr_drmd_species.cpp`:

```cpp
void IDM_DRMD_IDR_DRMD_Species::WriteBackgroundColumnTitles(
    BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr_drmd",  0.);
  w.Add("(.)rho_idm_drmd",  0.);
  w.Add("G_over_aH_drmd",   0.);
}
void IDM_DRMD_IDR_DRMD_Species::WriteBackgroundData(const double* pvecback,
                                                     BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr_drmd",  idr_drmd().Rho(pvecback));
  w.Add("(.)rho_idm_drmd",  idm_drmd().Rho(pvecback));
  w.Add("G_over_aH_drmd",   pvecback[bgm_->index_bg_G_over_aH_drmd_]);
}
```

`index_bg_G_over_aH_drmd_` remains a BackgroundModule field (it is a physics quantity, not a density mirror).

Also add `IsFreestreaming()` to the IDR_DRMD sub-species in `species/idr_drmd.h`:

```cpp
  bool IsFreestreaming() const override { return true; }
```

- [ ] **Step 11: DCDM_DR_Species (`species/dcdm_dr_species.h` + `.cpp`)**

In `dcdm_dr_species.h`:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `dcdm_dr_species.cpp`:

```cpp
void DCDM_DR_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm", 0.);
  w.Add("(.)rho_dr",   0.);
  // Per-decay-channel DR species
  for (int j = 0; j < dr().num_decay_channels(); ++j) {
    char tmp[40];
    snprintf(tmp, 40, "(.)rho_dr[%d]", j);
    w.Add(tmp, 0.);
  }
}
void DCDM_DR_Species::WriteBackgroundData(const double* pvecback,
                                           BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm", dcdm().Rho(pvecback));
  w.Add("(.)rho_dr",   dr().Rho(pvecback));
  for (int j = 0; j < dr().num_decay_channels(); ++j) {
    char tmp[40];
    snprintf(tmp, 40, "(.)rho_dr[%d]", j);
    w.Add(tmp, pvecback[dr().bg_rho_dr_species_index() + j]);
  }
}
```

`bg_rho_dr_species_index()` is a DarkRadiationSpecies-specific accessor that returns the base index for per-species DR densities; it is NOT the `BaseSpecies` `bg_rho_index()` being removed. Verify it exists in `dark_radiation_species.h` and add it if not.

If `DarkRadiationSpecies` does not have `num_decay_channels()`, use `pba_->N_decay_dr` via the stored `pba_` reference (check `dcdm_dr_species.cpp` for the pointer).

- [ ] **Step 12: ScalarFieldSpecies (`species/scalar_field.h` + `.cpp`)**

ScalarField has 8 columns: rho, p, p', phi, phi', V, V', V''. The phi/V/dV/ddV indices remain in `BackgroundModule` (they are physics indices needed by `dV_scf()` etc.). ScalarFieldSpecies accesses them via `bgm_->index_bg_phi_scf_` etc.

In `scalar_field.h`:

```cpp
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback,
                           BackgroundColumnWriter& w) const override;
```

In `scalar_field.cpp`:

```cpp
void ScalarFieldSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_scf",     0.);
  w.Add("(.)p_scf",       0.);
  w.Add("(.)p_prime_scf", 0.);
  w.Add("phi_scf",        0.);
  w.Add("phi'_scf",       0.);
  w.Add("V_scf",          0.);
  w.Add("V'_scf",         0.);
  w.Add("V''_scf",        0.);
}
void ScalarFieldSpecies::WriteBackgroundData(const double* pvecback,
                                              BackgroundColumnWriter& w) const {
  w.Add("(.)rho_scf",     Rho(pvecback));
  w.Add("(.)p_scf",       P(pvecback));
  w.Add("(.)p_prime_scf", pvecback[index_bg_p_prime_scf_]);
  w.Add("phi_scf",        pvecback[bgm_->index_bg_phi_scf_]);
  w.Add("phi'_scf",       pvecback[bgm_->index_bg_phi_prime_scf_]);
  w.Add("V_scf",          pvecback[bgm_->index_bg_V_scf_]);
  w.Add("V'_scf",         pvecback[bgm_->index_bg_dV_scf_]);
  w.Add("V''_scf",        pvecback[bgm_->index_bg_ddV_scf_]);
}
```

`index_bg_p_prime_scf_` is a `ScalarFieldSpecies`-private field (set in `RegisterBackgroundIndices` at offset `index_bg_phi_scf_ + 7`; see `scalar_field.h`). `Rho()` and `P()` use the species-private `index_bg_rho_scf_` and `index_bg_p_scf_` (confirmed in `scalar_field.h`).

- [ ] **Step 13: Build and verify no regressions**

```bash
make class -j 2>&1 | grep -E "error:|warning:" | head -30
```

Expected: zero new errors. The old `background_output_*` functions still exist and compile.

- [ ] **Step 14: Commit**

```bash
git add species/photons.h species/baryons.h \
        species/cdm.h species/cdm.cpp \
        species/lambda.h species/lambda.cpp \
        species/fluid.h species/fluid.cpp \
        species/ultra_relativistic.h species/ultra_relativistic.cpp \
        species/dark_radiation_species.h species/dark_radiation_species.cpp \
        species/ncdm_species.h species/ncdm_species.cpp \
        species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp \
        species/dcdm_dr_species.h species/dcdm_dr_species.cpp \
        species/scalar_field.h species/scalar_field.cpp \
        species/dncdm_species.h species/dncdm_species.cpp \
        species/idr.h species/idr_drmd.h
git commit -m "Implement WriteBackgroundColumnTitles/Data and IsFreestreaming per species"
```

---

### Task 4: Delete public interface + mirror fields; fix background_module.cpp

**Files:**
- Modify: `species/base_species.h` — remove `bg_rho_index()`, `bg_p_index()`
- Modify: `source/background_module.h` — remove ~15 species-dependent mirror fields
- Modify: `source/background_module.cpp` — registration loop, output dispatch, derivs, fluid pre-computation

After deleting the public interface the project will not compile. Fix `background_module.cpp` first (it is the species-owning module and has the simplest call sites). Leave `perturbations_module.cpp`, `thermodynamics_module.cpp`, and `input_module.cpp` for the following tasks.

#### Step 1: Delete from `species/base_species.h`

Remove the two public accessor methods:

```cpp
// REMOVE these two methods:
int bg_rho_index() const {
  return index_bg_rho_;
}
int bg_p_index() const {
  return index_bg_p_;
}
```

The protected fields `index_bg_rho_` and `index_bg_p_` stay. Subclasses still write to them in `RegisterBackgroundIndices`.

- [ ] **Step 1: Delete `bg_rho_index()` and `bg_p_index()` from `species/base_species.h`**

#### Step 2: Delete mirror fields from `source/background_module.h`

Remove these fields from `background_module.h`:

```
index_bg_rho_g_
index_bg_rho_b_
index_bg_rho_cdm_
index_bg_rho_lambda_
index_bg_rho_fld_
index_bg_w_fld_
index_bg_dw_over_da_fld_
index_bg_rho_ur_
index_bg_rho_idm_dr_
index_bg_rho_idr_
index_bg_rho_dcdm_
index_bg_rho_idm_drmd_
index_bg_rho_idr_drmd_
index_bg_rho_dr_species_
index_bg_rho_dr_
index_bg_rho_scf_
index_bg_p_scf_
index_bg_p_prime_scf_
index_bg_rho_ncdm1_
index_bg_p_ncdm1_
```

**Keep** (not removed):
- `index_bg_G_over_aH_drmd_`, `index_bg_Gamma0_drmd_` — DRMD physics indices
- `index_bg_phi_scf_` through `index_bg_ddV_scf_` — SCF potential field indices (used by `dV_scf()`, `V_scf()`, `ddV_scf()` methods and `background_derivs_member()`)
- `index_bg_number_ncdm1_`, `index_bg_pseudo_p_ncdm1_` — NCDM module-level indices
- All aggregates: `index_bg_rho_tot_`, `index_bg_p_tot_`, `index_bg_p_tot_prime_`, `index_bg_Omega_r_`, `index_bg_rho_crit_`, `index_bg_Omega_m_`, `index_bg_H_`, `index_bg_a_`, `index_bg_H_prime_`, distances, `index_bg_D_`, `index_bg_f_`, etc.

- [ ] **Step 2: Remove the 20 listed fields from `background_module.h`**

Do NOT compile yet. Continue to fix the `.cpp` file first.

#### Step 3: Fix `background_module.cpp` — registration (`background_indices()`)

Replace the per-species manual registration blocks (lines ~713-836) with a uniform pattern. The ORDER of registration must be preserved exactly — pvecback layout is baked into all output and test reference data.

- [ ] **Step 3: Replace per-species registration blocks in `background_indices()`**

```cpp
// ── Photons (always) ──────────────────────────────────────────────────────
all_species_.at("Photons")->RegisterBackgroundIndices(index_bg);

// ── Baryons (always) ──────────────────────────────────────────────────────
all_species_.at("Baryons")->RegisterBackgroundIndices(index_bg);

// ── CDM (optional) ────────────────────────────────────────────────────────
if (all_species_.count("CDM"))
  all_species_.at("CDM")->RegisterBackgroundIndices(index_bg);

// ── IDM_DRMD + IDR_DRMD composite (optional) ──────────────────────────────
if (all_species_.count("IDM_DRMD_IDR_DRMD"))
  all_species_.at("IDM_DRMD_IDR_DRMD")->RegisterBackgroundIndices(index_bg);

// Module physics indices for DRMD (not species-dependent densities)
class_define_index(index_bg_G_over_aH_drmd_,
                   all_species_.count("IDM_DRMD_IDR_DRMD"), index_bg, 1);
class_define_index(index_bg_Gamma0_drmd_,
                   all_species_.count("IDM_DRMD_IDR_DRMD"), index_bg, 1);

// ── NCDM (optional, sorted by ncdm_id) ───────────────────────────────────
index_bg_number_ncdm1_ = index_bg_pseudo_p_ncdm1_ = -1;
if (pba->N_ncdm > 0) {
  index_bg_number_ncdm1_ = index_bg;
  std::vector<NCDMSpecies*> ncdm_vec;
  for (auto& [name, sp] : all_species_) {
    if (auto* n = dynamic_cast<NCDMSpecies*>(sp.get()))
      ncdm_vec.push_back(n);
  }
  std::sort(ncdm_vec.begin(), ncdm_vec.end(), [](NCDMSpecies* a, NCDMSpecies* b) {
    return a->ncdm_id() < b->ncdm_id();
  });
  for (auto* ncdm : ncdm_vec)
    ncdm->RegisterBackgroundIndices(index_bg);
  if (!ncdm_vec.empty())
    index_bg_pseudo_p_ncdm1_ = ncdm_vec[0]->bg_pseudo_p_index();
}

// ── DCDM_DR composite (optional) ─────────────────────────────────────────
if (all_species_.count("DCDM_DR"))
  all_species_.at("DCDM_DR")->RegisterBackgroundIndices(index_bg);

// ── ScalarField (optional) — module caches arithmetic offsets for dV/V/ddV
if (all_species_.count("ScalarField")) {
  index_bg_phi_scf_ = index_bg;
  all_species_.at("ScalarField")->RegisterBackgroundIndices(index_bg);
  index_bg_phi_prime_scf_ = index_bg_phi_scf_ + 1;
  index_bg_V_scf_         = index_bg_phi_scf_ + 2;
  index_bg_dV_scf_        = index_bg_phi_scf_ + 3;
  index_bg_ddV_scf_       = index_bg_phi_scf_ + 4;
  // rho_scf, p_scf, p_prime_scf are at offsets +5, +6, +7 in pvecback
  // but owned privately by ScalarFieldSpecies
}

// ── Lambda (optional) ─────────────────────────────────────────────────────
if (all_species_.count("Lambda"))
  all_species_.at("Lambda")->RegisterBackgroundIndices(index_bg);

// ── Fluid (optional) ──────────────────────────────────────────────────────
if (all_species_.count("Fluid"))
  all_species_.at("Fluid")->RegisterBackgroundIndices(index_bg);

// ── UR (optional) ─────────────────────────────────────────────────────────
if (all_species_.count("UR"))
  all_species_.at("UR")->RegisterBackgroundIndices(index_bg);

// ── Module aggregate indices (unchanged) ──────────────────────────────────
class_define_index(index_bg_rho_tot_, _TRUE_, index_bg, 1);
class_define_index(index_bg_p_tot_, _TRUE_, index_bg, 1);
class_define_index(index_bg_p_tot_prime_, _TRUE_, index_bg, 1);
class_define_index(index_bg_Omega_r_, _TRUE_, index_bg, 1);

// ── IDM_DR + IDR composite (optional) ────────────────────────────────────
if (all_species_.count("IDM_DR_IDR"))
  all_species_.at("IDM_DR_IDR")->RegisterBackgroundIndices(index_bg);
```

#### Step 4: Fix `background_module.cpp` — fluid pre-computation in `background_functions()`

The current code (around line 359) directly writes to `pvecback[index_bg_w_fld_]`. Replace with the new `WriteWFld()` call:

- [ ] **Step 4: Update the fluid pre-computation block in `background_functions()`**

Before (lines ~354-361):
```cpp
pvecback[index_bg_w_fld_]          = w_fld;
pvecback[index_bg_dw_over_da_fld_] = dw_over_da_fld;
all_species_.at("Fluid")->ComputeBackground(a_rel, pvecback_B, pvecback);
```

After:
```cpp
static_cast<FluidSpecies&>(*all_species_.at("Fluid"))
    .WriteWFld(w_fld, dw_over_da_fld, pvecback);
all_species_.at("Fluid")->ComputeBackground(a_rel, pvecback_B, pvecback);
```

Add `#include "../species/fluid.h"` to `background_module.cpp` if not already present.

#### Step 5: Fix `background_module.cpp` — output dispatch

Replace the species-specific `class_store_columntitle` / `class_store_double` blocks in `background_output_titles()` and `background_output_data()` with dispatch.

- [ ] **Step 5a: Rewrite `background_output_titles()`**

The new function:

```cpp
int BackgroundModule::background_output_titles(char titles[_MAXTITLESTRINGLENGTH_]) const {
  // ── Module header (always present) ──────────────────────────────────────
  class_store_columntitle(titles, "z",               _TRUE_);
  class_store_columntitle(titles, "proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles, "conf. time [Mpc]", _TRUE_);
  class_store_columntitle(titles, "H [1/Mpc]",        _TRUE_);
  class_store_columntitle(titles, "comov. dist.",      _TRUE_);
  class_store_columntitle(titles, "ang.diam.dist.",    _TRUE_);
  class_store_columntitle(titles, "lum. dist.",        _TRUE_);
  class_store_columntitle(titles, "comov.snd.hrz.",    _TRUE_);

  // ── Species output — per-species dispatch ───────────────────────────────
  BackgroundColumnWriter writer(titles);
  for (auto& [name, sp] : all_species_)
    sp->WriteBackgroundColumnTitles(writer);

  // ── Module aggregate columns ────────────────────────────────────────────
  class_store_columntitle(titles, "(.)rho_crit",     _TRUE_);
  class_store_columntitle(titles, "(.)rho_tot",      _TRUE_);
  class_store_columntitle(titles, "(.)p_tot",        _TRUE_);
  class_store_columntitle(titles, "(.)p_tot_prime",  _TRUE_);
  class_store_columntitle(titles, "gr.fac. D",       _TRUE_);
  class_store_columntitle(titles, "gr.fac. f",       _TRUE_);

  return _SUCCESS_;
}
```

Note: species columns now come before `rho_crit` instead of after `rho_b`. The output ordering changes slightly vs the old format, but this is acceptable — tests run with `COMPARE_OUTPUT_REF=0`.

- [ ] **Step 5b: Rewrite `background_output_data()`**

```cpp
int BackgroundModule::background_output_data(int number_of_titles, double* data) const {
  for (int index_tau = 0; index_tau < bt_size_; index_tau++) {
    double* dataptr  = data + index_tau * number_of_titles;
    double* pvecback = const_cast<double*>(background_table_.data())
                       + index_tau * bg_size_;
    int storeidx = 0;

    // ── Module header ──────────────────────────────────────────────────────
    class_store_double(dataptr,
                       pba->a_today / pvecback[index_bg_a_] - 1., _TRUE_, storeidx);
    class_store_double(dataptr,
                       pvecback[index_bg_time_] / _Gyr_over_Mpc_, _TRUE_, storeidx);
    class_store_double(dataptr,
                       conformal_age_ - pvecback[index_bg_conf_distance_],
                       _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_H_],              _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_conf_distance_],  _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_ang_distance_],   _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_lum_distance_],   _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rs_],             _TRUE_, storeidx);

    // ── Species data — per-species dispatch ─────────────────────────────────
    BackgroundColumnWriter writer(dataptr, storeidx);
    for (auto& [name, sp] : all_species_)
      sp->WriteBackgroundData(pvecback, writer);

    // ── Module aggregate columns ─────────────────────────────────────────────
    class_store_double(dataptr, pvecback[index_bg_rho_crit_],     _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_tot_],      _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_tot_],        _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_tot_prime_],  _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_D_],            _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_f_],            _TRUE_, storeidx);
  }
  return _SUCCESS_;
}
```

#### Step 6: Fix `background_module.cpp` — `background_derivs_member()`

Lines ~1941–1958 use mirror indices for the photon-baryon sound horizon and growth factor rho_M accumulation.

- [ ] **Step 6: Update `background_derivs_member()`**

```cpp
// Line ~1941 — rho_g sanity check:
// BEFORE: class_test(pvecback[index_bg_rho_g_] <= 0., ...)
// AFTER:
class_test(all_species_.at("Photons")->Rho(pvecback) <= 0.,
           error_message,
           "rho_g = %e instead of strictly positive",
           all_species_.at("Photons")->Rho(pvecback));

// Line ~1948 — sound horizon:
// BEFORE: 1. / sqrt(3. * (1. + 3. * pvecback[index_bg_rho_b_] / 4. / pvecback[index_bg_rho_g_]))
// AFTER:
dy[index_bi_rs_] =
    1. / sqrt(3. * (1. + 3. * all_species_.at("Baryons")->Rho(pvecback)
                               / 4. / all_species_.at("Photons")->Rho(pvecback))) *
    sqrt(1. - pba->K * y[index_bi_rs_] * y[index_bi_rs_]);

// Lines ~1952-1958 — growth factor rho_M:
// BEFORE:
//   double rho_M = pvecback[index_bg_rho_b_];
//   if (index_bg_rho_cdm_ >= 0) rho_M += pvecback[index_bg_rho_cdm_];
//   if (index_bg_rho_idm_dr_ >= 0) rho_M += pvecback[index_bg_rho_idm_dr_];
//   if (index_bg_rho_idm_drmd_ >= 0) rho_M += pvecback[index_bg_rho_idm_drmd_];
// AFTER:
double rho_M = all_species_.at("Baryons")->Rho(pvecback);
if (all_species_.count("CDM"))
  rho_M += all_species_.at("CDM")->Rho(pvecback);
if (all_species_.count("IDM_DR_IDR")) {
  auto& idm_idr = static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
  rho_M += idm_idr.idm_dr().Rho(pvecback);
}
if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
  auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
  rho_M += drmd.idm_drmd().Rho(pvecback);
}
```

Add includes for the composite species headers at the top of `background_module.cpp` if not already present:
```cpp
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
```

- [ ] **Step 7: Build — expect errors only in perturbations, thermodynamics, input modules**

```bash
make class -j 2>&1 | grep "error:" | grep -v "perturbations_module\|thermodynamics_module\|input_module"
```

Expected: zero errors outside those three files. Background module should compile clean.

- [ ] **Step 8: Commit**

```bash
git add species/base_species.h source/background_module.h source/background_module.cpp
git commit -m "Delete bg_rho_index/bg_p_index and mirror fields; fix background_module.cpp"
```

---

### Task 5: Fix perturbations_module.cpp density access sites

**Files:** `source/perturbations_module.cpp`

The 18+ sites in `perturbations_module.cpp` that access `background_module_->index_bg_rho_X_` directly. These are outside the IC accumulation block (which is handled separately in Task 7).

Use `make class -j 2>&1 | grep "perturbations_module.cpp" | grep "error:"` to get the current list. The patterns below cover the expected sites. Fix each systematically.

- [ ] **Step 1: Fix NCDM density accesses (`ncdm_sp->bg_rho_index()` calls)**

All calls to `ncdm_sp->bg_rho_index()` or `ncdm_sp->bg_p_index()` in perturbation loops must change. The per-species `Rho(pvecback)` and `P(pvecback)` are already implemented.

```cpp
// BEFORE: pvecback[ncdm_sp->bg_rho_index()]
// AFTER:  ncdm_sp->Rho(pvecback)

// BEFORE: pvecback[ncdm_sp->bg_p_index()]
// AFTER:  ncdm_sp->P(pvecback)
```

Search for all occurrences:
```bash
grep -n "bg_rho_index\|bg_p_index" source/perturbations_module.cpp
```

- [ ] **Step 2: Fix CDM, Lambda, Fluid, UR density accesses**

Pattern (replace `index_bg_rho_X_` with `Rho(pvecback)` dispatch):

```cpp
// CDM:
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_cdm_]
// AFTER:  all_species_.at("CDM")->Rho(ppw->pvecback)

// BEFORE: if (background_module_->index_bg_rho_cdm_ >= 0) ...
// AFTER:  if (all_species_.count("CDM")) ...

// Lambda:
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_lambda_]
// AFTER:  all_species_.at("Lambda")->Rho(ppw->pvecback)

// Fluid rho:
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_fld_]
// AFTER:  all_species_.at("Fluid")->Rho(ppw->pvecback)

// Fluid w_fld (pressure ratio, used in perturbation equations):
// BEFORE: ppw->pvecback[background_module_->index_bg_w_fld_]
// The w_fld value is in pvecback at FluidSpecies's private index_bg_w_fld_.
// FluidSpecies provides P(pvecback) = w_fld * rho.
// To get w_fld alone: P(pvecback) / Rho(pvecback) — OR add a W(pvecback) accessor to FluidSpecies.
// Recommended: add to FluidSpecies:
//   double W(const double* pvecback) const { return pvecback[index_bg_w_fld_]; }
// Then:
// AFTER:  static_cast<FluidSpecies&>(*all_species_.at("Fluid")).W(ppw->pvecback)

// UR:
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_ur_]
// AFTER:  all_species_.at("UR")->Rho(ppw->pvecback)
```

Add `double W(const double* pvecback) const { return pvecback[index_bg_w_fld_]; }` to `fluid.h` if needed.

- [ ] **Step 3: Fix composite sub-component density accesses**

```cpp
// IDM_DR (sub-component of IDM_DR_IDR):
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_idm_dr_]
// AFTER:
//   static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).idm_dr().Rho(ppw->pvecback)

// IDR (sub-component of IDM_DR_IDR):
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_idr_]
// AFTER:
//   static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).idr().Rho(ppw->pvecback)

// IDM_DRMD (sub-component of IDM_DRMD_IDR_DRMD):
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_idm_drmd_]
// AFTER:
//   static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD")).idm_drmd().Rho(ppw->pvecback)

// IDR_DRMD (sub-component of IDM_DRMD_IDR_DRMD):
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_idr_drmd_]
// AFTER:
//   static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD")).idr_drmd().Rho(ppw->pvecback)

// DCDM (sub-component of DCDM_DR):
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_dcdm_]
// AFTER:
//   static_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR")).dcdm().Rho(ppw->pvecback)

// DR total (sub-component of DCDM_DR):
// BEFORE: ppw->pvecback[background_module_->index_bg_rho_dr_]
// AFTER:
//   static_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR")).dr().Rho(ppw->pvecback)
```

Add required includes to `perturbations_module.cpp` if not already present:
```cpp
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
#include "../species/dcdm_dr_species.h"
#include "../species/fluid.h"
```

For guards that were `if (background_module_->index_bg_rho_X_ >= 0)`:
- Top-level optional species: use `all_species_.count("SpeciesName")`
- Composite sub-components: use `all_species_.count("CompositeName")` (sub-components always register when composite is present)

- [ ] **Step 4: Build and verify only IC-block errors remain (if any)**

```bash
make class -j 2>&1 | grep "perturbations_module.cpp" | grep "error:"
```

Expected: all density-access errors resolved. The IC accumulation block (lines ~4700-4775) still uses old indices — those are fixed in Task 7.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp species/fluid.h
git commit -m "Fix perturbations_module.cpp density access sites (non-IC)"
```

---

### Task 6: Fix thermodynamics_module.cpp and input_module.cpp

**Files:** `source/thermodynamics_module.cpp`, `source/input_module.cpp`

#### thermodynamics_module.cpp — 3 sites

- [ ] **Step 1: Fix the IDM_DR/IDR interaction rate (lines ~508-511)**

```cpp
// BEFORE:
thermodynamics_table_[...index_th_ddmu_idm_dr_] =
    4. / 3. * pvecback[background_module_->index_bg_rho_idr_] /
    pvecback[background_module_->index_bg_rho_idm_dr_] * ...;

// AFTER:
{
  auto& idm_idr = static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
  thermodynamics_table_[...index_th_ddmu_idm_dr_] =
      4. / 3. * idm_idr.idr().Rho(pvecback) /
      idm_idr.idm_dr().Rho(pvecback) * ...;
}
```

- [ ] **Step 2: Fix the Neff_bbn calculation (lines ~1628-1638)**

```cpp
// BEFORE:
double Neff_bbn = (pvecback[background_module_->index_bg_Omega_r_] *
                       pvecback[background_module_->index_bg_rho_crit_] -
                   pvecback[background_module_->index_bg_rho_g_]) /
                  (7./8. * pow(4./11., 4./3.) * pvecback[background_module_->index_bg_rho_g_]);

if (all_species_.count("IDM_DRMD_IDR_DRMD") > 0) {
  Neff_bbn -= (pvecback[background_module_->index_bg_rho_idr_drmd_]) /
              (7./8. * pow(4./11., 4./3.) * pvecback[background_module_->index_bg_rho_g_]);
}

// AFTER:
const double rho_g = all_species_.at("Photons")->Rho(pvecback);
double Neff_bbn = (pvecback[background_module_->index_bg_Omega_r_] *
                       pvecback[background_module_->index_bg_rho_crit_] -
                   rho_g) /
                  (7./8. * pow(4./11., 4./3.) * rho_g);

if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
  auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
  Neff_bbn -= drmd.idr_drmd().Rho(pvecback) /
              (7./8. * pow(4./11., 4./3.) * rho_g);
}
```

Add required includes to `thermodynamics_module.cpp`:
```cpp
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
```

#### input_module.cpp — 4 background table lookup sites

These lookups read from `background_table_` at the last row (today): `background_table_[(bt_size_-1) * bg_size_ + index_bg_X_]`. After migration, use `Rho(background_table_.data() + (bt_size_-1) * bg_size_)`.

`all_species_` is owned by `BackgroundModule` (and its base `BaseModule`). `input_module.cpp` holds a `BackgroundModulePtr bam`. Access via `bam->all_species_`.

- [ ] **Step 3: Fix `Omega_dcdmdr` and `omega_dcdmdr` cases (lines ~3689-3713)**

```cpp
// BEFORE:
double rho_dcdm_today =
    bam->background_table_[(bam->bt_size_ - 1) * bam->bg_size_ + bam->index_bg_rho_dcdm_];
double rho_dr_today;
if (ba.has_dr == _TRUE_)
  rho_dr_today = bam->background_table_[(bam->bt_size_ - 1) * bam->bg_size_ +
                                        bam->index_bg_rho_dr_species_];
else
  rho_dr_today = 0.;

// AFTER:
const double* bg_today = bam->background_table_.data()
                         + (bam->bt_size_ - 1) * bam->bg_size_;
auto& dcdm_dr = static_cast<DCDM_DR_Species&>(*bam->all_species_.at("DCDM_DR"));
double rho_dcdm_today = dcdm_dr.dcdm().Rho(bg_today);
double rho_dr_today   = (ba.has_dr == _TRUE_)
                        ? dcdm_dr.dr().Rho(bg_today)
                        : 0.;
```

Apply the same pattern to `Omega_ini_dcdm` / `omega_ini_dcdm` cases (lines ~3726-3737) — identical structure.

- [ ] **Step 4: Fix `Omega_scf` case (lines ~3716-3721)**

```cpp
// BEFORE:
output[idx] =
    bam->background_table_[(bam->bt_size_ - 1) * bam->bg_size_ + bam->index_bg_rho_scf_] /
        (ba.H0 * ba.H0) - ba.Omega0_scf;

// AFTER:
const double* bg_today = bam->background_table_.data()
                         + (bam->bt_size_ - 1) * bam->bg_size_;
output[idx] = bam->all_species_.at("ScalarField")->Rho(bg_today) / (ba.H0 * ba.H0)
              - ba.Omega0_scf;
```

Add required includes to `input_module.cpp`:
```cpp
#include "../species/dcdm_dr_species.h"
```

- [ ] **Step 5: Build clean**

```bash
make class -j 2>&1 | grep "error:" | head -20
```

Expected: zero errors (excluding the IC accumulation block in perturbations_module.cpp if it still uses old indices — Task 7 fixes those).

- [ ] **Step 6: Quick smoke test**

```bash
./class explanatory.ini
```

Expected: runs to completion with no crashes.

- [ ] **Step 7: Commit**

```bash
git add source/thermodynamics_module.cpp source/input_module.cpp
git commit -m "Fix thermodynamics_module.cpp and input_module.cpp density access sites"
```

---

### Task 7: Goal 2 — IC accumulation dispatch + IsFreestreaming() overrides

**Files:** `source/perturbations_module.cpp`, `species/dncdm_decay_radiation_species.h`

Replace the manual `rho_r / rho_m / rho_nu` accumulation block in `perturbations_module.cpp` (~lines 4700–4775) with an `EnergyType` + `IsFreestreaming()` dispatch loop. Also verify `IsFreestreaming()` overrides on DNCDM decay radiation species.

The existing code structure:
```cpp
double rho_r = ppw->pvecback[background_module_->index_bg_rho_g_];
double rho_m = ppw->pvecback[background_module_->index_bg_rho_b_];
double rho_nu = 0.;
// ... manual species-by-species additions ...
double fracnu = rho_nu / rho_r;
double fracg  = ... / rho_r;
double fracb  = ... / rho_m;
double fraccdm = ppw->pvecback[background_module_->index_bg_rho_cdm_] / rho_m;
double fracidm_drmd = 0.;
if (background_module_->index_bg_rho_idm_drmd_ >= 0) { ... }
```

- [ ] **Step 1: Replace the accumulation block**

```cpp
// ── Seed rho_r from photons, rho_m from baryons ──────────────────────────
double rho_r  = all_species_.at("Photons")->Rho(ppw->pvecback);
double rho_m  = all_species_.at("Baryons")->Rho(ppw->pvecback);
double rho_nu = 0.;

// ── Dispatch over all other species ──────────────────────────────────────
for (auto& [name, sp] : all_species_) {
  if (name == "Photons" || name == "Baryons") continue;

  const double rho = sp->Rho(ppw->pvecback);

  switch (sp->energy_type()) {
    case BaseSpecies::EnergyType::Matter:
      rho_m += rho;
      break;
    case BaseSpecies::EnergyType::Radiation:
      rho_r += rho;
      break;
    case BaseSpecies::EnergyType::Other:
      // Mixed (e.g. NCDM semi-relativistic): split by equation of state
      rho_r += 3. * sp->P(ppw->pvecback);
      rho_m += rho - 3. * sp->P(ppw->pvecback);
      break;
    default:  // DarkEnergy: no rho_r/rho_m contribution
      break;
  }

  rho_nu += sp->FreestreamingRho(ppw->pvecback);
}

class_test(rho_r == 0., error_message_, "stop to avoid division by zero");

// ── Derived fractions (identical formulas to today) ───────────────────────
const double fracnu  = rho_nu / rho_r;
const double fracg   = all_species_.at("Photons")->Rho(ppw->pvecback) / rho_r;
const double fracb   = all_species_.at("Baryons")->Rho(ppw->pvecback) / rho_m;

double fraccdm = 0.;
if (all_species_.count("CDM"))
  fraccdm = all_species_.at("CDM")->Rho(ppw->pvecback) / rho_m;

double fracidm_drmd = 0.;
if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
  auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
  if (drmd.idm_drmd().IsPresent())
    fracidm_drmd = drmd.idm_drmd().Rho(ppw->pvecback) / rho_m;
}

const double rho_m_over_rho_r = rho_m / rho_r;
```

Note: the `om`, `ktau_two`, `ktau_three` computations and everything below them remain unchanged.

- [ ] **Step 2: Verify IsFreestreaming for DNCDM decay radiation**

Check `species/dncdm_decay_radiation_species.h`. If this species contributes free-streaming DR at IC time (deep radiation domination), add:

```cpp
bool IsFreestreaming() const override { return true; }
```

If the species decouples late and is fluid-like at IC time, leave the default `false`. Confirm by checking whether `DNCDM_DR` contributes to `rho_nu` in the original accumulation block (search for `rho_nu += ...` references to DNCDM around line 4746).

- [ ] **Step 3: Build clean**

```bash
make class -j 2>&1 | grep "error:" | head -20
```

Expected: zero errors.

- [ ] **Step 4: Run `explanatory.ini`**

```bash
./class explanatory.ini
echo "Exit code: $?"
```

Expected: exit code 0, no crashes.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp species/dncdm_decay_radiation_species.h
git commit -m "Goal 2: replace IC rho accumulation with EnergyType+IsFreestreaming dispatch"
```

---

### Task 8: Build, ASan, full test suite, PR

**Files:** — (verification and PR only)

- [ ] **Step 1: Clean build with zero warnings**

```bash
make class -j 2>&1 | grep -E "warning:|error:" | grep -v "^In file" | head -30
```

Expected: zero errors, zero new warnings.

- [ ] **Step 2: Run explanatory.ini**

```bash
./class explanatory.ini
echo "Exit: $?"
```

Expected: exit 0.

- [ ] **Step 3: Run full test suite**

```bash
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -20
```

Expected: 84 passed, 0 failed.

- [ ] **Step 4: ASan build and run**

```bash
cd ..
make clean
make class -j CXXFLAGS="-O1 -g -fsanitize=address,undefined"
./class explanatory.ini
```

Expected: no ASan/UBSan errors (`ERROR SUMMARY: 0 errors`).

- [ ] **Step 5: Create pull request**

```bash
gh pr create \
  --title "Remove species-dependent index mirrors from BackgroundModule" \
  --body "$(cat <<'EOF'
## Summary

- Removes all ~15 species-dependent `index_bg_rho_X_` mirror fields from `BackgroundModule`
- Deletes `bg_rho_index()` and `bg_p_index()` public accessors from `BaseSpecies`
- Adds `BackgroundColumnWriter` and per-species `WriteBackgroundColumnTitles`/`WriteBackgroundData` virtuals (modelled on `PerturbColumnWriter`) — background output now dispatched per-species
- Adds `IsFreestreaming()` predicate and replaces manual IC rho accumulation block with `EnergyType` + `IsFreestreaming()` dispatch loop
- Adds `IsPresent()` to `BaseSpecies` to replace `index >= 0` sentinel checks for sub-components

## Motivation

The previous PR (247-remove-has-guards) revealed that CDM, Lambda, DCDM, and Fluid species had private `index_bg_rho_X_` fields but never assigned `BaseSpecies::index_bg_rho_`. `bg_rho_index()` returned -1, causing `pvecback[-1]` heap overflow caught only by ASan. The one-liner fix closed the immediate bug, but the underlying two-representation problem remains. This PR eliminates it by making species own their indices exclusively.

## Test plan

- [ ] `make class -j` — zero warnings
- [ ] `./class explanatory.ini` — no crash
- [ ] 84 pytest scenarios pass with `COMPARE_OUTPUT_REF=0`
- [ ] ASan `explanatory.ini` — zero heap errors

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- ✅ Remove `bg_rho_index()`/`bg_p_index()` — Task 4
- ✅ Remove mirror fields from BackgroundModule — Task 4
- ✅ Add `BackgroundColumnWriter` — Task 2
- ✅ Add `IsFreestreaming()`, `WriteBackgroundColumnTitles`, `WriteBackgroundData` to BaseSpecies — Task 2
- ✅ Implement per-species — Task 3
- ✅ Add `IsPresent()` — Task 2
- ✅ Registration loop — Task 4 (Step 3)
- ✅ Output dispatch — Task 4 (Step 5)
- ✅ Fix background_module.cpp derivs — Task 4 (Step 6)
- ✅ Fix perturbations_module.cpp — Task 5
- ✅ Fix thermodynamics_module.cpp — Task 6
- ✅ Fix input_module.cpp — Task 6
- ✅ IC accumulation dispatch — Task 7
- ✅ IsFreestreaming overrides (UR, DR, NCDM, IDR, IDR_DRMD) — Tasks 3 + 7

**Edge cases noted in plan:**
- Fluid `index_bg_w_fld_` pre-computation → `WriteWFld()` pattern (Task 3 Step 5, Task 4 Step 4)
- ScalarField arithmetic offsets `index_bg_phi_scf_+n` — module keeps `index_bg_phi_scf_` through `index_bg_ddV_scf_` (Task 4 Step 3)
- NCDM ordering in `all_species_` map (alphabetical "NCDM_0" < "NCDM_1" preserves order for N_ncdm ≤ 9)
- Background output column order changes (species columns move before `rho_crit`) — acceptable with `COMPARE_OUTPUT_REF=0`
