# Perturbations Species Dispatch — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract species-specific logic from `perturbations_module.cpp` into species virtual methods, dispatched via `all_species_` loops, in three testable stages.

**Architecture:** `BaseSpecies` gains four virtual methods (`WriteOutputColumns`, `PrintVariables`, `FillSources`, `ApplyInitialConditions`). Two new context structs (`PerturbSourceContext`, `PerturbIcContext`) bundle shared pre-computed state. A `PerturbColumnWriter` helper guarantees title/data column ordering can never diverge. Module retains `index_tp_*` ownership.

**Tech Stack:** C++17, CMake/Makefile, pytest (`test_class.py`), classy Python wrapper.

---

## File Structure

| File | Change |
|------|--------|
| `species/perturb_source_context.h` | Add `PerturbColumnWriter`, update `PerturbSourceContext`, add `PerturbIcContext` |
| `species/base_species.h` | Add four virtual methods |
| `source/perturbations_module.h` | Add `SetSourceValue` accessor, expose `has_source_*` and `index_tp_*` needed by species |
| `source/perturbations_module.cpp` | Wire dispatch loops in 5 functions across 3 stages |
| `species/photons.h` / `photons.cpp` | Implement all 4 methods |
| `species/baryons.h` / `baryons.cpp` | Implement all 4 methods |
| `species/cdm.h` / `cdm.cpp` | Implement all 4 methods |
| `species/ultra_relativistic.h` / `ultra_relativistic.cpp` | Implement all 4 methods |
| `species/dark_radiation_species.h` / `dark_radiation_species.cpp` | Stage 1–3 (DR output + sources) |
| `species/dcdm.h` / `dcdm.cpp` | Stage 2–3 (DCDM sources + ICs) |
| `species/dcdm_dr_species.h` / `dcdm_dr_species.cpp` | Stage 1 output for composite |
| `species/ncdm_species.h` / `ncdm_species.cpp` | Implement all 4 methods |
| `species/dncdm_species.h` / `dncdm_species.cpp` | Stage 2–3 (mirrors NCDMSpecies) |
| `species/fluid.h` / `fluid.cpp` | Implement all 4 methods |
| `species/scalar_field.h` / `scalar_field.cpp` | Implement all 4 methods |
| `species/idm_dr_idr_species.h` / `idm_dr_idr_species.cpp` | Implement all 4 methods |
| `species/idm_drmd_idr_drmd_species.h` / `idm_drmd_idr_drmd_species.cpp` | Implement all 4 methods |

---

## STAGE 1 — Output columns and print variables

---

### Task 1: Infrastructure — `PerturbColumnWriter`, updated context structs, `BaseSpecies` additions

**Files:**
- Modify: `species/perturb_source_context.h`
- Modify: `species/base_species.h`
- Modify: `source/perturbations_module.h`

- [ ] **Step 1.1: Replace the body of `species/perturb_source_context.h`**

The file currently only contains `PerturbScalarContext`. Add the writer and the new context structs below it. Also add forward declarations at the top.

```cpp
#pragma once
#include <vector>

// Forward declarations
struct perturb_workspace;
struct perturb_vector;
struct precision;
class PerturbationsModule;

// ── existing PerturbScalarContext (unchanged) ────────────────────────────────
struct PerturbScalarContext {
  double k = 0., k2 = 0.;
  double a = 0., a2 = 0., a_prime_over_a = 0.;
  double metric_continuity = 0.;
  double metric_euler = 0.;
  double metric_shear = 0.;
  double metric_ufa_class = 0.;
  double cotKgen = 0., s2_squared = 1.;
  double delta_g = 0., theta_g = 0.;
  double shear_g = 0.;
  double delta_b = 0., theta_b = 0.;
  int idr_nature = 0;
  double R = 0.;
  double cb2 = 0.;
  double delta_p_b_over_rho_b = 0.;
  int gauge = 0;
};

// ── PerturbColumnWriter ──────────────────────────────────────────────────────
/**
 * Thin helper that handles both title-writing and data-writing in one pass.
 * Construct in title mode or one of two data modes; call Add() once per column.
 *
 * The two Add overloads let species use the same writer for:
 *   - WriteOutputColumns: Add(title, tp_index, active)  — value = tk[tp_index]
 *   - PrintVariables    : Add(title, value,    active)  — value passed directly
 */
class PerturbColumnWriter {
public:
  // Title mode (used from perturb_output_titles and perturb_prepare_k_output)
  explicit PerturbColumnWriter(char* titles)
    : titles_(titles) {}

  // Data mode for WriteOutputColumns: value = tk[tp_index]
  PerturbColumnWriter(double* dataptr, const double* tk, int& storeidx)
    : dataptr_(dataptr), tk_(tk), storeidx_(&storeidx) {}

  // Data mode for PrintVariables: value passed directly
  PerturbColumnWriter(double* dataptr, int& storeidx)
    : dataptr_(dataptr), storeidx_(&storeidx) {}

  bool IsTitleMode() const { return titles_ != nullptr; }

  // For WriteOutputColumns: species pass a tp_index; writer looks up tk[tp_index]
  void Add(const char* title, int tp_index, bool active);

  // For PrintVariables: species pass the value directly
  void Add(const char* title, double value, bool active);

private:
  char*         titles_   = nullptr;
  double*       dataptr_  = nullptr;
  const double* tk_       = nullptr;
  int*          storeidx_ = nullptr;
};

// ── PerturbSourceContext ─────────────────────────────────────────────────────
/**
 * Context passed to BaseSpecies::FillSources().
 * Bundles workspace, addressing info, and pre-computed N-body gauge corrections.
 */
struct PerturbSourceContext {
  PerturbationsModule* p_mod   = nullptr;
  perturb_workspace*   ppw     = nullptr;
  int index_md  = 0;
  int index_ic  = 0;
  int index_k   = 0;
  int index_tau = 0;
  double k             = 0.;
  double a_rel         = 0.;
  double a2_rel        = 0.;
  double a_prime_over_a = 0.;
  // Pre-computed N-body gauge corrections (0 when not has_Nbody_gauge_transfers)
  double theta_over_k2 = 0.;
  double theta_shift   = 0.;
};

// ── PerturbIcContext ─────────────────────────────────────────────────────────
/**
 * Context passed to BaseSpecies::ApplyInitialConditions().
 * The module pre-computes all shared IC ratios and seed values; species
 * just write to y[pv->index_pt_*].
 */
struct PerturbIcContext {
  // Density fractions
  double fracnu = 0., fracg = 0., fracb = 0., fraccdm = 0., fracidm_drmd = 0.;
  double rho_m_over_rho_r = 0., om = 0.;
  // (k*tau)^2, (k*tau)^3
  double ktau_two = 0., ktau_three = 0.;
  // Curvature factor: 1 - 3K/k^2
  double s2_squared = 1.;
  // IC seed values pre-computed by module (photon delta/theta set by module formulas)
  double delta_g_ic = 0., theta_g_ic = 0.;
  // Shared relativistic IC (= delta_g_ic for adiabatic, computed by module for isocurvatures)
  double delta_ur = 0., theta_ur = 0., shear_ur = 0., l3_ur = 0., delta_dr = 0.;
  // Synchronous metric IC
  double eta = 0.;
  // Gauge slip (only meaningful when gauge == newtonian, computed after species loop)
  double alpha = 0., alpha_prime = 0.;
  // Kinematics
  double k = 0., tau = 0., a = 0., a_prime_over_a = 0.;
  int index_ic = 0;
  int gauge = 0;
  // Pointers
  perturb_workspace*          ppw   = nullptr;
  const precision*            ppr   = nullptr;
  const PerturbationsModule*  p_mod = nullptr;
};
```

- [ ] **Step 1.2: Add `PerturbColumnWriter` implementation to a new `.cpp` file**

Create `species/perturb_column_writer.cpp`:

```cpp
#include "perturb_source_context.h"
#include "../include/common.h"   // for class_store_columntitle, class_store_double, _TRUE_, _FALSE_

void PerturbColumnWriter::Add(const char* title, int tp_index, bool active) {
  if (titles_) {
    class_store_columntitle(titles_, title, active ? _TRUE_ : _FALSE_);
  } else if (dataptr_ && tk_) {
    class_store_double(dataptr_, tk_[tp_index], active ? _TRUE_ : _FALSE_, *storeidx_);
  }
}

void PerturbColumnWriter::Add(const char* title, double value, bool active) {
  if (titles_) {
    class_store_columntitle(titles_, title, active ? _TRUE_ : _FALSE_);
  } else if (dataptr_) {
    class_store_double(dataptr_, value, active ? _TRUE_ : _FALSE_, *storeidx_);
  }
}
```

- [ ] **Step 1.3: Add new virtual methods to `species/base_species.h`**

After the existing `RhoPlusPShear` declaration (around line 193), add:

```cpp
  // ── Stage 1: Output ──────────────────────────────────────────────────────

  /**
   * Append this species' output columns (delta, theta) to writer.
   * Called from perturb_output_titles (title mode) and perturb_output_data (data mode).
   * Use writer.Add(title, mod.index_tp_XXX_, active) — writer looks up tk[tp_index].
   */
  virtual void WriteOutputColumns(PerturbColumnWriter& /*writer*/,
                                   const PerturbationsModule& /*mod*/,
                                   enum file_format /*fmt*/) const {}

  /**
   * Append this species' time-evolution variables to writer.
   * Called from perturb_prepare_k_output (title mode, tau=0/y=nullptr/ppw=nullptr)
   * and perturb_print_variables_member (data mode).
   * Guard computation behind: if (!writer.IsTitleMode()) { ... compute ... }
   * Use writer.Add(title, computed_value, active).
   */
  virtual void PrintVariables(PerturbColumnWriter& /*writer*/,
                               double /*tau*/,
                               const double* /*y*/,
                               const PerturbationsModule& /*mod*/,
                               const perturb_workspace* /*ppw*/) const {}

  // ── Stage 2: Source filling ───────────────────────────────────────────────

  /**
   * Write per-k, per-tau source values for this species into the source table.
   * Called in perturb_sources_member() after metric/temperature sources are set.
   * Use ctx.p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, ctx.p_mod->index_tp_XXX_,
   *                               ctx.index_tau, ctx.index_k, value).
   */
  virtual void FillSources(int /*index_md*/,
                             const double* /*y*/,
                             const double* /*dy*/,
                             int /*index_tau*/,
                             PerturbSourceContext& /*ctx*/) {}

  // ── Stage 3: Initial conditions ───────────────────────────────────────────

  /**
   * Write synchronous-gauge initial conditions for this species into y[].
   * y[] == ppw->pv->y. Use ctx.ppw->pv->index_pt_* for index lookup.
   * The module pre-computes ctx.delta_g_ic, ctx.delta_ur, etc.; species use them.
   * Guard each IC type: if (ctx.index_ic == ctx.p_mod->index_ic_ad_) { ... }
   * Newtonian gauge transformation is handled by the module after this loop.
   */
  virtual void ApplyInitialConditions(double* /*y*/,
                                       const PerturbIcContext& /*ctx*/) {}
```

Also add the required includes near the top of `base_species.h` (after the existing forward declarations):

```cpp
// PerturbColumnWriter is defined here
#include "perturb_source_context.h"
```

- [ ] **Step 1.4: Add `SetSourceValue` to `perturbations_module.h`**

Find the existing `GetPerturbs()` accessor block and add:

```cpp
  /** Write a value into the source table at (mode, ic, type, tau, k). */
  void SetSourceValue(int index_md, int index_ic, int index_tp,
                      int index_tau, int index_k, double value) {
    sources_[index_md][index_ic * tp_size_[index_md] + index_tp]
            [index_tau * k_size_[index_md] + index_k] = value;
  }
```

- [ ] **Step 1.5: Add `perturb_column_writer.cpp` to the build**

In `Makefile`, add `species/perturb_column_writer.cpp` to the list of compiled source files (same location as other species `.cpp` files). Also update `setup.py` / `pyproject.toml` if they enumerate source files.

- [ ] **Step 1.6: Build and confirm it compiles (no species implement the methods yet)**

```bash
make -j4 class 2>&1 | tail -20
```
Expected: builds cleanly. No species implement the new methods yet so all default no-ops.

- [ ] **Step 1.7: Commit**

```bash
git add species/perturb_source_context.h species/perturb_column_writer.cpp \
        species/base_species.h source/perturbations_module.h
git commit -m "feat: add PerturbColumnWriter, context structs, BaseSpecies dispatch interface"
```

---

### Task 2: Stage 1 — Photons `WriteOutputColumns` + `PrintVariables`

**Files:**
- Modify: `species/photons.h`
- Modify: `species/photons.cpp`

Photons contribute: `d_g`, `t_g` (transfer output) and `delta_g`, `theta_g`, `shear_g`, `pol0_g`, `pol1_g`, `pol2_g` (print variables).

- [ ] **Step 2.1: Add declarations to `species/photons.h`** (after existing PerturbTensorDerivs declaration)

```cpp
  void WriteOutputColumns(PerturbColumnWriter& writer, const PerturbationsModule& mod,
                           enum file_format fmt) const override;
  void PrintVariables(PerturbColumnWriter& writer, double tau, const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
```

- [ ] **Step 2.2: Implement in `species/photons.cpp`**

```cpp
#include "perturbations.h"   // for enum file_format, class_format, camb_format
// (photons.cpp already includes perturbations.h via photons.h)

void PhotonsSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                         const PerturbationsModule& mod,
                                         enum file_format fmt) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_)
      w.Add("d_g", mod.index_tp_delta_g_, true);
    if (ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_g", mod.index_tp_theta_g_, true);
  }
  // camb_format: photons not written separately (only d_g/k^2 via aggregate)
}

void PhotonsSpecies::PrintVariables(PerturbColumnWriter& w, double /*tau*/,
                                     const double* y,
                                     const PerturbationsModule& mod,
                                     const perturb_workspace* ppw) const {
  double delta_g = 0., theta_g = 0., shear_g = 0.;
  double pol0_g = 0., pol1_g = 0., pol2_g = 0.;

  if (!w.IsTitleMode()) {
    const ThermodynamicsModule* thermo = mod.GetThermodynamicsModule().get();
    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];
      if (ppw->approx[ppw->index_ap_tca] == (int)tca_on) {
        shear_g = ppw->tca_shear_g;
        pol0_g  = 2.5 * ppw->tca_shear_g;
        pol1_g  = 7./12. * 6./7. * ppaw_k(ppw) / ppw->pvecthermo[thermo->index_th_dkappa_] * ppw->tca_shear_g;
        pol2_g  = 0.5 * ppw->tca_shear_g;
      } else {
        shear_g = y[ppw->pv->index_pt_shear_g];
        pol0_g  = y[ppw->pv->index_pt_pol0_g];
        pol1_g  = y[ppw->pv->index_pt_pol1_g];
        pol2_g  = y[ppw->pv->index_pt_pol2_g];
      }
    } else {
      delta_g = ppw->rsa_delta_g;
      theta_g = ppw->rsa_theta_g;
    }
    // Synchronous → Newtonian gauge correction for print output
    if (mod.GetPerturbs()->gauge == synchronous) {
      const double H     = ppw->pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a     = ppw->pvecback[mod.GetBackgroundModule()->index_bg_a_];
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      delta_g -= 4. * H * a * alpha;
      theta_g += ppw->scalar_ctx.k2 * alpha;
    }
  }
  w.Add("delta_g", delta_g, true);
  w.Add("theta_g", theta_g, true);
  w.Add("shear_g", shear_g, true);
  w.Add("pol0_g",  pol0_g,  true);
  w.Add("pol1_g",  pol1_g,  true);
  w.Add("pol2_g",  pol2_g,  true);
}
```

> **Note:** `ppaw_k(ppw)` needs `k` — which is not on `ppw` directly. The print function uses `pppaw->k`. Since we pass only `ppw`, either pass `k` explicitly or retrieve it from `ppw->scalar_ctx.k`. Check what field is available; if not, add `double k` to `PrintVariables` signature or use `ppw->scalar_ctx.k`. For now, use `ppw->scalar_ctx.k` (set in `perturb_derivs_member` to the current k).

- [ ] **Step 2.3: Build and confirm**

```bash
make -j4 class 2>&1 | tail -5
```
Expected: clean build.

- [ ] **Step 2.4: Commit**

```bash
git add species/photons.h species/photons.cpp
git commit -m "feat(Stage1): PhotonsSpecies WriteOutputColumns + PrintVariables"
```

---

### Task 3: Stage 1 — Baryons and CDM

**Files:**
- Modify: `species/baryons.h`, `species/baryons.cpp`
- Modify: `species/cdm.h`, `species/cdm.cpp`

- [ ] **Step 3.1: `species/baryons.h`** — add declarations

```cpp
  void WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                           enum file_format fmt) const override;
  void PrintVariables(PerturbColumnWriter& w, double tau, const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
```

- [ ] **Step 3.2: `species/baryons.cpp`** — implement

```cpp
void BaryonsSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                         const PerturbationsModule& mod,
                                         enum file_format fmt) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_) w.Add("d_b", mod.index_tp_delta_b_, true);
    if (ppt->has_velocity_transfers == _TRUE_) w.Add("t_b", mod.index_tp_theta_b_, true);
  }
}

void BaryonsSpecies::PrintVariables(PerturbColumnWriter& w, double /*tau*/,
                                     const double* y, const PerturbationsModule& mod,
                                     const perturb_workspace* ppw) const {
  double delta_b = 0., theta_b = 0.;
  if (!w.IsTitleMode()) {
    delta_b = y[ppw->pv->index_pt_delta_b];
    theta_b = y[ppw->pv->index_pt_theta_b];
    if (mod.GetPerturbs()->gauge == synchronous) {
      const double H = ppw->pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a = ppw->pvecback[mod.GetBackgroundModule()->index_bg_a_];
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      delta_b -= 3. * H * a * alpha;
      theta_b += ppw->scalar_ctx.k2 * alpha;
    }
  }
  w.Add("delta_b", delta_b, true);
  w.Add("theta_b", theta_b, true);
}
```

- [ ] **Step 3.3: `species/cdm.h`** — add declarations

```cpp
  void WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                           enum file_format fmt) const override;
  void PrintVariables(PerturbColumnWriter& w, double tau, const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
```

- [ ] **Step 3.4: `species/cdm.cpp`** — implement

CDM theta is only output in Newtonian gauge (in synchronous gauge it is identically zero, but it IS printed as part of the verbose output with gauge correction applied).

```cpp
void CDMSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                     const PerturbationsModule& mod,
                                     enum file_format fmt) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    const background* pba = mod.GetBackground();
    if (ppt->has_density_transfers == _TRUE_)
      w.Add("d_cdm", mod.index_tp_delta_cdm_, pba->has_cdm == _TRUE_);
    if (ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_cdm", mod.index_tp_theta_cdm_,
            (pba->has_cdm == _TRUE_) && (ppt->gauge != synchronous));
  } else if (fmt == camb_format) {
    w.Add("-T_cdm/k2", mod.index_tp_delta_cdm_, mod.GetBackground()->has_cdm == _TRUE_);
  }
}

void CDMSpecies::PrintVariables(PerturbColumnWriter& w, double /*tau*/,
                                 const double* y, const PerturbationsModule& mod,
                                 const perturb_workspace* ppw) const {
  double delta_cdm = 0., theta_cdm = 0.;
  if (!w.IsTitleMode() && mod.GetBackground()->has_cdm == _TRUE_) {
    delta_cdm = y[ppw->pv->index_pt_delta_cdm];
    if (mod.GetPerturbs()->gauge == synchronous) {
      const double H = ppw->pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a = ppw->pvecback[mod.GetBackgroundModule()->index_bg_a_];
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      delta_cdm -= 3. * H * a * alpha;
      theta_cdm  = ppw->scalar_ctx.k2 * alpha;
    } else {
      theta_cdm = y[ppw->pv->index_pt_theta_cdm];
    }
  }
  const bool has = mod.GetBackground()->has_cdm == _TRUE_;
  w.Add("delta_cdm", delta_cdm, has);
  w.Add("theta_cdm", theta_cdm, has);
}
```

- [ ] **Step 3.5: Build and confirm**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 3.6: Commit**

```bash
git add species/baryons.h species/baryons.cpp species/cdm.h species/cdm.cpp
git commit -m "feat(Stage1): BaryonsSpecies, CDMSpecies WriteOutputColumns + PrintVariables"
```

---

### Task 4: Stage 1 — UltraRelativistic and DR/DCDM families

**Files:**
- Modify: `species/ultra_relativistic.h`, `species/ultra_relativistic.cpp`
- Modify: `species/dark_radiation_species.h`, `species/dark_radiation_species.cpp`
- Modify: `species/dcdm.h`, `species/dcdm.cpp`
- Modify: `species/dcdm_dr_species.h`, `species/dcdm_dr_species.cpp`

- [ ] **Step 4.1: `species/ultra_relativistic.h`** — add declarations

```cpp
  void WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                           enum file_format fmt) const override;
  void PrintVariables(PerturbColumnWriter& w, double tau, const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
```

- [ ] **Step 4.2: `species/ultra_relativistic.cpp`** — implement

```cpp
void UltraRelativisticSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                                    const PerturbationsModule& mod,
                                                    enum file_format fmt) const {
  const bool has = mod.GetBackground()->has_ur == _TRUE_;
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_) w.Add("d_ur", mod.index_tp_delta_ur_, has);
    if (ppt->has_velocity_transfers == _TRUE_) w.Add("t_ur", mod.index_tp_theta_ur_, has);
  } else if (fmt == camb_format) {
    w.Add("-T_ur/k2", mod.index_tp_delta_ur_, has);
  }
}

void UltraRelativisticSpecies::PrintVariables(PerturbColumnWriter& w, double /*tau*/,
                                               const double* y, const PerturbationsModule& mod,
                                               const perturb_workspace* ppw) const {
  double delta_ur = 0., theta_ur = 0., shear_ur = 0.;
  const bool has = mod.GetBackground()->has_ur == _TRUE_;
  if (!w.IsTitleMode() && has) {
    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
      delta_ur = y[ppw->pv->index_pt_delta_ur];
      theta_ur = y[ppw->pv->index_pt_theta_ur];
      shear_ur = y[ppw->pv->index_pt_shear_ur];
    } else {
      delta_ur = ppw->rsa_delta_ur;
      theta_ur = ppw->rsa_theta_ur;
    }
    if (mod.GetPerturbs()->gauge == synchronous) {
      const double H = ppw->pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a = ppw->pvecback[mod.GetBackgroundModule()->index_bg_a_];
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      delta_ur -= 4. * H * a * alpha;
      theta_ur += ppw->scalar_ctx.k2 * alpha;
    }
  }
  w.Add("delta_ur", delta_ur, has);
  w.Add("theta_ur", theta_ur, has);
  w.Add("shear_ur", shear_ur, has);
}
```

- [ ] **Step 4.3: DCDM+DR composite** — add `WriteOutputColumns` + `PrintVariables`

Look at `species/dcdm.h` and `species/dcdm_dr_species.h` to understand the class hierarchy. The `DcdmDrSpecies` composite holds both a `DCDMSpecies` and a `DarkRadiationSpecies`. Add `WriteOutputColumns`/`PrintVariables` to whichever class owns the `index_pt_delta_dcdm`/`index_pt_delta_dr` perturbation variables.

For `DCDMSpecies::WriteOutputColumns`:
```cpp
void DCDMSpecies::WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                                      enum file_format fmt) const {
  const bool has = mod.GetBackground()->has_dcdm == _TRUE_;
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_) w.Add("d_dcdm", mod.index_tp_delta_dcdm_, has);
    if (ppt->has_velocity_transfers == _TRUE_) w.Add("t_dcdm", mod.index_tp_theta_dcdm_, has);
  }
}
```

For `DarkRadiationSpecies::WriteOutputColumns`:
```cpp
void DarkRadiationSpecies::WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                                               enum file_format fmt) const {
  const bool has = mod.GetBackground()->has_dr == _TRUE_;
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_) w.Add("d_dr", mod.index_tp_delta_dr_, has);
    if (ppt->has_velocity_transfers == _TRUE_) w.Add("t_dr", mod.index_tp_theta_dr_, has);
  }
}
```

`PrintVariables` for DCDM:
```cpp
void DCDMSpecies::PrintVariables(PerturbColumnWriter& w, double /*tau*/, const double* y,
                                  const PerturbationsModule& mod, const perturb_workspace* ppw) const {
  double delta_dcdm = 0., theta_dcdm = 0.;
  const bool has = mod.GetBackground()->has_dcdm == _TRUE_;
  if (!w.IsTitleMode() && has) {
    delta_dcdm = y[ppw->pv->index_pt_delta_dcdm];
    theta_dcdm = y[ppw->pv->index_pt_theta_dcdm];
    if (mod.GetPerturbs()->gauge == synchronous) {
      const double H = ppw->pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a = ppw->pvecback[mod.GetBackgroundModule()->index_bg_a_];
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      delta_dcdm += alpha * (-a * mod.GetBackground()->Gamma_dcdm - 3. * H * a);
      theta_dcdm += ppw->scalar_ctx.k2 * alpha;
    }
  }
  w.Add("delta_dcdm", delta_dcdm, has);
  w.Add("theta_dcdm", theta_dcdm, has);
}
```

`PrintVariables` for `DarkRadiationSpecies` outputs per-decay-DR species using loop over `pba->N_decay_dr` (read the existing loop at lines 7665–7691 of `perturb_print_variables_member` and reproduce here).

- [ ] **Step 4.4: Build and confirm**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 4.5: Commit**

```bash
git add species/ultra_relativistic.h species/ultra_relativistic.cpp \
        species/dark_radiation_species.h species/dark_radiation_species.cpp \
        species/dcdm.h species/dcdm.cpp species/dcdm_dr_species.h species/dcdm_dr_species.cpp
git commit -m "feat(Stage1): UR, DR/DCDM WriteOutputColumns + PrintVariables"
```

---

### Task 5: Stage 1 — NCDM, IDM/IDR families, Fluid, ScalarField

**Files:**
- Modify: `species/ncdm_species.h`, `species/ncdm_species.cpp`
- Modify: `species/fluid.h`, `species/fluid.cpp`
- Modify: `species/scalar_field.h`, `species/scalar_field.cpp`
- Modify: `species/idm_dr_idr_species.h`, `species/idm_dr_idr_species.cpp`
- Modify: `species/idm_drmd_idr_drmd_species.h`, `species/idm_drmd_idr_drmd_species.cpp`

- [ ] **Step 5.1: `NCDMSpecies::WriteOutputColumns`**

NCDM needs `ncdm_id()` to look up the correct `index_tp_delta_ncdm1_ + n` and `index_tp_theta_ncdm1_ + n`.

```cpp
void NCDMSpecies::WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                                      enum file_format fmt) const {
  const int n = ncdm_id();
  const bool has = mod.GetBackground()->has_ncdm == _TRUE_;
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    char buf[64];
    if (ppt->has_density_transfers == _TRUE_) {
      snprintf(buf, sizeof(buf), "d_ncdm[%d]", n);
      w.Add(buf, mod.index_tp_delta_ncdm1_ + n, has);
    }
    if (ppt->has_velocity_transfers == _TRUE_) {
      snprintf(buf, sizeof(buf), "t_ncdm[%d]", n);
      w.Add(buf, mod.index_tp_theta_ncdm1_ + n, has);
    }
  } else if (fmt == camb_format) {
    w.Add("-T_ncdm/k2", mod.index_tp_delta_ncdm1_ + n, has);
  }
}
```

- [ ] **Step 5.2: `NCDMSpecies::PrintVariables`**

The print function stores delta_ncdm, theta_ncdm, shear_ncdm, delta_p_over_delta_rho_ncdm. These are pre-computed in the module into `ppw->delta_ncdm[n]` and `ppw->theta_ncdm[n]`. However in `perturb_print_variables_member` they are extracted from `ppw->pv->y` directly using `Delta()`, `Theta()`, `RhoPlusPShear()`, `DeltaP()`. Use those virtual methods:

```cpp
void NCDMSpecies::PrintVariables(PerturbColumnWriter& w, double /*tau*/, const double* y,
                                   const PerturbationsModule& mod,
                                   const perturb_workspace* ppw) const {
  const int n = ncdm_id();
  const bool has = mod.GetBackground()->has_ncdm == _TRUE_;
  double delta = 0., theta = 0., shear = 0., cs2 = 0.;
  if (!w.IsTitleMode() && has) {
    delta = Delta(ppw->pv, y, ppw->pvecback, ppw);
    theta = Theta(ppw->pv, y, ppw->pvecback, ppw);
    const double rho   = Rho(ppw->pvecback);
    const double p     = P(ppw->pvecback);
    shear = RhoPlusPShear(ppw->pv, y, ppw->pvecback, ppw) / (rho + p);
    const double delta_p = DeltaP(ppw->pv, y, ppw->pvecback, ppw);
    const double delta_rho = rho * delta;
    cs2 = (std::abs(delta_rho) > 1e-300) ? delta_p / delta_rho : 0.;
    // Synchronous → Newtonian gauge correction
    if (mod.GetPerturbs()->gauge == synchronous) {
      const double H = ppw->pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a = ppw->pvecback[mod.GetBackgroundModule()->index_bg_a_];
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      const double w_ncdm = p / rho;
      delta -= 3. * a * H * (1. + w_ncdm) * alpha;
      theta += ppw->scalar_ctx.k2 * alpha;
    }
  }
  char buf[64];
  snprintf(buf, sizeof(buf), "delta_ncdm[%d]", n); w.Add(buf, delta, has);
  snprintf(buf, sizeof(buf), "theta_ncdm[%d]", n); w.Add(buf, theta, has);
  snprintf(buf, sizeof(buf), "shear_ncdm[%d]", n); w.Add(buf, shear, has);
  snprintf(buf, sizeof(buf), "cs2_ncdm[%d]",   n); w.Add(buf, cs2,   has);
}
```

- [ ] **Step 5.3: IDM/IDR families**

`IDM_DR_IDR_Species::WriteOutputColumns`:
```cpp
void IDM_DR_IDR_Species::WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                                              enum file_format fmt) const {
  if (fmt != class_format) return;
  const perturbs* ppt = mod.GetPerturbs();
  const background* pba = mod.GetBackground();
  if (ppt->has_density_transfers == _TRUE_) {
    w.Add("d_idm_dr", mod.index_tp_delta_idm_dr_, pba->has_idm_dr == _TRUE_);
    w.Add("d_idr",    mod.index_tp_delta_idr_,    pba->has_idr    == _TRUE_);
  }
  if (ppt->has_velocity_transfers == _TRUE_) {
    w.Add("t_idm_dr", mod.index_tp_theta_idm_dr_, pba->has_idm_dr == _TRUE_);
    w.Add("t_idr",    mod.index_tp_theta_idr_,    pba->has_idr    == _TRUE_);
  }
}
```

`IDM_DR_IDR_Species::PrintVariables`: extract from `perturb_print_variables_member` lines 7572–7569 (IDR/IDM_DR section).

`IDM_DRMD_IDR_DRMD_Species` follows the same pattern but for `has_idr_drmd` / `has_idm_drmd`.

- [ ] **Step 5.4: `FluidSpecies::WriteOutputColumns`**

```cpp
void FluidSpecies::WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                                       enum file_format fmt) const {
  const bool has = mod.GetBackground()->has_fld == _TRUE_;
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_) w.Add("d_fld", mod.index_tp_delta_fld_, has);
    if (ppt->has_velocity_transfers == _TRUE_) w.Add("t_fld", mod.index_tp_theta_fld_, has);
  }
}
```

`FluidSpecies::PrintVariables`: read `ppw->delta_rho_fld`, `ppw->rho_plus_p_theta_fld`, `ppw->delta_p_fld` (as in perturb_print_variables_member lines 7896–7899).

- [ ] **Step 5.5: `ScalarFieldSpecies::WriteOutputColumns`**

```cpp
void ScalarFieldSpecies::WriteOutputColumns(PerturbColumnWriter& w, const PerturbationsModule& mod,
                                              enum file_format fmt) const {
  const bool has = mod.GetBackground()->has_scf == _TRUE_;
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->has_density_transfers == _TRUE_) w.Add("d_scf", mod.index_tp_delta_scf_, has);
    if (ppt->has_velocity_transfers == _TRUE_) w.Add("t__scf", mod.index_tp_theta_scf_, has);
  }
}
```

> Note the double underscore typo `"t__scf"` in the existing output_titles at line 377 — preserve it to maintain backward compatibility with users who depend on the column name.

`ScalarFieldSpecies::PrintVariables`: compute delta_scf and theta_scf as in perturb_print_variables_member lines 7675–7694.

- [ ] **Step 5.6: Build and confirm**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 5.7: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp \
        species/fluid.h species/fluid.cpp \
        species/scalar_field.h species/scalar_field.cpp \
        species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp
git commit -m "feat(Stage1): NCDM, IDM/IDR, Fluid, ScalarField WriteOutputColumns + PrintVariables"
```

---

### Task 6: Stage 1 — Wire module dispatch loops

**Files:**
- Modify: `source/perturbations_module.cpp`

This task replaces the per-species hardcoded blocks in four functions with dispatch loops.

- [ ] **Step 6.1: Replace species blocks in `perturb_output_titles`**

In `perturb_output_titles` (around line 325–378), replace the species-specific `class_store_columntitle` calls with:

```cpp
// After "k (h/Mpc)" title:
{
  PerturbColumnWriter w(titles);
  for (auto& [name, sp] : all_species_)
    sp->WriteOutputColumns(w, *this, output_format);
}
// Keep the aggregate titles (d_tot, phi, psi, h, eta, etc.) exactly as before:
if (ppt->has_density_transfers == _TRUE_) {
  class_store_columntitle(titles, "d_tot", _TRUE_);
  class_store_columntitle(titles, "phi", has_source_phi_);
  // ... etc unchanged ...
}
```

Also remove the per-species blocks (d_g, d_b, d_cdm, d_idm_dr, d_idm_drmd, d_fld, d_ur, d_idr, d_idr_drmd, the ncdm loop, d_dcdm, d_dr, d_scf, t_g, t_b, t_cdm, t_idm_dr, t_idm_drmd, t_fld, t_ur, t_idr, t_idr_drmd, the ncdm loop, t_dcdm, t_dr, t__scf).

- [ ] **Step 6.2: Replace species blocks in `perturb_output_data`**

In `perturb_output_data` (around line 233–285), replace per-species `class_store_double` calls with:

```cpp
if (ppt->has_density_transfers == _TRUE_) {
  {
    PerturbColumnWriter w(dataptr, tk, storeidx);
    for (auto& [name, sp] : all_species_)
      sp->WriteOutputColumns(w, *this, output_format);
  }
  // Keep aggregate stores (delta_tot, phi, psi, h, eta, etc.):
  class_store_double(dataptr, tk[index_tp_delta_tot_], has_source_delta_tot_, storeidx);
  // ... etc ...
}
if (ppt->has_velocity_transfers == _TRUE_) {
  {
    PerturbColumnWriter w(dataptr, tk, storeidx);
    for (auto& [name, sp] : all_species_)
      sp->WriteOutputColumns(w, *this, output_format);
  }
  // Keep theta_tot, theta_m, theta_cb:
  class_store_double(dataptr, tk[index_tp_theta_tot_], has_source_theta_tot_, storeidx);
  // ...
}
```

> **Important:** `WriteOutputColumns` must be called in the density block AND the velocity block — OR the species method guards internally on `has_density_transfers`/`has_velocity_transfers`. Choose one approach; the current design has species guard internally.  In that case the loop is called once before the aggregate stores:

```cpp
{
  PerturbColumnWriter w(dataptr, tk, storeidx);
  for (auto& [name, sp] : all_species_)
    sp->WriteOutputColumns(w, *this, output_format);
}
```

Check the existing column order: photons first, then baryons, then CDM, etc. The `all_species_` map is alphabetical. Verify the resulting column order still matches the titles in `perturb_output_titles` (they must be identical since both use the same `WriteOutputColumns`).

- [ ] **Step 6.3: Replace species blocks in `perturb_prepare_k_output`**

In `perturb_prepare_k_output` (around lines 2720–2807), replace per-species `class_store_columntitle` calls with:

```cpp
// Keep fixed columns: "tau [Mpc]", "a"
class_store_columntitle(scalar_titles_, "tau [Mpc]", _TRUE_);
class_store_columntitle(scalar_titles_, "a", _TRUE_);
// Dispatch to species (title mode, no y/ppw needed)
{
  PerturbColumnWriter w(scalar_titles_);
  for (auto& [name, sp] : all_species_)
    sp->PrintVariables(w, 0., nullptr, *this, nullptr);
}
// Keep fixed columns: psi, phi, perturbed recombination, etc.
class_store_columntitle(scalar_titles_, "psi", _TRUE_);
class_store_columntitle(scalar_titles_, "phi", _TRUE_);
// etc.
```

Remove the old per-species `class_store_columntitle` blocks.

- [ ] **Step 6.4: Replace species blocks in `perturb_print_variables_member`**

In `perturb_print_variables_member` (around lines 7821–7899), keep the pre-computation of `delta_g`, `theta_g`, `shear_g`, `pol0_g`, etc. (RSA/TCA/gauge corrections up to line 7805), then replace the storage loop with:

```cpp
class_store_double(dataptr, tau, _TRUE_, storeidx);
class_store_double(dataptr, pvecback[background_module_->index_bg_a_], _TRUE_, storeidx);
// Species dispatch for all per-species variables:
{
  PerturbColumnWriter w(dataptr, storeidx);
  for (auto& [name, sp] : all_species_)
    sp->PrintVariables(w, tau, y, *this, ppw);
}
// Keep fixed variables: psi, phi, perturbed recombination
class_store_double(dataptr, psi, _TRUE_, storeidx);
class_store_double(dataptr, phi, _TRUE_, storeidx);
class_store_double(dataptr, delta_temp, ppt->has_perturbed_recombination, storeidx);
class_store_double(dataptr, delta_chi,  ppt->has_perturbed_recombination, storeidx);
```

The pre-computation of `psi`, `phi`, `alpha`, `delta_temp`, `delta_chi` stays in the module. Species compute their own values from `ppw` and `y`.

- [ ] **Step 6.5: Build**

```bash
make -j4 class 2>&1 | tail -20
```

Fix any compile errors.

- [ ] **Step 6.6: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "feat(Stage1): wire WriteOutputColumns + PrintVariables dispatch loops in module"
```

---

### Task 7: Stage 1 — Test

- [ ] **Step 7.1: Build the Python wrapper**

```bash
pip install . 2>&1 | tail -10
```

- [ ] **Step 7.2: Run a baseline with `k_output_values`**

```python
# test_stage1.py
from classy import Class
cosmo = Class()
cosmo.set({'output': 'tCl,mPk', 'k_output_values': '0.01,0.1', 'P_k_max_1/Mpc': 2})
cosmo.compute()
cosmo.struct_cleanup()
```

Run: `python test_stage1.py`
Expected: no crash. Produces `output/` files with perturbation data.

- [ ] **Step 7.3: Compare column names with a reference run**

Build the reference CLASS (or use classyref if available):
```bash
# With reference binary in PATH:
class explanatory.ini  # produces reference output
```
Compare column headers: `head -1 output/xxx_scalars_k01_....dat` vs reference.
Expected: identical column names, same ordering.

- [ ] **Step 7.4: Run `TEST_LEVEL=0` pytest**

```bash
TEST_LEVEL=0 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
```
Expected: all tests pass (Stage 1 does not change numerical output — it only changes which code path produces the output data).

- [ ] **Step 7.5: If tests fail, diagnose**

Column count mismatch: a title is added by species AND kept in the module's aggregate block — find and remove duplicate.

Wrong values: compare data from old `class_store_double` path with new `WriteOutputColumns` path using print statements; verify `index_tp_*` values and `has_source_*` flags match.

- [ ] **Step 7.6: Commit test script**

```bash
git add test_stage1.py
git commit -m "test(Stage1): column output verification script"
```

---

## STAGE 2 — Source filling

---

### Task 8: Stage 2 — Photons and Baryons `FillSources`

**Files:**
- Modify: `species/photons.h`, `species/photons.cpp`
- Modify: `species/baryons.h`, `species/baryons.cpp`

First: update `PerturbSourceContext` to include `theta_over_k2` and `theta_shift` (already in spec — confirm they are in the struct from Task 1).

- [ ] **Step 8.1: Add declarations to `species/photons.h`**

```cpp
  void FillSources(int index_md, const double* y, const double* dy,
                   int index_tau, PerturbSourceContext& ctx) override;
```

- [ ] **Step 8.2: Implement `PhotonsSpecies::FillSources`**

Photons contribute `delta_g` and `theta_g` sources. Both have RSA-aware paths.

```cpp
void PhotonsSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                                   int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  const perturb_workspace* ppw = ctx.ppw;
  PerturbationsModule* mod = ctx.p_mod;

  if (mod->has_source_delta_g_ == _TRUE_) {
    double delta_g = (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on)
                   ? ppw->rsa_delta_g
                   : y[ppw->pv->index_pt_delta_g];
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_g_,
                        ctx.index_tau, ctx.index_k,
                        delta_g + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);
  }

  if (mod->has_source_theta_g_ == _TRUE_) {
    double theta_g = (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on)
                   ? ppw->rsa_theta_g
                   : y[ppw->pv->index_pt_theta_g];
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_g_,
                        ctx.index_tau, ctx.index_k,
                        theta_g + ctx.theta_shift);
  }
}
```

- [ ] **Step 8.3: Add declaration to `species/baryons.h`**

```cpp
  void FillSources(int index_md, const double* y, const double* dy,
                   int index_tau, PerturbSourceContext& ctx) override;
```

- [ ] **Step 8.4: Implement `BaryonsSpecies::FillSources`**

```cpp
void BaryonsSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                                   int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;

  if (mod->has_source_delta_b_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_b_,
                        ctx.index_tau, ctx.index_k,
                        y[ctx.ppw->pv->index_pt_delta_b]
                        + 3. * ctx.a_prime_over_a * ctx.theta_over_k2);

  if (mod->has_source_theta_b_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_b_,
                        ctx.index_tau, ctx.index_k,
                        y[ctx.ppw->pv->index_pt_theta_b] + ctx.theta_shift);
}
```

- [ ] **Step 8.5: Build**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 8.6: Commit**

```bash
git add species/photons.h species/photons.cpp species/baryons.h species/baryons.cpp
git commit -m "feat(Stage2): Photons, Baryons FillSources"
```

---

### Task 9: Stage 2 — CDM, UR, DR/DCDM `FillSources`

**Files:**
- Modify: `species/cdm.h`, `species/cdm.cpp`
- Modify: `species/ultra_relativistic.h`, `species/ultra_relativistic.cpp`
- Modify: `species/dcdm.h`, `species/dcdm.cpp`
- Modify: `species/dark_radiation_species.h`, `species/dark_radiation_species.cpp`

- [ ] **Step 9.1: `CDMSpecies::FillSources`**

CDM has `has_source_delta_cdm_` and `has_source_theta_cdm_`.

```cpp
void CDMSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                               int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;

  if (mod->has_source_delta_cdm_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_cdm_,
                        ctx.index_tau, ctx.index_k,
                        y[ctx.ppw->pv->index_pt_delta_cdm]
                        + 3. * ctx.a_prime_over_a * ctx.theta_over_k2);

  if (mod->has_source_theta_cdm_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_cdm_,
                        ctx.index_tau, ctx.index_k,
                        y[ctx.ppw->pv->index_pt_theta_cdm] + ctx.theta_shift);
}
```

- [ ] **Step 9.2: `UltraRelativisticSpecies::FillSources`**

UR has RSA path for both delta and theta.

```cpp
void UltraRelativisticSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                                             int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  const perturb_workspace* ppw = ctx.ppw;
  PerturbationsModule* mod = ctx.p_mod;

  if (mod->has_source_delta_ur_ == _TRUE_) {
    double delta_ur = (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off)
                    ? y[ppw->pv->index_pt_delta_ur] : ppw->rsa_delta_ur;
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_ur_,
                        ctx.index_tau, ctx.index_k,
                        delta_ur + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);
  }

  if (mod->has_source_theta_ur_ == _TRUE_) {
    double theta_ur = (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off)
                    ? y[ppw->pv->index_pt_theta_ur] : ppw->rsa_theta_ur;
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_ur_,
                        ctx.index_tau, ctx.index_k, theta_ur + ctx.theta_shift);
  }
}
```

- [ ] **Step 9.3: `DCDMSpecies::FillSources`**

DCDM uses a modified N-body correction: `(3*a'/a + a*Γ_dcdm)*theta_over_k2`.

```cpp
void DCDMSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                               int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;
  const background* pba = mod->GetBackground();

  if (mod->has_source_delta_dcdm_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_dcdm_,
                        ctx.index_tau, ctx.index_k,
                        y[ctx.ppw->pv->index_pt_delta_dcdm]
                        + (3. * ctx.a_prime_over_a + ctx.a_rel * pba->Gamma_dcdm)
                          * ctx.theta_over_k2);

  if (mod->has_source_theta_dcdm_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_dcdm_,
                        ctx.index_tau, ctx.index_k,
                        y[ctx.ppw->pv->index_pt_theta_dcdm] + ctx.theta_shift);
}
```

- [ ] **Step 9.4: `DarkRadiationSpecies::FillSources`**

DR uses `ppw->pv->index_pt_F0_dr_sum` and `F0_dr_sum + 1` for delta and theta. The normalization uses `r_dr = (a²/H0)² * rho_dr`.

```cpp
void DarkRadiationSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                                        int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;
  const background* pba = mod->GetBackground();
  const perturb_workspace* ppw = ctx.ppw;
  const double a2_over_H0_sq = (ctx.a2_rel / pba->H0) * (ctx.a2_rel / pba->H0);
  const double r_dr = a2_over_H0_sq * ppw->pvecback[mod->GetBackgroundModule()->index_bg_rho_dr_];

  if (mod->has_source_delta_dr_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_dr_,
                        ctx.index_tau, ctx.index_k,
                        y[ppw->pv->index_pt_F0_dr_sum] / r_dr
                        + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);

  if (mod->has_source_theta_dr_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_dr_,
                        ctx.index_tau, ctx.index_k,
                        3./4. * ctx.k * y[ppw->pv->index_pt_F0_dr_sum + 1] / r_dr
                        + ctx.theta_shift);
}
```

- [ ] **Step 9.5: Build**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 9.6: Commit**

```bash
git add species/cdm.h species/cdm.cpp \
        species/ultra_relativistic.h species/ultra_relativistic.cpp \
        species/dcdm.h species/dcdm.cpp \
        species/dark_radiation_species.h species/dark_radiation_species.cpp
git commit -m "feat(Stage2): CDM, UR, DR/DCDM FillSources"
```

---

### Task 10: Stage 2 — NCDM, IDM/IDR families, Fluid, ScalarField `FillSources`

**Files:**
- Modify: `species/ncdm_species.h`, `species/ncdm_species.cpp`
- Modify: `species/fluid.h`, `species/fluid.cpp`
- Modify: `species/scalar_field.h`, `species/scalar_field.cpp`
- Modify: `species/idm_dr_idr_species.h`, `species/idm_dr_idr_species.cpp`
- Modify: `species/idm_drmd_idr_drmd_species.h`, `species/idm_drmd_idr_drmd_species.cpp`

- [ ] **Step 10.1: `NCDMSpecies::FillSources`**

NCDM delta and theta use pre-computed `ppw->delta_ncdm[n]` and `ppw->theta_ncdm[n]`. The N-body correction uses `p/rho` ratio.

```cpp
void NCDMSpecies::FillSources(int index_md, const double* /*y*/, const double* /*dy*/,
                               int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;
  if (mod->has_source_delta_ncdm_ != _TRUE_) return;

  const int n = ncdm_id();
  const perturb_workspace* ppw = ctx.ppw;
  const double rho = ppw->pvecback[bg_rho_index()];
  const double p   = ppw->pvecback[bg_p_index()];

  mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_ncdm1_ + n,
                      ctx.index_tau, ctx.index_k,
                      ppw->delta_ncdm[n]
                      + 3. * ctx.a_prime_over_a * (1. + p/rho) * ctx.theta_over_k2);

  if (mod->has_source_theta_ncdm_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_ncdm1_ + n,
                        ctx.index_tau, ctx.index_k,
                        ppw->theta_ncdm[n] + ctx.theta_shift);
}
```

- [ ] **Step 10.2: `FluidSpecies::FillSources`**

```cpp
void FluidSpecies::FillSources(int index_md, const double* /*y*/, const double* /*dy*/,
                                int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;
  const perturb_workspace* ppw = ctx.ppw;

  if (mod->has_source_delta_fld_ == _TRUE_) {
    const double rho_fld = ppw->pvecback[mod->GetBackgroundModule()->index_bg_rho_fld_];
    const double w_fld   = ppw->pvecback[mod->GetBackgroundModule()->index_bg_w_fld_];
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_fld_,
                        ctx.index_tau, ctx.index_k,
                        ppw->delta_rho_fld / rho_fld
                        + 3. * ctx.a_prime_over_a * (1. + w_fld) * ctx.theta_over_k2);
  }

  if (mod->has_source_theta_fld_ == _TRUE_) {
    double w_fld, dw_over_da, integral;
    // Need background_w_fld to get w_fld at current a:
    mod->GetBackgroundModule()->background_w_fld(ctx.a_rel * mod->GetBackground()->a_today,
                                                   &w_fld, &dw_over_da, &integral);
    const double rho_fld = ppw->pvecback[mod->GetBackgroundModule()->index_bg_rho_fld_];
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_fld_,
                        ctx.index_tau, ctx.index_k,
                        ppw->rho_plus_p_theta_fld / ((1. + w_fld) * rho_fld)
                        + ctx.theta_shift);
  }
}
```

- [ ] **Step 10.3: `ScalarFieldSpecies::FillSources`**

SCF delta and theta are computed from `y[index_pt_phi_scf]` and `y[index_pt_phi_prime_scf]`. The expressions differ by gauge. Read lines 7162–7178 and 7305–7312 from `perturb_sources_member` for the exact formulas.

```cpp
void ScalarFieldSpecies::FillSources(int index_md, const double* y, const double* /*dy*/,
                                      int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;
  const perturb_workspace* ppw = ctx.ppw;
  const BackgroundModule* bgm  = mod->GetBackgroundModule().get();
  const perturbs* ppt = mod->GetPerturbs();

  if (mod->has_source_delta_scf_ == _TRUE_) {
    const double phi_prime = ppw->pvecback[bgm->index_bg_phi_prime_scf_];
    const double dV        = ppw->pvecback[bgm->index_bg_dV_scf_];
    double delta_rho_scf;
    if (ppt->gauge == synchronous) {
      delta_rho_scf = 1./3. * (1./ctx.a2_rel * phi_prime * y[ppw->pv->index_pt_phi_prime_scf]
                                + dV * y[ppw->pv->index_pt_phi_scf])
                      + 3. * ctx.a_prime_over_a
                        * (1. + ppw->pvecback[bgm->index_bg_p_scf_]
                               / ppw->pvecback[bgm->index_bg_rho_scf_])
                        * ctx.theta_over_k2;
    } else {
      delta_rho_scf = 1./3. * (1./ctx.a2_rel * phi_prime * y[ppw->pv->index_pt_phi_prime_scf]
                                + dV * y[ppw->pv->index_pt_phi_scf]
                                - 1./ctx.a2_rel * phi_prime * phi_prime
                                  * ppw->pvecmetric[ppw->index_mt_psi])
                      + 3. * ctx.a_prime_over_a
                        * (1. + ppw->pvecback[bgm->index_bg_p_scf_]
                               / ppw->pvecback[bgm->index_bg_rho_scf_])
                        * ctx.theta_over_k2;
    }
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_scf_,
                        ctx.index_tau, ctx.index_k,
                        delta_rho_scf / ppw->pvecback[bgm->index_bg_rho_scf_]);
  }

  if (mod->has_source_theta_scf_ == _TRUE_) {
    const double phi_prime = ppw->pvecback[bgm->index_bg_phi_prime_scf_];
    const double rho_plus_p_theta_scf =
      1./3. * ctx.k * ctx.k / ctx.a2_rel * phi_prime * y[ppw->pv->index_pt_phi_scf];
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_scf_,
                        ctx.index_tau, ctx.index_k,
                        rho_plus_p_theta_scf
                        / (ppw->pvecback[bgm->index_bg_rho_scf_]
                           + ppw->pvecback[bgm->index_bg_p_scf_])
                        + ctx.theta_shift);
  }
}
```

- [ ] **Step 10.4: IDM/IDR families**

`IDM_DR_IDR_Species::FillSources` covers IDR (with RSA path) and IDM_DR:

```cpp
void IDM_DR_IDR_Species::FillSources(int index_md, const double* y, const double* /*dy*/,
                                      int /*index_tau*/, PerturbSourceContext& ctx) {
  if (index_md != ctx.p_mod->index_md_scalars_) return;
  PerturbationsModule* mod = ctx.p_mod;
  const perturb_workspace* ppw = ctx.ppw;

  // IDR
  if (mod->has_source_delta_idr_ == _TRUE_) {
    double delta_idr = (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_off)
                     ? y[ppw->pv->index_pt_delta_idr] : ppw->rsa_delta_idr;
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_idr_,
                        ctx.index_tau, ctx.index_k,
                        delta_idr + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);
  }
  if (mod->has_source_theta_idr_ == _TRUE_) {
    double theta_idr = (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_off)
                     ? y[ppw->pv->index_pt_theta_idr] : ppw->rsa_theta_idr;
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_idr_,
                        ctx.index_tau, ctx.index_k, theta_idr + ctx.theta_shift);
  }
  // IDM_DR
  if (mod->has_source_delta_idm_dr_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_delta_idm_dr_,
                        ctx.index_tau, ctx.index_k,
                        y[ppw->pv->index_pt_delta_idm_dr]
                        + 3. * ctx.a_prime_over_a * ctx.theta_over_k2);
  if (mod->has_source_theta_idm_dr_ == _TRUE_)
    mod->SetSourceValue(ctx.index_md, ctx.index_ic, mod->index_tp_theta_idm_dr_,
                        ctx.index_tau, ctx.index_k,
                        y[ppw->pv->index_pt_theta_idm_dr] + ctx.theta_shift);
}
```

`IDM_DRMD_IDR_DRMD_Species::FillSources` follows the same pattern for `index_pt_delta_idr_drmd`, `index_pt_theta_idr_drmd`, `index_pt_delta_idm_drmd`, `index_pt_theta_idm_drmd`. IDR-DRMD uses a radiation N-body correction (factor 4); IDM-DRMD uses matter correction (factor 3).

> **Verify:** Check whether `has_source_theta_idr_drmd_` is set in `perturb_indices_of_perturbs` and whether the source was stored in the original `perturb_sources_member`. If missing, that is a pre-existing gap — add a `// TODO: verify` comment and do not fix it here.

- [ ] **Step 10.5: Build**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 10.6: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp \
        species/fluid.h species/fluid.cpp \
        species/scalar_field.h species/scalar_field.cpp \
        species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp
git commit -m "feat(Stage2): NCDM, IDM/IDR, Fluid, ScalarField FillSources"
```

---

### Task 11: Stage 2 — Wire `FillSources` in module

**Files:**
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 11.1: Build `PerturbSourceContext` before the source-filling loop**

In `perturb_sources_member`, after computing `theta_over_k2` and `theta_shift` (around line 7099–7107), construct the context:

```cpp
PerturbSourceContext src_ctx;
src_ctx.p_mod         = this;
src_ctx.ppw           = ppw;
src_ctx.index_md      = index_md;
src_ctx.index_ic      = index_ic;
src_ctx.index_k       = index_k;
src_ctx.index_tau     = index_tau;
src_ctx.k             = k;
src_ctx.a_rel         = a_rel;
src_ctx.a2_rel        = a2_rel;
src_ctx.a_prime_over_a = a_prime_over_a;
src_ctx.theta_over_k2 = theta_over_k2;
src_ctx.theta_shift   = theta_shift;
```

- [ ] **Step 11.2: Replace per-species source blocks with dispatch loop**

After the metric/temperature/aggregate sources (delta_tot, phi, psi, h, eta, delta_m, delta_cb, H_T_Nb_prime), add:

```cpp
// Species-specific sources (delta_g, theta_g, delta_b, theta_b, delta_cdm, ...)
for (auto& [name, sp] : all_species_)
  sp->FillSources(index_md, y, dy, index_tau, src_ctx);
```

Then delete all the old per-species blocks from `/* delta_g */` to `/* theta_ncdm1 */` (lines ~7125–7357 in the original).

- [ ] **Step 11.3: Build**

```bash
make -j4 class 2>&1 | tail -20
```

Fix any compile errors (missing includes, access to private members, etc.).

- [ ] **Step 11.4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "feat(Stage2): wire FillSources dispatch loop in perturb_sources_member"
```

---

### Task 12: Stage 2 — Test

- [ ] **Step 12.1: Build**

```bash
make -j4 class && pip install . 2>&1 | tail -10
```

- [ ] **Step 12.2: Run `TEST_LEVEL=0`**

```bash
TEST_LEVEL=0 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -40
```
Expected: all tests pass with differences below numerical tolerance.

- [ ] **Step 12.3: If tests fail — diagnostic procedure**

Identify which test scenario and which output differs:
```bash
TEST_LEVEL=0 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py -k "failing_scenario" -s 2>&1
```
Likely causes:
- Wrong N-body gauge correction: verify `theta_over_k2`/`theta_shift` are passed correctly in context.
- Wrong RSA branch: verify `ppw->approx[ppw->index_ap_rsa]` check matches original code.
- Missing source: verify `has_source_XXX_` flag is read from module and the species method is actually called.
- `SetSourceValue` off-by-one: verify `index_tp_* ` values match what was used in `_set_source_` macro.

- [ ] **Step 12.4: Commit fixes and re-run until passing**

---

## STAGE 3 — Initial conditions

---

### Task 13: Stage 3 — All species `ApplyInitialConditions`

**Files:**
- Modify: `species/photons.h`, `species/photons.cpp`
- Modify: `species/baryons.h`, `species/baryons.cpp`
- Modify: `species/cdm.h`, `species/cdm.cpp`
- Modify: `species/ultra_relativistic.h`, `species/ultra_relativistic.cpp`
- Modify: `species/idm_dr_idr_species.h`, `species/idm_dr_idr_species.cpp`
- Modify: `species/idm_drmd_idr_drmd_species.h`, `species/idm_drmd_idr_drmd_species.cpp`
- Modify: `species/dcdm.h`, `species/dcdm.cpp`
- Modify: `species/dark_radiation_species.h`, `species/dark_radiation_species.cpp`
- Modify: `species/ncdm_species.h`, `species/ncdm_species.cpp`
- Modify: `species/fluid.h`, `species/fluid.cpp`
- Modify: `species/scalar_field.h`, `species/scalar_field.cpp`

The module pre-computes all seed values (`ctx.delta_g_ic`, `ctx.theta_g_ic`, `ctx.delta_ur`, etc.) before calling species. Species write to `y[pv->index_pt_*]` using `ctx.ppw->pv`. Each species guards on `ctx.index_ic` via `ctx.p_mod->index_ic_ad_` etc.

- [ ] **Step 13.1: `PhotonsSpecies::ApplyInitialConditions`**

```cpp
void PhotonsSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;

  if (ctx.index_ic == mod->index_ic_ad_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_bi_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_nid_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_niv_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  }
  // tensor (index_ic_ten_): photon tensor IC set elsewhere (index_pt_gwdot, etc.)
}
```

- [ ] **Step 13.2: `BaryonsSpecies::ApplyInitialConditions`**

For all IC types, baryons track photons at leading order:

```cpp
void BaryonsSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;

  if (ctx.index_ic == mod->index_ic_ad_) {
    y[pv->index_pt_delta_b] = 3./4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[pv->index_pt_delta_b] = 3./4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_bi_) {
    // BI: baryon has isocurvature perturbation on top
    y[pv->index_pt_delta_b] = ctx.ppr->entropy_ini + 3./4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_nid_) {
    y[pv->index_pt_delta_b] = ctx.ppr->entropy_ini * ctx.fracnu / ctx.fracg / 8. * ctx.ktau_two;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  } else if (ctx.index_ic == mod->index_ic_niv_) {
    y[pv->index_pt_delta_b] = 3./4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  }
}
```

- [ ] **Step 13.3: `CDMSpecies::ApplyInitialConditions`**

CDM density tracks photons at 3/4, velocity = 0 in synchronous gauge:

```cpp
void CDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  // y[pv->index_pt_theta_cdm] = 0 already (zero-initialized) for synchronous gauge

  if (ctx.index_ic == mod->index_ic_ad_ || ctx.index_ic == mod->index_ic_bi_)
    y[pv->index_pt_delta_cdm] = 3./4. * ctx.delta_g_ic;
  else if (ctx.index_ic == mod->index_ic_cdi_)
    y[pv->index_pt_delta_cdm] = ctx.ppr->entropy_ini + 3./4. * ctx.delta_g_ic;
  else if (ctx.index_ic == mod->index_ic_nid_)
    y[pv->index_pt_delta_cdm] = -ctx.ppr->entropy_ini * ctx.fracnu * ctx.fracb / ctx.fracg / 80.
                                  * ctx.ktau_two * ctx.om * ctx.tau;
  else if (ctx.index_ic == mod->index_ic_niv_)
    y[pv->index_pt_delta_cdm] = -ctx.ppr->entropy_ini * 9./64. * ctx.fracnu * ctx.fracb / ctx.fracg
                                  * ctx.k * ctx.tau * ctx.om * ctx.tau;
}
```

- [ ] **Step 13.4: `UltraRelativisticSpecies::ApplyInitialConditions`**

UR uses `ctx.delta_ur`, `ctx.theta_ur`, `ctx.shear_ur`, `ctx.l3_ur` which the module computes:

```cpp
void UltraRelativisticSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  y[pv->index_pt_delta_ur] = ctx.delta_ur;
  y[pv->index_pt_theta_ur] = ctx.theta_ur;
  y[pv->index_pt_shear_ur] = ctx.shear_ur;
  y[pv->index_pt_l3_ur]    = ctx.l3_ur;
}
```

- [ ] **Step 13.5: IDM/IDR families `ApplyInitialConditions`**

IDM_DR (adiabatic only): `delta_idm_dr = 3/4 * delta_g`.
IDR: uses `delta_ur` / `theta_ur` (same as UR).
IDM_DRMD / IDR_DRMD: see existing code at lines 4997–5025.

```cpp
void IDM_DR_IDR_Species::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  if (ctx.index_ic != mod->index_ic_ad_) return;  // only adiabatic for IDM/IDR

  if (mod->GetBackground()->has_idm_dr == _TRUE_)
    y[pv->index_pt_delta_idm_dr] = 3./4. * ctx.delta_g_ic;
  if (mod->GetBackground()->has_idr == _TRUE_) {
    y[pv->index_pt_delta_idr] = ctx.delta_ur;
    y[pv->index_pt_theta_idr] = ctx.theta_ur;
  }
}
```

- [ ] **Step 13.6: `DCDMSpecies::ApplyInitialConditions`** (adiabatic only)

```cpp
void DCDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_) return;
  perturb_vector* pv = ctx.ppw->pv;
  y[pv->index_pt_delta_dcdm] = 3./4. * ctx.delta_g_ic;
  // theta_dcdm = 0 in synchronous gauge (already zero-initialized)
}
```

- [ ] **Step 13.7: `DarkRadiationSpecies::ApplyInitialConditions`** (adiabatic only)

```cpp
void DarkRadiationSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_) return;
  perturb_vector* pv = ctx.ppw->pv;
  y[pv->index_pt_F0_dr_sum] = ctx.delta_dr;
  // F0_dr_sum is rho-weighted sum; individual DR species ICs follow same formula
  // Read lines 5079 of perturb_initial_conditions for the exact dr IC
}
```

- [ ] **Step 13.8: `FluidSpecies::ApplyInitialConditions`**

```cpp
void FluidSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_) return;
  if (ctx.p_mod->GetBackground()->use_ppf == _TRUE_) return;  // PPF: Gamma_fld = 0
  perturb_vector* pv = ctx.ppw->pv;

  double w_fld, dw_over_da, integral;
  ctx.p_mod->GetBackgroundModule()->background_w_fld(ctx.a, &w_fld, &dw_over_da, &integral);

  y[pv->index_pt_delta_fld] =
    -ctx.ktau_two/4. * (1. + w_fld) * (4. - 3. * ctx.p_mod->GetBackground()->cs2_fld)
    / (4. - 6.*w_fld + 3.*ctx.p_mod->GetBackground()->cs2_fld)
    * ctx.ppr->curvature_ini * ctx.s2_squared;

  y[pv->index_pt_theta_fld] =
    -ctx.k * ctx.ktau_three / 4. * ctx.p_mod->GetBackground()->cs2_fld
    / (4. - 6.*w_fld + 3.*ctx.p_mod->GetBackground()->cs2_fld)
    * ctx.ppr->curvature_ini * ctx.s2_squared;
}
```

- [ ] **Step 13.9: `ScalarFieldSpecies::ApplyInitialConditions`**

SCF ICs are set to 0 (attractor not yet implemented):
```cpp
void ScalarFieldSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_) return;
  perturb_vector* pv = ctx.ppw->pv;
  y[pv->index_pt_phi_scf]       = 0.;
  y[pv->index_pt_phi_prime_scf] = 0.;
}
```

- [ ] **Step 13.10: `NCDMSpecies::ApplyInitialConditions`**

NCDM shares the UR ICs for all multipoles at leading order. The hierarchy (`index_pt_delta_ncdm` through `index_pt_shear_ncdm`) is set similarly to how UR is set.

Read lines 5397–5420 of `perturb_initial_conditions` for the NCDM hierarchy IC (it sets delta, theta, shear using delta_ur, theta_ur, shear_ur). Implement accordingly.

- [ ] **Step 13.11: Build**

```bash
make -j4 class 2>&1 | tail -5
```

- [ ] **Step 13.12: Commit**

```bash
git add species/photons.h species/photons.cpp \
        species/baryons.h species/baryons.cpp \
        species/cdm.h species/cdm.cpp \
        species/ultra_relativistic.h species/ultra_relativistic.cpp \
        species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp \
        species/dcdm.h species/dcdm.cpp \
        species/dark_radiation_species.h species/dark_radiation_species.cpp \
        species/ncdm_species.h species/ncdm_species.cpp \
        species/fluid.h species/fluid.cpp \
        species/scalar_field.h species/scalar_field.cpp
git commit -m "feat(Stage3): all species ApplyInitialConditions"
```

---

### Task 14: Stage 3 — Wire `ApplyInitialConditions` in module

**Files:**
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 14.1: Build `PerturbIcContext` in `perturb_initial_conditions`**

In `perturb_initial_conditions`, after the shared ratio computation block (lines ~4905–4943), construct the context. For each IC type, the module sets the seed values (`delta_g_ic`, `theta_g_ic`, `delta_ur`, etc.) before dispatching to species.

```cpp
PerturbIcContext ic_ctx;
ic_ctx.fracnu       = fracnu;
ic_ctx.fracg        = fracg;
ic_ctx.fracb        = fracb;
ic_ctx.fraccdm      = fraccdm;
ic_ctx.fracidm_drmd = fracidm_drmd;
ic_ctx.rho_m_over_rho_r = rho_m_over_rho_r;
ic_ctx.om           = om;
ic_ctx.ktau_two     = ktau_two;
ic_ctx.ktau_three   = ktau_three;
ic_ctx.s2_squared   = s2_squared;
ic_ctx.k            = k;
ic_ctx.tau          = tau;
ic_ctx.a            = a;
ic_ctx.a_prime_over_a = a_prime_over_a;
ic_ctx.index_ic     = index_ic;
ic_ctx.gauge        = ppt->gauge;
ic_ctx.ppw          = ppw;
ic_ctx.ppr          = ppr;
ic_ctx.p_mod        = this;
```

- [ ] **Step 14.2: Set seed values and dispatch for each IC type**

Replace the per-IC-type species blocks. For adiabatic (example):

```cpp
if ((ppt->has_ad == _TRUE_) && (index_ic == index_ic_ad_)) {
  // Module computes seed values (physics stays here)
  ic_ctx.delta_g_ic = - ktau_two/3. * (1. - om*tau/5.) * ppr->curvature_ini * s2_squared;
  ic_ctx.theta_g_ic = - k*ktau_three/36.
                       * (1. - 3.*(1.+5.*fracb-fracnu)/20./(1.-fracnu)*om*tau)
                       * ppr->curvature_ini * s2_squared;

  if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_) || (pba->has_dr == _TRUE_) || (pba->has_idr == _TRUE_)) {
    ic_ctx.delta_ur = ic_ctx.delta_g_ic;
    ic_ctx.theta_ur = - k*ktau_three/36./(4.*fracnu+15.)
                       * (4.*fracnu+11.+12.*s2_squared
                          - 3.*(8.*fracnu*fracnu+50.*fracnu+275.)/20./(2.*fracnu+15.)*tau*om)
                       * ppr->curvature_ini * s2_squared;
    ic_ctx.shear_ur = ktau_two/(45.+12.*fracnu) * (3.*s2_squared-1.)
                       * (1.+(4.*fracnu-5.)/4./(2.*fracnu+15.)*tau*om) * ppr->curvature_ini;
    ic_ctx.l3_ur    = ktau_three*2./7./(12.*fracnu+45.) * ppr->curvature_ini;
    if (pba->has_dr == _TRUE_) ic_ctx.delta_dr = ic_ctx.delta_ur;
  }
  ic_ctx.eta = ppr->curvature_ini * (1. - ktau_two/12./(15.+4.*fracnu)
                 * (5. + 4.*s2_squared*fracnu
                    - (16.*fracnu*fracnu+280.*fracnu+325)/10./(2.*fracnu+15.)*tau*om));

  // Dispatch to species
  for (auto& [name, sp] : all_species_)
    sp->ApplyInitialConditions(ppw->pv->y, ic_ctx);

  // Set synchronous metric IC (not species-specific)
  if (ppt->gauge == synchronous)
    ppw->pv->y[ppw->pv->index_pt_eta] = ic_ctx.eta;
}
```

Repeat this pattern for CDI, BI, NID, NIV — setting the appropriate `ic_ctx.delta_g_ic`, `ic_ctx.theta_g_ic`, `ic_ctx.delta_ur` for each IC type before dispatching.

- [ ] **Step 14.3: Keep Newtonian gauge transformation block unchanged**

The Newtonian gauge transformation block (lines ~5252–5400) stays in the module — it computes `alpha` from the now-set state and applies corrections. This block is NOT delegated to species.

- [ ] **Step 14.4: Build**

```bash
make -j4 class 2>&1 | tail -20
```

- [ ] **Step 14.5: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "feat(Stage3): wire ApplyInitialConditions dispatch in perturb_initial_conditions"
```

---

### Task 15: Stage 3 — Test

- [ ] **Step 15.1: Build**

```bash
make -j4 class && pip install . 2>&1 | tail -10
```

- [ ] **Step 15.2: Run `TEST_LEVEL=1`**

```bash
TEST_LEVEL=1 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -40
```
Expected: all tests pass.

- [ ] **Step 15.3: If tests fail — IC diagnostic**

IC bugs usually show up as wrong large-scale power (low-k). The most common mistake:
- Wrong `delta_g_ic` used by baryons (3/4 factor missing)
- Isocurvature mode applied to wrong species
- `ppw->pv->index_pt_*` is `-1` (unregistered) and writing to `y[-1]` is UB — check the perturbation vector initialization

Add `printf("IC: delta_b = %e, expected %e\n", y[pv->index_pt_delta_b], ...)` guards to narrow down.

- [ ] **Step 15.4: Run `TEST_LEVEL=2` (final validation)**

```bash
TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -40
```

This takes ~1 hour. Run only after `TEST_LEVEL=1` passes.

- [ ] **Step 15.5: Commit test results and notes**

```bash
git commit --allow-empty -m "test(Stage3): TEST_LEVEL=1 passing, Stage3 complete"
```
