# Thermodynamics Module Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove all `if(all_species_.count("IDM_DR_IDR") > 0)` and `if(all_species_.count("IDM_DRMD_IDR_DRMD") > 0)` guards from `ThermodynamicsModule`, replacing them with registry loops over six virtual hooks on `BaseSpecies`.

**Architecture:** A new `ThermoTableContext` struct bundles the thermodynamics table + metadata for multi-row hook calls. `IDM_DR_IDR_Species` implements five hooks; `IDM_DRMD_IDR_DRMD_Species` implements one. `ThermodynamicsModule` becomes species-agnostic. `ThermodynamicsModule` public index members (`index_th_dmu_idm_dr_` etc.) are preserved for backward compatibility with `perturbations_module.cpp`; they are synced from the species after `RegisterThermodynamicsIndices`.

**Tech Stack:** C++17, CLASS cosmology code.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `species/thermo_context.h` | `ThermoTableContext` struct definition |
| Modify | `species/base_species.h` | Add `#include "thermo_context.h"` + 6 virtual hooks |
| Modify | `species/idm_dr_idr_species.h` | Add `const thermo& th_`, 9 index members, 5 hook declarations |
| Modify | `species/idm_dr_idr_species.cpp` | Implement 5 hooks; update `ReadIDRIni`/`ReadIDMDRIni` constructor calls; update `AddCouplingDerivs` to use own index members |
| Modify | `species/idm_drmd_idr_drmd_species.h` | Add `NeffBbnCorrection` declaration |
| Modify | `species/idm_drmd_idr_drmd_species.cpp` | Implement `NeffBbnCorrection` |
| Modify | `source/thermodynamics_module.h` | Remove 9 IDM_DR_IDR index members from public section |
| Modify | `source/thermodynamics_module.cpp` | Six sites: replace guarded blocks with registry loops; construct `ThermoTableContext`; include `idm_dr_idr_species.h` for index sync |

---

### Task 1: Create `ThermoTableContext` and add hooks to `BaseSpecies`

**Files:**
- Create: `species/thermo_context.h`
- Modify: `species/base_species.h`

- [ ] **Step 1: Create `species/thermo_context.h`**

```cpp
#pragma once
#include <vector>
#include "common.h"

class BackgroundModule;
struct precision;

/**
 * Context passed to thermodynamics hooks in BaseSpecies.
 * Bundles the thermodynamics table, coordinate arrays, and related metadata
 * so species can perform multi-row operations (spline, integrate, etc.)
 * without parameter explosion.
 */
struct ThermoTableContext {
  std::vector<double>&        table;            /**< thermodynamics_table_ (writable) */
  const std::vector<double>&  tau_table;        /**< conformal time for each row */
  const std::vector<double>&  z_table;          /**< redshift for each row */
  int                         th_size;          /**< number of columns (stride) */
  const BackgroundModule*     background_module;
  const precision*            ppr;
  ErrorMsg&                   error_message;
};
```

- [ ] **Step 2: Add `#include "thermo_context.h"` and the 6 hooks to `species/base_species.h`**

Add after the existing `#include "perturb_context.h"` line:
```cpp
#include "thermo_context.h"
```

Add after the existing `BackgroundDerivs` hook (inside the `// ── Background` section), a new `// ── Thermodynamics` section:

```cpp
  // ── Thermodynamics ────────────────────────────────────────────────────────

  /**
   * Claim consecutive slots in pvecthermo. Called once during thermodynamics_indices().
   * Implementation must assign its index_th_* members and increment index_th accordingly.
   */
  virtual void RegisterThermodynamicsIndices(int& index_th) {}

  /**
   * Fill per-row thermodynamics entries at redshift z.
   * thermo_row = thermodynamics_table_.data() + index_tau * th_size.
   * Called inside the main loop over index_tau in thermodynamics_init().
   */
  virtual void FillThermodynamicsRow(int index_tau, double z,
      const double* pvecback, double* thermo_row) {}

  /**
   * Post-process the completed thermodynamics table.
   * Called once after array_smooth(). Responsible for optical depths,
   * sign restoration, visibility functions, and temperature evolution.
   * Returns _SUCCESS_ or sets ctx.error_message and returns _FAILURE_.
   */
  virtual int PostProcessThermodynamicsTable(ThermoTableContext& ctx) { return _SUCCESS_; }

  /**
   * Fill pvecthermo entries for z >= z_table[tt_size-1] (high-z extrapolation).
   * Called inside thermodynamics_at_z() when z is above the table range.
   */
  virtual void FillHighZExtrapolation(double z,
      const double* pvecback, double* pvecthermo) const {}

  /**
   * Compute the IDR free-streaming conformal time.
   * Sets tau_out to the result and returns _SUCCESS_, or returns _FAILURE_.
   * tau_out = -1 means "not applicable for this species".
   * Called after tau_free_streaming_ is set in thermodynamics_init().
   *
   * Note: returns int (not double) so that class_call() error propagation works correctly.
   */
  virtual int GetIdrFreeStreamingTau(ThermoTableContext& ctx, double& tau_out) const {
    tau_out = -1.;
    return _SUCCESS_;
  }

  /**
   * Return the correction to subtract from Neff_bbn in thermodynamics_helium_from_bbn().
   * pvecback is evaluated at the BBN redshift.
   */
  virtual double NeffBbnCorrection(const double* pvecback) const { return 0.; }
```

- [ ] **Step 3: Build to verify compilation**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build, no errors.

- [ ] **Step 4: Commit**

```bash
git add species/thermo_context.h species/base_species.h
git commit -m "refactor: add ThermoTableContext and thermodynamics hooks to BaseSpecies"
```

---

### Task 2: Implement `RegisterThermodynamicsIndices` in `IDM_DR_IDR_Species`

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`
- Modify: `source/thermodynamics_module.h`
- Modify: `source/thermodynamics_module.cpp`

- [ ] **Step 1: Add 9 index members and hook declaration to `species/idm_dr_idr_species.h`**

Inside the `IDM_DR_IDR_Species` class, add in the `public` or `private` section (these must be accessible for the sync in `thermodynamics_module.cpp` — make them `public`):

```cpp
  void RegisterThermodynamicsIndices(int& index_th) override;

  // Set by RegisterThermodynamicsIndices; -1 = not registered
  int index_th_dmu_idm_dr_   = -1;
  int index_th_ddmu_idm_dr_  = -1;
  int index_th_dddmu_idm_dr_ = -1;
  int index_th_tau_idm_dr_   = -1;
  int index_th_tau_idr_      = -1;
  int index_th_g_idm_dr_     = -1;
  int index_th_cidm_dr2_     = -1;
  int index_th_Tidm_dr_      = -1;
  int index_th_dmu_idr_      = -1;
```

- [ ] **Step 2: Implement `RegisterThermodynamicsIndices` in `species/idm_dr_idr_species.cpp`**

Add after the constructor definition:

```cpp
void IDM_DR_IDR_Species::RegisterThermodynamicsIndices(int& index_th) {
  index_th_dmu_idm_dr_   = index_th++;
  index_th_ddmu_idm_dr_  = index_th++;
  index_th_dddmu_idm_dr_ = index_th++;
  index_th_tau_idm_dr_   = index_th++;
  index_th_tau_idr_      = index_th++;
  index_th_g_idm_dr_     = index_th++;
  index_th_cidm_dr2_     = index_th++;
  index_th_Tidm_dr_      = index_th++;
  index_th_dmu_idr_      = index_th++;
}
```

- [ ] **Step 3: Remove the 9 IDM_DR_IDR index members from `source/thermodynamics_module.h`**

Find and remove these lines from the public section of `ThermodynamicsModule`:

```cpp
  int index_th_dmu_idm_dr_;    /**< scattering rate of idr with idm_dr (i.e. idr opacity to idm_dr scattering) (units 1/Mpc) */
  int index_th_ddmu_idm_dr_;   /**< derivative of this scattering rate */
  int index_th_dddmu_idm_dr_;  /**< second derivative of this scattering rate */
  int index_th_dmu_idr_;       /**< idr self-interaction rate */
  int index_th_tau_idm_dr_;    /**< optical depth of idm_dr (due to interactions with idr) */
  int index_th_tau_idr_;       /**< optical depth of idr (due to self-interactions) */
  int index_th_g_idm_dr_;      /**< visibility function of idm_idr */
  int index_th_cidm_dr2_;      /**< interacting dark matter squared sound speed \f$ c_{dm}^2 \f$ */
  int index_th_Tidm_dr_;       /**< temperature of DM interacting with DR \f$ T_{idm_dr} \f$ */
```

- [ ] **Step 4: Wire `thermodynamics_indices()` in `source/thermodynamics_module.cpp`**

At the top of `thermodynamics_module.cpp`, add the include (after existing includes):
```cpp
#include "idm_dr_idr_species.h"
```

In `ThermodynamicsModule::thermodynamics_indices()`, find the block:

```cpp
  if(all_species_.count("IDM_DR_IDR") > 0){
    index_th_dmu_idm_dr_ = index;
    index++;
    index_th_ddmu_idm_dr_ = index;
    index++;
    index_th_dddmu_idm_dr_ = index;
    index++;
    index_th_tau_idm_dr_ = index;
    index++;
    index_th_tau_idr_ = index;
    index++;
    index_th_g_idm_dr_ = index;
    index++;
    index_th_cidm_dr2_ = index;
    index++;
    index_th_Tidm_dr_ = index;
    index++;
    index_th_dmu_idr_ = index;
    index++;
  }
```

Replace with:

```cpp
  for (auto& [name, sp] : all_species_)
    sp->RegisterThermodynamicsIndices(index);
```

- [ ] **Step 5: Update `AddCouplingDerivs` in `idm_dr_idr_species.cpp` to use own index members**

Find in `AddCouplingDerivs`:
```cpp
  auto* pth_mod  = ppaw.perturbations_module->GetThermodynamicsModule().get();
```
and every reference to `pth_mod->index_th_dmu_idm_dr_` and `pth_mod->index_th_cidm_dr2_`. Replace them to use the species' own members:

```cpp
  // Remove: auto* pth_mod = ppaw.perturbations_module->GetThermodynamicsModule().get();
  // (keep pth_mod only if needed for other accesses; if only used for indices, remove entirely)

  const double dmu_idm_dr = pvecthermo[index_th_dmu_idm_dr_];
```

And replace `pth_mod->index_th_cidm_dr2_` with `index_th_cidm_dr2_`, and `pth_mod->GetThermodynamics()->b_idr` / `pth_mod->GetThermodynamics()->nindex_idm_dr` with `th_.b_idr` / `th_.nindex_idm_dr`.

Note: `th_` member will be added in Task 3. For now only do the index replacements that use the new index members. The `GetThermodynamics()` calls referencing `b_idr` and `nindex_idm_dr` will be fixed in Task 3 once `th_` is available.

Specifically, replace in `AddCouplingDerivs`:

```cpp
  auto* pth_mod  = ppaw.perturbations_module->GetThermodynamicsModule().get();
  const double* pvecback   = ppw->pvecback;
  const double* pvecthermo = ppw->pvecthermo;
  auto* ppt = ppaw.perturbations_module->GetPerturbs();

  const double dmu_idm_dr = pvecthermo[pth_mod->index_th_dmu_idm_dr_];
```

with:

```cpp
  const double* pvecback   = ppw->pvecback;
  const double* pvecthermo = ppw->pvecthermo;
  auto* ppt = ppaw.perturbations_module->GetPerturbs();
  auto* pth_mod  = ppaw.perturbations_module->GetThermodynamicsModule().get();

  const double dmu_idm_dr = pvecthermo[index_th_dmu_idm_dr_];
```

And replace `pvecthermo[pth_mod->index_th_cidm_dr2_]` with `pvecthermo[index_th_cidm_dr2_]` (two occurrences).

Keep `pth_mod` for now since it's still used for `pth_mod->GetThermodynamics()->b_idr` — that will be cleaned up in Task 3.

- [ ] **Step 6: Build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build.

- [ ] **Step 7: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        source/thermodynamics_module.h source/thermodynamics_module.cpp
git commit -m "refactor: implement RegisterThermodynamicsIndices for IDM_DR_IDR species"
```

---

### Task 3: Add `const thermo& th_` and implement `FillThermodynamicsRow`

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`
- Modify: `source/thermodynamics_module.cpp`

- [ ] **Step 1: Add `const thermo& th_` member and hook declaration to `idm_dr_idr_species.h`**

Add to the constructor declaration:
```cpp
  explicit IDM_DR_IDR_Species(const background& pba, const thermo& th);
```

Add to private members:
```cpp
  const thermo& th_;
```

Add hook declaration:
```cpp
  void FillThermodynamicsRow(int index_tau, double z,
      const double* pvecback, double* thermo_row) override;
```

- [ ] **Step 2: Update constructor and `ReadIDRIni`/`ReadIDMDRIni` in `idm_dr_idr_species.cpp`**

Update constructor definition:
```cpp
IDM_DR_IDR_Species::IDM_DR_IDR_Species(const background& pba, const thermo& th)
  : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other)
  , pba_(pba)
  , th_(th)
{
  auto idm = std::make_unique<IDM_DRSpecies>(pba);
  auto idr = std::make_unique<IDRSpecies>(pba);
  idm_dr_ = idm.get();
  idr_    = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}
```

In both `ReadIDRIni` and `ReadIDMDRIni`, find:
```cpp
    instances.push_back(std::make_unique<IDM_DR_IDR_Species>(*pba));
```
and replace with:
```cpp
    instances.push_back(std::make_unique<IDM_DR_IDR_Species>(*pba, context.th));
```
(There are two such lines, one in each method — replace both.)

- [ ] **Step 3: Implement `FillThermodynamicsRow` in `idm_dr_idr_species.cpp`**

```cpp
void IDM_DR_IDR_Species::FillThermodynamicsRow(int /*index_tau*/, double z,
    const double* pvecback, double* thermo_row) {
  /* idr opacity to idm_dr scattering */
  thermo_row[index_th_dmu_idm_dr_] =
    th_.a_idm_dr * pow((1. + z)/1.e7, th_.nindex_idm_dr)
    * pba_.Omega0_idm_dr * pow(pba_.h, 2);

  /* [Sinv * dmu_idm_dr] stored temporarily in ddmu_idm_dr */
  thermo_row[index_th_ddmu_idm_dr_] =
    4./3. * pvecback[background_module_->index_bg_rho_idr_]
          / pvecback[background_module_->index_bg_rho_idm_dr_]
    * thermo_row[index_th_dmu_idm_dr_];

  /* idr self-interaction rate */
  thermo_row[index_th_dmu_idr_] =
    th_.b_idr * pow((1. + z)/1.e7, th_.nindex_idm_dr)
    * pba_.Omega0_idr * pow(pba_.h, 2);
}
```

- [ ] **Step 4: Finish cleaning up `AddCouplingDerivs` — replace `pth_mod->GetThermodynamics()->...` with `th_.*`**

In `AddCouplingDerivs`, replace:
```cpp
    const double dmu_idr = pth_mod->GetThermodynamics()->b_idr / pth_mod->GetThermodynamics()->a_idm_dr
                         * pba_.Omega0_idr / pba_.Omega0_idm_dr * dmu_idm_dr;
```
with:
```cpp
    const double dmu_idr = th_.b_idr / th_.a_idm_dr
                         * pba_.Omega0_idr / pba_.Omega0_idm_dr * dmu_idm_dr;
```

And replace:
```cpp
    const double tca_shear_idm_dr = 0.5 * 8./15. / dmu_idm_dr / ppt->alpha_idm_dr[0] * (theta_idm_dr + ctx.metric_shear);
    const double tca_slip_idm_dr = (pth_mod->GetThermodynamics()->nindex_idm_dr - 2./(1.+Sinv)) * ctx.a_prime_over_a * (theta_idm_dr - theta_idr)
```
with:
```cpp
    const double tca_shear_idm_dr = 0.5 * 8./15. / dmu_idm_dr / ppt->alpha_idm_dr[0] * (theta_idm_dr + ctx.metric_shear);
    const double tca_slip_idm_dr = (th_.nindex_idm_dr - 2./(1.+Sinv)) * ctx.a_prime_over_a * (theta_idm_dr - theta_idr)
```

Now remove the `pth_mod` variable entirely from `AddCouplingDerivs` (it is no longer used after these replacements).

- [ ] **Step 5: Wire `FillThermodynamicsRow` in `thermodynamics_module.cpp`**

In `ThermodynamicsModule::thermodynamics_init()`, find the main loop:
```cpp
  for (index_tau = 0; index_tau < tt_size_; index_tau++) {

    class_call(background_module_->background_at_tau(...), ...);

    R = 3./4.*pvecback[...index_bg_rho_b_]/pvecback[...index_bg_rho_g_];

    thermodynamics_table_[index_tau*th_size_+index_th_ddkappa_] = -1./R*...;

    if(all_species_.count("IDM_DR_IDR") > 0) {
      // ... three lines filling dmu_idm_dr, ddmu_idm_dr, dmu_idr
    }
  }
```

Remove the `if(all_species_.count("IDM_DR_IDR") > 0)` block inside the loop and replace with:

```cpp
    for (auto& [name, sp] : all_species_)
      sp->FillThermodynamicsRow(index_tau, z_table_[index_tau], pvecback.data(),
                                 thermodynamics_table_.data() + index_tau*th_size_);
```

- [ ] **Step 6: Build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build.

- [ ] **Step 7: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        source/thermodynamics_module.cpp
git commit -m "refactor: implement FillThermodynamicsRow for IDM_DR_IDR species"
```

---

### Task 4: Implement `PostProcessThermodynamicsTable` and construct `ThermoTableContext`

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`
- Modify: `source/thermodynamics_module.cpp`

- [ ] **Step 1: Add hook declaration to `idm_dr_idr_species.h`**

```cpp
  int PostProcessThermodynamicsTable(ThermoTableContext& ctx) override;
```

- [ ] **Step 2: Add `#include "arrays.h"` to `idm_dr_idr_species.cpp`**

Add after existing includes:
```cpp
#include "arrays.h"
```

- [ ] **Step 3: Implement `PostProcessThermodynamicsTable` in `idm_dr_idr_species.cpp`**

```cpp
int IDM_DR_IDR_Species::PostProcessThermodynamicsTable(ThermoTableContext& ctx) {
  const int tt_size = (int)ctx.z_table.size();
  double* table     = ctx.table.data();
  ErrorMsg& errmsg  = ctx.error_message;

  // ── Phase 1: optical depth integration ────────────────────────────────────

  /* spline [Sinv*dmu_idm_dr] (currently in ddmu column) → dddmu (workspace) */
  class_call(array_spline_table_line_to_line(
      ctx.tau_table.data(), tt_size, table, ctx.th_size,
      index_th_ddmu_idm_dr_, index_th_dddmu_idm_dr_,
      _SPLINE_EST_DERIV_, errmsg),
    errmsg, errmsg);

  /* integrate [Sinv*dmu_idm_dr] → -tau_idm_dr (sign flipped later) */
  class_call(array_integrate_spline_table_line_to_line(
      ctx.tau_table.data(), tt_size, table, ctx.th_size,
      index_th_ddmu_idm_dr_, index_th_dddmu_idm_dr_, index_th_tau_idm_dr_,
      errmsg),
    errmsg, errmsg);

  /* spline dmu_idm_dr → dddmu (workspace, reused) */
  class_call(array_spline_table_line_to_line(
      ctx.tau_table.data(), tt_size, table, ctx.th_size,
      index_th_dmu_idm_dr_, index_th_dddmu_idm_dr_,
      _SPLINE_EST_DERIV_, errmsg),
    errmsg, errmsg);

  /* integrate dmu_idm_dr → -tau_idr (sign flipped later) */
  class_call(array_integrate_spline_table_line_to_line(
      ctx.tau_table.data(), tt_size, table, ctx.th_size,
      index_th_dmu_idm_dr_, index_th_dddmu_idm_dr_, index_th_tau_idr_,
      errmsg),
    errmsg, errmsg);

  // ── Phase 2: sign restore + visibility function ────────────────────────────

  for (int index_tau = 0; index_tau < tt_size; index_tau++) {
    double* row = table + index_tau * ctx.th_size;
    row[index_th_tau_idm_dr_] *= -1.;
    row[index_th_tau_idr_]    *= -1.;
    row[index_th_g_idm_dr_]   = row[index_th_ddmu_idm_dr_]
                               * exp(-row[index_th_tau_idm_dr_]);
  }

  // ── Phase 3: true spline derivatives of dmu_idm_dr ────────────────────────

  /* second derivative of dmu_idm_dr w.r.t. tau */
  class_call(array_spline_table_line_to_line(
      ctx.tau_table.data(), tt_size, table, ctx.th_size,
      index_th_dmu_idm_dr_, index_th_dddmu_idm_dr_,
      _SPLINE_EST_DERIV_, errmsg),
    errmsg, errmsg);

  /* first derivative of dmu_idm_dr w.r.t. tau (overwrites ddmu with true value) */
  class_call(array_derive_spline_table_line_to_line(
      ctx.tau_table.data(), tt_size, table, ctx.th_size,
      index_th_dmu_idm_dr_, index_th_dddmu_idm_dr_, index_th_ddmu_idm_dr_,
      errmsg),
    errmsg, errmsg);

  // ── Phase 4: T_idm_dr and c_s^2 evolution ────────────────────────────────

  double Gamma_heat_idm_dr, dTdz_idm_dr, T_idm_dr, T_idr, dz, T_adia, z_adia, tau, z;
  int last_index_back = ctx.background_module->bg_size_ - 1;
  std::vector<double> pvecback(ctx.background_module->bg_size_);

  /* (A) initial value at maximum z (minimum tau) */
  z = ctx.z_table[tt_size - 1];

  class_call(ctx.background_module->background_tau_of_z(z, &tau),
    ctx.background_module->error_message_, errmsg);
  class_call(ctx.background_module->background_at_tau(
      tau, pba_.short_info, pba_.inter_normal, &last_index_back, pvecback.data()),
    ctx.background_module->error_message_, errmsg);

  Gamma_heat_idm_dr = 2.*pba_.Omega0_idr*pow(pba_.h,2)*th_.a_idm_dr
    *pow((1.+z),(th_.nindex_idm_dr+1.))/pow(1.e7,th_.nindex_idm_dr);

  if (Gamma_heat_idm_dr > 1.e-3*pvecback[ctx.background_module->index_bg_a_]
                                 *pvecback[ctx.background_module->index_bg_H_]) {
    T_idr = pba_.T_idr*(1.+z);
    T_idm_dr = T_idr;
    dTdz_idm_dr = pba_.T_idr;
  } else {
    T_idr = pba_.T_idr*(1.+z);
    T_idm_dr = Gamma_heat_idm_dr
      / (pvecback[ctx.background_module->index_bg_a_]*pvecback[ctx.background_module->index_bg_H_])
      / (1. + Gamma_heat_idm_dr/(pvecback[ctx.background_module->index_bg_a_]*pvecback[ctx.background_module->index_bg_H_]))
      * T_idr;
    dTdz_idm_dr = 2.*T_idm_dr
      - Gamma_heat_idm_dr/pvecback[ctx.background_module->index_bg_H_]*(T_idr - T_idm_dr);
  }

  table[(tt_size - 1)*ctx.th_size + index_th_Tidm_dr_]  = T_idm_dr;
  table[(tt_size - 1)*ctx.th_size + index_th_cidm_dr2_] =
    _k_B_*T_idm_dr/_eV_/th_.m_idm*(1. + dTdz_idm_dr/3./T_idm_dr);

  T_adia = T_idm_dr;
  z_adia = z;

  /* (B) iterate from high-z to low-z */
  for (int index_tau = tt_size - 2; index_tau >= 0; index_tau--) {

    /* (B1) tight coupling: Gamma >> H */
    if (Gamma_heat_idm_dr > 1.e3*pvecback[ctx.background_module->index_bg_a_]
                                  *pvecback[ctx.background_module->index_bg_H_]) {
      z = ctx.z_table[index_tau];
      T_idr = pba_.T_idr*(1.+z);
      T_idm_dr = T_idr;
      Gamma_heat_idm_dr = 2.*pba_.Omega0_idr*pow(pba_.h,2)*th_.a_idm_dr
        *pow((1.+z),(th_.nindex_idm_dr+1.))/pow(1.e7,th_.nindex_idm_dr);
      class_call(ctx.background_module->background_tau_of_z(z, &tau),
        ctx.background_module->error_message_, errmsg);
      class_call(ctx.background_module->background_at_tau(
          tau, pba_.short_info, pba_.inter_normal, &last_index_back, pvecback.data()),
        ctx.background_module->error_message_, errmsg);
      dTdz_idm_dr = pba_.T_idr;
    }

    /* (B2) intermediate: integrate dT/dz */
    else if (Gamma_heat_idm_dr > 1.e-3*pvecback[ctx.background_module->index_bg_a_]
                                        *pvecback[ctx.background_module->index_bg_H_]) {
      dz = ctx.z_table[index_tau + 1] - ctx.z_table[index_tau];

      /* (B2a) non-stiff: single Euler step */
      if (dz < pvecback[ctx.background_module->index_bg_H_]/Gamma_heat_idm_dr/10.) {
        z = ctx.z_table[index_tau];
        T_idr = pba_.T_idr*(1.+z);
        T_idm_dr -= dTdz_idm_dr*dz;
        Gamma_heat_idm_dr = 2.*pba_.Omega0_idr*pow(pba_.h,2)*th_.a_idm_dr
          *pow((1.+z),(th_.nindex_idm_dr+1.))/pow(1.e7,th_.nindex_idm_dr);
        class_call(ctx.background_module->background_tau_of_z(z, &tau),
          ctx.background_module->error_message_, errmsg);
        class_call(ctx.background_module->background_at_tau(
            tau, pba_.short_info, pba_.inter_normal, &last_index_back, pvecback.data()),
          ctx.background_module->error_message_, errmsg);
        dTdz_idm_dr = 2.*pvecback[ctx.background_module->index_bg_a_]*T_idm_dr
          - Gamma_heat_idm_dr/pvecback[ctx.background_module->index_bg_H_]*(T_idr - T_idm_dr);
      }

      /* (B2b) stiff: sub-stepped Euler */
      else {
        int N_sub_steps = (int)(dz/(pvecback[ctx.background_module->index_bg_H_]/Gamma_heat_idm_dr/10.)) + 1;
        double dz_sub_step = dz/N_sub_steps;
        for (int n = 0; n < N_sub_steps; n++) {
          z -= dz_sub_step;
          if (n == (N_sub_steps - 1)) z = ctx.z_table[index_tau];
          T_idr = pba_.T_idr*(1.+z);
          T_idm_dr -= dTdz_idm_dr*dz_sub_step;
          Gamma_heat_idm_dr = 2.*pba_.Omega0_idr*pow(pba_.h,2)*th_.a_idm_dr
            *pow((1.+z),(th_.nindex_idm_dr+1.))/pow(1.e7,th_.nindex_idm_dr);
          class_call(ctx.background_module->background_tau_of_z(z, &tau),
            ctx.background_module->error_message_, errmsg);
          class_call(ctx.background_module->background_at_tau(
              tau, pba_.short_info, pba_.inter_normal, &last_index_back, pvecback.data()),
            ctx.background_module->error_message_, errmsg);
          dTdz_idm_dr = 2.*pvecback[ctx.background_module->index_bg_a_]*T_idm_dr
            - Gamma_heat_idm_dr/pvecback[ctx.background_module->index_bg_H_]*(T_idr - T_idm_dr);
        }
      }

      /* update T_adia, z_adia after any B2 step (B2a or B2b) */
      T_adia = T_idm_dr;
      z_adia = z;
    }

    /* (B3) decoupled: T_idm_dr ~ a^{-2} */
    else {
      z = ctx.z_table[index_tau];
      T_idr = pba_.T_idr*(1.+z);
      T_idm_dr = T_adia * pow((1.+z)/(1.+z_adia), 2);
      Gamma_heat_idm_dr = 2.*pba_.Omega0_idr*pow(pba_.h,2)*th_.a_idm_dr
        *pow((1.+z),(th_.nindex_idm_dr+1.))/pow(1.e7,th_.nindex_idm_dr);
      class_call(ctx.background_module->background_tau_of_z(z, &tau),
        ctx.background_module->error_message_, errmsg);
      class_call(ctx.background_module->background_at_tau(
          tau, pba_.short_info, pba_.inter_normal, &last_index_back, pvecback.data()),
        ctx.background_module->error_message_, errmsg);
      dTdz_idm_dr = 2./(1+z)*T_idm_dr;
    }

    table[index_tau*ctx.th_size + index_th_Tidm_dr_]  = T_idm_dr;
    table[index_tau*ctx.th_size + index_th_cidm_dr2_] =
      _k_B_*T_idm_dr/_eV_/th_.m_idm*(1. + dTdz_idm_dr/3./T_idm_dr);
  }

  return _SUCCESS_;
}
```

- [ ] **Step 4: Construct `ThermoTableContext` in `thermodynamics_init()` and call `PostProcessThermodynamicsTable`**

In `thermodynamics_module.cpp`, right after the `tau_table` is populated (the `for (index_tau = 0; ...)` loop that fills `tau_table` via `background_tau_of_z`, around line 449–453), add:

```cpp
  ThermoTableContext thermo_ctx{thermodynamics_table_, tau_table, z_table_,
                                 th_size_, background_module_.get(), ppr, error_message_};
```

Then find the `array_smooth` call:
```cpp
  class_call(array_smooth(thermodynamics_table_.data(),
                          th_size_,
                          tt_size_,
                          index_th_rate_,
                          ppr->thermo_rate_smoothing_radius,
                          error_message_), ...);
```

After it, add:
```cpp
  for (auto& [name, sp] : all_species_)
    class_call(sp->PostProcessThermodynamicsTable(thermo_ctx),
               thermo_ctx.error_message, error_message_);
```

- [ ] **Step 5: Remove the three old IDM_DR_IDR blocks from `thermodynamics_init()`**

Remove the following three blocks:

**Block A** (post-loop, lines ~528–580) — the `if(all_species_.count("IDM_DR_IDR") > 0)` block containing the four `array_spline_table_line_to_line` / `array_integrate_spline_table_line_to_line` calls.

**Block B** (inside the final forward loop, lines ~773–786) — the `if(all_species_.count("IDM_DR_IDR") > 0)` block that restores signs and computes `g_idm_dr`. Remove only the `if` block; keep the surrounding loop logic.

**Block C** (lines ~800–970) — the `if(all_species_.count("IDM_DR_IDR") > 0)` block for T_idm_dr evolution.

- [ ] **Step 6: Build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build.

- [ ] **Step 7: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        source/thermodynamics_module.cpp
git commit -m "refactor: implement PostProcessThermodynamicsTable for IDM_DR_IDR species"
```

---

### Task 5: Implement `FillHighZExtrapolation`

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`
- Modify: `source/thermodynamics_module.cpp`

- [ ] **Step 1: Add hook declaration to `idm_dr_idr_species.h`**

```cpp
  void FillHighZExtrapolation(double z,
      const double* pvecback, double* pvecthermo) const override;
```

- [ ] **Step 2: Implement `FillHighZExtrapolation` in `idm_dr_idr_species.cpp`**

The method needs access to `thermodynamics_table_` for the last two rows (for `tau_idm_dr` / `tau_idr` linear extrapolation) and `z_table_`. Pass these via the species storing a pointer to the module, or — simpler — cache what's needed. But `FillHighZExtrapolation` is `const` and doesn't receive the table context. 

Solution: cache the table's raw pointer and key dimensions at the end of `PostProcessThermodynamicsTable`. `thermodynamics_table_` is a member of `ThermodynamicsModule` and is not resized after `PostProcess`, so its pointer remains valid.

In `idm_dr_idr_species.h`, add private members:
```cpp
  // Cached after PostProcessThermodynamicsTable for use in FillHighZExtrapolation.
  // Points into thermodynamics_table_ which lives in ThermodynamicsModule.
  const double* thermo_table_data_ = nullptr; // thermodynamics_table_.data()
  int           tt_size_cached_    = 0;
  int           th_size_cached_    = 0;
  double        z_last_            = 0.;
  double        z_second_last_     = 0.;
```

At the end of `PostProcessThermodynamicsTable`, before `return _SUCCESS_`:
```cpp
  // Cache for FillHighZExtrapolation
  thermo_table_data_ = ctx.table.data();
  tt_size_cached_    = tt_size;
  th_size_cached_    = ctx.th_size;
  z_last_            = ctx.z_table[tt_size - 1];
  z_second_last_     = ctx.z_table[tt_size - 2];
```

Now implement `FillHighZExtrapolation`:
```cpp
void IDM_DR_IDR_Species::FillHighZExtrapolation(double z,
    const double* pvecback, double* pvecthermo) const {
  const double* row_last   = thermo_table_data_ + (tt_size_cached_ - 1) * th_size_cached_;
  const double* row_second = thermo_table_data_ + (tt_size_cached_ - 2) * th_size_cached_;

  /* dmu_idm_dr and its derivatives */
  pvecthermo[index_th_dmu_idm_dr_] =
    th_.a_idm_dr * pow((1. + z)/1.e7, th_.nindex_idm_dr)
    * pba_.Omega0_idm_dr * pow(pba_.h, 2);

  pvecthermo[index_th_ddmu_idm_dr_] =
    -pvecback[background_module_->index_bg_H_] * th_.nindex_idm_dr/(1 + z)
    * pvecthermo[index_th_dmu_idm_dr_];

  pvecthermo[index_th_dddmu_idm_dr_] =
    (pvecback[background_module_->index_bg_H_]*pvecback[background_module_->index_bg_H_]/(1.+z)
     - pvecback[background_module_->index_bg_H_prime_])
    * th_.nindex_idm_dr/(1.+z) * pvecthermo[index_th_dmu_idm_dr_];

  /* dmu_idr self-interaction rate */
  pvecthermo[index_th_dmu_idr_] =
    th_.b_idr * pow((1. + z)/1.e7, th_.nindex_idm_dr)
    * pba_.Omega0_idr * pow(pba_.h, 2);

  /* linear extrapolation of tau_idm_dr and tau_idr from last two table rows */
  pvecthermo[index_th_tau_idm_dr_] =
    row_last[index_th_tau_idm_dr_]
    + (row_last[index_th_tau_idm_dr_] - row_second[index_th_tau_idm_dr_])
      * (z - z_last_) / (z_last_ - z_second_last_);

  pvecthermo[index_th_tau_idr_] =
    row_last[index_th_tau_idr_]
    + (row_last[index_th_tau_idr_] - row_second[index_th_tau_idr_])
      * (z - z_last_) / (z_last_ - z_second_last_);

  /* visibility function: constant extrapolation from last row */
  pvecthermo[index_th_g_idm_dr_] = row_last[index_th_g_idm_dr_];

  /* IDM sound speed and temperature */
  pvecthermo[index_th_cidm_dr2_] = 4.*_k_B_*pba_.T_idr*(1.+z)/_eV_/3./th_.m_idm;
  pvecthermo[index_th_Tidm_dr_]  = pba_.T_idr*(1.+z);
}
```

- [ ] **Step 3: Wire `FillHighZExtrapolation` in `thermodynamics_at_z()`**

In `thermodynamics_module.cpp`, inside `thermodynamics_at_z()`, find the high-z branch (the `if (z >= z_table_[tt_size_ - 1])` block). At the end of this block, just before the closing `}`, remove the existing:

```cpp
    if(all_species_.count("IDM_DR_IDR") > 0){
      // ... (the block filling pvecthermo[index_th_dmu_idm_dr_] etc.)
    }
```

and replace with:

```cpp
    for (auto& [name, sp] : all_species_)
      sp->FillHighZExtrapolation(z, pvecback, pvecthermo);
```

- [ ] **Step 4: Build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build.

- [ ] **Step 5: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        source/thermodynamics_module.cpp
git commit -m "refactor: implement FillHighZExtrapolation for IDM_DR_IDR species"
```

---

### Task 6: Implement `GetIdrFreeStreamingTau`

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`
- Modify: `source/thermodynamics_module.cpp`

- [ ] **Step 1: Add hook declaration to `idm_dr_idr_species.h`**

```cpp
  int GetIdrFreeStreamingTau(ThermoTableContext& ctx, double& tau_out) const override;
```

- [ ] **Step 2: Implement `GetIdrFreeStreamingTau` in `idm_dr_idr_species.cpp`**

```cpp
int IDM_DR_IDR_Species::GetIdrFreeStreamingTau(ThermoTableContext& ctx, double& tau_out) const {
  const int tt_size   = (int)ctx.z_table.size();
  const double* table = ctx.table.data();
  double tau;
  int index_tau;

  if (th_.nindex_idm_dr >= 2) {
    /* start one below top of table, search backwards (decreasing index = later times) */
    index_tau = tt_size - 2;
    class_call(ctx.background_module->background_tau_of_z(ctx.z_table[index_tau], &tau),
               ctx.background_module->error_message_, ctx.error_message);

    while ((1./table[index_tau*ctx.th_size + index_th_dmu_idm_dr_]/tau
            < ctx.ppr->idr_streaming_trigger_tau_c_over_tau)
           && index_tau > 0) {
      index_tau--;
      class_call(ctx.background_module->background_tau_of_z(ctx.z_table[index_tau], &tau),
                 ctx.background_module->error_message_, ctx.error_message);
    }
  } else {
    /* start at bottom of table, search forwards (increasing index = earlier times) */
    index_tau = 0;
    class_call(ctx.background_module->background_tau_of_z(ctx.z_table[index_tau], &tau),
               ctx.background_module->error_message_, ctx.error_message);

    while ((1./table[index_tau*ctx.th_size + index_th_dmu_idm_dr_]/tau
            < ctx.ppr->idr_streaming_trigger_tau_c_over_tau)
           && index_tau < tt_size - 1) {
      index_tau++;
      class_call(ctx.background_module->background_tau_of_z(ctx.z_table[index_tau], &tau),
                 ctx.background_module->error_message_, ctx.error_message);
    }
  }

  tau_out = tau;
  return _SUCCESS_;
}

- [ ] **Step 3: Wire `GetIdrFreeStreamingTau` in `thermodynamics_init()`**

Find in `thermodynamics_module.cpp` the block starting with:
```cpp
  /** - Find interacting dark radiation free-streaming time */
  int index_tau_fs = index_tau;
  double tau_idm_dr_fs=0.;
  if(all_species_.count("IDM_DR_IDR") > 0) {
    // ... the search loop and tau_idr_free_streaming_ = tau; assignments
  }
```

Replace the entire `if(all_species_.count("IDM_DR_IDR") > 0)` block with:

```cpp
  for (auto& [name, sp] : all_species_) {
    double t = -1.;
    class_call(sp->GetIdrFreeStreamingTau(thermo_ctx, t),
               thermo_ctx.error_message, error_message_);
    if (t > 0.) tau_idr_free_streaming_ = t;
  }
```

- [ ] **Step 4: Build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build.

- [ ] **Step 5: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp \
        source/thermodynamics_module.cpp
git commit -m "refactor: implement GetIdrFreeStreamingTau for IDM_DR_IDR species"
```

---

### Task 7: Implement `NeffBbnCorrection` for `IDM_DRMD_IDR_DRMD_Species`

**Files:**
- Modify: `species/idm_drmd_idr_drmd_species.h`
- Modify: `species/idm_drmd_idr_drmd_species.cpp`
- Modify: `source/thermodynamics_module.cpp`

- [ ] **Step 1: Add hook declaration to `idm_drmd_idr_drmd_species.h`**

```cpp
  double NeffBbnCorrection(const double* pvecback) const override;
```

- [ ] **Step 2: Implement `NeffBbnCorrection` in `idm_drmd_idr_drmd_species.cpp`**

```cpp
double IDM_DRMD_IDR_DRMD_Species::NeffBbnCorrection(const double* pvecback) const {
  return pvecback[background_module_->index_bg_rho_idr_drmd_]
         / (7./8. * pow(4./11., 4./3.) * pvecback[background_module_->index_bg_rho_g_]);
}
```

- [ ] **Step 3: Wire `NeffBbnCorrection` in `thermodynamics_helium_from_bbn()`**

Find the block:
```cpp
  /**DRMD**/
  if(all_species_.count("IDM_DRMD_IDR_DRMD") > 0){
    Neff_bbn -= (pvecback[background_module_->index_bg_rho_idr_drmd_])/(7./8.*pow(4./11.,4./3.)*pvecback[background_module_->index_bg_rho_g_]);
  }
```

Replace with:
```cpp
  for (auto& [name, sp] : all_species_)
    Neff_bbn -= sp->NeffBbnCorrection(pvecback.data());
```

Note: `thermodynamics_helium_from_bbn()` uses a local `pvecback` — verify the variable name is `pvecback` (not `pvecback.data()`). Adjust as needed to match the local variable type in that function.

- [ ] **Step 4: Build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build.

- [ ] **Step 5: Commit**

```bash
git add species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp \
        source/thermodynamics_module.cpp
git commit -m "refactor: implement NeffBbnCorrection for IDM_DRMD_IDR_DRMD species"
```

---

### Task 8: Final Verification

**Files:** None modified — run only.

- [ ] **Step 1: Full build**

```bash
cd /Users/au192734/Projects/class_gemini && make -j8 class
```
Expected: Clean build, zero warnings introduced by this refactor.

- [ ] **Step 2: Run standard regression tests**

```bash
cd /Users/au192734/Projects/class_gemini && python -m pytest python/test_class.py -v 2>&1 | tail -20
```
Expected: All tests pass (same as before refactor).

- [ ] **Step 3: Verify IDM_DR_IDR scenario produces identical output**

Create a temporary ini file `/tmp/test_idm.ini`:
```ini
Omega_b = 0.022
Omega_cdm = 0.1
h = 0.67
N_idr = 0.34
Omega_idm_dr = 0.001
a_idm_dr = 1e-4
nindex_idm_dr = 4
m_idm = 1e8
b_idr = 0.5
output = tCl
```

Checkout the commit just before this refactor, run CLASS, save output. Then checkout current HEAD, run CLASS again, diff the outputs:

```bash
# On current HEAD (after refactor):
./class /tmp/test_idm.ini
mv output/test_idm_cl.dat /tmp/cl_after.dat

# Verify cl_after.dat exists and is non-empty
wc -l /tmp/cl_after.dat
```

If you have the pre-refactor output, compare with:
```bash
diff /tmp/cl_before.dat /tmp/cl_after.dat
```
Expected: No differences (or only floating-point rounding at the last digit).

- [ ] **Step 4: Confirm no species-specific guards remain in `thermodynamics_module.cpp`**

```bash
grep -n 'IDM_DR_IDR\|IDM_DRMD_IDR_DRMD' source/thermodynamics_module.cpp
```
Expected: No output (zero matches).

- [ ] **Step 5: Final commit**

```bash
git commit --allow-empty -m "refactor: complete thermodynamics module species refactor"
```
