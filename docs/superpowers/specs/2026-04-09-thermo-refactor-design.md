# Spec: Thermodynamics Module Refactor (Registry-Driven Physics)

**Status:** Draft
**Date:** 2026-04-09
**Goal:** Refactor `ThermodynamicsModule` to use the species registry (`all_species_`) for all species-specific thermodynamics operations, removing hardcoded `if(all_species_.count("IDM_DR_IDR") > 0)` and `if(all_species_.count("IDM_DRMD_IDR_DRMD") > 0)` guards.

## 1. Scope

Standard baryon-photon physics (RECFAST/HyRec recombination, reionization, kappa, g, Tb, cb2) stays in `ThermodynamicsModule` unchanged. Only the dark-species extensions are migrated:

| Species | Blocks removed from ThermodynamicsModule |
|---|---|
| `IDM_DR_IDR` | Index registration, per-row rate filling, optical depth integration, sign restore + visibility, T_idm_dr evolution, high-z extrapolation, IDR free-streaming time |
| `IDM_DRMD_IDR_DRMD` | Neff_bbn correction in `thermodynamics_helium_from_bbn` |

## 2. New Type: `ThermoTableContext`

A new header `species/thermo_context.h` (following the `perturb_context.h` precedent) defines:

```cpp
#pragma once
#include <vector>
#include "common.h"

class BackgroundModule;
struct precision;

struct ThermoTableContext {
  std::vector<double>&        table;            // thermodynamics_table_
  const std::vector<double>&  tau_table;        // local tau_table in thermodynamics_init
  const std::vector<double>&  z_table;          // z_table_
  int                         th_size;
  const BackgroundModule*     background_module;
  const precision*            ppr;
  ErrorMsg&                   error_message;
};
```

`tt_size` is obtained as `z_table.size()` wherever needed. `ThermodynamicsModule` constructs one `ThermoTableContext` instance after `tau_table` is populated and reuses it for all hook sites that need table access.

## 3. BaseSpecies Hooks

Six virtual methods added to `species/base_species.h`, all no-ops / zero by default:

```cpp
// ── Thermodynamics ────────────────────────────────────────────────────────

/**
 * Claim consecutive slots in pvecthermo. Called once during thermodynamics_indices().
 * Implementation must assign its index_th_* members and increment index_th.
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
 * Called once after array_smooth(). Responsible for: optical depth
 * integration, sign restoration, visibility functions, and temperature
 * evolution (e.g. T_idm_dr).
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
 * Compute and return the IDR free-streaming conformal time, or -1 if not applicable.
 * Called after tau_free_streaming_ is set in thermodynamics_init().
 */
virtual double GetIdrFreeStreamingTau(const ThermoTableContext& ctx) const { return -1.; }

/**
 * Return a correction to subtract from Neff_bbn in thermodynamics_helium_from_bbn().
 * pvecback is evaluated at the BBN redshift.
 */
virtual double NeffBbnCorrection(const double* pvecback) const { return 0.; }
```

`base_species.h` gains `#include "thermo_context.h"`.

## 4. IDM_DR_IDR_Species

### 4.1 Constructor

Add `const thermo& th` parameter (following the existing `const background& pba` pattern). `ReadIni` passes `context.th`.

```cpp
IDM_DR_IDR_Species(const background& pba, const thermo& th);
```

### 4.2 New private members

```cpp
const thermo& th_;

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

### 4.3 Hook implementations

**`RegisterThermodynamicsIndices`:** Claims 9 consecutive slots and stores them in the members above. Moves the existing `if(all_species_.count("IDM_DR_IDR") > 0)` block from `thermodynamics_indices()`.

**`FillThermodynamicsRow`:** Fills:
- `thermo_row[index_th_dmu_idm_dr_]` = IDR opacity to IDM scattering
- `thermo_row[index_th_ddmu_idm_dr_]` = `[Sinv * dmu_idm_dr]` (stored temporarily; overwritten in PostProcess)
- `thermo_row[index_th_dmu_idr_]` = IDR self-interaction rate

Uses `th_.a_idm_dr`, `th_.nindex_idm_dr`, `th_.b_idr`, `pba_.Omega0_idm_dr`, `pba_.Omega0_idr`, `pba_.h` and `background_module_->index_bg_rho_idr_`, `background_module_->index_bg_rho_idm_dr_`.

**`PostProcessThermodynamicsTable`:** All executed in sequence within one method:
1. Spline `ddmu_idm_dr` (Sinv*dmu) → `dddmu_idm_dr`; integrate → `tau_idm_dr`
2. Spline `dmu_idm_dr` → `dddmu_idm_dr`; integrate → `tau_idr`
3. For each row: flip sign of `tau_idm_dr` and `tau_idr`; compute `g_idm_dr = Sinv*dmu * exp(-tau_idm_dr)`
4. Spline `dmu_idm_dr` → `dddmu_idm_dr` (true second derivative); derive → true `ddmu_idm_dr`
5. Forward-Euler integration for `T_idm_dr` and `cidm_dr2` (with sub-stepping for stiff regime)

All array utility calls (`array_spline_table_line_to_line`, `array_integrate_spline_table_line_to_line`, `array_derive_spline_table_line_to_line`) are called directly from the species `.cpp` via `#include "arrays.h"`.

**`FillHighZExtrapolation`:** Fills all 9 pvecthermo slots for z ≥ z_initial. Moves the existing `if(all_species_.count("IDM_DR_IDR") > 0)` block from `thermodynamics_at_z()`.

**`GetIdrFreeStreamingTau`:** Runs the IDR free-streaming search loop (currently in `thermodynamics_init` after `tau_free_streaming_` is set). Returns the found `tau`; returns `-1.` if the search is not applicable.

`NeffBbnCorrection` is not overridden (inherits the 0-returning default).

## 5. IDM_DRMD_IDR_DRMD_Species

One override:

```cpp
double NeffBbnCorrection(const double* pvecback) const override {
  return pvecback[background_module_->index_bg_rho_idr_drmd_]
         / (7./8. * pow(4./11., 4./3.) * pvecback[background_module_->index_bg_rho_g_]);
}
```

`background_module_` is available via the inherited `SetBackgroundModule` hook. All other hooks inherit default no-ops.

## 6. ThermodynamicsModule Changes

### 6.1 `thermodynamics_indices()`

After the core index block:
```cpp
for (auto& [name, sp] : all_species_)
  sp->RegisterThermodynamicsIndices(index);
```
Remove the `if(all_species_.count("IDM_DR_IDR") > 0)` block.

### 6.2 Main loop (Loop 1) in `thermodynamics_init()`

At the end of each row iteration:
```cpp
for (auto& [name, sp] : all_species_)
  sp->FillThermodynamicsRow(index_tau, z_table_[index_tau], pvecback.data(),
                             thermodynamics_table_.data() + index_tau*th_size_);
```
Remove the `if(all_species_.count("IDM_DR_IDR") > 0)` block inside the loop.

Also remove the two post-loop IDM_DR_IDR blocks (optical depth spline+integrate section, and the sign restore + visibility block inside the final row loop).

### 6.3 Context construction

`ThermoTableContext` is constructed once immediately after `tau_table` is populated (~line 453), and kept alive through the rest of `thermodynamics_init()`:

```cpp
ThermoTableContext ctx{thermodynamics_table_, tau_table, z_table_,
                       th_size_, background_module_.get(), ppr, error_message_};
```

### 6.4 After `array_smooth()`

```cpp
for (auto& [name, sp] : all_species_)
  class_call(sp->PostProcessThermodynamicsTable(ctx), ctx.error_message, error_message_);
```

### 6.5 Free-streaming search (renumbered from 6.4)

Uses the same `ctx`. After `tau_free_streaming_` is set:
```cpp
for (auto& [name, sp] : all_species_) {
  double t = sp->GetIdrFreeStreamingTau(ctx);
  if (t > 0.) tau_idr_free_streaming_ = t;
}
```
Remove the `if(all_species_.count("IDM_DR_IDR") > 0)` block.

### 6.6 `thermodynamics_at_z()`, high-z branch

At end of the high-z block:
```cpp
for (auto& [name, sp] : all_species_)
  sp->FillHighZExtrapolation(z, pvecback, pvecthermo);
```
Remove the `if(all_species_.count("IDM_DR_IDR") > 0)` block.

### 6.7 `thermodynamics_helium_from_bbn()`

Replace the `if(all_species_.count("IDM_DRMD_IDR_DRMD") > 0)` block with:
```cpp
for (auto& [name, sp] : all_species_)
  Neff_bbn -= sp->NeffBbnCorrection(pvecback.data());
```

## 7. Verification Plan

1. **Build:** `make -j8 class`
2. **Regression:** `pytest -v python/test_class.py` covering standard ΛCDM and IDM_DR_IDR scenarios
3. **Consistency:** Compare `thermodynamics_table` output before and after for an IDM_DR_IDR scenario (e.g. using `explanatory.ini` with IDM parameters) to confirm bit-level agreement
