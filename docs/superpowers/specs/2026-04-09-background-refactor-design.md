# Spec: Phase 3 - Background Module Refactor (Registry-Driven Physics)

**Status:** Draft
**Date:** 2026-04-09
**Goal:** Refactor `BackgroundModule` to use the species registry (`all_species_`) for index registration, initial conditions, and ODE derivatives, removing hardcoded `if (pba->has_xxx)` guards.

## 1. Architecture

The refactor follows a "Physics Hook" pattern. `BackgroundModule` will no longer "know" about specific species. Instead, it will iterate through the `all_species_` map and call standardized virtual methods on each `BaseSpecies` instance.

### 1.1. BaseSpecies Hooks

Add the following virtual methods to `BaseSpecies`:

```cpp
/**
 * Register background indices for this species.
 * @param index_bg Running index in the background vector.
 */
virtual void RegisterBackgroundIndices(int& index_bg) {}

/**
 * Register integration indices for this species.
 * @param index_bi Running index in the integration vector.
 */
virtual void RegisterIntegrationIndices(int& index_bi) {}

/**
 * Set initial conditions for background ODE variables.
 * @param a_rel Scale factor relative to today (a/a_today).
 * @param y Vector of variables to be integrated.
 */
virtual void SetBackgroundInitialConditions(double a_rel, double* y) {}

/**
 * Compute derivatives for background ODE variables.
 * @param tau Conformal time.
 * @param y Vector of integrated variables.
 * @param dy Vector of derivatives.
 * @param pvecback Vector of all background quantities.
 */
virtual void BackgroundDerivs(double tau, double* y, double* dy, double* pvecback) {}
```

## 2. Implementation Strategy

### 2.1. BackgroundModule::background_functions

Already partially refactored in Phase 2. It uses `sp->ComputeBackground()` and an `accumulate` helper.
**Change:** Ensure all species (including `Fluid` and `ScalarField`) are correctly integrated into this loop or have their dependencies (like `w_fld`) handled before the loop.

### 2.2. BackgroundModule::background_indices

**Current State:** Contains hardcoded blocks for every species to assign `index_bg_rho_xxx_` and `index_bi_rho_xxx_`.
**Refactor:** 
1. Call `sp->RegisterBackgroundIndices(index_bg)` in a loop for all species.
2. Call `sp->RegisterIntegrationIndices(index_bi)` in a loop for all species.
3. Move species-specific logic (like the special `ScalarField` or `Fluid` index assignments) into their respective classes.

### 2.3. BackgroundModule::background_initial_conditions

**Current State:** Hardcoded blocks for `Omega_rad` calculation and species-specific ICs.
**Refactor:**
1. Dynamically calculate `Omega_rad` by iterating through `all_species_` and checking `sp->energy_type() == Radiation`.
2. Call `sp->SetBackgroundInitialConditions()` in a loop.

### 2.4. BackgroundModule::background_derivs_member

**Current State:** Hardcoded `rho_M` calculation and species-specific derivative calls.
**Refactor:**
1. Dynamically calculate `rho_M` by iterating through `all_species_` and checking `sp->energy_type() == Matter`.
2. Call `sp->BackgroundDerivs()` in a loop (already exists for some species, but needs full coverage).

## 3. Species-Specific Migrations

### 3.1. ScalarField
- Move potential calculations (`V_scf`, `dV_scf`, etc.) into `ScalarFieldSpecies`.
- Move attractor initial condition logic into `ScalarFieldSpecies::SetBackgroundInitialConditions`.

### 3.2. Fluid
- Move `background_w_fld` logic into `FluidSpecies`.
- Handle `EDE` (Early Dark Energy) specific checks within the species.

### 3.3. IDM_DRMD / IDR_DRMD
- Move `background_idm_drmd` logic into the composite class or children.

## 4. Verification Plan

1. **Build:** Ensure `make class` passes.
2. **Regression:** Run `pytest -v -m test_scenario python/test_class.py` (TEST_LEVEL=1).
3. **Consistency:** Compare `background_table` output for a complex scenario (e.g., `explanatory.ini`) before and after refactoring to ensure bit-wise (or near bit-wise) equality.
