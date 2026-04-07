# Species Refactoring: CompositeSpecies and Uniform Dispatch

**Date:** 2026-04-06  
**Context:** Finalising the species-as-plugins refactoring started in PRs #194, #198, #211, #215.  
**Goal:** The main modules (background, perturbations, thermodynamics) must never check `has_X` to dispatch species logic. They loop over `all_species_` uniformly. Anything species-dependent lives in the species class.

---

## 1. Architecture Overview

Three independent but coordinated changes:

1. **`CompositeSpecies`** — a new `BaseSpecies` subclass that owns N child species and acts as the single `all_species_` entry for a physically coupled sector. The composite mediates coupling: children write non-interacting (free-streaming) terms into `dy`; the composite then writes coupling terms via `AddCouplingDerivs`.

2. **NCDM into `all_species_`** — each NCDM flavor becomes an individual `"NCDM_0"`, `"NCDM_1"`, … entry in `all_species_`. The `ncdm_species_` vector and the ~30 loops that use it are removed.

3. **Uniform species loop** — all `has_`-guarded per-species calls in the main modules are replaced by a single `for (auto& [name, sp] : all_species_)` loop. The `std::map` alphabetical order is acceptable; only PPF (FluidSpecies) requires a deferred second pass.

Three concrete composites are created:
- `DCDM_DR_Species`
- `IDM_DR_IDR_Species`
- `IDM_DRMD_IDR_DRMD_Species`

---

## 2. `CompositeSpecies` Class

```cpp
class CompositeSpecies : public BaseSpecies {
public:
  CompositeSpecies(std::string name, EnergyType energy_type);

  // Registration: delegates to children in insertion order
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void RegisterPerturbationIndices(perturb_vector*, const precision*,
                                   int& index_pt, const perturb_workspace*, int gauge) override;

  // Background: sums over children
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  // Perturbation derivs: two-phase (children free-streaming, then coupling)
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  // Stress-energy: weighted sum so Rho()*Delta() is correct for Einstein equations
  double Delta(const perturb_vector*, const double* y,
               const double* pvecback, const perturb_workspace*) const override;
  double Theta(...) const override;
  double DeltaP(...) const override;
  double RhoPlusPShear(...) const override;

protected:
  std::vector<std::unique_ptr<BaseSpecies>> children_;

  // Subclasses override to add coupling terms after children have run
  virtual void AddCouplingDerivs(double tau, const double* y, double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {}
};
```

**Key invariants:**
- `PerturbDerivs` calls each child's `PerturbDerivs`, then calls `AddCouplingDerivs`.
- `Delta()` returns `Σ(ρᵢ·δᵢ) / Σρᵢ` so that `Rho() * Delta()` gives the correct `δρ_total` for the Einstein equations.
- Children are `unique_ptr`-owned by the composite and never appear directly in `all_species_`.
- `EnergyType::Other` is used for composites spanning both matter and radiation (e.g. DCDM+DR). The background accumulation loop already handles `Other` as `rho_m += rho - 3p; rho_r += 3p`.

---

## 3. NCDM Migration

**In `InputModule`:**
```cpp
for (int n = 0; n < pba->N_ncdm; n++) {
  std::string name = "NCDM_" + std::to_string(n);
  all_species_[name] = std::make_unique<NCDMSpecies>(n, ncdm_, pba, nullptr);
}
```

- `ncdm_species_` vector in `BaseModule` is removed.
- Each `NCDMSpecies` continues to own its `ncdm_id_` and look up properties from the shared `NonColdDarkMatter` object.
- The ~30 `for (auto* ncdm_sp : ncdm_species_)` loops in `perturbations_module.cpp` and elsewhere collapse into the generic `for (auto& [name, sp] : all_species_)` loop.

**Rationale:** NCDM flavors are physically independent — they have no coupling to each other. Wrapping them in a composite would be wrong. They belong as individual entries in `all_species_` from the start; `ncdm_species_` is a historical workaround.

---

## 4. Uniform Species Loop

Replace all `has_`-guarded per-species calls in the main modules with:

```cpp
// Index registration
for (auto& [name, sp] : all_species_)
  sp->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, gauge);

// Derivs — first pass (all non-deferred species)
for (auto& [name, sp] : all_species_)
  if (!sp->RequiresDeferredPerturbDerivs())
    sp->PerturbDerivs(tau, y, dy, ppaw);

// Derivs — deferred pass (PPF fluid must go last)
for (auto& [name, sp] : all_species_)
  if (sp->RequiresDeferredPerturbDerivs())
    sp->PerturbDerivs(tau, y, dy, ppaw);
```

`BaseSpecies::RequiresDeferredPerturbDerivs()` defaults to `false`. `FluidSpecies` overrides it to `true`.

`FluidSpecies` also requires special treatment in the background `ComputeBackground` loop (it currently has a `if (name == "Fluid") continue;` skip because it needs `w_fld` setup first). A `virtual bool RequiresDeferredBackground() const` method handles this analogously.

**`has_` flags:** Remain in the `background` struct for validation/error-checking (e.g. the synchronous gauge CDM requirement check) but are no longer used to gate species method calls anywhere.

**`std::map` order:** Alphabetical order is acceptable. Species derivs only read from `y`, so call order is irrelevant for correctness. PPF is the sole exception, handled by the deferred pass above.

---

## 5. Concrete Composite Implementations

### `DCDM_DR_Species`

**Children:** `DCDMSpecies`, `DarkRadiationSpecies`.

**Background coupling:** `DCDM_DR_Species` overrides `BackgroundDerivs`, calls `CompositeSpecies::BackgroundDerivs` (which delegates to children), then adds the DR source from DCDM decay:
```cpp
dy[bi_dr_species_ + 0] += a * Gamma_dcdm * pvecback[index_bg_rho_dcdm];
```
No separate virtual `AddCouplingBackgroundDerivs` hook is needed — concrete composites simply override `BackgroundDerivs` and call the base first. This replaces the `has_dcdm` guard currently inside `DarkRadiationSpecies::BackgroundDerivs`.

**`AddCouplingDerivs`:** DR free-streaming hierarchy is written by `DarkRadiationSpecies::PerturbDerivs` without DCDM-sourced terms. The composite then adds:
```cpp
dy[base_dr + 0] += rprime_dr * (y[pv->index_pt_delta_dcdm] + metric_euler/k2);
dy[base_dr + 1] += rprime_dr / k * y[pv->index_pt_theta_dcdm];
```
DCDM has no coupling terms flowing back from DR; `DCDMSpecies::PerturbDerivs` requires no change.

**Joint input:** `Omega0_dcdmdr` (combined DCDM+DR energy density) is resolved in `InputModule` when constructing the composite, splitting it between children.

---

### `IDM_DR_IDR_Species`

**Children:** `IDM_DRSpecies`, `IDRSpecies`.

**`AddCouplingDerivs`:** Momentum-transfer terms that currently live inside each child behind `has_idm_dr` guards:
```cpp
// IDR receives momentum from IDM_DR
dy[pv->index_pt_theta_idr] +=
    dmu_idm_dr * (y[pv->index_pt_theta_idm_dr] - y[pv->index_pt_theta_idr]);

// IDM_DR receives back-reaction
dy[pv->index_pt_theta_idm_dr] -=
    R_idr * dmu_idm_dr * (y[pv->index_pt_theta_idm_dr] - y[pv->index_pt_theta_idr]);
```
TCA slip and shear terms (currently in `PerturbScalarContext`) also move here.

**`PerturbScalarContext` fields removed:** `delta_idr`, `theta_idr`, `theta_idm_dr`, `R_idr`, `tca_shear_idm_dr`, `tca_slip_idm_dr`. The composite reads directly from `y` using children's stored `index_pt_*` members.

---

### `IDM_DRMD_IDR_DRMD_Species`

Same pattern as `IDM_DR_IDR_Species`. Coupling terms from both children's `PerturbDerivs` that check `has_idr_drmd`/`has_idm_drmd` migrate to `AddCouplingDerivs`. Fields `delta_idr_drmd`, `theta_idr_drmd`, `delta_idm_drmd`, `theta_idm_drmd` are removed from `PerturbScalarContext`.

---

## 6. `PerturbScalarContext` After Refactoring

**Removed** (owned by composites):
- `delta_idr`, `theta_idr`, `shear_idr`
- `delta_idm_dr`, `theta_idm_dr`
- `delta_idr_drmd`, `theta_idr_drmd`
- `delta_idm_drmd`, `theta_idm_drmd`
- `R_idr`, `tca_shear_idm_dr`, `tca_slip_idm_dr`
- `idr_nature` (owned by `IDRSpecies` / `IDM_DR_IDR_Species`)

**Kept** (used across multiple unrelated species):
- `k`, `k2`, `a`, `a2`, `a_prime_over_a`
- `metric_continuity`, `metric_euler`, `metric_shear`, `metric_ufa_class`
- `cotKgen`, `s2_squared`
- `delta_g`, `theta_g`, `shear_g` (photons: read by baryons for TCA)
- `delta_b`, `theta_b` (baryons: read by photons for TCA)
- `R`, `cb2`, `delta_p_b_over_rho_b` (photon-baryon TCA)
- `gauge`

---

## 7. GitHub Workflow

1. Create a GitHub issue listing the remaining work (this spec).
2. Work in branch `{issue_number}-composite-species`.
3. Open a PR when ready, referencing the issue.
