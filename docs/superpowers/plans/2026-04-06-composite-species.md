# CompositeSpecies and Uniform Species Dispatch — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finalise the species-as-plugins refactoring by introducing `CompositeSpecies` for coupled sectors (DCDM+DR, IDM_DR+IDR, IDM_DRMD+IDR_DRMD), migrating NCDM to individual `all_species_` entries, and eliminating all `has_`-guarded per-species dispatch from the main modules.

**Architecture:** `CompositeSpecies : BaseSpecies` owns child species and mediates coupling — children write free-streaming terms into `dy`, the composite's `AddCouplingDerivs` writes coupling terms. All modules iterate over `all_species_` uniformly; `ncdm_species_` vector is removed.

**Tech Stack:** C++17, make build system (`make class`), Python test suite (`python/test_class.py`).

**Spec:** `docs/superpowers/specs/2026-04-06-species-composite-design.md`

---

## File Map

| Action | File | Purpose |
|--------|------|---------|
| Create | `species/composite_species.h` | `CompositeSpecies` base class |
| Create | `species/composite_species.cpp` | `CompositeSpecies` implementations |
| Create | `species/dcdm_dr_species.h` | `DCDM_DR_Species` composite |
| Create | `species/dcdm_dr_species.cpp` | `DCDM_DR_Species` implementations |
| Create | `species/idm_dr_idr_species.h` | `IDM_DR_IDR_Species` composite |
| Create | `species/idm_dr_idr_species.cpp` | `IDM_DR_IDR_Species` implementations |
| Create | `species/idm_drmd_idr_drmd_species.h` | `IDM_DRMD_IDR_DRMD_Species` composite |
| Create | `species/idm_drmd_idr_drmd_species.cpp` | `IDM_DRMD_IDR_DRMD_Species` implementations |
| Modify | `species/base_species.h` | Add `RequiresDeferredPerturbDerivs()`, `RequiresDeferredBackground()` |
| Modify | `species/perturb_context.h` | Remove coupling fields owned by composites |
| Modify | `species/dark_radiation_species.cpp` | Remove `has_dcdm` coupling in `BackgroundDerivs` and `PerturbDerivs` |
| Modify | `species/dcdm.cpp` | No coupling terms to remove (DCDM does not read DR) |
| Modify | `species/interacting_species.cpp` | Remove `has_idm_dr`/`has_idr_drmd`/`has_idm_drmd` coupling terms from IDR, IDM_DR, IDR_DRMD, IDM_DRMD `PerturbDerivs` |
| Modify | `source/input_module.cpp` | Register composites; keep NCDM per-flavour registration (already done) |
| Modify | `source/base_module.h` | Remove `ncdm_species_` vector and `PopulateNCDMVector()` |
| Modify | `source/background_module.h` | Update index members for composite-owned species |
| Modify | `source/background_module.cpp` | Replace `has_`-guarded dispatch with uniform loop; update index capture from composites |
| Modify | `source/perturbations_module.cpp` | Replace `has_`-guards + `ncdm_species_` loops with uniform loop |
| Modify | `source/thermodynamics_module.cpp` | Replace `has_`-guards with uniform loop |
| Modify | `Makefile` | Add new `.opp` files to `SPECIES_OPP` |

---

## Task 1: Create GitHub Issue and Branch

**Files:** none (git/GitHub ops only)

- [ ] **Step 1: Create the GitHub issue**

```bash
gh issue create \
  --title "Finalise species extraction: CompositeSpecies, uniform dispatch, has_-guard removal" \
  --body "$(cat docs/superpowers/specs/2026-04-06-species-composite-design.md)"
```

Note the issue number printed (e.g. `216`). Use it in the next step.

- [ ] **Step 2: Create the branch**

```bash
git checkout -b 216-composite-species   # replace 216 with actual issue number
```

---

## Task 2: Add `RequiresDeferredPerturbDerivs` and `RequiresDeferredBackground` to BaseSpecies

`FluidSpecies` (PPF dark energy) must run its perturbation derivs after all other species because it reads the accumulated stress-energy. It is also skipped in the background `ComputeBackground` loop because it needs `w_fld` setup first. Two virtual flags handle this.

**Files:**
- Modify: `species/base_species.h`

- [ ] **Step 1: Add the two virtual methods to `BaseSpecies`**

In `species/base_species.h`, inside the `public:` section after `DpDloga`, add:

```cpp
  /**
   * Returns true if this species' PerturbDerivs must run AFTER all other
   * species in a second pass. Used for PPF fluid (FluidSpecies).
   */
  virtual bool RequiresDeferredPerturbDerivs() const { return false; }

  /**
   * Returns true if this species' ComputeBackground must be deferred.
   * Used for FluidSpecies which needs w_fld evaluated before it can run.
   */
  virtual bool RequiresDeferredBackground() const { return false; }
```

- [ ] **Step 2: Override in `FluidSpecies`**

In `species/fluid.h`, add to the public section:

```cpp
  bool RequiresDeferredPerturbDerivs() const override { return true; }
  bool RequiresDeferredBackground() const override { return true; }
```

- [ ] **Step 3: Build to verify no compilation errors**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly (no new errors).

- [ ] **Step 4: Commit**

```bash
git add species/base_species.h species/fluid.h
git commit -m "Add RequiresDeferredPerturbDerivs/Background virtual flags to BaseSpecies"
```

---

## Task 3: Implement `CompositeSpecies` Base Class

**Files:**
- Create: `species/composite_species.h`
- Create: `species/composite_species.cpp`
- Modify: `Makefile`

- [ ] **Step 1: Write `species/composite_species.h`**

```cpp
#pragma once
#include "base_species.h"
#include <vector>
#include <memory>

/**
 * CompositeSpecies: a BaseSpecies that owns N child species and acts as a
 * single all_species_ entry for a physically coupled sector.
 *
 * Background methods sum over children.
 * PerturbDerivs runs a two-phase dispatch:
 *   1. Each child's PerturbDerivs (free-streaming terms)
 *   2. AddCouplingDerivs (coupling terms — override in concrete subclasses)
 *
 * Delta/Theta/DeltaP/RhoPlusPShear return the rho-weighted sum so that
 *   Rho() * Delta() == sum_i(rho_i * delta_i)
 * which is what the Einstein equations need.
 */
class CompositeSpecies : public BaseSpecies {
public:
  CompositeSpecies(std::string name, EnergyType energy_type)
    : BaseSpecies(std::move(name), energy_type) {}

  // ── Registration ────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void RegisterVectorPerturbationIndices(perturb_vector* pv, int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;
  void RegisterTensorPerturbationIndices(perturb_vector* pv, int& index_pt,
                                          const perturb_workspace* ppw,
                                          int gauge) override;

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void ComputeBackground(double a_rel, const double* pvecback_B,
                         double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy,
                        const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  // ── Perturbations ────────────────────────────────────────────────────────
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void PerturbVectorDerivs(double tau, const double* y, double* dy,
                            const perturb_parameters_and_workspace& ppaw) override;
  void PerturbTensorDerivs(double tau, const double* y, double* dy,
                            const perturb_parameters_and_workspace& ppaw) override;

  // Rho-weighted sums: Rho()*Delta() == sum_i(rho_i*delta_i)
  double Delta(const perturb_vector* pv, const double* y,
               const double* pvecback, const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv, const double* y,
               const double* pvecback, const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv, const double* y,
                const double* pvecback, const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv, const double* y,
                       const double* pvecback, const perturb_workspace* ppw) const override;

protected:
  std::vector<std::unique_ptr<BaseSpecies>> children_;

  /**
   * Override in concrete subclasses to add coupling terms to dy after
   * all children have written their free-streaming contributions.
   * Default: no-op.
   */
  virtual void AddCouplingDerivs(double tau, const double* y, double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {}
};
```

- [ ] **Step 2: Write `species/composite_species.cpp`**

```cpp
#include "composite_species.h"

void CompositeSpecies::RegisterBackgroundIndices(int& index_bg) {
  for (auto& child : children_)
    child->RegisterBackgroundIndices(index_bg);
}

void CompositeSpecies::RegisterIntegrationIndices(int& index_bi) {
  for (auto& child : children_)
    child->RegisterIntegrationIndices(index_bi);
}

void CompositeSpecies::RegisterPerturbationIndices(
    perturb_vector* pv, const precision* ppr, int& index_pt,
    const perturb_workspace* ppw, int gauge) {
  for (auto& child : children_)
    child->RegisterPerturbationIndices(pv, ppr, index_pt, ppw, gauge);
}

void CompositeSpecies::RegisterVectorPerturbationIndices(
    perturb_vector* pv, int& index_pt, const perturb_workspace* ppw, int gauge) {
  for (auto& child : children_)
    child->RegisterVectorPerturbationIndices(pv, index_pt, ppw, gauge);
}

void CompositeSpecies::RegisterTensorPerturbationIndices(
    perturb_vector* pv, int& index_pt, const perturb_workspace* ppw, int gauge) {
  for (auto& child : children_)
    child->RegisterTensorPerturbationIndices(pv, index_pt, ppw, gauge);
}

void CompositeSpecies::SetBackgroundModule(const BackgroundModule* bgm) {
  for (auto& child : children_)
    child->SetBackgroundModule(bgm);
}

void CompositeSpecies::ComputeBackground(double a_rel,
                                          const double* pvecback_B,
                                          double* pvecback) {
  for (auto& child : children_)
    child->ComputeBackground(a_rel, pvecback_B, pvecback);
}

void CompositeSpecies::BackgroundDerivs(double tau, const double* y,
                                         double* dy, const double* pvecback) {
  for (auto& child : children_)
    child->BackgroundDerivs(tau, y, dy, pvecback);
}

double CompositeSpecies::Rho(const double* pvecback) const {
  double rho = 0.;
  for (const auto& child : children_) rho += child->Rho(pvecback);
  return rho;
}

double CompositeSpecies::P(const double* pvecback) const {
  double p = 0.;
  for (const auto& child : children_) p += child->P(pvecback);
  return p;
}

double CompositeSpecies::DpDloga(const double* pvecback) const {
  double dp = 0.;
  for (const auto& child : children_) dp += child->DpDloga(pvecback);
  return dp;
}

void CompositeSpecies::PerturbDerivs(double tau, const double* y, double* dy,
                                      const perturb_parameters_and_workspace& ppaw) {
  for (auto& child : children_)
    child->PerturbDerivs(tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);
}

void CompositeSpecies::PerturbVectorDerivs(double tau, const double* y, double* dy,
                                            const perturb_parameters_and_workspace& ppaw) {
  for (auto& child : children_)
    child->PerturbVectorDerivs(tau, y, dy, ppaw);
}

void CompositeSpecies::PerturbTensorDerivs(double tau, const double* y, double* dy,
                                            const perturb_parameters_and_workspace& ppaw) {
  for (auto& child : children_)
    child->PerturbTensorDerivs(tau, y, dy, ppaw);
}

double CompositeSpecies::Delta(const perturb_vector* pv, const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const {
  double rho_delta = 0., rho_total = 0.;
  for (const auto& child : children_) {
    const double rho = child->Rho(pvecback);
    rho_delta += rho * child->Delta(pv, y, pvecback, ppw);
    rho_total += rho;
  }
  return (rho_total > 0.) ? rho_delta / rho_total : 0.;
}

double CompositeSpecies::Theta(const perturb_vector* pv, const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const {
  double rho_theta = 0., rho_total = 0.;
  for (const auto& child : children_) {
    const double rho = child->Rho(pvecback);
    rho_theta += rho * child->Theta(pv, y, pvecback, ppw);
    rho_total += rho;
  }
  return (rho_total > 0.) ? rho_theta / rho_total : 0.;
}

double CompositeSpecies::DeltaP(const perturb_vector* pv, const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  double dp = 0.;
  for (const auto& child : children_)
    dp += child->DeltaP(pv, y, pvecback, ppw);
  return dp;
}

double CompositeSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  double s = 0.;
  for (const auto& child : children_)
    s += child->RhoPlusPShear(pv, y, pvecback, ppw);
  return s;
}
```

- [ ] **Step 3: Add `composite_species.opp` to `Makefile`**

In `Makefile`, find the line:
```
SPECIES_OPP = cdm.opp photons.opp baryons.opp lambda.opp ultra_relativistic.opp fluid.opp dcdm.opp dark_radiation_species.opp ncdm_species.opp scalar_field.opp interacting_species.opp
```
Change it to:
```
SPECIES_OPP = cdm.opp photons.opp baryons.opp lambda.opp ultra_relativistic.opp fluid.opp dcdm.opp dark_radiation_species.opp ncdm_species.opp scalar_field.opp interacting_species.opp composite_species.opp
```

- [ ] **Step 4: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly.

- [ ] **Step 5: Commit**

```bash
git add species/composite_species.h species/composite_species.cpp Makefile
git commit -m "Add CompositeSpecies base class with two-phase PerturbDerivs dispatch"
```

---

## Task 4: Create `DCDM_DR_Species` Composite

**Context:** `DCDMSpecies::PerturbDerivs` has no coupling terms pointing to DR (matter doesn't respond to radiation at this level). `DarkRadiationSpecies::PerturbDerivs` currently reads `delta_dcdm`/`theta_dcdm` from `y[pv->index_pt_delta_dcdm]` behind a `has_dcdm` guard — these are the coupling terms. `DarkRadiationSpecies::BackgroundDerivs` similarly has a `has_dcdm` source term.

The plan:
1. Strip the DCDM-dependent code from `DarkRadiationSpecies::PerturbDerivs` and `BackgroundDerivs`.
2. Create `DCDM_DR_Species` composite that adds those terms in `AddCouplingDerivs` and overrides `BackgroundDerivs`.

**Files:**
- Modify: `species/dark_radiation_species.cpp`
- Create: `species/dcdm_dr_species.h`
- Create: `species/dcdm_dr_species.cpp`
- Modify: `Makefile`

- [ ] **Step 1: Strip the DCDM-sourced block from `DarkRadiationSpecies::PerturbDerivs`**

In `species/dark_radiation_species.cpp`, in `DarkRadiationSpecies::PerturbDerivs`, the entire block:
```cpp
  int index_dr = 0;

  // DCDM-sourced DR
  if (pba_->has_dcdm == _TRUE_) {
    ...
    ++index_dr;
  }
```
…currently computes the full DR Boltzmann hierarchy multiplied by `r_dr` and adds `rprime_dr * (delta_dcdm + ...)` / `rprime_dr / k * theta_dcdm` coupling sources. Replace the block so that `delta_dcdm` and `theta_dcdm` are always zero (the free-streaming part), i.e. remove the `rprime_dr` terms:

```cpp
  int index_dr = 0;

  // Free-streaming DR hierarchy (coupling source from DCDM added by DCDM_DR_Species::AddCouplingDerivs)
  if (pba_->N_decay_dr > 0) {
    const int base = pv->index_pt_F0_dr_species + index_dr * (pv->l_max_dr + 1);
    double r_dr = pvecback[bgm->index_bg_rho_dr_species_]
                  * (a * a) * (a * a) / (pba_->H0 * pba_->H0);

    // l=0: no rprime_dr source here
    dy[base + 0] = -k * y[base + 1] - 4./3. * metric_continuity * r_dr;
    // l=1: no rprime_dr source here
    dy[base + 1] = k/3. * y[base + 0]
                   - 2./3. * k * y[base + 2] * s2_squared
                   + 4. * metric_euler / (3. * k) * r_dr;
    // l=2
    dy[base + 2] = 8./15. * (3./4. * k * y[base + 1] + metric_shear * r_dr)
                   - 3./5. * k * s_l[3] / s_l[2] * y[base + 3];
    // l=3
    {
      int l = 3;
      dy[base + l] = k / (2.*l + 1.) * (l * s_l[l] * s_l[2] * y[base + l - 1]
                                          - (l + 1.) * s_l[l + 1] * y[base + l + 1]);
    }
    // l=4..l_max_dr-1
    for (int l = 4; l < pv->l_max_dr; ++l)
      dy[base + l] = k / (2.*l + 1.) * (l * s_l[l] * y[base + l - 1]
                                          - (l + 1.) * s_l[l + 1] * y[base + l + 1]);
    // l=l_max_dr
    {
      int l = pv->l_max_dr;
      dy[base + l] = k * (s_l[l] * y[base + l - 1] - (1. + l) * cotKgen * y[base + l]);
    }
    // Accumulate into sum
    for (int l = 0; l <= pv->l_max_dr; ++l)
      dy[pv->index_pt_F0_dr_sum + l] += dy[base + l];

    ++index_dr;
  }
```

Also remove the `if (!pba_->has_dr) return;` early-return guard at the top of the function — the composite only calls this when DR is present.

- [ ] **Step 2: Strip the DCDM source from `DarkRadiationSpecies::BackgroundDerivs`**

In `species/dark_radiation_species.cpp`, in `BackgroundDerivs`, remove the `has_dcdm` branch:

```cpp
  for (int n = 0; n < pba_->N_decay_dr; ++n) {
    // Dilution only; DCDM decay source added by DCDM_DR_Species::BackgroundDerivs override
    dy[index_bi_rho_dr_species_ + n] =
        -4. * a * H * y[index_bi_rho_dr_species_ + n];
  }
```

- [ ] **Step 3: Write `species/dcdm_dr_species.h`**

```cpp
#pragma once
#include "composite_species.h"
#include "dcdm.h"
#include "dark_radiation_species.h"
#include "background.h"

class BackgroundModule;

/**
 * DCDM_DR_Species: composite for Decaying Cold Dark Matter + its decay radiation.
 *
 * EnergyType::Other so the background loop splits rho correctly:
 *   rho_m += rho - 3p  (= rho_dcdm)
 *   rho_r += 3p        (= rho_dr)
 *
 * BackgroundDerivs override adds the DCDM→DR decay source on top of children.
 * AddCouplingDerivs adds the DR Boltzmann l=0,1 coupling source terms from DCDM.
 */
class DCDM_DR_Species : public CompositeSpecies {
public:
  DCDM_DR_Species(std::shared_ptr<DarkRadiation> dr,
                  const background* pba,
                  const BackgroundModule* bgm);

  void SetBackgroundModule(const BackgroundModule* bgm) override;

  // Override to add DCDM→DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy,
                        const double* pvecback) override;

  // Typed accessors so BackgroundModule can capture the child indices
  DCDMSpecies&             dcdm() { return *dcdm_; }
  DarkRadiationSpecies&    dr()   { return *dr_sp_; }
  const DCDMSpecies&       dcdm() const { return *dcdm_; }
  const DarkRadiationSpecies& dr() const { return *dr_sp_; }

protected:
  void AddCouplingDerivs(double tau, const double* y, double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

private:
  DCDMSpecies*          dcdm_ = nullptr;   // non-owning pointer into children_
  DarkRadiationSpecies* dr_sp_ = nullptr;  // non-owning pointer into children_
  const background*     pba_;
  const BackgroundModule* bgm_ = nullptr;
};
```

- [ ] **Step 4: Write `species/dcdm_dr_species.cpp`**

```cpp
#include "dcdm_dr_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include <cmath>

DCDM_DR_Species::DCDM_DR_Species(std::shared_ptr<DarkRadiation> dr,
                                   const background* pba,
                                   const BackgroundModule* bgm)
  : CompositeSpecies("DCDM_DR", BaseSpecies::EnergyType::Other)
  , pba_(pba), bgm_(bgm)
{
  auto dcdm = std::make_unique<DCDMSpecies>(*pba);
  auto dr_sp = std::make_unique<DarkRadiationSpecies>(dr, pba, bgm);
  dcdm_ = dcdm.get();
  dr_sp_ = dr_sp.get();
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(dr_sp));
}

void DCDM_DR_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

void DCDM_DR_Species::BackgroundDerivs(double tau, const double* y,
                                        double* dy, const double* pvecback) {
  // Children handle their own dilution terms
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);

  // DCDM→DR decay source (first DR channel)
  const double a = pvecback[bgm_->index_bg_a_];
  dy[dr_sp_->bi_rho_dr_species_index()] +=
      a * pba_->Gamma_dcdm * pvecback[dcdm_->bg_rho_index()];
}

void DCDM_DR_Species::AddCouplingDerivs(double tau, const double* y, double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  if (pv->index_pt_delta_dcdm < 0 || pv->index_pt_F0_dr_species < 0) return;

  auto* bgm = ppaw.perturbations_module->GetBackgroundModule().get();
  const double* pvecback = ppw->pvecback;
  const double a = ctx.a;
  const double k = ctx.k;

  const int base = pv->index_pt_F0_dr_species;  // first DR channel, index_dr=0
  const double r_dr = pvecback[bgm->index_bg_rho_dr_species_]
                      * (a * a) * (a * a) / (pba_->H0 * pba_->H0);
  const double rprime_dr = pba_->Gamma_dcdm
                           * pvecback[dcdm_->bg_rho_index()]
                           * std::pow(a, 5) / (pba_->H0 * pba_->H0);

  const double delta_dcdm = y[pv->index_pt_delta_dcdm];
  const double theta_dcdm = y[pv->index_pt_theta_dcdm];

  // Add DCDM source to DR l=0 and l=1
  dy[base + 0] += rprime_dr * (delta_dcdm + ctx.metric_euler / (k * k));
  dy[base + 1] += rprime_dr / k * theta_dcdm;

  // Keep sum slots consistent
  dy[pv->index_pt_F0_dr_sum + 0] += rprime_dr * (delta_dcdm + ctx.metric_euler / (k * k));
  dy[pv->index_pt_F0_dr_sum + 1] += rprime_dr / k * theta_dcdm;
}
```

- [ ] **Step 5: Add to Makefile**

In `Makefile`, update `SPECIES_OPP` to include `dcdm_dr_species.opp`:
```
SPECIES_OPP = cdm.opp photons.opp baryons.opp lambda.opp ultra_relativistic.opp fluid.opp dcdm.opp dark_radiation_species.opp ncdm_species.opp scalar_field.opp interacting_species.opp composite_species.opp dcdm_dr_species.opp
```

- [ ] **Step 6: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly.

- [ ] **Step 7: Commit**

```bash
git add species/dark_radiation_species.cpp species/dcdm_dr_species.h species/dcdm_dr_species.cpp Makefile
git commit -m "Add DCDM_DR_Species composite; strip DCDM coupling from DarkRadiationSpecies"
```

---

## Task 5: Create `IDM_DR_IDR_Species` Composite

**Context:** `IDM_DRSpecies::PerturbDerivs` reads `ctx.R_idr` and `ctx.theta_idr` (IDR state). `IDRSpecies::PerturbDerivs` reads `ctx.theta_idm_dr`, `ctx.delta_idm_dr`, `ctx.R_idr`, `ctx.tca_shear_idm_dr`, `ctx.tca_slip_idm_dr` (IDM_DR state). All these become direct reads from `y` in the composite.

The plan:
1. Strip all IDM_DR↔IDR coupling terms from both individual species.
2. Create the composite with `AddCouplingDerivs`.

**Files:**
- Modify: `species/interacting_species.cpp`
- Create: `species/idm_dr_idr_species.h`
- Create: `species/idm_dr_idr_species.cpp`
- Modify: `Makefile`

- [ ] **Step 1: Strip coupling from `IDM_DRSpecies::PerturbDerivs`**

In `species/interacting_species.cpp`, `IDM_DRSpecies::PerturbDerivs` currently has two branches depending on `tca_idm_dr_off`. Replace the full function body with just the free-streaming (non-TCA, non-coupling) part:

```cpp
void IDM_DRSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                   const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  // idm_dr density: same as CDM continuity equation
  dy[pv->index_pt_delta_idm_dr] = -(y[pv->index_pt_theta_idm_dr] + ctx.metric_continuity);

  // idm_dr velocity: Hubble friction only (coupling to IDR added by IDM_DR_IDR_Species)
  dy[pv->index_pt_theta_idm_dr] = -ctx.a_prime_over_a * y[pv->index_pt_theta_idm_dr]
                                  + ctx.metric_euler;
}
```

- [ ] **Step 2: Strip coupling from `IDRSpecies::PerturbDerivs`**

`IDRSpecies::PerturbDerivs` is long. Keep the structural skeleton (RSA check, delta continuity, free-streaming hierarchy) but remove all IDM_DR references.

Replace the function body in `species/interacting_species.cpp`:

```cpp
void IDRSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  if (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_on) return;

  auto* ppt = ppaw.perturbations_module->GetPerturbs();
  const double* s_l = ppw->s_l;
  const double k = ctx.k;
  const double cotKgen = ctx.cotKgen;

  // idr density
  dy[pv->index_pt_delta_idr] = -4./3.*(y[pv->index_pt_theta_idr] + ctx.metric_continuity);

  if (pba_.idr_nature == idr_free_streaming) {
    // idr velocity (free-streaming; IDM_DR coupling added by IDM_DR_IDR_Species)
    dy[pv->index_pt_theta_idr] = ctx.k2*(y[pv->index_pt_delta_idr]/4.
                                         - ctx.s2_squared*y[pv->index_pt_shear_idr])
                                 + ctx.metric_euler;

    // idr shear
    int l = 2;
    dy[pv->index_pt_shear_idr] = 0.5*(8./15.*(y[pv->index_pt_theta_idr] + ctx.metric_shear)
                                       - 3./5.*k*s_l[3]/s_l[2]*y[pv->index_pt_shear_idr+1]);

    // idr l=3
    l = 3;
    dy[pv->index_pt_l3_idr] = k/(2.*l+1.)*(l*2.*s_l[l]*s_l[2]*y[pv->index_pt_shear_idr]
                                             - (l+1.)*s_l[l+1]*y[pv->index_pt_l3_idr+1]);

    // idr l>3
    for (l = 4; l < pv->l_max_idr; l++)
      dy[pv->index_pt_delta_idr+l] = k/(2.*l+1)*(l*s_l[l]*y[pv->index_pt_delta_idr+l-1]
                                                   - (l+1.)*s_l[l+1]*y[pv->index_pt_delta_idr+l+1]);

    // idr lmax
    l = pv->l_max_idr;
    dy[pv->index_pt_delta_idr+l] = k*(s_l[l]*y[pv->index_pt_delta_idr+l-1]
                                       - (1.+l)*cotKgen*y[pv->index_pt_delta_idr+l]);
  } else {
    // Fluid idr velocity
    dy[pv->index_pt_theta_idr] = ctx.k2/4. * y[pv->index_pt_delta_idr] + ctx.metric_euler;
  }
}
```

Note: `pba_.idr_nature` does not exist yet — `idr_nature` is currently in `PerturbScalarContext`. Add `int idr_nature` as a member to `IDRSpecies` (stored from `pba_` at construction time, or read from `ppt` via `ppaw`). The simplest approach: keep `idr_nature` in `PerturbScalarContext` for now (it is not coupling state, it is a model parameter). Replace `pba_.idr_nature` above with `ctx.idr_nature` and leave `idr_nature` in `PerturbScalarContext`.

Also: remove the `alpha_idm_dr`/`beta_idr` Compton-scattering collision terms from all the hierarchy equations (they depend on `dmu_idm_dr` which is coupling state). These move to `AddCouplingDerivs`.

- [ ] **Step 3: Write `species/idm_dr_idr_species.h`**

```cpp
#pragma once
#include "composite_species.h"
#include "idm_dr.h"
#include "idr.h"
#include "background.h"

/**
 * IDM_DR_IDR_Species: composite for interacting dark matter + interacting dark radiation.
 *
 * Children handle free-streaming terms; this composite's AddCouplingDerivs
 * adds the momentum-exchange and TCA terms that couple IDM_DR to IDR.
 */
class IDM_DR_IDR_Species : public CompositeSpecies {
public:
  explicit IDM_DR_IDR_Species(const background& pba);

  IDM_DRSpecies& idm_dr() { return *idm_dr_; }
  IDRSpecies&    idr()    { return *idr_; }
  const IDM_DRSpecies& idm_dr() const { return *idm_dr_; }
  const IDRSpecies&    idr()    const { return *idr_; }

protected:
  void AddCouplingDerivs(double tau, const double* y, double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

private:
  IDM_DRSpecies* idm_dr_ = nullptr;
  IDRSpecies*    idr_    = nullptr;
  const background& pba_;
};
```

- [ ] **Step 4: Write `species/idm_dr_idr_species.cpp`**

```cpp
#include "idm_dr_idr_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

IDM_DR_IDR_Species::IDM_DR_IDR_Species(const background& pba)
  : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other)
  , pba_(pba)
{
  auto idm = std::make_unique<IDM_DRSpecies>(pba);
  auto idr = std::make_unique<IDRSpecies>(pba);
  idm_dr_ = idm.get();
  idr_    = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}

void IDM_DR_IDR_Species::AddCouplingDerivs(double tau, const double* y, double* dy,
                                            const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  if (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_on) return;

  auto* pth_mod = ppaw.perturbations_module->GetThermodynamicsModule().get();
  auto* bgm     = ppaw.perturbations_module->GetBackgroundModule().get();
  const double* pvecback  = ppw->pvecback;
  const double* pvecthermo = ppw->pvecthermo;
  auto* ppt = ppaw.perturbations_module->GetPerturbs();

  const double dmu_idm_dr = pvecthermo[pth_mod->index_th_dmu_idm_dr_];
  const double Sinv       = 4./3. * pvecback[idr_->bg_rho_index()]
                            / pvecback[idm_dr_->bg_rho_index()];  // R_idr

  const double theta_idm_dr = y[pv->index_pt_theta_idm_dr];
  const double theta_idr    = (pv->index_pt_theta_idr >= 0)
                              ? y[pv->index_pt_theta_idr] : 0.;

  if (ppw->approx[ppw->index_ap_tca_idm_dr] == (int)tca_idm_dr_off) {
    const double dmu_idr = pth_mod->GetThermodynamics()->b_idr
                         / pth_mod->GetThermodynamics()->a_idm_dr
                         * pba_.Omega0_idr / pba_.Omega0_idm_dr * dmu_idm_dr;
    const double* s_l = ppw->s_l;

    // IDM_DR velocity coupling
    dy[pv->index_pt_theta_idm_dr] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr)
                                    - ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_]
                                      * y[pv->index_pt_delta_idm_dr];

    // IDR velocity coupling
    if (ctx.idr_nature == idr_free_streaming) {
      dy[pv->index_pt_theta_idr] += dmu_idm_dr * (theta_idm_dr - theta_idr);

      // IDR Compton collision terms in the higher multipoles
      const int l_max = pv->l_max_idr;
      for (int l = 2; l <= l_max; l++) {
        dy[pv->index_pt_delta_idr + l] -=
            (ppt->alpha_idm_dr[l-2] * dmu_idm_dr + ppt->beta_idr[l-2] * dmu_idr)
            * y[pv->index_pt_delta_idr + l];
      }
    }
  } else {
    // TCA on: replace both velocity equations with joint TCA expressions
    const double tca_shear_idm_dr =
        0.5 * 8./15. / dmu_idm_dr / ppt->alpha_idm_dr[0]
        * (theta_idm_dr + ctx.metric_shear);
    const double tca_slip_idm_dr =
        (pth_mod->GetThermodynamics()->nindex_idm_dr - 2./(1.+Sinv))
        * ctx.a_prime_over_a * (theta_idm_dr - theta_idr)
        + 1./(1.+Sinv) / dmu_idm_dr
        * (-ctx.a_prime_over_a * theta_idm_dr
           + ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_]
             * y[pv->index_pt_delta_idm_dr]
           + ctx.k2 * Sinv * (y[pv->index_pt_delta_idr]/4. - tca_shear_idm_dr));

    dy[pv->index_pt_theta_idm_dr] =
        1./(1.+Sinv) * (-ctx.a_prime_over_a*theta_idm_dr
                        + ctx.k2*pvecthermo[pth_mod->index_th_cidm_dr2_]
                          *y[pv->index_pt_delta_idm_dr]
                        + ctx.k2*Sinv*(y[pv->index_pt_delta_idr]/4. - tca_shear_idm_dr))
        + ctx.metric_euler + Sinv/(1.+Sinv)*tca_slip_idm_dr;

    dy[pv->index_pt_theta_idr] =
        1./(1.+Sinv) * (-ctx.a_prime_over_a*theta_idm_dr
                        + ctx.k2*pvecthermo[pth_mod->index_th_cidm_dr2_]
                          *y[pv->index_pt_delta_idm_dr]
                        + ctx.k2*Sinv*(y[pv->index_pt_delta_idr]/4. - tca_shear_idm_dr))
        + ctx.metric_euler - 1./(1.+Sinv)*tca_slip_idm_dr;
  }
}
```

- [ ] **Step 5: Add to Makefile**

```
SPECIES_OPP = ... composite_species.opp dcdm_dr_species.opp idm_dr_idr_species.opp
```

- [ ] **Step 6: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly.

- [ ] **Step 7: Commit**

```bash
git add species/interacting_species.cpp species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp Makefile
git commit -m "Add IDM_DR_IDR_Species composite; strip IDM_DR/IDR coupling from individual species"
```

---

## Task 6: Create `IDM_DRMD_IDR_DRMD_Species` Composite

Same pattern as Task 5 for the DRMD model variant.

**Files:**
- Modify: `species/interacting_species.cpp` (IDM_DRMDSpecies, IDR_DRMDSpecies)
- Create: `species/idm_drmd_idr_drmd_species.h`
- Create: `species/idm_drmd_idr_drmd_species.cpp`
- Modify: `Makefile`

- [ ] **Step 1: Strip coupling from `IDM_DRMDSpecies::PerturbDerivs`**

In `interacting_species.cpp`, `IDM_DRMDSpecies::PerturbDerivs` currently checks `pba_.has_idr_drmd` and calls `background_idm_drmd`. Remove this block — leave only the free-streaming part:

```cpp
void IDM_DRMDSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                     const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[pv->index_pt_delta_idm_drmd] = -(y[pv->index_pt_theta_idm_drmd] + ctx.metric_continuity);
  // Free-streaming velocity: coupling to IDR_DRMD added by IDM_DRMD_IDR_DRMD_Species
  dy[pv->index_pt_theta_idm_drmd] = -ctx.a_prime_over_a * y[pv->index_pt_theta_idm_drmd]
                                    + ctx.metric_euler;
}
```

- [ ] **Step 2: Strip coupling from `IDR_DRMDSpecies::PerturbDerivs`**

Similarly remove the `has_idm_drmd` block:

```cpp
void IDR_DRMDSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                     const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[pv->index_pt_delta_idr_drmd] = -4./3.*(y[pv->index_pt_theta_idr_drmd] + ctx.metric_continuity);
  // Free-streaming velocity: coupling to IDM_DRMD added by IDM_DRMD_IDR_DRMD_Species
  dy[pv->index_pt_theta_idr_drmd] = 0.25 * ctx.k2 * y[pv->index_pt_delta_idr_drmd]
                                   + ctx.metric_euler;
}
```

- [ ] **Step 3: Write `species/idm_drmd_idr_drmd_species.h`**

```cpp
#pragma once
#include "composite_species.h"
#include "idm_drmd.h"
#include "idr_drmd.h"
#include "background.h"

class IDM_DRMD_IDR_DRMD_Species : public CompositeSpecies {
public:
  explicit IDM_DRMD_IDR_DRMD_Species(const background& pba);

  IDM_DRMDSpecies& idm_drmd() { return *idm_drmd_; }
  IDR_DRMDSpecies& idr_drmd() { return *idr_drmd_; }
  const IDM_DRMDSpecies& idm_drmd() const { return *idm_drmd_; }
  const IDR_DRMDSpecies& idr_drmd() const { return *idr_drmd_; }

protected:
  void AddCouplingDerivs(double tau, const double* y, double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

private:
  IDM_DRMDSpecies* idm_drmd_ = nullptr;
  IDR_DRMDSpecies* idr_drmd_ = nullptr;
  const background& pba_;
};
```

- [ ] **Step 4: Write `species/idm_drmd_idr_drmd_species.cpp`**

```cpp
#include "idm_drmd_idr_drmd_species.h"
#include "background_module.h"
#include "perturbations_module.h"

IDM_DRMD_IDR_DRMD_Species::IDM_DRMD_IDR_DRMD_Species(const background& pba)
  : CompositeSpecies("IDM_DRMD_IDR_DRMD", BaseSpecies::EnergyType::Other)
  , pba_(pba)
{
  auto idm = std::make_unique<IDM_DRMDSpecies>(pba);
  auto idr = std::make_unique<IDR_DRMDSpecies>(pba);
  idm_drmd_ = idm.get();
  idr_drmd_ = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}

void IDM_DRMD_IDR_DRMD_Species::AddCouplingDerivs(
    double tau, const double* y, double* dy,
    const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*   ppw = ppaw.ppw;
  const perturb_vector*      pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  auto* bgm     = ppaw.perturbations_module->GetBackgroundModule().get();
  const double* pvecback = ppw->pvecback;

  double Rint, csp2, Gint;
  class_call(bgm->background_idm_drmd(
      ctx.a,
      pvecback[idm_drmd_->bg_rho_index()] / pvecback[idr_drmd_->bg_rho_index()],
      &Rint, &csp2, &Gint),
    bgm->error_message_, bgm->error_message_);

  const double theta_idm = y[pv->index_pt_theta_idm_drmd];
  const double theta_idr = y[pv->index_pt_theta_idr_drmd];
  const double delta_idr = y[pv->index_pt_delta_idr_drmd];

  if (ppw->approx[ppw->index_ap_tca_idm_drmd] == (int)tca_idm_drmd_on) {
    double GdDelta_idm = 3.*csp2*(ctx.a_prime_over_a*theta_idm + ctx.k2*delta_idr/4.);
    dy[pv->index_pt_theta_idm_drmd] = 0.25*ctx.k2*delta_idr + ctx.metric_euler
                                     - GdDelta_idm * Rint;
    double GdDelta_idr = 3.*csp2*(ctx.a_prime_over_a*theta_idm + ctx.k2*delta_idr/4.);
    dy[pv->index_pt_theta_idr_drmd] -= GdDelta_idr * Rint;
  } else {
    dy[pv->index_pt_theta_idm_drmd] += Gint*(theta_idr - theta_idm);
    dy[pv->index_pt_theta_idr_drmd] -= Gint * Rint * (theta_idr - theta_idm);
  }
}
```

- [ ] **Step 5: Add to Makefile**

```
SPECIES_OPP = ... idm_dr_idr_species.opp idm_drmd_idr_drmd_species.opp
```

- [ ] **Step 6: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly.

- [ ] **Step 7: Commit**

```bash
git add species/interacting_species.cpp species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp Makefile
git commit -m "Add IDM_DRMD_IDR_DRMD_Species composite; strip DRMD coupling from individual species"
```

---

## Task 7: Wire Composites into `InputModule` and `BackgroundModule`

**Files:**
- Modify: `source/input_module.cpp`
- Modify: `source/background_module.cpp`
- Modify: `source/background_module.h`

- [ ] **Step 1: Update `InputModule::ConstructSpecies`**

In `source/input_module.cpp`, add the new composite headers at the top of `ConstructSpecies` (or in the includes section):

```cpp
#include "../species/dcdm_dr_species.h"
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
```

Replace the individual DCDM/DR registrations:

```cpp
  // OLD:
  // if (pba->has_dcdm == _TRUE_) all_species_["DCDM"] = ...;
  // if (pba->has_dr   == _TRUE_) all_species_["DR"]   = ...;
  //
  // NEW:
  if (pba->has_dcdm == _TRUE_ || pba->has_dr == _TRUE_) {
    all_species_["DCDM_DR"] = std::make_unique<DCDM_DR_Species>(dr_, pba, nullptr);
  }
```

Replace IDM_DR/IDR:

```cpp
  // OLD:
  // if (pba->has_idm_dr == _TRUE_) all_species_["IDM_DR"] = ...;
  // if (pba->has_idr    == _TRUE_) all_species_["IDR"]    = ...;
  //
  // NEW:
  if (pba->has_idm_dr == _TRUE_ || pba->has_idr == _TRUE_) {
    all_species_["IDM_DR_IDR"] = std::make_unique<IDM_DR_IDR_Species>(*pba);
  }
```

Replace IDM_DRMD/IDR_DRMD:

```cpp
  // OLD:
  // if (pba->has_idm_drmd == _TRUE_) all_species_["IDM_DRMD"] = ...;
  // if (pba->has_idr_drmd == _TRUE_) all_species_["IDR_DRMD"] = ...;
  //
  // NEW:
  if (pba->has_idm_drmd == _TRUE_ || pba->has_idr_drmd == _TRUE_) {
    all_species_["IDM_DRMD_IDR_DRMD"] = std::make_unique<IDM_DRMD_IDR_DRMD_Species>(*pba);
  }
```

- [ ] **Step 2: Update `BackgroundModule` to capture indices from composites**

In `source/background_module.h`, add accessors for the composite child indices. The members `index_bg_rho_dcdm_`, `index_bg_rho_dr_`, etc. must be set from the composite's children after registration.

In `source/background_module.cpp`, in `background_indices()`, replace the old individual species blocks:

```cpp
  // DCDM + DR composite
  if (pba->has_dcdm == _TRUE_ || pba->has_dr == _TRUE_) {
    auto& dcdm_dr = static_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
    all_species_.at("DCDM_DR")->RegisterBackgroundIndices(index_bg);
    index_bg_rho_dcdm_       = dcdm_dr.dcdm().bg_rho_index();
    index_bg_rho_dr_species_ = dcdm_dr.dr().bg_rho_dr_species_index();
    index_bg_rho_dr_         = dcdm_dr.dr().bg_rho_index();
  } else {
    index_bg_rho_dcdm_ = index_bg_rho_dr_ = index_bg_rho_dr_species_ = -1;
  }

  // IDM_DR + IDR composite
  if (pba->has_idm_dr == _TRUE_ || pba->has_idr == _TRUE_) {
    auto& idm_idr = static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
    all_species_.at("IDM_DR_IDR")->RegisterBackgroundIndices(index_bg);
    index_bg_rho_idm_dr_ = idm_idr.idm_dr().bg_rho_index();
    index_bg_rho_idr_    = idm_idr.idr().bg_rho_index();
  } else {
    index_bg_rho_idm_dr_ = index_bg_rho_idr_ = -1;
  }

  // IDM_DRMD + IDR_DRMD composite
  if (pba->has_idm_drmd == _TRUE_ || pba->has_idr_drmd == _TRUE_) {
    auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
    all_species_.at("IDM_DRMD_IDR_DRMD")->RegisterBackgroundIndices(index_bg);
    index_bg_rho_idm_drmd_ = drmd.idm_drmd().bg_rho_index();
    index_bg_rho_idr_drmd_ = drmd.idr_drmd().bg_rho_index();
  } else {
    index_bg_rho_idm_drmd_ = index_bg_rho_idr_drmd_ = -1;
  }
```

Do the same for the integration-index equivalents (`index_bi_rho_dcdm_`, `index_bi_rho_dr_species_`) using the child accessors `dcdm().bi_rho_index()` and `dr().bi_rho_dr_species_index()`.

Also add the composite includes at the top of `background_module.cpp`:
```cpp
#include "../species/dcdm_dr_species.h"
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
```

- [ ] **Step 3: Update `BackgroundModule::background_derivs_member`**

In `background_module.cpp`, in `background_derivs_member`, the existing code has a special `if (pba->has_dr) { ... }` block that handles the DR ODE inline (bypassing the species loop). Replace this with the composite's `BackgroundDerivs` which now handles both DCDM and DR. Ensure the generic `all_species_` loop no longer skips "DR" (since "DR" no longer exists as a standalone key).

Search for and remove:
```cpp
    if (name == "DR") continue;
```

- [ ] **Step 4: Build**

```bash
make class -j4 2>&1 | tail -30
```

Expected: compiles cleanly. Fix any linker/include errors.

- [ ] **Step 5: Run a quick sanity check**

```bash
./class explanatory.ini 2>&1 | tail -5
```

Expected: CLASS runs without crash for the default (no DCDM/IDR) configuration.

- [ ] **Step 6: Commit**

```bash
git add source/input_module.cpp source/background_module.cpp source/background_module.h
git commit -m "Wire DCDM_DR, IDM_DR_IDR, IDM_DRMD_IDR_DRMD composites into InputModule and BackgroundModule"
```

---

## Task 8: Replace `ncdm_species_` Loops in `background_module.cpp`

**Files:**
- Modify: `source/background_module.cpp`

`ncdm_species_` is still used in `background_module.cpp` for index registration and derivs loops. Replace each with a loop that iterates `all_species_` and dynamic-casts.

- [ ] **Step 1: Add a helper lambda at the top of each function that needs it**

Wherever `ncdm_species_` is looped, replace with:

```cpp
// Iterate all NCDM species in ncdm_id order
std::vector<NCDMSpecies*> ncdm_vec;
for (auto& [name, sp] : all_species_) {
  if (auto* n = dynamic_cast<NCDMSpecies*>(sp.get())) ncdm_vec.push_back(n);
}
std::sort(ncdm_vec.begin(), ncdm_vec.end(),
          [](NCDMSpecies* a, NCDMSpecies* b){ return a->ncdm_id() < b->ncdm_id(); });
```

Then use `ncdm_vec` wherever `ncdm_species_` was used.

- [ ] **Step 2: Update `background_indices()` NCDM block**

Replace:
```cpp
for (auto* ncdm : ncdm_species_) {
  ncdm->RegisterBackgroundIndices(index_bg);
}
if (!ncdm_species_.empty()) {
  index_bg_rho_ncdm1_      = ncdm_species_[0]->bg_rho_index();
  ...
}
```

With the lambda pattern above (ncdm_vec), then:
```cpp
for (auto* ncdm : ncdm_vec) {
  ncdm->RegisterBackgroundIndices(index_bg);
}
if (!ncdm_vec.empty()) {
  index_bg_rho_ncdm1_      = ncdm_vec[0]->bg_rho_index();
  index_bg_p_ncdm1_        = ncdm_vec[0]->bg_p_index();
  index_bg_pseudo_p_ncdm1_ = ncdm_vec[0]->bg_pseudo_p_index();
}
```

- [ ] **Step 3: Replace remaining `ncdm_species_` loops**

For each of the remaining ~5 occurrences in `background_module.cpp` (integration indices, background derivs for decay_dr, etc.), apply the same lambda pattern.

- [ ] **Step 4: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly.

- [ ] **Step 5: Commit**

```bash
git add source/background_module.cpp
git commit -m "Replace ncdm_species_ loops in background_module with all_species_ dynamic_cast iteration"
```

---

## Task 9: Replace `ncdm_species_` Loops in `perturbations_module.cpp`

There are ~30 occurrences. The pattern is the same lambda used in Task 8.

**Files:**
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 1: Add the NCDM helper lambda to each function that uses `ncdm_species_`**

For each function in `perturbations_module.cpp` that contains `for (auto* ncdm_sp : ncdm_species_)`, add the sorted lambda at the top of that function's relevant scope and use `ncdm_vec` in place of `ncdm_species_`.

Use this as a search-replace guide:
```bash
grep -n "ncdm_species_" source/perturbations_module.cpp
```

There are approximately these call sites:
- `perturb_indices` (index registration)
- `perturb_initial_conditions` (several loops)
- `perturb_vector_init` (registration)
- `perturb_derivs_member` (ODE RHS)
- `perturb_sources` (output)
- `perturb_store_source_functions`

Apply the lambda in each.

- [ ] **Step 2: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly.

- [ ] **Step 3: Run tests**

```bash
cd /Users/au192734/Projects/class_claude && python python/test_class.py -v 2>&1 | tail -30
```

Expected: all existing tests pass (NCDM behaviour unchanged).

- [ ] **Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "Replace ncdm_species_ loops in perturbations_module with all_species_ dynamic_cast iteration"
```

---

## Task 10: Remove `ncdm_species_` from `BaseModule`

**Files:**
- Modify: `source/base_module.h`
- Modify: `source/input_module.cpp` (remove `PopulateNCDMVector`)

- [ ] **Step 1: Remove `ncdm_species_` and `PopulateNCDMVector` from `base_module.h`**

Delete:
```cpp
  std::vector<NCDMSpecies*> ncdm_species_;
```
and:
```cpp
  void PopulateNCDMVector();
```
and the call `PopulateNCDMVector();` in the constructor.

- [ ] **Step 2: Remove `BaseModule::PopulateNCDMVector` definition from `input_module.cpp`**

Delete the entire `void BaseModule::PopulateNCDMVector() { ... }` function body (lines ~4112–4130).

- [ ] **Step 3: Build**

```bash
make class -j4 2>&1 | tail -20
```

Fix any remaining `ncdm_species_` references (there may be a few in other modules — find with `grep -rn ncdm_species_ source/`).

- [ ] **Step 4: Commit**

```bash
git add source/base_module.h source/input_module.cpp
git commit -m "Remove ncdm_species_ vector and PopulateNCDMVector from BaseModule"
```

---

## Task 11: Uniform Species Loop in `background_module.cpp`

Replace remaining `has_`-guarded calls with loops over `all_species_`.

**Files:**
- Modify: `source/background_module.cpp`

- [ ] **Step 1: Replace `background_indices()` with uniform loop**

After the fixed indices (a, H, H', photons, baryons), the remaining species index registration becomes:

```cpp
  for (auto& [name, sp] : all_species_) {
    // Fixed species already registered above; composites handle children internally
    if (name == "Photons" || name == "Baryons") continue;
    sp->RegisterBackgroundIndices(index_bg);
  }
```

Remove all the individual `if (pba->has_X)` blocks that call `.at("X")->RegisterBackgroundIndices(...)`. The `index_bg_rho_X_` captures from composites (done in Task 7) stay — they just happen after the composite's `RegisterBackgroundIndices` call inside the loop.

Note: because these captures require knowing which species is which composite, move the composite-specific index capture code into typed `SetBackgroundModule`/`RegisterBackgroundIndices` overrides in each concrete composite, or keep the explicit captures just after the generic loop in `background_indices()`. The latter (keep explicit captures after the loop) is simpler.

- [ ] **Step 2: Replace `ComputeBackground` dispatch with uniform two-pass loop**

Find the existing `ComputeBackground` loop (already partially uniform). Make it fully uniform:

```cpp
  // First pass: all non-deferred species
  for (const auto& [name, sp] : all_species_) {
    if (sp->RequiresDeferredBackground()) continue;
    sp->ComputeBackground(a_rel, pvecback_B, pvecback);
    accumulate(*sp);
  }
  // Second pass: deferred (FluidSpecies / PPF)
  for (const auto& [name, sp] : all_species_) {
    if (!sp->RequiresDeferredBackground()) continue;
    sp->ComputeBackground(a_rel, pvecback_B, pvecback);
    accumulate(*sp);
  }
```

Remove the old `if (name == "Fluid") continue;` skip.

- [ ] **Step 3: Replace `BackgroundDerivs` dispatch**

The existing `background_derivs_member` has explicit `has_idm_dr`, `has_idm_drmd` calls followed by the generic loop. Replace with the generic loop (composites now handle all coupling internally):

```cpp
  for (const auto& [name, sp] : all_species_) {
    sp->BackgroundDerivs(tau, y, dy, pvecback);
  }
```

Remove the explicit `if (pba->has_idm_dr)` / `if (pba->has_idm_drmd)` blocks and the inline DR ODE code (now in `DCDM_DR_Species::BackgroundDerivs`).

- [ ] **Step 4: Build and run**

```bash
make class -j4 2>&1 | tail -20
./class explanatory.ini 2>&1 | tail -5
```

Expected: clean build, CLASS runs without crash.

- [ ] **Step 5: Commit**

```bash
git add source/background_module.cpp
git commit -m "Uniform species loop in background_module: replace has_-guards with all_species_ iteration"
```

---

## Task 12: Uniform Species Loop in `perturbations_module.cpp`

Replace the remaining 342 `has_`-guarded calls. The most important are `RegisterPerturbationIndices` and `PerturbDerivs`.

**Files:**
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 1: Replace `RegisterPerturbationIndices` dispatch (scalar)**

Find the block in `perturb_indices` that explicitly calls each species' registration. Replace:

```cpp
  // Photons and Baryons first (special l_max setup needed)
  ppv->l_max_g     = ppr->l_max_g;
  ppv->l_max_pol_g = ppr->l_max_pol_g;
  all_species_.at("Photons")->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);
  all_species_.at("Baryons")->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);

  // All remaining species (composites handle children internally)
  for (auto& [name, sp] : all_species_) {
    if (name == "Photons" || name == "Baryons") continue;
    // Set any l_max fields needed before registration
    if (name == "UR") ppv->l_max_ur = ppr->l_max_ur;
    sp->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);
  }
```

The `ppv->l_max_dr`, `ppv->N_ncdm`, `ppv->l_max_ncdm`, etc. setups that currently precede the individual calls must move into the respective species' `RegisterPerturbationIndices` methods (or remain as pre-loop setup). Check each species' registration for any `ppv->l_max_X` setup that the main code currently does before calling the species — move those into the species or composite as needed.

- [ ] **Step 2: Replace `PerturbDerivs` dispatch**

Find `perturb_derivs_member` and the explicit `if (pba->has_X)` derivs calls. Replace with the two-pass loop:

```cpp
  for (auto& [name, sp] : all_species_) {
    if (!sp->RequiresDeferredPerturbDerivs())
      sp->PerturbDerivs(tau, y, dy, ppaw);
  }
  for (auto& [name, sp] : all_species_) {
    if (sp->RequiresDeferredPerturbDerivs())
      sp->PerturbDerivs(tau, y, dy, ppaw);
  }
```

- [ ] **Step 3: Replace remaining `has_`-guards for output/title/diagnostic code**

The ~60 remaining guards in `perturb_sources`, `perturb_store_source_functions`, column-title registration etc. are mostly of the form:

```cpp
class_store_columntitle(titles, "d_cdm", pba->has_cdm);
```

These are not species dispatch — they're output formatting. A clean approach: add a virtual `std::string perturbation_output_title() const` to `BaseSpecies` (default empty), override in each species to return the relevant title, and loop over `all_species_` to build the title list. However, this is a large change; if time is limited, leave output-formatting guards as-is (they're harmless — they just need the `has_` flags to stay in `background`). Mark these with `// TODO: migrate to species virtual method` comments.

- [ ] **Step 4: Build and test**

```bash
make class -j4 2>&1 | tail -20
python python/test_class.py -v 2>&1 | tail -30
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "Uniform species loop in perturbations_module: replace has_-guards and ncdm loops"
```

---

## Task 13: Uniform Loop in `thermodynamics_module.cpp`

**Files:**
- Modify: `source/thermodynamics_module.cpp`

The 18 `has_`-guards in thermodynamics mostly guard IDM_DR-specific recombination physics (drag terms, heat transfer). These cannot be replaced by a simple species loop because thermodynamics doesn't call `BaseSpecies` virtual methods — it does its own physics computation that happens to be conditional on species presence.

- [ ] **Step 1: Audit each `has_` guard**

```bash
grep -n "pba->has_\|pba_->has_" source/thermodynamics_module.cpp
```

For each occurrence, determine if it gates a species virtual method call (replaceable with a loop) or gates thermodynamics-internal physics (keep as-is but replace `pba->has_idm_dr` with `all_species_.count("IDM_DR_IDR") > 0`).

- [ ] **Step 2: Replace species-presence checks with `all_species_` existence checks**

Change guards like:
```cpp
if (pba->has_idm_dr == _TRUE_) { ... }
```
to:
```cpp
if (all_species_.count("IDM_DR_IDR") > 0) { ... }
```

This eliminates the dependency on `has_` flags while keeping the thermodynamics logic unchanged.

- [ ] **Step 3: Build and test**

```bash
make class -j4 2>&1 | tail -20
python python/test_class.py -v 2>&1 | tail -30
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add source/thermodynamics_module.cpp
git commit -m "Replace has_-guards in thermodynamics_module with all_species_.count() checks"
```

---

## Task 14: Clean Up `PerturbScalarContext`

**Files:**
- Modify: `species/perturb_context.h`
- Modify: `source/perturbations_module.cpp` (remove population of removed fields)

- [ ] **Step 1: Remove coupling fields from `PerturbScalarContext`**

In `species/perturb_context.h`, remove the following fields (now owned by composites):

```cpp
  // REMOVE these:
  double delta_idr = 0., theta_idr = 0., shear_idr = 0.;
  double delta_idm_dr = 0., theta_idm_dr = 0.;
  double delta_idr_drmd = 0., theta_idr_drmd = 0.;
  double delta_idm_drmd = 0., theta_idm_drmd = 0.;
  double R_idr = 0.;
  double tca_shear_idm_dr = 0.;
  double tca_slip_idm_dr = 0.;
```

Keep `idr_nature` (it is a model parameter read by IDRSpecies, not coupling state).

- [ ] **Step 2: Remove population of removed fields in `perturbations_module.cpp`**

Find the block in `perturb_derivs_member` (around line 6290) that sets these fields:
```cpp
ppw->scalar_ctx.delta_idr     = delta_idr;
ppw->scalar_ctx.theta_idr     = theta_idr;
ppw->scalar_ctx.delta_idm_dr  = y[ppw->pv->index_pt_delta_idm_dr];
...
ppw->scalar_ctx.tca_shear_idm_dr = ...;
ppw->scalar_ctx.tca_slip_idm_dr  = ...;
```

Delete all of these assignments.

- [ ] **Step 3: Build**

```bash
make class -j4 2>&1 | tail -20
```

Expected: compiles cleanly. If any species still reads a removed field, fix it to read from `y` directly.

- [ ] **Step 4: Run full test suite**

```bash
python python/test_class.py -v 2>&1 | tail -40
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add species/perturb_context.h source/perturbations_module.cpp
git commit -m "Remove coupling fields from PerturbScalarContext; composites read from y directly"
```

---

## Task 15: Integration Test with Interacting Species

Verify the composites produce correct physics by running CLASS with DCDM, IDM_DR+IDR, and NCDM inputs.

**Files:** none (testing only)

- [ ] **Step 1: Test DCDM+DR**

```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini \
  Omega_dcdmdr=0.12 Gamma_dcdm=10.0 output=tCl,pCl,lCl 2>&1 | tail -10
```

Expected: runs to completion, produces output files.

- [ ] **Step 2: Test IDM_DR+IDR**

```bash
./class explanatory.ini \
  Omega_idm_dr=0.12 xi_idr=0.3 a_idm_dr=1e-4 output=tCl 2>&1 | tail -10
```

Expected: runs to completion.

- [ ] **Step 3: Test NCDM**

```bash
./class explanatory.ini \
  N_ncdm=2 m_ncdm="0.06,0.03" deg_ncdm="2.0,1.0" output=tCl 2>&1 | tail -10
```

Expected: runs to completion, results match pre-refactoring output.

- [ ] **Step 4: Run full Python test suite**

```bash
python python/test_class.py -v 2>&1 | tee test_output.txt | tail -40
```

Expected: all tests pass. If any fail, investigate the delta in physics output.

---

## Task 16: Create Pull Request

**Files:** none (GitHub ops)

- [ ] **Step 1: Push the branch**

```bash
git push -u origin 216-composite-species   # use actual issue number
```

- [ ] **Step 2: Create the PR**

```bash
gh pr create \
  --title "CompositeSpecies for coupled sectors; uniform species dispatch; has_-guard removal" \
  --body "$(cat <<'EOF'
## Summary

- Introduces `CompositeSpecies : BaseSpecies` with two-phase `PerturbDerivs` (children write free-streaming terms; composite adds coupling terms via `AddCouplingDerivs`)
- Three concrete composites: `DCDM_DR_Species`, `IDM_DR_IDR_Species`, `IDM_DRMD_IDR_DRMD_Species`
- NCDM flavours were already in `all_species_` as individual entries; removes `ncdm_species_` vector and `PopulateNCDMVector`
- Background, perturbations and thermodynamics modules replaced `has_`-guarded per-species dispatch with uniform `for (auto& [name, sp] : all_species_)` loops
- Removes coupling fields (`delta_idr`, `theta_idm_dr`, etc.) from `PerturbScalarContext`

Closes #216

## Test plan

- [ ] `make class` builds cleanly with no warnings
- [ ] `python/test_class.py` passes all existing tests
- [ ] Manual run with `Omega_dcdmdr=0.12 Gamma_dcdm=10.0` produces valid Cl spectrum
- [ ] Manual run with IDM_DR+IDR inputs produces valid output
- [ ] Manual run with `N_ncdm=2` produces valid output
EOF
)"
```
