# Rewrite `perturb_total_stress_energy` as a Species-Dispatch Loop

**Date:** 2026-04-22
**Owner:** Thomas Tram
**Status:** Design approved; plan pending.

## Problem

`PerturbationsModule::perturb_total_stress_energy` (`source/perturbations_module.cpp:~6090-6560`) is ~400 lines that should be roughly 20. The species-dispatch refactor replaced C pointer indirection with virtual calls but preserved the per-species block structure of the original C code. Composites already expose correct aggregate `Rho`/`Delta`/`Theta`/`DeltaP`/`RhoPlusPShear` (see `species/composite_species.cpp:61-163`), but the module bypasses them and reaches into composite children via `static_cast<XSpecies&>(*all_species_.at("X"))` (~10 sites in this function alone). The function also pulls fluid equation-of-state computation, PPF physics, IDR-under-TCA shear preamble, and DNCDM rescaling into itself — all of which are species-local physics that leaked out.

## Goals

1. Collapse `perturb_total_stress_energy` to a two-pass loop over `all_species_` where each pass is ~7 lines.
2. Move species-local physics into the species classes. After this PR, the only name-based reach-in remaining in the function is the one required to skip PPF fluid in the loop.
3. Preserve current numerical output on the full reference test suite (`TEST_LEVEL=2 COMPARE_OUTPUT_REF=1`).

## Non-goals

- Eliminating the other ~20 `static_cast<XSpecies&>(...)` downcast sites elsewhere in the module (background_module, thermodynamics_module, etc.). Follow-up PRs.
- Fixing the IDM_DR matter-tally asymmetry. Preserved in this PR via `IDM_DRSpecies::IsMatterSpecies() { return false; }` with a follow-up issue.
- Refactoring the `init()`/`free()` module lifecycle, raw `new`/`delete` in cosmology.cpp, or other parallel tidy-ups.
- Any change to `InputModule::ConstructSpecies` wiring.

## Design

### Architectural principle (reference)

`all_species_` holds gravity-only-coupled sectors. A sector is either a standalone `BaseSpecies` or a `CompositeSpecies` whose children couple to each other non-gravitationally. Module-side stress-energy accumulation is symmetric over sectors: each sector's `Rho`/`Delta`/`Theta`/`DeltaP`/`RhoPlusPShear` already exposes aggregate values, so the module needs no knowledge of sector internals.

### New virtuals on `BaseSpecies`

Five additions (and no new hook types beyond these):

```cpp
virtual bool IsMatterSpecies() const { return energy_type_ == EnergyType::Matter; }

virtual double MatterRho(const double* pvecback) const {
  return IsMatterSpecies() ? Rho(pvecback) : 0.;
}
virtual double MatterRhoDelta(const perturb_vector* pv, const double* y,
                              const double* pvecback, const perturb_workspace* ppw) const {
  return IsMatterSpecies() ? Rho(pvecback) * Delta(pv, y, pvecback, ppw) : 0.;
}
virtual double MatterRhoPlusPTheta(const perturb_vector* pv, const double* y,
                                   const double* pvecback, const perturb_workspace* ppw) const {
  return IsMatterSpecies() ? (Rho(pvecback) + P(pvecback)) * Theta(pv, y, pvecback, ppw) : 0.;
}
virtual double MatterRhoPlusP(const double* pvecback) const {
  return IsMatterSpecies() ? Rho(pvecback) + P(pvecback) : 0.;
}
```

Rationale: simple matter species inherit the defaults (no duplicated Rho/Delta logic). Composites with mixed matter+radiation children override all four to sum over matter children only. `IsMatterSpecies` is a per-sector gate for the loop.

### Overrides

- `IDM_DRSpecies::IsMatterSpecies() override { return false; }` — preserves current behavior where IDM_DR is excluded from the matter tally. One-line comment points to the follow-up issue for the asymmetry.
- `CompositeSpecies::IsMatterSpecies()`: returns true iff any child is matter.
- `CompositeSpecies::MatterRho` / `MatterRhoDelta` / `MatterRhoPlusPTheta` / `MatterRhoPlusP`: sum over children. Default implementation delegates to each child's corresponding `Matter*` method — which, for a non-matter child, returns 0 (so the sum automatically drops radiation).

No `PrepareStressEnergy` hook. Species compute what they need inside their own methods.

### Species-side relocations

#### `FluidSpecies` owns the fluid equation of state

- Move `BackgroundModule::background_w_fld` into `FluidSpecies::ComputeWFld(a, &w, &dw_over_da, &integral)`. `BackgroundModule` call sites (during background integration) dispatch via `fluid_species_->ComputeWFld(...)`.
- `FluidSpecies::Delta` / `Theta` / `DeltaP`: implement non-PPF fluid math using `pvecback[index_bg_w_fld_]` / `pvecback[index_bg_dw_over_da_fld_]` (indices already owned by the species). Currently `DeltaP` returns 0; after the move it returns the correctly gauge-transformed δp (the body currently at `perturbations_module.cpp:6444-6446`).
- New method `FluidSpecies::ComputePpf(k, a, a_prime_over_a, y, ppw)`: holds the PPF body (the `use_ppf == _TRUE_` branch, currently at `perturbations_module.cpp:6448+`). The module calls it explicitly after the main loop.

#### `IDRSpecies` owns its TCA shear

- Private `const` helper `double IDRSpecies::TcaShearIdr(pv, y, ppw) const` evaluating the TCA formula (currently at `perturbations_module.cpp:6165-6167`).
- `IDRSpecies::RhoPlusPShear` and `IDM_DR_IDR_Species::PrintVariables` (`species/idm_dr_idr_species.cpp:152`) both call `TcaShearIdr` when the shear register is absent.
- Delete `ppw->tca_shear_idm_dr` field; delete the corresponding preamble block in the module.
- New hook `BaseSpecies::SetThermodynamicsModule(const ThermodynamicsModule*)` (mirrors existing `SetBackgroundModule`). `IDRSpecies` stores the pointer to read `pvecthermo[index_th_dmu_idm_dr_]`.

#### `NCDMSpecies` / `DNCDMSpecies`: no `ncdm_id` exposure, no ppw arrays

- `NCDMSpecies::Delta` / `Theta` / `RhoPlusPShear`: return raw values from `y[]` (today's unrescaled path).
- `DNCDMSpecies::Delta` / `Theta` / `RhoPlusPShear` override: return rescaled values. `RescaledNCDMPerturbations` (currently a method on `PerturbationsModule`) moves into `DNCDMSpecies` as a private helper.
- Delete `ppw->delta_ncdm` / `ppw->theta_ncdm` / `ppw->shear_ncdm` arrays. All consumers update to call `species->Delta(...)` etc. directly.
- `ncdm_id` as a `BaseSpecies`-reachable concept goes away; if it survives at all it's a private member of `NCDMSpecies` used for generating output column titles.
- Delete `dynamic_cast<NCDMSpecies*>` / `dynamic_cast<DNCDMSpecies*>` inside `perturbations_module.cpp:6372-6376` and anywhere else.
- `ncdm_species_sorted_` cache is no longer needed in the stress-energy path. Retain it only if a separate path still uses it; otherwise delete.

### The rewritten function

```cpp
int PerturbationsModule::perturb_total_stress_energy(
    double k, double tau, double* y, perturb_workspace* ppw) {
  const double a = /* ... */;
  const double a_prime_over_a = /* ... */;
  const double k2 = k * k;

  // (a) Approximation-dependent preamble for scalar_ctx variables that feed
  //     species (photons shear_g under TCA, baryon δp). Unchanged from today.
  compute_scalar_ctx_preamble(k, k2, a, a_prime_over_a, y, ppw);

  // (b) Pass 1: shear (needed for ψ in Newtonian gauge).
  ppw->rho_plus_p_shear = 0.;
  for (auto& sp : all_species_)
    ppw->rho_plus_p_shear += sp->RhoPlusPShear(ppw->pv, y, ppw->pvecback, ppw);

  // (c) Derive ψ from total shear (Newtonian); no-op synchronous. ScalarField
  //     reads it via scalar_ctx inside its own Delta method.
  ppw->scalar_ctx.psi = derive_psi_from_shear(ppw);

  // (d) Pass 2: everything else.
  ppw->delta_rho = ppw->rho_plus_p_theta = ppw->delta_p = ppw->rho_plus_p_tot = 0.;
  double delta_rho_m = 0., rho_m = 0., rho_plus_p_theta_m = 0., rho_plus_p_m = 0.;

  for (auto& sp : all_species_) {
    if (sp->name() == "Fluid" && pba->use_ppf == _TRUE_) continue;  // PPF: handled below

    const double rho = sp->Rho(ppw->pvecback);
    const double rho_plus_p = rho + sp->P(ppw->pvecback);
    ppw->delta_rho        += rho        * sp->Delta(ppw->pv, y, ppw->pvecback, ppw);
    ppw->rho_plus_p_theta += rho_plus_p * sp->Theta(ppw->pv, y, ppw->pvecback, ppw);
    ppw->delta_p          += sp->DeltaP(ppw->pv, y, ppw->pvecback, ppw);
    ppw->rho_plus_p_tot   += rho_plus_p;

    if (sp->IsMatterSpecies()) {
      delta_rho_m        += sp->MatterRhoDelta(ppw->pv, y, ppw->pvecback, ppw);
      rho_m              += sp->MatterRho(ppw->pvecback);
      rho_plus_p_theta_m += sp->MatterRhoPlusPTheta(ppw->pv, y, ppw->pvecback, ppw);
      rho_plus_p_m       += sp->MatterRhoPlusP(ppw->pvecback);
    }
  }

  // (e) PPF fluid: one species allowed, depends on everything else; kept
  //     special in the module, body lives in FluidSpecies::ComputePpf.
  if (all_species_.count("Fluid") && pba->use_ppf == _TRUE_) {
    static_cast<FluidSpecies&>(*all_species_.at("Fluid"))
        .ComputePpf(k, a, a_prime_over_a, y, ppw);
  }

  // (f) cb / m wrap-up.
  if (has_source_delta_m_ && has_source_delta_cb_)
    ppw->delta_cb = delta_rho_m / rho_m;
  if ((has_source_delta_m_ || has_source_theta_m_) &&
      (has_source_delta_cb_ || has_source_theta_cb_))
    ppw->theta_cb = rho_plus_p_theta_m / rho_plus_p_m;

  return _SUCCESS_;
}
```

The `static_cast<FluidSpecies&>` at (e) is the only downcast that remains in this function — an explicit PPF special case.

## Staging plan

Each commit builds cleanly and passes `./class explanatory.ini`. Full reference suite runs at the end; bit-identical output is the target, fallback to existing tolerance if order-of-accumulation roundoff drift appears.

1. **Add `IsMatterSpecies` + four `Matter*` virtuals to `BaseSpecies` with defaults.** Add `CompositeSpecies` overrides that sum over matter children. Add `IDM_DRSpecies::IsMatterSpecies() { return false; }`. No behavior change; new virtuals unused.
2. **Move `background_w_fld` into `FluidSpecies::ComputeWFld`.** Update `BackgroundModule` call sites to dispatch. No behavior change.
3. **Move PPF body into `FluidSpecies::ComputePpf`. Move non-PPF fluid math into `FluidSpecies::Delta`/`Theta`/`DeltaP`.** Module still calls them from the existing fluid block. No behavior change.
4. **Move `RescaledNCDMPerturbations` into `DNCDMSpecies`; override `Delta`/`Theta`/`RhoPlusPShear`. Update all consumers of `ppw->delta_ncdm/theta_ncdm/shear_ncdm` to call species methods; delete the arrays.** Delete `ncdm_id` exposure outside `NCDMSpecies`. Largest commit.
5. **Move IDR TCA shear into `IDRSpecies::TcaShearIdr` helper. Add `BaseSpecies::SetThermodynamicsModule`. Delete `ppw->tca_shear_idm_dr`.** Update the two consumers.
6. **Rewrite `perturb_total_stress_energy` as two-pass loop + PPF call.** Delete per-species blocks at 6247-6425 and NCDM loop at 6360-6409. The 400-line function collapses here.
7. **Run `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` reference suite.** Document any non-bit-identical drift in the PR description.

Each commit is independently bisect-friendly.

## Verification

- `make clean && make class -j` clean build.
- `./class explanatory.ini` succeeds after every commit.
- `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py` after commits 1-6.
- `cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py` at the end. Target: bit-identical. Acceptable: within existing reference tolerance, drift noted in PR.
- All three build systems updated if new files added: `Makefile`, `setup.py`, `CLASS.xcodeproj/project.pbxproj`.

## Risks

- **`RescaledNCDMPerturbations` body may reach module-private state.** If it does, commit 4 grows to port that state onto `DNCDMSpecies` or `NCDMBaseSpecies`.
- **`ppw->delta_ncdm[n]` consumer count unknown.** Likely sites: `perturb_print_variables_member`, `perturb_sources_member`, `perturb_vector_init`. If consumer count exceeds ~5 sites, mid-PR check-in before continuing commit 4.
- **ψ ordering assumption.** Current code derives ψ implicitly via `perturb_einstein`. The two-pass design assumes no species reads ψ *during* shear accumulation. Verify before committing to two passes; if violated, design moves to three passes or recombines.
- **Order of accumulation changing last-bit roundoff.** Not a correctness concern; may require a follow-up tolerance update in the reference comparison.

## Follow-ups (explicitly not in this PR)

- IDM_DR matter-tally asymmetry (open issue from commit 1 comment).
- Remaining `static_cast<XSpecies&>(...)` downcasts in `background_module`, `thermodynamics_module`, `input_module`.
- `init()` / `free()` lifecycle removal across modules.
- `cosmology.cpp` raw `new` → `std::make_unique`.
