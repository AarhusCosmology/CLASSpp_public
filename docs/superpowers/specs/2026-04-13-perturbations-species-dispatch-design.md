# Perturbations Species Dispatch — Design Spec

**Date:** 2026-04-13  
**Repo:** AarhusCosmology/CLASSpp  
**Branch target:** master (incremental PRs per stage)

---

## Goal

Extract species-specific logic from `perturbations_module.cpp` into the species classes under `species/`. The four dispatch zones are handled in three ordered stages, each tested before the next begins.

This round deliberately does **not** move `index_tp_*` ownership to species — that is deferred to a future round because it also requires simultaneous changes to the transfer, spectra, and Python wrapper modules.

---

## Stages

| Stage | Functions refactored | Test |
|-------|---------------------|------|
| 1 | `perturb_output_titles`, `perturb_output_data`, `perturb_print_variables_member` | `k_output_values` output + verbose run |
| 2 | `perturb_sources_member` | `TEST_LEVEL=0` pytest |
| 3 | `perturb_initial_conditions` | `TEST_LEVEL=1` pytest |

---

## Architecture

### Index ownership

`PerturbationsModule` retains all `index_tp_*` members and their assignment in `perturb_indices_of_perturbs()`. Species read these back through the module reference passed in context. Index delegation to species is out of scope for this round.

### Module call sites

Each refactored function gains a species dispatch loop over `all_species_` (a `std::map`, iterated alphabetically). Species-specific `if (pba->has_XXX)` blocks for individual species are removed and replaced by the loop. Aggregate/metric sources (`index_tp_phi_`, `index_tp_delta_tot_`, temperature source terms, etc.) remain in the module.

---

## New Infrastructure

### `PerturbColumnWriter` (Stage 1)

Thin helper in `species/perturb_source_context.h`. Constructed in either title-writing mode or data-writing mode. Species call `writer.Add(title, tp_index, active)` once per column — eliminating the current fragility where `class_store_columntitle` and `class_store_double` must be manually kept in sync across two separate functions.

```cpp
class PerturbColumnWriter {
public:
  explicit PerturbColumnWriter(char* titles);
  PerturbColumnWriter(double* dataptr, const double* tk, int& storeidx);

  // Call once per column, in the same order for titles and data.
  void Add(const char* title, int tp_index, bool active);

private:
  char*         titles_   = nullptr;
  double*       dataptr_  = nullptr;
  const double* tk_       = nullptr;
  int*          storeidx_ = nullptr;
};
```

### `PerturbSourceContext` (Stage 2)

POD struct in `species/perturb_source_context.h`. Bundles the workspace and addressing info needed by `FillSources`. Species write to the source table via `p_mod->SetSourceValue(...)`.

```cpp
struct PerturbSourceContext {
  PerturbationsModule* p_mod;
  perturb_workspace*   ppw;
  int index_md, index_ic, index_k, index_tau;
  double k;
  double a, a_prime_over_a;
};
```

`SetSourceValue` is a public accessor on `PerturbationsModule`:
```cpp
void SetSourceValue(int index_md, int index_ic, int index_tp,
                    int index_tau, int index_k, double value);
```

### `PerturbIcContext` (Stage 3)

POD struct in `species/perturb_source_context.h`. Bundles all shared IC ratios pre-computed by the module. Species write directly into `y[]` (`ppw->pv->y`).

```cpp
struct PerturbIcContext {
  // Shared ratios pre-computed by the module
  double fracnu, fracg, fracb, fraccdm, fracidm_drmd;
  double rho_m_over_rho_r, om;
  double ktau_two, ktau_three;
  double s2_squared;          // curvature factor: 1 - 3K/k^2
  // Pre-computed relativistic perturbations (needed for coupled ICs)
  double delta_ur, theta_ur, shear_ur, l3_ur, delta_dr, eta;
  double alpha, alpha_prime;  // synchronous gauge slip
  // Kinematics
  double k, tau, a, a_prime_over_a;
  int index_ic;               // which IC type (ad, cdi, bi, nid, niv, ten)
  int gauge;
  perturb_workspace* ppw;     // for writing ppw->pv->y[index_pt_*]
};
```

---

## `BaseSpecies` Interface Additions

All methods have empty default implementations. No existing species break at the point of adding the interface; species are converted one by one within each stage.

```cpp
// Stage 1 — output
virtual void WriteOutputColumns(PerturbColumnWriter& writer,
                                 const PerturbationsModule& mod,
                                 enum file_format fmt,
                                 BaseSpecies::TransferColumnSection section = BaseSpecies::TransferColumnSection::all) const {}

virtual void PrintVariables(PerturbColumnWriter& writer,
                             double tau, const double* y,
                             const PerturbationsModule& mod,
                             const perturb_workspace* ppw) const {}

// Stage 2 — source filling
virtual void FillSources(const double* y, const double* dy,
                          PerturbSourceContext& ctx) {}

// Stage 3 — initial conditions
virtual void ApplyInitialConditions(double* y,
                                     const PerturbIcContext& ctx) {}
```

---

## Species Covered

All species in `all_species_`: Photons, Baryons, CDM, UltraRelativistic, DarkRadiation/DCDM composite, FluidSpecies, ScalarField, NCDMSpecies, IDM-DR/IDR composite, IDM-DRMD/IDR-DRMD composite, DNCDM/DecayRadiation composite.

Lambda has no perturbation variables and needs no implementation.

---

## What Stays in the Module

- Temperature and polarization source terms (`index_tp_t0_`, `index_tp_t1_`, `index_tp_t2_`, `index_tp_p_`)
- Metric source terms (`index_tp_phi_`, `index_tp_psi_`, `index_tp_h_`, `index_tp_eta_`, etc.)
- Aggregate sources (`index_tp_delta_tot_`, `index_tp_delta_m_`, `index_tp_theta_tot_`, etc.)
- All `index_tp_*` assignment in `perturb_indices_of_perturbs()`
- All `has_source_*` flag assignment
- Shared IC physics: ratio computation (`fracnu`, `om`, etc.), neutrino hierarchy values, gauge transformation

Species `FillSources()` implementations only own species-local transfer entries
such as `delta_i`/`theta_i`. Module-owned CMB temperature/polarization, metric,
and aggregate `index_tp_*` slots remain written exclusively in
`perturbations_module.cpp`.

---

## Testing Protocol

### Stage 1
Run with `k_output_values = '0.01, 0.1'` and compare ASCII perturbation output against a reference CLASS run (column names and values must match). `perturb_print_variables_member` is triggered automatically by the same `k_output_values` condition, so both functions are exercised in the same run.

### Stage 2
```bash
TEST_LEVEL=0 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py
```
All scenarios must pass before proceeding to Stage 3.

### Stage 3
```bash
TEST_LEVEL=1 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py
```

---

## Out of Scope

- `index_tp_*` ownership delegation to species (next round, requires transfer + spectra + Python wrapper changes)
- `perturb_vector_init` / perturbation index registration (already done in prior PRs)
- Non-scalar modes for species that only have scalar perturbations
