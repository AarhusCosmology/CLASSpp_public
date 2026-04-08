# DNCDM Composite Species Design

**Date:** 2026-04-08
**Context:** Extracting Decaying Non-Cold Dark Matter (DNCDM) from `NCDMSpecies` and migrating it to the `CompositeSpecies` pattern to support a 1-to-1 relationship between DNCDM flavors and their decay radiation, enabling future momentum-dependent decay products.

## 1. Architecture Overview

1. **`DNCDMSpecies`**: Extract the `decay_dr` type from `NCDMSpecies` into a dedicated `DNCDMSpecies` class. `NCDMSpecies` will revert to handling only standard, non-decaying NCDM.
2. **`DNCDM_DecayRadiationSpecies`**: Create a new species class specifically for the decay products of a *single* DNCDM flavor. This replaces the shared `DarkRadiationSpecies` for DNCDM decay products, making it trivial to swap in a momentum-dependent DR model later.
3. **`DNCDM_DR_Species` Composite**: For each DNCDM flavor requested, instantiate a `DNCDM_DR_Species` composite (e.g., `DNCDM_DR_0`, `DNCDM_DR_1`). The composite owns exactly one `DNCDMSpecies` and one `DNCDM_DecayRadiationSpecies`. The composite's `BackgroundDerivs` and `AddCouplingDerivs` handle the energy and momentum transfer.

## 2. DNCDMSpecies

- Inherits from `BaseSpecies`.
- Handles the background distribution function `lnf` and its momentum derivatives `dlnfdlnq`, which are specific to decaying NCDM.
- In `PerturbDerivs`, handles the free-streaming Boltzmann hierarchy for the massive decaying particle.

## 3. DNCDM_DecayRadiationSpecies

- Inherits from `BaseSpecies`.
- Initially implements a standard, momentum-independent massless radiation hierarchy (`F0`).
- Provides a strict 1-to-1 mapping for a single DNCDM flavor's decay product.

## 4. DNCDM_DR_Species Composite

- Inherits from `CompositeSpecies`.
- **Background Coupling**: In `BackgroundDerivs`, adds the decay source term `+ a * Gamma * M_1 * n_dncdm` to the `DNCDM_DecayRadiationSpecies` energy density, removing energy from the `DNCDMSpecies` via the distribution function evolution.
- **Perturbation Coupling**: In `AddCouplingDerivs`, adds the collision terms from the DNCDM hierarchy into the `DNCDM_DecayRadiationSpecies` multipoles (`l=0` and `l=1`).

## 5. Main Modules Cleanup

- **`InputModule`**: Create individual `DNCDM_DR_Species` composites instead of setting `has_ncdm_decay_dr` and appending them to `ncdm_species_`.
- **`BackgroundModule` / `PerturbationsModule`**: Remove `has_ncdm_decay_dr` explicit checks. Since DNCDM is now a composite, the uniform `all_species_` loop introduced in PR #217 will automatically handle index registration and ODE derivatives.
- Remove contiguous `index_bg_lnf_ncdm_decay_dr1_` and `index_bi_rho_dr_species_` array assumptions for DNCDM in `background_module.cpp`. Each composite will own its indices.