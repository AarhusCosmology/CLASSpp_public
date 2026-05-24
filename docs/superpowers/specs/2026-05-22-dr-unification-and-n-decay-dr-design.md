# DR species unification + remove `pba->N_decay_dr` — design (PR B of #268 follow-ups)

## Goal

Finish the dark-radiation (DR) migration and delete `pba->N_decay_dr`. Concretely:

1. Collapse the legacy multi-channel `DarkRadiationSpecies` to a single self-contained DR channel, so
   `N_decay_dr` evaporates (each DR species is intrinsically one channel).
2. Unify the two DR implementations into **one generic `DarkRadiationSpecies`** class, deleting
   `DNCDM_DecayRadiationSpecies`.
3. Finish the perturbation-layout migration: remove the parallel `perturb_vector` DR fields and move
   the deferred `DCDM_DR` switch-copy / IC logic into the species hook.
4. Migrate the 3 idm_dr/idr/drmd dispatch guards inherited from PR C, via per-sub-species presence
   accessors on the IDM composites.

This is **not** behaviour-neutral at the bit level — collapsing the vestigial/"sum" slots shrinks the
perturbation vector, shifting downstream indices. Physics is unchanged; verification is the ~0.1%
tolerance scenario suite (never bit-identical; handle Cℓ zero-crossings), per project guidance.

## Background: where `N_decay_dr` lives today

`pba->N_decay_dr = (Omega0_dcdmdr>0?1:0) + n_dncdm` (set in `input_module.cpp:343`) — the total count
of DR-producing decay channels. 19 uses across `background.h` (the field), `input_module` (setter),
`background_module:701` (cached mirror `index_bg_rho_dr_`), `perturbations_module:3871/4933` (the
deferred DCDM_DR switch-copy + IC loops), `dcdm_dr_species` (output column loops), and
`dark_radiation_species.cpp` (10 — the legacy multi-channel layout/loops).

The DR subsystem is half-migrated:
- **`DNCDM_DecayRadiationSpecies`** (dncdm→dr) is already modern: each decaying ncdm owns its own
  `DNCDM_DR_Species` wrapping a single-hierarchy DR species (`idx_F0`, own `index_bg_rho_`/
  `index_bi_rho_`); it never touches `N_decay_dr`.
- **`DarkRadiationSpecies`** (dcdm→dr, inside `DCDM_DR_Species`) is the legacy holdout: although
  dcdm→dr is one channel, it lays out `N_decay_dr`×(l_max_dr+1) per-species slots + a sum hierarchy,
  evolves only channel 0 (`PerturbDerivs` has no per-channel loop — the higher channels are IC'd from
  zero-ρ and never evolved → dead weight), and still drives the parallel `perturb_vector` fields
  (`index_pt_F0_dr_sum/species`, `l_max_dr`).

**Verified equivalence of the two DR implementations:** background evolution is dilution-only
(`dρ = -4aH·ρ`) in both; the free-streaming Boltzmann hierarchy (l=0,1,2 with metric sources, l=3
special, l=4..lmax-1 generic, lmax truncation) and `Delta`/`Theta`/`DeltaP`/`RhoPlusPShear` are
byte-for-byte identical. The decay **source** is injected by the parent composite, not the DR species
(`DCDM_DR::BackgroundDerivs` adds `a·Γ_dcdm·ρ_dcdm`; `DNCDM_DR::BackgroundDerivs` adds `a·Γ·M·n`).
The *only* current difference is output: column naming + the synchronous-gauge correction in
`PrintVariables`.

## Architecture decisions

### 1. One generic single-channel `DarkRadiationSpecies`

`DarkRadiationSpecies` becomes a single-channel DR species (mirroring today's
`DNCDM_DecayRadiationSpecies`): one `idx_F0` hierarchy, one `index_bg_rho_`/`index_bi_rho_`, no sum,
no per-species array, no `N_decay_dr`. It holds all shared evolution (dilution, the Boltzmann
hierarchy, ICs, `Delta`/`Theta`/`DeltaP`/`RhoPlusPShear`, `Rho`/`P`/`DpDloga`).

`DNCDM_DecayRadiationSpecies` is **deleted**. Both `DCDM_DR_Species` and `DNCDM_DR_Species`
instantiate `DarkRadiationSpecies` (a composite needing N decay-radiation products simply owns N
instances). DR is **not shared** between decay channels — the old shared-pool optimization is dropped.

### 2. Generic synchronous-gauge output correction (enables the single class)

The hardcoded `decay_corr = a·Γ_dcdm·ρ_dcdm/ρ_dr` is exactly `ρ̇_dr/ρ_dr` minus dilution. Since the
parent composite already computes the full `ρ̇_dr` in `BackgroundDerivs` (dilution + decay source), the
output correction `(ρ̇_dr/ρ_dr)·α` is computed **generically** with no `Γ_dcdm` reference, so no
decay-channel-specific code remains in `DarkRadiationSpecies`.

**Recommended wiring (finalise in the plan):** the parent composite supplies the DR species with its
decay-source term the same way it injects it into `BackgroundDerivs` — e.g. a `rho_dot` provider /
small hook set at construction (defaulting to pure dilution `-4aH·ρ`). The DR species then forms
`(ρ̇_dr/ρ_dr)·α` itself. Alternatives considered: (a) the composite drives the DR perturbation output
and adds the source-correction (splits output across two classes); (b) cache `ρ̇_dr` in a background
slot during integration. Option chosen at plan time; all keep `DarkRadiationSpecies` channel-agnostic.

### 3. Output naming per instance

Column names derive from the DR species' instance name rather than hardcoded `delta_ncdm[n]`/`d_dr`.
The multi-channel `[n]` indexing in the current `PrintVariables` disappears (single channel). The
**total** `d_dr`/`t_dr` transfer (`index_tp_delta_dr_`/`index_tp_theta_dr_`) becomes a sum over the DR
species, consistent with the species-dispatch stress-energy loop (#261/#262) — the materialized "sum"
hierarchy is therefore removed.

### 4. Finish the layout migration

- Remove the parallel `perturb_vector` fields `index_pt_F0_dr_sum`, `index_pt_F0_dr_species`, and the
  `perturb_vector` copies of `l_max_dr`/`l_max_dr_col` (the precision `l_max_dr` config stays; the
  hierarchy length lives on the layout as `l_max`).
- Migrate the deferred `DCDM_DR` switch-copy block (`perturbations_module:3866`) and IC block
  (`:4931`) into `DarkRadiationSpecies::CopyPerturbationsAcrossSwitch` + the existing
  `ApplyInitialConditions` (the hook pattern PR A established for the other species). With single
  channels these become simple per-hierarchy copies; the two consecutive `if (count("DCDM_DR"))`
  blocks collapse out of the module.
- `background_module:701`'s `index_bg_rho_dr_` mirror and the `N_decay_dr` loops in
  `perturbations_module`/`dcdm_dr_species` go away or source their count from the (now singular) DR
  species.

### 5. Inherited-from-C IDM composite accessors (separable subsystem)

Add presence accessors capturing the `pba` flags at construction:
- `IDM_DR_IDR_Species`: `bool has_idm_dr() const`, `bool has_idr() const`
- `IDM_DRMD_IDR_DRMD_Species`: `bool has_idm_drmd() const`, `bool has_idr_drmd() const`

Then migrate the 3 guards off `pba->has_*`: `background_module:532`
(`has_idr_drmd && has_idm_drmd` → composite `.has_idr_drmd() && .has_idm_drmd()`), the
`nonlinear_module` `has_idm_dr` Halofit warning, and the `output_module` `has_idm_dr` thermo-file
header (→ the composite's `has_idm_dr()`). These were *not* swappable for `count("IDM_DR_IDR")` in
PR C because the composites are created on `has_idm_dr || has_idr`.

## Components touched

| File | Change |
|---|---|
| `species/dark_radiation_species.{h,cpp}` | Collapse to single channel; becomes the sole DR class holding all shared evolution; generic gauge correction; `CopyPerturbationsAcrossSwitch`; per-instance output naming |
| `species/dncdm_decay_radiation_species.{h,cpp}` | **Deleted** (remove from Makefile / setup.py / pbxproj) |
| `species/dncdm_dr_species.{h,cpp}` | Use `DarkRadiationSpecies` instead of `DNCDM_DecayRadiationSpecies` |
| `species/dcdm_dr_species.{h,cpp}` | Single-channel DR; supply decay-source for the gauge correction; drop `N_decay_dr` output loops |
| `source/background.h` | Remove `int N_decay_dr` |
| `source/input_module.cpp` | Remove the `N_decay_dr` setter (and its comment) |
| `source/background_module.cpp` | Remove the `index_bg_rho_dr_` mirror `+ N_decay_dr`; migrate guard `:532` to composite accessors |
| `source/perturbations_module.cpp` | Remove the `N_decay_dr` switch-copy/IC loops (now in the species hook); remove parallel-field plumbing |
| `source/perturbations.h` | Remove `index_pt_F0_dr_sum/species`, `l_max_dr`/`l_max_dr_col` perturb_vector fields |
| `species/idm_dr_idr_species.h`, `species/idm_drmd_idr_drmd_species.h` | Add presence accessors |
| `source/nonlinear_module.cpp`, `source/output_module.cpp` | Migrate the 2 `has_idm_dr` guards to composite accessors |

## Build systems

Deleting `dncdm_decay_radiation_species.{h,cpp}` requires updating **all three** build systems:
`Makefile` (SPECIES_OPP), `setup.py` (source_files), `CLASS.xcodeproj/project.pbxproj` (user verifies
Xcode). No new files are added.

## Suggested commit structure (separable chunks within the PR)

1. **Single-channel + unify DR**: collapse `DarkRadiationSpecies`, generic gauge correction, delete
   `DNCDM_DecayRadiationSpecies`, rewire both composites. Remove `N_decay_dr`.
2. **Finish layout migration**: parallel `perturb_vector` fields, switch-copy/IC hook, total-DR output.
3. **IDM composite accessors + 3 guards** (independent subsystem).

Chunks 1–2 are coupled (DR); chunk 3 is independent and could be split into its own PR if review gets
heavy.

## The shared-sum is buggy in multi-channel cases (verification consequence)

Tracing the legacy `index_pt_F0_dr_sum` end-to-end shows it is not merely vestigial — it
double-counts the decay source whenever there is more than the single dcdm channel:

- **dcdm-only:** `DarkRadiationSpecies` allocates a separate "sum" hierarchy alongside the
  "species[0]" hierarchy; both evolve identically and the collision source is added to both. The sum
  is a redundant copy of the one channel ⇒ the migration is **physically identical**.
- **dncdm-only:** there is no `DarkRadiationSpecies`, so `perturbations_module.cpp:3396` repoints
  `index_pt_F0_dr_sum` at the *first* dncdm's own `idx_F0`; `DNCDM_DR::AddCouplingDerivs`
  (`dncdm_dr_species.cpp:261-262`) then adds the collision term to **both** `base = my.dr.idx_F0` and
  `index_pt_F0_dr_sum`, which for the first/only dncdm is the **same slot** (source added twice).
- **combined dcdm+dncdm:** the stress-energy total double-counts each dncdm's collision (once via
  `DCDM_DR`'s sum-based DR, once via each `DNCDM_DR`'s own hierarchy).

This is a latent bug (combined was never tested; not guaranteed to have ever worked). The per-species
design removes the shared sum and is correct by construction, so dncdm/combined results **will
change** — matching the existing committed reference there would be *wrong*.

## Verification

- `make class -j` clean; `make classy`.
- **dcdm-only is the strict behavior-preserving baseline** — ~0.1% tolerance match (NOT bit-identical;
  zero-crossing-aware), incl. `test_dcdm_dr_matches_reference`.
- **dncdm-only and combined: correct-by-construction (option a)** — the migration fixes the
  double-count, so we do not match the old reference; sanity-check magnitude/shape (each DR species is
  standard free-streaming + its own decay source) and the author eyeballs `dncdm_dr` / combined
  output. (If ever used in a paper, test against known limits then.) Regenerate the `dncdm_dr` /
  combined references afterward.
- Full 84-scenario suite (`TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 -m test_scenario`) for no crashes.
- The `idm_dr`/`idr-only`/`drmd-partial` reference tests confirm the IDM guard migration is neutral.

## Open implementation details (resolve in the plan)

- Exact wiring of the generic `ρ̇_dr/ρ_dr` gauge correction (see decision 2).
- Per-instance DR output column names (what `DCDM_DR`'s and `DNCDM_DR`'s columns are called) and how
  the total `d_dr`/`t_dr` is summed across DR species.
- Whether `l_max_dr_col` (collision-term truncation) is still referenced anywhere after the
  perturb_vector copy is removed.
