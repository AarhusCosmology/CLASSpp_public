# Encapsulate Species Perturbation Layouts (and Remove `ncdm_id`)

## Background

The earlier spec at `docs/superpowers/specs/2026-05-03-remove-ncdm-id-design.md` proposed moving NCDM perturbation indices from the flat per-`pv` arrays `pv->index_ncdm_[ncdm_id_]/l_max_ncdm[ncdm_id_]/q_size_ncdm[ncdm_id_]` onto member fields of `NCDMBaseSpecies` (`pt_l_max_/pt_q_size_/pt_index_per_q_`). That implementation crashed under multi-threaded evolution: `BaseSpecies` instances are shared across OpenMP threads, but per-thread `RegisterPerturbationIndices` calls each tried to overwrite the same instance's `pt_*` members. The result was a write-write-read race producing wrong indices, out-of-bounds writes into `y`, and either heap corruption or NaN derivatives that stalled the integrator.

The bug exposed a broader architectural problem. `perturb_vector` today has many hardcoded fields per species type (`index_pt_phi`, `index_pt_theta_b`, `index_pt_delta_g`, `index_pt_eta`, `l_max_g`, etc.). Each field assumes one instance per species type. Even setting aside `ncdm_id`, the codebase cannot today instantiate (say) two `ScalarFieldSpecies` instances side-by-side without clobbering shared `pv->index_pt_phi_scf` slots.

This spec replaces the previous one. The goal is broader: encapsulate every species's perturbation layout so each species owns its layout type, multi-instance support is automatic, and the integer `ncdm_id` (and its flat-array siblings) are gone for good.

## Goals

1. Each species defines its own polymorphic layout struct holding its perturbation indices.
2. `perturb_vector` holds a `std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>>` parallel to `all_species_`. Per-thread isolation is automatic because each thread allocates its own `pv` inside the OpenMP region.
3. All species's `Register*PerturbationIndices` populate their own layout slot. All read methods (`PerturbDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`, `CopyPerturbationsAcrossSwitch`, `MarkUsedInSources`) take the layout as argument.
4. Hardcoded per-species fields on `perturb_vector` (`index_pt_phi`, `index_pt_delta_g`, `index_ncdm_`, `l_max_ncdm`, etc.) are deleted. Only true non-species state (metric `eta`/`phi`, perturbed-recombination, integration storage) stays as bare fields.
5. The integer `ncdm_id` and its associated members (`pba->N_ncdm`, `Omega0_ncdm_tot`, `N_decay_dr`, `PerturbationsModule::ncdm_species_sorted_`, `pv->index_pt_psi0_ncdm1`) are deleted.
6. NCDM-family output column titles use instance names (`(.)rho_<instance>`, `delta_<instance>`). DNCDM_DR composite is renamed to `<dncdm_instance_name>` (and DR child to `<dncdm_instance_name>_DR`).
7. Shooter slot registration: each species owns its shooter slots via `RegisterShootingIndices`/`ComputeShootingResidual`/`ComputeShootingGuess` hooks. Cosmological-level `theta_s` shooting remains module-level.
8. Public API change: `PerturbationsModule::RescaledNCDMPerturbations(int n_ncdm, ...)` → `(BaseSpecies*, ...)`.

## Non-goals

- Generalising every multi-instance species through input parsing. The architectural support exists after this work but the input parser still accepts only one `ScalarField`, one `Fluid`, etc. Adding multi-instance input syntax is a follow-up.
- Migrating non-NCDM approximation-switch copy blocks (UR ufa, IDR rsa/tca, IDM_DRMD tca) into the new `CopyPerturbationsAcrossSwitch` pattern. The hook is general; the wholesale migration is a follow-up.
- `pba->has_ncdm` / `pba->has_ncdm_decay_dr` has-guards — opportunistic cleanup track.
- Rewriting `Omega0_ncdm_tot` consumers (5 sites) to avoid the `dynamic_cast<NCDMBaseSpecies*>` filter — flagged with `TODO(architecture)`, redesigned alongside `Omega0_b/Omega0_cdm` later.

## Architecture

### `BaseSpecies::PerturbLayout`

```cpp
class BaseSpecies {
 public:
  // Polymorphic layout base — one subclass per concrete species. Each layout
  // holds the indices/sizes/etc. its species needs to address its slots in
  // perturb_vector::y. The layout's lifetime equals its owning perturb_vector's.
  struct PerturbLayout {
    virtual ~PerturbLayout() = default;
  };

  // Each concrete species creates its own subclass. Called once per pv during
  // perturb_vector_init, before Register*PerturbationIndices.
  virtual std::unique_ptr<PerturbLayout> CreatePerturbLayout() const = 0;
  // ... see API surface below
};
```

Each concrete species defines a subclass:

```cpp
struct NCDMPerturbLayout : BaseSpecies::PerturbLayout {
  int l_max  = -1;
  int q_size = -1;
  std::vector<int> index_per_q;          // absolute offsets into pv->y
  int total_size() const { return q_size * (l_max + 1); }
};

struct ScalarFieldPerturbLayout : BaseSpecies::PerturbLayout {
  int idx_phi       = -1;
  int idx_phi_prime = -1;
};

struct PhotonsPerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1, idx_theta = -1, idx_shear = -1, idx_l3 = -1;
  int idx_pol0 = -1,  idx_pol1 = -1,  idx_pol2 = -1, idx_pol3 = -1;
  int l_max = -1;
  int l_max_pol = -1;
};
```

Concrete species methods downcast once at the top:

```cpp
void NCDMSpecies::PerturbDerivs(const PerturbLayout& base, /* ... */) {
  const auto& layout = static_cast<const NCDMPerturbLayout&>(base);
  // use layout.q_size, layout.l_max, layout.index_per_q[iq]
}
```

The downcast is `static_cast` (not `dynamic_cast`): the layout was returned by `this->CreatePerturbLayout()`, so the type is guaranteed.

### `perturb_vector`

```cpp
struct perturb_vector {
  // Per-species layouts, parallel to all_species_ (lex-key order).
  std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> species_layouts;

  // Non-species state: metric perturbations, perturbed recombination,
  // integration storage. These are NOT multi-instance and stay as bare fields.
  int index_pt_eta                                 = -1;  // synchronous metric
  int index_pt_phi                                 = -1;  // newtonian metric
  int index_pt_perturbed_recombination_delta_temp  = -1;
  int index_pt_perturbed_recombination_delta_chi   = -1;
  int index_pt_hv_prime                            = -1;  // vector mode metric
  int index_pt_V                                   = -1;  // vector mode (newtonian)
  int index_pt_gw                                  = -1;  // tensor metric
  int index_pt_gwdot                               = -1;

  int pt_size = 0;                                  // total slots = sum of layout sizes + bare fields
  std::vector<double> y_storage, dy_storage;
  double* y  = nullptr;
  double* dy = nullptr;
  std::vector<int> used_in_sources_storage;
  int* used_in_sources = nullptr;
};
```

All `index_pt_delta_g/theta_g/shear_g/l3_g/pol0_g/...` (photon), `index_pt_delta_b/theta_b` (baryon), `index_pt_delta_cdm/theta_cdm` (CDM), `index_pt_phi_scf/phi_prime_scf` (scalar field), `index_pt_delta_dcdm/delta_dr/...` (decaying CDM, dark radiation), `index_pt_delta_ur/theta_ur/shear_ur/l3_ur` (UR, plus `l_max_ur`), `index_pt_delta_idr/theta_idr/shear_idr/l3_idr` (plus `l_max_idr`), `index_pt_delta_idm_dr/theta_idm_dr` (IDM_DR), `index_pt_delta_idm_drmd/theta_idm_drmd/delta_idr_drmd/theta_idr_drmd` (IDM_DRMD/IDR_DRMD), `index_ncdm_/l_max_ncdm/q_size_ncdm/N_ncdm/index_pt_psi0_ncdm1/l_max_ncdm_storage/q_size_ncdm_storage` (NCDM family) are **deleted** from `perturb_vector`. They live in their species's layout subclass.

### `SpeciesCollection` indexing

`SpeciesCollection` gains `operator[](size_t i) const` returning `BaseSpecies*` (no key wrapper). The keyed-lookup API (`at(key)`, `count(key)`, `find(key)`, `photons()`, `baryons()`) stays. Existing key-iteration `for (auto& [name, sp] : all_species_)` continues to work; new code can use `for (size_t i = 0; i < all_species_.size(); ++i) { all_species_[i]->Foo(); }`.

### Module dispatch

```cpp
// e.g. in perturb_derivs_member, scalar mode:
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->PerturbDerivs(*ppw->pv->species_layouts[i],
                                 tau, y, dy, *pppaw);
}
```

`O(1)` direct vector indexing, no map lookup. The "must be parallel" invariant is held by `perturb_vector_init` (where `species_layouts.push_back(sp->CreatePerturbLayout())` is called once per species in `all_species_` iteration order).

### Composite species (nested composition)

`DNCDM_DR_Species`, `DCDM_DR_Species`, `IDM_DR_IDR_Species`, `IDM_DRMD_IDR_DRMD_Species` are `CompositeSpecies` whose children are NOT first-class members of `all_species_`. The composite's layout subclass owns its children's layouts as nested members:

```cpp
struct DNCDM_DR_PerturbLayout : BaseSpecies::PerturbLayout {
  NCDMPerturbLayout dncdm;
  DRPerturbLayout   dr;
};

void DNCDM_DR_Species::RegisterPerturbationIndices(
    PerturbLayout& base, const precision* ppr, int& index_pt,
    const perturb_workspace* ppw, int gauge) {
  auto& my = static_cast<DNCDM_DR_PerturbLayout&>(base);
  dncdm_->RegisterPerturbationIndices(my.dncdm, ppr, index_pt, ppw, gauge);
  dr_->RegisterPerturbationIndices(my.dr, ppr, index_pt, ppw, gauge);
}

void DNCDM_DR_Species::PerturbDerivs(const PerturbLayout& base, /* ... */) {
  const auto& my = static_cast<const DNCDM_DR_PerturbLayout&>(base);
  dncdm_->PerturbDerivs(my.dncdm, /* ... */);
  dr_->PerturbDerivs(my.dr, /* ... */);
}
```

Children's layouts live ONLY inside the composite's layout. Children never appear in `pv->species_layouts`. The composite is fully responsible for routing layout pieces to its children.

### Bookkeeping hooks

```cpp
class BaseSpecies {
  // No StashPerturbationLayout — old pv still holds its species_layouts during
  // a switch, so "previous layout" = old_pv->species_layouts[i].

  virtual void CopyPerturbationsAcrossSwitch(
      const PerturbLayout& old_layout,
      const PerturbLayout& new_layout,
      const double* old_y, double* new_y,
      const PerturbSwitchContext& ctx) const {}

  virtual void MarkUsedInSources(
      const PerturbLayout& layout,
      int* used_in_sources) const {}
};

struct PerturbSwitchContext {
  double k        = 0.;
  double a        = 0.;
  double a_today  = 1.;
  const double* pvecback = nullptr;
};
```

Module wiring in `perturb_vector_init`'s "switching approximation" branch:

```cpp
PerturbSwitchContext ctx{ k, a_now, pba->a_today, ppw->pvecback };
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->CopyPerturbationsAcrossSwitch(
      *ppw->pv->species_layouts[i],   // old
      *new_pv->species_layouts[i],    // new
      ppw->pv->y, new_pv->y,
      ctx);
}
// after y_storage allocated and used_in_sources defaulted to TRUE:
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->MarkUsedInSources(*new_pv->species_layouts[i],
                                     new_pv->used_in_sources);
}
```

The 8 NCDM-shaped switch-copy blocks (`ppv->y[ppv->index_ncdm_[n][q]+l] = ppw->pv->y[...]`), the analogous UR `ufa`/IDR `tca/rsa`/IDM_DRMD `tca` slot-by-slot copy blocks, and the NCDM `used_in_sources` masking block in `perturb_vector_init` all collapse into the two loops above. Per-species copy logic (slot-by-slot or NCDM fluid-approximation collapse) lives entirely in each species's `CopyPerturbationsAcrossSwitch` override.

The non-NCDM switch-copy blocks (UR ufa, IDR rsa/tca, IDM_DRMD tca) are functionally slot-by-slot copies — porting them to the species hook is mechanical. This spec performs the NCDM-family migrations (the only ones with non-trivial collapse logic, in NCDM's fluid approximation case); the UR/IDR/IDM_DRMD migrations are deferred to a follow-up because the slot-by-slot ones are pure mechanical refactor with no behaviour change.

## API surface on `BaseSpecies`

```cpp
class BaseSpecies {
 public:
  struct PerturbLayout { virtual ~PerturbLayout() = default; };
  virtual std::unique_ptr<PerturbLayout> CreatePerturbLayout() const = 0;

  virtual void RegisterPerturbationIndices(
      PerturbLayout& layout, const precision* ppr, int& index_pt,
      const perturb_workspace* ppw, int gauge) {}
  virtual void RegisterVectorPerturbationIndices(
      PerturbLayout& layout, int& index_pt,
      const perturb_workspace* ppw, int gauge) {}
  virtual void RegisterTensorPerturbationIndices(
      PerturbLayout& layout, int& index_pt,
      const perturb_workspace* ppw, int gauge) {}

  virtual void PerturbDerivs(
      const PerturbLayout& layout, double tau,
      const double* y, double* dy,
      const perturb_parameters_and_workspace& ppaw) = 0;
  virtual void PerturbVectorDerivs(
      const PerturbLayout&, double, const double*, double*,
      const perturb_parameters_and_workspace&) {}
  virtual void PerturbTensorDerivs(
      const PerturbLayout&, double, const double*, double*,
      const perturb_parameters_and_workspace&) {}

  virtual double Delta(
      const PerturbLayout& layout, const double* y,
      const double* pvecback, const perturb_workspace* ppw) const = 0;
  virtual double Theta(
      const PerturbLayout&, const double*, const double*,
      const perturb_workspace*) const = 0;
  virtual double DeltaP(
      const PerturbLayout&, const double*, const double*,
      const perturb_workspace*) const = 0;
  virtual double RhoPlusPShear(
      const PerturbLayout&, const double*, const double*,
      const perturb_workspace*) const = 0;

  virtual void ApplyInitialConditions(
      const PerturbLayout& layout, double* y,
      const PerturbIcContext& ctx) {}
  virtual void FillSources(
      const PerturbLayout& layout, const double* y, const double* dy,
      PerturbSourceContext& ctx) {}

  virtual void CopyPerturbationsAcrossSwitch(
      const PerturbLayout& old_layout, const PerturbLayout& new_layout,
      const double* old_y, double* new_y,
      const PerturbSwitchContext& ctx) const {}
  virtual void MarkUsedInSources(
      const PerturbLayout& layout, int* used_in_sources) const {}

  virtual void ContributeTensorGwSource(
      double a, double a_today, const double* y,
      perturb_workspace* ppw,
      const PerturbLayout& layout) const {}
  virtual void WriteTensorOutputColumnTitles(
      char* tensor_titles) const {}

  // Output / transfer (no PerturbLayout dependency — transfer indices are
  // species-instance state set ONCE, no per-thread race):
  virtual void RegisterTransferIndices(
      int& index_type, const struct perturbs* ppt) {}
  virtual void WriteOutputColumns(
      PerturbColumnWriter& w, const PerturbationsModule& mod,
      enum file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const {}
  virtual void PrintVariables(
      PerturbColumnWriter& w, double tau, const double* y,
      const PerturbationsModule& mod, const perturb_workspace* ppw) const {}

  // Shooter (cosmology-side, executed during InputModule shooting; not
  // per-thread — runs once per parameter-search iteration):
  virtual void RegisterShootingIndices(
      int& index_sh, FileContent* pfc, const background& pba) {}
  virtual void ComputeShootingResidual(
      const ShootingContext& ctx, double* residual) const {}
  virtual void ComputeShootingGuess(
      const ShootingContext& ctx, double* guess) const {}
};
```

## Migration strategy

The migration uses a **dual-write** pattern: every species's `Register*PerturbationIndices` populates BOTH its new layout fields AND the legacy `pv->index_pt_*` field (for migrating species) or `pv->index_ncdm_[n]` array (for NCDM-family). Module readers that haven't yet been migrated continue to read the legacy fields; migrated readers use the layout. Once every reader of a given field is migrated, the legacy field is deleted in a final cleanup commit.

### Per-species migration order

1. **Lambda** — no perturbation slots. Trivial layout (empty). Establishes the pattern with minimum risk.
2. **Fluid** — single `idx_delta`, `idx_theta`. Simple.
3. **CDM** — single `idx_delta` (synchronous gauge: theta is metric-derived, no `idx_theta`).
4. **DCDM, DR** — single `idx_delta`, `idx_theta` each (DR is photon-like radiation with multipoles; the layout has `l_max`, multipole indices).
5. **ScalarField** — `idx_phi`, `idx_phi_prime`.
6. **UR** — `idx_delta/theta/shear/l3`, `l_max`. The `ufa` switch participates in `CopyPerturbationsAcrossSwitch` (slot-by-slot copy when ufa is on; full hierarchy otherwise).
7. **IDR** — same shape as UR, plus `idr_nature` (free-streaming vs fluid) gating which fields are used.
8. **IDM_DR, IDM_DRMD** — `idx_delta/theta` plus tca interaction terms.
9. **Baryons** — `idx_delta/theta` (and `idx_pol0_b` etc. for tensor mode).
10. **Photons** — most fields: scalar (delta/theta/shear/l3/pol0/pol1/pol2/pol3, `l_max_g`, `l_max_pol_g`), tensor (subset), tight-coupling approximation logic.
11. **NCDM** — full hierarchy + fluid approximation collapse path (the non-trivial switch case). Replaces flat `pv->index_ncdm_[n]/l_max_ncdm[n]/q_size_ncdm[n]` arrays.
12. **DNCDM** — full hierarchy, no fluid approximation. Per-species `dlnf0_dlnq` accessor differs from NCDM (uses `pvecback[bg_dlnfdlnq_index_ + iq]`).
13. **NCDMInteractingSpecies** — inherits from NCDMSpecies, no override needed.
14. **DCDM_DR composite** — nested composition of `DCDMSpecies` + `DRSpecies` layouts.
15. **DNCDM_DR composite** — nested composition of `DNCDMSpecies` + `DRSpecies` layouts.
16. **IDM_DR_IDR composite, IDM_DRMD_IDR_DRMD composite** — nested composition.

For each species the migration is:
1. Add `XYZPerturbLayout` subclass next to the species declaration.
2. Add `CreatePerturbLayout()` override returning `std::make_unique<XYZPerturbLayout>()`.
3. Have `Register*PerturbationIndices` take `(PerturbLayout&, ...)` and populate the layout (continue dual-writing to the legacy `pv->*` field).
4. Migrate the species's reader methods (`PerturbDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`) to take and read from the layout.
5. Migrate module-side reads of the species's `pv->*` field to read from `pv->species_layouts[i]` via downcast.
6. Build, run regression scenarios with `output/scen_*.dat` cleared first, verify class exits 0, diff against baseline. Bit-identical for non-NCDM-family scenarios; NCDM-family covered in steps 11–13.

After all species migrate, delete the now-unused legacy fields from `perturb_vector` (single cleanup commit per category: photons, baryons, NCDM family, etc.).

### `ncdm_id` removal

After all NCDM-family species (NCDM, DNCDM, NCDMInteracting, DNCDM_DR composite) migrate to per-species layouts, the integer `ncdm_id` becomes unused. Final removal:

- Delete `ncdm_id_` field, `SetNcdmId(int)` virtual, `ncdm_id()` accessor on every NCDM-family species.
- Delete the `ctx.ncdm_id_next` counter from `SpeciesBuildContext`. Delete the `(*ctx.ncdm_id_next)++` increments from each NCDM-family `CreateAll`.
- Delete `pba->N_ncdm`. Replace consumers (5 sites) with sums over `all_species_` filtered by `dynamic_cast<NCDMBaseSpecies*>`. Mark each substitution `// TODO(architecture)` (deferred follow-up).
- Delete `pba->Omega0_ncdm_tot`. Replace consumers similarly.
- Delete `pba->N_decay_dr`. Replace consumers similarly.
- Delete `PerturbationsModule::ncdm_species_sorted_` member and its construction in `perturbations_init`.
- `PerturbationsModule::RescaledNCDMPerturbations(int n_ncdm, ...)` → `(BaseSpecies*, ...)`. Caller passes the species pointer it already has in scope.

### Output column titles

This change is NCDM-family only — other species (Photons, Baryons, CDM, ScalarField, …) currently have unindexed titles (`delta_b`, `theta_b`, `phi_scf`, `phi_prime_scf`, …) and stay unchanged. When the input parser eventually accepts multi-instance non-NCDM species, those titles will need a similar rename (out of scope here).

The `[%d]` indexing in NCDM-family titles is replaced by instance names. With `nu1.type = ncdm_standard`, the instance is named `nu1`; titles use that name directly:

- `(.)rho_ncdm[0]` → `(.)rho_nu1`
- `delta_ncdm[0]` → `delta_nu1`
- `theta_ncdm[0]` → `theta_nu1`
- `shear_ncdm[0]` → `shear_nu1` (tensor mode)
- `(.)number_ncdm[0]` → `(.)number_nu1`
- `(.)p_ncdm[0]` → `(.)p_nu1`
- `(.)pseudo_p_ncdm[0]` → `(.)pseudo_p_nu1`
- `(.)dlnfdlnq_ncdm[0][q]` → `(.)dlnfdlnq_nu1[q]`

For DNCDM, equivalent: `dncdm1.m = ...` produces `(.)rho_dncdm1`. The DNCDM_DR composite key (formerly `DNCDM_DR_<ncdm_id>`) becomes `<dncdm_instance_name>` and the DR child becomes `<dncdm_instance_name>_DR`.

This is a **breaking change for downstream parsers** that match column titles by the `[%d]` regex. Documented in the PR description.

### Shooter

Shooter / root-finding runs inside `InputModule` as part of input-parameter resolution. It is **not parallelised** — shooting happens before any per-k-mode evolution, and `InputModule` runs single-threaded. So the shooter API doesn't need per-thread state; species can store their shooter slot indices as plain instance members (set during `RegisterShootingIndices`).

Today `InputModule::input_class_root_finder` has hardcoded cases for `theta_s`, `Omega_dcdmdr` (DCDM_DR), `Omega_dncdm_decay_dr` (DNCDM_DR), `phi_ini` / `phi_prime_ini` (ScalarField). With this work, each species owns its shooter slots:

```cpp
struct ShootingContext {
  FileContent* pfc;
  const background* pba;
  // Currently-fitted parameter values during root finding:
  const double* current_values;
  // Indices of this species's slots in current_values:
  // populated by RegisterShootingIndices.
};

class BaseSpecies {
  virtual void RegisterShootingIndices(int& index_sh,
                                        FileContent* pfc,
                                        const background& pba) {}
  virtual void ComputeShootingResidual(const ShootingContext& ctx,
                                        double* residual) const {}
  virtual void ComputeShootingGuess(const ShootingContext& ctx,
                                     double* guess) const {}
};
```

Each species registers as many slots as it needs (e.g., DNCDM_DR registers 1 — `Omega_dncdm_decay_dr`). InputModule's shooter loops over `all_species_` calling `RegisterShootingIndices` first to allocate, then `ComputeShootingResidual` / `ComputeShootingGuess` during root-finding iterations. The cosmology-level `theta_s` slot remains module-allocated (it is not species-owned).

`DNCDM_DR_Species`, `DCDM_DR_Species`, and `ScalarFieldSpecies` are the species providing shooter overrides for this work.

## Testing strategy

Five regression scenarios at `test/scenarios/`:
- `ncdm_single.ini` — one NCDM, full Cl.
- `ncdm_multi_unsorted.ini` — two NCDMs (`nu_a`, `nu_b`) named so construction order ≠ lex order. Detects accidental order dependence.
- `ncdm_self_interacting.ini` — one self-interacting NCDM, full Cl.
- `dncdm_dr.ini` — one DNCDM_DR, background-only (full Cl with DNCDM is too slow for the regression cycle).
- `ncdm_dncdm_idmdr_combined.ini` — NCDM + DNCDM_DR + IDM_DR_IDR, background-only.

A baseline of these scenarios was captured at commit `2f64cffe` and stored in `/tmp/ncdm_baseline/`.

### Regression cycle (mandatory order)

After every behaviour-touching commit:

1. `rm -f output/scen_*.dat` — clear stale outputs from previous successful runs.
2. `make class -j` — clean build.
3. `for f in test/scenarios/*.ini; do ./class "$f" || { echo "FAIL: $f"; exit 1; }; done` — explicitly check exit code per scenario.
4. Diff each generated `.dat` against `/tmp/ncdm_baseline/`. Bit-identical (`max_rel_diff = 0.00e+00`) expected for all dual-write commits and all migrations except column-title renames.

Single-thread runs (with `OMP_NUM_THREADS=1`) are NOT a substitute. Every regression cycle uses the default thread count to exercise the per-thread `pv` invariants. Tasks that knowingly diverge bit-identicality (only the column-title rename and the eventual int→species-pointer dispatch in `RescaledNCDMPerturbations`) document the expected diff in the task description.

After all migrations land, run:
```bash
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py
```

This is the existing test_scenario suite under `python/`. It exercises a broader range of input combinations than the 5 regression scenarios.

### Final audit checks

- `grep -rn "ncdm_id" source/ species/ include/` returns 0 matches.
- `grep -rn "index_pt_phi\b\|index_pt_delta_g\b" source/perturbations_module.h` returns 0 matches (deleted).
- `grep -rn "pba->N_ncdm\|pba->Omega0_ncdm_tot\|pba->N_decay_dr" source/ species/` returns 0 matches.
- `grep -rn "ncdm_species_sorted_\b\|index_pt_psi0_ncdm1\b" source/` returns 0 matches.

## Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| **Multi-threading data race resurfaces.** The whole point of the refactor is to fix this. Verified by the regression cycle running multi-threaded. | Every regression run uses default `OMP_NUM_THREADS`. Single-thread is a debugging tool only. |
| **Iteration-order shift breaks output column ordering.** With `species_layouts` parallel to lex-sorted `all_species_`, multi-NCDM scenarios get column reordering vs construction-order baseline. | `ncdm_multi_unsorted` scenario is the canary — the column-title rename in §5 (output titles) is the change that knowingly diverges bit-identicality there. After that point, `nu_a` precedes `nu_b` in column order regardless of construction order. |
| **Breaking change to downstream parsers** that match `[%d]` in column titles. | Documented in the PR. Migration guide lines: `delta_ncdm[0]` → `delta_<instance>`. |
| **`background.h` ABI change.** Removing `pba->N_ncdm/Omega0_ncdm_tot/N_decay_dr` breaks the C API exposed via `class.h`. | Python wrapper rebuild covered by `python -m pytest` step. C-API consumers external to the repo flagged in PR. |
| **Polymorphic `unique_ptr<PerturbLayout>` overhead** — virtual destructor + downcast per access. | Downcast is `static_cast` (zero overhead). Virtual destructor fires only on `pv` destruction (once per switch). Hot-path is `layout.field_x` after the cast — same cost as today's direct member access. |
| **Composite layout-routing mistakes.** A composite that fails to forward all child layout pieces to its children silently leaves the child's slots un-initialised. | Each composite has a regression scenario that exercises its evolution path (`dncdm_dr.ini`, `ncdm_dncdm_idmdr_combined.ini`). NaN propagation breaks the integrator quickly. |
| **`StashPerturbationLayout` removal mishandled.** With both old and new `pv->species_layouts` available during a switch, `CopyPerturbationsAcrossSwitch` reads from the old one. If the old `pv` is destroyed too early, the read is on freed memory. | `pv` is moved into `ppw->pv` only AFTER the copy loop runs. The old pv is destroyed at the end of `perturb_vector_init`. Standard `delete` ordering. |
| **Stale `output/scen_*.dat` masking real failures.** This was the bug pattern in the previous attempt. | `rm -f output/scen_*.dat` is the first step of every regression cycle. Class exit code is explicitly checked. |
