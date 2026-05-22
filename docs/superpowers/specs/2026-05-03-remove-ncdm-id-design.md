# Remove `ncdm_id` — Design

> **Status:** spec for a single PR; implementation plan to follow under `docs/superpowers/plans/`.
> **Repo:** `AarhusCosmology/CLASSpp`
> **Date:** 2026-05-03

## Goal

Eliminate the integer `ncdm_id` and every flat-array sibling it indexes into. After this PR, no file in `source/`, `species/`, or `include/` contains the symbol `ncdm_id`. NCDM-family quantities are owned by their species instances; module code iterates `all_species_` and dispatches via virtual hooks instead of reaching into shared per-id arrays.

## Background

CLASSpp's species refactor has reached the point where `all_species_` is a sorted-vector container keyed by the explicit instance-name string, composite species (DCDM_DR, IDM_DR_IDR, IDM_DRMD_IDR_DRMD, DNCDM_DR) are first-class sectors, and most module-side species dispatch goes through the `BaseSpecies` virtual interface. The `ncdm_id` integer is the last load-bearing remnant of the legacy "index everything by an int handle" pattern. It still wires together:

1. **Per-instance perturbation storage** — `pv->index_ncdm_[id]` (a `std::map<int, std::vector<int>>`), `pv->l_max_ncdm[id]`, `pv->q_size_ncdm[id]`, plus `pv->index_pt_psi0_ncdm1`, `pv->N_ncdm`.
2. **Output column titles** — `(.)rho_ncdm[%d]`, `delta_ncdm[%d]`, `lnf_dncdm[%d][%d]`, etc.
3. **Deterministic ordering** — `[](a,b){ return a->ncdm_id() < b->ncdm_id(); }` comparators feeding `ncdm_species_sorted_`, the DR registration loop, and the tensor-mode title emission.
4. **DNCDM_DR composite naming** — `"DNCDM_DR_" + std::to_string(ncdm_id)` and `"DNCDM_DecayRadiation_" + std::to_string(ncdm_id)`.
5. **Background "ncdm1" first-slot indices** — `index_bg_number_ncdm1_`, `index_bg_pseudo_p_ncdm1_` (set from the `ncdm_id == 0` species).
6. **Input-module shooter** — `target_values[idx + dncdm_id]` (Omega/Neff/deg targets ordered positionally by ncdm_id).
7. **`RescaledNCDMPerturbations(int n_ncdm, ...)`** — public API that takes an integer ncdm_id and looks up the matching DNCDM.
8. **Background struct** — `pba->N_ncdm`, `pba->Omega0_ncdm_tot`, `pba->N_decay_dr` (count + running sum siblings of the per-id pattern; the per-id parameter arrays themselves were removed in earlier work).

This PR is scope-A: full removal in one PR. The result is a uniform "loop over `all_species_`, dispatch via virtual hook, each species owns its own state" pattern with no NCDM-shaped flat arrays anywhere in the codebase.

## Non-goals

- Removing `pba->has_ncdm` / `pba->has_ncdm_decay_dr` — these are has-guards on the opportunistic-cleanup track, not part of this PR.
- Migrating the *non-NCDM* approximation-switch y-copy blocks (UR ufa, IDR rsa/tca, IDM_DRMD tca) into the new `CopyPerturbationsAcrossSwitch` pattern. Those blocks stay as-is in this PR; the hook is introduced in a generalisable way so the follow-up PR is mechanical.
- Redesigning `theta_s` shooting as a "global cosmology" species target — it stays on `InputModule`.

## Architectural endpoint

After this PR:

- The symbol `ncdm_id` does not appear in `source/`, `species/`, or `include/`.
- The `perturb_vector` struct fields `index_ncdm_`, `l_max_ncdm`, `q_size_ncdm`, `l_max_ncdm_storage`, `q_size_ncdm_storage`, `index_pt_psi0_ncdm1`, `N_ncdm` are deleted.
- The `background` struct fields `N_ncdm`, `Omega0_ncdm_tot`, `N_decay_dr` are deleted.
- `PerturbationsModule::ncdm_species_sorted_` is deleted.
- The module-side `dynamic_cast<NCDMSpecies*>` / `dynamic_cast<DNCDMSpecies*>` filtering loops are deleted; module code iterates `all_species_` and dispatches via virtual hooks. (One temporary exception: `Omega0_ncdm_tot` consumers — flagged below.)
- Each NCDM-family species owns its own perturbation indices, output columns, transfer-function indices, and shooter indices.
- The `DNCDM_DR_Species` composite uses its DNCDM child's instance name as both its `BaseSpecies::name()` and its `SpeciesCollection` insertion key. The DR child's name is `<instance>_DR`.

## Approach (hybrid — Approach 3 from brainstorm)

Per-species storage everywhere. For physics-shaped loops that already have a natural virtual on `BaseSpecies` (source/stress accumulation, IC application, derivative dispatch, output emission, tensor sources), push into the species. For the bookkeeping inside `perturb_vector_init` (the switch-time y-copying and `used_in_sources` masking), introduce generic `CopyPerturbationsAcrossSwitch` and `MarkUsedInSources` hooks on `BaseSpecies` that default to a no-op.

This satisfies the "loop + dispatch, never type-pick" architectural rule without bloating `BaseSpecies` with NCDM-only hooks. The bookkeeping hooks generalise to UR/IDR/IDM_DRMD switches in a follow-up PR.

## Component design

### 1. Per-species perturbation storage

Each NCDM-family species (`NCDMSpecies`, `DNCDMSpecies`, `NCDMInteractingSpecies`, `DNCDM_DecayRadiationSpecies`) gains its own perturbation storage on `NCDMBaseSpecies`:

```cpp
int pt_l_max_   = -1;             // truncation multipole; respects fluid-approx switch
int pt_q_size_  = -1;             // number of momentum bins; respects fluid-approx switch
std::vector<int> pt_index_per_q_; // base offset into y/dy for each momentum bin

// Snapshot of the previous registration, for CopyPerturbationsAcrossSwitch:
int pt_l_max_previous_   = -1;
int pt_q_size_previous_  = -1;
std::vector<int> pt_index_per_q_previous_;
```

Public accessors: `pt_l_max()`, `pt_q_size()`, `pt_index_at(int q)`, `pt_total_size()` (= `pt_q_size_ * (pt_l_max_ + 1)`).

These are written exclusively by the species's own `RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr, int& index_pt, const perturb_workspace* ppw, int gauge)`, which:
- decides fluid-approx-on vs -off `(l_max, q_size)`,
- advances `index_pt` by `pt_total_size()`,
- stores the per-q offsets into `pt_index_per_q_`,
- no longer touches `pv->index_ncdm_`, `pv->l_max_ncdm`, `pv->q_size_ncdm`, `pv->index_pt_psi0_ncdm1`, or `pv->N_ncdm`.

The `DNCDM_DR_Species` composite reaches into its DNCDM child via the existing `dncdm()` accessor (`dncdm_.pt_index_at(q) + l`) instead of the shared map.

`NCDMInteractingSpecies` does not override `RegisterPerturbationIndices` — it inherits the standard NCDM slot layout from `NCDMSpecies` and only adds collision terms inside `PerturbDerivs`.

**Known-unknown:** `pv->index_pt_psi0_ncdm1` (the legacy "first NCDM slot" sentinel) — set when the first NCDM species registers, then read in some print/output path not yet pinned down. The implementation plan includes a step to grep all read sites and either rewrite each through the per-species accessor or delete it if dead.

### 2. Loop migration

#### (a) Physics loops — extend existing `BaseSpecies` virtuals

| Site (today, `perturbations_module.cpp`) | After |
|---|---|
| Source accumulation NCDM block (~line 6390) | `BaseSpecies::AccumulateSources` is the dispatch hook; NCDM-family override reads its own `pt_index_at(q)`. Module-side loop is `for (auto& [_, sp] : all_species_) sp->AccumulateSources(...)` with no NCDM-special block. |
| Tensor source title emission (~line 2540) | NCDM-family override of `WriteTensorOutputColumns(PerturbColumnWriter&)` emits `delta_<instance>`, `theta_<instance>`, `shear_<instance>`. |
| Tensor `PerturbDerivs` for NCDM (~line 7455) | New `BaseSpecies::PerturbDerivsTensor(...)` hook (default no-op); NCDM-family override moves the per-momentum-bin loop into the species and reads its own slots. The existing scalar `PerturbDerivs` is unchanged. |
| `class_test` initial-time-too-late warning (~line 2555), embeds `ncdm_sp->ncdm_id()` | Embed `ncdm_sp->name()` (instance name) instead. |

#### (b) Bookkeeping loops — new generic `BaseSpecies` hooks

Three hooks on `BaseSpecies`, all defaulting to no-op:

```cpp
virtual void StashPerturbationLayout();        // copies pt_*_ → pt_*_previous_
virtual void CopyPerturbationsAcrossSwitch(    // reads old_y at pt_*_previous_, writes new_y at pt_*_
    const double* old_y, double* new_y) const;
virtual void MarkUsedInSources(int* used) const; // marks species' slots in used_in_sources
```

Module flow at every approximation switch (replaces the NCDM-shaped block of lines ~3870–4620):

```cpp
for (auto& [_, sp] : all_species_) sp->StashPerturbationLayout();
auto* new_pv = new perturb_vector();
int index_pt = 0;
for (auto& [_, sp] : all_species_) sp->RegisterPerturbationIndices(new_pv, ppr, index_pt, ppw, gauge);
new_pv->y_storage.resize(index_pt); new_pv->y = new_pv->y_storage.data();
for (auto& [_, sp] : all_species_) sp->CopyPerturbationsAcrossSwitch(old_pv->y, new_pv->y);
for (auto& [_, sp] : all_species_) sp->MarkUsedInSources(new_pv->used_in_sources);
delete old_pv; ppw->pv = new_pv;
```

The protocol is wrapped in a non-virtual `BaseSpecies::Migrate*` method (or a free function in `perturb_vector_init`) that calls the three virtual hooks in the fixed order, so call sites cannot reorder them.

NCDM-family overrides handle the ncdmfa-on (truncate to `l ≤ 2`, `q_size = 1`) and ncdmfa-off (rebuild from truncated state via the existing rescaling logic) transitions internally. The other transition types (UR ufa, IDR rsa/tca, IDM_DRMD tca) keep their current module-side blocks in this PR and migrate in a follow-up.

#### (c) Identity / API changes

- `RescaledNCDMPerturbations(int n_ncdm, ...)` becomes `RescaledNCDMPerturbations(BaseSpecies* sp, ...)`. The single call site already iterates species.
- `ncdm_species_sorted_` is deleted entirely. Iteration is `for (auto& [_, sp] : all_species_)` with virtual dispatch; `all_species_` is already deterministic (instance-name lexical sorted-vector).

### 3. Shooter index registration

Modeled on `RegisterIntegrationIndices(int& index_bi)`. New hooks on `BaseSpecies` (all default to no-op):

```cpp
virtual void RegisterShootingIndices(
    int& index_sh,
    FileContent* pfc,
    const background& ba);

virtual void ComputeShootingResidual(
    double* output,
    const double* xvalues,
    const ShootingContext& ctx) const;

virtual void ComputeShootingGuess(
    double* xguess,
    double* dxdy,
    const ShootingContext& ctx) const;
```

`ShootingContext` bundles what today's `input_try_unknown_parameters` / `input_get_guess` cases need (cosmology pointer, computed background table for residual evaluation, verbose flag).

Per-species storage mirrors the `index_bi_*` pattern: each species holds its own `int sh_*_index_ = -1` plus the captured target value, populated during `RegisterShootingIndices`, read during residual + guess.

Flow inside `InputModule`:

```cpp
int index_sh = 0;
input_module_->RegisterCosmologicalShootingIndices(index_sh);            // theta_s only (see below)
for (auto& [_, sp] : all_species_) sp->RegisterShootingIndices(index_sh, &file_content_, background_);
// unknown_parameters_size = index_sh; allocate xvalues/output as today.

// During solver iteration:
input_module_->ComputeCosmologicalShootingResidual(output, xvalues, ctx);
for (auto& [_, sp] : all_species_) sp->ComputeShootingResidual(output, xvalues, ctx);
// (analogous for ComputeShootingGuess)
```

After this PR, `theta_s` is the **only** target left on `RegisterCosmologicalShootingIndices` — every other current target is species-keyed and moves into a species below.

What moves out of `InputModule`:

- `omega_dncdmdr` / `Omega_dncdmdr` / `Neff_ini_dncdm` / `deg_ncdm_decay_dr` / `omega_ini_dncdm` / `Omega_ini_dncdm` cases (positional `dncdm_id`-indexed loops at lines ~3625–3665 and ~3771–3845 of `input_module.cpp`) → into `DNCDM_DR_Species::ComputeShootingResidual` / `ComputeShootingGuess`. Each instance computes only its own slot.
- `omega_dcdmdr` / `Omega_dcdmdr` / `omega_ini_dcdm` / `Omega_ini_dcdm` cases (single slot) → into `DCDM_DR_Species`.
- `Omega_scf` → into `ScalarFieldSpecies`.

The `pfzw->target_name`, `pfzw->target_size`, and global enum-dispatch switches are deleted. Species translate user-facing input into their own internal target choice during `RegisterShootingIndices`.

### 4. Background-struct cleanup

Delete from `background.h`:

- `int N_ncdm`
- `double Omega0_ncdm_tot`
- `int N_decay_dr`

Read-site rewrites (~17 total):

| Pattern today | Rewrite |
|---|---|
| `if (pba->N_ncdm > 0) { ... }` as a presence guard around code that iterates NCDM-family species | Delete the guard; iteration is naturally empty if no NCDM-family species are present. |
| `pba->N_ncdm` used as an array allocation size (`perturbations_module.cpp:1089, 1105, 6903–6905`) | Replaced by per-species transfer-function registration (see "Transfer functions" below) and locally-scoped scratch vectors sized by the per-species loop. |
| `pba->Omega0_ncdm_tot` (5 sites: `nonlinear_module`, `thermodynamics_module`, `input_module`, perturbations_module RSA path) | Each consumer does `double sum = 0; for (auto& [_, sp] : all_species_) if (auto* nsp = dynamic_cast<NCDMBaseSpecies*>(sp.get())) sum += nsp->GetOmega0()`. **This is the only spot in the PR that temporarily breaks the "no species-type picking" rule.** Flagged with a `TODO(architecture): redesign once Omega0_b/Omega0_cdm get the same treatment` comment. |
| `pba->N_decay_dr` | Same treatment — count of `DNCDM_DR_Species` instances, computed on demand. |

`InputModule` writes (`pba->N_ncdm = n_ncdm;`, `pba->Omega0_ncdm_tot = ...`, `input_module.cpp:346–364`) are deleted.

### 5. Output column titles + DNCDM_DR composite naming

| Today | After |
|---|---|
| `BaseSpecies::name()` of DNCDM_DR composite = `"DNCDM_DR_" + std::to_string(ncdm_id)` | `name()` = the DNCDM child's instance name (e.g. `"nu_decay1"`) |
| `BaseSpecies::name()` of DR child = `"DNCDM_DecayRadiation_" + std::to_string(ncdm_id)` | `name()` = `"<instance>_DR"` (e.g. `"nu_decay1_DR"`) |
| `SpeciesCollection` insertion key for the composite | The instance name (e.g. `"nu_decay1"`) |
| Sub-species accessed via `composite.dncdm()` / `composite.dr()` | unchanged |

Background-output column titles (each species's own `WriteBackgroundColumnTitles`):

| Today | After |
|---|---|
| `(.)rho_ncdm[0]`, `(.)number_ncdm[2]`, `(.)p_ncdm[1]` | `(.)rho_<instance>`, `(.)number_<instance>`, `(.)p_<instance>` |
| `lnf_dncdm[0][3]`, `dlnfdlnq_dncdm[0][3]` | `lnf_<instance>[3]`, `dlnfdlnq_<instance>[3]` (inner `[3]` is momentum-bin index, kept) |
| `(.)rho_dr_species` (DNCDM DR child) | `(.)rho_<instance>_DR` |

Perturbation-output columns (scalar + tensor, each species's `WriteOutputColumns` / `PrintVariables`):

| Today | After |
|---|---|
| `delta_ncdm[0]`, `theta_ncdm[0]`, `shear_ncdm[0]` | `delta_<instance>`, `theta_<instance>`, `shear_<instance>` |

Transfer functions (`perturbations_module.cpp:1089/1105` flat block-of-N):

- The flat `class_define_index(index_tp_delta_ncdm1_, ..., index_type, pba->N_ncdm)` block is deleted.
- A new `BaseSpecies::RegisterTransferIndices(int& index_type)` hook (default no-op) lets each species register zero or more transfer-function slots and store its own `int tp_delta_index_ = -1`, `int tp_theta_index_ = -1`, etc.
- NCDM-family overrides this hook to register one delta + one theta (+ tensor versions) per instance.
- Transfer-output emission iterates `all_species_` and dispatches into the species — same pattern as background-output columns.

Verbose budget prints in `background_module.cpp` (line ~559 — `"-> N_eff = %g for ncdm species %d"`, line ~2115ff) — the format strings change `%d` for ncdm_id to `%s` for instance name. Already inside per-species printer methods (`PrintNeffInfo`, `PrintMassInfo`, `PrintOmegaInfo`).

## Testing strategy

Build + smoke gates after every implementation task:

```bash
make class -j 2>&1 | tail -20
./class explanatory.ini 2>&1 | tail -5
```

Full regression gate before PR:

```bash
cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py
```

Targeted scenarios (.ini files) covering paths the smoke test does not:

- 1 plain NCDM species (`m_ncdm = 0.06`)
- ≥2 NCDM species with non-monotonic instance names (e.g. `nu_b.type = ncdm_standard`, `nu_a.type = ncdm_standard`) — verifies iteration is by instance-name lexical order, not legacy ncdm_id construction order
- 1 DNCDM_DR
- 1 NCDM self-interacting (G_eff > 0)
- A scenario combining NCDM + DNCDM_DR + IDM_DR_IDR (cross-composite dispatch)

For each scenario: build a numerical-regression baseline on `master` *before* the branch starts, then diff post-branch outputs against it. Acceptable diff window: bit-identical for `(.)rho_*` columns; physical quantities (Cl's, P(k)) within `1e-10` relative tolerance.

## Risks

1. **Iteration-order change (highest risk).** Today's NCDM iteration order = construction order = `ncdm_id`. After: `all_species_` order = instance-name lexical. For .ini files that name species in a different order than they appear (`nu_b` before `nu_a`), output column order changes; if any physics depends on iteration order at ulp scale, values drift. Mitigation: the multi-NCDM regression scenarios above; if drift appears, identify the order-dependence and either fix it or sort by a canonical criterion at the call site.

2. **Output column titles renamed.** Anything downstream that greps for `(.)rho_ncdm[0]`-style names breaks — Python wrapper, user analysis scripts, external test harnesses. Mitigation: scan `python/` for hard-coded NCDM column-name patterns before committing; update wrapper accessors. CHANGELOG entry calling out the rename.

3. **`background.h` ABI change.** Removing `N_ncdm` / `Omega0_ncdm_tot` / `N_decay_dr` from `struct background` is a binary-compat break for anything compiled against the old header — including the Python wrapper if built against an older header. Mitigation: ensure `setup.py` rebuilds against current headers; verify Cython side does not `cdef` these fields.

4. **Approximation-switch y-copy correctness.** The new `StashPerturbationLayout` → `RegisterPerturbationIndices(new)` → `CopyPerturbationsAcrossSwitch` protocol must run in exactly that order. If a species's `pt_*_previous_` is read after `pt_*_` has been overwritten, values are corrupted silently. Mitigation: wrap the protocol in a non-virtual method that calls the three virtual hooks in fixed order.

5. **`DNCDM_DR` composite-name rename ABI.** Anything that looks up `all_species_.at("DNCDM_DR_0")` breaks. Mitigation: after the rename, grep the codebase for any literal `"DNCDM_DR_"` string; the only legitimate construction site should be `dncdm_dr_species.cpp` itself.

6. **`pv->index_pt_psi0_ncdm1`** (the legacy "first NCDM slot" sentinel). Known-unknown from §1. The implementation plan includes a grep step to find every read site; expectation is they are either dead code or trivially replaceable by "the first NCDM-family species in `all_species_`" found via the iteration loop.

## PR-size note

This is a ~1500–2500 line diff across ~15 files — bigger than recent `#266` and `#264` PRs. If splitting becomes necessary the natural seam is:

- (a) per-species perturbation storage + column titles + composite-naming + `ncdm_species_sorted_` deletion + `RescaledNCDMPerturbations` API change,
- (b) shooter index registration + `background` cleanup.

Both halves are independently mergeable; (b) does not depend on (a) for compilation, only for the architectural endpoint.

## Build-system reminder

New `.cpp`/`.h` files (if any added — none currently planned, but `ShootingContext` may warrant a small header) must be added to all three build systems:

- `Makefile` (SPECIES_OPP / SOURCE_OPP lists)
- `setup.py` (species/source `source_files` lists)
- `CLASS.xcodeproj/project.pbxproj`
