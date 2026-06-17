# Species-registered transfer-source indices (#309)

- **Issue:** AarhusCosmology/CLASSpp#309 — "v4: species-registered source types (remove per-species `has_source_*`/`index_tp_*` from `PerturbationsModule`)"
- **Date:** 2026-06-16
- **Status:** Design approved; ready for implementation plan.
- **Predecessor:** PR #317 dissolved the per-species `dynamic_cast` *contribution* loops and deferred the source-slot layer here.

## Problem

`source/perturbations_module.h:82-172` hard-codes a `index_tp_delta_X_` / `index_tp_theta_X_` member **and** a `has_source_delta_X_` / `has_source_theta_X_` flag for every species that can emit a transfer function (`dcdm`, `fld`, `scf`, `dr`, `ur`, `idr`, `idm_dr`, `idr_drmd`, `idm_drmd`, `ncdm1`, …). Setting those up takes three module-owned mechanisms that all encode species identity in the module:

1. **Stringly-typed gating** — `has_source_delta_X_` is set via `all_species_.count("DCDM_DR")`, `all_species_.count("Fluid")`, … (`perturbations_module.cpp:957-1013`).
2. **A fixed `class_define_index` block** — one line per species type (`perturbations_module.cpp:1068-1117`).
3. **Multi-instance `dynamic_cast` slot loops** — `SetSourceSlot` assignment for NCDM and for DR-emitting composites (`perturbations_module.cpp:639-664`), the layer #308 explicitly deferred to this issue.

Consequently "add a species without touching the main code" holds only for species that want no transfer-function output. This is the last large structural gap in the species-as-plugins design.

## Key scoping finding

The per-species `index_tp_delta_X_` / `index_tp_theta_X_` source types are read in **exactly two kinds of place**:

- the owning species' `FillSources` (writes the source value), and
- the owning species' `WriteOutputColumns` (emits the `d_X` / `t_X` Tk/vTk column).

The transfer, output, spectra, nonlinear, and lensing modules consume **only module-owned** source types — `delta_m`, `delta_cb`, `theta_m`, `theta_cb`, `phi`, `phi_prime`, `phi_plus_psi`, `psi`, the CMB temperature terms `t0`/`t1`/`t2`, and polarization `p` (verified: `transfer_module.cpp:436-450, 1032-1105`; no per-species `index_tp_*` reads in output/spectra/nonlinear/lensing). Therefore #309 is **contained to `source/perturbations_module.*` plus `species/`**; no transfer/output/spectra changes are required, and the module-owned types stay exactly as they are.

## Architectural principle: this is the *background* index pattern, not the perturbation one

Three index spaces exist on a species:

| Space | Running index | k-mode dependent? | Where state lives | Hook |
|-------|---------------|-------------------|-------------------|------|
| Background | `index_bg` | No (set once) | **Plain species member** (`index_bg_rho_`) | `RegisterBackgroundIndices(int&)` |
| Perturbation / y-vector | `index_pt` | **Yes** (approximation switches & `l_max(k)` change the active variable set) | `perturb_vector::species_layouts` (per-k, off the const/shared species) | `RegisterPerturbationIndices(PerturbLayout&, …)` |
| **Transfer source (this issue)** | `index_tp` | No (set once in `perturb_indices_of_perturbs`; `tp_size_[index_md]` identical for every k) | **Plain species member** | `RegisterTransferSourceIndices(int&, …)` (new) |

Transfer-source indices are computed once before the k-loop and never vary per k, exactly like background indices. They therefore belong as **plain members on the species**, written once at registration and read read-only thereafter — including inside the parallel k-loop, which is safe for the same reason `index_bg_rho_` is (write-once, then shared-immutable). They must **not** be pushed onto `PerturbLayout`; that machinery exists only because perturbation indices genuinely differ per k.

## Design

### 1. New hook on `BaseSpecies`

A default-no-op virtual placed beside the other register hooks in `species/base_species.h`:

```cpp
/**
 * Register this species' transfer-source (index_tp) slots, set once in
 * perturb_indices_of_perturbs() before the k-loop. Mirrors
 * RegisterBackgroundIndices: the species bumps the running index_tp by as
 * many slots as it needs and caches each slot in its own member.
 * Default: no-op (species that emit no transfer functions).
 */
virtual void RegisterTransferSourceIndices(int& index_tp,
                                           const SourceRequestContext& ctx) {}
```

Unlike `RegisterPerturbationIndices`, it takes **no `PerturbLayout&`** — that absence is the marker that this is k-independent, member-stored state. Unlike `RegisterBackgroundIndices`, it is **not** pure-virtual: many species emit no sources (Lambda, CDM in synchronous gauge for θ), so a no-op default avoids empty stubs.

### 2. `SourceRequestContext`

A small by-value struct (new, e.g. in `species/perturb_source_context.h`) carrying only the *output requests* a species cannot infer from itself:

```cpp
struct SourceRequestContext {
  bool wants_density;   // ppt->has_density_transfers
  bool wants_velocity;  // ppt->has_velocity_transfers
  int  gauge;           // possible_gauges; CDM registers theta only when != synchronous
};
```

(No `has_cl_number_count` field: every source that flag gates — `delta_m`, `theta_m`, `theta_cb` — is module-owned, not species-registered, so no species' registration depends on it.)

It is plain const input, not per-k state, so it pulls nothing onto the layout. (Gating is retained deliberately: computing δ/θ for non-trivial species — NCDM momentum-grid integrals especially — is a real cost worth skipping when no transfer output is requested.)

### 3. Call site — one flat loop

In `PerturbationsModule::perturb_indices_of_perturbs()`, scalar-mode branch, the per-species `has_source_*` gating block, the per-species `class_define_index` lines, and **both** `SetSourceSlot` `dynamic_cast` loops are replaced by:

```cpp
// --- module-owned head (stay as members; transfer reads these by name) ---
class_define_index(index_tp_t0_,        has_source_t_,        index_type, 1);
class_define_index(index_tp_t1_,        has_source_t_,        index_type, 1);
class_define_index(index_tp_delta_m_,   has_source_delta_m_,  index_type, 1);
class_define_index(index_tp_delta_cb_,  has_source_delta_cb_, index_type, 1);
class_define_index(index_tp_delta_tot_, has_source_delta_tot_,index_type, 1);
class_define_index(index_tp_theta_m_,   has_source_theta_m_,  index_type, 1);
class_define_index(index_tp_theta_cb_,  has_source_theta_cb_, index_type, 1);
class_define_index(index_tp_theta_tot_, has_source_theta_tot_,index_type, 1);

// --- per-species transfer sources: one clean loop, lex order = output order ---
const SourceRequestContext src_ctx{ppt->has_density_transfers,
                                   ppt->has_velocity_transfers,
                                   ppt->gauge};
for (auto& [name, sp] : all_species_)
    sp->RegisterTransferSourceIndices(index_type, src_ctx);

// --- module-owned tail (metric potentials) ---
class_define_index(index_tp_phi_, has_source_phi_, index_type, 1);
... (phi_prime, phi_plus_psi, psi, h, h_prime, eta, eta_prime, H_T_Nb_prime, k2gamma_Nb) ...

tp_size_[index_md] = index_type;
```

This is the flat loop `background_indices()` never achieved — possible here because per-species source types do not interleave with positionally-dependent module indices, and full-dissolve removes the shared base-block addressing (`index_tp_*_ncdm1_ + slot`) that forced background's `dynamic_cast`/`first_ncdm` caching.

The `SetSourceSlot` ordering loops at `perturbations_module.cpp:639-664` are deleted entirely; registration order (lex `all_species_` order) replaces them.

### 4. Member storage

Each source-writing **leaf** species gains members initialized to `-1` (matching `index_bg_rho_`):

```cpp
int index_tp_delta_ = -1;
int index_tp_theta_ = -1;
```

`class_define_index` writes the member only when its condition is true, leaving `-1` otherwise — so an absent/unrequested source keeps `-1`, exactly the existing convention. `FillSources` / `WriteOutputColumns` switch from `mod.index_tp_delta_X_` to the species' own `index_tp_delta_`. The `active` argument they already pass to `writer.Add(...)` is now simply `index_tp_delta_ >= 0`.

### 5. Per-species migration

Twelve source-writing units, in three structural shapes:

**(a) Single-instance leaves** — register own δ/θ into own members; swap reads in `FillSources`/`WriteOutputColumns`:

| Species | Sources | Gating |
|---|---|---|
| `PhotonsSpecies` | δ_g, θ_g | density / velocity |
| `BaryonsSpecies` | δ_b, θ_b | density / velocity |
| `CDMSpecies` | δ_cdm, θ_cdm | density / (velocity **and** gauge≠synchronous) |
| `FluidSpecies` | δ_fld, θ_fld | density / velocity |
| `ScalarFieldSpecies` | δ_scf, θ_scf | density / velocity |
| `UltraRelativisticSpecies` | δ_ur, θ_ur | density / velocity |

**(b) Interacting composites** — register two channels each:

| Species | Sources |
|---|---|
| `IDM_DR_IDR` | δ/θ for `idm_dr` and δ/θ for `idr` |
| `IDM_DRMD_IDR_DRMD` | δ/θ for `idm_drmd` and δ/θ for `idr_drmd` |

**(c) Multi-instance / DR-emitting** — each instance/child registers its **own** slot; composites forward to children like `CompositeSpecies::RegisterBackgroundIndices` already does, and read the child's cached member in `FillSources` (the existing `dr_sp_->source_slot()` access pattern, now `dr_sp_->index_tp_delta_`):

| Species | Sources | Note |
|---|---|---|
| `NCDMSpecies` (×N) | δ_ncdm, θ_ncdm | each instance its own pair; the `index_tp_*_ncdm1_ + slot` base-block is gone |
| `DCDM_DR` | δ_dcdm, θ_dcdm (dcdm channel) + δ_dr, θ_dr (dr child) | composite forwards to its dcdm/dr children |
| `DNCDM_DR` | δ_dr, θ_dr (dr child) | the dncdm child **registers a slot it never writes** (no `FillSources`), preserving today's behavior where DNCDM consumes an NCDM slot silently |

Note the DCDM column/source split: the standalone `DCDMSpecies` child owns the `d_dcdm`/`t_dcdm` **column** (`dcdm.cpp` `WriteOutputColumns`) while the `DCDM_DR` composite owns the δ_dcdm **source value** (`dcdm_dr_species.cpp` `FillSources`). Both must address the same slot, so the δ_dcdm/θ_dcdm members live on the child and the composite reads them through its child pointer — consistent with how the DR child slot is already shared.

### 6. Deleted vs. kept

**Deleted from `perturbations_module.h` / `.cpp`:**
- every per-species `index_tp_delta_X_` / `index_tp_theta_X_` member;
- every per-species `has_source_delta_X_` / `has_source_theta_X_` flag;
- the stringly-typed `all_species_.count("…")` source-gating block (`:957-1013`);
- the per-species `class_define_index` block (`:1075-1107`);
- both `SetSourceSlot` `dynamic_cast` loops (`:639-664`);
- `SetSourceSlot` / `source_slot_` / `source_slot()` on `NCDMSpecies` and `DarkRadiationSpecies`.

`NcdmFamily(all_species_).size()` and `DrSpeciesCount(...)` are no longer needed to **size** the source block (each instance self-registers), removing those call sites from `perturb_indices_of_perturbs`.

**Kept as members** (module-owned; transfer reads them by name): `index_tp_t0_/t1_/t2_/p_`, `index_tp_delta_m_/cb_/tot_`, `index_tp_theta_m_/cb_/tot_`, and all metric/potential types (`phi`, `phi_prime`, `phi_plus_psi`, `psi`, `h`, `h_prime`, `eta`, `eta_prime`, `H_T_Nb_prime`, `k2gamma_Nb`) plus their `has_source_*` flags.

## Behavior preservation & testing

Tk/vTk **column order** is already produced by the species-iteration column writer (`perturbations_module.cpp:290-389`), not by the `class_define_index` order, so reordering the internal `index_tp` allocation does not move output columns. Output is therefore expected **byte-identical**; the internal source-table column order changes but is never observed outside the module.

Acceptance:
- Full scenario suite with `classyref` installed, `TEST_LEVEL=2 -m test_scenario`: all pass.
- Named multi-species guards green: `dcdm_dr`, `tensor_massive_ncdm`, `dncdm_dr_computes`, plus NCDM and IDM_DR / IDM_DRMD scenarios.
- `COMPARE_OUTPUT_REF=1` against the master baseline: no output deltas (or, if any appear, explained at the ~0.1 % tolerance with explicit column-name comparison — never blind max-rel-diff, per project guidance on `Cl^TE` zero-crossings).
- Clean `pip install . --no-build-isolation` with no new warnings.

## Out of scope (sibling follow-ups)

- **`background_indices()` cleanup** — the background analog of this issue: it suffers the same stringly `count()/at()` gating and NCDM/DNCDM `dynamic_cast` loops, blocked by the same shared base-block addressing (`index_bg_number_ncdm1_`, `index_bg_phi_scf_`). Applying this clean-loop dissolve there is a separate change.
- **Residual `HasNcdm` / `NcdmFamily` family-detection casts** still tracked by #308 — they serve N_eff / budget accounting, not source sizing, and are untouched here.
