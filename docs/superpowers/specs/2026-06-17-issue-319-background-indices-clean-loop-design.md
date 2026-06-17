# Issue #319 — clean-loop `background_indices()` + a post-accumulation hook

**Date:** 2026-06-17
**Issue:** [#319](https://github.com/AarhusCosmology/CLASSpp/issues/319) — apply the #309 clean-loop treatment to `BackgroundModule::background_indices()`
**Tech Stack:** C++17 (`source/background_module.*`, species plugins in `species/`); Python `classy` Cython wrapper; `pytest`; build via `pip install . --no-build-isolation`; reference comparison via `classyref` + the `TEST_LEVEL`/`COMPARE_OUTPUT_REF` env knobs in `python/test_class.py`.

## Goal

`BackgroundModule::background_indices()` is the background-index analog of #309's transfer-source registration, but it never became a flat `for (auto& [name, sp] : all_species_) sp->RegisterBackgroundIndices(index_bg);` loop. Instead it is a ~115-line sequence of stringly-typed presence gates (`if (all_species_.count("CDM")) all_species_.at("CDM")->RegisterBackgroundIndices(...)`), two `dynamic_cast` loops (NCDM, DNCDM_DR), and module-owned indices positionally interleaved with specific species.

This issue dissolves all of that, applying the same treatment #309 achieved for transfer sources: **every species registers its own background indices in a single flat loop**, module-owned aggregates live in a head/tail around that loop, and the structural reason a species index has stayed module-side (the DRMD `G_over_aH`) is fixed at the root rather than papered over.

`RegisterBackgroundIndices(int& index_bg)` is *already* the per-species virtual (pure virtual on `BaseSpecies`, implemented by every species). So the core work is entirely at the call site; the species-side registration interface already exists.

## Investigation: why these indices stayed on the module

The per-step background fill is `BackgroundModule::background_functions()`. It runs in two implicit phases:

1. **Pre-Friedmann (accumulation).** A flat loop calls `sp->ComputeBackground(a_rel, pvecback_B, pvecback)` for every species; each writes its own density/pressure slots. An `accumulate` lambda sums `rho_tot`, `p_tot`, `rho_r`, `rho_m`, and `dp_dloga`.
2. **Friedmann.** `H = sqrt(rho_tot - K/a²)` (the *only* place the Friedmann equation is assumed).

`G_over_aH_drmd = Gint / (H · a_rel)` needs `H`, which is only known *after* the accumulation loop. There is **no post-Friedmann (post-H) species hook**, so the write was parked in a hardcoded module block after `H` (line 402-411 in the current tree) and the index was left on the module next to the write. The scalar field hits the identical wall: its `p_tot_prime` contribution also needs `H`, so it too is a hardcoded post-H module block (`ComputePPrimeAndWrite`, line 396-400). The `dp_dloga` accumulator exists *solely* to defer the `a·H` multiplication out of the accumulation phase — its only consumer is `p_tot_prime = a·H·dp_dloga`.

So the root cause is a missing post-accumulation lifecycle hook. The codebase already has phase-ordering hooks for the analogous pre-H case (`RequiresDeferredBackground()` for Fluid, `RequiresDeferredPerturbDerivs()` for PPF fluid), so adding a post-accumulation hook fits the existing design vocabulary.

**Precedent for relocating physics into species:** `FluidSpecies::ComputeWFld` and `FluidSpecies::ComputePpf` already hold what used to be `BackgroundModule::background_w_fld` / PPF code, with thin delegating module wrappers kept for cross-module callers (thermodynamics, perturbations, nonlinear). The DRMD relocation in this issue mirrors that exactly.

## Design

### A. Flat registration loop (core #309 treatment)

`background_indices()` becomes:

```
index_bg = 0
// head (module-owned, always)
a, H, H_prime          → bg_size_short_ = index_bg
// single flat loop — every species owns its registration
for (auto& [name, sp] : all_species_) sp->RegisterBackgroundIndices(index_bg);
// tail (module-owned aggregates)
rho_tot, p_tot, p_tot_prime, Omega_r   → bg_size_normal_ = index_bg
// long vector (module-owned)
rho_crit, Omega_m, conf_distance, ang_distance, lum_distance, time, rs, D, f, …
```

The stringly-typed gates for CDM / IDM_DRMD_IDR_DRMD / ScalarField / Lambda / Fluid / UR / IDM_DR_IDR and both `dynamic_cast` loops (NCDM, DNCDM_DR) all collapse into the single loop. Photons/baryons are no longer special-cased to the front; they register in lexicographic position like everything else.

**Index *order* in the bg-vector changes** (now lexicographic by species name). All access is by named index (`pvecback[index_bg_*]`), so model output is unaffected. The ODE state vector uses a *separate* running index (`index_bi`, unchanged by this issue).

### B. Dead-index deletions

These module members are `class_define_index`'d (reserving a `pvecback` column) but never written *and* never read — delete them and their declarations:

- `index_bg_number_ncdm1_` — the NCDM species own and write their own `index_bg_number_`; the cached module offset and its "others contiguous" comment are a vestige of the old C layout, read nowhere.
- `index_bg_pseudo_p_ncdm1_` — same.
- `index_bg_Gamma0_drmd_` — a reserved bg-table column, never written or read. (Not to be confused with the live scalar `double Gamma0_drmd_`, which is the DRMD interaction rate and is handled in §E.)

### C. Scalar-field index ownership

`ScalarFieldSpecies` already registers and owns `index_bg_phi_scf_`, `index_bg_phi_prime_scf_`, `index_bg_V_scf_`, `index_bg_dV_scf_`, `index_bg_ddV_scf_`, `index_bg_p_prime_scf_` in its own `RegisterBackgroundIndices`/members. The `BackgroundModule` keeps duplicate copies of the first five purely so `scalar_field.cpp` can read them via `bgm->index_bg_*_scf_`. **Delete the five module copies** and re-route every `bgm->index_bg_*_scf_` / `mod.GetBackgroundModule()->index_bg_*_scf_` read in `scalar_field.cpp` to the species' own members (all such reads are inside `ScalarFieldSpecies` methods, so `this->` is in scope).

### D. Two post-accumulation operations, cleanly split

The single post-H module block is replaced by a post-accumulation pass over the flat species map. The two operations are kept distinct (no return-value bolted onto a writer hook):

**D1. `PPrime()` replaces `DpDloga()` — a post-H pressure-derivative accumulator.**

- Delete the `DpDloga()` pure-virtual on `BaseSpecies` and all ~20 overrides, the `dp_dloga` local, and the `p_tot_prime = a·H·dp_dloga` formula plus the scf `ComputePPrimeAndWrite` exception block.
- Add to `BaseSpecies`:
  ```cpp
  // Post-Friedmann. This species' conformal-time pressure derivative p'.
  // pvecback_B = ODE state (IN). Default 0 (matter / Lambda).
  virtual double PPrime(double a, double H,
                        const double* pvecback_B, const double* pvecback) const { return 0.; }
  ```
  Default `0.` lets the matter/Λ species (CDM, baryons, Lambda, IDM_DR, IDM_DRMD, DCDM) drop their override entirely. Radiation/NCDM/Fluid/DR/IDR override returning `a·H·(old DpDloga value)`. `CompositeSpecies::PPrime` sums children.
- **scf folds in like every other species:** its `PPrime` computes `phi_prime·(-phi_prime·H/a - 2/3·dV)` reading `phi_prime` from `pvecback_B[index_bi_phi_prime_scf_]` and `dV` via `dV_scf(phi)` — no special-casing in the module.
- `background_functions` pre-H accumulate loop sums only `rho_tot`/`p_tot`/`rho_r`/`rho_m`. After `H` is computed, a post-accumulation loop does `for (sp) p_tot_prime += sp->PPrime(a, H, pvecback_B, pvecback);` then writes `pvecback[index_bg_p_tot_prime_] = p_tot_prime`.

**D2. `FinalizeBackground()` — void; species write their own post-H *output* slots.**

```cpp
// Post-Friedmann. Write this species' H-dependent owned slots. pvecback = OUT.
virtual void FinalizeBackground(double a, double H,
                                const double* pvecback_B, double* pvecback) {}
```
- `ScalarFieldSpecies::FinalizeBackground` writes its `index_bg_p_prime_scf_` output column (the same `p'` expression as its `PPrime`; the recomputation is three multiplies).
- `IDM_DRMD_IDR_DRMD_Species::FinalizeBackground` writes its now-species-owned `index_bg_G_over_aH_drmd_` = `Gint / (H · a_rel)`.
- Runs in the same post-accumulation loop pass as D1.

**`pvecback` discipline.** `H` is passed explicitly so neither hook re-reads it from `pvecback`. scf reads its ODE inputs from `pvecback_B` (the input vector), keeping `pvecback` write-only for scf. *Honesty note:* DRMD's `Gint` needs the idm/idr density ratio, which lives in `pvecback` from the accumulation pass, so `FinalizeBackground` reads those two cross-species densities — that is reading other species' finalized inputs, not re-reading its own writes. Kept and flagged rather than pretended pure.

### E. DRMD physics relocation (mirrors Fluid)

`index_bg_G_over_aH_drmd_` moves onto `IDM_DRMD_IDR_DRMD_Species` (the composite is the natural owner — it is a property of the idm↔idr interaction and already exposes `idm_drmd()`/`idr_drmd()` children), registered in its `RegisterBackgroundIndices` and exposed via `bg_G_over_aH_index()`. To let the species *compute* the value:

- Move `background_idm_drmd`'s body (the `Rint`/`csp2`/`Gint` computation) into `IDM_DRMD_IDR_DRMD_Species::ComputeIdmDrmd(...)`, mirroring `FluidSpecies::ComputeWFld`.
- Move the IC-derived constants it depends on — `Gamma0_drmd_`, `z_stop_` (and the output diagnostic `f_idr_drmd_`) — into the species, computed in its `SetBackgroundInitialConditions` hook (it already exposes `z_stop()` / `G_over_aH_drmd()` input accessors).
- `BackgroundModule::background_idm_drmd` becomes a thin delegating wrapper (like `background_w_fld`), so the two `perturbations_module` callers and the species' own perturbation-side callers are untouched.

Module read/write sites that move:

- `background_functions` per-step write already has the `drmd` species reference in hand → moves into `FinalizeBackground` (§D2).
- `idm_drmd_idr_drmd_species.cpp` self-read of `bgm_->index_bg_G_over_aH_drmd_` becomes its own member.
- **The decoupling-redshift scan moves too** (§F1) — `z_dec_drmd_` and `G_over_aH_tmp_` are not module diagnostics, they are DRMD physics; they leave the module entirely.

### F1. `ProcessBackgroundTable()` — post-table-fill species analysis

The decoupling-redshift scan (`background_solve`, line 886-894) reads the *completed* `G_over_aH` column across all table rows to find where `G_over_aH ≈ 1`, storing `z_dec_drmd_`; `G_over_aH_tmp_` is its running scan-scratch and the only reader of `z_dec_drmd_` is a verbose printf. This is a table-scope, DRMD-specific reduction — it does not belong in the module's generic table loop. Add a third background lifecycle hook:

```cpp
// Called once after background_solve has filled the entire background table.
// Lets a species run table-scope analysis over its own columns.
virtual void ProcessBackgroundTable(const double* background_table, int n_rows,
                                    int row_stride, const double* z_table) {}
```

`IDM_DRMD_IDR_DRMD_Species::ProcessBackgroundTable` performs the scan over its own `bg_G_over_aH_index()` column and stores `z_dec_drmd_` / `G_over_aH_tmp_` as **species members**. The module calls `for (sp) sp->ProcessBackgroundTable(...)` once after the table is built (a separate, cheap single pass rather than piggybacking the existing loop — clean separation). The verbose banner reads `drmd.z_dec_drmd()` via accessor. Since `z_dec_drmd_` feeds only a printf, this relocation has **zero effect on numerical output** (trivially byte-identical).

## Components touched

| File | Change |
|------|--------|
| `source/background_module.cpp` | `background_indices()` → head + flat loop + tail; delete dead members; pre-H accumulate drops `dp_dloga`; add post-accumulation loop (`PPrime` sum + `FinalizeBackground`); add post-table-fill `ProcessBackgroundTable` loop; decoupling-redshift scan removed from the table loop; `background_idm_drmd` → thin wrapper; DRMD IC block delegates to species; verbose banner reads `drmd.z_dec_drmd()` |
| `source/background_module.h` | Remove `index_bg_number_ncdm1_`, `index_bg_pseudo_p_ncdm1_`, `index_bg_Gamma0_drmd_`, 5 `index_bg_*_scf_`, `z_dec_drmd_`, `G_over_aH_tmp_`; move `Gamma0_drmd_`/`z_stop_`/`f_idr_drmd_` ownership to species (keep wrapper decl) |
| `species/base_species.h` | Remove `DpDloga()` pure-virtual; add `PPrime()` (default 0), `FinalizeBackground()` and `ProcessBackgroundTable()` (default no-op) |
| `species/scalar_field.{h,cpp}` | Drop `DpDloga`; route bg-index reads to own members; implement `PPrime` (folds in old exception); `FinalizeBackground` writes `p_prime_scf`; delete `ComputePPrimeAndWrite` |
| `species/idm_drmd_idr_drmd_species.{h,cpp}` | Own `index_bg_G_over_aH_drmd_` + `bg_G_over_aH_index()`; `ComputeIdmDrmd` physics method; hold `Gamma0_drmd_`/`z_stop_`/`f_idr_drmd_` from `SetBackgroundInitialConditions`; `FinalizeBackground`; own `z_dec_drmd_`/`G_over_aH_tmp_` + `z_dec_drmd()` accessor; `ProcessBackgroundTable` decoupling scan |
| All other `species/*.{h,cpp}` | Remove `DpDloga` override; add `PPrime` override only where non-zero (UR, photons, ncdm, fluid, dr, idr, idr_drmd, dncdm, composite) |

## Verification

Mirrors #309. The local `classyref` reference was regenerated from current master on 2026-06-17, so `COMPARE_OUTPUT_REF` is a clean gate again.

- **Characterization golden test** (new, `python/`): `get_background()` columns/values across ~7 cosmologies (ΛCDM, ncdm, scf, idm_drmd, idm_dr, dncdm_dr, curvature). Assert **exact equality on every column except `p_tot_prime`**, which is checked at the project tolerance (see below).
- **Scenario suite:** `COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 pytest -m test_scenario test_class.py` against the fresh `classyref`. Expect all green (no newtonian stale-ref fails).
- **Targeted sync + newtonian** scenarios for scf, idm_drmd, ncdm.

**Bit-identical scope (honest):** the index relocations (A/B/C), `FinalizeBackground` writes, and the DRMD physics relocation are byte-identical. The `PPrime` refactor is **not** bit-identical for `p_tot_prime`: it distributes the `a·H` factor across the per-species sum (`Σ aH·dᵢ` vs `aH·Σ dᵢ`) and folds scf into the sum rather than adding it last — both ULP-level reorderings. Per project policy (bit-identical is the wrong bar; verify at ~0.1% and handle reductions carefully — see `feedback_no_bit_identical_requirement`, `feedback_vectorization_reduction_drift`), verify TT/P(k) stay <0.1% and regenerate the characterization baseline for `p_tot_prime` rather than asserting exact equality. Every other background column stays exactly equal.

## Out of scope

- Fully relocating `background_w_fld` / `background_idm_drmd` *callers* in other modules — the thin wrapper keeps them working; this issue only moves the physics body + IC constants.
- The residual #308 family-detection casts elsewhere (`HasNcdm`, etc.) — #308's domain. (This issue does dissolve the NCDM/DNCDM casts *in `background_indices()`*.)
- `#310` (species scratch out of `perturb_workspace`) — closely related but separate.
