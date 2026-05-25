# Newtonian-gauge transformation via species dispatch — design

**Date:** 2026-05-24
**Status:** Design (pending plan)
**Area:** `source/perturbations_module.cpp` — scalar-mode initial conditions, `if (ppt->gauge == newtonian)` block (currently lines 4633–4935).
**Context:** Continuation of the species-as-plugins refactoring (PR #268 follow-up family). The TODO at `perturbations_module.cpp:4919` explicitly anticipates this work and names the target hook `PerturbSynchronousToNewtonian`.

## Goal

Replace the ~280-line Newtonian-gauge initial-condition block — which currently performs ~15 `all_species_.count("X")` dispatches with concrete-type downcasts — with species-polymorphic dispatch, so that adding a new species requires **no** edits to this block. This is the species-picking anti-pattern the architecture is removing (module code must loop + dispatch, never downcast).

This is a research code, so we improve the design rather than mechanically relocating the existing blocks: we derive the gauge shift from the background continuity equation (a single source of truth shared with `BackgroundDerivs`), fix a latent quirk, and fix a confirmed sign bug.

## The unifying observation

The synchronous→Newtonian gauge transformation of a density contrast is purely *kinematic*: it depends only on how the species' **background** density evolves. For a time shift `α = (h' + 6η')/(2k²)` (Ma & Bertschinger), every fluid-like species transforms by the same two lines:

```
delta += (ρ̇/ρ) · alpha          // ρ̇ ≡ dρ̄/dτ  (conformal-time derivative of background density)
theta += k² · alpha              // shear, l3, and all higher moments are gauge-invariant
```

where `ρ̇/ρ = -3ℋ(1+w)` by the continuity equation (`ℋ ≡ a'/a = a·H`), with extra terms for non-standard species (e.g. decay).

This reproduces **every** inline block today (`ℋ = a_prime_over_a`):

| species | current inline shift | `ρ̇/ρ` | agree? |
|---|---|---|---|
| photons (w=⅓) | `δ -= 4ℋα` | `-4ℋ` | ✓ |
| baryons, cdm, idm_dr, idm_drmd (w=0) | `δ -= 3ℋα` | `-3ℋ` | ✓ |
| idr_drmd (w=⅓) | `δ -= 4ℋα` | `-4ℋ` | ✓ |
| dcdm (decaying) | `δ += (-3ℋ-aΓ)α` | `-3ℋ-aΓ` | ✓ |
| UR / IDR (radiation) | re-seed = `δ_sync − 4ℋα` | `-4ℋ` | ✓ |
| DR (decaying radiation) | re-seed = `(δ_sync + (-4ℋ+decay)α)·r_dr` | `-4ℋ` + decay sink | ✓ (override: decay term + `r_dr` rescale) |
| **fluid (non-PPF)** | `δ += 3(1+w)ℋα` | `-3(1+w)ℋ` | ✗ **opposite sign — bug** |

The DCDM case is the proof: `DCDMSpecies::BackgroundDerivs` (dcdm.cpp:37) computes `dρ/dτ = -a(3H+Γ)ρ`, i.e. `ρ̇/ρ = -3ℋ-aΓ`, which is *identical* to its inline gauge coefficient. The same physics that integrates the background also drives the gauge transform.

The default `ρ̇/ρ = -3ℋ(Rho+P)/Rho` reproduces the exact analytic derivative for every analytic species (photons `a⁻⁴` → `-4ℋ`; matter `a⁻³` → `-3ℋ`).

## What the block does today (three operations)

The Newtonian block decomposes into three logically distinct operations:

1. **Accumulate `alpha`** (4657–4724): build `delta_tot` / `velocity_tot` from per-species density/velocity, compute `alpha`, then set the metric `phi = eta − ℋα`.
2. **Gauge-shift each species' own variables** (4726–4838): the ~10 downcast blocks doing `delta -= n·ℋα`, `theta += k²α`.
3. **Re-seed the relativistic species** (4843–4934): UR/IDR/NCDM/DR get their `y[]` written from the shifted `delta_ur`/`theta_ur`/`shear_ur`/`l3_ur`/`delta_dr`.

Key fact enabling the redesign: `dispatch_species_ic()` (block b) runs `ApplyInitialConditions` for **all** species in **both** gauges, so the *synchronous* IC is already in `y[]` for UR/IDR/NCDM/DR by the time we reach the Newtonian block. Operation 3's "re-seed" is therefore just an in-place shift of values already present — identical to operation 2. (Every special species — NCDM, ScalarField, IDM_DR_IDR, DCDM_DR, DNCDM_DR — already has a layout-based `ApplyInitialConditions`, confirmed.)

## Architecture

### New per-species physics accessor

```cpp
// BaseSpecies — log conformal-time derivative of the background density, ρ̇/ρ.
// Single source of truth for both background evolution and the gauge transform.
// Default: continuity relation. Overridden by species with extra source/sink terms.
virtual double RhoDotOverRho(const double* pvecback, double a_prime_over_a) const {
  return -3. * a_prime_over_a * (Rho(pvecback) + P(pvecback)) / Rho(pvecback);
}
```

- **Override — DCDMSpecies:** `-3·a_prime_over_a − a·Γ` (the decay sink already in its `BackgroundDerivs`).
- (No other simple species needs to override; DR/NCDM handle their own transform — see below.)

`a_prime_over_a` is passed explicitly because `BaseSpecies` does not store the background `H`/`a` indices. Its sole caller is the base transform default.

### New per-species transform hook

```cpp
// BaseSpecies — transform this species' own perturbation variables from
// synchronous to Newtonian gauge in place. Called by the module AFTER alpha is
// known (ctx.alpha filled). Reuses PerturbIcContext (already reserves alpha/alpha_prime
// "filled by module AFTER the species loop for the Newtonian gauge transformation").
// Default: the universal fluid-like rule.
virtual void PerturbSynchronousToNewtonian(const PerturbLayout& layout,
                                           double* y,
                                           const PerturbIcContext& ctx) { /* universal δ/θ shift */ }
```

The base default applies `delta += RhoDotOverRho(...)·alpha` and `theta += k²·alpha` to the species' own `idx_delta`/`idx_theta` (guarded by `>= 0`). This covers: **photons, baryons, cdm, idm_dr, idm_drmd, idr_drmd, dcdm, UR, IDR, and fluid** (fluid now via the corrected sign — see Behavior changes).

**Overrides:**
- **NCDMBaseSpecies** (and DNCDM): per-q shift of the F0/F1 moments `∝ dln f0/dln q`; shear/l3 (F2/F3) gauge-invariant. Mirrors the current per-q block (4898–4914) but as an in-place shift of the already-seeded synchronous moments.
- **DarkRadiationSpecies:** decay-corrected `δ_dr` shift and the `r_dr` rescaling of the stored F-hierarchy (current 4922–4934). DR lives inside `DCDM_DR_Species` / `DNCDM_DR_Species`; the composite delegates.
- **ScalarFieldSpecies:** shifts `φ` and `φ′` using `index_bg_phi_prime_scf_` (current 4792–4810); not a fluid `δ`, so it does not use `RhoDotOverRho`.
- **Composites** (`DCDM_DR_Species`, `IDM_DR_IDR_Species`, `IDM_DRMD_IDR_DRMD_Species`, `DNCDM_DR_Species`): delegate to each member's hook (same delegation pattern as `CopyPerturbationsAcrossSwitch`). Members remain usable standalone.
  - **`IDM_DR_IDR_Species` requires care (LIVE in Newtonian gauge):** in tight coupling the idm_dr velocity is *locked* to the idr velocity rather than taking its own `k²α` (current 4870–4876: `theta_idm_dr = theta_ur`). The composite override reproduces this lock. `idr` itself is universal radiation. (The `class_test` gating idr_drmd does **not** gate idr.)

### Operation 1 — module-level `alpha` reduction (no downcasts)

The metric/`alpha` computation stays in the module (it owns the metric), but becomes a clean loop mirroring `perturb_total_stress_energy` Pass 2 (perturbations_module.cpp:5915–5916):

```
delta_rho_ic = 0;  rho_plus_p_theta_ic = 0;
for each species with energy_type != DarkEnergy:        // radiation + matter + "Other" (semi-rel NCDM)
    rho = sp.Rho(pvecback);  rho_plus_p = rho + sp.P(pvecback);
    delta_rho_ic        += rho        · sp.Delta(layout, pv, y, pvecback, ppw);   // reads synchronous y[]
    rho_plus_p_theta_ic += rho_plus_p · sp.Theta(layout, pv, y, pvecback, ppw);
delta_tot    = delta_rho_ic        / (rho_r + rho_m);
velocity_tot = rho_plus_p_theta_ic / (rho_r + rho_m);
alpha = (eta + 1.5·ℋ²/(k²·s2²)·(delta_tot + 3·ℋ/k²·velocity_tot)) / ℋ;
y[index_pt_phi] = eta − ℋ·alpha;
```

This **fixes the `fraccdm` quirk** automatically: today `fraccdm` is set only when `count("CDM")` is true (4365), so DCDM-/IDM_DR-only models silently drop cold matter from `alpha`. Summing `Rho·Delta` over all species includes them. The `(rho_r + rho_m)` denominator already equals the sum of all non-dark-energy `Rho` (the rho_r/rho_m construction loop at 4322–4350 splits "Other" into r+m, which sum back to `rho`), so the normalization is consistent.

**Ordering:** Operation 1 reads synchronous `y[]`, so it must run before the transform dispatch (which shifts `y[]`). The module computes `alpha`, sets `ctx.alpha` (and `ctx.alpha_prime`), then runs the second per-species loop calling `PerturbSynchronousToNewtonian`. The `delta_ur`/`theta_ur`/… locals and block (e) disappear.

## Data flow

```
block b: dispatch_species_ic() → ApplyInitialConditions writes SYNCHRONOUS IC into y[] (all species, both gauges)
            │
ppt->gauge == newtonian:
   op1 (module loop): Rho/P/Delta/Theta over species → delta_tot, velocity_tot → alpha → metric phi
            │  ctx.alpha = alpha;  ctx.alpha_prime = alpha_prime;
   op2+3 (dispatch):  for each species: sp->PerturbSynchronousToNewtonian(layout, y, ctx)
                         base default: δ += (ρ̇/ρ)α ; θ += k²α
                         overrides:    NCDM (per-q), DR (decay+rescale), ScalarField (φ,φ′), composites (delegate)
```

## Behavior changes (both confirmed)

1. **Fluid sign — bug fix.** The non-PPF fluid's `δ += 3(1+w)ℋα` is opposite to `ρ̇/ρ` and to every other species. Confirmed against the newest public CLASS, which uses the opposite (correct) sign; ours is the bug. The universal rule yields the correct `δ -= 3(1+w)ℋα`. Only affects the `use_ppf == _FALSE_` path (niche, which is why it went unnoticed).

2. **NCDM in `alpha` (minor).** Today NCDM enters `alpha` lumped into `fracnu·delta_ur` (the analytic shared relativistic IC). The redesign uses `Rho_ncdm·Delta_ncdm()` (the actual integrated moment ≈ `delta_ur`). Tiny numerical difference, arguably more correct, expected well within the ~0.1% tolerance suite.

## Error handling / guards

- The `class_test(pba->has_idr_drmd == _TRUE_)` at 4634 (Newtonian gauge not tested for the DRMD implementation) **stays** — it is a validation guard, not a dispatch guard. IDM_DRMD/IDR_DRMD branches remain unreachable in Newtonian gauge; their composite hook can still be implemented (universal) for consistency, but is dead here.
- All hook bodies guard each index with `idx >= 0` (species/mode may not have registered that slot), matching existing layout-based methods.

## Components touched

- `species/base_species.h` — add `RhoDotOverRho` and `PerturbSynchronousToNewtonian` virtuals (with universal default).
- `species/dcdm.{h,cpp}` — override `RhoDotOverRho`.
- `species/ncdm_species.{h,cpp}` — override `PerturbSynchronousToNewtonian` (per-q).
- `species/dark_radiation_species.{h,cpp}` — override (decay + `r_dr` rescale).
- `species/scalar_field.{h,cpp}` — override (φ, φ′).
- `species/dcdm_dr_species.{h,cpp}`, `species/idm_dr_idr_species.{h,cpp}`, `species/idm_drmd_idr_drmd_species.{h,cpp}`, `species/dncdm_dr_species.{h,cpp}` — composite delegation (IDM_DR_IDR also reproduces the idm_dr↔idr velocity lock).
- `source/perturbations_module.cpp` — replace the Newtonian block (4633–4935): op1 reduction loop + metric `phi` + `ctx.alpha` set + transform dispatch loop. Delete the inline shift blocks and block (e).

**No new files** → no build-system list changes (Makefile / setup.py / pbxproj unaffected; only existing files edited).

## Testing strategy

Behavior-preserving except the two documented changes; verified with the ~0.1% tolerance suite (atol+rtol, zero-crossing-aware — not bit-identical), per project convention.

- **Newtonian-gauge coverage across species** (the live ones): ΛCDM (photons/baryons/cdm/UR), NCDM, DCDM_DR, IDM_DR_IDR — run in Newtonian gauge, compare TT/EE/etc. to the synchronous-gauge reference within tolerance (cross-gauge agreement is the physical invariant).
- **Fluid sign fix:** a non-PPF fluid scenario (`use_ppf=no`, w0/wa) in Newtonian gauge — expected to *change* vs the old buggy output and to *agree* with public CLASS / with the synchronous-gauge result within tolerance. This is the acceptance test for change (1).
- **`fraccdm` fix:** a DCDM_DR-only (no literal CDM) Newtonian-gauge run — confirm `alpha` now includes the cold-matter contribution (cross-gauge agreement improves or holds).
- **Full 84-scenario grid** + DCDM strict baseline (`test_dcdm_dr_matches_reference`) to catch synchronous-gauge regressions (synchronous path must be untouched).
- TDD per species: assert the new hook reproduces the old inline shift for each species before deleting the inline block (a unit-level check that `RhoDotOverRho·α` equals the hand-written coefficient).
```
