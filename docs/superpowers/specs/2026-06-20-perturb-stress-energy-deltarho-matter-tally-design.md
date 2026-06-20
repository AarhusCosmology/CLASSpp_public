# Perturbations stress-energy: `DeltaRho`/`RhoPlusPTheta` interface + generic matter tally

Date: 2026-06-20
Issue: none filed (v4 cleanup programme; follow-on to #325 dispatch cache)

## Problem

`PerturbationsModule::perturb_total_stress_energy` (the scalar block,
`perturbations_module.cpp:4795–5030`, ~233 lines) is a hot function — it runs
inside every RHS evaluation via `perturb_einstein` — yet it is large, repeats
work, and does per-RHS type/string dispatch:

1. **Photons and baryons are special-cased** with a dedicated arithmetic block
   that initialises `delta_rho / rho_plus_p_theta / rho_plus_p_shear / delta_p /
   rho_plus_p_tot`, duplicating what a generic species loop already does.
2. **Two passes over `all_species_`.** Pass 1 accumulates `rho_plus_p_shear`
   alone, skipping photons/baryons by `sp->name()` **string compare**; pass 2
   accumulates everything else.
3. **The matter tally double-computes.** For every matter species it calls
   `MatterRhoDelta` / `MatterRhoPlusPTheta`, which recompute `Delta` / `Theta`
   a **second time** (already computed for the Einstein totals), gated per
   species by `IsMatterSpecies()` / `IsColdMatterSpecies()` virtual calls.
4. **A per-RHS string map lookup** (`all_species_.count("IDM_DR_IDR")`,
   `4870`) to set `scalar_ctx.idr_nature`, whose value is constant per run.

It also surfaced a **physics bug** (see Key findings #4).

### Out of scope

Tensor-mode `dynamic_cast<NCDMSpecies>` / `<DNCDM_DR_Species>` picking
(`5108–5119`) is a separate cleanup, not addressed here. The vector block is
unchanged.

## Key findings (drive the design)

1. **`Delta()` is almost always `ρ·δ`.** Every external consumer — both
   stress-energy loops, the IC `delta_tot` loop (`3960`), the matter tally, and
   all four composites' internal combination — immediately multiplies `Delta`
   by `Rho`. The *only* fractional-δ consumer is NCDM source output
   (`ncdm_species.cpp:285`). So the natural interface primitive is `DeltaRho()`
   (`= ρδ`), with the rare fractional consumer computing `DeltaRho()/Rho()`
   inline. (User: "Just do DeltaRho()/Rho(): this is completely logical.")

2. **`Theta()` is the same story for momentum.** The Einstein equations need
   `(ρ+P)θ`, not `θ`. Every consumer multiplies `Theta` by `ρ+P`; the only
   fractional-θ consumer is NCDM source output (`ncdm_species.cpp:297`). So the
   primitive is `RhoPlusPTheta()` (`= (ρ+P)θ`), with the rare consumer computing
   `RhoPlusPTheta()/(Rho()+P())` inline.

3. **The matter tally needs no per-species `Matter*` machinery.** The matter
   density can be extracted from quantities already computed in the loop using
   the **`ρ − 3P`** identity the *background* already uses
   (`background_module.cpp:344`, `rho_m += rho - 3p` for `EnergyType::Other`):
   - background matter density:  `ρ − 3P`     (0 for radiation: `ρ − 3·(ρ/3)`)
   - matter density perturbation: `δρ − 3δP`  (0 for radiation: `δρ − 3·(δρ/3)`)

   `δP` is exactly `DeltaP`, already computed in the loop, and `δρ` is
   `DeltaRho` — so the tally is **pure arithmetic on already-computed values,
   zero extra virtual calls.** This deletes the entire `Matter*` family
   (`MatterRho`, `MatterRhoDelta`, `MatterRhoPlusP`, `MatterRhoPlusPTheta`,
   their out-of-line defs in `base_species.cpp`, and all four composite
   overrides).

   - **Exactness.** `δρ − 3δP` recovers the matter contrast *exactly* when the
     matter part is pressureless (`δP=0`) and the rest is radiation-like
     (`δP=δρ/3`). This holds for cold species and for the **DCDM_DR** composite
     (cold dcdm + dark radiation ⇒ `δρ−3δP = δρ_dcdm` exactly). It is an
     approximation for **warm** species (standalone NCDM/DNCDM and the
     DNCDM_DR composite), where `δP` is an independent variable — there the
     contrast is pressure-subtracted. This is accepted intentionally: `delta_m`
     is an inherently lossy, gauge-dependent proxy (structure formation is
     sourced by the full `δρ`/potential, in comoving gauge), and the new
     definition makes the perturbation tally **consistent with the background**.
   - **`delta_cb` (cold tally) is unchanged**: cold species have `P=0`, so
     `ρ−3P=ρ` and `δρ−3δP=δρ`.

4. **The matter tally is real and load-bearing, and currently has a bug.** The
   four tallied quantities each feed an output (gated by `has_source_delta_m_ ||
   has_source_theta_m_`, i.e. only when P(k)/transfers are requested):

   | quantity  | set   | source slot              | consumer |
   |-----------|-------|--------------------------|----------|
   | `delta_m` | 5018  | `index_tp_delta_m_`      | T_m (`transfer:1055`), **P_m(k)** (`nonlinear:1579`), nl-corrections |
   | `delta_cb`| 4990  | `index_tp_delta_cb_`     | **P_cb(k)** (`nonlinear:1582`) |
   | `theta_m` | 5027  | `index_tp_theta_m_`      | velocity/RSD transfers (`transfer:1058/61/64`) |
   | `theta_cb`| 4993  | `index_tp_theta_cb_`     | cb velocity transfer (`transfer:444/465`) |

   `theta_m` is needed even for `delta_m`-only output: the gauge-invariant
   `delta_m` is built from it at `4718` (`delta_m += 3aH·theta_m/k²`).

   **The bug:** IDM_DR overrides `IsMatterSpecies()→false`
   (`idm_dr.h:123`, a documented asymmetry), so it is excluded from `delta_m` —
   while the background `rho_m` *includes* it (`EnergyType::Matter`). In a
   cosmology where IDM_DR is the bulk of the dark matter, `P_m(k)` is built from
   the *subdominant* matter only, with the wrong normalization. The new generic
   rule (every non-dark-energy species contributes `δρ−3δP`) has no
   `IsMatterSpecies` hook to special-case on, so **IDM_DR enters the tally
   naturally and the bug is fixed.** (User: "it is time to fix that asymmetry.")

5. **Photon shear / baryon pressure context is photon/baryon physics.** The
   `scalar_ctx.shear_g` (TCA/RSA-corrected photon shear) and
   `scalar_ctx.delta_p_b_over_rho_b` prep blocks are only consumed by
   `PhotonsSpecies::RhoPlusPShear` and `BaryonsSpecies::DeltaP` respectively
   (verified by grep; no other reader, not read in `PerturbDerivs`, and the
   source-filling `delta_p_b_over_rho_b` at `6006` is an independent local in
   the RHS). They can move into those species and the two `scalar_ctx` fields be
   deleted, making the module loop fully generic. Photons/baryons already reach
   thermo via `GetThermodynamicsModule()` and gauge via `GetPerturbs()`; they
   gain stored `thm_`/`ppt_` pointers for the no-`ppaw` methods.

6. **`PhotonsSpecies::RhoPlusPShear` is currently dead in the scalar path**
   (pass 1 skips photons; the live value is the inline `4904`). Unifying the
   loop makes it live — which is exactly why finding #5's shear formula must
   live there and reproduce the inline value in every approximation branch.

7. **PPF fluid is not in `active_species`** (its δ/θ slots are unregistered),
   yet its background `ρ+P` must still enter `rho_plus_p_tot` — today picked up
   because the loop iterates all of `all_species_`. The new loop iterates
   `active_species`, so the PPF block must add `f->Rho()+f->P()` explicitly.
   Lambda (the only other non-active species in play) has `ρ+P=0`.

## Goals

- Reduce the scalar block by ~80 % and remove all per-RHS string/name/type
  dispatch and the second `Delta`/`Theta` computation.
- One lex-order loop over `active_species` (from #325), with the matter tally
  fused in via cached booleans and already-computed values.
- Replace the `Delta`/`Theta` virtuals with `DeltaRho`/`RhoPlusPTheta`; delete
  the `Matter*` family and the `IsMatterSpecies` override chain.
- Fix the IDM_DR P(k) bug as a natural consequence.

## Design

### 1. `BaseSpecies` interface

```cpp
// Pure virtuals — absolute stress-energy contributions:
virtual double DeltaRho(const PerturbLayout&, const perturb_vector*, const double* y,
                        const double* pvecback, const perturb_workspace*) const = 0;   // = ρδ
virtual double RhoPlusPTheta(const PerturbLayout&, const perturb_vector*, const double* y,
                             const double* pvecback, const perturb_workspace*) const = 0; // = (ρ+P)θ
// DeltaP, RhoPlusPShear: unchanged (already absolute).
```

- **Removed:** `Delta()` and `Theta()` virtuals. The two fractional consumers
  (NCDM `FillSources`) compute `DeltaRho()/Rho()` and
  `RhoPlusPTheta()/(Rho()+P())` inline.
- **Removed:** `MatterRho`, `MatterRhoDelta`, `MatterRhoPlusP`,
  `MatterRhoPlusPTheta` (+ `base_species.cpp` out-of-line defs + composite
  overrides in `composite_species.cpp`, `dcdm_dr_species.cpp`,
  `dncdm_dr_species.cpp`, etc.).
- **Removed:** `IsMatterSpecies()` (base default + overrides in
  `ncdm_base_species.h:77`, `idm_dr.h:123`, any composite) — after verifying the
  matter tally is its only remaining consumer.
- **Kept, redefined:** `IsColdMatterSpecies()` default becomes
  `energy_type_ == EnergyType::Matter` (was `IsMatterSpecies()`). NCDM/DNCDM
  override `→ false` stays; composite override (child scan) stays. Used once per
  pv to set the cold/warm cache bit.

Per-species conversions are mechanical, e.g.:
```cpp
// CDM:      DeltaRho = pvecback[index_bg_rho_] * y[layout.idx_delta];
//           RhoPlusPTheta = pvecback[index_bg_rho_] * y[layout.idx_theta];   // P=0
// Photons:  DeltaRho = (idx_delta>=0) ? rho_g*y[idx_delta] : 0.;
//           RhoPlusPTheta = (idx_theta>=0) ? 4./3.*rho_g*y[idx_theta] : 0.;
// Composite DCDM_DR (and the other three): the guarded weighted-average Delta/Theta
//           overrides collapse to sums:
//           DeltaRho      = dcdm_->DeltaRho(my.dcdm,…) + dr_sp_->DeltaRho(my.dr,…);
//           RhoPlusPTheta = dcdm_->RhoPlusPTheta(my.dcdm,…) + dr_sp_->RhoPlusPTheta(my.dr,…);
```

### 2. Photon/baryon context into the species (finding #5/#6)

Delete `scalar_ctx.shear_g` and `scalar_ctx.delta_p_b_over_rho_b`. Add stored
`thm_`/`ppt_` to `PhotonsSpecies` and `BaryonsSpecies`
(`SetThermodynamicsModule`/`SetPerturbs` overrides — verify the module calls
both on every species).

```cpp
double PhotonsSpecies::RhoPlusPShear(const PerturbLayout& L, …, const double* pvecback,
                                     const perturb_workspace* ppw) const {
  const double rho_g = pvecback[index_bg_rho_];
  if (L.idx_shear >= 0) return 4./3.*rho_g*y[L.idx_shear];          // no approximation
  if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) return 0.;    // RSA: shear neglected
  if (ppt_->gauge == possible_gauges::newtonian)                    // TCA, 1st order
    return 4./3.*rho_g*(16./45./ppw->pvecthermo[thm_->index_th_dkappa_]*y[L.idx_theta]);
  return 0.;  // sync TCA: re-set later in perturb_einstein (unchanged)
}

double BaryonsSpecies::DeltaP(const PerturbLayout& L, const perturb_vector* pv, const double* y,
                              const double* pvecback, const perturb_workspace* ppw) const {
  const double rho_b = pvecback[index_bg_rho_];
  if (ppt_->has_perturbed_recombination && ppw->approx[ppw->index_ap_tca]==(int)tca_off)
    return rho_b*ppw->pvecthermo[thm_->index_th_wb_]*
           (y[L.idx_delta]+y[pv->index_pt_perturbed_recombination_delta_temp]);
  return rho_b*ppw->pvecthermo[thm_->index_th_cb2_]*y[L.idx_delta];
}
```

Values are byte-for-byte the formulas relocated from the module prep block.
`idr_nature` becomes a module member `idr_nature_` resolved once at construction.

### 3. The cache (extend `active_species`)

```cpp
struct ActiveSpecies {
  BaseSpecies*                species;
  BaseSpecies::PerturbLayout* layout;
  bool not_dark_energy;   // energy_type() != EnergyType::DarkEnergy   (tally inclusion)
  bool is_cold;           // IsColdMatterSpecies()                      (delta_cb bucket)
};
```

Filled once per pv where `active_species` is built
(`perturbations_module.cpp:2968` and the vector/tensor twins). Radiation is
included by `not_dark_energy` but contributes 0 via `ρ−3P`/`δρ−3δP`, so its
bucket assignment is irrelevant.

### 4. The new scalar block

```cpp
if (_scalars_) {
  const double* pb = ppw->pvecback.data();
  const perturb_vector* pv = ppw->pv.get();

  ppw->scalar_ctx.k=k; ppw->scalar_ctx.k2=k2; ppw->scalar_ctx.a=a; ppw->scalar_ctx.a2=a2;
  ppw->scalar_ctx.gauge = static_cast<int>(ppt->gauge);
  ppw->scalar_ctx.idr_nature = idr_nature_;          // cached — no map lookup

  ppw->delta_rho=0.; ppw->rho_plus_p_theta=0.; ppw->rho_plus_p_shear=0.;
  ppw->delta_p=0.; ppw->rho_plus_p_tot=0.;
  struct Tally { double drho=0,rho=0,rppt=0,rpp=0; } cold, warm;
  const bool tally = has_source_delta_m_ || has_source_theta_m_;

  for (const auto& e : pv->active_species) {
    const auto& L     = *e.layout;
    const double rho  = e.species->Rho(pb);
    const double p    = e.species->P(pb);
    const double rpp  = rho + p;
    const double drho = e.species->DeltaRho(L, pv, y, pb, ppw);
    const double dp   = e.species->DeltaP(L, pv, y, pb, ppw);
    const double rppt = e.species->RhoPlusPTheta(L, pv, y, pb, ppw);

    ppw->delta_rho        += drho;
    ppw->rho_plus_p_theta += rppt;
    ppw->delta_p          += dp;
    ppw->rho_plus_p_shear += e.species->RhoPlusPShear(L, pv, y, pb, ppw);
    ppw->rho_plus_p_tot   += rpp;

    if (tally && e.not_dark_energy) {
      const double m_rho  = rho - 3.*p;                       // 0 for radiation
      const double m_drho = drho - 3.*dp;                     // 0 for radiation
      const double m_rppt = (rpp > 0.) ? m_rho*rppt/rpp : 0.; // ≈ (ρ−3P)·θ
      Tally& t = e.is_cold ? cold : warm;
      t.drho += m_drho; t.rho += m_rho; t.rppt += m_rppt; t.rpp += m_rho;
    }
  }

  if (has_source_delta_m_ && has_source_delta_cb_) ppw->delta_cb = cold.drho/cold.rho;
  if ((has_source_delta_m_||has_source_theta_m_)&&(has_source_delta_cb_||has_source_theta_cb_))
    ppw->theta_cb = cold.rppt/cold.rpp;

  if (auto* f = ppf_fluid()) {                          // module-owned closure, last
    f->ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
    ppw->delta_rho        += ppw->delta_rho_fld;
    ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
    ppw->delta_p          += ppw->delta_p_fld;
    ppw->rho_plus_p_tot   += f->Rho(pb) + f->P(pb);     // PPF excluded from active_species
  }

  if (has_source_delta_m_)
    ppw->delta_m = (cold.drho+warm.drho)/(cold.rho+warm.rho);
  if (has_source_delta_m_||has_source_theta_m_)
    ppw->theta_m = (cold.rppt+warm.rppt)/(cold.rpp+warm.rpp);
}
```

~50 lines vs 233. Each species' `DeltaRho`/`RhoPlusPTheta`/`DeltaP`/`Rho`/`P` is
computed **once**; the matter tally is arithmetic on those values.

### 5. Callers outside the function

- IC `delta_tot` loop (`3960–3978`): `rho*sp->Delta` → `sp->DeltaRho`;
  `rho_plus_p*sp->Theta` → `sp->RhoPlusPTheta` (still skips `DarkEnergy`).
- NCDM `FillSources` (`ncdm_species.cpp:285/297`): fractional δ/θ via
  `DeltaRho()/Rho()` and `RhoPlusPTheta()/(Rho()+P())`.

## Behavior changes (intentional — regenerate reference, document)

- **Warm NCDM/DNCDM `delta_m`/`theta_m`** become pressure-subtracted (`ρ−3P`
  weighting), now consistent with the background. Magnitude ~`3w·f_ν`
  (~0.001–0.1 % at z=0; larger at high-z mP(k)).
- **IDM_DR** now enters `delta_m`/`delta_cb` (bug fix); large change for
  IDM_DR-dominated cosmologies.
- **`delta_cb` for purely-cold species, and all non-matter-output runs,
  unchanged.**

## Correctness argument

- Totals (`delta_rho/rho_plus_p_theta/delta_p/rho_plus_p_shear/rho_plus_p_tot`)
  are the same sums as today, reordered into lex order over `active_species`
  (photons/baryons no longer summed first). Pure-lex reordering perturbs FP
  reductions at the ULP level (amplified by the ODE) — accepted under the
  ≤0.1 % bar ([[feedback_no_bit_identical_requirement]],
  [[feedback_vectorization_reduction_drift]]).
- Photon/baryon become ordinary loop members: `ρ+P=4/3ρ_g`/`ρ_b`,
  `RhoPlusPShear` reproduces the old inline `4904` in every approximation
  branch, baryon shear = 0, `DeltaP` reproduces `4905`.
- PPF `ρ+P` re-added explicitly; Lambda (`ρ+P=0`) and other non-active species
  carry no `rho_plus_p_tot` — **verify** no non-active species has `ρ+P≠0`
  besides PPF.
- Matter tally: `ρ−3P`/`δρ−3δP` is exact for cold and DCDM_DR; the warm/IDM_DR
  changes are the intended behavior changes above.

## Verification

1. Build clean (CMake) Debug (asserts) + Release.
2. 4-scenario diff — PPF fluid / non-PPF fluid / massive NCDM / tensors —
   TT & mP(k) ≤0.1 % (TE zero-crossing handled, never blind max-rel-diff).
3. Dedicated IDM_DR and DNCDM scenarios: confirm the `delta_m`/P(k) change is
   the expected direction/magnitude (IDM_DR now *included*), not a regression.
4. Full `TEST_LEVEL=2` against regenerated `classyref`; regenerate goldens.
5. Re-profile `class_profiled` (Planck-2018 single-NCDM) vs the #325 baseline.

## Risks

- **Photon/baryon `thm_`/`ppt_` not wired** → assert non-null in the new
  methods; verify the module calls `SetThermodynamicsModule`/`SetPerturbs` on
  every species.
- **`IsMatterSpecies` has a non-tally consumer** → grep before deletion; if any
  remains, keep it but stop using it in the tally.
- **`rho_plus_p_tot` drops a non-active species with `ρ+P≠0`** → covered by the
  PPF + non-PPF-fluid scenarios and the explicit verification step.
- **Composite `IsColdMatterSpecies` misclassifies the cold/warm bucket**
  (DCDM_DR→cold, DNCDM_DR→warm) → verify against the existing child-scan
  override.
