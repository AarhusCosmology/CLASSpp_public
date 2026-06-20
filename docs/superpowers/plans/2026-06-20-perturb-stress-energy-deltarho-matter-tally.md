# perturb_total_stress_energy: DeltaRho/RhoPlusPTheta + generic matter tally — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse `PerturbationsModule::perturb_total_stress_energy` (scalar block) into one lex-order loop over `active_species`, replacing the `Delta()`/`Theta()` virtuals with the absolute primitives `DeltaRho()`/`RhoPlusPTheta()`, deleting the `Matter*` tally family in favour of a generic `(ρ−3P)`/`(δρ−3δP)` matter tally, and fixing the IDM_DR matter-tally exclusion bug.

**Architecture:** Four ordered changes, each independently buildable. (1) Interface rename `Delta→DeltaRho`, `Theta→RhoPlusPTheta` across `BaseSpecies` and all species, with consumers updated — byte-identical for ΛCDM, ULP for warm/scf (division cancels). (2) Move photon-shear/baryon-pressure context into the species; cache `idr_nature_`. (3) Rewrite the scalar block as one loop with a fused `(ρ−3P)` matter tally and an extended `active_species` cache. (4) Delete the now-dead `Matter*` family and `IsMatterSpecies`.

**Tech Stack:** C++17, CMake + scikit-build-core, `class` CLI + `classy`/`classyref` Python modules, `python/test_class.py` harness (`TEST_LEVEL`, `COMPARE_OUTPUT_REF`).

## Global Constraints

- **Never `git add -A`** in this repo (in-source CMake/Xcode artifacts leak in). Stage named files only.
- **Verification bar:** ≤0.1 % on TT and mP(k), with TE/EE zero-crossing handling — never a blind max-rel-diff. Byte-identical is *not* required ([[feedback_no_bit_identical_requirement]]).
- **Headers are C++-only:** no `#ifdef __cplusplus`, no `extern "C"` ([[feedback_cpp_only_no_c_guards]]).
- **`cclassy.pxd` is generated** — never hand-edit ([[feedback_cclassy_pxd_generated]]).
- **Module code loops + dispatches; no new downcasts.** Typed photon/baryon accessors are the only allowed exceptions ([[feedback_no_species_picking_in_modules]]).
- Branch already created: `perturb-stress-energy-deltarho-matter-tally` (spec committed there).
- Build (CLI): `cmake -S . -B build/cmake && cmake --build build/cmake --parallel` → `build/cmake/class`.
- Build (python): `pip install .` (rebuilds `classy`).
- Spec: `docs/superpowers/specs/2026-06-20-perturb-stress-energy-deltarho-matter-tally-design.md`.

---

## File Structure

| File | Responsibility | Tasks |
|------|----------------|-------|
| `species/base_species.h` | virtual interface: add `DeltaRho`/`RhoPlusPTheta`, remove `Delta`/`Theta`; (T4) remove `Matter*`/`IsMatterSpecies`, redefine `IsColdMatterSpecies` | 1, 4 |
| `species/base_species.cpp` | out-of-line `MatterRhoDelta`/`MatterRhoPlusPTheta` → use new primitives; (T4) delete | 1, 4 |
| `species/{cdm,baryons,photons,dcdm,ultra_relativistic,dark_radiation_species,fluid,scalar_field,ncdm_species,dncdm_species,interacting_species,lambda}.{h,cpp}` | per-species `DeltaRho`/`RhoPlusPTheta` bodies | 1 |
| `species/{dcdm_dr,dncdm_dr,idm_dr_idr,idm_drmd_idr_drmd}_species.cpp` | composite `DeltaRho`/`RhoPlusPTheta` = child sums; (T4) delete `MatterRho*` overrides | 1, 4 |
| `species/photons.{h,cpp}`, `species/baryons.{h,cpp}` | (T2) store `thm_`/`ppt_`; absorb shear / δp_b | 2 |
| `species/ncdm_base_species.h`, `species/idm_dr.h` | (T4) drop `IsMatterSpecies` override; keep/adjust `IsColdMatterSpecies` | 4 |
| `source/perturbations.h` | (T3) extend `ActiveSpecies` with `not_dark_energy`/`is_cold` | 3 |
| `source/perturbations_module.{h,cpp}` | call-site updates (T1); context deletion + `idr_nature_` (T2); function rewrite + cache fill (T3) | 1, 2, 3 |

---

## Task 1: Interface — `DeltaRho`/`RhoPlusPTheta` replace `Delta`/`Theta`

**Files:**
- Modify: `species/base_species.h`, `species/base_species.cpp`
- Modify: all species `.h`/`.cpp` (bodies below + table)
- Modify: `source/perturbations_module.cpp` (call sites), `species/ncdm_species.cpp` (FillSources)
- Test: `python/test_class.py` + CLI golden diff

**Interfaces:**
- Produces: `virtual double BaseSpecies::DeltaRho(const PerturbLayout&, const perturb_vector*, const double* y, const double* pvecback, const perturb_workspace*) const = 0;` (`= ρδ`) and `virtual double BaseSpecies::RhoPlusPTheta(...) const = 0;` (same signature, `= (ρ+P)θ`). `Delta`/`Theta` no longer exist.

- [ ] **Step 1: Edit `base_species.h` — swap the pure virtuals**

Replace the `Delta` declaration (`base_species.h:341-345`) and `Theta` declaration (`:347-352`) with:

```cpp
  /** Density contribution to the stress-energy tensor, δρ = ρ·δ. */
  virtual double DeltaRho(const PerturbLayout& layout,
                          const perturb_vector* pv,
                          const double* y,
                          const double* pvecback,
                          const perturb_workspace* ppw) const = 0;

  /** Momentum density (ρ+P)·θ. */
  virtual double RhoPlusPTheta(const PerturbLayout& layout,
                               const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const = 0;
```

- [ ] **Step 2: Edit `base_species.cpp` — matter helpers use the new primitives**

These are deleted in Task 4 but must compile now. Replace bodies:

```cpp
double BaseSpecies::MatterRhoDelta(const perturb_vector* pv, const double* y,
                                   const double* pvecback, const perturb_workspace* ppw) const {
  if (!IsMatterSpecies()) return 0.;
  return DeltaRho(*pv->species_layouts.at(collection_index_), pv, y, pvecback, ppw);
}
double BaseSpecies::MatterRhoPlusPTheta(const perturb_vector* pv, const double* y,
                                        const double* pvecback, const perturb_workspace* ppw) const {
  if (!IsMatterSpecies()) return 0.;
  return RhoPlusPTheta(*pv->species_layouts.at(collection_index_), pv, y, pvecback, ppw);
}
```

- [ ] **Step 3: Convert the trivial species (rename method + scale body)**

For each row, rename the method (declaration in `.h`, definition in `.cpp`) `Delta→DeltaRho`, `Theta→RhoPlusPTheta`, and replace the `return` expression. `pvb` ≡ the `pvecback` parameter (un-comment it where currently `/*pvecback*/`).

| Species (file) | `DeltaRho` returns | `RhoPlusPTheta` returns |
|---|---|---|
| CDM (`cdm.cpp`) | `pvb[index_bg_rho_] * y[layout.idx_delta]` | `(layout.idx_theta >= 0) ? pvb[index_bg_rho_] * y[layout.idx_theta] : 0.` |
| Baryons (`baryons.cpp`) | `pvb[index_bg_rho_] * y[layout.idx_delta]` | `pvb[index_bg_rho_] * y[layout.idx_theta]` |
| DCDM (`dcdm.cpp`) | `pvb[index_bg_rho_] * y[layout.idx_delta]` | `pvb[index_bg_rho_] * y[layout.idx_theta]` |
| IDM_DR (`interacting_species.cpp:52`) | `pvb[index_bg_rho_] * y[layout.idx_delta]` | `pvb[index_bg_rho_] * y[layout.idx_theta]` |
| IDM_DRMD (`interacting_species.cpp:358`) | `pvb[index_bg_rho_] * y[layout.idx_delta]` | `pvb[index_bg_rho_] * y[layout.idx_theta]` |
| Photons (`photons.cpp`) | `(layout.idx_delta >= 0) ? pvb[index_bg_rho_]*y[layout.idx_delta] : 0.` | `(layout.idx_theta >= 0) ? 4./3.*pvb[index_bg_rho_]*y[layout.idx_theta] : 0.` |
| UR (`ultra_relativistic.cpp`) | `(layout.idx_delta >= 0) ? Rho(pvb)*y[layout.idx_delta] : 0.` | `(layout.idx_theta >= 0) ? (Rho(pvb)+P(pvb))*y[layout.idx_theta] : 0.` |
| IDR (`interacting_species.cpp:206`) | `(layout.idx_delta >= 0) ? Rho(pvb)*y[layout.idx_delta] : 0.` | `(layout.idx_theta >= 0) ? (Rho(pvb)+P(pvb))*y[layout.idx_theta] : 0.` |
| IDR_DRMD (`interacting_species.cpp:437`) | `(layout.idx_delta >= 0) ? Rho(pvb)*y[layout.idx_delta] : 0.` | `(layout.idx_theta >= 0) ? (Rho(pvb)+P(pvb))*y[layout.idx_theta] : 0.` |
| Fluid (`fluid.cpp`) | `(layout.idx_delta >= 0) ? Rho(pvb)*y[layout.idx_delta] : 0.` | `(layout.idx_theta >= 0) ? (Rho(pvb)+P(pvb))*y[layout.idx_theta] : 0.` |
| Lambda (`lambda.h`) | `0.` | `0.` |

(Photon uses the historical `4./3.*ρ_g`; UR/IDR/IDR_DRMD use `(Rho+P)` to match the historical generic-loop `rho_plus_p` exactly.)

- [ ] **Step 4: Convert `dark_radiation_species.cpp` (division cancels)**

```cpp
double DarkRadiationSpecies::DeltaRho(const BaseSpecies::PerturbLayout& base, const perturb_vector*,
                                      const double* y, const double* pvecback, const perturb_workspace*) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0 || pvecback[index_bg_rho_] <= 0.) return 0.;
  double a = pvecback[bgm_->index_bg_a_], a2 = a*a;
  double rho_dr_over_f = (pba_->H0/a2)*(pba_->H0/a2);
  return rho_dr_over_f * y[layout.idx_F0];                       // was … / pvecback[index_bg_rho_]
}
double DarkRadiationSpecies::RhoPlusPTheta(const BaseSpecies::PerturbLayout& base, const perturb_vector*,
                                           const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0 || pvecback[index_bg_rho_] <= 0.) return 0.;
  double a = pvecback[bgm_->index_bg_a_], a2 = a*a;
  double rho_dr_over_f = (pba_->H0/a2)*(pba_->H0/a2);
  double k = ppw->scalar_ctx.k;
  return k * rho_dr_over_f * y[layout.idx_F0 + 1];               // (ρ+P)θ; the /ρ and 4/3·(3/4) cancel
}
```

- [ ] **Step 5: Convert `scalar_field.cpp` (division cancels)**

`DeltaRho`: identical to current `Delta` (`:374-397`) but return `delta_rho` instead of `delta_rho / rho` (drop the `rho==0` early-return and the final division; keep the Newtonian `psi` branch). `RhoPlusPTheta`: identical to current `Theta` (`:399-413`) but `return (1./3.) * k2/a2 * phi_prime * y[layout.idx_phi];` (drop `rho_plus_p` and its division/guard).

- [ ] **Step 6: Convert `ncdm_species.cpp` (division cancels, two branches)**

```cpp
double NCDMSpecies::DeltaRho(const BaseSpecies::PerturbLayout& base, const perturb_vector*,
                             const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.index_per_q.empty()) return 0.;
  if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on)
    return pvecback[index_bg_rho_] * y[layout.index_per_q[0]];
  const double a = ppw->scalar_ctx.a;
  double rho_delta_ncdm = 0.0;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q = q_[iq], epsilon = std::sqrt(q*q + std::pow(M_*a, 2));
    rho_delta_ncdm += q*q*epsilon * w_[iq] * y[layout.index_per_q[iq]];
  }
  return rho_delta_ncdm * factor_ * std::pow(pba_->a_today/a, 4);   // was … / pvecback[index_bg_rho_]
}
double NCDMSpecies::RhoPlusPTheta(const BaseSpecies::PerturbLayout& base, const perturb_vector*,
                                  const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.index_per_q.empty()) return 0.;
  if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on)
    return (pvecback[index_bg_rho_] + pvecback[index_bg_p_]) * y[layout.index_per_q[0] + 1];
  const double a = ppw->scalar_ctx.a;
  double rho_plus_p_theta_ncdm = 0.0;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q = q_[iq];
    rho_plus_p_theta_ncdm += q*q*q * w_[iq] * y[layout.index_per_q[iq] + 1];
  }
  return rho_plus_p_theta_ncdm * ppw->scalar_ctx.k * factor_ * std::pow(pba_->a_today/a, 4);  // was … / (ρ+p)
}
```

- [ ] **Step 7: Convert `dncdm_species.cpp` (scale the rescaled fraction)**

```cpp
double DNCDMSpecies::DeltaRho(const BaseSpecies::PerturbLayout& base, const perturb_vector*,
                              const double* /*y*/, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty()) return 0.;
  auto [d, t, s] = RescaledPerturbations(layout, ppw->scalar_ctx.a, ppw->scalar_ctx.k, ppw);
  return Rho(pvecback) * d;
}
double DNCDMSpecies::RhoPlusPTheta(const BaseSpecies::PerturbLayout& base, const perturb_vector*,
                                   const double* /*y*/, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty()) return 0.;
  auto [d, t, s] = RescaledPerturbations(layout, ppw->scalar_ctx.a, ppw->scalar_ctx.k, ppw);
  return (Rho(pvecback) + P(pvecback)) * t;
}
```

- [ ] **Step 8: Convert the four composites (child sums)**

`dcdm_dr_species.cpp` (`Delta` `:112-126`, `Theta` `:128-142`) → :

```cpp
double DCDM_DR_Species::DeltaRho(const BaseSpecies::PerturbLayout& base, const perturb_vector* pv,
                                 const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  return dcdm_->DeltaRho(my.dcdm, pv, y, pvecback, ppw) + dr_sp_->DeltaRho(my.dr, pv, y, pvecback, ppw);
}
double DCDM_DR_Species::RhoPlusPTheta(const BaseSpecies::PerturbLayout& base, const perturb_vector* pv,
                                      const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  return dcdm_->RhoPlusPTheta(my.dcdm, pv, y, pvecback, ppw) + dr_sp_->RhoPlusPTheta(my.dr, pv, y, pvecback, ppw);
}
```

Apply the identical pattern to the other three (rename + child sum), using each one's child accessors / sub-layout fields:

| Composite (file) | children → sub-layouts |
|---|---|
| `dncdm_dr_species.cpp` (`:116`,`:?`) | `dncdm_->…(my.dncdm,…) + dr_sp_->…(my.dr,…)` |
| `idm_dr_idr_species.cpp` (`:279`) | `idm_dr_->…(my.idm_dr,…) + idr_->…(my.idr,…)` |
| `idm_drmd_idr_drmd_species.cpp` (`:291`) | `idm_drmd_->…(my.idm_drmd,…) + idr_drmd_->…(my.idr_drmd,…)` |

Also update each composite's `MatterRhoDelta`/`MatterRhoPlusPTheta` overrides (e.g. `dcdm_dr_species.cpp:152/163`): replace `child_->Rho(pvb)*child_->Delta(...)` with `child_->DeltaRho(...)` and `(child_->Rho+child_->P)*child_->Theta(...)` with `child_->RhoPlusPTheta(...)`. (These overrides are deleted in Task 4, but must compile now.)

- [ ] **Step 9: Update the module call sites**

In `perturbations_module.cpp`:

- IC `delta_tot` loop (`:3971-3975`):
```cpp
      delta_rho_ic        += sp->DeltaRho(layout, ppw->pv.get(), ppw->pv->y.data(), ppw->pvecback.data(), ppw);
      rho_plus_p_theta_ic += sp->RhoPlusPTheta(layout, ppw->pv.get(), ppw->pv->y.data(), ppw->pvecback.data(), ppw);
```
  (drop the now-unused `rho`/`rho_plus_p` locals there).
- Generic pass-2 loop (`:4960-4962`):
```cpp
        ppw->delta_rho        += sp->DeltaRho(layout, ppw->pv.get(), y, ppw->pvecback.data(), ppw);
        ppw->rho_plus_p_theta += sp->RhoPlusPTheta(layout, ppw->pv.get(), y, ppw->pvecback.data(), ppw);
```
  (`rho_plus_p` is still needed for `rho_plus_p_tot`; keep it.)
- Photon/baryon dedicated block (`:4893-4905`): replace the four `PH.Delta/PH.Theta/BA.Delta/BA.Theta` calls with direct layout reads (byte-identical, transitional — block is removed in Task 3):
```cpp
      const double delta_g = (g_tse_lay.idx_delta >= 0) ? y[g_tse_lay.idx_delta] : 0.;
      const double theta_g = (g_tse_lay.idx_theta >= 0) ? y[g_tse_lay.idx_theta] : 0.;
      const double delta_b = y[b_tse_pre_lay.idx_delta];
      const double theta_b = y[b_tse_pre_lay.idx_theta];
```

- [ ] **Step 10: Update NCDM `FillSources` fractional consumers (`ncdm_species.cpp:285/297`)**

```cpp
          DeltaRho(layout, pv, y, pvecback, ppw) / pvecback[index_bg_rho_] +
              3. * ctx.a_prime_over_a * (1. + w) * ctx.theta_over_k2);
…
          RhoPlusPTheta(layout, pv, y, pvecback, ppw) / (pvecback[index_bg_rho_] + pvecback[index_bg_p_])
              + ctx.theta_shift);
```

- [ ] **Step 11: Build**

Run: `cmake --build build/cmake --parallel 2>&1 | tail -20`
Expected: links cleanly; no references to `->Delta(`/`->Theta(` remain (`grep -rn '\->Delta(\|\.Delta(\|->Theta(\|\.Theta(' source species` returns nothing).

- [ ] **Step 12: Golden diff — ΛCDM byte-identical, NCDM ≤0.1 %**

Run the CLI on a ΛCDM and an NCDM ini, diff outputs vs a baseline built from `master` (`git stash`/worktree). ΛCDM (`base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini` with `m_ncdm` removed, or `explanatory.ini` trimmed): expect **identical** `cl`/`pk`. NCDM ini: expect ≤0.1 % TT & mP(k).
```bash
build/cmake/class <lcdm.ini>; diff -q output/<lcdm>_cl.dat ref/<lcdm>_cl.dat
```

- [ ] **Step 13: Commit**

```bash
git add species/*.h species/*.cpp source/perturbations_module.cpp
git commit -m "perturb: DeltaRho/RhoPlusPTheta replace Delta/Theta in BaseSpecies"
```

---

## Task 2: Move photon-shear / baryon-pressure context into the species; cache `idr_nature_`

**Files:**
- Modify: `species/photons.h`, `species/photons.cpp`, `species/baryons.h`, `species/baryons.cpp`
- Modify: `source/perturbations_module.h` (`idr_nature_` member), `source/perturbations_module.cpp`
- Modify: `species/perturb_source_context.h` is unaffected; `scalar_ctx` struct holder (find: `grep -rn "double shear_g\|delta_p_b_over_rho_b" source species`)
- Test: CLI golden diff (byte-identical)

**Interfaces:**
- Consumes: `DeltaRho`/`RhoPlusPTheta` (Task 1).
- Produces: `scalar_ctx.shear_g` and `scalar_ctx.delta_p_b_over_rho_b` no longer exist; `PerturbationsModule::idr_nature_` (int) member.

- [ ] **Step 1: Give `PhotonsSpecies` stored `thm_`/`ppt_`**

In `photons.h`, add members + overrides (model on `ncdm_base_species.h:147-152`):
```cpp
  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override { thm_ = thm; }
  void SetPerturbs(const perturbs* ppt) override { ppt_ = ppt; }
 private:
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_ = nullptr;
```
Same for `baryons.h`.

- [ ] **Step 2: `PhotonsSpecies::RhoPlusPShear` self-computes the TCA/RSA shear**

Replace `photons.cpp:437-449` body:
```cpp
  const auto& layout = static_cast<const PerturbLayout&>(base);
  const double rho_g = pvecback[index_bg_rho_];
  if (layout.idx_shear >= 0) return 4./3.*rho_g*y[layout.idx_shear];          // no approximation
  if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) return 0.;             // RSA: shear neglected
  if (ppt_->gauge == possible_gauges::newtonian)                             // TCA, 1st order
    return 4./3.*rho_g*(16./45./ppw->pvecthermo[thm_->index_th_dkappa_]*y[layout.idx_theta]);
  return 0.;  // synchronous TCA: re-set later in perturb_einstein
```

- [ ] **Step 3: `BaryonsSpecies::DeltaP` self-computes δp_b**

Replace `baryons.cpp:101-107` body:
```cpp
  const auto& layout = static_cast<const PerturbLayout&>(base);
  const double rho_b = pvecback[index_bg_rho_];
  if (ppt_->has_perturbed_recombination && ppw->approx[ppw->index_ap_tca] == (int)tca_off)
    return rho_b * ppw->pvecthermo[thm_->index_th_wb_] *
           (y[layout.idx_delta] + y[pv->index_pt_perturbed_recombination_delta_temp]);
  return rho_b * ppw->pvecthermo[thm_->index_th_cb2_] * y[layout.idx_delta];
```
(Un-comment the `pv` and `y` parameters in the signature.)

- [ ] **Step 4: Add `idr_nature_` module member, resolve once**

In `perturbations_module.h` add `int idr_nature_ = idr_free_streaming;`. In the constructor (after `all_species_` is available; near where other species params are read), set:
```cpp
  if (all_species_.count("IDM_DR_IDR"))
    idr_nature_ = static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).idr().idr_nature();
```

- [ ] **Step 5: Point the current function's dedicated block at the species / member**

In `perturb_total_stress_energy` (still the old structure):
- Delete the shear-prep block (`:4806-4844`) and the δp_b-prep block (`:4846-4859`).
- Delete `ppw->scalar_ctx.shear_g = …` (`:4863`) and `ppw->scalar_ctx.delta_p_b_over_rho_b = …` (`:4864`).
- Replace the idr_nature block (`:4870-4878`) with `ppw->scalar_ctx.idr_nature = idr_nature_;`.
- In the dedicated arithmetic, source photon shear and baryon δp from the species methods:
```cpp
      ppw->rho_plus_p_shear = all_species_.photons().RhoPlusPShear(*ppw->pv->photon_layout,
                                  ppw->pv.get(), y, ppw->pvecback.data(), ppw);
      ppw->delta_p = 1./3.*rho_g*delta_g +
                     all_species_.baryons().DeltaP(*ppw->pv->baryon_layout, ppw->pv.get(), y,
                                                   ppw->pvecback.data(), ppw);
```

- [ ] **Step 6: Delete the `scalar_ctx` fields**

Remove `shear_g` and `delta_p_b_over_rho_b` from the `scalar_ctx` struct definition (found via `grep -rn "double shear_g" source species`). Build will flag any missed reader.

- [ ] **Step 7: Build + confirm no stragglers**

Run: `cmake --build build/cmake --parallel 2>&1 | tail -20`
Run: `grep -rn "scalar_ctx.shear_g\|scalar_ctx.delta_p_b_over_rho_b" source species` → empty.

- [ ] **Step 8: Golden diff (byte-identical)**

Run the photon-TCA-sensitive scenarios: a synchronous-gauge ΛCDM, a Newtonian-gauge run, and a `perturbed_recombination` run. Expect **byte-identical** `cl` (formulas only relocated).

- [ ] **Step 9: Commit**

```bash
git add species/photons.h species/photons.cpp species/baryons.h species/baryons.cpp \
        source/perturbations_module.h source/perturbations_module.cpp <scalar_ctx-header>
git commit -m "perturb: move photon shear / baryon delta_p into species; cache idr_nature"
```

---

## Task 3: Rewrite the scalar block as one loop + generic `(ρ−3P)` matter tally

**Files:**
- Modify: `source/perturbations.h` (`ActiveSpecies` struct)
- Modify: `source/perturbations_module.cpp` (`active_species` build sites + the scalar block)
- Test: CLI golden diff (ΛCDM ≤0.1 %; NCDM/IDM_DR expected physics) + regenerate reference

**Interfaces:**
- Consumes: `DeltaRho`/`RhoPlusPTheta`/`DeltaP`/`RhoPlusPShear`, the self-computing photon/baryon methods, `idr_nature_`.
- Produces: `ActiveSpecies{ species, layout, bool not_dark_energy, bool is_cold }`.

- [ ] **Step 1: Extend `ActiveSpecies` (`perturbations.h:259-263`)**

```cpp
  struct ActiveSpecies {
    BaseSpecies*                species;
    BaseSpecies::PerturbLayout* layout;
    bool                        not_dark_energy;  // energy_type() != DarkEnergy  (matter-tally inclusion)
    bool                        is_cold;          // IsColdMatterSpecies()        (delta_cb bucket)
  };
```

- [ ] **Step 2: Fill the bits at the three build sites (`perturbations_module.cpp:2968, 3006, 3036`)**

```cpp
        if (index_pt > before)
          ppv->active_species.push_back({entry.get(), ppv->species_layouts[i].get(),
                                         entry->energy_type() != BaseSpecies::EnergyType::DarkEnergy,
                                         entry->IsColdMatterSpecies()});
```

- [ ] **Step 3: Replace the scalar block (`perturbations_module.cpp:4795-5030`)**

```cpp
  if (_scalars_) {
    const double* pb = ppw->pvecback.data();
    const perturb_vector* pv = ppw->pv.get();

    ppw->scalar_ctx.k=k; ppw->scalar_ctx.k2=k2; ppw->scalar_ctx.a=a; ppw->scalar_ctx.a2=a2;
    ppw->scalar_ctx.gauge = static_cast<int>(ppt->gauge);
    ppw->scalar_ctx.idr_nature = idr_nature_;

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
        const double m_rho  = rho - 3.*p;
        const double m_drho = drho - 3.*dp;
        const double m_rppt = (rpp > 0.) ? m_rho*rppt/rpp : 0.;
        Tally& t = e.is_cold ? cold : warm;
        t.drho += m_drho; t.rho += m_rho; t.rppt += m_rppt; t.rpp += m_rho;
      }
    }

    if (has_source_delta_m_ && has_source_delta_cb_) ppw->delta_cb = cold.drho/cold.rho;
    if ((has_source_delta_m_||has_source_theta_m_)&&(has_source_delta_cb_||has_source_theta_cb_))
      ppw->theta_cb = cold.rppt/cold.rpp;

    if (auto* f = ppf_fluid()) {
      f->ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
      ppw->delta_rho        += ppw->delta_rho_fld;
      ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
      ppw->delta_p          += ppw->delta_p_fld;
      ppw->rho_plus_p_tot   += f->Rho(pb) + f->P(pb);
    }

    if (has_source_delta_m_)
      ppw->delta_m = (cold.drho+warm.drho)/(cold.rho+warm.rho);
    if (has_source_delta_m_||has_source_theta_m_)
      ppw->theta_m = (cold.rppt+warm.rppt)/(cold.rpp+warm.rpp);
  }
```

Leave the `_vectors_` and `_tensors_` blocks untouched.

- [ ] **Step 4: Build**

Run: `cmake --build build/cmake --parallel 2>&1 | tail -20` → clean.

- [ ] **Step 5: Verify the structural change is ULP-only on ΛCDM**

Run a ΛCDM scenario (no NCDM, no IDM_DR): expect ≤0.1 % TT & mP(k) vs the `master` baseline (cold matter `ρ−3P=ρ`, so only the lex-reorder ULP differs).

- [ ] **Step 6: Verify the intended physics changes**

- NCDM scenario: `delta_m`/P_m(k) shifts by ~`3w·f_ν` (small, now pressure-subtracted). Confirm sign/magnitude, not a blow-up.
- **IDM_DR scenario** (`grep -l idm_dr *.ini` or craft one with `Omega_idm_dr` dominant): confirm IDM_DR now **appears** in P_m(k) — compare against `master` (where it was absent) and sanity-check the large-scale P_m(k) tracks the IDM_DR density. This is the bug fix.
- PPF-fluid scenario: `rho_plus_p_tot` correct (verify Ω closure / H unaffected; `delta_rho_fld` path intact).

- [ ] **Step 7: Verify no non-active species drops a `ρ+P` term**

Confirm `rho_plus_p_tot` matches `master` (≤ULP) across ΛCDM + non-PPF-fluid + PPF-fluid: any non-active species with `ρ+P≠0` other than PPF would show here.

- [ ] **Step 8: Commit**

```bash
git add source/perturbations.h source/perturbations_module.cpp
git commit -m "perturb: single-loop stress-energy with generic (rho-3P) matter tally; fix IDM_DR P(k)"
```

---

## Task 4: Delete the dead `Matter*` family and `IsMatterSpecies`

**Files:**
- Modify: `species/base_species.h`, `species/base_species.cpp`
- Modify: `species/ncdm_base_species.h`, `species/idm_dr.h`, composite `*.cpp`/`*.h`
- Test: build + golden diff (byte-identical — pure dead-code removal)

**Interfaces:**
- Consumes: nothing (these symbols are unused after Task 3).
- Produces: `IsColdMatterSpecies()` base default becomes `energy_type_ == EnergyType::Matter`.

- [ ] **Step 1: Confirm the symbols are dead**

Run: `grep -rn "MatterRhoDelta\|MatterRhoPlusPTheta\|MatterRho\b\|MatterRhoPlusP\b\|IsMatterSpecies" source species`
Expected: matches only in *declarations/overrides* (no call sites in `source/`). If any call site remains, stop and route it through Task 3's tally instead.

- [ ] **Step 2: Delete from `base_species.h`**

Remove `MatterRho` (`:582-584`), `MatterRhoDelta` (`:592-595`), `MatterRhoPlusPTheta` (`:598-601`), `MatterRhoPlusP` (`:604-606`), and `IsMatterSpecies` (`:508-510`). Change `IsColdMatterSpecies` default body (`:614-616`) to:
```cpp
  virtual bool IsColdMatterSpecies() const { return energy_type_ == EnergyType::Matter; }
```

- [ ] **Step 3: Delete `base_species.cpp` out-of-line defs**

Remove `MatterRhoDelta` and `MatterRhoPlusPTheta` (the whole file body from Task 1 Step 2).

- [ ] **Step 4: Delete overrides in species**

- `ncdm_base_species.h`: delete `IsMatterSpecies` override (`:77-79`); **keep** `IsColdMatterSpecies()→false` (`:85-87`).
- `idm_dr.h`: delete the `IsMatterSpecies()→false` override (`:118-125`) — this is the asymmetry fix; IDM_DR now inherits `IsColdMatterSpecies()→true` (it is `EnergyType::Matter`).
- Composites: delete `MatterRho`/`MatterRhoPlusP` overrides (`composite_species.cpp:105/112`) and `MatterRhoDelta`/`MatterRhoPlusPTheta` overrides (`dcdm_dr_species.cpp:144-166`, and the equivalents in `dncdm_dr_species.cpp`, `idm_dr_idr_species.cpp`, `idm_drmd_idr_drmd_species.cpp`); delete their declarations in the matching `.h`. **Verify** each composite's `IsColdMatterSpecies` override still classifies correctly (DCDM_DR→cold via dcdm child; DNCDM_DR→warm via dncdm child) — if a composite relied on the base default that called `IsMatterSpecies`, give it an explicit child-scan override.

- [ ] **Step 5: Build + grep clean**

Run: `cmake --build build/cmake --parallel 2>&1 | tail -20` → clean.
Run: `grep -rn "IsMatterSpecies\|MatterRhoDelta\|MatterRhoPlusPTheta" source species` → empty.

- [ ] **Step 6: Golden diff (byte-identical to end of Task 3)**

Re-run the Task 3 scenarios; outputs must be **identical** to the post-Task-3 binary (this task removes only unreferenced code).

- [ ] **Step 7: Commit**

```bash
git add species/base_species.h species/base_species.cpp species/ncdm_base_species.h \
        species/idm_dr.h species/composite_species.cpp species/*dr_species.cpp species/*dr_species.h \
        species/*drmd_species.cpp species/*drmd_species.h
git commit -m "perturb: delete dead Matter* family + IsMatterSpecies; IsColdMatterSpecies on energy_type"
```

---

## Task 5: Full regression + reference regeneration + profile

**Files:** none (verification only)

- [ ] **Step 1: Rebuild `classy`**

Run: `pip install . 2>&1 | tail -5`

- [ ] **Step 2: Full `TEST_LEVEL=2` suite (pre-regeneration, expect the intended diffs)**

Run: `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest python/test_class.py -q 2>&1 | tail -30`
Expected: green except the documented intentional changes (warm-NCDM/DNCDM `delta_m`; IDM_DR scenarios) flagged by `COMPARE_OUTPUT_REF`. Inspect each flagged scenario and confirm it matches the spec's "Behavior changes".

- [ ] **Step 3: Regenerate `classyref` + goldens**

Rebuild `classyref` from this branch (the module-name swap per the CMake `CLASS_PYTHON_MODULE_NAME` flow) and regenerate any committed golden outputs. Then re-run Step 2 — expect fully green.

- [ ] **Step 4: Re-profile**

Run `class_profiled` on the Planck-2018 single-NCDM benchmark (single-threaded) and record the `perturbations` median vs the #325 baseline. Expect a reduction (no per-RHS string lookup, no second `Delta`/`Theta`, fused tally).

- [ ] **Step 5: Final commit (goldens/reference, if tracked)**

```bash
git add <regenerated golden files>
git commit -m "perturb: regenerate goldens/classyref for (rho-3P) matter tally + IDM_DR P(k) fix"
```

---

## Self-Review

**Spec coverage:** §1 interface → T1; §2 context relocation → T2; §3 cache → T3 Steps 1-2; §4 function → T3 Step 3; matter-tally `(ρ−3P)` → T3 Step 3; `Matter*`/`IsMatterSpecies` deletion → T4; IDM_DR fix → T3 Step 6 / T4 Step 4; IC `delta_tot` + NCDM FillSources → T1 Steps 9-10; behavior-change verification → T3 Steps 5-6, T5 Step 2; risks (thm_/ppt_ wiring, rho_plus_p_tot drop, composite cold/warm) → T2 Step 1, T3 Step 7, T4 Step 4. All covered.

**Placeholder scan:** the `<scalar_ctx-header>` / `<lcdm.ini>` / `<regenerated golden files>` tokens are deliberate "resolve-by-grep/ls" placeholders, each paired with the exact `grep`/`ls` to resolve them — not vague gaps.

**Type consistency:** `DeltaRho`/`RhoPlusPTheta` use the identical 5-arg signature of the removed `Delta`/`Theta` throughout; `ActiveSpecies` field names `not_dark_energy`/`is_cold` match between the struct (T3 S1), the fill sites (T3 S2), and the loop (T3 S3); `idr_nature_` consistent across T2 S4 and T3 S3.
