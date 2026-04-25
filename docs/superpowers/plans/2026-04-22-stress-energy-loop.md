# `perturb_total_stress_energy` Species-Dispatch Rewrite — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse `PerturbationsModule::perturb_total_stress_energy` (~400 lines) into a two-pass loop over `all_species_` (~20 lines) by moving species-local physics (fluid EoS, PPF, IDR TCA shear, DNCDM rescaling) into the species classes and adding `IsMatterSpecies` + four `Matter*` virtuals to `BaseSpecies`.

**Architecture:** Composites already expose correct aggregate `Rho`/`Delta`/`Theta`/`DeltaP`/`RhoPlusPShear` — the module just needs to use them. Each non-PPF sector participates in a symmetric loop. PPF is the single explicit special case (one allowed, depends on all other species, called after the main loop).

**Tech Stack:** C++17, `make class -j` for build, pytest regression suite in `python/test_class.py` for verification. No unit-test framework for C++ — verification is bit-identical reference comparison at the end.

**Spec:** [docs/superpowers/specs/2026-04-22-stress-energy-loop-design.md](../specs/2026-04-22-stress-energy-loop-design.md)

---

## Conventions

- After every code change: `make class -j 2>&1 | tail -10` and `./class explanatory.ini 2>&1 | tail -3`. Both must succeed.
- Fast smoke: `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py -k "test_scenario_00"`. Runs one scenario quickly.
- Reference gate (end only): `cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py`.
- Every new `.cpp`/`.h` must land in `Makefile`, `setup.py`, **and** `CLASS.xcodeproj/project.pbxproj`. No new files are expected in this plan; if a task adds one, add the build-system update as an explicit step.
- Commits use `git -c commit.gpgsign=false commit ...` if gpg signing is not set up. Follow the repo's existing commit-message style (short prefix like `background:`, `species:`, `perturbations:` matching the touched area).

---

## Task 1: Create branch

**Files:** none (git ops).

- [ ] **Step 1: Create issue** *(optional — if user wants GitHub issue tracking)*

```bash
cd /Users/au192734/Projects/class_claude
gh issue create \
  --title "Rewrite perturb_total_stress_energy as a species-dispatch loop" \
  --body "See docs/superpowers/specs/2026-04-22-stress-energy-loop-design.md for design, docs/superpowers/plans/2026-04-22-stress-energy-loop.md for plan."
```

Note the issue number. If you skip the issue, use `stress-energy-loop` as the branch name suffix.

- [ ] **Step 2: Create branch**

```bash
git checkout -b <issue-number>-stress-energy-loop master
```

- [ ] **Step 3: Baseline check**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build; `All parameters and all quantities computed successfully`. If not, something is broken on master before you start — stop and investigate.

---

## Task 2: Add `IsMatterSpecies` + four `Matter*` virtuals to `BaseSpecies`

**Files:**
- Modify: `species/base_species.h`
- Modify: `species/composite_species.h`
- Modify: `species/composite_species.cpp`
- Modify: `species/idm_dr.h`

- [ ] **Step 1: Add the five virtuals to `BaseSpecies`**

In `species/base_species.h`, find the section just before the `protected:` block (around line 301, right after `ApplyInitialConditions`). Add:

```cpp
  // ── Matter tally ──────────────────────────────────────────────────────────

  /**
   * True iff this sector participates in the matter tally (delta_m, theta_m,
   * and the delta_cb / theta_cb passes that cap out in perturb_total_stress_energy).
   * Default: based on EnergyType::Matter. Overrides: IDM_DR returns false
   * (preserves a pre-existing asymmetry in the matter tally — see follow-up issue).
   * Composites return true iff any child is matter.
   */
  virtual bool IsMatterSpecies() const {
    return energy_type_ == EnergyType::Matter;
  }

  /** Rho contribution to the matter tally. Default: Rho() if IsMatterSpecies, else 0. */
  virtual double MatterRho(const double* pvecback) const {
    return IsMatterSpecies() ? Rho(pvecback) : 0.;
  }

  /** Rho * Delta contribution to delta_rho_m. Default: Rho*Delta if IsMatterSpecies, else 0. */
  virtual double MatterRhoDelta(const perturb_vector* pv,
                                const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const {
    return IsMatterSpecies() ? Rho(pvecback) * Delta(pv, y, pvecback, ppw) : 0.;
  }

  /** (Rho+P) * Theta contribution to rho_plus_p_theta_m. */
  virtual double MatterRhoPlusPTheta(const perturb_vector* pv,
                                     const double* y,
                                     const double* pvecback,
                                     const perturb_workspace* ppw) const {
    return IsMatterSpecies() ? (Rho(pvecback) + P(pvecback)) * Theta(pv, y, pvecback, ppw) : 0.;
  }

  /** Rho+P contribution to rho_plus_p_m. */
  virtual double MatterRhoPlusP(const double* pvecback) const {
    return IsMatterSpecies() ? Rho(pvecback) + P(pvecback) : 0.;
  }
```

- [ ] **Step 2: Add composite overrides — declare in `species/composite_species.h`**

Find the `// ── Perturbations ─` block and after `RhoPlusPShear` (around line 86), add:

```cpp
  // ── Matter tally ────────────────────────────────────────────────────────
  bool IsMatterSpecies() const override;
  double MatterRho(const double* pvecback) const override;
  double MatterRhoDelta(const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const override;
  double MatterRhoPlusPTheta(const perturb_vector* pv,
                             const double* y,
                             const double* pvecback,
                             const perturb_workspace* ppw) const override;
  double MatterRhoPlusP(const double* pvecback) const override;
```

- [ ] **Step 3: Add composite implementations — `species/composite_species.cpp`**

Append at the end of the file (after the existing `RhoPlusPShear` implementation):

```cpp
bool CompositeSpecies::IsMatterSpecies() const {
  for (const auto& child : children_)
    if (child->IsMatterSpecies()) return true;
  return false;
}

double CompositeSpecies::MatterRho(const double* pvecback) const {
  double r = 0.;
  for (const auto& child : children_)
    r += child->MatterRho(pvecback);
  return r;
}

double CompositeSpecies::MatterRhoDelta(const perturb_vector* pv,
                                       const double* y,
                                       const double* pvecback,
                                       const perturb_workspace* ppw) const {
  double rd = 0.;
  for (const auto& child : children_)
    rd += child->MatterRhoDelta(pv, y, pvecback, ppw);
  return rd;
}

double CompositeSpecies::MatterRhoPlusPTheta(const perturb_vector* pv,
                                             const double* y,
                                             const double* pvecback,
                                             const perturb_workspace* ppw) const {
  double rpt = 0.;
  for (const auto& child : children_)
    rpt += child->MatterRhoPlusPTheta(pv, y, pvecback, ppw);
  return rpt;
}

double CompositeSpecies::MatterRhoPlusP(const double* pvecback) const {
  double rp = 0.;
  for (const auto& child : children_)
    rp += child->MatterRhoPlusP(pvecback);
  return rp;
}
```

- [ ] **Step 4: Add `IDM_DR` override in `species/idm_dr.h`**

Inside the `IDM_DRSpecies` class body (it already has `EnergyType::Matter`), add after the constructor:

```cpp
  /**
   * IDM_DR is excluded from the matter tally by current convention. This is
   * an asymmetry (IDM_DRMD, DCDM, NCDM are all included) that predates this
   * refactor — see follow-up issue.
   */
  bool IsMatterSpecies() const override { return false; }
```

- [ ] **Step 5: Build and verify**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build (possibly one or two unused-parameter warnings at worst), successful smoke test. No numerical change — the new methods are defined but not yet called anywhere.

- [ ] **Step 6: Commit**

```bash
git add species/base_species.h species/composite_species.h species/composite_species.cpp species/idm_dr.h
git -c commit.gpgsign=false commit -m "species: add IsMatterSpecies + Matter* virtuals to BaseSpecies

No behavior change. Adds the hooks that perturb_total_stress_energy will
dispatch through once rewritten as a loop. Composites override to sum over
matter children; IDM_DR opts out of the matter tally (pre-existing asymmetry,
tracked separately).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: Move `background_w_fld` body into `FluidSpecies::ComputeWFld`

`BackgroundModule::background_w_fld` becomes a thin delegator to preserve the 9+ call sites; the implementation moves into `FluidSpecies`.

**Files:**
- Modify: `species/fluid.h`
- Modify: `species/fluid.cpp`
- Modify: `source/background_module.h`
- Modify: `source/background_module.cpp`

- [ ] **Step 1: Add `bgm_` member and `SetBackgroundModule` to `FluidSpecies`**

In `species/fluid.h`, inside `FluidSpecies` class, add after the constructor (alongside the existing methods):

```cpp
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }

  /**
   * Fluid equation-of-state evaluation. Owns w_fld(a), dw_fld/da(a), and the
   * analytic integral used by background integration. Previously lived on
   * BackgroundModule; moved here because it's purely fluid physics.
   */
  int ComputeWFld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) const;
```

And in the `private:` section:

```cpp
  const BackgroundModule* bgm_ = nullptr;
```

Forward-declare `class BackgroundModule;` at the top of the file if not already done (check the existing fluid.h for the pattern — `dcdm_dr_species.h` is the reference).

- [ ] **Step 2: Move `background_w_fld` body into `FluidSpecies::ComputeWFld`**

Open `source/background_module.cpp`. The body of `BackgroundModule::background_w_fld` is at lines 458–560 (starting at `int BackgroundModule::background_w_fld(...)` and ending at `return _SUCCESS_;`).

Copy that body into `species/fluid.cpp` as `FluidSpecies::ComputeWFld`. The body uses `pba->...` fields and `all_species_.count(...)` calls. On `FluidSpecies`, `pba_` is already the `const background&` reference, and `bgm_->all_species_` provides the species collection.

In `species/fluid.cpp`, append:

```cpp
int FluidSpecies::ComputeWFld(double a,
                              double* w_fld,
                              double* dw_over_da_fld,
                              double* integral_fld) const {
  double Omega_ede          = 0.;
  double dOmega_ede_over_da = 0.;
  double a_eq               = 0.0;

  // w(a)
  switch (pba_.fluid_equation_of_state) {
    case CLP:
      *w_fld = pba_.w0_fld + pba_.wa_fld * (1. - a / pba_.a_today);
      break;
    case EDE: {
      Omega_ede = (pba_.Omega0_fld - pba_.Omega_EDE * (1. - pow(a, -3. * pba_.w0_fld))) /
                      (pba_.Omega0_fld + (1. - pba_.Omega0_fld) * pow(a, 3. * pba_.w0_fld)) +
                  pba_.Omega_EDE * (1. - pow(a, -3. * pba_.w0_fld));

      dOmega_ede_over_da =
          -pba_.Omega_EDE * 3. * pba_.w0_fld * pow(a, -3. * pba_.w0_fld - 1.) /
              (pba_.Omega0_fld + (1. - pba_.Omega0_fld) * pow(a, 3. * pba_.w0_fld)) -
          (pba_.Omega0_fld - pba_.Omega_EDE * (1. - pow(a, -3. * pba_.w0_fld))) *
              (1. - pba_.Omega0_fld) * 3. * pba_.w0_fld * pow(a, 3. * pba_.w0_fld - 1.) /
              pow(pba_.Omega0_fld + (1. - pba_.Omega0_fld) * pow(a, 3. * pba_.w0_fld), 2) +
          pba_.Omega_EDE * 3. * pba_.w0_fld * pow(a, -3. * pba_.w0_fld - 1.);

      double Omega_r =
          pba_.Omega0_g *
          (1. + 3.046 * 7. / 8. * pow(4. / 11., 4. / 3.));
      double Omega_m = pba_.Omega0_b;
      const auto& all_species = bgm_->all_species_;
      if (all_species.count("CDM"))               Omega_m += pba_.Omega0_cdm;
      if (all_species.count("IDM_DR_IDR"))        Omega_m += pba_.Omega0_idm_dr;
      if (all_species.count("IDM_DRMD_IDR_DRMD")) Omega_m += pba_.Omega0_idm_drmd;
      if (all_species.count("DCDM_DR"))
        class_stop(bgm_->error_message_,
                   "Early Dark Energy not compatible with decaying Dark Matter because "
                   "we omitted to code the calculation of a_eq in that case");
      a_eq = Omega_r / Omega_m;

      *w_fld = -dOmega_ede_over_da * a / Omega_ede / 3. / (1. - Omega_ede) +
               a_eq / 3. / (a + a_eq);
      break;
    }
  }

  // dw/da
  switch (pba_.fluid_equation_of_state) {
    case CLP:
      *dw_over_da_fld = -pba_.wa_fld / pba_.a_today;
      break;
    case EDE: {
      double d2Omega_ede_over_da2 = 0.;
      *dw_over_da_fld = -d2Omega_ede_over_da2 * a / 3. / (1. - Omega_ede) / Omega_ede -
                        dOmega_ede_over_da / 3. / (1. - Omega_ede) / Omega_ede +
                        dOmega_ede_over_da * dOmega_ede_over_da * a / 3. / (1. - Omega_ede) /
                            (1. - Omega_ede) / Omega_ede +
                        a_eq / 3. / (a + a_eq) / (a + a_eq);
      break;
    }
  }

  // Analytic integral: \int_a^{a_today} 3 (1 + w_fld(a')) / a' da'.
  // Used only in background initial conditions. If a future w(a) has no
  // closed form, set *integral_fld = 0 and integrate numerically at IC time.
  switch (pba_.fluid_equation_of_state) {
    case CLP:
      *integral_fld = 3. * ((1. + pba_.w0_fld + pba_.wa_fld) * log(pba_.a_today / a) +
                            pba_.wa_fld * (a / pba_.a_today - 1.));
      break;
    case EDE:
      class_stop(bgm_->error_message_,
                 "EDE implementation not finished: to finish it, read the comments in background.c "
                 "just before this line\n");
      break;
  }

  return _SUCCESS_;
}
```

- [ ] **Step 3: Reduce `BackgroundModule::background_w_fld` to a delegator**

Replace the body of `BackgroundModule::background_w_fld` in `source/background_module.cpp` (currently lines 458–560) with:

```cpp
int BackgroundModule::background_w_fld(double a,
                                       double* w_fld,
                                       double* dw_over_da_fld,
                                       double* integral_fld) const {
  return static_cast<FluidSpecies&>(*all_species_.at("Fluid"))
      .ComputeWFld(a, w_fld, dw_over_da_fld, integral_fld);
}
```

Add `#include "../species/fluid.h"` at the top of `source/background_module.cpp` if not already present.

- [ ] **Step 4: Wire `SetBackgroundModule` on the fluid instance**

Find the place where other species receive `SetBackgroundModule` calls — grep:

```bash
grep -n "SetBackgroundModule" source/*.cpp source/*.h species/*.cpp species/*.h
```

Expected: `BackgroundModule` already calls `SetBackgroundModule` on each species during construction (often in `BackgroundModule::BackgroundModule` or a helper). Ensure `FluidSpecies` is included in that sweep. If the existing sweep is a loop over `all_species_`, no change is needed beyond Step 1 (which added the override).

- [ ] **Step 5: Build and verify**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py -k "test_scenario_00" 2>&1 | tail -10
cd ..
```

Expected: clean build, smoke-test passes, test_scenario_00 passes. Numerics are identical (same code, just moved).

- [ ] **Step 6: Commit**

```bash
git add species/fluid.h species/fluid.cpp source/background_module.cpp source/background_module.h
git -c commit.gpgsign=false commit -m "fluid: move background_w_fld body into FluidSpecies::ComputeWFld

BackgroundModule::background_w_fld becomes a thin delegator to the species.
All call sites untouched.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 4: Move non-PPF fluid math into `FluidSpecies::Delta`/`Theta`/`DeltaP`

Make `FluidSpecies` return correct perturbation contributions for the non-PPF case so the main loop can dispatch through it.

**Files:**
- Modify: `species/fluid.h`
- Modify: `species/fluid.cpp`
- Modify: `source/perturbations_module.cpp` (only the fluid block inside `perturb_total_stress_energy`)

- [ ] **Step 1: Inspect current `FluidSpecies::Delta/Theta/DeltaP`**

Open `species/fluid.cpp` around line 223. Current:

```cpp
double FluidSpecies::Delta(const perturb_vector* pv, ...) const {
  return (pv->index_pt_delta_fld >= 0) ? y[pv->index_pt_delta_fld] : 0.;
}
double FluidSpecies::Theta(const perturb_vector* pv, ...) const {
  return (pv->index_pt_theta_fld >= 0) ? y[pv->index_pt_theta_fld] : 0.;
}
double FluidSpecies::DeltaP(...) const {
  return 0.;  // wrong for non-PPF
}
```

`Delta` and `Theta` return the raw `y[]` values. For non-PPF, that gives the correct `delta_fld` (contribution to `delta_rho` is `rho * delta_fld`) and `theta_fld` (contribution to `rho_plus_p_theta` needs `(rho + p) * theta_fld = (1 + w) * rho * theta_fld`). But `Theta` currently returns raw `y[theta_fld]` without multiplying by `(1+w)` elsewhere — the main loop currently does `rho_plus_p_theta_fld = (1 + w_fld) * Rho * y[theta_fld]`. Since the main loop accumulates `rho_plus_p * Theta`, if `Theta` returns raw `y[theta_fld]` then `(rho+p)*Theta = (1+w)*rho*y[theta_fld]` — consistent. OK, `Theta` can stay as-is.

The only method that currently lies is `DeltaP`: it must return the non-PPF δp.

- [ ] **Step 2: Rewrite `FluidSpecies::DeltaP` for the non-PPF path**

In `species/fluid.cpp`, replace:

```cpp
double FluidSpecies::DeltaP(const perturb_vector* /*pv*/,
                            const double* /*y*/,
                            const double* /*pvecback*/,
                            const perturb_workspace* /*ppw*/) const {
  return 0.;
}
```

with:

```cpp
double FluidSpecies::DeltaP(const perturb_vector* pv,
                            const double* y,
                            const double* pvecback,
                            const perturb_workspace* ppw) const {
  // PPF uses a dedicated path (ComputePpf) and this method is not called
  // for PPF — the module skips FluidSpecies in the main loop when use_ppf.
  if (pv->index_pt_delta_fld < 0 || pv->index_pt_theta_fld < 0) return 0.;

  const double k2             = ppw->scalar_ctx.k2;
  const double a              = ppw->scalar_ctx.a;
  const double a_prime_over_a = pvecback[bgm_->index_bg_H_] * a;  // aH

  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld);
  const double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

  const double rho               = Rho(pvecback);
  const double delta_rho_fld     = rho * y[pv->index_pt_delta_fld];
  const double rho_plus_p_theta_fld =
      (1. + w_fld) * rho * y[pv->index_pt_theta_fld];
  const double ca2_fld = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;

  return pba_.cs2_fld * delta_rho_fld +
         (pba_.cs2_fld - ca2_fld) * (3. * a_prime_over_a * rho_plus_p_theta_fld / k2);
}
```

- [ ] **Step 3: Verify `scalar_ctx.a`, `scalar_ctx.k2`, `bgm_->index_bg_H_` are available**

```bash
grep -n "scalar_ctx.a\b\|scalar_ctx.k2\b" source/*.cpp species/*.cpp species/*.h | head
grep -n "index_bg_H_" source/*.h source/*.cpp | head
```

Expected: `scalar_ctx.a` and `scalar_ctx.k2` are already populated in `perturb_total_stress_energy` at lines 6186–6187; `BackgroundModule::index_bg_H_` is a public member (grep confirms). If `scalar_ctx.a` is not populated before the fluid block today, add it in the preamble.

- [ ] **Step 4: Do NOT yet change the module-side fluid block**

Leave the fluid block at `source/perturbations_module.cpp:6430–6553` unchanged for this task. `FluidSpecies::DeltaP` now returns the correct value, but the module still writes to `ppw->delta_rho_fld` etc. directly. The old path and new path coexist without conflict because nothing in the module currently calls `FluidSpecies::DeltaP` for non-PPF — that happens in Task 8.

- [ ] **Step 5: Build and verify**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py -k "test_scenario_00" 2>&1 | tail -10
cd ..
```

Expected: build clean, smoke test passes, numerics unchanged (the new `DeltaP` implementation is dead code for now).

- [ ] **Step 6: Commit**

```bash
git add species/fluid.cpp species/fluid.h
git -c commit.gpgsign=false commit -m "fluid: implement non-PPF DeltaP correctly inside FluidSpecies

Currently unused by the module. Task 8 will dispatch the main loop through
FluidSpecies::Delta/Theta/DeltaP for the non-PPF case.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 5: Move PPF body into `FluidSpecies::ComputePpf`

**Files:**
- Modify: `species/fluid.h`
- Modify: `species/fluid.cpp`
- Modify: `source/perturbations_module.cpp` (fluid block only)

- [ ] **Step 1: Declare `ComputePpf` in `species/fluid.h`**

Add inside the `FluidSpecies` class, after the existing perturbation methods:

```cpp
  /**
   * Compute PPF fluid contribution to ppw (writes ppw->delta_rho_fld,
   * rho_plus_p_theta_fld, delta_p_fld, S_fld, Gamma_prime_fld, and adds
   * to ppw->delta_rho / rho_plus_p_theta / delta_p / rho_plus_p_tot).
   * Only called when pba->use_ppf == _TRUE_. Moved from PerturbationsModule
   * — PPF is fluid-specific physics that depends on the rest of the universe.
   */
  void ComputePpf(double k, double a, double a_prime_over_a,
                  const double* y, perturb_workspace* ppw) const;
```

- [ ] **Step 2: Copy PPF body into `FluidSpecies::ComputePpf`**

In `source/perturbations_module.cpp`, the PPF body is at lines 6448–6546 (between `else {` and the closing `}` just before `ppw->delta_rho += ppw->delta_rho_fld;`).

In `species/fluid.cpp`, append:

```cpp
void FluidSpecies::ComputePpf(double k,
                              double a,
                              double a_prime_over_a,
                              const double* y,
                              perturb_workspace* ppw) const {
  const double a2 = a * a;
  const double k2 = k * k;

  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld);
  const double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

  const double s2sq               = ppw->s_l[2] * ppw->s_l[2];
  const double c_gamma_k_H_square = pow(pba_.c_gamma_over_c_fld * k / a_prime_over_a, 2) *
                                    pba_.cs2_fld;
  // <<< paste lines 6455-6546 verbatim here, replacing:
  //     - pba->               with pba_.
  //     - all_species_.at("Fluid")->Rho(ppw->pvecback)   with Rho(ppw->pvecback)
  //     - background_module_->index_bg_H_              with bgm_->index_bg_H_
  //     - background_module_->index_bg_H_prime_        with bgm_->index_bg_H_prime_
  //     - background_module_->index_bg_rho_tot_        with bgm_->index_bg_rho_tot_
  //     - background_module_->index_bg_p_tot_          with bgm_->index_bg_p_tot_
  //     - background_module_->index_bg_p_tot_prime_    with bgm_->index_bg_p_tot_prime_
  //     - ppr->                                        with p_mod_->GetPrecision()->  (see next step)
  //     Do NOT include the trailing `ppw->delta_rho += ppw->delta_rho_fld;` block
  //     at lines 6548-6552 — that stays in the module for now, called
  //     after this function.

  // After paste:
  // (nothing else — the function's outputs are the ppw side-channel writes.)
}
```

**Precision access:** `ppr` in the module is a `BaseModule` protected member. In `FluidSpecies`, the only way to reach precision is via a stored pointer. Add `SetPrecision(const precision* ppr) { ppr_ = ppr; }` to `FluidSpecies` (mirror the `SetBackgroundModule` pattern) and a `const precision* ppr_ = nullptr;` member. Wire it alongside `SetBackgroundModule` at species-construction time. Replace `ppr->` with `ppr_->` in the pasted body.

Alternatively (simpler): pass `ppr` as an argument to `ComputePpf` from the module.

Choose the simpler: pass `ppr` as an argument. Update the declaration:

```cpp
  void ComputePpf(double k, double a, double a_prime_over_a,
                  const precision* ppr,
                  const double* y, perturb_workspace* ppw) const;
```

And use `ppr->` inside verbatim.

- [ ] **Step 3: Replace the module-side PPF body with a call**

In `source/perturbations_module.cpp`, edit the fluid block. Old body (lines 6430–6553, roughly):

```cpp
if (all_species_.count("Fluid")) {
  double w_fld, dw_over_da_fld, integral_fld;
  class_call(background_module_->background_w_fld(...), ...);
  double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

  if (pba->use_ppf == _FALSE_) {
    // non-PPF body (lines 6438-6447)
  }
  else {
    // PPF body (lines 6448-6546)
  }

  ppw->delta_rho        += ppw->delta_rho_fld;
  ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
  ppw->delta_p          += ppw->delta_p_fld;
  ppw->rho_plus_p_tot   += (1. + w_fld) * all_species_.at("Fluid")->Rho(ppw->pvecback);
}
```

Replace with:

```cpp
if (all_species_.count("Fluid")) {
  auto& fluid = static_cast<FluidSpecies&>(*all_species_.at("Fluid"));
  double w_fld, dw_over_da_fld, integral_fld;
  class_call(fluid.ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld),
             background_module_->error_message_, error_message_);

  if (pba->use_ppf == _FALSE_) {
    ppw->delta_rho_fld        = fluid.Rho(ppw->pvecback) * y[ppw->pv->index_pt_delta_fld];
    ppw->rho_plus_p_theta_fld = (1. + w_fld) * fluid.Rho(ppw->pvecback) *
                                y[ppw->pv->index_pt_theta_fld];
    const double w_prime_fld  = dw_over_da_fld * a_prime_over_a * a;
    const double ca2_fld      = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;
    ppw->delta_p_fld = pba->cs2_fld * ppw->delta_rho_fld +
                       (pba->cs2_fld - ca2_fld) *
                           (3 * a_prime_over_a * ppw->rho_plus_p_theta_fld / k / k);
  }
  else {
    fluid.ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
  }

  ppw->delta_rho        += ppw->delta_rho_fld;
  ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
  ppw->delta_p          += ppw->delta_p_fld;
  ppw->rho_plus_p_tot   += (1. + w_fld) * fluid.Rho(ppw->pvecback);
}
```

The non-PPF inline body stays for now (Task 8 removes it). Only PPF has moved.

- [ ] **Step 4: Build and verify**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py 2>&1 | tail -15
cd ..
```

Expected: clean build, full fast suite passes. Bit-identical.

- [ ] **Step 5: Commit**

```bash
git add species/fluid.h species/fluid.cpp source/perturbations_module.cpp
git -c commit.gpgsign=false commit -m "fluid: move PPF body into FluidSpecies::ComputePpf

Module calls it from the fluid block. Non-PPF inline body stays until
the main-loop rewrite in a later commit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 6: DNCDM rescaling to species; delete `ppw->delta_ncdm[n]` / `theta_ncdm[n]` / `shear_ncdm[n]`

Four sub-commits. Each one builds and smoke-tests.

**Files:**
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`
- Modify: `species/ncdm_species.cpp`
- Modify: `source/perturbations_module.cpp`
- Modify: `source/perturbations_module.h`
- Modify: `source/perturbations.h` (if `ppw->delta_ncdm` is declared there)

### Task 6a: Move `RescaledNCDMPerturbations` into `DNCDMSpecies` as a helper

- [ ] **Step 1: Declare `DNCDMSpecies::RescaledPerturbations`**

In `species/dncdm_species.h`, add to the class body (private section):

```cpp
 private:
  /**
   * Returns rescaled (delta, theta, shear) for this decaying NCDM flavor.
   * Used by this species' Delta/Theta/RhoPlusPShear overrides when lnf
   * is near the precision floor.
   */
  std::tuple<double, double, double> RescaledPerturbations(double a, double k,
                                                           const perturb_workspace* ppw) const;
```

- [ ] **Step 2: Copy the body from `PerturbationsModule::RescaledNCDMPerturbations`**

Source at `source/perturbations_module.cpp:8233–8282`. Paste into `species/dncdm_species.cpp`:

```cpp
std::tuple<double, double, double> DNCDMSpecies::RescaledPerturbations(
    double a, double k, const perturb_workspace* ppw) const {
  double rho_scaled              = 0.;
  double rho_plus_p_scaled       = 0.;
  double rho_delta_scaled        = 0.;
  double rho_plus_p_theta_scaled = 0.;
  double rho_plus_p_shear_scaled = 0.;

  const double* lnf_array = ppw->pvecback + bg_lnf_index();
  const double lnN        = GetRescalingFactor(lnf_array);

  for (int index_q = 0; index_q < q_size(); index_q++) {
    const int index_pt = ppw->pv->index_ncdm_[ncdm_id()][index_q];
    const double dq    = dq()[index_q];
    const double lnf   = lnf_array[index_q];
    const double q     = q()[index_q];
    const double q2    = q * q;
    const double epsilon = sqrt(q2 + a * a * M() * M());

    rho_scaled              += dq * pow(q, 2) * epsilon * exp(lnN + lnf);
    rho_plus_p_scaled       += dq * pow(q, 2) * (epsilon + q2 / 3. / epsilon) *
                               exp(lnN + lnf);
    rho_delta_scaled        += dq * pow(q, 2) * epsilon * exp(lnN + lnf) *
                               ppw->pv->y[index_pt];
    rho_plus_p_theta_scaled += dq * pow(q, 3) * exp(lnN + lnf) *
                               ppw->pv->y[index_pt + 1];
    rho_plus_p_shear_scaled += dq * pow(q, 4) / epsilon * exp(lnN + lnf) *
                               ppw->pv->y[index_pt + 2];
  }
  rho_plus_p_theta_scaled *= k;
  rho_plus_p_shear_scaled *= 2. / 3.;

  const double delta = rho_delta_scaled / rho_scaled;
  const double theta = rho_plus_p_theta_scaled / rho_plus_p_scaled;
  const double shear = rho_plus_p_shear_scaled / rho_plus_p_scaled;

  return {delta, theta, shear};
}
```

Add `#include <tuple>` at the top of `dncdm_species.cpp` if not present.

- [ ] **Step 3: Reduce `PerturbationsModule::RescaledNCDMPerturbations` to a delegator**

Edit `source/perturbations_module.cpp:8233–8282` to:

```cpp
std::tuple<double, double, double> PerturbationsModule::RescaledNCDMPerturbations(
    int n_ncdm, double a, double k, perturb_workspace* ppw) {
  for (auto& sp : all_species_) {
    if (auto* composite = dynamic_cast<DNCDM_DR_Species*>(sp.get())) {
      if (composite->dncdm().ncdm_id() == n_ncdm) {
        return composite->dncdm().RescaledPerturbations(a, k, ppw);
      }
    }
  }
  throw std::runtime_error("RescaledNCDMPerturbations: invalid ncdm_id");
}
```

- [ ] **Step 4: Build and smoke test**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
```

- [ ] **Step 5: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp source/perturbations_module.cpp
git -c commit.gpgsign=false commit -m "dncdm: move RescaledNCDMPerturbations into DNCDMSpecies

PerturbationsModule method becomes a thin lookup that delegates. No
behavior change yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 6b: Override `Delta`/`Theta`/`RhoPlusPShear` on `DNCDMSpecies`

- [ ] **Step 1: Add overrides to `DNCDMSpecies`**

In `species/dncdm_species.h`, add to public section:

```cpp
  double Delta(const perturb_vector* pv, const double* y,
               const double* pvecback, const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv, const double* y,
               const double* pvecback, const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv, const double* y,
                       const double* pvecback, const perturb_workspace* ppw) const override;

  double MatterRhoDelta(const perturb_vector* pv, const double* y,
                        const double* pvecback, const perturb_workspace* ppw) const override;
  double MatterRhoPlusPTheta(const perturb_vector* pv, const double* y,
                             const double* pvecback, const perturb_workspace* ppw) const override;
```

- [ ] **Step 2: Implement in `species/dncdm_species.cpp`**

```cpp
double DNCDMSpecies::Delta(const perturb_vector* /*pv*/, const double* /*y*/,
                           const double* pvecback, const perturb_workspace* ppw) const {
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(a, k, ppw);
  return d;
}
double DNCDMSpecies::Theta(const perturb_vector* /*pv*/, const double* /*y*/,
                           const double* pvecback, const perturb_workspace* ppw) const {
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(a, k, ppw);
  return t;
}
double DNCDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/,
                                   const double* pvecback, const perturb_workspace* ppw) const {
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(a, k, ppw);
  return (Rho(pvecback) + P(pvecback)) * s;
}

double DNCDMSpecies::MatterRhoDelta(const perturb_vector* pv, const double* y,
                                    const double* pvecback,
                                    const perturb_workspace* ppw) const {
  return IsMatterSpecies() ? Rho(pvecback) * Delta(pv, y, pvecback, ppw) : 0.;
}
double DNCDMSpecies::MatterRhoPlusPTheta(const perturb_vector* pv, const double* y,
                                         const double* pvecback,
                                         const perturb_workspace* ppw) const {
  return IsMatterSpecies() ? (Rho(pvecback) + P(pvecback)) * Theta(pv, y, pvecback, ppw) : 0.;
}
```

*Note:* `Delta`/`Theta`/`RhoPlusPShear` each recompute the rescaling. If profiling shows this is a hot spot, a per-ppw cache becomes the follow-up. Do not optimize now.

- [ ] **Step 3: Build and smoke test**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py 2>&1 | tail -15
cd ..
```

Expected: clean build, tests pass. Numerics unchanged because the module still uses the old `RescaledNCDMPerturbations` path.

- [ ] **Step 4: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp
git -c commit.gpgsign=false commit -m "dncdm: override Delta/Theta/RhoPlusPShear to return rescaled values

Used by Task 8's new main loop. Not called yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 6c: Replace all consumers of `ppw->delta_ncdm[n]` / `theta_ncdm[n]` / `shear_ncdm[n]`

- [ ] **Step 1: Find all consumers**

```bash
grep -n "ppw->delta_ncdm\|ppw->theta_ncdm\|ppw->shear_ncdm" source/*.cpp source/*.h species/*.cpp species/*.h
```

Expected sites (as of master):

| File | Line | Op | Context |
|------|------|----|---------|
| `source/perturbations_module.cpp` | ~6382–6390 | write | stress-energy function |
| `source/perturbations_module.cpp` | ~6401, 6405 | read | stress-energy matter tally |
| `species/ncdm_species.cpp` | ~220, 231 | read | NCDMSpecies::FillSources |

If the grep finds additional sites, add steps to update them. Do not proceed until the full list is enumerated.

- [ ] **Step 2: Update `species/ncdm_species.cpp:220` (`FillSources` — delta_ncdm[n] read)**

Find the block (around line 215–225) that uses `ppw->delta_ncdm[n]`. Replace `ppw->delta_ncdm[n]` with `Delta(pv, y, ppw->pvecback, ppw)` and `ppw->theta_ncdm[n]` with `Theta(pv, y, ppw->pvecback, ppw)`. For NCDMSpecies these return the same unrescaled values; for DNCDMSpecies they return rescaled — matching prior behavior.

Example (replace old):
```cpp
set_source(...,  ppw->delta_ncdm[n] + 3. * ctx.a_prime_over_a * (1. + w) * ...);
```
with:
```cpp
set_source(..., Delta(pv, y, ppw->pvecback, ppw) +
               3. * ctx.a_prime_over_a * (1. + w) * ...);
```

(Follow the exact signature of the surrounding code — `pv` and `y` are in scope via the context object or passed in.)

Verify `pv` and `y` availability in the surrounding method signature. If not, pass them or reach through ctx.

- [ ] **Step 3: Build and smoke test**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py 2>&1 | tail -15
cd ..
```

Expected: bit-identical.

- [ ] **Step 4: Commit**

```bash
git add species/ncdm_species.cpp
git -c commit.gpgsign=false commit -m "ncdm: replace ppw->delta_ncdm[n]/theta_ncdm[n] reads with Delta()/Theta()

Preparation for removing the ppw arrays.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 6d: Delete `ppw->delta_ncdm` / `theta_ncdm` / `shear_ncdm`

(The stress-energy function is the only remaining consumer; Task 8 will remove it. For a bisect-friendly history, defer the field deletion until after Task 8 lands. Skip 6d here and add it to Task 8.)

---

## Task 7: Move IDR TCA shear into `IDRSpecies`; delete `ppw->tca_shear_idm_dr`

**Files:**
- Modify: `species/base_species.h`
- Modify: `species/idr.h`
- Modify: `species/interacting_species.cpp`
- Modify: `species/idm_dr_idr_species.cpp`
- Modify: `source/perturbations_module.cpp`
- Modify: `source/perturbations.h`

- [ ] **Step 1: Add `SetThermodynamicsModule` and `SetPerturbs` hooks to `BaseSpecies`**

Two mutator hooks are needed here, mirroring the existing `SetBackgroundModule`. `alpha_idm_dr[0]` lives on the `perturbs` struct (`source/perturbations.h:199`), so IDR needs `ppt` too.

In `species/base_species.h`, add near the existing `SetBackgroundModule` (around line 55):

```cpp
// At top of file (alongside existing forward declarations):
class ThermodynamicsModule;
struct perturbs;  // (may already be forward-declared)

// Inside BaseSpecies:
virtual void SetThermodynamicsModule(const ThermodynamicsModule* /*thm*/) {}
virtual void SetPerturbs(const perturbs* /*ppt*/) {}
```

- [ ] **Step 2: Add `thm_` / `ppt_` members and `TcaShearIdr` helper to `IDRSpecies`**

In `species/idr.h`, inside the class body:

```cpp
  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override { thm_ = thm; }
  void SetPerturbs(const perturbs* ppt) override { ppt_ = ppt; }

  /**
   * IDR shear under TCA with IDM_DR. Returns 0 outside the TCA regime.
   * Reads pvecthermo[index_th_dmu_idm_dr_] for the scattering rate and
   * ppt_->alpha_idm_dr[0] for the angular collisional factor.
   */
  double TcaShearIdr(const perturb_vector* pv, const double* y,
                     const perturb_workspace* ppw) const;

 private:
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_ = nullptr;
```

Add `#include "thermodynamics_module.h"` in `species/idr.h` (or forward-declare to match existing patterns).

- [ ] **Step 3: Implement `TcaShearIdr` and update `RhoPlusPShear`**

In `species/interacting_species.cpp`, add:

```cpp
double IDRSpecies::TcaShearIdr(const perturb_vector* pv, const double* y,
                               const perturb_workspace* ppw) const {
  if (ppw->scalar_ctx.idr_nature != idr_free_streaming) return 0.;
  if (pv->index_pt_shear_idr >= 0) return 0.;            // shear integrated → not TCA
  if (ppw->approx[ppw->index_ap_rsa_idr] != (int) rsa_idr_off) return 0.;
  if (ppw->approx[ppw->index_ap_tca_idm_dr] != (int) tca_idm_dr_on) return 0.;
  if (!pba_.has_idm_dr) return 0.;
  if (ppw->scalar_ctx.gauge != newtonian) return 0.;     // synchronous: set in perturb_einstein

  return 0.5 * (8. / 15. /
                ppw->pvecthermo[thm_->index_th_dmu_idm_dr_] /
                ppt_->alpha_idm_dr[0] *
                y[pv->index_pt_theta_idr]);
}
```

Modify `IDRSpecies::RhoPlusPShear` at `species/interacting_species.cpp:129-140`:

```cpp
double IDRSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  if (ppw->scalar_ctx.idr_nature != idr_free_streaming) return 0.;
  const double shear_idr = (pv->index_pt_shear_idr >= 0)
                               ? y[pv->index_pt_shear_idr]
                               : TcaShearIdr(pv, y, ppw);
  return 4. / 3. * pvecback[index_bg_rho_] * shear_idr;
}
```

- [ ] **Step 4: Update `IDM_DR_IDR_Species::PrintVariables` (species/idm_dr_idr_species.cpp:152)**

Replace:
```cpp
shear_idr = ppw->tca_shear_idm_dr;
```
with:
```cpp
shear_idr = idr_->TcaShearIdr(pv, y, ppw);
```

Make `TcaShearIdr` public or add a public wrapper on `IDRSpecies`. Friend-class or public — whichever is cleaner.

- [ ] **Step 5: Delete the preamble block in the module**

In `source/perturbations_module.cpp`, delete lines 6157–6178 (the IDR under-TCA shear preamble block). The new code in `IDRSpecies::TcaShearIdr` replaces its effect.

- [ ] **Step 6: Check `tca_shear_idm_dr` still needed**

```bash
grep -n "tca_shear_idm_dr" source/*.cpp source/*.h species/*.cpp species/*.h
```

Remaining site of concern: `source/perturbations_module.cpp:4175` writes `ppv->y[ppv->index_pt_shear_idr] = ppw->tca_shear_idm_dr;` during `perturb_vector_init`. This reads the TCA shear at vector initialization. Replace with a `TcaShearIdr(...)` call on the IDR sub-species (reached via the IDM_DR_IDR composite). Confirm the call-site has `pv`, `y`, `ppw` in scope.

Once this site is updated, the `ppw->tca_shear_idm_dr` field is dead.

- [ ] **Step 7: Delete the field**

In `source/perturbations.h:391`, delete:
```cpp
double tca_shear_idm_dr;
```

- [ ] **Step 8: Build and test**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py 2>&1 | tail -15
cd ..
```

Expected: bit-identical. If tests fail at an IDM_DR_IDR scenario, the TCA formula has not been faithfully preserved — diff the removed and new code carefully.

- [ ] **Step 9: Commit**

```bash
git add species/base_species.h species/idr.h species/interacting_species.cpp species/idm_dr_idr_species.cpp source/perturbations_module.cpp source/perturbations.h
git -c commit.gpgsign=false commit -m "idr: move TCA shear computation into IDRSpecies

Delete ppw->tca_shear_idm_dr side channel. Add BaseSpecies::SetThermodynamicsModule
hook so IDRSpecies can read pvecthermo[index_th_dmu_idm_dr_].

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 8: Rewrite `perturb_total_stress_energy` as two-pass loop + PPF call

This is the big commit. The old per-species blocks (6198–6553) and the NCDM dedicated loop (6360–6409) all get deleted and replaced with two loops.

**Files:**
- Modify: `source/perturbations_module.cpp` (function body)
- Modify: `source/perturbations.h` (remove `ppw->delta_ncdm`, `theta_ncdm`, `shear_ncdm` arrays)

- [ ] **Step 1: Verify the ψ-ordering assumption**

```bash
grep -nE "scalar_ctx\.psi\b|->psi\b" source/*.cpp species/*.cpp | head
```

Expected: no species method currently reads `ψ` during stress-energy accumulation. If some species does (e.g. `ScalarFieldSpecies::Delta`), confirm the current code computes ψ before calling Delta on that species — the "must be below all species with non-zero shear" comment at line 6411. If so, the two-pass design is correct. If ψ is read inside `RhoPlusPShear` for any species, the design falls to three passes.

- [ ] **Step 2: Inspect how ψ is derived today**

```bash
grep -n "index_pt_phi\|scalar_ctx.psi\|perturb_einstein" source/perturbations_module.cpp | head -30
```

Expected: `perturb_einstein` (a separate function) derives ψ from the Einstein equations using accumulated `rho_plus_p_shear`. In the current code, ψ is NOT derived inside `perturb_total_stress_energy` — it's derived later in `perturb_einstein`. Consequence: the two-pass is conceptually:

```
Pass 1: shear
  -> perturb_einstein will use this shear to derive ψ (happens outside this function)
Pass 2: delta_rho, theta, delta_p, rho_plus_p_tot
```

The scalar field's `Delta` that depends on ψ is the one exception — today it reads `ppw->pvecmetric[ppw->index_mt_alpha]` or similar. Confirm `ScalarFieldSpecies::Delta` in `species/scalar_field.cpp` does not reach for ψ directly; it reads from `scalar_ctx` or `pvecmetric`.

If the investigation shows ψ is *not* needed during pass 2 (because ScalarField uses pvecmetric/scalar_ctx quantities already populated), skip the explicit "derive ψ" step — pass 1 just accumulates shear and exits; pass 2 runs directly.

- [ ] **Step 3: Write the new function body**

Open `source/perturbations_module.cpp`. Locate `perturb_total_stress_energy` (around line 6060). Replace the body from just after the preamble setup (approximately line 6191) through the end of the fluid block (line 6553) with:

```cpp
    // (b) Pass 1: shear (feeds perturb_einstein's ψ derivation downstream)
    ppw->rho_plus_p_shear = 0.;
    for (auto& sp : all_species_) {
      if (sp->name() == "Fluid" && pba->use_ppf == _TRUE_) continue;  // PPF handled below
      ppw->rho_plus_p_shear += sp->RhoPlusPShear(ppw->pv, y, ppw->pvecback, ppw);
    }

    // (c) Pass 2: delta_rho, theta, delta_p, rho_plus_p_tot, matter tally
    ppw->delta_rho        = 0.;
    ppw->rho_plus_p_theta = 0.;
    ppw->delta_p          = 0.;
    ppw->rho_plus_p_tot   = 0.;
    double delta_rho_m = 0., rho_m = 0., rho_plus_p_theta_m = 0., rho_plus_p_m = 0.;

    for (auto& sp : all_species_) {
      if (sp->name() == "Fluid" && pba->use_ppf == _TRUE_) continue;  // PPF handled below

      const double rho        = sp->Rho(ppw->pvecback);
      const double rho_plus_p = rho + sp->P(ppw->pvecback);
      ppw->delta_rho        += rho        * sp->Delta(ppw->pv, y, ppw->pvecback, ppw);
      ppw->rho_plus_p_theta += rho_plus_p * sp->Theta(ppw->pv, y, ppw->pvecback, ppw);
      ppw->delta_p          += sp->DeltaP(ppw->pv, y, ppw->pvecback, ppw);
      ppw->rho_plus_p_tot   += rho_plus_p;

      if (sp->IsMatterSpecies()) {
        delta_rho_m        += sp->MatterRhoDelta(ppw->pv, y, ppw->pvecback, ppw);
        rho_m              += sp->MatterRho(ppw->pvecback);
        rho_plus_p_theta_m += sp->MatterRhoPlusPTheta(ppw->pv, y, ppw->pvecback, ppw);
        rho_plus_p_m       += sp->MatterRhoPlusP(ppw->pvecback);
      }
    }

    // (d) PPF fluid: one species only, depends on everything else; kept special.
    if (all_species_.count("Fluid") && pba->use_ppf == _TRUE_) {
      auto& fluid = static_cast<FluidSpecies&>(*all_species_.at("Fluid"));
      double w_fld, dw_over_da_fld, integral_fld;
      class_call(fluid.ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld),
                 background_module_->error_message_, error_message_);
      fluid.ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
      ppw->delta_rho        += ppw->delta_rho_fld;
      ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
      ppw->delta_p          += ppw->delta_p_fld;
      ppw->rho_plus_p_tot   += (1. + w_fld) * fluid.Rho(ppw->pvecback);
    }

    // (e) delta_cb / theta_cb wrap-up — unchanged logic, just reorganized.
    if ((has_source_delta_m_ == _TRUE_) && (has_source_delta_cb_ == _TRUE_))
      ppw->delta_cb = delta_rho_m / rho_m;
    if (((has_source_delta_m_ == _TRUE_) || (has_source_theta_m_ == _TRUE_)) &&
        ((has_source_delta_cb_ == _TRUE_) || (has_source_theta_cb_ == _TRUE_)))
      ppw->theta_cb = rho_plus_p_theta_m / rho_plus_p_m;
```

**Important:** the code at lines 6557–6570 (storing `ppw->delta_m` / `ppw->theta_m` for downstream use) must remain. Preserve that block.

- [ ] **Step 4: Confirm `scalar_ctx.a` and `scalar_ctx.k2` are populated before the loops**

The existing preamble populates these at `perturbations_module.cpp:6184–6189`. Keep that as-is.

- [ ] **Step 5: Delete the now-dead arrays**

```bash
grep -n "delta_ncdm\[\|theta_ncdm\[\|shear_ncdm\[" source/*.cpp source/*.h species/*.cpp species/*.h | grep -v "^source/perturbations_module.cpp:7"
```

Remaining consumers should be only the output-path at 7337–7345 (local vector, not ppw) and title strings. If `ppw->delta_ncdm`, `ppw->theta_ncdm`, `ppw->shear_ncdm` no longer have any consumer, delete their declarations in `source/perturbations.h` (and wherever they are allocated in `perturb_workspace_init`).

```bash
grep -n "delta_ncdm\b\|theta_ncdm\b\|shear_ncdm\b" source/perturbations.h source/perturbations_module.cpp
```

Delete the field declarations and their allocations. The vector `q_size_ncdm` elsewhere is a different array — don't touch it.

- [ ] **Step 6: Build, smoke test, fast suite**

```bash
make clean && make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

Expected: all tests pass. If any scenario differs, the refactor has changed numerics — bisect against the previous commit and find the accumulation order or guard condition that was dropped.

- [ ] **Step 7: Commit**

```bash
git add source/perturbations_module.cpp source/perturbations.h
git -c commit.gpgsign=false commit -m "perturbations: collapse perturb_total_stress_energy to a species-dispatch loop

The function goes from ~400 lines of per-species blocks to a two-pass loop
over all_species_, with PPF fluid as the single explicit special case after
the main loop. All species-local physics (fluid EoS, PPF, IDR TCA shear,
DNCDM rescaling) now lives inside the respective species. No new abstractions
— just calls through the existing BaseSpecies interface.

Delete now-unused ppw->delta_ncdm / theta_ncdm / shear_ncdm arrays.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 9: Full reference suite + PR

- [ ] **Step 1: Run the full reference comparison**

```bash
cd /Users/au192734/Projects/class_claude
make clean && make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
cd python
TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tee /tmp/refsuite.log | tail -50
cd ..
```

Expected: all scenarios pass, bit-identical against the reference. If any scenario fails:
- Inspect `/tmp/refsuite.log` for the delta.
- If the delta is within the built-in tolerance but not bit-identical, it's almost certainly an order-of-accumulation change in pass 2. Document in the PR body; optionally add a tolerance note. Do not block on it.
- If the delta is outside tolerance, it's a bug. Bisect the plan commits: the breaking commit is the one to revisit.

- [ ] **Step 2: Clean up any stragglers**

```bash
grep -n "RescaledNCDMPerturbations\|tca_shear_idm_dr\|delta_ncdm\[" source/*.cpp source/*.h species/*.cpp species/*.h
```

Expected: only legitimate residuals (e.g. `tca_shear_idm_dr` as a local var in `idm_dr_idr_species.cpp:261`, NCDM output titles in 7337+). Any remaining module-level `RescaledNCDMPerturbations` delegator and `ppw->tca_shear_idm_dr` usage should be gone.

- [ ] **Step 3: Final build-system check**

No new files were added in this plan. Confirm:

```bash
git diff master --stat
```

Expected: only `.h`/`.cpp` modifications in `source/`, `species/`; no new files. If any file was added, update `Makefile`, `setup.py`, `CLASS.xcodeproj/project.pbxproj` accordingly.

- [ ] **Step 4: Push branch and open PR**

```bash
git push -u origin <branch-name>
gh pr create \
  --title "Rewrite perturb_total_stress_energy as a species-dispatch loop" \
  --body "$(cat <<'EOF'
## Summary

Collapses `PerturbationsModule::perturb_total_stress_energy` from ~400 lines of per-species blocks to a two-pass loop over `all_species_`. Species-local physics (fluid equation of state, PPF, IDR TCA shear, DNCDM rescaling) moves into the respective species classes. PPF stays as the one explicit special case in the module — one species allowed, depends on everything else, called after the main loop.

## Changes

- `BaseSpecies`: add `IsMatterSpecies()` and four `Matter*` virtuals with energy-type-driven defaults.
- `CompositeSpecies`: override `Matter*` to sum matter children.
- `IDM_DRSpecies`: `IsMatterSpecies() override { return false; }` (preserves pre-existing asymmetry in the matter tally; follow-up issue tracks the asymmetry itself).
- `FluidSpecies`: owns `ComputeWFld` (body moved from `BackgroundModule::background_w_fld`, which becomes a delegator), `ComputePpf` (PPF body), and a correct `DeltaP` for the non-PPF case.
- `DNCDMSpecies`: owns `RescaledPerturbations` (moved from `PerturbationsModule::RescaledNCDMPerturbations`); overrides `Delta/Theta/RhoPlusPShear` to return rescaled values.
- `IDRSpecies`: owns `TcaShearIdr` (moved from the stress-energy preamble); `ppw->tca_shear_idm_dr` side channel deleted.
- `ppw->delta_ncdm[n]` / `theta_ncdm[n]` / `shear_ncdm[n]` arrays deleted.
- New `BaseSpecies::SetThermodynamicsModule` hook (mirrors the existing `SetBackgroundModule`).

## Design doc

See [`docs/superpowers/specs/2026-04-22-stress-energy-loop-design.md`](docs/superpowers/specs/2026-04-22-stress-energy-loop-design.md).

## Test plan

- [x] `make clean && make class -j` clean
- [x] `./class explanatory.ini` succeeds
- [x] `TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py` passes
- [x] `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py` passes (bit-identical / within tolerance; see commit history for any documented drift)

## Follow-ups (not in this PR)

- The IDM_DR matter-tally asymmetry (why is IDM_DR excluded but IDM_DRMD / DCDM / NCDM included?).
- Remaining `static_cast<XSpecies&>(*all_species_.at("X"))` downcasts in `background_module`, `thermodynamics_module`, `input_module`.
- `init()`/`free()` lifecycle removal across modules.
- `cosmology.cpp` raw `new` → `std::make_unique`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Self-review summary

- Every task has exact file paths and code blocks where code changes are required.
- Each commit builds and passes at minimum the `explanatory.ini` smoke test; most also pass `TEST_LEVEL=1`.
- The big commit (Task 8) lands after all preparatory commits so that a bisect identifying Task 8 as the break point immediately tells you the rewrite introduced drift — not some earlier plumbing.
- No `TBD`, no "see above," no unspecified "handle edge cases." Where a verbatim copy of existing code is needed (PPF body, w_fld body), the plan calls out the exact line range and the variable renames required.
- The one place with flexible wording is Task 8 Step 2 — "investigate whether ψ is needed during pass 2." This is unavoidable without reading `perturb_einstein` in depth; if the investigation reveals three passes are needed, the function body in Step 3 expands by one loop. Noted as a risk in the spec.
