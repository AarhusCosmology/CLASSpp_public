# Newtonian-gauge transformation via species dispatch — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the ~280-line `if (ppt->gauge == newtonian)` initial-condition block in `source/perturbations_module.cpp` (currently lines 4633–4935) with species-polymorphic dispatch, so adding a species requires no edits to this block.

**Architecture:** Derive the density-contrast gauge shift from the background continuity equation via a new `BaseSpecies::RhoDotOverRho` accessor (single source of truth, shared with `BackgroundDerivs`). The transformation becomes two universal lines — `delta += (ρ̇/ρ)·α`, `theta += k²·α` — applied through a new `PerturbSynchronousToNewtonian` hook (base default + a protected helper); special species (NCDM, DarkRadiation, ScalarField) and composites override. The `alpha` accumulation stays module-level but becomes a downcast-free `Rho·Delta` reduction, which also fixes the `fraccdm` quirk.

**Tech Stack:** C++17, the CLASSpp species-plugin architecture, `make class` / `make classy`, pytest harness in `python/test_class.py` (cross-gauge + reference comparisons), `test/scenarios/compare_tol.py`.

**Design doc:** `docs/superpowers/specs/2026-05-24-newtonian-gauge-dispatch-design.md`

**Behavior changes (all intentional, all verified):**
1. **Fluid sign fix** — non-PPF fluid `delta` shift flips to match `ρ̇/ρ` and public CLASS (current `+= 3(1+w)ℋα` is a confirmed bug).
2. **`fraccdm` fix** — every matter species contributes `Rho·Delta` to `alpha` (not just literal "CDM").
3. **`alpha` now includes IDR and full NCDM** — the reduction sums all non-dark-energy species (current code lumped NCDM into `fracnu·delta_ur` and omitted IDR). Verified within tolerance.
4. **DNCDM_DR / combined-DR Newtonian re-seed** — each composite now re-seeds its own dark radiation (current block e is `count("DCDM_DR")`-only and shares a `delta_dr` local across channels). Only affects the currently-untested DNCDM/combined Newtonian paths; DCDM-only is preserved exactly.

---

## File Structure

| File | Responsibility / change |
|---|---|
| `species/base_species.h` | New `RhoDotOverRho` + `PerturbSynchronousToNewtonian` virtuals; protected `ApplyFluidLikeNewtonianShift` helper. |
| `species/dcdm.{h,cpp}` | Override `RhoDotOverRho` (decay sink); universal `PerturbSynchronousToNewtonian`. |
| `species/photons.{h,cpp}`, `baryons.{h,cpp}`, `cdm.{h,cpp}`, `fluid.{h,cpp}`, `ultra_relativistic.{h,cpp}` | Universal `PerturbSynchronousToNewtonian`. |
| `species/idm_dr.{h,cpp}`, `idr.{h,cpp}`, `idm_drmd.{h,cpp}`, `idr_drmd.{h,cpp}` | Universal `PerturbSynchronousToNewtonian` (composite members). |
| `species/scalar_field.{h,cpp}` | Override `PerturbSynchronousToNewtonian` (φ, φ′). |
| `species/ncdm_species.{h,cpp}` | `NCDMBaseSpecies::PerturbSynchronousToNewtonian` (per-q) + protected virtual `Dlnf0Dlnq`. |
| `species/dncdm_species.{h,cpp}` | Override `Dlnf0Dlnq` (time-dependent). |
| `species/dark_radiation_species.{h,cpp}` | Public `PerturbNewtonianReseed(layout, y, ctx, decay_corr)` helper (re-seed F-hierarchy). |
| `species/dcdm_dr_species.{h,cpp}`, `dncdm_dr_species.{h,cpp}`, `idm_dr_idr_species.{h,cpp}`, `idm_drmd_idr_drmd_species.{h,cpp}` | `PerturbSynchronousToNewtonian` delegating to members; composites compute their DR `decay_corr`. |
| `source/perturbations_module.cpp` | Replace block (d) shifts + block (e) re-seed with the dispatch loop (Task 7); replace the `alpha` accumulation with the downcast-free reduction (Task 8). |
| `test/scenarios/*.ini` | New Newtonian-gauge test scenarios + baselines. |

No new source files → no build-system list changes (Makefile / setup.py / pbxproj untouched).

---

### Task 1: Branch + capture baselines from current code

**Files:**
- Create: `test/scenarios/gauge_lcdm.ini`, `gauge_ncdm.ini`, `gauge_dcdm.ini`, `gauge_idmdr.ini`, `gauge_fluid.ini`, `gauge_scf.ini`

- [ ] **Step 1: Create a feature branch**

```bash
cd /Users/au192734/Projects/class_claude
git checkout -b newtonian-gauge-dispatch
```

- [ ] **Step 2: Create the test scenarios** (one base `.ini` per species family; each gets run in both gauges by Step 4)

Create `test/scenarios/gauge_lcdm.ini`:
```ini
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
N_ur = 3.044
```

Create `test/scenarios/gauge_ncdm.ini`:
```ini
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
N_ur = 2.0308
N_ncdm = 1
m_ncdm = 0.06
```

Create `test/scenarios/gauge_dcdm.ini`:
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_dcdmdr = 0.12
Gamma_dcdm = 10.0
```

Create `test/scenarios/gauge_idmdr.ini`:
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_idm_dr = 0.12
xi_idr = 0.3
a_idm_dr = 1e-4
```

Create `test/scenarios/gauge_fluid.ini` (non-PPF fluid — the bug-fix target):
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_Lambda = 0
fluid_equation_of_state = CLP
w0_fld = -0.9
wa_fld = 0.1
cs2_fld = 1.0
use_ppf = no
```

Create `test/scenarios/gauge_scf.ini`:
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_Lambda = 0
Omega_fld = 0
Omega_scf = -1
attractor_ic_scf = yes
scf_parameters = 10.0, 0.0, 0.0, 0.0, 0.0
```

- [ ] **Step 3: Build the current (baseline) binary**

Run: `make class -j4`
Expected: `gcc ... -o class` with no errors; `./class` exists.

- [ ] **Step 4: Generate baseline outputs in both gauges**

```bash
mkdir -p /tmp/gauge_base
for s in gauge_lcdm gauge_ncdm gauge_dcdm gauge_idmdr gauge_fluid gauge_scf; do
  for g in synchronous newtonian; do
    d=/tmp/gauge_base/${s}_${g}; mkdir -p "$d"
    printf 'gauge = %s\nroot = %s/out_\n' "$g" "$d" > /tmp/run.ini
    cat test/scenarios/${s}.ini >> /tmp/run.ini
    ./class /tmp/run.ini
  done
done
ls /tmp/gauge_base/gauge_lcdm_newtonian/
```
Expected: each directory contains `out_cl.dat` (and `out_pk.dat` where mPk requested). These are the frozen baselines for later comparison.

- [ ] **Step 5: Record current synchronous-vs-Newtonian agreement** (so we know the starting point, especially the fluid)

```bash
for s in gauge_lcdm gauge_ncdm gauge_dcdm gauge_idmdr gauge_fluid gauge_scf; do
  echo "== $s : synchronous vs newtonian =="
  python test/scenarios/compare_tol.py /tmp/gauge_base/${s}_synchronous /tmp/gauge_base/${s}_newtonian "out_*.dat"
done
```
Expected: most `OK`; **`gauge_fluid` likely FAILs or shows elevated `worst_vs_colpeak`** — that is the pre-existing fluid sign bug. Note the numbers; after the fix this scenario must improve to `OK`.

- [ ] **Step 6: Commit the scenarios**

```bash
git add test/scenarios/gauge_*.ini docs/superpowers/specs/2026-05-24-newtonian-gauge-dispatch-design.md docs/superpowers/plans/2026-05-24-newtonian-gauge-dispatch.md
git commit -m "test: add Newtonian-gauge scenarios + commit gauge-dispatch design/plan"
```

---

### Task 2: Add `RhoDotOverRho` accessor (base default + DCDM override)

**Files:**
- Modify: `species/base_species.h` (after `CopyPerturbationsAcrossSwitch`, ~line 455)
- Modify: `species/dcdm.h`, `species/dcdm.cpp`

- [ ] **Step 1: Add the virtual to `BaseSpecies`**

In `species/base_species.h`, immediately after the `CopyPerturbationsAcrossSwitch` declaration (the `const {}` ending ~line 455), insert:

```cpp
  // ── Synchronous → Newtonian gauge transformation ──────────────────────────
  /**
   * Log conformal-time derivative of the background density, ρ̇/ρ (ρ̇ ≡ dρ̄/dτ).
   * Single source of truth for the density-contrast gauge shift; mirrors the
   * continuity equation used in BackgroundDerivs. Default: ρ̇/ρ = -3ℋ(Rho+P)/Rho.
   * Species with extra source/sink terms (e.g. decay) override this.
   * @param a_prime_over_a  ℋ = a'/a (base does not store H/a indices).
   */
  virtual double RhoDotOverRho(const double* pvecback, double a_prime_over_a) const {
    return -3. * a_prime_over_a * (Rho(pvecback) + P(pvecback)) / Rho(pvecback);
  }
```

- [ ] **Step 2: Declare the DCDM override**

In `species/dcdm.h`, next to the existing `double DpDloga(...) const override;` (~line 45), add:

```cpp
  double RhoDotOverRho(const double* pvecback, double a_prime_over_a) const override;
```

- [ ] **Step 3: Implement the DCDM override**

In `species/dcdm.cpp`, after `DpDloga` (~line 49), add:

```cpp
double DCDMSpecies::RhoDotOverRho(const double* pvecback, double a_prime_over_a) const {
  // ρ̇/ρ = -a(3H + Γ) = -3ℋ - aΓ  (same physics as BackgroundDerivs).
  const double a = pvecback[bgm_->index_bg_a_];
  return -3. * a_prime_over_a - a * pba_.Gamma_dcdm;
}
```

- [ ] **Step 4: Build**

Run: `make class -j4`
Expected: compiles clean (method is not yet called → no behavior change).

- [ ] **Step 5: Commit**

```bash
git add species/base_species.h species/dcdm.h species/dcdm.cpp
git commit -m "feat(species): add RhoDotOverRho continuity accessor (base + DCDM decay override)"
```

---

### Task 3: Add the transform hook + helper to `BaseSpecies`

**Files:**
- Modify: `species/base_species.h`

- [ ] **Step 1: Add the public virtual**

In `species/base_species.h`, directly below the `RhoDotOverRho` virtual added in Task 2, insert:

```cpp
  /**
   * Transform this species' own perturbation variables from synchronous to
   * Newtonian gauge, in place in y[]. Called by the module after alpha is known
   * (ctx.alpha / ctx.alpha_prime filled). The synchronous IC is already in y[]
   * (ApplyInitialConditions runs in both gauges), so this is a pure shift/re-seed.
   * Default: no-op (species with no perturbed variables, e.g. cosmological constant).
   */
  virtual void PerturbSynchronousToNewtonian(const PerturbLayout& /*layout*/,
                                             double* /*y*/,
                                             const PerturbIcContext& /*ctx*/) {}
```

- [ ] **Step 2: Add the protected helper**

In `species/base_species.h`, inside the existing `protected:` section (after the `BaseSpecies(std::string, EnergyType)` constructor, ~line 545), insert:

```cpp
  /**
   * Universal fluid-like gauge shift: delta += (ρ̇/ρ)·alpha ; theta += k²·alpha
   * (shear, l3 and higher moments are gauge-invariant). Fluid-like species call
   * this from PerturbSynchronousToNewtonian, supplying their own delta/theta slots.
   */
  void ApplyFluidLikeNewtonianShift(double* y, int idx_delta, int idx_theta,
                                    const double* pvecback, const PerturbIcContext& ctx) const {
    if (idx_delta >= 0)
      y[idx_delta] += RhoDotOverRho(pvecback, ctx.a_prime_over_a) * ctx.alpha;
    if (idx_theta >= 0)
      y[idx_theta] += ctx.k * ctx.k * ctx.alpha;
  }
```

- [ ] **Step 3: Build**

Run: `make class -j4`
Expected: compiles clean (still uncalled).

- [ ] **Step 4: Commit**

```bash
git add species/base_species.h
git commit -m "feat(species): add PerturbSynchronousToNewtonian hook + fluid-like shift helper"
```

---

### Task 4: Universal override for standalone radiation/matter species

**Files:**
- Modify: `species/photons.{h,cpp}`, `baryons.{h,cpp}`, `cdm.{h,cpp}`, `fluid.{h,cpp}`, `ultra_relativistic.{h,cpp}`

For each species the `.h` gets one declaration and the `.cpp` one definition. The body is identical except the layout type.

- [ ] **Step 1: Photons** — declare in `species/photons.h` next to its other layout-based overrides:

```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/photons.cpp`:
```cpp
void PhotonsSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                   double* y,
                                                   const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```

- [ ] **Step 2: Baryons** — declare in `species/baryons.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/baryons.cpp`:
```cpp
void BaryonsSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                   double* y,
                                                   const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```

- [ ] **Step 3: CDM** — declare in `species/cdm.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/cdm.cpp`:
```cpp
void CDMSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                               double* y,
                                               const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```
(CDM `theta` starts at 0 — `y_storage` is zero-initialised and CDM's `ApplyInitialConditions` does not set `theta` — so `+=` reproduces the old `theta = k²α`.)

- [ ] **Step 4: Fluid** — declare in `species/fluid.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/fluid.cpp`:
```cpp
void FluidSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                 double* y,
                                                 const PerturbIcContext& ctx) {
  // Non-PPF only: in PPF mode idx_delta/idx_theta are unregistered (idx == -1) and
  // the helper is a no-op. delta uses the universal ρ̇/ρ shift = -3(1+w)ℋα, which
  // corrects the historical opposite-sign bug.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```

- [ ] **Step 5: UltraRelativistic** — declare in `species/ultra_relativistic.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/ultra_relativistic.cpp`:
```cpp
void UltraRelativisticSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                             double* y,
                                                             const PerturbIcContext& ctx) {
  // delta/theta shift; shear (idx_shear) and l3 (idx_l3) are gauge-invariant.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```

- [ ] **Step 6: Build**

Run: `make class -j4`
Expected: compiles clean (still uncalled by the module).

- [ ] **Step 7: Commit**

```bash
git add species/photons.* species/baryons.* species/cdm.* species/fluid.* species/ultra_relativistic.*
git commit -m "feat(species): universal Newtonian-gauge transform for photons/baryons/cdm/fluid/UR"
```

---

### Task 5: ScalarField and NCDM overrides

**Files:**
- Modify: `species/scalar_field.{h,cpp}`, `species/ncdm_species.{h,cpp}`, `species/dncdm_species.{h,cpp}`

- [ ] **Step 1: ScalarField — declare** in `species/scalar_field.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```

- [ ] **Step 2: ScalarField — define** in `species/scalar_field.cpp` (mirrors the old inline block at perturbations_module.cpp:4792–4810; synchronous IC is φ=φ′=0, so `+=` == overwrite):
```cpp
void ScalarFieldSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                       double* y,
                                                       const PerturbIcContext& ctx) {
  const auto& l            = static_cast<const PerturbLayout&>(base);
  const double* pvecback   = ctx.ppw->pvecback;
  const BackgroundModule* bgm = ctx.p_mod->GetBackgroundModule();
  const double phi_prime   = pvecback[bgm->index_bg_phi_prime_scf_];
  const double phi_scf     = pvecback[bgm->index_bg_phi_scf_];
  if (l.idx_phi >= 0)
    y[l.idx_phi] += ctx.alpha * phi_prime;
  if (l.idx_phi_prime >= 0)
    y[l.idx_phi_prime] += (-2. * ctx.a_prime_over_a * ctx.alpha * phi_prime -
                           ctx.a * ctx.a * bgm->dV_scf(phi_scf) * ctx.alpha +
                           phi_prime * ctx.alpha_prime);
}
```
(Confirm the `BackgroundModule` accessor name; the old block used `background_module_->index_bg_phi_prime_scf_`, `index_bg_phi_scf_`, `dV_scf(...)`. `ctx.p_mod->GetBackgroundModule()` exists — it is used in `cdm.cpp:239`.)

- [ ] **Step 3: NCDM — add the `Dlnf0Dlnq` protected virtual.** In `species/ncdm_species.h`, in the `protected:` area of `NCDMBaseSpecies`, add:
```cpp
  /** dln f0 / dln q at the q-node. Default: the static tabulated value (stable NCDM).
      DNCDMSpecies overrides with the time-dependent background value. */
  virtual double Dlnf0Dlnq(int index_q, const double* /*pvecback*/) const {
    return dlnf0_dlnq()[index_q];
  }
```
And declare the transform override on `NCDMBaseSpecies` (public):
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```

- [ ] **Step 4: NCDM — define the transform** in `species/ncdm_species.cpp` (mirrors perturbations_module.cpp:4898–4914; re-seed per-q from the gauge-shifted free-streaming radiation IC, `delta_ur` shifting as radiation `-4ℋα`):
```cpp
void NCDMBaseSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;
  const double* pvecback = ctx.ppw->pvecback;
  const double delta_ur  = ctx.delta_ur - 4. * ctx.a_prime_over_a * ctx.alpha;  // radiation shift
  const double theta_ur  = ctx.theta_ur + ctx.k * ctx.k * ctx.alpha;
  const int lmax = layout.l_max;
  for (int index_q = 0; index_q < layout.q_size; ++index_q) {
    const int idx           = layout.index_per_q[index_q];
    const double q          = q()[index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * M() * M());
    const double dlnf0_dlnq = Dlnf0Dlnq(index_q, pvecback);
    y[idx + 0] = -0.25 * delta_ur * dlnf0_dlnq;
    if (lmax >= 1)
      y[idx + 1] = -epsilon / 3. / q / ctx.k * theta_ur * dlnf0_dlnq;
    if (lmax >= 2)
      y[idx + 2] = -0.5 * ctx.shear_ur * dlnf0_dlnq;  // shear/l3 gauge-invariant (unshifted)
    if (lmax >= 3)
      y[idx + 3] = -0.25 * ctx.l3_ur * dlnf0_dlnq;
  }
}
```
(Use the accessors `q()`, `M()`, `dlnf0_dlnq()` — confirmed used by the old module block via the `NCDMBaseSpecies*`.)

- [ ] **Step 5: DNCDM — override `Dlnf0Dlnq`.** In `species/dncdm_species.h` (protected), add:
```cpp
  double Dlnf0Dlnq(int index_q, const double* pvecback) const override {
    return pvecback[bg_dlnfdlnq_index() + index_q];  // time-dependent for decaying NCDM
  }
```

- [ ] **Step 6: Build**

Run: `make class -j4`
Expected: compiles clean.

- [ ] **Step 7: Commit**

```bash
git add species/scalar_field.* species/ncdm_species.* species/dncdm_species.*
git commit -m "feat(species): Newtonian-gauge transform for ScalarField + NCDM/DNCDM (per-q)"
```

---

### Task 6: DarkRadiation re-seed helper + composite delegation

**Files:**
- Modify: `species/dark_radiation_species.{h,cpp}`, `species/dcdm_dr_species.{h,cpp}`, `species/dncdm_dr_species.{h,cpp}`, `species/idm_dr.{h,cpp}`, `species/idr.{h,cpp}`, `species/idm_dr_idr_species.{h,cpp}`, `species/idm_drmd.{h,cpp}`, `species/idr_drmd.{h,cpp}`, `species/idm_drmd_idr_drmd_species.{h,cpp}`

- [ ] **Step 1: DarkRadiation re-seed helper — declare** in `species/dark_radiation_species.h` (public):
```cpp
  /** Re-seed the DR multipole hierarchy from the gauge-shifted IC. decay_corr is
      supplied by the owning composite (= aΓ·ρ_parent/ρ_dr; 0 when ρ_dr == 0). */
  void PerturbNewtonianReseed(const PerturbLayout& layout, double* y,
                              const PerturbIcContext& ctx, double decay_corr) const;
```

- [ ] **Step 2: DarkRadiation re-seed helper — define** in `species/dark_radiation_species.cpp` (mirrors `ApplyInitialConditions` at lines 161–182 + the old block e at perturbations_module.cpp:4922–4934, with shifted `delta_dr`/`theta_ur`):
```cpp
void DarkRadiationSpecies::PerturbNewtonianReseed(const PerturbLayout& layout, double* y,
                                                  const PerturbIcContext& ctx,
                                                  double decay_corr) const {
  if (layout.idx_F0 < 0)
    return;
  const double* pvecback = ctx.ppw->pvecback;
  const double r_dr      = std::pow(std::pow(ctx.a / pba_->a_today, 2) / pba_->H0, 2) *
                           pvecback[index_bg_rho_];
  const double delta_dr  = ctx.delta_dr + (-4. * ctx.a_prime_over_a + decay_corr) * ctx.alpha;
  const double theta_ur  = ctx.theta_ur + ctx.k * ctx.k * ctx.alpha;
  y[layout.idx_F0 + 0] = delta_dr * r_dr;
  if (layout.l_max >= 1)
    y[layout.idx_F0 + 1] = 4. / (3. * ctx.k) * theta_ur * r_dr;
  if (layout.l_max >= 2)
    y[layout.idx_F0 + 2] = 2. * ctx.shear_ur * r_dr;  // shear/l3 gauge-invariant
  if (layout.l_max >= 3)
    y[layout.idx_F0 + 3] = ctx.l3_ur * r_dr;
}
```

- [ ] **Step 3: DCDM_DR composite — declare** in `species/dcdm_dr_species.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
**Define** in `species/dcdm_dr_species.cpp`:
```cpp
void DCDM_DR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dcdm_->PerturbSynchronousToNewtonian(my.dcdm, y, ctx);  // matter+decay via RhoDotOverRho
  const double* pvecback  = ctx.ppw->pvecback;
  const double rho_dr     = dr_sp_->Rho(pvecback);
  const double decay_corr = (rho_dr > 0.) ?
      ctx.a * pba_->Gamma_dcdm * dcdm_->Rho(pvecback) / rho_dr : 0.;
  dr_sp_->PerturbNewtonianReseed(my.dr, y, ctx, decay_corr);
}
```

- [ ] **Step 4: DNCDM_DR composite — declare** in `species/dncdm_dr_species.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
**Define** in `species/dncdm_dr_species.cpp` (decay rate/parent density from the DNCDM child — fixes the old `count("DCDM_DR")`-only re-seed):
```cpp
void DNCDM_DR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                     double* y,
                                                     const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dncdm_->PerturbSynchronousToNewtonian(my.dncdm, y, ctx);  // NCDMBaseSpecies per-q transform
  const double* pvecback  = ctx.ppw->pvecback;
  const double rho_dr     = dr_sp_->Rho(pvecback);
  const double decay_corr = (rho_dr > 0.) ?
      ctx.a * dncdm_->Gamma() * dncdm_->Rho(pvecback) / rho_dr : 0.;
  dr_sp_->PerturbNewtonianReseed(my.dr, y, ctx, decay_corr);
}
```

- [ ] **Step 5: IDM_DR + IDR members — universal overrides.** Declare in `species/idm_dr.h` and `species/idr.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/idm_dr.cpp`:
```cpp
void IDM_DRSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                  double* y,
                                                  const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```
Define in `species/idr.cpp`:
```cpp
void IDRSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                               double* y,
                                               const PerturbIcContext& ctx) {
  // delta/theta shift; shear/l3 gauge-invariant. idm_dr's synchronous theta is
  // already theta_ur (tight-coupling lock), so the universal += reproduces it.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```

- [ ] **Step 6: IDM_DR_IDR composite — declare** in `species/idm_dr_idr_species.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
**Define** in `species/idm_dr_idr_species.cpp`:
```cpp
void IDM_DR_IDR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                       double* y,
                                                       const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  idm_dr_->PerturbSynchronousToNewtonian(my.idm_dr, y, ctx);
  idr_->PerturbSynchronousToNewtonian(my.idr, y, ctx);
}
```

- [ ] **Step 7: IDM_DRMD + IDR_DRMD members — universal overrides** (dead in Newtonian gauge per the `has_idr_drmd` `class_test`, implemented for consistency). Declare in `species/idm_drmd.h` and `species/idr_drmd.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
Define in `species/idm_drmd.cpp`:
```cpp
void IDM_DRMDSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```
Define in `species/idr_drmd.cpp`:
```cpp
void IDR_DRMDSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}
```

- [ ] **Step 8: IDM_DRMD_IDR_DRMD composite — declare** in `species/idm_drmd_idr_drmd_species.h`:
```cpp
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
```
**Define** in `species/idm_drmd_idr_drmd_species.cpp`:
```cpp
void IDM_DRMD_IDR_DRMD_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                              double* y,
                                                              const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  idm_drmd_->PerturbSynchronousToNewtonian(my.idm_drmd, y, ctx);
  idr_drmd_->PerturbSynchronousToNewtonian(my.idr_drmd, y, ctx);
}
```

- [ ] **Step 9: Build**

Run: `make class -j4`
Expected: compiles clean (all overrides still uncalled by the module).

- [ ] **Step 10: Commit**

```bash
git add species/dark_radiation_species.* species/dcdm_dr_species.* species/dncdm_dr_species.* \
        species/idm_dr.* species/idr.* species/idm_dr_idr_species.* \
        species/idm_drmd.* species/idr_drmd.* species/idm_drmd_idr_drmd_species.*
git commit -m "feat(species): DR re-seed helper + composite Newtonian-gauge delegation"
```

---

### Task 7: Switch the module to dispatch (block d shifts + block e re-seed)

**Files:**
- Modify: `source/perturbations_module.cpp:4724-4935`

- [ ] **Step 1: Replace the inline shifts and re-seed with the dispatch loop.**

In `source/perturbations_module.cpp`, the block currently runs:
```cpp
      ppw->pv->y[ppw->pv->index_pt_phi] = eta - a_prime_over_a * alpha;
      ppw->pv->y[g_lay_ic.idx_delta] -= 4. * a_prime_over_a * alpha;
      ... (all inline shift blocks, lines 4726-4838) ...
    } /* end of gauge transformation to newtonian gauge */

    /** - (e) ... */
    if (ppt->gauge == newtonian) {
      ... (block e re-seed, lines 4843-4934) ...
    }
```
Replace everything from line 4726 (`ppw->pv->y[g_lay_ic.idx_delta] -= ...`) through the end of block (e) (line 4935, the closing `}` of `if (ppt->gauge == newtonian)`) with the single dispatch loop:
```cpp
      /** - (d/e) Per-species synchronous->Newtonian gauge transformation. Each
          species shifts its own variables in place (the synchronous IC is already
          in y[] from ApplyInitialConditions, which runs in both gauges). This
          replaces the former inline per-species shift blocks and relativistic re-seed. */
      ic_ctx.alpha       = alpha;
      ic_ctx.alpha_prime = 0.;  // historical: alpha_prime was hardcoded 0 in the scalar-field block
      {
        size_t i = 0;
        for (auto& sp : all_species_) {
          sp->PerturbSynchronousToNewtonian(*ppw->pv->species_layouts[i], ppw->pv->y, ic_ctx);
          ++i;
        }
      }
    } /* end of gauge transformation to newtonian gauge */
```
Keep the metric line `ppw->pv->y[ppw->pv->index_pt_phi] = eta - a_prime_over_a * alpha;` (line 4724) immediately before this. **Do not** touch the `alpha` accumulation above it (lines 4657–4722) — that is Task 8. The local `delta_ur`/`theta_ur`/… and the `g_lay_ic`/`b_lay_ic` reads in the accumulation stay for now.

- [ ] **Step 2: Build**

Run: `make class -j4`
Expected: compiles clean. (If unused-variable warnings appear for `g_lay_ic`/`b_lay_ic`, ignore — they are still read by the accumulation in 4657–4722.)

- [ ] **Step 3: Regenerate Newtonian outputs and compare to baseline (non-fluid must be ~identical)**

```bash
mkdir -p /tmp/gauge_new
for s in gauge_lcdm gauge_ncdm gauge_dcdm gauge_idmdr gauge_scf; do
  d=/tmp/gauge_new/${s}_newtonian; mkdir -p "$d"
  printf 'gauge = newtonian\nroot = %s/out_\n' "$d" > /tmp/run.ini
  cat test/scenarios/${s}.ini >> /tmp/run.ini
  ./class /tmp/run.ini
  echo "== $s : new-newtonian vs baseline-newtonian =="
  python test/scenarios/compare_tol.py /tmp/gauge_base/${s}_newtonian "$d" "out_*.dat"
done
```
Expected: `gauge_lcdm`, `gauge_dcdm`, `gauge_scf` → all `OK` (behavior-identical). `gauge_ncdm`, `gauge_idmdr` → `OK` (the alpha accumulation is unchanged in this task; only the per-species transform moved). Any `FAIL` here means a transform override does not reproduce its old inline block — debug that species before proceeding.

- [ ] **Step 4: Confirm the fluid sign fix via cross-gauge agreement**

```bash
d=/tmp/gauge_new/gauge_fluid_newtonian; mkdir -p "$d"
printf 'gauge = newtonian\nroot = %s/out_\n' "$d" > /tmp/run.ini
cat test/scenarios/gauge_fluid.ini >> /tmp/run.ini
./class /tmp/run.ini
echo "== fluid: NEW newtonian vs synchronous (should now agree) =="
python test/scenarios/compare_tol.py /tmp/gauge_base/gauge_fluid_synchronous "$d" "out_*.dat"
```
Expected: `OK` (or markedly improved `worst_vs_colpeak` vs Task 1 Step 5). This demonstrates the corrected fluid sign restores gauge invariance.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "refactor(perturbations): dispatch Newtonian-gauge transform per species; fix fluid sign"
```

---

### Task 8: Replace the `alpha` accumulation with the downcast-free reduction

**Files:**
- Modify: `source/perturbations_module.cpp:4657-4722`

- [ ] **Step 1: Verify which species are dark energy** (so the reduction skips exactly the ones the old `delta_tot` omitted — fluid, scalar field, Λ).

Run: `grep -rn "EnergyType::" species/fluid.cpp species/scalar_field.h species/lambda.cpp species/*.cpp | grep -iE "BaseSpecies\(|CompositeSpecies\("`
Result (verified during execution): Fluid and Lambda are `EnergyType::DarkEnergy` (skipped). **ScalarField is `EnergyType::Other`, NOT DarkEnergy** — and per the user's decision it is *included* in the reduction (all species participate; the scalar field simply contributes 0 at IC because its synchronous field perturbation is δφ=δφ′=0). The composites (DCDM_DR/DNCDM_DR/IDM_DR_IDR/IDM_DRMD_IDR_DRMD) and NCDM are `Other` and are included (they aggregate their children via `Rho`/`Delta`). So the reduction skips only `DarkEnergy`; do NOT add ScalarField to a skip set. See the execution-deviations section for the `scalar_ctx` setup this requires.

- [ ] **Step 2: Replace the accumulation.** In `source/perturbations_module.cpp`, replace the block from line 4657 (`if (all_species_.count("CDM")) {` — the `delta_cdm` gathering) through line 4722 (`alpha = (...) / a_prime_over_a;`) with:

```cpp
      /* delta_tot = δρ/ρ_c and velocity_tot = (ρ+p)θ/ρ_c with ρ_c = rho_r + rho_m,
         summed over every non-dark-energy species reading the just-set synchronous
         y[]. Downcast-free analogue of perturb_total_stress_energy; includes all
         matter species (the old fraccdm channel counted only literal "CDM"). */
      double delta_rho_ic = 0., rho_plus_p_theta_ic = 0.;
      {
        size_t i = 0;
        for (const auto& sp : all_species_) {
          if (sp->energy_type() == BaseSpecies::EnergyType::DarkEnergy) {
            ++i;
            continue;
          }
          const double rho        = sp->Rho(ppw->pvecback);
          const double rho_plus_p = rho + sp->P(ppw->pvecback);
          const auto& layout      = *ppw->pv->species_layouts[i];
          delta_rho_ic +=
              rho * sp->Delta(layout, ppw->pv, ppw->pv->y, ppw->pvecback, ppw);
          rho_plus_p_theta_ic +=
              rho_plus_p * sp->Theta(layout, ppw->pv, ppw->pv->y, ppw->pvecback, ppw);
          ++i;
        }
      }
      double delta_tot    = delta_rho_ic / (rho_r + rho_m);
      double velocity_tot = rho_plus_p_theta_ic / (rho_r + rho_m);

      alpha = (eta + 3. / 2. * a_prime_over_a * a_prime_over_a / k / k / s2_squared *
                         (delta_tot + 3. * a_prime_over_a / k / k * velocity_tot)) /
              a_prime_over_a;
```

- [ ] **Step 3: Remove now-unused locals.** Deleting the accumulation drops the last readers of `delta_cdm`, `fracg`, `fracnu`, `fracb`, `fraccdm`, `fracidm_drmd`, `g_lay_ic`, `b_lay_ic`. `fracnu`/`fracg`/`fracb`/`fraccdm` are also copied into `ic_ctx` (lines 4399–4403) and read by species `ApplyInitialConditions` (e.g. `cdm.cpp:145`), so **keep their computation** (lines 4355–4374) and the `ic_ctx` assignments. Only delete the `double delta_cdm` declaration at line 4298 if the compiler reports it unused after this change. Build will tell you (Step 4); remove whatever it flags as unused.

- [ ] **Step 4: Build**

Run: `make class -j4`
Expected: compiles clean. Remove any variable the compiler reports as unused-and-set (`delta_cdm`, possibly `g_lay_ic`/`b_lay_ic` if they were only used here).

- [ ] **Step 5: Regenerate and compare — cross-gauge invariance is the arbiter**

```bash
for s in gauge_lcdm gauge_ncdm gauge_dcdm gauge_idmdr gauge_fluid gauge_scf; do
  d=/tmp/gauge_new2/${s}_newtonian; mkdir -p "$d"
  printf 'gauge = newtonian\nroot = %s/out_\n' "$d" > /tmp/run.ini
  cat test/scenarios/${s}.ini >> /tmp/run.ini
  ./class /tmp/run.ini
  echo "== $s : newtonian (post-alpha-redesign) vs SYNCHRONOUS baseline =="
  python test/scenarios/compare_tol.py /tmp/gauge_base/${s}_synchronous "$d" "out_*.dat"
done
```
Expected: all `OK`. Cross-gauge agreement (Newtonian vs synchronous) is the physical invariant; the `fraccdm`/IDR/NCDM changes are improvements, so agreement should hold or improve. If `gauge_idmdr` or `gauge_ncdm` shows a small elevated `worst_vs_colpeak` (≤ 1e-3) that is the documented IDR/NCDM-in-alpha change — acceptable within the 0.1% tolerance. Anything larger: investigate.

- [ ] **Step 6: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "refactor(perturbations): downcast-free alpha reduction; fix fraccdm cold-matter omission"
```

---

### Task 9: Full verification + cleanup

**Files:**
- Modify: `source/perturbations_module.cpp` (dead-code removal only)

- [ ] **Step 1: Remove residual dead locals.** Confirm `delta_ur`, `theta_ur`, `shear_ur`, `l3_ur`, `delta_dr`, `delta_cdm` declared at perturbations_module.cpp:4298–4300 are no longer read anywhere in the function (they were consumed by the old inline blocks and the `dispatch_species_ic` lambda still assigns the locals from `ic_ctx`). The `dispatch_species_ic` lambda assigns `delta_ur = ic_ctx.delta_ur; ...` (lines 4420–4425) — if nothing reads those locals now, delete those six assignment lines and the six declarations. Keep `eta` (still used at line 4628 synchronous + 4724 metric phi) and `alpha`/`alpha_prime`.

Run: `grep -n "delta_ur\|theta_ur\|shear_ur\|l3_ur\|delta_dr\|delta_cdm" source/perturbations_module.cpp` and confirm no remaining reads in the IC function before deleting.

- [ ] **Step 2: Build clean**

Run: `make class -j4 2>&1 | grep -i warning`
Expected: no new warnings about the touched code.

- [ ] **Step 3: Build the python wrapper**

Run: `make classy`
Expected: builds `classy` without error.

- [ ] **Step 4: Run the targeted reference tests** (compare against `classyref`; require it installed)

Run:
```bash
cd python && python -m pytest test_class.py -k "test_dcdm_dr_matches_reference or test_idm_dr_idr_perturbations_match_reference or test_rs_drag_matches_reference or test_theta_s_shooting_matches_reference or test_idr_without_idm_dr_computes or test_drmd_without_idr_drmd_computes" -v
```
Expected: all PASS (these run in default/synchronous gauge — unaffected by the Newtonian changes — and confirm no synchronous-path regression).

- [ ] **Step 5: Run the cross-gauge comparison suite**

Run:
```bash
cd python && COMPARE_OUTPUT_GAUGE=1 TEST_LEVEL=1 python -m pytest test_class.py -k "test_scenario" -q
```
Expected: PASS within `COMPARE_CL_RELATIVE_ERROR_GAUGE` (1.5e-2) and `COMPARE_PK_RELATIVE_ERROR_GAUGE` (5e-2). This exercises the parametrized grid in both gauges. Any `fail/` PDFs indicate a gauge mismatch to investigate.

- [ ] **Step 6: Run the scenario .ini grid (both gauges, tight tolerance)**

Run:
```bash
cd /Users/au192734/Projects/class_claude
for s in $(ls test/scenarios/*.ini | xargs -n1 basename | sed 's/.ini//'); do
  for g in synchronous newtonian; do
    d=/tmp/grid/${s}_${g}; mkdir -p "$d"
    printf 'gauge = %s\nroot = %s/out_\n' "$g" "$d" > /tmp/run.ini
    cat test/scenarios/${s}.ini >> /tmp/run.ini
    ./class /tmp/run.ini 2>/dev/null || echo "SKIP $s ($g) — scenario not valid in this gauge"
  done
done
echo "grid complete"
```
Expected: all scenarios run (DRMD scenarios will `SKIP`/error in Newtonian by design — the `has_idr_drmd` `class_test`). No crashes.

- [ ] **Step 7: Final commit**

```bash
git add source/perturbations_module.cpp
git commit -m "refactor(perturbations): remove dead gauge-transform locals"
```

---

## Self-Review

**Spec coverage:**
- New `RhoDotOverRho` accessor (default + DCDM) → Task 2. ✓
- New `PerturbSynchronousToNewtonian` hook + helper → Task 3. ✓
- Universal species (photons/baryons/cdm/fluid/UR/idm_dr/idr/idm_drmd/idr_drmd/dcdm) → Tasks 4, 6. ✓
- ScalarField + NCDM/DNCDM overrides → Task 5. ✓
- DarkRadiation re-seed + composites (DCDM_DR, DNCDM_DR, IDM_DR_IDR, IDM_DRMD) → Task 6. ✓
- Op2+op3 dispatch in module + block (e) removal → Task 7. ✓
- Op1 downcast-free reduction + `fraccdm` fix → Task 8. ✓
- Metric `phi` stays module-owned → kept in Task 7 (line 4724 untouched). ✓
- Behavior changes 1–4 → fluid sign (Task 4/7 Step 4), fraccdm (Task 8), IDR/NCDM-in-alpha (Task 8 Step 5), DNCDM/combined DR re-seed (Task 6 Step 4). ✓
- Testing: cross-gauge + reference + scenario grid → Tasks 1, 7, 8, 9. ✓
- `class_test(has_idr_drmd)` guard stays → not touched (Task 7 replaces only lines ≥ 4657; the guard at 4634 is left in place). ✓

**Type consistency:** `PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout&, double*, const PerturbIcContext&)` and `RhoDotOverRho(const double*, double)` and `ApplyFluidLikeNewtonianShift(double*, int, int, const double*, const PerturbIcContext&)` and `DarkRadiationSpecies::PerturbNewtonianReseed(const PerturbLayout&, double*, const PerturbIcContext&, double)` and `NCDMBaseSpecies::Dlnf0Dlnq(int, const double*)` are used identically across all tasks. Member names (`dcdm_`, `dr_sp_`, `pba_`, `dncdm_`, `idm_dr_`, `idr_`, `idm_drmd_`, `idr_drmd_`) match the headers. ✓

**Placeholder scan:** No TBD/TODO; every code step shows full code; every test step shows the command and expected result. The two "confirm the accessor name" notes (Task 5 Step 2 `GetBackgroundModule`, Task 8 Step 1 energy-type) are explicit verification steps with the fallback action stated, not placeholders. ✓

---

## Execution deviations (recorded post-implementation)

Five intentional deviations from the plan as written, all verified:

1. **`Dlnf0Dlnq` reused, not added (Task 5).** The codebase already had a pure-virtual `NCDMBaseSpecies::GetDlnf0DlnqForTensor(iq, pvecback)` (NCDM→static, DNCDM→time-dependent), semantically identical to the planned new `Dlnf0Dlnq`. Renamed it `GetDlnf0Dlnq` (mode-agnostic) and reused it in the NCDM per-q transform — no duplicate virtual, no separate DNCDM edit. The NCDM transform lives on **`NCDMBaseSpecies` (`ncdm_base_species.{h,cpp}`)**, not `ncdm_species.*`.

2. **IDM/IDR member impls live in `species/interacting_species.cpp`** (Task 6), not separate `idm_dr.cpp`/`idr.cpp`/`idm_drmd.cpp`/`idr_drmd.cpp` — all four member overrides went there.

3. **DCDM needed its own override (commit 304b9bde).** The plan's task steps added `cdm` (Task 4) but not `dcdm`; the DCDM_DR composite delegates to `dcdm_->PerturbSynchronousToNewtonian`, which would have hit the base no-op. Added `DCDMSpecies::PerturbSynchronousToNewtonian` (universal helper; its decay-aware `RhoDotOverRho` reproduces `(-3ℋ-aΓ)α`).

4. **ScalarField is `EnergyType::Other`, included with a synchronous `scalar_ctx` (Task 8, user decision).** All species participate in alpha. ScalarField contributes 0 at IC (δφ=δφ′=0) but its `Delta()` selects branch from `ppw->scalar_ctx.gauge` and the Newtonian branch reads workspace state invalid at IC. Fix: before the reduction, force `ppw->scalar_ctx.{gauge=synchronous, a=a, a2=a², k=k, k2=k²}` (the kinematic fields every species' Delta/Theta reads at IC). A first attempt set only `a2`/`k2`/`gauge` and crashed `gauge_ncdm` (NCDM::Delta/Theta read the non-squared `scalar_ctx.a`/`k` → `factor=(a_today/0)⁴`→NaN→singular matrix); fixed by also setting `a` and `k`.

5. **Three PRE-EXISTING cross-gauge failures surfaced, out of scope.** Newtonian-vs-synchronous at 0.1% tol: NCDM ~0.7% (FA evolution, not IC), DCDM fast-decay (Γ=10) Pk ~67%, scalar-field ~100%. All present on `master`, unchanged by this refactor (behavior-delta vs old inline is sub-ppm for every scenario). Candidates for separate follow-up investigation.
