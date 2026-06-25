# Fluid / PPF Class Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract `PpfFluid : public FluidSpecies` so the PPF closure stops being a `use_ppf_` bool threaded through `FluidSpecies`, and flip module-side PPF detection from `use_ppf()` to `dynamic_cast<PpfFluid*>`.

**Architecture:** Pure relocation. `FluidSpecies` keeps all shared background + the true-fluid δ/θ perturbations; `PpfFluid` overrides only the perturbation methods that differ and adds `ComputePpf`. No arithmetic or summation order changes, so output is expected **bit-identical**. Done in build-green steps: extend the safety net → scaffold the empty subclass → move PPF behavior into overrides (temporary duplication) → cut over construction + module detection → delete the now-dead code from `FluidSpecies`.

**Tech Stack:** C++17, CMake (`build/cmake`), Python `classy` extension (scikit-build-core), pytest characterization goldens (`.npz`).

**Spec:** `docs/superpowers/specs/2026-06-25-fluid-ppf-class-split-design.md`

## Global Constraints

- **C++-only headers** — plain C++, no `extern "C"`, no `#ifdef __cplusplus` ([[feedback_cpp_only_no_c_guards]]).
- **Never `git add -A`** — stage explicit paths; in-source CMake/Xcode artifacts must not be swept in ([[feedback_never_git_add_all]]).
- **Never hand-edit `python/cclassy.pxd`** — it is regenerated from headers at build by `generate_wrapper.py`. (No Python code references `use_ppf()`/`c_gamma_over_c_fld()`/`ComputePpf`, so removing them is safe.)
- **Verification bar = bit-identical.** This is a pure relocation; any output difference is a red flag to investigate, not to rebaseline ([[feedback_no_bit_identical_requirement]] still applies for the *lens*: use ≤0.1% + `Cl^TE` zero-crossing handling, never blind max-rel-diff).
- **No species-type picking in modules** — the new module-side `dynamic_cast<PpfFluid*>` is a one-line capability/type query (≤1 fluid), not a behavioural downcast ([[feedback_no_species_picking_in_modules]]).
- **Compile gate:** `make class`. **Behavioral gate:** `make classy-pip-dev` then `pytest`. **Suite:** `make test`.

> **TDD note for a relocation:** there is no new behavior to test-drive. The discipline here is *characterization*: Task 1 locks current output into goldens, and every later task must keep them **green**. "Red" only ever means a regression.

---

### Task 1: Extend the characterization safety net (non-PPF + EDE)

The committed goldens cover lcdm + the PPF `"fluid"` case but not a non-PPF fluid or EDE — the two `FluidSpecies` paths the split most affects. Lock them in **before** touching code, generated from the current (pre-refactor) tree.

**Files:**
- Modify: `python/gen_transfer_golden.py` (the `CASES` dict, ~line 19-27)
- Modify: `python/gen_background_golden.py` (the `CASES` dict, ~line 14)
- Create (generated): `python/transfer_golden/fluid_nonppf.npz`, `python/transfer_golden/fluid_ede.npz`, `python/background_golden/fluid_ede.npz`

**Interfaces:**
- Produces: two new transfer cases `fluid_nonppf`, `fluid_ede` and one new background case `fluid_ede`, consumed automatically by the parametrized `test_transfer_columns.py` / `test_background_columns.py`.

- [ ] **Step 1: Add the cases to both generators**

In `python/gen_transfer_golden.py`, inside `CASES`, after the existing `"fluid"` entry:

```python
    "fluid_nonppf": {
        "Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.1, "use_ppf": "no",
        "N_ur": 3.044,
    },
    "fluid_ede": {
        "Omega_Lambda": 0.0, "fluid_equation_of_state": "EDE",
        "w0_fld": -0.9, "Omega_EDE": 0.01, "N_ur": 3.044,
    },
```

In `python/gen_background_golden.py`, inside `CASES`, after the existing `"fluid"` entry:

```python
    "fluid_ede": {
        "Omega_Lambda": 0.0, "fluid_equation_of_state": "EDE",
        "w0_fld": -0.9, "Omega_EDE": 0.01, "N_ur": 3.044,
    },
```

- [ ] **Step 2: Build classy from the current tree**

Run: `make classy-pip-dev`
Expected: build succeeds, `classy` importable.

- [ ] **Step 3: Generate the new goldens (pre-refactor physics)**

Run: `python python/gen_transfer_golden.py && python python/gen_background_golden.py`
Expected: writes `python/transfer_golden/fluid_nonppf.npz`, `python/transfer_golden/fluid_ede.npz`, `python/background_golden/fluid_ede.npz` (plus rewrites existing ones identically).

- [ ] **Step 4: Confirm the full column suite passes (incl. new cases)**

Run: `pytest python/test_transfer_columns.py python/test_background_columns.py -v`
Expected: PASS, including `fluid_nonppf`, `fluid_ede`.

- [ ] **Step 5: Commit (stage explicit paths only)**

```bash
git add python/gen_transfer_golden.py python/gen_background_golden.py \
        python/transfer_golden/fluid_nonppf.npz python/transfer_golden/fluid_ede.npz \
        python/background_golden/fluid_ede.npz
git commit -m "test: characterize non-PPF and EDE fluid before PPF class split

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DSiFdA1yNDHQB67g3CX6HJ"
```

---

### Task 2: Scaffold `PpfFluid` (empty subclass) + open up `FluidSpecies`

Create the new class as a do-nothing subclass that compiles and is wired into the build but is not yet constructed anywhere. Promote the `FluidSpecies` members the PPF code will need.

**Files:**
- Create: `species/ppf_fluid.h`, `species/ppf_fluid.cpp`
- Modify: `species/fluid.h` (move 5 members `private` → `protected`)
- Modify: `CMakeLists.txt:83` (add `species/ppf_fluid.cpp`)

**Interfaces:**
- Produces: `class PpfFluid : public FluidSpecies` with a constructor `PpfFluid(const background& pba, double omega0_fld, equation_of_state fluid_eos, double w0_fld, double wa_fld, double cs2_fld, double Omega_EDE, double c_gamma_over_c_fld)`. Consumed by Tasks 4-5.

- [ ] **Step 1: Promote the members PPF needs from `private` to `protected`**

In `species/fluid.h`, move these five out of the `private:` block (line ~177-188) into a new `protected:` section (keep the rest private):

```cpp
 protected:
  const BackgroundModule* bgm_ = nullptr;
  int index_bg_rho_fld_        = -1;
  int index_bg_w_fld_          = -1;
  int index_bg_dw_over_da_fld_ = -1;
  double cs2_fld_              = 1.;
```

(`pba_`, `Omega0_fld_`, `fluid_eos_`, `w0_fld_`, `wa_fld_`, `Omega_EDE_`, the `index_bi_/index_tp_` slots, and — for now — `use_ppf_`/`c_gamma_over_c_fld_`/`idx_Gamma` stay as they are.)

- [ ] **Step 2: Create `species/ppf_fluid.h`**

```cpp
#pragma once

#include "fluid.h"

/**
 * PPF (Parametrised Post-Friedmann) dark-energy fluid. Same background as
 * FluidSpecies; the perturbations use the Gamma closure (Hu 0808.3125) instead
 * of delta/theta. At most one PPF fluid may exist (PPF is defined relative to
 * the whole universe).
 */
class PpfFluid : public FluidSpecies {
 public:
  PpfFluid(const background& pba,
           double omega0_fld,
           equation_of_state fluid_eos,
           double w0_fld,
           double wa_fld,
           double cs2_fld,
           double Omega_EDE,
           double c_gamma_over_c_fld);
};
```

- [ ] **Step 3: Create `species/ppf_fluid.cpp`**

```cpp
#include "ppf_fluid.h"

PpfFluid::PpfFluid(const background& pba,
                   double omega0_fld,
                   equation_of_state fluid_eos,
                   double w0_fld,
                   double wa_fld,
                   double cs2_fld,
                   double Omega_EDE,
                   double c_gamma_over_c_fld)
    : FluidSpecies(pba, omega0_fld, fluid_eos, w0_fld, wa_fld, cs2_fld, Omega_EDE,
                   /*use_ppf=*/true, c_gamma_over_c_fld) {}
```

(For now it forwards `use_ppf=true` to the still-existing base ctor param; Task 5 removes that param.)

- [ ] **Step 4: Wire into CMake**

In `CMakeLists.txt`, add after line 83 (`species/fluid.cpp`):

```cmake
  species/ppf_fluid.cpp
```

- [ ] **Step 5: Build**

Run: `make class`
Expected: compiles cleanly; `PpfFluid` exists but is unconstructed.

- [ ] **Step 6: Commit**

```bash
git add species/ppf_fluid.h species/ppf_fluid.cpp species/fluid.h CMakeLists.txt
git commit -m "fluid: scaffold empty PpfFluid subclass + protected members

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DSiFdA1yNDHQB67g3CX6HJ"
```

---

### Task 3: Move PPF perturbation behavior into `PpfFluid` overrides

Add the overrides + `ComputePpf` to `PpfFluid`, copied from the `if (use_ppf_)` true-branches of `FluidSpecies`. `FluidSpecies` is left untouched (temporary duplication — resolved in Task 5). `PpfFluid` reuses the inherited `FluidSpecies::PerturbLayout` (which still owns `idx_Gamma`), so no layout clash yet.

**Files:**
- Modify: `species/ppf_fluid.h`, `species/ppf_fluid.cpp`

**Interfaces:**
- Consumes: `FluidSpecies::ComputeWFld`, `Rho`, `W`, the protected `cs2_fld_`/`bgm_`/`index_bg_*`, and the inherited `FluidSpecies::PerturbLayout` (with `idx_Gamma`).
- Produces: `PpfFluid` overrides of `RegisterPerturbationIndices`, `PerturbDerivs`, `FillSources`, `PrintVariables`, `CopyPerturbationsAcrossSwitch`, plus `void ComputePpf(double k, double a, double a_prime_over_a, const precision* ppr, const double* y, perturb_workspace* ppw) const`.

- [ ] **Step 1: Declare the overrides in `species/ppf_fluid.h`**

Add inside `class PpfFluid` (after the constructor):

```cpp
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  void ComputePpf(double k,
                  double a,
                  double a_prime_over_a,
                  const precision* ppr,
                  const double* y,
                  perturb_workspace* ppw) const;
```

- [ ] **Step 2: Implement the overrides in `species/ppf_fluid.cpp`**

These are the PPF branches lifted verbatim from `FluidSpecies`. `PerturbLayout` here resolves to the inherited `FluidSpecies::PerturbLayout`.

```cpp
void PpfFluid::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                           perturb_vector* /*pv*/,
                                           const precision* /*ppr*/,
                                           int& index_pt,
                                           const perturb_workspace* /*ppw*/,
                                           int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  class_define_index(layout.idx_Gamma, _TRUE_, index_pt, 1);
}

void PpfFluid::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                             double /*tau*/,
                             const double* /*y*/,
                             double* dy,
                             const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  dy[layout.idx_Gamma] = ppaw.ppw->Gamma_prime_fld;
}

void PpfFluid::FillSources(const BaseSpecies::PerturbLayout& /*base*/,
                           const double* /*y*/,
                           const double* /*dy*/,
                           PerturbSourceContext& ctx) const {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const double* pvecback     = ppw->pvecback.data();
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  if (index_tp_delta_ >= 0) {
    const double w_fld     = W(pvecback);
    const double delta_fld = ppw->delta_rho_fld / Rho(pvecback);
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp_delta_, ctx.index_tau,
                          ctx.index_k,
                          delta_fld + 3. * ctx.a_prime_over_a * (1. + w_fld) * ctx.theta_over_k2);
  }
  if (index_tp_theta_ >= 0) {
    const double w_fld     = W(pvecback);
    const double theta_fld = ppw->rho_plus_p_theta_fld / (1. + w_fld) / Rho(pvecback);
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp_theta_, ctx.index_tau,
                          ctx.index_k, theta_fld + ctx.theta_shift);
  }
}

void PpfFluid::PrintVariables(PerturbColumnWriter& w,
                              double /*tau*/,
                              const double* /*y*/,
                              const PerturbationsModule& /*mod*/,
                              const perturb_workspace* ppw) const {
  double delta_rho_fld = 0., rho_plus_p_theta_fld = 0., delta_p_fld = 0.;
  if (!w.IsTitleMode()) {
    delta_rho_fld        = ppw->delta_rho_fld;
    rho_plus_p_theta_fld = ppw->rho_plus_p_theta_fld;
    delta_p_fld          = ppw->delta_p_fld;
  }
  w.Add("delta_rho_fld", delta_rho_fld, true);
  w.Add("rho_plus_p_theta_fld", rho_plus_p_theta_fld, true);
  w.Add("delta_p_fld", delta_p_fld, true);
}

void PpfFluid::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                             const BaseSpecies::PerturbLayout& new_base,
                                             const double* old_y,
                                             double* new_y,
                                             const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_Gamma >= 0 && new_l.idx_Gamma >= 0)
    new_y[new_l.idx_Gamma] = old_y[old_l.idx_Gamma];
}
```

- [ ] **Step 3: Move `ComputePpf` into `PpfFluid`**

Copy the body of `FluidSpecies::ComputePpf` verbatim from `species/fluid.cpp:422-537` into `species/ppf_fluid.cpp`, changing only the method qualifier `FluidSpecies::ComputePpf` → `PpfFluid::ComputePpf`. The `static_cast<const PerturbLayout&>` at the old line 449 now resolves to the inherited `FluidSpecies::PerturbLayout` (still valid). Leave the `FluidSpecies` copy in place for now.

Add the needed includes at the top of `species/ppf_fluid.cpp`:

```cpp
#include "ppf_fluid.h"

#include "background_module.h"
#include "perturbations_module.h"
```

- [ ] **Step 4: Build**

Run: `make class`
Expected: compiles. `PpfFluid` now has full PPF behavior but is still unconstructed (no behavioral change yet).

- [ ] **Step 5: Commit**

```bash
git add species/ppf_fluid.h species/ppf_fluid.cpp
git commit -m "fluid: implement PPF behavior on PpfFluid (overrides + ComputePpf)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DSiFdA1yNDHQB67g3CX6HJ"
```

---

### Task 4: Cut over — construct `PpfFluid` + module detection by type

Flip the factory to mint `PpfFluid` for `use_ppf`, and the module to detect PPF via `dynamic_cast`. After this task PPF physics runs through `PpfFluid`'s overrides. **This is the critical behavioral gate** — the goldens must stay green.

**Files:**
- Modify: `species/fluid.cpp:628-637` (`CreateAll` construction branch)
- Modify: `species/fluid.cpp` top (add `#include "ppf_fluid.h"`)
- Modify: `source/perturbations_module.h:27,40-41` (pointer + accessor type)
- Modify: `source/perturbations_module.cpp:93-97` (resolution) + add `#include "ppf_fluid.h"`

**Interfaces:**
- Consumes: `PpfFluid` ctor (Task 2), `PpfFluid` overrides (Task 3).
- Produces: a `PpfFluid*` reachable via `PerturbationsModule::ppf_fluid()`.

- [ ] **Step 1: Branch the factory construction**

In `species/fluid.cpp`, add `#include "ppf_fluid.h"` near the top, then replace the `result.push_back({...})` at lines 628-637 with:

```cpp
  std::unique_ptr<BaseSpecies> sp;
  if (use_ppf) {
    sp = std::make_unique<PpfFluid>(*ctx.pba, omega0_fld, fluid_eos, w0_fld, wa_fld,
                                    cs2_fld, Omega_EDE, c_gamma_over_c_fld);
  }
  else {
    sp = std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld, fluid_eos, w0_fld, wa_fld,
                                        cs2_fld, Omega_EDE, /*use_ppf=*/false,
                                        c_gamma_over_c_fld);
  }
  result.push_back({"Fluid", std::move(sp)});
```

(The `FluidSpecies` ctor still takes `use_ppf`/`c_gamma_over_c_fld` at this point; Task 5 trims them.)

- [ ] **Step 2: Change the module pointer + accessor type**

In `source/perturbations_module.h`, change line 27 and the accessor at 40-41 from `FluidSpecies*` to `PpfFluid*`:

```cpp
  PpfFluid* ppf_fluid = nullptr;   // the single PPF fluid, if present
```
```cpp
  PpfFluid* ppf_fluid() const {
    return resolved_.ppf_fluid;
  }
```

Add `#include "ppf_fluid.h"` to `source/perturbations_module.h` (or the `.cpp` if the header only needs a forward declaration — a forward `class PpfFluid;` in the header plus the include in the `.cpp` is cleaner).

- [ ] **Step 3: Replace the resolution with a type query**

In `source/perturbations_module.cpp`, add `#include "ppf_fluid.h"`, then replace lines 93-97:

```cpp
  if (auto* p = all_species_.find("Fluid"))
    if (auto* f = dynamic_cast<PpfFluid*>(p->get()))
      resolved_.ppf_fluid = f;
```

- [ ] **Step 4: Build**

Run: `make class`
Expected: compiles. PPF now dispatches through `PpfFluid`; the `FluidSpecies` `use_ppf_` branches are dead for constructed objects.

- [ ] **Step 5: Behavioral gate — goldens must stay green**

Run: `make classy-pip-dev && pytest python/test_transfer_columns.py python/test_background_columns.py -v`
Expected: PASS for all cases — crucially `fluid` (PPF, now via `PpfFluid`), `fluid_nonppf`, `fluid_ede`, `lcdm`. Background bit-identical; transfer within the goldens' tolerances. Any failure here means the move changed behavior — STOP and diff (per superpowers:systematic-debugging), do not rebaseline.

- [ ] **Step 6: Commit**

```bash
git add species/fluid.cpp source/perturbations_module.h source/perturbations_module.cpp
git commit -m "fluid: route PPF through PpfFluid; detect it by dynamic_cast

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DSiFdA1yNDHQB67g3CX6HJ"
```

---

### Task 5: Remove the dead PPF code from `FluidSpecies` + finalize shape

Delete what's now unreachable in `FluidSpecies`, move `idx_Gamma` ownership to `PpfFluid`, and trim the ctor. After this `FluidSpecies` is the pure true fluid.

**Files:**
- Modify: `species/fluid.h`, `species/fluid.cpp`
- Modify: `species/ppf_fluid.h` (own `PerturbLayout` + `CreatePerturbLayout`)
- Modify: `species/fluid.cpp:628-636` (`FluidSpecies` ctor call drops 2 args)
- Modify: `source/background.h:77`, `source/input_module.cpp:831-832` (comments)

**Interfaces:**
- Produces: `FluidSpecies` ctor `FluidSpecies(const background&, double omega0_fld, equation_of_state, double w0_fld, double wa_fld, double cs2_fld, double Omega_EDE)` (no `use_ppf`/`c_gamma_over_c_fld`). `PpfFluid::PerturbLayout : FluidSpecies::PerturbLayout { int idx_Gamma; }`.

- [ ] **Step 1: Slim `FluidSpecies` declarations (`species/fluid.h`)**

- Remove `idx_Gamma` from `FluidSpecies::PerturbLayout` (line 25).
- Remove the `bool use_ppf` and `double c_gamma_over_c_fld` ctor params (lines 39-40).
- Remove accessors `use_ppf()` (61-63) and `c_gamma_over_c_fld()` (64-66).
- Remove the `ComputePpf` declaration (102-113).
- Remove members `use_ppf_` (187) and `c_gamma_over_c_fld_` (188).

- [ ] **Step 2: Slim `FluidSpecies` definitions (`species/fluid.cpp`)**

- ctor (14-26): drop the `use_ppf`/`c_gamma_over_c_fld` params and their initialisers.
- `RegisterPerturbationIndices` (116-122): delete the `if(!use_ppf_)`/`else` — keep only:
  ```cpp
    auto& layout = static_cast<PerturbLayout&>(base);
    class_define_index(layout.idx_delta, _TRUE_, index_pt, 1);
    class_define_index(layout.idx_theta, _TRUE_, index_pt, 1);
  ```
- `PerturbDerivs` (140-157): delete the `if(!use_ppf_)` wrapper and the `else { dy[idx_Gamma]=… }` — keep only the δ/θ body.
- `FillSources` (177, 191): drop the ternaries — `const double delta_fld = y[layout.idx_delta];` and `const double theta_fld = y[layout.idx_theta];`.
- `ApplyInitialConditions` (208-209): delete `if (use_ppf_) return;` (the `idx_delta<0||idx_theta<0` guard at 210 stays).
- `PrintVariables` (245-259): delete the `if (use_ppf_)` branch, keep only the `else` (StressEnergy-derived) body.
- `CopyPerturbationsAcrossSwitch` (564-565): delete the `idx_Gamma` copy lines.
- Delete the whole `FluidSpecies::ComputePpf` definition (422-537).

- [ ] **Step 3: Give `PpfFluid` its own layout (`species/ppf_fluid.h`)**

Add inside `class PpfFluid` and provide the factory override:

```cpp
  struct PerturbLayout : FluidSpecies::PerturbLayout {
    int idx_Gamma = -1;
  };
  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }
```

Now `PpfFluid`'s overrides + `ComputePpf` cast to *this* `PerturbLayout` (they already say `static_cast<...PerturbLayout&>`, which now binds to `PpfFluid::PerturbLayout`). The inherited self-zeroing `StressEnergy`/`ApplyInitialConditions` cast to `FluidSpecies::PerturbLayout` — valid upcast, and `idx_delta`/`idx_theta` are `-1` for PPF.

- [ ] **Step 4: Trim the `PpfFluid` ctor forwarding (`species/ppf_fluid.cpp`)**

The base ctor no longer takes `use_ppf`. Update the forwarding:

```cpp
    : FluidSpecies(pba, omega0_fld, fluid_eos, w0_fld, wa_fld, cs2_fld, Omega_EDE),
      c_gamma_over_c_fld_(c_gamma_over_c_fld) {}
```

Add the member to `species/ppf_fluid.h` (`private:`):

```cpp
 private:
  double c_gamma_over_c_fld_ = 0.4;
```

`ComputePpf` reads `c_gamma_over_c_fld_` — now its own member.

- [ ] **Step 5: Fix the factory `FluidSpecies` branch (`species/fluid.cpp`)**

In the `else` branch from Task 4, drop the two trailing args:

```cpp
    sp = std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld, fluid_eos, w0_fld, wa_fld,
                                        cs2_fld, Omega_EDE);
```

- [ ] **Step 6: Update stale comments**

- `source/background.h:77`: change the comment so `use_ppf`/`c_gamma_over_c_fld` are described as living on `PpfFluid`, the rest on `FluidSpecies`.
- `source/input_module.cpp:831-832`: same — note PPF-only params live on `PpfFluid`.

- [ ] **Step 7: Build**

Run: `make class`
Expected: compiles; no `use_ppf_` symbol remains in `FluidSpecies`.

- [ ] **Step 8: Full verification**

```bash
make classy-pip-dev
pytest python/test_transfer_columns.py python/test_background_columns.py -v
make test
```
Expected: all PASS, every fluid case bit-identical/within-tol vs Task 1 goldens; ctest suite (`test-parser`, `test-bisection`, `test-photons`, dotsyntax/ncdm scenarios) green.

- [ ] **Step 9: Grep for leftovers**

Run: `grep -rn 'use_ppf' species/fluid.h species/fluid.cpp`
Expected: no matches (the flag is fully gone from `FluidSpecies`; remaining `use_ppf` references should only be the input *string key* in `CreateAll` and `PpfFluid` plumbing).

- [ ] **Step 10: Commit**

```bash
git add species/fluid.h species/fluid.cpp species/ppf_fluid.h species/ppf_fluid.cpp \
        source/background.h source/input_module.cpp
git commit -m "fluid: remove dead use_ppf branches; FluidSpecies is now the true fluid

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01DSiFdA1yNDHQB67g3CX6HJ"
```

---

## Self-Review

**Spec coverage:**
- §2 hierarchy (`PpfFluid : FluidSpecies`, EoS enum stays) → Tasks 2-5. ✓
- §3 split inventory (every `use_ppf_` site) → Task 5 Step 2 enumerates each by line. ✓
- §3 inherited safety net (`ApplyInitialConditions`/`StressEnergy`/`PerturbSynchronousToNewtonian` self-zero) → Task 5 Step 3 note; not overridden. ✓
- §3 member access (private→protected) → Task 2 Step 1. ✓
- §3 file layout (`species/ppf_fluid.{h,cpp}`) → Task 2. ✓
- §4 module `dynamic_cast<PpfFluid*>` → Task 4 Steps 2-3. ✓
- §5 factory branch → Task 4 Step 1 + Task 5 Step 5. ✓
- §6 comments → Task 5 Step 6. ✓
- §7 verification (bit-identical, scenarios both paths) → Task 1 goldens + Task 4 Step 5 + Task 5 Step 8. ✓
- §8 (multi-fluid) → explicitly out of scope; no task. ✓

**Placeholder scan:** no TBD/TODO; `ComputePpf` is a verbatim move with an exact source range (`fluid.cpp:422-537`) rather than re-pasted 115 lines. All other code steps show full code.

**Type consistency:** `PpfFluid` ctor signature identical in Task 2 (decl), Task 4 (call), Task 5 (forwarding). `ppf_fluid` typed `PpfFluid*` consistently in Task 4. `PerturbLayout`/`idx_Gamma` ownership transfers cleanly: inherited base layout in Task 3 → own subclass layout in Task 5, both valid for the same `static_cast` text.
