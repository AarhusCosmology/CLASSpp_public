# Fluid / PPF class split (Tier A)

Date: 2026-06-25
Scope decision: **Tier A — class split only.** Multi-fluid (N DE fluids) is
explicitly deferred; see §8. Builds on merged #308 (the module-side
`ppf_fluid()` accessor and the fluid-in-generic-loop refactor already exist).

## 1. Motivation & scope

Today `FluidSpecies` carries a construction-fixed `use_ppf_` bool and threads it
through the perturbation code as `if (use_ppf_) … else …` branches, plus a block
of methods that return 0 / no-op in PPF mode. The two schemes (true fluid:
δ/θ + cs²; PPF: Γ + modified-gravity closure) share the entire background but
diverge in perturbations — the textbook case for inheritance.

**Why only Tier A (single fluid, no N-fluid machinery).** A cost/benefit pass
established that multiple DE fluids are *background-indistinguishable* from one
combined fluid carrying `w_eff(a) = ΣPᵢ/Σρᵢ` (since `H(a)` depends only on
`ρ_DE(a)`). The sole irreducible gain is per-component sound speed `cs²` (distinct
clustering scales) and mixing closure schemes (PPF phantom-crossing + a clustering
true fluid) — both niche-research, weakly constrained. The class split, by
contrast, is worth doing on its own: it is the user's original motivation and
stands alone at single-fluid. So we do the split now and defer the genuine
multi-fluid tax (effective-w aggregate + per-instance closure + dot-naming) to
§8, exactly where the #308 spec left it ("when the right shape is knowable").

This is a **pure relocation**: no arithmetic and no summation order changes, so
output is expected **bit-identical** (see §7).

## 2. Class hierarchy

```
FluidSpecies            (concrete, the true fluid + ALL shared background)
  background: rho_fld ODE, w(a) via ComputeWFld, EoS enum, Rho/P/PPrime
  perturb:    delta/theta evolution (cs2_fld)
        |
        +-- PpfFluid : public FluidSpecies
              holds:     c_gamma_over_c_fld_
              layout:    PerturbLayout : FluidSpecies::PerturbLayout { idx_Gamma }
              overrides: CreatePerturbLayout, RegisterPerturbationIndices,
                         PerturbDerivs, FillSources, PrintVariables,
                         CopyPerturbationsAcrossSwitch
              adds:      ComputePpf
              inherits (deliberately, self-zeroing):
                         ApplyInitialConditions, StressEnergy,
                         PerturbSynchronousToNewtonian, RegisterTransferSourceIndices,
                         WriteOutputColumns, all background methods
```

The EoS parametrisation (`fluid_eos_`) stays a per-instance enum on
`FluidSpecies`. "Proper EDE / beyond-w(a)" physics is explicitly **not** a fluid —
it escapes to its own `BaseSpecies` species (forced by `BaseSpecies` pure
virtuals), so no abstract `FluidBase` is warranted for a hierarchy of two.

## 3. The split inventory (file:line, against current `species/fluid.{h,cpp}`)

### FluidSpecies loses its PPF burden — becomes the pure true fluid
- **ctor** (`fluid.cpp:14-26`, `fluid.h:32-40`): drop `use_ppf` and
  `c_gamma_over_c_fld` params; drop members `use_ppf_` (`fluid.h:187`),
  `c_gamma_over_c_fld_` (`fluid.h:188`) and accessors `use_ppf()` (`fluid.h:61`),
  `c_gamma_over_c_fld()` (`fluid.h:64`).
- **PerturbLayout** (`fluid.h:22-26`): keep `idx_delta`, `idx_theta`; **remove
  `idx_Gamma`** (moves to `PpfFluid::PerturbLayout`).
- **RegisterPerturbationIndices** (`fluid.cpp:116-122`): drop the `if(!use_ppf_)`
  — always register `idx_delta`/`idx_theta`.
- **PerturbDerivs** (`fluid.cpp:140-157`): drop the branch — keep only the
  δ/θ evolution.
- **FillSources** (`fluid.cpp:177,191`): drop the `use_ppf_ ? … : …` ternaries —
  `delta_fld = y[idx_delta]`, `theta_fld = y[idx_theta]`.
- **ApplyInitialConditions** (`fluid.cpp:208-209`): drop `if(use_ppf_) return;`.
  The existing guard `if (idx_delta<0||idx_theta<0) return;` (`fluid.cpp:210`)
  stays and makes this method a safe no-op when inherited by PPF.
- **PrintVariables** (`fluid.cpp:245-259`): drop the branch — keep only the
  `StressEnergy`-derived path.
- **CopyPerturbationsAcrossSwitch** (`fluid.cpp:560-565`): drop the `idx_Gamma`
  copy (moves to the PPF override); keep δ/θ copies.
- **ComputePpf** (`fluid.h:108`, `fluid.cpp:422-537`): **moves to PpfFluid.**

### Unchanged on FluidSpecies, inherited by PpfFluid (the safety net)
`StressEnergy` (`fluid.cpp:266-316`), `PerturbSynchronousToNewtonian`
(`fluid.cpp:541-549`), `RegisterTransferSourceIndices`, `WriteOutputColumns`
(`fluid.cpp:224-236`), and all background methods. Each is guarded on
`idx_delta>=0`/`idx_theta>=0` or operates only on background, so for a PPF fluid
(which never registers δ/θ) they **self-zero / no-op** with no override — this is
the load-bearing #308 finding.

### PpfFluid (new) — `species/ppf_fluid.{h,cpp}`
- `ctor(... same background params ..., double c_gamma_over_c_fld)` forwards the
  background params to the `FluidSpecies` base ctor and stores
  `c_gamma_over_c_fld_`.
- `struct PerturbLayout : FluidSpecies::PerturbLayout { int idx_Gamma = -1; };`
  + `CreatePerturbLayout()` override returning it. Inherited methods that cast to
  `FluidSpecies::PerturbLayout` (e.g. the self-zeroing `StressEnergy`) remain
  valid via upcast.
- Overrides: `RegisterPerturbationIndices` (register `idx_Gamma` only),
  `PerturbDerivs` (`dy[idx_Gamma] = ppw->Gamma_prime_fld`), `FillSources` (publish
  `ppw->delta_rho_fld`/`rho_plus_p_theta_fld`), `PrintVariables` (read the
  published `_fld` workspace fields), `CopyPerturbationsAcrossSwitch` (copy
  `idx_Gamma`).
- `ComputePpf` moved here verbatim; its casts to `PerturbLayout` now resolve to
  `PpfFluid::PerturbLayout`.

### Member access
`ComputePpf` and the overrides read FluidSpecies internals currently `private`.
Promote these from `private` → `protected`: `bgm_`, `index_bg_rho_fld_`,
`index_bg_w_fld_`, `index_bg_dw_over_da_fld_`, `cs2_fld_`, and
`index_tp_delta_`/`index_tp_theta_` (the latter two read by `PpfFluid::FillSources`).
Everything else stays `private`. (`Rho`/`P`/`W`/`ComputeWFld` are already public.)
`c_gamma_over_c_fld_` is *not* promoted — it leaves `FluidSpecies` entirely and
becomes a `private` member of `PpfFluid` (§3, PpfFluid).

## 4. Module side (`source/perturbations_module.{h,cpp}`)

- `perturbations_module.h:27,40-41`: change the stored pointer + accessor type
  from `FluidSpecies* ppf_fluid` to `PpfFluid* ppf_fluid`.
- `perturbations_module.cpp:93-97`: replace the `find("Fluid")` +
  `static_cast<FluidSpecies*>` + `f->use_ppf()` resolution with a **type query**:
  ```cpp
  if (auto* p = all_species_.find("Fluid"))
    if (auto* f = dynamic_cast<PpfFluid*>(p->get()))
      resolved_.ppf_fluid = f;
  ```
  The type now *is* the capability — `use_ppf()` is gone. Aligns with
  "no species-type picking in modules" (this is a capability/type query, not a
  behavioural downcast).
- Sites `:555` (`!ppf_fluid()`) and `:4911` (`if (auto* f = ppf_fluid())`) are
  **unchanged** — they already speak through the accessor.
- Add `#include "ppf_fluid.h"` where `PpfFluid` is named.

## 5. Input (`species/fluid.cpp::CreateAll`)

Parsing of `use_ppf` (`fluid.cpp:617-624`) and `c_gamma_over_c_fld`
(`fluid.cpp:625-626`) is unchanged. Only the construction branches on the flag:

```cpp
auto sp = use_ppf
  ? std::unique_ptr<BaseSpecies>(std::make_unique<PpfFluid>(*ctx.pba, omega0_fld,
        fluid_eos, w0_fld, wa_fld, cs2_fld, Omega_EDE, c_gamma_over_c_fld))
  : std::unique_ptr<BaseSpecies>(std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld,
        fluid_eos, w0_fld, wa_fld, cs2_fld, Omega_EDE));
result.push_back({"Fluid", std::move(sp)});
```

Both are keyed `"Fluid"` (single instance), so `kAllSpeciesFactories`, the
closure machinery, and `background_w_fld`'s `static_cast<FluidSpecies&>(at("Fluid"))`
(`background_module.cpp:442`) are all **untouched** — a `PpfFluid` *is* a
`FluidSpecies`, so the background cast and `ComputeWFld` call still resolve.

## 6. Comment hygiene
- `background.h:77` and `input_module.cpp:831-832`: update the comments that list
  `use_ppf`/`c_gamma_over_c_fld` as living on `FluidSpecies` to reflect that the
  PPF-only ones now live on `PpfFluid`.

## 7. Verification

Because this is a pure relocation (no arithmetic reordered, branches become
virtual dispatch over identical expressions), the **expectation is bit-identical**
output for every scenario — any difference is a red flag to investigate, not to
rebaseline away. Per repo norms ([[feedback_no_bit_identical_requirement]]),
still verify with the ≤0.1% lens and `Cl^TE` zero-crossing handling rather than
blind max-rel-diff.

Scenarios, both gauges (synchronous + newtonian), lensed TT + P(k):
1. **ΛCDM** (no fluid) — exercises the unchanged path.
2. **PPF CPL** (`use_ppf=yes`, `wa_fld≠0`) — exercises `PpfFluid`.
3. **non-PPF CPL** (`use_ppf=no`) — exercises the slimmed `FluidSpecies`.
4. ~~**EDE** (`fluid_equation_of_state=EDE`)~~ — *not testable*: the `EDE` branch of
   `ComputeWFld` is an unfinished `class_stop` stub (fluid.cpp), so no EDE golden can
   be generated. CLP is the only working EoS; the split is EoS-agnostic regardless.

Verification was performed as a HEAD-vs-master A/B (scale-relative, zero-crossing-safe
≤0.1%), **not** against the committed goldens — those are pre-existing stale vs the
ffast-math build (#338). Result: `fluid` (PPF→PpfFluid) 5.99e-08; `fluid_nonppf`
(true→FluidSpecies) 5.43e-04; true-fluid `dy[]` arithmetic character-identical to master.

Gate: `python/test_background_columns.py`, `python/test_transfer_columns.py`, and
`TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` against the existing `classyref`. If
bit-identical holds, **no golden/classyref regeneration is needed** (unlike #308).

## 8. Deferred: full multi-fluid (Tier C) — recorded so the design isn't lost

Not implemented now. When a genuine multi-`cs²` DE model needs it, the worked-out
design is:

- **Class:** unchanged — `FluidSpecies` / `PpfFluid` already mint cleanly as
  `unique_ptr<BaseSpecies>`; the generic stress-energy loop already iterates all
  species. `ComputePpf` already self-subtracts its own ρ+p, so it composes with N
  true fluids (they correctly remain in the rest-of-universe total).
- **Input:** named typed dot-notation, conforming to the existing NCDM mechanism
  (`pfc.instances_with("type", "fluid")` / `"ppf_fluid"`; cf.
  `input_module.cpp:594`, `ncdm_species.cpp:205`). The `.type` value *names the
  class*, so `use_ppf` survives only as a legacy flat-key translation. Output
  columns keyed by instance `name()` (like NCDM); legacy `"Fluid"` keeps `…_fld`.
  Enforce **≤1 PPF instance** in the resolver. ≥2 `DarkEnergy` species coexisting
  is then possible.
- **Effective w:** replace `background_w_fld` with a background-module service
  returning `{w_eff, dw_eff/da}` = `Σ Pᵢ/Σ ρᵢ` over `EnergyType::DarkEnergy`
  species (Fluid+Λ; `ScalarField` is `Other`, excluded — matches today). The
  `integral` out-param leaves this service and becomes per-instance (each
  `FluidSpecies::SetBackgroundInitialConditions`). Consumers: halofit
  (`nonlinear_module.cpp:1965,2887`), thermo (`thermodynamics_module.cpp:2910`),
  perturbation-IC checks (`perturbations_module.cpp:546-547`).
- **Closure (Option A, per-instance, omitted-Ω):** add `closure_key` to
  `SpeciesBuildContext`; a closure prepass picks the single instance omitting its
  Ω (validate exactly ≤1 across all closure-capable species). `ConstructSpecies`
  lowers from per-entry to per-instance: pass 1 builds `key != closure_key`
  (non-fluid factories unaffected), pass 2 builds only `closure_key` with the
  override; `CreateAll` reads phase from `omega0_closure_override.has_value()`.
- **Pk_equal** (deferred #308 sites): guard to single-fluid (error on multi).
