# Type-3 momentum-transfer coupling (CDM ↔ scalar field) — design

**Date:** 2026-06-26
**Reference physics:** Pourtsidou & Tram, *Reconciling CMB and structure growth measurements
with dark energy interactions*, arXiv:1604.04222v2 (the "Type-3", pure-momentum-transfer
coupled-quintessence model). LaTeX source in `arXiv-1604.04222v2/paperv2.tex`.
**Reference implementation:** the original (old-C CLASS) lives on
`github.com/ThomasTram/class` branch `CoupledQuintessence` (`source/perturbations.c`,
`source/background.c`). The equations there are the **validated** reference (they produced the
published figures) and the spec follows their *code-faithful* forms, which are algebraically
expanded versions of the paper's equations. In that code the coupling parameter β is called
`scf_veta` ("veta"); `factorveta` is a vestigial `1.0` debug hook.

## 1. Purpose and success criteria

Implement the Type-3 pure-momentum-transfer coupling between cold dark matter and a quintessence
scalar field in the CLASSpp **species system**. The primary goal is to **stress-test and improve
the species system** by driving a genuine bidirectional perturbation coupling between two
first-class species through the existing `CompositeSpecies` / `AddCouplingDerivs` machinery, and
to land the reusable infrastructure that any "something-coupled-to-CDM" model needs.

**Success = reproduce the qualitative features of the paper's Figure 2** (`cl_ratio_beta.pdf`,
`pk_ratio_beta.pdf` — the `T3_ph` ratio plots over a wide range of large negative β): the
large-scale ISW-only signature in `Cl^TT/Cl^TT(β=0)`, and the `P(k)` growth suppression that
turns up again below β ≈ −10². The strong-drag (large |β|) regime is therefore **in scope**.

**Non-goals / out of scope (first cut):**
- **Newtonian gauge coupling.** Postponed — it requires a fresh derivation not present in the
  paper or the original code. Synchronous gauge only. (See §7.)
- The MCMC / parameter-inference results of the paper.
- Energy-exchange (Pourtsidou-Tram Types 1 & 2, `Q ≠ 0`) — see constraint ③ in §2.

## 2. Species-system constraints surfaced (and their resolutions)

The user explicitly asked to flag where the current design makes a CDM + scalar-field composite
hard. Three were found:

- **① (Blocking) Synchronous-gauge CDM has no velocity variable.** `CDMSpecies` sets
  `idx_theta = -1` in synchronous gauge (`cdm.cpp:58-62`); CDM *defines* the gauge by being at
  rest. The θ=0 assumption is hard-coded in three places: `RegisterPerturbationIndices`
  (`cdm.cpp:58-62`), `PerturbDerivs` synchronous branch (`cdm.cpp:79-81`), and `PrintVariables`
  (`cdm.cpp:194-196`); the transfer-source path gates θ on `gauge != synchronous`
  (`cdm.cpp:38-42`). **Resolution:** an opt-in construction flag on `CDMSpecies` (see §4.2),
  used directly — *not* a subclass. A subclass would have to override three methods and reach
  `private index_tp_theta_`, plus spawn a new `kTypeName`/factory type. This is **general
  infrastructure** for any X-coupled-to-CDM model, kept opt-in so β=0 / non-coupled runs are
  byte-identical (the original registered `theta_cdm` unconditionally, paying the cost on every
  synchronous run).

- **② (Minor) Composite stress-energy must become coupling-aware.** The scalar field's momentum
  `(ρ+p)θ_φ` depends on the *sibling's* `θ_cdm`. Composites today sum independent children's
  `StressEnergy`. **Resolution:** the composite adds the one cross-term after summing, exactly
  symmetric to the existing `AddCouplingDerivs`. No new infrastructure — it extends an existing
  pattern.

- **③ (Future, non-blocking) Composites cannot modify a child's *background*.**
  `AddCouplingDerivs` is perturbation-only; `CompositeSpecies::BackgroundDerivs` just sums
  children. Type-3 sidesteps this because β is a self-contained reparametrization the scalar
  field owns (its ρ, P, and KG). A true energy-exchange model (`Q ≠ 0`) would need a
  background-coupling hook the composite currently lacks. **Flagged for the roadmap; not built
  here.**

## 3. Decisions (locked during brainstorming)

| Decision | Choice |
|---|---|
| Sequencing | **One integrated branch** (infra + physics land together). |
| Fidelity target | **Figure 2** (`T3_ph` ratio plots), large negative β / strong drag in scope. |
| Gauge | **Synchronous only**; Newtonian postponed; **fail loudly** if Newtonian + β≠0. |
| Abstraction level | **Faithful Type-3 + reusable seams** (no speculative general framework). |
| CDM θ in synchronous | Opt-in **construction flag**, used directly (no subclass). |
| β on the scalar field | **Parametrize** the existing `ScalarFieldSpecies` (default 0), no subclass. |
| Potential V(φ) | **Injectable** — extracted into an injectable bundle (see §4.1); step 1. |
| Coupling owner | Dedicated `Type3Species : CompositeSpecies` (the `IDM_DR_IDR` precedent). |

## 4. Components

Three species pieces plus a potential refactor, mirroring `IDM_DR_IDR`. The phrase "who owns β"
matters: β modifies the scalar field's *own* equations (owned by the scalar-field child) and the
*coupling* (owned by the composite, which reads β back from the child via an accessor).

### 4.1 Step 1 — Injectable scalar-field potential (standalone refactor)

`V/dV/ddV` are called only in `ComputeBackground` (`scalar_field.cpp:106-108`) and the attractor
IC (`scalar_field.cpp:85`); the perturbation code reads `V''` from `pvecback`
(`scalar_field.cpp:199`). So the potential is **background-only** — never in the per-k hot loop —
and `std::function` indirection costs nothing measurable.

- Introduce a small value type bundling the three callables, signature **(A)**:
  ```cpp
  struct ScalarFieldPotential {
    std::function<double(double phi, const std::vector<double>& params)> V, dV, ddV;
  };
  ```
- `ScalarFieldSpecies` continues to own `scf_parameters_` and `scf_tuning_index_`; the shooter
  tunes `params[tuning_index]` exactly as today (`scalar_field.cpp:485`,`499-626`). The injected
  functions supply the *form*; the species owns the *numbers* and the shooting knob.
- **Default** bundle = today's 4-parameter built-in `V = exp(-λφ)·((φ-B)^α + A)`
  (`scalar_field.cpp:26-58`), wrapped as lambdas. Existing `.ini` files stay **byte-identical**.
- Constructor gains an optional `ScalarFieldPotential` argument (defaulted). Internal `V_scf`
  etc. become calls into the bundle, passing `scf_parameters_`.

This independently removes the "potential-picking inside the species" smell: future models bring
their own potential without touching `ScalarFieldSpecies`.

### 4.2 `CDMSpecies` — opt-in coupled velocity (general infra)

Add a construction flag, e.g. `bool coupled_` (default `false`), reached via constructor /
`CreateAll`. When `false`, every path is byte-identical to today. When `true`:
- `RegisterPerturbationIndices`: also allocate `idx_theta` in synchronous gauge.
- `PerturbDerivs` synchronous branch: `dy[idx_delta] = -(θ + metric_continuity)` and
  `dy[idx_theta] = -(a'/a)·θ` (+ `metric_euler`, ≡0 in synchronous). This is **free-streaming
  only** — the momentum-transfer source is added by the composite.
- `RegisterTransferSourceIndices`: register the θ transfer-source in synchronous too.
- `ApplyInitialConditions`: θ_cdm = 0 (adiabatic).
- `PrintVariables`: output the real θ instead of the forced 0 (`cdm.cpp:194-196`).
- `StressEnergy` (`cdm.cpp:96`) already reads θ via the `idx_theta >= 0` sentinel → `ρ·θ` falls
  out for free; `CopyPerturbationsAcrossSwitch` / `PerturbSynchronousToNewtonian` are already
  θ-sentinel-generic.

### 4.3 `ScalarFieldSpecies` — gains β (default 0)

β is a constructor parameter (default 0 → byte-identical). It owns everything β does to the
field's *own* equations (background and the free-streaming KG / stress-energy); see §5.

### 4.4 `Type3Species : CompositeSpecies` — owns the coupling

Children (constructed directly, like `IDM_DR_IDR_Species`): one **coupled** `CDMSpecies`, one
**β** `ScalarFieldSpecies`. Nested `PerturbLayout { CDMSpecies::PerturbLayout cdm;
ScalarFieldSpecies::PerturbLayout scf; }`. Reads β back from the scalar-field child via an
accessor. Implements:
- `AddCouplingDerivs` → the φ-KG θ_cdm source **and** the full O(β) CDM-Euler momentum-transfer
  source (§5).
- `StressEnergy` override → sum children, then add the `−(2β/3)·Z̄²·θ_cdm` cross-term (§5).
- Standard composite forwarding for layout registration, IC, source filling, output (mirror
  `IDM_DR_IDR_Species`).

**Data flow per RHS call (synchronous):** module fills `PerturbScalarContext` →
`Type3Species::PerturbDerivs` calls the scalar-field child (own modified KG) and the CDM child
(free-streaming δ, θ) → `AddCouplingDerivs` reads `δφ, δφ', θ_cdm` + background `φ'_bg, dV,
ρ_cdm` and writes coupling into `dy[idx_phi_prime]` and `dy[idx_theta_cdm]`. For the Einstein
equations, `Type3Species::StressEnergy` returns summed-plus-cross-term totals.

## 5. Equations (synchronous gauge, code-faithful)

Conventions: `δφ = y[idx_phi]`, `δφ' = y[idx_phi_prime]`, `θc = y[idx_theta_cdm]`,
`Z̄ = −φ'_bg/a`, shared denominator `D ≡ 3·ρ_cdm·a² − 2β·Z̄²`. All forms reduce to current code
at β=0.

### 5.1 Background — `ScalarFieldSpecies` (owns β)
- `ρ_φ = [(1−2β)·φ'²/(2a²) + V] / 3`
- `P_φ = [(1−2β)·φ'²/(2a²) − V] / 3`
- KG: `φ'' = −a·(2H·φ' + a·dV/(1−2β))` — `1/(1−2β)` on `dV`
- `P'_φ = φ'·(−(1−2β)·φ'·H/a − (2/3)·dV)` — `(1−2β)` on the kinetic term only

### 5.2 Perturbations — child-written (free-streaming)
- Scalar field's own modified KG (only change vs `scalar_field.cpp:208-210` is `k² → k²/(1−2β)`):
  ```
  dy[idx_phi]       = δφ'
  dy[idx_phi_prime] = −2(a'/a)·δφ' − metric_continuity·φ'_bg − (k²/(1−2β) + a²·V'')·δφ
  ```
- Coupled CDM continuity: `dy[idx_delta_cdm] = −(θc + metric_continuity)`
- Coupled CDM Euler free part: `dy[idx_theta_cdm] = −(a'/a)·θc`

### 5.3 Perturbations — composite-written (coupling), via `AddCouplingDerivs`
- φ-KG source: `dy[idx_phi_prime] += 2·(β/(1−2β))·φ'_bg·θc`
- CDM Euler momentum-transfer source (orig `perturbations.c:6916-6925`):
  ```
  dy[idx_theta_cdm] +=
      − k²·2β·[(1−2β)·Z̄²·δφ' + dV·δφ] / [(1−2β)·D]
      + 4β·dV·φ'_bg·(k²·δφ/φ'_bg − 2β·θc) / [(1−2β)·D]
      − (6β·(a'/a)·Z̄² + 4β·dV·φ'_bg)·θc / D
  ```

### 5.4 Stress-energy (feeds the Einstein equations)
- Scalar-field child: kinetic terms of `δρ_φ, δp_φ` pick up `(1−2β)` (its own β-mod, not a
  cross-term); `(ρ+p)θ_φ = (1/3)(k²/a²)·φ'_bg·δφ` is **unchanged**.
- Composite `StressEnergy` override, after summing children:
  `(ρ+p)θ_total += −(2β/3)·Z̄²·θc`.

### 5.5 Physical / numerical notes
- θc initial condition = 0 (adiabatic). The momentum source is the only thing that moves it, so
  **β=0 ⟹ θc ≡ 0 ⟹ existing synchronous runs unchanged**; no residual-gauge fix needed.
- For β<0 (target regime), `D = 3ρ_cdm·a² + |2β|·Z̄² > 0` always — no singularity.
- Add a `class_test` for the physicality bound **β < 1/2** (ghost/strong-coupling above), and a
  guard against `D → 0` (relevant only for the β=+0.1 branch).

## 6. Construction, input, shooting, guards

- **Input parameter:** β read from a new key, proposed **`scf_veta`** (continuity with the
  original; maps to β). Absent or `0` ⟹ no composite (plain CDM + plain scalar field). Present
  and nonzero ⟹ build the `Type3Species` composite.
- **Factory:** add `Type3Species::CreateAll` to `kAllSpeciesFactories` (`all_species.h`). It
  reads `scf_veta`; if active, it reads the scalar-field inputs (`Omega_scf`, `scf_parameters`,
  `scf_tuning_index`, `attractor_ic_scf`, `scf_shooting_parameter` — reuse the logic in
  `ScalarFieldSpecies::CreateAll`) and constructs the β scalar-field child (with the 1EXP
  potential lambdas, `params=[V₀, λ]`, `tuning_index=0`) and the coupled CDM child.
- **Ω-budget routing:** when the coupling is active, the **cdm** budget is routed to the
  composite's CDM child and the standalone-CDM slot is zeroed, following the
  `ReadCoupledOmegaBudget` / `IDM_DR_IDR` precedent (so `CDMSpecies::CreateAll`, which reads
  `ctx.omega_budget->cdm` at `cdm.cpp:248-253`, does not also build a CDM). *Implementation note:
  confirm the exact budget edit against the input module.*
- **Shooting:** V₀/Ω_φ closure is unchanged — the composite's scalar-field child carries the
  existing `needs_shooting_` / `shooting_target_` / `ComputeShootingGuess` path (β is an input,
  not a shot parameter). θ_s shooting is the standard input-module path. The composite forwards
  `GetShootingTargets` from its scalar-field child.
- **Gauge guard:** `class_test` erroring clearly ("Type-3 coupling is implemented in synchronous
  gauge only") when `gauge == newtonian` and β≠0. The CDM velocity flag itself still works in
  both gauges (Newtonian CDM already carries θ); only the *coupling source* is synchronous-only.

## 7. Testing & validation

Per repo conventions: no bit-identical requirement for physics changes (~0.1% tolerance, handle
`Cl^TE` zero-crossings); use the classyref / `COMPARE_OUTPUT_REF` workflow for the
refactor/no-op gates; don't gate refactors on committed goldens that are stale under ffast-math —
use HEAD-vs-master A/B with a scale-relative metric.

1. **Step 1 (potential injection) — byte-identical.** Default-potential build must reproduce
   master for all existing scalar-field runs (`COMPARE_OUTPUT_REF` / classyref). Strong gate; a
   pure refactor.
2. **Step 2 (β param, default 0) — byte-identical** for β=0 / all existing runs.
3. **Step 3 (CDM coupled flag, default off) — byte-identical** for all non-coupled runs.
4. **Step 4 (Type-3 composite):**
   - **β=0 coupling-active sanity:** with the composite built but β=0, θc must stay 0 and outputs
     must match plain CDM + scalar field to ~0.1% (the coupling is a no-op at β=0).
   - **β≠0 against Figure 2:** for a range of negative β, form the ratios to the β=0 run and check
     the qualitative features — large-scale-only `Cl^TT` ISW signature, `P(k)` suppression with
     the turn-up below β ≈ −10² (reproduce `cl_ratio_beta` / `pk_ratio_beta`).
   - **Internal consistency:** Ω_φ closure via shooting converges; the strong-drag (large |β|)
     system integrates stably with the default ndf15 stiff evolver.
5. **Unit-style coverage** for the coupling terms (in the spirit of `photons_formula_test` /
   `species_type_name_test`): assert `AddCouplingDerivs` and the `StressEnergy` cross-term vanish
   at β=0, and spot-check the coupling expressions against hand-evaluated values.

## 8. Build order (single branch)

1. **Potential injection** refactor (§4.1) — verify byte-identical.
2. **`CDMSpecies` opt-in coupled flag** (§4.2) — off by default, verify byte-identical.
3. **`ScalarFieldSpecies` β parameter** (§4.3, §5.1, §5.2 scalar-field part) — default 0, verify
   byte-identical.
4. **`Type3Species` composite** (§4.4, §5.3, §5.4) + factory/budget/guards (§6) — the coupling.
5. **Validation** against Figure 2 (§7).

## 9. Follow-ups (tracked, not built here)

- Newtonian-gauge coupling (derive from the action / gauge-transform the synchronous equations;
  unlocks the synchronous↔Newtonian agreement cross-check as an internal ground truth).
- A composite **background-coupling** hook (constraint ③) to enable energy-exchange (Type 1/2,
  `Q ≠ 0`) models.
- Wrapper (`classy`) exposure of `scf_veta` and `theta_cdm` output, if not automatic.
