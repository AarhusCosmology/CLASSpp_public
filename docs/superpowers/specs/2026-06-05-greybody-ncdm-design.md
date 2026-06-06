# Grey-body NCDM species — design

**Date:** 2026-06-05
**Status:** Approved design, pre-implementation
**Supersedes:** PR #94 ("Grey body", branch `grey-body`, open since 2020) — re-implemented
from scratch against the post-refactor species system rather than rebased.

## Motivation

CLASS's non-cold-dark-matter (NCDM) machinery assumes a Fermi–Dirac phase-space
distribution (PSD). Non-thermally-produced dark matter (freeze-in, decay, etc.) has a
momentum distribution that departs from thermal. The **grey-body** parameterization
replaces the fixed FD form with a 3-parameter family:

```
f0(q) = 2/(2π)³ · (q/q0)^(α-1) · 1/(exp(α·x·q) + 1)
```

evaluated in log-space for numerical stability:

```
log f0 = (α-1)·log q − (α-1)·log q0 − log(exp(α·x·q) + 1) + log(2/(2π)³)
```

The family is specified either **directly** by (α, q0, x), or **by its velocity moments**
(model-independent): given the dimensionless moments M₂, M₃ and the ratio
`r = M₂·M₄/M₃²`, solve for (α, x, q0). The moment route is the scientifically useful one —
it lets a model be characterized by a few low-order moments without committing to a
microphysical production mechanism.

This design re-implements the feature of PR #94 cleanly. PR #94 predates the entire
species refactor and cannot be rebased; we carry forward its physics and drop its cruft.

## What PR #94 contained (reference)

Core physics (carried forward):
- Grey-body log-space PSD in `background_ncdm_distribution`.
- `greybodyparams` helper: moment → (α, x, q0) inverse solve via bisection on
  `r_M(α) = r`, using Riemann ζ extensively. Small-argument Taylor fallbacks for the
  `J_α`/`H_α` auxiliary functions near the α = −1, −2 poles.
- A hand-rolled `std::riemann_zeta` (Borwein η alternating-series) injected into
  `namespace std` under `#ifdef __clang__`, because libc++ lacks `std::riemann_zeta`.
- `greybody_moments()` computing M₂/M₃/M₄, exposed as global derived parameters
  `M_2`/`M_3`/`M_4` via `background_module` + `classy`.
- Two non-thermal quadrature methods: `qm_GB_Laguerre` (generalized Gauss–Laguerre
  tuned to the `q^(α-1) e^(-αxq)` shape) and `qm_trapz_log` (log-spaced trapezoid with
  analytic qmin/qmax determination).
- A `m_ncdm_M2` "mass × M₂" input convenience (mass = m_M2 / M₂).

Dropped (cruft):
- `if (1==1)` block dumping `ncdm_distribution_output_*.dat` every run.
- `std::cout << "Exception J/H"` debug prints; commented-out code.
- Committed `tools/.non_cold_dark_matter.cpp.swp`.
- `Pk_r_plot.py`, `cl_r_plot.py`, `plotting.py`; Xcode project edits.
- Raw `malloc`/`free` of `m_ncdm_M2` into an uninitialized pointer.
- The `#ifdef __clang__` / `namespace std` injection (replaced — see below).

## Confirmed technical constraint

`std::riemann_zeta` is **still absent from Apple clang's libc++ as of 2026-06**
(verified by compile test: Apple clang 21.0.0, `-std=c++17`). So a vendored ζ is required.
**Decision:** vendor the Borwein η-series ζ in *our* namespace, used uniformly on **all**
compilers (not gated on `__clang__`), so GCC and clang produce identical `f0` — important
for reproducibility and CI. **Location:** keep it private to the grey-body translation unit
(anonymous namespace) for now; promote to a shared `tools/` header only if a second
consumer appears.

## Architecture

### Class
New `species/greybody_ncdm_species.{h,cpp}`:

```cpp
class GreyBodyNCDMSpecies : public NCDMSpecies { ... };
```

**Inherits `NCDMSpecies`** (not `NCDMBaseSpecies`). The perturbation hierarchy, background
indices, momenta integration, sources, and ICs are identical to standard NCDM; only the
PSD `f0(q)` and the quadrature sampling differ. Inheriting from `NCDMSpecies` means we
override a handful of small hooks and inherit everything else.

Registered in `species/all_species.h`:
```cpp
SpeciesFactoryEntry{"NCDMGreyBody", &GreyBodyNCDMSpecies::CreateAll},
```
selected by dot-syntax `type = ncdm_greybody`. `CreateAll` mirrors
`NCDMSpecies::CreateAll`, iterating `instances_with("type", "ncdm_greybody")`.

### PSD made overridable (minimal base refactor)
Today `NCDMBaseSpecies::DistributionFunction` is a `static` private callback that branches
file-table-interpolation vs. Fermi–Dirac. Extract only the analytic branch into a protected
virtual:

```cpp
// NCDMBaseSpecies
virtual double EvaluatePsdAnalytic(double q) const;   // base impl: FD with ksi
```

The static `DistributionFunction` callback keeps the shared file path and calls
`p->sp->EvaluatePsdAnalytic(q)` for the analytic case. `GreyBodyNCDMSpecies` overrides it
with the log-space grey-body expression. This relocation is byte-identical for every
existing species (the FD expression is unchanged; only its call site moves).

### Construction ordering (deferred-init constructor)
`InitQuadrature` currently runs inside the base constructor, so an overridden
`EvaluatePsdAnalytic`/`FillQuadratureParams` would NOT be dispatched during construction
(virtual-call-in-constructor resolves to the base). To let the grey-body PSD drive its own
sampling in a single clean pass — with no throwaway FD quadrature polluting `deg_`/`factor_`
— we split base construction into two phases, confined to `ncdm_base_species.*` and
`ncdm_species.*`:

- `NCDMBaseSpecies`: `ReadParametersByInstance` no longer calls `InitQuadrature`. The
  quadrature + `M_` setup moves to a protected `BuildQuadratureAndMass(settings)`. The
  normal base constructor calls `ReadParametersByInstance` then `BuildQuadratureAndMass`
  (unchanged net behavior). A new protected tagged constructor
  `NCDMBaseSpecies(..., DeferInit)` runs `ReadParametersByInstance` only.
- `NCDMSpecies`: its mass/Omega closure moves from the constructor body into a protected
  `ResolveMassOmegaClosure(settings)`. The public constructor is unchanged net behavior. A
  protected tagged constructor `NCDMSpecies(..., DeferInit)` forwards to the base
  `DeferInit` constructor and runs no closure.
- `GreyBodyNCDMSpecies` constructor delegates to `NCDMSpecies(..., DeferInit)`, then in its
  body reads grey-body input, sets `alpha_/x_times_alpha_/alpham1_logq0_` and the default
  quadrature strategy, and finally calls `BuildQuadratureAndMass(settings)` (now dispatching
  the grey-body PSD) followed by `ResolveMassOmegaClosure(settings)`.

`DNCDMSpecies` (uses the normal `NCDMBaseSpecies` constructor) and `NCDMInteractingSpecies`
(uses the normal `NCDMSpecies` constructor) are untouched and byte-stable.

### Quadrature integration (species-agnostic)
Port `qm_GB_Laguerre` and `qm_trapz_log` into `tools/quadrature.cpp`. To keep
`quadrature.cpp` free of species types, the three grey-body parameters
(`alpha`, `x_times_alpha`, `alpham1_logq0`) become optional fields on the
`DistributionParams` struct (defaulted to FD-equivalents). They are filled via a virtual:

```cpp
// NCDMBaseSpecies
virtual void FillQuadratureParams(DistributionParams& p) const {}   // base: no-op
```

overridden only by `GreyBodyNCDMSpecies`. The quadrature routine reads the parameters from
the struct — no downcasting, no species-type picking. The grey-body species defaults its
sampling method through:

```cpp
virtual int DefaultQuadratureStrategy() const;   // base: qm_auto; greybody: GB method
```

### Moment inverse solve
Reworked `greybodyparams` as a self-contained helper in `greybody_ncdm_species.cpp`
(anonymous namespace):
- `riemann_zeta(s)` — Borwein η-series, our namespace, all compilers.
- Root-find α from `r_M(α) = r` using **`tools/bisection.h`** `bisect_value`, keeping the
  PR's bracket-search to bound the root before bisection.
- Closed forms for `x`, `q0`; store `alpha`, `x·alpha`, `(α−1)·log q0`.
- Retain Taylor fallbacks for `J_α`/`H_α` near α = −1, −2.

### Derived moments
`GreyBodyMoments()` method on the species (the PR's `greybody_moments`, fixed to use the
species' own quadrature rather than hardcoded index 0). Exposed **per-species** by
overriding `BaseSpecies::GetParam`:

```cpp
std::optional<double> GetParam(const std::string& name) const override;
// returns M_2 / M_3 / M_4 / alpha / q0 / x on match
```

This rides the existing `BackgroundModule::GetSpeciesParam(key, param)` →
`species->GetParam(param)` path (already used by IDM_DR), so `classy` retrieves
`GetSpeciesParam(b"<instance>", b"M_2")` keyed by the user's instance name. The
`classy.pyx` derived-parameter wiring is the **final** task (needs a per-species
derived-name convention). In moment mode the moments verify the inverse solve reproduced
the inputs; M₄ is a genuine prediction.

## Input schema (dot-syntax, explicit mode)

```
gb.type             = ncdm_greybody
gb.parameterization = direct | moments

# mass: exactly one of
gb.m    = 0.06        # eV, direct
gb.m_M2 = ...         # mass·M2 convenience; mass = m_M2 / M2  (moments mode)

# direct mode:
gb.alpha = ...   gb.q0 = ...   gb.x = ...

# moments mode:
gb.r = ...       gb.M2 = ...   gb.M3 = ...

# optional, otherwise GB default:
gb.quadrature_strategy = ...
gb.momenta_bins = ...   gb.max_q = ...   gb.deg = ...   gb.T = ...   gb.Omega = ...
```

Validation: the chosen `parameterization` requires its full field-group; a clear error on
missing/extra fields or on a "mass too low" condition (the PR threw at m < 0.01 eV in
moment mode — retained as a guarded check, not a hard global throw).

Legacy `psd_parameters = [flag,…]` `.ini` files are **not** auto-translated: the only such
files live on the unmerged `grey-body` branch, so there is nothing in the wild to preserve.

## Out of scope
- No global `background_module` M_2/M_3/M_4 fields (replaced by per-species `GetParam`).
- No plotting scripts, no Xcode project changes.
- Promotion of ζ to a shared header (deferred until a second consumer needs it).

## Testing (TDD)
1. **Unit — ζ:** `riemann_zeta` vs known values (ζ(2)=π²/6, ζ(3)≈1.2020569, ζ(4)=π⁴/90).
2. **Unit — round trip:** feed (r, M₂, M₃) → solve (α,x,q0) → `GreyBodyMoments()` returns
   the same M₂, M₃ and a self-consistent M₄.
3. **Physics — FD limit:** direct mode with α→1, x→1, q0→1 reduces to Fermi–Dirac;
   cross-check Ω_ncdm and N_eff against an equivalent `type=ncdm_standard` species to
   ~0.1%.
4. **Integration:** `test/scenarios/ncdm_greybody.ini` (both modes) runs end-to-end; a
   dot-syntax example `.ini`.
5. **Regression:** existing NCDM scenarios byte-stable (the PSD refactor is a pure
   relocation for non-grey-body species).

## Implementation order (for the plan)
1. Base refactor: extract `EvaluatePsdAnalytic`; verify all NCDM tests byte-stable.
2. Vendored ζ + unit test.
3. Moment inverse-solve helper + round-trip unit test (uses `tools/bisection.h`).
4. `GreyBodyNCDMSpecies` skeleton: `CreateAll`, registration, input parsing, both modes,
   PSD override. FD-limit physics test.
5. Quadrature: `DistributionParams` GB fields + `FillQuadratureParams` +
   `DefaultQuadratureStrategy`; port `qm_GB_Laguerre` / `qm_trapz_log`.
6. `GreyBodyMoments()` + `GetParam` override; `classy.pyx` per-species derived wiring.
7. Scenario `.ini`s + docs.
