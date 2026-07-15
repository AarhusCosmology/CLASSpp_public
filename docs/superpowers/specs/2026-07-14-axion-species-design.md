# Axion support for CLASSpp — design

**Date:** 2026-07-14
**Status:** approved autonomously (user on vacation; requested end-to-end delivery)
**Reference:** Vivian Poulin's AxiCLASS (github.com/PoulinV/AxiCLASS), papers arXiv:1806.10608,
arXiv:1811.04083, arXiv:1905.12618. A shallow clone of AxiCLASS lives in the session scratchpad
for formula extraction; line references below are to its `source/*.c`.

## Goal

Add axion-like scalar-field physics to CLASSpp, covering the two regimes AxiCLASS serves:

1. **Exact treatment** — an axion potential for the existing `ScalarFieldSpecies`, evolving the
   Klein-Gordon equation throughout. Valid whenever the field oscillates slowly enough for the
   ODE to resolve (axion dark energy, late oscillation onset, m up to ~10³ H₀).
2. **Effective EDE fluid** — AxiCLASS's `pheno_axion` parametrization: a fluid whose equation of
   state transitions from w_i = −1 to w_f = (n−1)/(n+1) at scale factor a_c, with the
   Poulin-et-al k- and a-dependent effective sound speed. This is the numerically cheap
   parametrization used in published EDE analyses, valid deep into the oscillatory regime where
   the exact ODE is intractable (e.g. a_c ~ 10⁻³·⁵, m/H₀ ~ 10⁵⁺).

## Approaches considered

- **(A) Axion potential only (full KG).** Minimal, but EDE science regime unreachable
  (resolving ~10⁵–10⁶ field oscillations per k-mode).
- **(B) Full AxiCLASS parity: KG→fluid runtime switching** (`scf_evolve_as_fluid`,
  m/H threshold). Requires a new approximation-switching seam in *background* integration,
  which CLASSpp does not have (perturbations have one; background does not). Too invasive to
  build and validate unattended; also the AxiCLASS implementation carries several
  numerical-security hacks (`n_axion_security`, `security_small_Omega_scf`) that signal
  fragility.
- **(C) = (A) + pheno-axion EDE fluid.** Chosen. Both components are independently useful,
  independently testable, fit the species-plugin architecture with no module changes, and
  overlap in a validation window (moderate a_c) where they can be cross-checked. (B) becomes a
  follow-up issue.

## Component 1: axion potential for `ScalarFieldSpecies`

### Physics

V(φ) = m² f² (1 − cos(φ/f))ⁿ in CLASS internal units: φ and f in reduced-Planck-mass units,
m in 1/Mpc, so V carries 1/Mpc² like the existing potentials (cf. AxiCLASS `V_p_scf`,
background.c ≈ line 3680). With u ≡ 1 − cos(φ/f) and using sin² = u(2−u), cos = 1−u:

- V   = m² f² uⁿ
- V′  = m² f n uⁿ⁻¹ sin(φ/f)
- V″  = m² n uⁿ⁻¹ [ (n−1)(2−u) + (1−u) ]

The V″ form is exactly equivalent to the naive second derivative but contains no negative
powers of u, so it is regular at the minimum u → 0 for all n ≥ 1 (AxiCLASS's version
"bugs sometimes for n=1" — do **not** port their formula). n may be non-integer, n ≥ 1.

### Interface

New factory `AxionScalarFieldPotential()` in `species/scalar_field_potential.{h,cpp}` returning
a `ScalarFieldPotential` bundle with params layout **[m, f, n, Θ_ini]** (documented in the
header). The bundle's `shooting_guess` uses the frozen-field approximation
(ρ_today ≈ V(φ_ini) when the field has not yet thawed):

- m_guess = sqrt( 3 H₀² Ω / ( f² (1 − cos Θ_ini)ⁿ ) ), with the matching analytic dxdy
  (dΩ/dm = 2 m f² (1−cos Θ_ini)ⁿ / (3 H₀²)).

This guess is exact for a still-frozen field and a serviceable Newton start otherwise.

### Input vocabulary (flat keys, AxiCLASS-compatible names)

| key | meaning | default |
|---|---|---|
| `scf_potential` | `axion` selects the new bundle; absent keeps the historical default potential; any other value errors | historical |
| `m_axion` | m in 1/Mpc — optional Newton seed for the shooting (replaces the frozen-field guess) | frozen-field guess |
| `f_axion` | f in reduced-Planck units | required |
| `n_axion` | exponent n ≥ 1 | 1 |
| `Theta_initial_scf` | Θ_ini ∈ (0, π); φ_ini = Θ_ini · f | required |
| `Omega_scf` / closure | shooting target for m (tuning index 0) | one of the two required |

Wiring: `ScalarFieldSpecies::CreateAll` reads `scf_potential`; when `axion`, it builds
scf_parameters = [m, f, n, Θ_ini], forces tuning_index = 0, forces **frozen ICs**
(φ_ini = Θ_ini·f, φ′_ini = 0, attractor off) and throws if the user explicitly requests
`attractor_ic_scf = yes` (a (1−cos) potential has no exponential attractor) or passes
`scf_parameters` (superseded by the named keys). m is **always resolved by shooting** against
Omega_scf (explicit value, or the budget-closure override — closure mode arms the same
target). There is no "specify m, derive Ω" mode: the budget architecture has no place for a
floating Ω, and closure mode IS the floating-Ω mode. The Type-3 composite keeps its hardwired
1EXP potential; axion+Type3 is out of scope.

### What does NOT change

No changes to `ScalarFieldSpecies` background/perturbation equations — the KG machinery,
including the Newtonian-gauge q-variable (#348), is potential-agnostic. Perturbations use
V′ and V″ through the existing indices.

## Component 2: pheno-axion EDE fluid

### Physics (re-derived; port physics, not AxiCLASS bugs)

Let Δw = w_f − w_i, r = 3Δw/ν (ν = `nu_fld`, transition rapidity), x(a) = (a_c/a)^r:

- **w(a)** = w_i + Δw / (1 + x)                       [background.c ≈ 770]
- **dw/da** = Δw · (r/a) · x / (1+x)²                 (AxiCLASS's expression at ≈ line 815 has a
  C comma-operator bug and mixes (1+w_f) with Δw — derive fresh, verify by finite difference)
- **∫ₐ¹ 3(1+w)/a′ da′** = 3(1+w_i) ln(1/a) + ν · [ ln(1 + a_c^−r) − ln(1 + (a/a_c)^r) ]
  (substitute u = (a/a_c)^r; the Δw part integrates to (3Δw/r) ln(1+u) = ν ln(1+u)).
  Verify against numerical quadrature in a unit test; AxiCLASS's equivalent closed form is at
  ≈ line 850. ρ_fld(a) = ρ_fld,0 · exp(integral). Note: AxiCLASS's input.c 4138/4154 version
  of this exponent agrees with the derivation above only for w_i = −1 (it hardwires mixed
  (1+w_f)/(w_f−w_i) conventions); the fresh derivation is authoritative, the quadrature test
  is the referee.
- **Effective sound speed** (fluid rest frame), from 1806.10608:
  c_s²(k,a) = [2 a² (n−1) ϖ² + k²] / [2 a² (n+1) ϖ² + k²], with
  ϖ(a) = ω_axion · a^(−3(n−1)/(n+1))  [perturbations.c ≈ 7802]. Limits: k→∞ ⇒ 1; k→0 ⇒
  (n−1)/(n+1) = w_f.
- **Derived mapping** [background.c 982–1023, Eqs. 27/28/30 of 1806.10608]: from
  (a_c, Ω_fld_ac, Θ_i, n) derive m_fld, α_fld and ω_axion = H₀ · m_fld · (1−cos Θ_i)^((n−1)/2) · G(a_c,n),
  where G uses Euler Γ functions — use `std::tgamma` (no GSL dependency). Port verbatim,
  including the p = 1/2 (a_c < a_eq) vs 2/3 branch and E(a_c).
- **Density normalization** [input.c 4113–4169, fully algebraic — no shooting]:
  - `Omega_fld_ac` (in ρ_crit,0 units) ⇒ Ω0_fld = Ω_fld_ac / exp(integral(a_c)).
  - `fraction_fld_ac` f ⇒ Ω_fld_ac = Ω_tot(a_c) · f/(1−f) with Ω_tot(a_c) from the standard
    radiation+matter(+Λ) scalings — extract the exact expression from input.c.

### Interface

- `enum equation_of_state { CLP, EDE, PhenoAxion }` (source/background.h; string value
  `pheno_axion` in `fluid_equation_of_state` parsing — keep AxiCLASS's .ini vocabulary).
- New class `AxionEDEFluid : public FluidSpecies` in `species/axion_ede_fluid.{h,cpp}`
  (precedent: `PpfFluid` is likewise a FluidSpecies subclass built by
  `FluidSpecies::CreateAll`). It owns a_c, n, ν, w_i, w_f, Θ_i, ω_axion and overrides:
  - `ComputeWFld` — sigmoid w(a), analytic dw/da, closed-form integral;
  - a **new virtual sound-speed hook** on FluidSpecies: `virtual double Cs2(double k, double a) const`
    returning `cs2_fld_` in the base class; the subclass returns the GDM formula. All fluid
    perturbation code paths (PerturbDerivs, ApplyInitialConditions, StressEnergy if it uses
    cs2, and PpfFluid where applicable) must be swept to route through `Cs2(k, a)` —
    a plan task enumerates the call sites.
  - `ApplyInitialConditions` — AxiCLASS's EDE-adapted adiabatic ICs
    (perturbations.c ≈ 5784–5807: the (32 + 6 c_s² + 12 w_f) denominator variant); these are
    regular as w → −1 because δ, θ ∝ (1+w).
- `FluidSpecies::CreateAll` constructs `AxionEDEFluid` when eos = pheno_axion. PPF is not
  meaningful for this species (w > −1 for all a > 0; the GDM c_s² must enter the true fluid
  equations): if the user sets `use_ppf = yes` together with pheno_axion, **throw**; default
  use_ppf for pheno_axion is no. Any base-class "w must be > −1 without PPF" validation must
  evaluate w at finite a (w_i = −1 is only the a → 0 asymptote).
- Budget: Ω0_fld is derived, so the closure species must be Lambda (the default); throw a
  clear error if the user simultaneously fixes Omega_fld/Omega_Lambda inconsistently or passes
  more than one of {Omega_fld, Omega_fld_ac, fraction_fld_ac} (mirror input.c 4124–4128).

### Input vocabulary

| key | meaning | default |
|---|---|---|
| `fluid_equation_of_state` | `pheno_axion` | — |
| `n_pheno_axion` | n (sets w_f = (n−1)/(n+1)) | required |
| `log10_axion_ac` or `a_c` | transition scale factor (exactly one) | required |
| `fraction_fld_ac` or `Omega_fld_ac` or `Omega_fld` | density normalization (exactly one) | required |
| `Theta_initial_fld` | Θ_i ∈ (0, π), enters ω_axion only | required (AxiCLASS errors without it) |
| `nu_fld` | transition rapidity ν | 1 (AxiCLASS input.c:7742) |
| `w_fld_i`, `w_fld_f` | optional overrides of −1 and (n−1)/(n+1) | derived |

Where the table says "AxiCLASS default (extract)", the implementer reads the default from
AxiCLASS input.c and records it in the .ini documentation.

### Validation rules

n ≥ 1; 0 < a_c < 1; 0 < fraction_fld_ac < 1; Θ_i ∈ (0, π); ν > 0; mutually exclusive key
groups enforced with clear messages.

## Testing

Unit tests (gtest, alongside existing `species/*_test.cpp`):
1. Axion potential: V′ and V″ against central finite differences of V over a φ grid including
   the minimum, for n ∈ {1, 2, 3, 2.5}; V″ regularity at φ = 0.
2. Shooting guess: for a frozen-field configuration, guess reproduces Ω to first order.
3. w(a): sigmoid limits (w → w_i for a ≪ a_c, → w_f for a ≫ a_c); dw/da vs finite difference;
   closed-form integral vs adaptive quadrature of 3(1+w)/a over several (a_c, n, ν).
4. c_s²(k,a): k→0 and k→∞ limits; monotonicity in k.
5. Ω0_fld ↔ Ω_fld_ac round-trip through the closed form.

Integration tests (.ini level, following the repo's existing pattern; `git add -f` for
gitignored .ini fixtures):
6. scf-axion smoke run (background + Cls) with m ~ few·H₀ (axion DE regime) — runs to
   completion, Ω_scf hit by shooting to the standard tolerance.
7. pheno-axion EDE run (a_c = 10⁻³·⁵, f = 0.1, n = 3) — runs to completion; f_EDE(a_c)
   reproduces the input fraction to ~1%; Cl TT differs from ΛCDM (feature actually does
   something).
8. Decoupling: f → 10⁻⁸ recovers the ΛCDM Cls to < 0.01% (scale-aware comparison, no blind
   max-rel-diff; Cl^TE zero-crossing rule applies).
9. Cross-validation (physics sanity, loose tolerance): late-transition case (a_c ~ 0.3,
   small m) run both as scf-axion (exact KG) and pheno-axion fluid; background ρ(a) asymptotes
   (frozen plateau, a^(−3(1+w_f)) tail) agree at the ~10% level.

Best-effort (not a merge gate): build the cloned AxiCLASS (old C CLASS 2.x; may not compile on
current macOS/clang — skip gracefully) and compare Cl TT for one matched pheno-axion
configuration at ~1% tolerance.

Regression: full existing test suite must stay green; both features are strictly opt-in
(no default-path output change), so `COMPARE_OUTPUT_REF` against classyref must be unaffected.

## Post-validation amendments (2026-07-14, after integration testing)

Integration validation (plan Task 7) found two pre-existing module checks that were
incompatible with pheno_axion; both fixes are species-owned virtuals, keeping physics
decisions out of module code:

1. **Phantom-divide guard.** `perturb_init`'s historical inline test
   `(1+w(0))(1+w(1)) ≤ 0` fired on the sigmoid's exact w(0) = −1 asymptote (a false
   positive: w > −1 strictly for all a > 0). Replaced by
   `virtual bool FluidSpecies::ReachesPhantomDivide()` — base class reproduces the
   historical endpoint test exactly (CLP/EDE behavior unchanged, including blocking exact
   w = −1 contact, whose true-fluid equations divide by 1+w); `AxionEDEFluid` returns
   false. This implements the "evaluate w at finite a" requirement stated in Component 2's
   interface section.
2. **HyRec is not supported with pheno_axion in v1.** HyRec reconstructs its internal
   dark-energy density from one CPL (w0, wa) pair applied to Ω_Λ+Ω_fld (hyrec/history.c);
   no such pair represents Λ plus a frozen→dilution sigmoid, and the old tangent-at-a=1
   fill collapsed HyRec's recombination solver (z_rec < 500). Replaced by
   `virtual bool FluidSpecies::HyrecCplApproximation(double* w0, double* wa)` — base class
   returns the historical tangent (behavior unchanged for CLP/EDE); `AxionEDEFluid`
   refuses, and thermodynamics throws an actionable error directing to
   `recombination = RECFAST` (the default, which reads the true background table and is
   unaffected). Follow-up issue: feed HyRec the true H(z) — this would also correct the
   pre-existing approximation for scalar-field dark-energy runs, whose density HyRec's
   internal model ignores entirely.
   Consequently the "recombination-code cross-check" in Testing is amended to: RECFAST
   runs end-to-end; HyRec fails fast with the actionable message.

## Non-goals / follow-ups (file as issues)

- KG→fluid runtime switching for the scalar field (`scf_evolve_as_fluid` parity) — needs a
  background approximation-switch seam; design separately.
- axionCAMB-style perturbation treatment (`scf_evolve_like_axionCAMB`).
- Axion potential inside the Type-3 composite.
- `phi_2n` / `axionquad` potentials (trivial to add later via the same bundle seam if wanted).
- Derived-parameter output of the f_EDE peak (a_peak, f_ede_peak) via the wrapper.

## Deliverable

One feature branch (`axion-species`), two logical commits (component 1, component 2), one PR.
The pyx/pxd wrapper is untouched unless a struct member changes (none planned — all new state
lives in species classes; `cclassy.pxd` is generated, never hand-edited).
