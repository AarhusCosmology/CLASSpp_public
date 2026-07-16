# Thermal Axion Hot Dark Matter (`ncdm_axion`) — Design

**Date:** 2026-07-16
**Status:** Approved autonomously (Thomas on vacation; pre-authorized full pipeline: brainstorm → spec → plan → subagent implementation).
**References:** arXiv:1307.0615 (Archidiacono, Hannestad, Mirizzi, Raffelt, Wong — "Axion hot dark matter bounds after Planck 2013"); hep-ph/0504059 (Hannestad, Mirizzi, Raffelt); arXiv:0803.1585, 1004.0695 (Hannestad et al.).

## Goal

Add thermal axions as a hot-dark-matter species to CLASSpp. A thermal axion is a
relic that decoupled from the SM plasma while relativistic — exactly like a massive
neutrino, except it is a **boson** (Bose-Einstein statistics) with **one** internal
degree of freedom (real scalar) and a decoupling-dependent temperature. PR #370
added the axion as dark energy (scf potential + pheno-axion EDE fluid); this spec
covers the complementary HDM regime (eV-scale mass, free-streaming relic).

## Physics

Phase-space distribution (CLASS reduced momentum q = p/(T_a·a), PSD in CLASS
normalization where ncdm_standard's Fermi-Dirac at deg=1 is `2/(2π)³ · 1/(e^q+1)`,
i.e. deg=1 ≡ 2 fermionic dof for ν+ν̄):

```
f0_axion(q) = 1/(2π)³ · 1/(e^q − 1)          (deg = 1 ≡ ONE bosonic dof)
```

Temperature today from instantaneous-decoupling entropy conservation:

```
T_a/T_γ = (43/11 / g*S(T_D))^(1/3)
```

with `g*S(T_D)` the entropy degrees of freedom at axion decoupling
(43/11 = 2 + (7/8)·6·(4/11) is today's g*S; g*S = 10.75 reproduces the neutrino
value (4/11)^(1/3) ≈ 0.7138). For hadronic axions decoupling around the QCD epoch,
g*S ~ 15–70 → T_a/T_γ ~ 0.38–0.64.

Everything downstream is standard NCDM physics and is inherited: ΔN_eff in the
massless limit is (8/7)·(1/2)·(T_a/T_ν)⁴ per dof; Ω_a h² = m_a n_a/ρ_c with
n_a = (ζ₃/π²)(T_a)³ per dof; free-streaming suppression of P(k); Boltzmann
hierarchy driven by dlnf0/dlnq.

## Approaches considered

**A. Subclass `NCDMSpecies`, override the PSD (chosen).** Exactly the
`GreyBodyNCDMSpecies` precedent: DeferInit constructor → read params → set up PSD →
`BuildQuadratureAndMass` + `ResolveMassOmegaClosure`. ~120 LOC + tests. Statistics
enters the shared machinery *only* via the quadrature products (`q_`, `w_`,
`dlnf0_dlnq_`, `M_`), so all background/perturbation/FA/tensor/source code is
inherited unmodified. This is also what Thomas asked for.

**B. Statistics switch on `ncdm_standard`** (`statistics = bose|fermi` input).
Fewer files, but muddies deg semantics (FD deg=1 ≡ 2 dof; the axion needs 1),
weakens the species-as-plugins design (`type = ncdm_axion` in an ini documents
intent; a flag combination does not), and adds a branch to every ncdm instance.
Rejected.

**C. PSD file only (no code).** `ncdm_standard` + `use_psd_file` with a tabulated
Bose-Einstein f0 works today but has poor UX, spline/extrapolation error, and no
first-class parametrization. Rejected as the deliverable — but used as the
independent A/B validation path in tests, since it exercises identical downstream
code from an independently constructed PSD.

## Design

### New species

- **Files:** `species/axion_ncdm_species.h`, `species/axion_ncdm_species.cpp`
- **Class:** `AxionNCDMSpecies : public NCDMSpecies` (mirrors `GreyBodyNCDMSpecies`)
- **Type name:** `kTypeName = "ncdm_axion"` (family convention: `ncdm_standard`,
  `ncdm_greybody`, `ncdm_decay_dr`, `ncdm_self_interacting`)
- **PSD override:**
  ```cpp
  double AxionNCDMSpecies::EvaluatePsdAnalytic(double q) const {
    return 1.0 / pow(2. * _PI_, 3) / std::expm1(q);   // Bose-Einstein, 1 dof
  }
  ```
  `std::expm1` for small-q accuracy. The 1/q divergence at q→0 is integrable
  (measure q²·f0 → q) and no caller evaluates q ≤ 0: Gauss-Kronrod nodes are
  interior to (0,1)-mapped domains, Laguerre nodes are > 0, and the dlnf0_dlnq
  finite-difference stencil is bounded by dq ≤ (0.5−ε)·q at the first node.
  Under -ffast-math, never materialize the q=0 infinity (no eager evaluation,
  no guards that compute 1/expm1(0)).
- **Constructor** (DeferInit pattern, exactly like greybody):
  1. `NCDMSpecies(pfc, name, settings, pba, bgm, DeferInit{})` — base reads the
     shared fields (`m`, `deg`, `Omega`/`omega`, quadrature knobs…).
  2. Read `T` and `gstar_dec` via `SpeciesInput::get<double>` (presence-aware).
     Require **exactly one**; if `gstar_dec`: `T_ = cbrt((43./11.)/gstar_dec)`.
  3. Validation (throw `std::invalid_argument` with the instance name):
     - both or neither of `T`/`gstar_dec` given,
     - `gstar_dec <= 43./11.` (axion cannot be hotter than the photons it
       decoupled from; also guards div-by-zero/negative cube root),
     - `ksi_ != 0` (chemical potential undefined for this species; μ>0 BE
       diverges),
     - `use_psd_file` set (contradiction — the PSD *is* the model; use
       `ncdm_standard` for file PSDs).
  4. `BuildQuadratureAndMass(settings)` then `ResolveMassOmegaClosure(settings)`
     (dispatches the overridden PSD into quadrature; inherits the m ↔ Omega
     closure, including "both m and Omega ⇒ rescale deg" semantics).
- **`CreateAll`**: iterate `ctx.pfc->instances_with("type", kTypeName)`, mark
  `.type` consumed, construct per instance — verbatim greybody pattern.
  Dot-syntax only; **no** legacy `N_ncdm`-style counting.
- **No other overrides.** `DefaultQuadratureStrategy` stays `qm_auto` (adaptive
  Gauss-Kronrod/Laguerre selection handles any exponentially decaying PSD; the
  FD-tuned candidate rule simply won't be selected when it doesn't converge).
  `GetParam` not needed (nothing new to expose; `T` and `deg` are base state).

### Input interface

```ini
ax.type      = ncdm_axion
ax.m         = 1.0          # eV  (or ax.Omega / ax.omega — inherited closure)
ax.T         = 0.39         # T_a/T_cmb …
ax.gstar_dec = 30           # … XOR g*S at decoupling → T = (43/11/g*S)^(1/3)
ax.deg       = 1            # optional; internal bosonic dof (default 1, real scalar)
# inherited & allowed: quadrature_strategy, momenta_bins, momenta_bins_bg, max_q
# rejected: ksi, use_psd_file/psd_filename
```

Note the deliberate semantic difference from `ncdm_standard`: there deg=1 means one
*neutrino family* (2 fermionic dof); here deg=1 means one *real scalar* (1 bosonic
dof). deg=2 gives a complex scalar / generic thermal boson, so the species doubles
as a generic hot-boson model. Massless-limit default is inherited: if both m and
Omega are absent, m = 1e-5 eV (≈ dark radiation).

### Registration

- `species/all_species.h`: include header + one `SpeciesFactoryEntry` row.
- `species/species_type_name_test.cpp`: add the kTypeName assert and extend the
  pinned `expected` set (test fails otherwise — it pins the registry exactly).

### What is inherited with zero changes

Background momenta integration, scalar + tensor Boltzmann hierarchies, fluid
approximation and switch copying, transfer sources, N_eff reporting/printing,
a_ini relativistic search, halofit mass warning, output columns, shooting.
`GetNeff` normalizes to one standard massless neutrino, `factor_` and
`rho_nu_rel_` are statistics-agnostic (verified by reading the base).

## Error handling

All validation throws `std::invalid_argument` naming the instance, matching
greybody. No silent fallbacks; in particular there is **no default temperature**
(the base's 0.71611 neutrino default is wrong for an axion, so T/gstar_dec is
mandatory).

## Testing

New target `test-axion-ncdm` (`species/axion_ncdm_test.cpp`), registered in **all
three** places: `add_executable` + the `foreach(_t IN ITEMS …)` list in
CMakeLists.txt, and `TEST_TARGETS` in the Makefile (CI only builds what is named
there).

Unit assertions (analytic, tolerance ~1e-4 unless noted — quadrature-limited,
never bit-identical):

1. **Factory:** `type = ncdm_axion` instance constructs via `CreateAll`;
   registry pinned in species_type_name_test.
2. **Bose/Fermi ratio:** massless `ncdm_axion` vs massless `ncdm_standard` at the
   same T, deg=1: `GetNeff` ratio = (π⁴/15)/(2·7π⁴/120) = **4/7** — pins both the
   statistics and the 1-vs-2 dof normalization in one number.
3. **ΔN_eff:** massless axion at T: GetNeff(0) = (8/7)·(1/2)·(T/(4/11)^(1/3))⁴.
4. **Ω_a h²:** m = 1 eV, T = 0.39: matches m·n_a/ρ_crit computed from first
   principles (ζ₃, physical constants) in the test, tolerance 1e-3.
5. **Non-relativistic limit:** rho/n → m (unit-converted) for M ≫ 1 at z = 0.
6. **gstar_dec mapping:** g*S = 10.75 → T = (43/11/10.75)^(1/3) ≈ 0.71377.
7. **Error cases:** the four constructor throws.
8. **Quadrature convergence:** qm_auto rho(z=0) vs high-resolution manual grid
   (`momenta_bins` large) agree within tol_ncdm_bg-driven tolerance; also
   dlnf0_dlnq at the smallest node is near the BE small-q limit (→ −1 as q→0).
9. **PSD-file A/B (the independent path):** write a dense BE table to a temp
   file, build `ncdm_standard` with `use_psd_file=1`, same T and m, and compare
   rho at several z and Omega0 against `ncdm_axion` (tolerance ~1e-4, set by
   spline error). Since every downstream consumer sees only (q, w, dlnf0_dlnq, M),
   background agreement here implies perturbation agreement by construction.
10. **Pipeline smoke test:** build a small `Cosmology` with an ncdm_axion
    instance (background at least, following dcdm_wdm_test.cpp's pattern) and
    assert finite, sane output — proves no assert/NaN trips end-to-end wiring.

Full-run validation (manual, not CI): one `./class` run with
`test/dotsyntax_axion.ini` at default precision; eyeball Ω budget and N_eff in
the log. Committed goldens are known-stale under -ffast-math (#338) — do not gate
on them.

## Documentation

- `explanatory.ini`: short block in the ncdm section (§6) documenting
  `ncdm_axion` dot-syntax, deg semantics, and the T/gstar_dec choice, in the
  style of the `ncdm_decay_dr` note.
- `test/dotsyntax_axion.ini`: minimal working example (likely gitignored —
  stage with `git add -f` like the other dot-syntax inis).

## Non-goals

- **m_a → T_D hadronic freeze-out mapping** (the chiral-PT π+π ↔ π+a rate of
  1307.0615): model-dependent, and post-2021 lattice/EFT work has superseded the
  1993 Chang-Choi rates — hardcoding them would bake in stale physics. Users
  supply T or g*S(T_D) from their preferred production calculation.
- Chemical potential for the axion (rejected input).
- Coupling to / unification with the PR #370 EDE axion (different regime; the
  thermal axion is a free-streaming relic, not a rolling field).
- Legacy `N_ncdm`-family counted input; python wrapper changes (none needed —
  no new module members, cclassy.pxd untouched and auto-generated anyway).

## Decisions made autonomously (would otherwise have been questions)

1. Parametrize by (m, T | gstar_dec), not by m alone via the hadronic-axion
   relation — see Non-goals.
2. Type name `ncdm_axion` over `thermal_axion`/`ncdm_bose` — NCDM family naming
   convention; deg makes it a generic thermal boson anyway.
3. Mandatory explicit temperature (no neutrino-default fallback) — silently
   inheriting 0.71611 would be a wrong-physics trap.
4. deg=1 ≡ 1 bosonic dof (documented divergence from ncdm_standard's deg=1 ≡ 2
   fermionic dof).
