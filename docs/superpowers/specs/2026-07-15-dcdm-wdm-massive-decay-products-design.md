# DCDM → massive decay products (dcdm_wdm) — design

**Date:** 2026-07-15
**Reference:** arXiv:2606.14849 (Bencke, Lee, Kamionkowski, "Constraints with CMB lensing on
dark matter decays to massive decay products"), local copy in `arXiv-2606.14849v1/`.
**Mode:** autonomous (user on vacation, explicitly delegated design decisions).

## 1. Goal

Implement dark matter decaying to two *massive* decay products in CLASSpp:

- A cold parent χ (fraction *f* of the dark matter) decays with rate Γ into two identical
  daughters, each of mass m_d = ε·m_χ/2 with ε ∈ [0, 1).
- Each daughter is born with fixed physical momentum p = (m_χ/2)·√(1−ε²), i.e. kick velocity
  v_kick = c·√(1−ε²); each daughter's injection energy is exactly m_χ/2.
- Daughters form a warm species whose phase-space distribution f₀(q, τ) *builds up* over time:
  a decay at conformal time τ′ injects particles at comoving momentum q = a(τ′)·p (monochromatic
  injection at the moving cutoff).

The reference paper solves an integral equation for the daughter perturbations (CLASSIER-DDM).
We instead evolve the standard per-momentum-bin Boltzmann hierarchy with a time-dependent
injection source (per Thomas: the integral-equation method is neither wanted nor necessary).
This is the ODE-based analogue of Abellán–Murgia–Poulin (arXiv:2102.12498), generalized to two
massive daughters (which is *simpler*: a single daughter species, no massless channel).

Validation deliverable: a notebook in `notebooks/` describing the feature and recreating the
physics panels of the paper's Fig. 3 (ΔP(k)/P, ΔC_ℓ^TT/C, ΔC_L^φφ/C for several lifetimes at
f = 0.1), plus limit checks.

## 2. Physics

### 2.1 Parameters (per dot-syntax instance, type `dcdm_wdm`)

| Input | Meaning |
|---|---|
| `Gamma` \| `log10Gamma` \| `lifetime` \| `log10lifetime` | decay rate, exactly one (same units/conversions as `ncdm_decay_dr`: Gamma in km/s/Mpc, lifetime in years) |
| `epsilon` XOR `vkick` | mass retention ε, or kick velocity v/c = √(1−ε²) (numerically better for small kicks) |
| `Omega_ini` \| `omega_ini` XOR `Omega_dcdmwdm` \| `omega_dcdmwdm` | parent initial abundance (CLASS `Omega_ini_dcdm` convention) or combined sector density today |
| `momenta_bins` (default 96) | number of q bins (background grid = perturbation grid) |
| `q_min_ratio` (default 1e−4) | grid lower edge as fraction of q_kick (covers injection from a = 1e−4) |
| `kernel_width` (default 1.0) | injection-kernel width in units of the local bin spacing (in ln q) |
| `l_max` (default `ppr->l_max_ncdm`) | daughter hierarchy truncation |

The paper's *f* is not a direct input: users set stable `Omega_cdm` plus this instance's
abundance, exactly as for dcdm_dr; the notebook provides the (f, ω_dm) → inputs mapping.

### 2.2 Internal normalization

Only ratios matter (paper: results independent of m_χ). We fix the daughter temperature unit
T_wdm = T_cmb and *choose* the dimensionless kick momentum q_kick ≡ p/T₀ = 10. Then the
dimensionless daughter mass is M = q_kick·ε/√(1−ε²) (from p/m_d = √(1−ε²)/ε), and the cutoff
sweeps q_cut(τ) = a(τ)·q_kick, reaching q_kick today (a = 1). ε → 0 gives M → 0 (massless
daughters ≡ dark radiation, a validation limit); ε → 1 gives M → ∞ (reject ε above 1−1e−9 with
a clear error; users wanting v_kick ≲ 1e−4 should input `vkick`, computed exactly).

The distribution amplitude is arbitrary up to the `factor_` normalization of the NCDM machinery
(ρ = factor/a⁴ · Σᵢ dqᵢ qᵢ² εᵢ fᵢ with ε(q,a) = √(q² + a²M²)). We choose factor_ = H₀² so the
integrated fᵢ are Ω-like magnitudes; no physical constant enters the amplitude because the
injection normalization (below) is anchored to ρ_dcdm.

### 2.3 Background

Parent: unchanged `DCDMSpecies` ODE, dρ_dcdm/dτ = −3ℋρ_dcdm − aΓρ_dcdm.

Daughter: per-bin ODE variables fᵢ (distribution value at qᵢ) and gᵢ ≡ (∂f₀/∂ln q)|_{qᵢ},
both starting at 0:

    dfᵢ/dτ = Jᵢ(τ),        dgᵢ/dτ = (∂J/∂ln q)|_{qᵢ}(τ)

with the *normalized smooth kernel* source (u ≡ ln q):

    Jᵢ(τ) = A(τ) · Gᵢ(τ) / D(τ)
    Gᵢ(τ) = exp( −(uᵢ − u_cut(τ))² / (2σᵢ²) ),   σᵢ = kernel_width · Δuᵢ(local)
    D(τ)  = Σⱼ dqⱼ qⱼ² εⱼ(a) Gⱼ(τ)
    A(τ)  = a Γ ρ_dcdm(τ) · a⁴ / factor_

This makes the energy-injection sum rule **exact by construction** for any grid and kernel:
(factor/a⁴)·Σᵢ dqᵢ qᵢ² εᵢ Jᵢ = a Γ ρ_dcdm, i.e. ρ̇_wdm + 3ℋ(ρ_wdm + P_wdm) = a Γ ρ_dcdm
(total parent energy m_χ per decay goes to the daughters — each carries exactly m_χ/2). The
smearing (delta → kernel of width ~1 bin) is the only discretization approximation; it converges
away with `momenta_bins`. ∂J/∂ln q treats D as q-independent: ∂G/∂u = −(u−u_cut)/σ²·G, analytic.
Injection is suppressed (J ≡ 0) while q_cut < grid lower edge (energy error ~ decayed fraction at
a = q_min_ratio, ≲ 1e−9 for the paper's Γ range). The grid extends ~5% above q_kick so the kernel
stays resolved through a = 1.

Continuum limit check: the stationary solution reproduces f̃₀(q) = Γ a N_χ /(4π q³ ℋ)|_{τ_q(q)}
(paper Eq. 6) up to kernel smearing; a golden test compares ρ_wdm(a) against the exact
injection-time integral ρ_wdm(a) = ∫ dτ′ a′Γρ_dcdm(a′)·(a′/a)³·√(ε² + (1−ε²)(a′/a)²) computed
by quadrature from the background table.

pvecback slots owned by the daughter: scalars ρ, P, pseudo-P (for PPrime), n; per-bin arrays
gᵢ and Jᵢ (needed by the perturbations, interpolated from the background table — both smooth by
construction, no splining of a step). w_bg_ᵢ = fᵢ·dqᵢ is refreshed each ComputeBackground (the
DNCDM trick) so the inherited `ComputeMomenta` produces n/ρ/P/pseudo-P unchanged.

### 2.4 Perturbations (synchronous gauge, scalars)

We evolve **unnormalized** per-bin multipoles ψ_ℓ(qᵢ) ≡ (Δf)_ℓ(qᵢ) — not the CLASS Ψ = Δf/f₀ —
because f₀ starts at zero (Ψ undefined pre-injection) while the Δf equation is closed, linear,
and regular everywhere. From ∂Δf/∂τ + ikμ(q/ε)Δf + (∂f₀/∂ln q)[metric] = (∂f₀/∂τ)δ_dcdm
(paper Eq. 8; Aoyama et al. 2014):

    dψ₀/dτ = −(qk/ε)ψ₁ + (metric_continuity/3)·gᵢ            + Jᵢ·δ_dcdm      ← source: composite
    dψ₁/dτ = (qk/3ε)(ψ₀ − 2s₂ψ₂) − (ε·metric_euler/(3qk))·gᵢ
    dψ₂/dτ = (qk/5ε)(2s₂ψ₁ − 3s₃ψ₃) − s₂·(2/15)·metric_shear·gᵢ
    ℓ ≥ 3: standard ncdm couplings; standard truncation at l_max.

These are exactly the CLASS ncdm hierarchy equations with dln f₀/dln q replaced by the
*unnormalized* gᵢ = ∂f₀/∂ln q, plus the injection source on ℓ = 0. The steep moving edge of
f₀ lives in gᵢ, so the metric terms automatically reproduce the exact injection jump conditions
(Ψ₀ ← δ_χ − ḣ/6ℋ etc.) in the continuum limit — no state jumps, no evolver restarts, no
integration-interval surgery; ndf15 resolves each bin's injection window natively.

Stress-energy (per paper Eq. 11, no f₀ weights — ψ carries them):

    δρ = (factor/a⁴) Σ dq q² ε ψ₀            δP = (factor/3a⁴) Σ dq q⁴/ε ψ₀
    (ρ̄+P̄)θ = (factor/a⁴) k Σ dq q³ ψ₁       (ρ̄+P̄)σ = (2factor/3a⁴) Σ dq q⁴/ε ψ₂

Initial conditions: ψ ≡ 0 (daughter unpopulated at τ_ini). For robustness against short
lifetimes where some injection precedes a mode's start time: ψ₀,ᵢ = fᵢ(τ_ini)·δ_dcdm_ic (adiabatic;
daughters injected superhorizon share the parent's density contrast), ψ_{ℓ≥1} = 0.

Scope guards:
- **Synchronous gauge only** (paper's gauge; Type-3 precedent). Input-time error otherwise.
  The ℓ=1 parent-velocity source is omitted (θ_dcdm ≡ 0 in synchronous) and documented.
- **No fluid approximation** for the daughter (DNCDM precedent: full hierarchy always).
- **Tensors:** daughter registers no tensor slots (zero contribution; documented). Parent dcdm
  already contributes nothing to tensors.

Classification: daughter ClustersAsMatter = true, IsColdMatterSpecies = false (inherited from
NCDMBaseSpecies) → enters δ_m/P(k) as warm matter via the generic ρ−3P tally; parent stays cold
matter. Composite EnergyType::Other (background ρ_m += ρ−3P, ρ_r += 3P split, as DCDM_DR).
Daughter IsFreestreaming = false, GetOmega0 = 0 (decay product starts at zero), NeutrinoOmega0 =
0, Neff contribution 0 (it is empty at IC time).

## 3. Architecture

Mirrors the `DNCDM_DR_Species` composite pattern (dot-syntax instances, shooting, generic
composite child loops — including the #358 generic-layout requirement).

New files:
- `species/wdm_decay_product.{h,cpp}` — `WdmDecayProductSpecies : NCDMBaseSpecies`
  (DeferInit constructor path; builds its own log-spaced q grid + trapezoidal dqᵢ instead of the
  PSD-driven quadrature; owns fᵢ/gᵢ bi-slots, per-bin bg columns, hierarchy, stress-energy,
  transfer sources d_wdm/t_wdm with N-body gauge corrections via RhoDotOverRho override
  ρ̇/ρ = (−3ℋ(ρ+P) + aΓρ_dcdm)/ρ).
- `species/dcdm_wdm_species.{h,cpp}` — `DCDM_WDM_Species : CompositeSpecies`, children
  {DCDMSpecies, WdmDecayProductSpecies}. Owns: input parsing (CreateAll over instances with
  `type = dcdm_wdm`), gauge guard, background injection wiring (composite ComputeBackground
  writes Jᵢ columns after the child loop, using ρ_dcdm from child 0; BackgroundDerivs adds
  dfᵢ/dτ = Jᵢ, dgᵢ/dτ), AddCouplingDerivs (dψ₀,ᵢ += Jᵢ·δ_dcdm), shooting hooks
  (target/unknown pairs {Omega_ini ↔ Omega_dcdmwdm}, DCDM_DR pattern), background output
  columns `(.)rho_dcdm_<inst>`, `(.)rho_wdm_<inst>` (no collision with a coexisting DCDM_DR).
- `species/dcdm_wdm_test.cpp` — unit/golden tests (below). Registered in **both**
  `CMakeLists.txt` and Makefile `TEST_TARGETS` (recurring trap).
- Factory: `all_species.h` gets `SpeciesFactoryEntry{DCDM_WDM_Species::kTypeName = "dcdm_wdm",
  &DCDM_WDM_Species::CreateAll}`.

No changes to evolvers, perturbation module dispatch, or the Einstein-equation tally: everything
flows through the existing species interfaces. `cclassy.pxd` is auto-generated — never hand-edit.

## 4. Alternatives considered

1. **Integral-equation method (CLASSIER-DDM)** — rejected up front by Thomas; alien to the
   CLASS ODE architecture.
2. **Exact delta-function injection via integration-interval boundaries** (each bin switches on
   at τ_q with jump ICs): exact per-bin ICs, but requires ~`momenta_bins` evolver restarts per
   k-mode (ndf15 loses multistep history at each), plus invasive interval machinery in the
   perturbation module. Rejected: cost and intrusiveness for no accuracy gain over the smeared
   kernel at equal bin count (the delta is *already* coarse-grained by momentum discretization).
3. **Viscous-fluid daughter approximation** (Abellán et al. production mode): cheap enough for
   MCMC but approximate and unvalidatable without an exact reference — it needs *this* feature
   first. Deferred as a natural follow-up.

## 5. Testing & validation

Unit/golden (`dcdm_wdm_test`, fast, CI):
1. Background energy conservation: ρ_wdm(a) vs the exact injection-time quadrature (§2.3),
   tolerance ~1e−3 (kernel/grid resolution), plus the exact sum-rule identity at the ODE level.
2. ε → 0 limit: background ρ_wdm tracks ρ_dr of an equivalent `dcdm_dr` run (sub-0.1%).
3. ε → 1⁻ (vkick = 1e−4): observables ≈ ΛCDM with identical ω_m budget.
4. Gauge guard: newtonian + dcdm_wdm → clear input error.
5. Composite layout invariant (mirrors composite_layout_test expectations for the new type).

Integration: full C_ℓ run vs master must be unchanged when the species is absent
(COMPARE_OUTPUT_REF); verification of changed-physics runs uses scale-relative metrics with
~0.1% tolerance and Cl^TE zero-crossing handling — never blind max-rel-diff, and never a
bit-identical requirement.

Notebook `notebooks/ddm-massive-decay-products.ipynb`:
- Feature description, parameter mapping (f, Γ, ε) → inputs, internal-normalization notes.
- Limit checks: f → 0 ≡ ΛCDM; ε → 0 ≡ dcdm_dr; vkick → 0 ≡ ΛCDM.
- Recreate paper Fig. 3 panels: ΔP(k)/P, ΔC_ℓ^TT/C_ℓ, ΔC_L^φφ/C_L for Γ⁻¹ ∈ {0.1, 1, 10} Gyr at
  f = 0.1 with v_kick at the paper's quoted 1σ bounds (read off Fig. 2; exact reproduction of
  the bounds themselves — the generalized Fisher pipeline — is out of scope), with the expected
  free-streaming scales k_fs and ℓ_fs overplotted (paper Eq. 1).
- Momentum-grid convergence check (bins 48/96/192).

Success criteria: suppression sets in at the predicted k_fs; C_L^φφ suppression ~5–10% for the
paper's 1σ models; limits reproduce ΛCDM/dcdm_dr within tolerance; energy conservation holds.

## 6. Risks / open points

- Kernel smearing broadens the f₀ edge → slightly smeared free-streaming cutoff. Controlled by
  `momenta_bins`; quantified by the convergence test.
- gᵢ near the edge is large (steep lnf drop) — but always multiplied into integrals with weight
  fᵢ; the analytic-ODE construction (no numerical differentiation of a step) keeps it stable.
- ndf15 must resolve each bin's injection window: adds O(bins × few) steps per integration —
  acceptable for a validation-grade feature; the fluid follow-up is the MCMC path.
- Very short lifetimes (Γ⁻¹ ≪ 0.01 Gyr, decay before z ~ 1e4) are out of the supported window
  (grid lower edge); guard with a warning keyed to Γ vs H(a = q_min_ratio).
