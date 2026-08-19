# Design: Inverse decays + quantum statistics for decaying warm dark matter (WDM → DR)

**Date:** 2026-07-24
**Status:** draft — written autonomously (Thomas asked for brainstorm→plan→implement with
subagents; limited feedback available). Review welcome at any point; nothing here is
irreversible.
**Physics reference:** arXiv:2011.01502 (Barenboim et al.; LaTeX source in
`arXiv-2011.01502v2/`). Prior art: Emil Brinch Holm's PhD thesis ch. 5–6
(`Emil_PhD_thesis.pdf`) and the reference branch `origin/PhD2024-EBH`.
Research extracts live in the session scratchpad (`research-thesis.md`,
`research-branch.md`, `research-plugins.md`, `formalism-notes.md`).

## 1. Goal

Extend the existing decaying-NCDM model (`ncdm_decay_dr`: massive parent WDM → massless
dark radiation, currently decay-only with a bulk integrated DR fluid) to the full
ν_H ↔ ν_l + φ system of arXiv:2011.01502:

- **momentum-resolved daughters**: a massless fermion PSD f_l(q) and a massless boson
  PSD f_φ(q), each with its own grid and Boltzmann hierarchy;
- **inverse decays** (l + φ → H), which couple different momentum regions of all three
  distributions;
- **quantum statistics** (Pauli blocking of H and l, Bose enhancement of φ).

Both at the background level (evolving f̄_i(q,τ)) and at the linear perturbation level
(collision terms in the Ψ_{i,ℓ}(q) hierarchies, synchronous gauge).

Physics motivations (thesis ch. 6): (i) the data-preferred DWDM corner is fast
*relativistic* decay, exactly where inverse decays + QS matter; (ii) the unresolved
γ³-vs-γ⁵ neutrino-lifetime transport-rate dispute can be settled by measuring the
anisotropic-stress suppression during the decay/inverse-decay equilibrium phase — Emil's
explicit recommendation is that this is answerable with QS off, i.e. before the hardest
regime is even needed.

## 2. What killed the previous attempt (2020–2024), distilled

1. **Discrete non-conservation.** The three reduced 1-D collision integrals were
   discretized *independently* (per-species trapezoid grids, off-grid f by linear
   interpolation). The continuum conservation proofs rely on reordering a double
   integral under kinematic Θ-functions (paper App. C); independent quadratures break
   that pairing, so number and energy leak, accumulating over ~10⁴ steps.
2. **Moving-endpoint discontinuities.** Kinematic limits (ε₁±q₁)/2, |a²m²/4q − q| sweep
   across a *static* grid; splicing them into a trapezoid makes the RHS jump each time an
   endpoint crosses a node. The last 15 branch commits chased exactly this.
3. **The calibration band-aid.** Rescaling daughter derivatives by A+B·q (3×3 solve per
   step) restored conservation *post hoc* but visibly deleted physics (the lasing
   fermion-energy spike) and has a degenerate branch that zeroes the parent RHS.
4. **Neutrino lasing non-convergence.** With QS on, the Bose runaway
   ∂_τ f_φ ∝ +f_φ ∫(f_H − f_l) populates the *lowest sampled* boson bin fastest; with
   leaking conservation the runaway is unbounded and the result never converges in
   q_min. (Barenboim et al. avoided it only by using q_min ~ 0.1.)
5. Stiffness at equilibrium onset/exit; ndf15's numjac fought the moving RHS structure;
   rkdp45 was ported as a workaround.

## 3. Approaches considered

**A. Re-port Emil's scheme into the plugin architecture** (independent quadratures +
A+B·q calibration). Fast to write, known territory — but inherits failure modes 1–4
verbatim. The thesis itself says the calibration is suspect. Rejected.

**B. Shared discrete transition network (chosen).** Discretize the *joint* transition
kernel once per RHS evaluation and derive every dec/inv/qs term in all three equations
from the same discrete transitions. Details in §4. Conservation and detailed balance
hold by construction; the moving-endpoint problem disappears because nodes move smoothly
inside the band instead of being spliced into a fixed quadrature. Slightly more code
(gather/scatter kernel), slightly diffusive in momentum (two-bin deposits) — both
acceptable, and the kernel is a small, isolated, unit-testable component.

**C. Spectral/moment expansion of the PSDs.** Poor fit for the spiky, non-thermal
distributions (mono-energetic injection lines, low-q lasing peak). Rejected.

## 4. The core numerical design: conservative transition network

### 4.1 Continuum structure being discretized

All three background equations are band integrals of one phase-space factor
(paper eq. A.30-A.32 / thesis eq. 6.2), with ε₁ = ε₂ + ε₃ enforced kinematically:

    Λ(q₁,q₂,q₃) = f̄_l(q₂) f̄_φ(q₃) [1−f̄_H(q₁)] − f̄_H(q₁)[1−f̄_l(q₂)][1+f̄_φ(q₃)]
                = f̄_l f̄_φ − f̄_H + f̄_H (f̄_l − f̄_φ)          ← cubic terms cancel exactly

    df̄_H/dτ(q₁) = +K/(ε₁q₁) ∫_{band(q₁)} dq₂ Λ,          K ≡ a² m_H² 𝔤²/4π = a² m_H Γ
    df̄_l/dτ(q₂) = −K/q₂²   ∫ dq₁ (q₁/ε₁) Λ
    df̄_φ/dτ(q₃) = −2K/q₃²  ∫ dq₁ (q₁/ε₁) Λ               (×2 = dof ratio g_H/g_φ)

with band(q₁) = [(ε₁−q₁)/2, (ε₁+q₁)/2] (massless daughters; flat in daughter energy).
Λ = 0 identically for FD/BE distributions with common T and μ_H = μ_l + μ_φ.

### 4.2 Discrete network

Grids: parent {q₁ᵢ} (existing DNCDM background grid), daughters {q₂ⱼ}, {q₃ₖ}
(static, log-spaced with explicit q_min; see §6 grids).

Per RHS evaluation, for each parent bin i:

1. Place S Gauss–Legendre nodes q₂*ₛ across band(q₁ᵢ) (S ≈ 4–8, precision parameter).
   Node positions are smooth functions of a — nothing is ever spliced.
   Partner momentum q₃*ₛ = ε₁ᵢ − q₂*ₛ exactly.
2. Gather f̄_l(q₂*ₛ) and f̄_φ(q₃*ₛ) by linear interpolation on their grids, with weights
   w. Evaluate Λₛ once.
3. Scatter the transition rate into all three derivative arrays with the *same* weights
   (gather-with-w / scatter-with-wᵀ):
   - parent bin i:  +(node weight)·Λₛ · K/(ε₁ᵢq₁ᵢ) · (parent measure)
   - daughter bins (j, j+1): two-bin split λ, 1−λ with λ q₂ⱼ + (1−λ) q₂ⱼ₊₁ = q₂*ₛ
   - sibling bins (k, k+1): same for q₃*ₛ.

Because ε = q for massless daughters, the two-bin linear split conserves *number and
energy simultaneously and exactly*; total energy per transition is exact because
q₂*ₛ + q₃*ₛ = ε₁ᵢ by construction. Working in per-bin particle-number variables
(N-weights: dN ∝ q² f dq absorbed into the deposit bookkeeping) makes both conservation
laws hold to machine precision **for any grid resolution** — this replaces the A+B·q
calibration entirely.

Consequences:
- **Detailed balance is exact on the grid**: Λₛ = 0 for discrete FD/BE equilibria
  (ε₁ = q₂* + q₃* holds exactly), so the stiff quasi-equilibrium phase has a genuine
  discrete fixed point. Residual interpolation error only redistributes *along* the
  conservation manifold (it cannot create/destroy N or E) and converges away with grid
  refinement.
- **Lasing becomes bounded and physical**: the runaway saturates when the available
  inversion ∫(f_H − f_l) is consumed — which exact conservation enforces. f_φ in the
  lowest bin may grow huge (transient condensate — that is real physics), but the
  observables ρ_φ, N_φ, and the fermion-energy spike are capped by the conserved
  budgets. Grid convergence is then judged on observables, not on f(q_min). A
  `bose_occupancy` diagnostic column will flag condensate-like behavior.
- **No endpoint discontinuities**: the RHS is continuous in τ (kinks only when a Gauss
  node crosses a grid point — benign, C⁰).

### 4.3 Perturbations on the same network

Collision terms 𝒞^{(1)}_ℓ (thesis eqs. 6.4–6.8) use the same transitions with angular
factors P_ℓ(cos α*ₛ), P_ℓ(cos β*ₛ), P_ℓ(cos γ*ₛ) evaluated by recurrence at the exact
node kinematics (branch lesson: never interpolate P_ℓ tables). Gather Ψ at nodes,
scatter with the same weights, per ℓ up to a truncation `l_max_dncdm_col`.

Exactness inherited per transition:
- ℓ=0: number & energy conservation of the perturbed operator — same argument as
  background with f̄ → f̄Ψ linearization.
- ℓ=1: momentum conservation — q⃗₁ = q⃗₂ + q⃗₃ projects to
  q₁ = q₂*cos α* + q₃*cos β* exactly at each node, and the two-bin deposits preserve
  first moments.

The parent hierarchy has *no* pure-decay collision term (it cancels identically in Ψ
variables — paper §4.3); inverse decays are what first source it. Decay-only runs must
therefore reproduce the existing dncdm_dr perturbations exactly — a strong regression
anchor.

## 5. Architecture (species-as-plugins)

New, parallel to the existing path — the current `ncdm_decay_dr` (parent + bulk DR)
remains untouched when the new flags are off. No behavior change for existing users.

- **`tools/decay_transition_kernel.{h,cpp}`** — the isolated core: owns the three grid
  views, node placement, gather/scatter weights, Λ evaluation, and the per-ℓ angular
  factors. Pure functions of (a, m, Γ, grids, f-arrays); no CLASS plumbing. Unit-tested
  standalone to machine-precision conservation (this is the component everything else
  trusts).
- **`species/dr_psd_species.{h,cpp}`** — ONE daughter class (`DrPsdSpecies`,
  kTypeName `dr_psd`), statistics = fermion|boson as a parameter (branch "sibling"
  lesson: code one daughter equation). A massless NCDM-like species: per-q background
  PSD evolved as ODEs (f directly, positivity floor — not ln f, since gains dominate),
  per-q massless hierarchy in perturbations, dlnf̄/dlnq computed live from a ln f–ln q
  spline (no frozen-PSD special cases). Standalone it is a free-streaming massless
  PSD species ("all species are equal").
- **`species/dncdm_inv_species.{h,cpp}`** — composite (parent `DNCDMSpecies` + fermion
  `DrPsdSpecies` + boson `DrPsdSpecies`) owning the kernel. Background coupling via
  `BackgroundDerivs` override; perturbation coupling via the existing
  `AddCouplingDerivs` seam. Follows the composite layout contract
  (generic `child_layouts`, per the #358/#359 lesson). Parent's own decay loss comes
  through the kernel too (never double-counted with DNCDM's diagonal ln f decay term —
  the parent child runs in "collision-owned" mode when inverse decays are on: f
  variables, kernel-supplied RHS).
- **Input** (dot-notation on the dncdm instance):
  - `dncdmX.inverse_decays = yes|no` (default no)
  - `dncdmX.quantum_statistics = yes|no` (default no; class_test: requires
    inverse_decays — qs terms are only meaningful with both dec+inv present)
  - daughter grid controls: `dncdmX.dr_q_min/dr_q_max/dr_N_q` (+ `_pt` variants),
    log-spaced trapezoid default; parent grid reuses existing DNCDM controls
  - `l_max_dncdm_col` precision parameter (collision truncation)
  - optional initial daughter abundances (default zero)
  - Γ input unchanged (`dncdmX.Gamma`/lifetime forms); all collision prefactors reduce
    to K = a²m_H Γ, so no separate coupling 𝔤 input is needed.
- **Guards** (error-severity conventions: class_test_severe only for structure):
  synchronous gauge only, no tensors with inverse decays on, no ncdm fluid
  approximation for the three collision-coupled species (all `class_test` at input
  time, following the type3 precedent).
- **Evolver**: background integrates all three PSDs inside the existing background ODE
  system (ndf15 default; RHS is now smooth so numjac should cope — the union sparsity
  of the kernel couplings is declared static). Perturbations: expected to prefer
  `evolver=2` (rkdp45), as for dcdm_wdm (14× lesson); both must work.
- **Speed lever (later stage, optional)**: once the parent EoS drops below a threshold,
  inverse+QS are kinematically off and the model reduces to the existing decay-only
  path; a switch-to-decay-only approximation (reusing the dncdm switch-copy machinery)
  can be added behind a flag after correctness is established. Not in v1.

## 6. Grids

- Parent: existing DNCDM grid machinery.
- Daughters: log-spaced trapezoid with explicit q_min (thesis compromise:
  bg q_min=1e-3, q_max=1e2, ~100 bins; pt q_min=1e-2, q_max=25, ~15 bins as defaults).
  Gauss–Laguerre and Kainulainen were both tried and rejected historically.
- The #384 injection-adapted quantile grid is a candidate refinement for the daughter
  grids (injection support is analytically predictable in the decay-only late phase);
  deferred until conservation-exact baseline works — with exact conservation the grid
  only affects resolution, not budgets, so this becomes a pure accuracy knob.
- Separate background/perturbation grids with a log-spline bridge (branch lesson 7),
  `N_bg ≥ N_pt` enforced.

## 7. Validation ladder (in order; each step gates the next)

1. **Kernel unit tests** (new test executable, registered in BOTH CMakeLists.txt and
   Makefile TEST_TARGETS): machine-precision number/energy conservation for random
   PSDs and random (a, m, Γ); Λ ≡ 0 on discrete FD/BE equilibria (detailed balance);
   ℓ=1 momentum identity; RHS continuity as band edges cross grid nodes.
2. **Decay-only regression**: `inverse_decays=no` runs must match current master
   `ncdm_decay_dr` background AND C_ℓ/P(k) to ~0.1% (never bit-identical; TE
   zero-crossing handling per repo conventions). Also: `inverse_decays=yes` with tiny Γ
   → ΛCDM+ν limit.
3. **Background physics**: dec+inv reaches quasi-equilibrium; final PSDs fit
   1/(exp((q−μ)/T)±1) with μ_H = μ_l + μ_φ; non-relativistic limit shuts inv+qs off
   automatically; comoving sector energy redshifts exactly as radiation during
   equilibrium; comparison against Barenboim et al. fig. 2 morphology (dec / dec+inv /
   dec+inv+qs energy-density panels).
4. **Lasing run** (m=0.5 eV, large Γ, qs on): observables (ρ_i, N_i, fermion-energy
   spike) must be bounded and q_min-convergent even as f_φ(q_min) grows; conservation
   columns stay at machine precision.
5. **Perturbation identities**: discrete ℓ=0/ℓ=1 conservation of the collision
   operator; parent δ,θ,σ converge to daughters' during equilibrium; decay-only
   perturbations identical to master.
6. **The physics question (stretch goal)**: measure a⁴Π_{νφ} during equilibrium
   (thesis eq. 6.19) with qs off — does anisotropic stress get suppressed, and with
   which effective γ-power? This is deliverable-grade output for the γ³/γ⁵ dispute.

## 8. Out of scope (v1)

Massive daughters (m_l > 0) — kinematics kept general in the kernel API but only the
massless path implemented/tested. Newtonian gauge, tensors, fluid approximation for the
coupled trio, DCDM-approximation speed switch, 2↔2 scattering channels, MCMC/emulator
work. Multiple simultaneous inverse-decay families (allowed by layout, untested).
