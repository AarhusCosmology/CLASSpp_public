# `dcdm_wdm` injection-adapted momentum grid — design

**Date:** 2026-07-18
**Species:** `dcdm_wdm` (dark matter → two massive daughters, arXiv:2606.14849; PR #371).
**Related:** `docs/superpowers/specs/2026-07-15-dcdm-wdm-massive-decay-products-design.md` (original
feature), `notebooks/dcdm-wdm-massive-decay-products.ipynb` §5 (convergence).
**Scope:** the daughter's momentum-grid *placement and spacing* only. No change to the injection
sum-rule, the perturbation hierarchy, the stress-energy, or the composite/shooting machinery.

> **Status (2026-07-23, v3 BUILT).** The erf deposit described below (§4.3) is landed and
> verified: `ff014f78` (smooth erf energy deposit replaces CIC; `ndf15` guard removed),
> `5addd728` (static-sparsity floor on the injection source — fixes an `ndf15` numjac
> pattern-lock blowup the erf deposit alone did not close), `1586d250` (completes `∂J/∂ln q`
> with the `-(3+q²/ε²)J` measure term — fixes a flat ~−7.6% low-k `P_m` deficit the one-term
> derivative silently produced). Verification, in order:
> - **Hang-regression gate** (the exact sweep that spun 8h on the withdrawn CIC deposit,
>   `.superpowers/sdd/task-2-report.md`): **8/8** green on `evolver=2` (RKDP45, 1.6–6.6 s each)
>   and, after the sparsity-floor fix, **6/6** green on the default `evolver=1` (ndf15,
>   13.9–81.7 s each) — both evolvers now complete every previously-hanging configuration.
> - **Low-k ΛCDM anchor** (`task-3-report.md` §5): `P_m(k)` at `k≤0.01 h/Mpc` matches a
>   ΛCDM control to **−0.04%** (`R−1`, modelA, Γ⁻¹=1 Gyr/v=0.011) / −0.07% (modelB) / −0.01%
>   (v→0 CDM-limit null test) — all inside the `<1%`/`<0.3%` gates; before the dJ/dlnq fix
>   these were ≈−7.6% (see `1586d250`'s root cause).
> - **Old-vs-new convergence comparison** (redone on this fixed build; the earlier `311337bc`
>   numbers were stale, generated before `1586d250`): truth = this branch's grid at 192 bins,
>   independently validated by the two checks above. At a 1% `P_m(k)` target the
>   injection-adapted grid needs **16 bins** vs **70–73** for master's uniform-`ln q` grid
>   (**4.4–4.6×** fewer); old-192 vs new-192 max`|ΔP_m/P_m|` is `1.3–1.7e-3`, concentrated away
>   from low k (`~1.1e-4` at `k=0.01 h/Mpc`) — consistent with master's flat **+1.9%**
>   `⟨β²⟩` deposit-placement bias (`task-3-diagnosis.md`) rather than a shared-continuum-limit
>   problem.
>
> Historical findings that remain load-bearing: **edges clamp, they do not drop** (the drop failed
> the over-time golden, 6.7% at z=9 — `err ≈ q_edge_tol/F(z)`), the **over-time energy golden is
> the decisive test** (`< 5e−3` at z ∈ {9, 3, 0}), and the erf deposit's width
> `σ = clamp(Δu_loc, 0.03, 0.25)` plus the `1e−12` static-sparsity floor are both required for
> `ndf15` (§4.3/§4.4) — the withdrawn order-1 CIC deposit (`93b9bf00`..`2340f3f3`) banned
> `ndf15` outright and, even under RKDP45, aliased against the background table at
> `momenta_bins ≳ 48` (spline-interpolated `index_bg_inj_` columns ringing, the perturbation RHS
> turning noisy, adaptive step control collapsing to micro-steps — a crawl indistinguishable
> from a hang).

## 1. Problem

The WDM daughter's phase-space distribution is *built up* by decay: a decay at scale factor `a′`
injects two daughters at comoving momentum `q = a′·q_kick` (`q_kick = 10` by convention), each
carrying energy `m_χ/2`. Nothing is injected above `q_kick`, so `f₀` has **compact support**
`[0, q_kick]`.

The grid (`species/wdm_decay_product.cpp:123-132`) is **uniform in `u = ln q`** over a **fixed
4-decade** window `[q_min_ratio·q_kick, ~q_kick] = [1e-3, 10]`, ~24 bins/decade at the default
`momenta_bins = 96`. But the injection weight per `d ln q` is a **peak**, not a plateau:

    w(q) ∝ a′^{3/2} · e^{−Γ·t(a′)} |_{a′ = q/q_kick}   (matter era),

peaking at `Γ·t = 1`, i.e. `t = 1/Γ` (the lifetime), with a shallow `q^{3/2}` rise below and an
exponential cutoff above. Which part of the fixed window is populated slides with `Γ`: for
`Γ⁻¹ = 0.1 Gyr` the peak sits at `q ≈ 0.3` and the whole top decade is empty; for `10 Gyr` the peak
is at `q ≈ 7` and the bottom is empty. Either way most of the uniform grid's bins do little work.
Notebook §5: 96 bins give only 0.4% convergence vs 192 (`|P/P₁₉₂−1|`), roughly first-order.

## 2. Goal

Replace the uniform-in-`ln q` grid with a **non-uniform, injection-adapted grid** whose bins carry
an **equal number of injected decays** — i.e. bins concentrated where `f₀` actually lives, none
wasted on empty momentum ranges. Same deterministic `momenta_bins` count → substantially better
convergence (to be quantified by the notebook re-run); users bank the saving by lowering
`momenta_bins`. This is the "smarter sampling" idea the earlier uniform-edge sketch was a weaker
approximation of.

## 3. Key insight — equidistribute the decays

Decay is a Poisson process in cosmic time: `dN_decay/dt ∝ e^{−Γt}`. With `g ≡ Γ·t`,

    dN/dg ∝ e^{−g}   (exact, any era)   ⟹   cumulative fraction  F(g) = 1 − e^{−g}.

`F` is the **decayed fraction**. Sampling uniformly in `F` therefore puts an **equal number of
decays — equivalently equal injected energy** (every decay injects `m_χ`) — in every bin. For the
non-relativistic daughters comoving number is conserved after injection, so equal-number-in-`F`
also means the grid is matched to today's `f₀(q)`. Two clean properties:

- **The equidistribution is cosmology-free** (pure Poisson-in-time). The cosmology enters *only*
  through the map `g ↦ a′` that places each decay-time quantile at a momentum. Evaluating that map
  in a **fixed fiducial `H(t)`** makes the grid a function of **`Γ` alone** — invariant across
  Ω_m, ω_b, H₀, … at fixed `Γ`, so grid-discretization error is a smooth function of `Γ` only
  (clean MCMC / finite-difference derivatives), and no `pba`/shooting state is needed at construction.
- **It auto-centers on `q_peak`** (`Γt = 1`) and **subsumes both edges** — the window is just
  `F ∈ [F_lo, F_hi]` — with no separate "where's the peak / where to cut" logic.

This is the spirit of thermal-ncdm quadrature (place nodes by the physical weight) but **structurally
different**: `get_qsampling` builds a *Gaussian quadrature* (nodes+weights for polynomial exactness
against the Fermi-Dirac weight). The daughter grid is a set of **evolved Boltzmann-hierarchy bins**
with midpoint cell widths `dq_i` and a moving injection kernel — it needs grid *points*, not Gauss
nodes. So we build the quantile grid analytically rather than reuse `get_qsampling`.

## 4. Design

### 4.1 Fiducial `H(t)`

Fixed placement-only ΛCDM constants (documented as *not* physics):

    Ω_m,fid = 0.31,   Ω_r,fid ≈ 9.2e−5  (T_cmb = 2.7255, N_eff = 3.044),   h_fid = 0.674,
    H₀,fid = h_fid · 1e5 / _c_   (code units, 1/Mpc, matching Gamma_).

Exact radiation+matter cosmic time (Λ negligible at the relevant epochs), verified against both
limits:

    t(a) = (2 / (3 H₀ Ω_m²)) · [ √(Ω_r + Ω_m a)·(Ω_m a − 2Ω_r) + 2 Ω_r^{3/2} ]

monotone in `a`, so `a′(g) = t⁻¹(g/Γ)` by a bracketed Newton/bisection (limiting forms
`a ≈ [(3/2)(g/Γ)H₀√Ω_m]^{2/3}` matter, `a ≈ [2(g/Γ)H₀√Ω_r]^{1/2}` radiation as seeds; the
rad/matter crossover is at `a_eq = Ω_r/Ω_m ≈ 3e−4`).

### 4.2 Analytic quantile grid

Let `N = momenta_bins`. Energy budget `q_edge_tol` (default `1e−3`) sets the trimmed span:

    F_lo = q_edge_tol,
    F_hi = min( 1 − q_edge_tol,  F(a′=1) ),        F(a′=1) = 1 − e^{−Γ·t_fid(1)}   (today's cap)

Cell **edges** uniform in `F`, midpoints at half-cells:

    F_edge[j]  = F_lo + j·(F_hi − F_lo)/N,          j = 0..N
    F_mid[i]   = ½(F_edge[i] + F_edge[i+1]),         i = 0..N−1
    g = −ln(1 − F);   a′ = t⁻¹(g/Γ);   q = a′·q_kick

giving grid points `q_[i] = q(F_mid[i])` and midpoint weights
`dq_[i] = q(F_edge[i+1]) − q(F_edge[i])`. `u_[i] = ln q_[i]` and per-bin `Δu_[i] = u_edge[i+1] −
u_edge[i]` (the local `ln q` spacing) are stored for the kernel. `q_bg_ = q_` (background =
perturbation grid, unchanged coupling). The exact energy sum-rule (`D = Σ dq_j q_j² ε_j G_j`,
`wdm_decay_product.cpp:229`) already tolerates arbitrary `dq_j`, so it holds verbatim.

`F(a′=1)` cap means: **fast decays** (F(1)≈1) → top quantile at `1−q_edge_tol` → `q_hi < q_kick`
(top trimmed; the §4.3 deposit clamps the off-grid upper tail into the top bin — retained, not dropped); **slow decays** (F(1) < 1−q_edge_tol) → top
quantile at `a′=1` → `q_hi = q_kick` (no top trim). Symmetric on the low side by `F_lo`.

### 4.3 Conservative injection deposit — cell-integrated Gaussian (erf) [v3, current]

> Replaces the order-1 CIC of v2 (kept below for history). Design forces, in order:
> (1) **exact energy conservation** — `Σ Wᵢ = 1` identically, on and off the grid;
> (2) **smooth (`C^∞`) source in time** — both `ndf15` and RKDP45 integrate it, and the tabulated
>     `index_bg_inj_` columns are representable by the background table
>     (`back_integration_stepsize = 7e−3` in `ln a`) at any `momenta_bins`;
> (3) **bounded placement bias on any grid** — the original local-`σ` Gaussian failed the golden at
>     39% because `σ` followed the bin width onto the wide low-`F` tail (`σ` up to ~1.5).

**Width.** `σ(u_cut) = clamp(Δu_loc(u_cut), kSigmaMin = 0.03, kSigmaMax = 0.25)` where `Δu_loc` is
the **analytic** local quantile-bin width (smooth, no table lookup):

    dF/du = Γ e^{−Γ t_fid(a)} / H_fid(a),   a = e^{u_cut}/kQKick
    Δu_loc(u_cut) = [(F_hi − F_lo)/N] · H_fid(a) e^{+Γ t_fid(a)} / Γ

with `H_fid(a) = H0_fid √(Ω_m a⁻³ + Ω_r a⁻⁴)` (same fiducial constants as §4.1). Evaluate in log
form and clamp *before* exponentiating (`Γ t_fid` can exceed 700 for fast decays; never materialize
`inf` under `-ffast-math`). The floor `0.03` ≈ 4 background-table samples per `σ` (representability);
the cap `0.25` bounds the relativistic-case placement bias (`~3σ² ≲ 0.19` in `u`, and only `O(v²)`
in energy for the non-relativistic daughters). In the peak region `Δu_loc ∝ 1/N`, so the deposit
sharpens with resolution until the floor — convergence with `N` is preserved where the mass lives.
The min/max clamp corners are two isolated `C⁰` kinks in `σ(t)`; acceptable (regime-switch-like);
smooth them only if an evolver measurably stumbles there.

**Deposit.** Let `Φ` be the standard normal CDF, `c_j = Φ((u_edge[j] − u_cut)/σ)` at the `N+1` cell
edges. Energy weights:

    Wᵢ      = c_{i+1} − c_i                      (interior cells)
    W₀     += c_0                                 (below-grid mass clamped into the bottom bin)
    W_{N−1}+= 1 − c_N                             (above-grid mass clamped into the top bin)
    ⇒  Σᵢ Wᵢ = 1   identically, for every u_cut (telescoping) — no D-renormalization

    A   = a·Γ·ρ_dcdm · a⁴ / factor
    Jᵢ  = A · Wᵢ / (dqᵢ qᵢ² εᵢ)     ⇒   Σᵢ dqᵢ qᵢ² εᵢ Jᵢ = A   (sum-rule EXACT)

All weights in `[0,1]`; `f` stays physical. The energy centroid is `u_cut` up to `O(Δu²)` cell
effects where cells are resolved (`Δu ≲ σ`) and up to sub-bin placement (`≤ Δu/2`, unavoidable for
any one-point-per-bin scheme) inside wide tail cells — in both regimes the daughters' non-relativistic
`ε ≈ a·M` makes the energy error `O(v²)`. The **decisive** correctness check remains the over-time
golden (§6.2), which this deposit must pass at `< 5e−3`.

**`∂J/∂ln q` (perturbations only) — CORRECTED 2026-07-23 (low-k field diagnosis).** The deposit's
continuum source as a function of the grid coordinate is `J(u) = A·φ_σ(u−u_cut)/(q³ε)` (the
energy-fraction deposit divides by the measure), so its full derivative has TWO terms:

    ∂J/∂ln q = [A·φ′_σ/(q³ε)]  −  (3 + q²/ε²)·J

Discretely, with the normal pdf `p(x) = e^{−x²/(2σ²)}/(σ√(2π))` at the cell edges:

    dJdlnqᵢ = A · [p(u_edge[i+1] − u_cut) − p(u_edge[i] − u_cut)] / (dqᵢ qᵢ² εᵢ)  −  (3 + qᵢ²/εᵢ²)·Jᵢ

The first (edge-difference) term telescopes to `≈ 0` under the `Σ dqᵢ qᵢ² εᵢ (·)` moment; the
second supplies the integration-by-parts value `−Σ dqᵢ qᵢ² εᵢ (3 + q²/ε²) Jᵢ ≈ −3A` (non-rel)
that the continuum demands. **Omitting the second term kills the daughter's metric-driven growth**
(the `g`-channel's ρ-weighted moment drives the fluid-limit continuity equation): measured effect
= a flat **−7.6%** P_m deficit at k ≤ 0.01 h/Mpc vs a ΛCDM control, surviving `v_kick → 0` and
bin-count changes (see `.superpowers/sdd/task-3-lowk-diagnosis.md`). Master's old kernel never had
this bug — its D-normalized `J = A·G(q)/D` has no `1/(q³ε)` factor, so its pointwise
`J·(u_cut−u)/σ²` is already the full derivative (the same `q³`-tilt that biases its momentum
placement supplies its `g`-moment correctly; the two builds' defects are complementary).
No edge special-casing; exponentials clamped master-style (`x²/(2σ²) ≥ 60 → 0`). The decisive
unit-level check is the **moment identity** `Σ dq q²ε·dJdlnq = −Σ dq q²ε·(3+q²/ε²)·J` for interior
cutoffs; the decisive field check is `P_m(model)/P_m(ΛCDM control) ≈ 1` at k ≤ 0.01.

**Evolvers.** The source is `C^∞` in time (up to the two isolated `σ`-clamp kinks): **both `ndf15`
(default) and RKDP45 integrate it**. The constructor guard rejecting `evolver=1`
(commit `5539964c`) is **removed**, along with its unit test and the test-fixture `evolver=2`
forcing. RKDP45 remains the faster choice for production (~14×, non-stiff hierarchy) — a
recommendation for the notebook, not a constraint in code.

**Static-sparsity floor (added 2026-07-23 after the Task-2 field gate).** `ndf15`'s numjac derives
the Jacobian sparsity pattern from **exact zeros** in a dense Jacobian
(`tools/evolver_ndf15.cpp`: `fabs(dFdy) != 0`), locks it permanently after `trust_sparse`
consecutive repeats, and thereafter attributes grouped finite-difference increments through the
frozen pattern. A *moving compact* source violates that static-sparsity contract: early on the
deposit is clamped into bin 0 only, a too-narrow pattern gets locked, and when the footprint later
sweeps into the grid the grouped Jacobian is corrupted — background `rho_dcdm` blows up
(`rho_crit <= 0` within 0.2 s; field gate, bins 32/64/96, both test models; master's Gaussian
kernel has the **same pre-existing defect** at bins=32). Fix (species-side, contract-conforming):
floor the deposit weights, renormalized so `Σ W = 1` stays exact —
`Wᵢ ← (Wᵢ + f)/(1 + N·f)`, `f = kSparsityFloor = 1e−12` — and seed `∂J/∂ln q` with the same
magnitude, so every injection row is *structurally* coupled to `ρ_dcdm` at all times at a
physically invisible level (~1e−12 of the instantaneous deposit). This also keeps the tabulated
`index_bg_inj_` columns structurally nonzero, protecting the *perturbation* system's
`δ_dcdm`-column pattern in `AddCouplingDerivs` from the same lock. An evolver-side fix (periodic
pattern re-verification) was considered and rejected here: it touches the shared, perf-tuned
default evolver for all species — worth a separate issue, not a side effect of this branch.

### 4.3-v2 (historical) Conservative energy-CIC injection deposit — superseded by v3

> **Revision (2026-07-19).** The first design used a **Gaussian kernel of local width**
> `σ = kernel_width·Δu_local` plus smooth onset/shutoff gates. Implementation exposed a
> discretization defect the design missed: covariant sector energy conservation is exact for *any*
> grid (the injection **rate** telescopes — see below), but the golden also checks that `ρ_wdm(a)`
> matches the ideal **monochromatic-kick** history, and a finite-width Gaussian *smears* the kick in
> `u = ln q`. Because the energy measure `dq·q²·ε ∝ q³` (relativistic) weights the symmetric-in-`u`
> Gaussian asymmetrically, the deposited energy's centroid shifts by `≈ 3σ²` — a **systematic
> O(σ²) momentum bias**. Daughters land at the wrong `q` and redshift wrongly. The old uniform grid
> had `σ ≈ 0.1` so the bias sat under the `5e−3` tolerance; the quantile grid's low-`F` (early-time)
> bins are up to ~15× wider (`Δu` to ~1.5), blowing `σ` up to ~0.3–1.5 and the z=9 energy golden to
> **39%**. The fix replaces the smearing kernel with a **conservative order-1 (linear) CIC deposit**
> whose energy centroid is exact by construction. (The `C⁰` linear source **requires** an explicit
> evolver — `evolver=2`/RKDP45; the stiff `ndf15` default hangs on its kinks, see §4.4/§8.)

> **Revision (2026-07-23).** Two more findings while landing (see top-of-doc Status): the edge deposit
> **clamps** the off-grid share into the edge bin instead of dropping it (the drop failed the
> over-time golden at high z), and the run must use **RKDP45** (`ndf15` hangs on the `C⁰` source). The
> interior CIC and the exact in-grid sum-rule are unchanged.

**Deposit.** The kick is monochromatic at `u_cut = ln(a·q_kick)`. Deposit **energy fractions** onto
the grid with a **linear (order-1) CIC** stencil `Wᵢ` in `u = ln q`, using the **actual bracketing
bin centres** so it is exactly partition-of-unity and first-moment-exact on any grid — `u_cut`
between `u_k, u_{k+1}` → `W_k = (u_{k+1}−u_cut)/(u_{k+1}−u_k)`, `W_{k+1} = (u_cut−u_k)/(u_{k+1}−u_k)`,
all other `Wᵢ = 0`:

    Σᵢ Wᵢ = 1        Σᵢ Wᵢ (uᵢ − u_cut) = 0        ⇒  energy centroid at u_cut, exactly
    A   = a·Γ·ρ_dcdm · a⁴ / factor
    Jᵢ  = A · Wᵢ / (dqᵢ qᵢ² εᵢ)     ⇒   Σᵢ dqᵢ qᵢ² εᵢ Jᵢ = A   (sum-rule EXACT, no D-normalization)

The exact sum-rule is automatic (`Σᵢ dqᵢ qᵢ² εᵢ Jᵢ = A·Σᵢ Wᵢ = A`); the two nonzero weights are
always in `[0,1]` (`f` stays physical). First-moment exactness kills the systematic bias; the only
residual is the **non-systematic** sub-bin variance `O(Δu²)`, which does not accumulate coherently.
(First-moment exactness is an **interior** property; at the grid edges the deposit clamps into the
one edge bin — `Σ Wᵢ = 1` is preserved, i.e. energy is conserved, but the centroid is the edge bin,
not `u_cut`. Harmless here: the daughters are non-relativistic so `ε ≈ a·M` barely depends on `q`.)

- **Why order-1 only (higher order deferred).** Order-1 is the *uniquely robust* choice on this grid:
  it interpolates between the actual bracketing bins, so it stays first-moment-exact even on the wide
  low-`F` tail (`Δu` up to ~1.5). Higher B-spline orders (TSC/cubic) would smooth the source, but on a
  wide bin a *moment-corrected* stencil is non-smooth (the correction jumps as the stencil gains/loses
  a bin) and an *uncorrected* one has `O(Δuᵖ)` first-moment error that diverges for `Δu>1` (bias
  returns). Higher order would require first bounding the bin width (`Δu_max` cap, blending the
  F-quantile with uniform-`ln q`) so no bin is wide.

  **Re-evaluated 2026-07-23.** One historical blocker for higher order — *edge* energy conservation
  (a wide stencil hanging partly off-grid breaks partition-of-unity) — is now **resolved in
  principle**: the same clamp/renormalize we adopted at order 1 (renormalize the in-grid stencil
  weights to sum to 1) generalizes to any order. And a smoother source (order-2 → `C¹` `J`, `C⁰`
  `∂J/∂ln q`) would very likely let even `ndf15` integrate without hanging. **We still defer it**,
  deliberately: (a) `ndf15` buys nothing here — the daughter hierarchy is non-stiff, so RKDP is
  already correct *and* ~14× faster; smoothing the source only unlocks the slower path; (b) the
  "user leaves the default `ndf15` and silently hangs" robustness gap is better closed by a small
  evolver guard (§4.4) than by a deposit rewrite; (c) higher order buys *smoothness* but not
  *accuracy* in this regime — the daughters are non-relativistic (`ε ≈ a·M`, q-independent), so
  first-moment placement barely affects energy (this is exactly why clamping into one edge bin still
  conserves to `0.4%`). The interior wide-bin tension (smooth XOR first-moment-exact) is untouched by
  the clamp and only *bites* once daughters are relativistic (`ε ≈ q`, large `v_kick`). **That
  relativistic-daughter case is the real trigger** to revisit — and there one wants a `Δu_max` cap
  **and** higher-order together. Until then, YAGNI. `kernel_width` is removed with no replacement
  input.
- **Edges — clamp, don't drop (revised 2026-07-23).** When `u_cut` falls below the lowest bin centre
  (or above the highest) the whole off-grid share is **clamped into that edge bin** (`W_edge = 1`), so
  `Σᵢ Wᵢ = 1` **always** and no injected energy is lost. This *replaces* both the onset gate and the
  high-side shutoff. The earlier design instead *dropped* the off-grid share (one-sided ramp,
  `Σᵢ Wᵢ < 1`, bounded by `q_edge_tol`); that removed the earliest `q_edge_tol` of decays, which is a
  large fraction of the small cumulative decayed-fraction at high z, and **failed the over-time golden
  (6.7% at z=9)**. Clamping conserves energy from the first decay (golden `3.7e−3` at z=9). The deposit
  stays `C⁰` (`W_edge` is flat at 1 past the edge, then ramps down entering the interior). It is
  polynomial — **no `exp()` in the kernel**, so the old `-ffast-math`/exp-clamp-60 concern is gone.
- **`∂J/∂ln q` (perturbations only).** Sources the per-bin `g ≡ ∂f/∂ln q` hierarchy variable
  (background injection passes `nullptr`). Use the analytic `∂/∂ln q` of the linear-CIC profile
  (a `±1/Δu` step across the two stencil bins), keeping `g` consistent with `∂f/∂ln q`. Verified
  finite/well-behaved by the perturbed smoke test; a smoother source is the deferred-higher-order
  lever if needed.

### 4.4 Guards, inputs, backward-compat

- **Negligible-decay fallback:** if `F(a′=1) ≤ 2·q_edge_tol` (span collapses; daughter essentially
  unpopulated), fall back to a minimal valid uniform-`ln q` grid over `[q_kick/30, q_kick]` so the
  hierarchy/kernel stay well-formed. Documented; observationally irrelevant (empty species).
- **Count validation:** `momenta_bins ∈ [8,4096]` (the linear deposit's 2-bin stencil imposes no
  further floor). The old `> 3·kernel_width + 4` check is removed with `kernel_width`.
- **`kernel_width`** — **removed**, no replacement input (the order-1 deposit has no width knob).
  Higher CIC order is deferred (§4.3) and would reintroduce an input then.
- **`q_edge_tol`** — new input, default `1e−3`, range `[1e−12, 1e−2]`. It sets the trimmed **grid
  F-span** (which decay-quantiles get their own bins). With the clamp (§4.3) the trimmed tail is
  **lumped into the edge bin, not dropped**, so it is no longer an energy-budget knob — it is a
  grid-resolution knob (how finely the extreme-`F` decays are resolved in `q`); total energy is
  conserved for any value.
- **Evolver — no restriction (v3).** The erf deposit is smooth, so **both evolvers work**; the
  `5539964c` constructor guard, its unit test, and the test-fixture `evolver=2` forcing are removed.
  The over-time golden must pass under the **default `ndf15`** (this is the hang-regression test);
  RKDP45 stays the documented fast path (~14×) for production runs.
- **`q_min_ratio`** — retained as an **optional hard floor** on `q_lo` only (`q_lo ≥
  q_min_ratio·q_kick`); default unset. `max_q`/uniform-grid semantics are gone (the grid is
  injection-adapted). This is a research feature validated by the notebook, not a pinned API.
- The original §6 **short-lifetime limitation dissolves** — the window follows the decay.

## 5. Files touched

- `species/wdm_decay_product.cpp` / `.h` (v3 delta over the landed CIC state): replace
  `CicWeights` with the §4.3 erf weights (+ the analytic `Δu_loc`/`σ` helper); replace the
  `FillInjection` weight call and the three hand-written `dJdlnq` edge branches with the
  edge-pdf-difference form; delete the `ndf15` constructor guard. The grid build
  (`BuildInjectionAdaptedGrid`), inputs (`q_edge_tol`, `q_min_ratio` floor), and everything
  downstream (`ComputeBackground`, hierarchy, stress-energy) are unchanged.
- `species/dcdm_wdm_test.cpp` (v3 delta): drop the evolver-guard test and the `evolver=2` fixture
  forcing; replace the CIC exact-centroid (1e−9) assert with the σ-tolerance interior-centroid
  check; keep sum-rule exactness, edge clamps, `∂J/∂ln q` sign structure, grid tests, and the
  decisive over-time golden — now run under the **default `ndf15`** (§6.2).

`cclassy.pxd` is auto-generated — never hand-edit. `q_edge_tol` flows via the existing
`SpeciesInput`/dot-syntax path; no wrapper edits.

## 6. Testing & validation

Unit/golden (`dcdm_wdm_test`, fast, CI — registered in **both** `CMakeLists.txt` and Makefile
`TEST_TARGETS`):

1. **Equidistribution (construction).** The per-bin decayed-fraction width `ΔF_i` is uniform to
   machine precision and `q_[i]` equal the analytic quantile map `q(F_mid[i])`; `q_[i]` monotone;
   grid spans `[q_lo, q_hi] ⊆ (0, q_kick]`.
2. **Energy budget + over-time conservation (decisive).** With the edge clamp (§4.3) the realized
   sum-rule fraction `Σ Wᵢ = 1` at **every** `a` (below, inside, and above the grid) — no injected
   energy is dropped — and the instantaneous sum-rule `Σ dqᵢ qᵢ² εᵢ Jᵢ = A` holds exactly. The
   **decisive** check is the over-time golden: `ρ_wdm(a)` at z ∈ {9, 3, 0} matches the analytic
   monochromatic-kick injection integral to `< 5e−3`. **v2 CIC result:** `3.7e−3 / 1.0e−3 / 2.7e−4`
   at z=9/3/0 (the local-σ Gaussian failed at 39%; the pre-clamp order-1 *drop* failed at 6.7% at
   z=9). The v3 erf deposit must re-pass the same bar — under the default `ndf15` (§6.6).
   > **Note (2026-07-23).** The spec originally proposed a `Δu_max` bin cap as the fallback if the
   > golden leaked. It leaked — but the cause was **not** the interior placement residual the cap
   > addresses; it was the **edge drop** (`Σ Wᵢ = 0` below the grid at high z; `Σ Wᵢ = 1` and zero
   > energy loss for every *in-grid* cutoff). A `Δu_max` cap would not have fixed it. The clamp did.
   > The `Δu_max` cap remains a (still-unused) lever only for the *interior* residual, relevant to a
   > future relativistic-daughter case (§4.3 re-eval, §8).
3. **Cosmology invariance.** Two runs differing only in Ω_m (or h) produce **identical** `q_[i]`
   (machine precision) at fixed `Γ` — the defining property of the fiducial map.
4. **Edge placement across lifetimes.** For `Γ⁻¹ ∈ {0.1, 1, 10} Gyr`, `q_peak` lands at `Γt_fid=1`,
   `q_hi = q_kick` for the slow case (no top trim) and `< q_kick` for the fast case (top trimmed; the
   deposit clamps the off-grid weight into the top bin — retained, not dropped).
5. **Deposit correctness (v3).** Energy sum-rule exact (`Σ dqᵢ qᵢ² εᵢ Jᵢ = A` to 1e−12 — holds for
   the erf weights identically); all weights in `[0,1]`; past either edge the whole share is
   **clamped into the edge bin** (`Σ Wᵢ = 1`, energy retained — verified at, half-a-bin below, and
   well below/above the edge); **interior** energy centroid `|Σ (dqᵢ qᵢ² εᵢ Jᵢ) uᵢ / A − u_cut| <
   0.5·σ` at probes in the well-resolved mid-`F` region (the CIC's 1e−9 exactness is a property of
   the linear stencil, not a physical requirement — the golden is the physical requirement);
   `J` peaks at the cutoff bin, `∂J/∂ln q ≥ 0` below the cutoff and `≤ 0` above. The **perturbed
   smoke test** confirms the smooth injection source drives finite, sensible `Cₗ`/`P(k)`.
6. **ndf15 hang-regression (v3, new).** The full background golden run (§6.2) executes under the
   **default evolver (no `evolver` key)** and completes — the constructor accepts `evolver=1`
   again, and the smooth source must not thrash the stiff step control.
7. Existing construction / gauge-guard / composite-layout / `ε→0` / `vkick→0` tests pass unchanged.

Integration, comparison harness & notebook:
- Runs without a `dcdm_wdm` species are byte-for-byte unaffected (`COMPARE_OUTPUT_REF` / classyref).
  `dcdm_wdm` runs shift (better-placed grid) — verify with scale-relative metrics at ~0.1%, Cl^TE
  zero-crossing handling; **never** blind max-rel-diff or a bit-identical requirement.
- **Old-vs-new convergence comparison (the previously-failed deliverable).** Old grid = a master
  build installed as module `classyref` (from the existing `.worktrees/master-oldgrid` worktree,
  `pip install . --config-settings=cmake.define.CLASS_PYTHON_MODULE_NAME=classyref` with the
  project-name switch — the established reference workflow). New grid = the branch `classy`. The
  driver (`notebooks/dcdm_wdm_convergence.py`) must be **hang-proof**: each (model × bins) run in a
  separate subprocess with a hard timeout (~30 min), results appended to the output `.npz` after
  *every* run (a hang or crash loses one point, not the sweep), timestamped progress to a log file.
  Sweep `bins ∈ {16, 24, 32, 48, 96, 192}`, reference = 192, two models (`Γ⁻¹ = 0.1` Gyr `v=0.03`;
  `Γ⁻¹ = 1` Gyr `v=0.011`), `evolver=2` on both builds (fair and fast). Verify old-192 vs new-192
  agree (consistency of the continuum limit) before trusting the reference. Deliverables: the
  bins-for-target table (`0.01/0.003/0.001` in `max|ΔP/P|`) and the convergence figure.
- Re-run notebook §5 convergence with a low-bin sweep (e.g. 16/24/32/48/96): expect the
  injection-adapted grid to reach today's 96-bin accuracy at materially fewer bins. Regenerate
  `notebooks/dcdm_wdm_fig3_cache.npz` and the Fig. 3 panels; update §5 (grid is now injection-adapted,
  not a fixed window) and the §6 short-lifetime note (auto-followed).

Success criteria: cosmology-invariant grid at fixed Γ; energy budget honored both ends; Fig. 3
physics reproduced at fewer bins than the uniform 96; all limits and the composite/shooting paths
unchanged.

## 7. Out of scope / rejected

- **Observable-weighted quantiles** (blend the injection number with a velocity/pressure `q⁴/ε`
  weight for relativistic daughters, or the perturbation free-streaming integrand). Equal injected
  *number* is exact for the non-relativistic regime the paper studies and is the clean default;
  reweighting is a future refinement if a relativistic-daughter case demands it.
- **Reusing `get_qsampling`** — it is Gaussian quadrature against a static PSD; wrong structure for
  an evolved-bin grid with a moving kernel (§3).
- **`(0, q_kick]`-style compactifying remaps** — `f₀` has compact support with no infinite high-`q`
  tail to tame; solves a problem this species does not have.

## 8. Risks / open points

- **`t⁻¹(g/Γ)` robustness:** monotone and smooth; bracket in `[a_min, 1]`, seed with the limiting
  forms, fall back to bisection. Built `N` times once at construction — cost is irrelevant.
- **High-`q` free-streaming resolution:** equal-number bins are sparse in the exponential upper
  tail, but that tail carries negligible `f₀`, so its weight in `δρ/θ/σ` (and any phase recurrence)
  is negligible — the same reason weight-adapted thermal-ncdm grids are stable. Confirm via the
  notebook `P(k)` cutoff shape; add an observable-weighted term (§7) only if it shows under-resolution.
- **Low-`F` tail — RESOLVED, but not the way this section predicted.** The over-time golden **did**
  fail (6.7% at z=9). The cause was **not** the interior `O(Δu²)` sub-bin-variance residual named here
  — that residual loses *zero* energy (`Σ Wᵢ = 1` for every in-grid cutoff, confirmed by probe). The
  cause was the **edge drop**: below the grid `Σ Wᵢ = 0`, so the earliest `q_edge_tol` of decays was
  discarded, which is a large fraction of the small cumulative decayed-fraction at high z
  (`err ≈ q_edge_tol / F(z)`; z=0 predicted `1e−3/0.174 = 0.59%` vs measured `0.57%`). Fix = **clamp,
  don't drop** (§4.3), *not* the `Δu_max` cap. The interior residual is genuinely negligible here; the
  `Δu_max` cap stays an unused lever for a future relativistic-daughter case.
- **Source smoothness / evolver — RESOLVED by v3, and worse than v2 believed.** The v2 claim "RKDP
  handles the `C⁰` source" was falsified in the field: an `evolver=2` sweep against the CIC build
  spun 8+ hours at 100% CPU (2026-07-23, process sampled + killed). Root cause was not the evolver
  choice but **table representability**: per-bin CIC tents of width `Δu_loc < back_integration_stepsize`
  alias in the tabulated `index_bg_inj_` columns; the spline-ringing injection feeds
  `AddCouplingDerivs`, the perturbation RHS turns noisy, and any adaptive step control collapses.
  The v3 erf deposit with `σ ≥ kSigmaMin = 0.03` (~4 table samples) is smooth *and* representable —
  both evolvers work; the guard is removed. Higher-order CIC (TSC/cubic) is moot: it smooths `J` in
  *q*, but its per-bin time-support still tracks the local bin width, so the aliasing returns at
  high `N`. The `σ` floor is the correct mechanism, not stencil order.
- **Notebook baseline churn:** cached Fig. 3 and convergence numbers shift; regeneration is part of
  the deliverable, not a regression.
- **Fiducial mismatch far from ΛCDM:** bins are placed slightly sub-optimally (never *wrong* — the
  physics is exact for any grid, and the clamp conserves energy for any edge placement); the realized
  edge quantile drifts from `q_edge_tol` only by the weak `√(Ω_m,fid/Ω_m,true)` factor.
