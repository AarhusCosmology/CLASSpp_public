# A reduced perturbation model for ν_H ↔ ν_l + φ

**Status:** design, 2026-08-12; **implemented** and shipped as
`<instance>.dr_reduced_moments = N` (the "moments" rung of `dr_representation`).
At the time of writing, its one approximation (§2.5) was **unvalidated** — see M4,
which is a null result, not a supporting one. That is still the caveat to carry: the
rung is trusted as an intermediate between the exact and RTA rungs, and is what the
RTA's coefficients are fitted against, not as an independent oracle.
**Goal:** CMB and P(k) at a cost that permits parameter estimation, with the exact
background retained and every coefficient of the approximation read off that background.

---

## 0. What stays exact

The background is untouched: all three species keep their momentum grids, the
`DecayTransitionKernel` network, quantum statistics, and detailed balance. It already
costs a small fraction of a run (`evolver_background = 3` made it 22×) and it is the
thing that fixes N_eff, the decay epoch and the daughter PSDs.

Everything below concerns the *perturbations only*, and every coefficient it uses is a
functional of the exact background — never a fitted number.

---

## 1. What the data says

Thirteen measurements from `~/dncdm-harness/hpc_cmb` plus two local runs. They are what the
design is built on, so they are stated before the scheme. One of them (M4) is a null
result that removes the empirical support for the design's central approximation; it is
stated at the same weight as the rest.

### M0 — a pure fluid RTA fails; the parent must stay momentum-resolved

Fitting a k-independent 3×3 matrix `A(a)` in `dΠ_i/dτ|_coll = Σ_j A_ij Π_j` on two
k-modes and predicting the third gives held-out errors of 20–100% before the decay
epoch and >100% during it. The collision genuinely resolves the parent's momentum
profile, and the massive parent phase-mixes (free streaming is `k q/ε`, q-dependent), so
no fixed q-shape represents it. **The parent's hierarchy is not reducible.** It is also
not the cost: 16 bins × 18 multipoles = 288 variables against the daughters' 7524.

### M1 — the collision lives in ~1.5 decades of scale factor

`|dΠ/dτ|_coll` against the total `|d(a⁴Π)/dτ|`, windowed RMS in ln a:

| Γ, m | window where the ratio exceeds 10⁻³ | peak ratio |
|---|---|---|
| 1e7, 0.06 | a = 1.0e-3 … 4e-2 | ~8 |
| 1e9, 0.3 | a = 3e-4 … 1e-2 | ~60 |

Outside the window the ratio falls by ten or more orders of magnitude. Before it, the
sector sits near detailed balance and an adiabatic perturbation is nearly collisionless;
after it, the parent is extinct (ρ_H/ρ_tot ~ 1e-60).

A peak ratio ≫ 1 means the collision is large but nearly cancels against free streaming
— the tight-coupling signature, and the reason the system is stiff.

### M2 — three quarters to nineteen twentieths of the steps are outside that window

One row per evolver step in the `k_output` dumps:

| cell | k | total steps | inside window | outside |
|---|---|---|---|---|
| G1e5_m0.06 | 1e-1 | 81 732 | 10 310 | 87% |
| G1e7_m0.06 | 1e-1 | 149 122 | 18 568 | 88% |
| G1e7_m0.3 | 1e-1 | 91 284 | 5 490 | 94% |
| G1e9_m0.3 | 1e-1 | 376 499 | 147 078 | 61% |

Those steps are paid for the stiff eigenvalue −aΓ acting on a parent that no longer
exists. This corroborates the earlier "freezing at a = 0.3 is bit-identical and 4.2×"
measurement, and says the switch can be far more aggressive than a = 0.3.

### M3 — the daughters are massless, so their integrated hierarchy is *exactly* closed

For a massless species the free-streaming operator is `∂_τ δf + i k μ δf = C`: the
transport term carries no q at all. Therefore

* any q-moment `∫dq q^n δf(q, ℓ)` obeys the same closed ℓ-hierarchy under free
  streaming, for any n and any q-profile, with **no** closure error; and
* a perturbation that starts separable, `δf(q,ℓ) = Φ(q) A_ℓ(τ)`, stays separable.

So collapsing the daughters to a few q-moments is exact wherever the collision is off,
which by M1/M2 is most of the run. Only inside the window is any approximation made.

### M4 — the daughter closure is UNVALIDATED, and the obvious test has no power

This finding replaces an earlier, wrong one. It is the reason the design below is a
proposal rather than a result.

The test attempted: model the daughters' ℓ = 2 collision term as the exact
energy-weighted Legendre transfer evaluated on the parent's q-resolved δf_H, with **no
daughter momentum resolution**, one constant fitted per window on k = 1e-2 and evaluated
out-of-sample on k = 1e-1. Two local runs with `CLASS_DNCDM_PERT_PSD=2`.

* **At Γ = 1e5 it appears to work** — residuals of 0.1–3% of the total `d(a⁴Π)/dτ`, worst
  band 21% for φ. **That number is meaningless**, for two compounding reasons.
* **The weight is unidentifiable.** Five candidate weights were compared by asking whether
  `dΠ_coll / (w · δf_H)` is a smooth background quantity — sign consistency and the
  interquartile spread of ln|ratio| within a window, which needs no absolute normalisation.
  The candidates were `q²F_2` (Legendre), `q²` (no beaming), `q⁴/ε` (i.e. Π_H itself),
  `q²F_2/ε`, and a deliberately meaningless **uniform** weight. At Γ = 1e5 *and* at Γ = 1e7
  all five score the same — sign 0.73–1.00, IQR 1.0–2.4, i.e. a factor 3–11 swing in the
  ratio inside a single window — and **the nonsense weight is often the best of them**
  (sign 0.95–0.97 at Γ = 1e7 against the Legendre weight's 0.81–0.99). No fixed q-projection
  of δf_H is the daughters' source, and this family of tests cannot identify one.
* **Out of sample the model simply fails at Γ = 1e7**: the residual equals the RMS of the
  collision term itself. The in-sample correlation is strong (|r| = 0.84–0.998) but its
  **sign flips** between windows and between k-modes — two different projections of a
  rotating vector, not a physical proportionality.
* Rank note, since it is easy to over-read: over the *whole* window the parent's
  δf_H(q, ℓ=2) is 94% rank-1 at Γ = 1e5 and only 63.8% at Γ = 1e7, which is why a
  whole-window constant survives at the lower Γ. Within narrow sub-windows both are 92–98%
  rank-1. So the degeneracy is real but it is not the whole story: the parent-only model is
  wrong at *both* Γ, and at Γ = 1e5 a per-window fitted constant was simply absorbing it.

The reading: the daughters' collision term is *not* dominated by the decay source. The
inverse-decay and quantum-statistics legs are large — consistent with the measured peak
occupations f̄_φ ≈ 140 and f̄_l ≈ 1.01 — and they are exactly the legs a reduced daughter
has to approximate — the ℓ = 2 collision term is **not a functional of the parent alone**.
**The load-bearing approximation in §2 has no supporting measurement, and one contrary one.**

What this does *not* touch: M1, M2, M3, M5, M6 and the F_ℓ identities are all either
emitted columns, counts, theorems, or verified algebra. The architecture stands; its one
approximation is untested, and §4 test 2 is now the test that decides it, not a
confirmation of something already known.

### M7 — the parent is never the whole source, and is a minority late

The perturbation operator is linear in (F_H, F_l, F_φ), so the split is exact: run
`ApplyPerturbationOperator` three times with one leg live and the others zeroed. Driven on
real Γ = 1e7 states (background + δf pulled out of a CLASS run), the reconstruction closes
to 1e-16. Share of |dE_ℓ=2/dτ| on the fermion daughter:

| a | parent | fermion | boson |
|---|---|---|---|
| 2.6e-3 | 50% | 45% | 5% |
| 7.6e-3 | 42% | 54% | 4% |
| 1.3e-2 | 29% | 53% | 18% |
| 2.2e-2 | 0.5% | 91% | 9% |

Roughly half parent / half daughters early, and **daughter-dominated late** (the boson is
similar). This is the mechanism behind M4: no weighting of δf_H can reproduce a quantity
that is mostly built from δf_l and δf_φ. It also says the RTA/loss leg of §2.5 is not a
correction to be got roughly right — late in the window it *is* the answer.

### M8 — a moment/group daughter basis converges slowly on the collision rate

Same probe, with the daughters' δf replaced by its projection onto a small basis (weighted
least squares in the q³dq measure, parent kept exact, no fitted constant). Median relative
error of the ℓ=2 collision moment rate over the window, at k = 1e-2:

| basis | DOF | dE_H | dE_l | dE_φ |
|---|---|---|---|---|
| amp1 `{f̄}` — pure RTA | 1 | 9.8 | 1.8 | 15.5 |
| eq2 `{−∂f̄/∂lnq, −q∂f̄/∂lnq}` | 2 | 7.8 | 2.6 | 12.0 |
| 8 groups | 8 | 4.2 | 2.4 | 2.5 |
| 16 groups | 16 | 1.8 | 1.2 | 0.9 |
| 32 groups | 32 | 0.9 | 0.7 | 0.4 |
| 100 groups | 100 | 1.3 | 0.33 | 0.16 |
| 134 groups | 134 | 0.48 | 0.19 | 0.06 |

Roughly 1/N, not spectral, and not monotone (g50's dE_H is worse than g32's) — the residual
is fine q-structure that a group projection samples erratically. **On this metric the
two-moment design of §2.1 is dead by two orders of magnitude.**

⚠ **But this metric is not calibrated to the observables, and the one anchor available says
it is wildly pessimistic.** g100 is roughly the resolution of the r2 rung, and it scores
16–133% here — while the actual r1-vs-r2 C_ℓ^TT difference is **0.27%** (M6). So an O(1)
error in the ℓ=2 collision moment rate is worth O(0.3%) in C_ℓ. Most of the collision-rate
error evidently cancels in the observable, which is consistent with M5: the structure it is
sensitive to is injection-time phase that the metric integrates away.

**Do not design against dE/dτ either.** The decisive measurement is the observable's own
sensitivity to daughter perturbation resolution, which `dr_bg_refine` already produces with
no new code — sweep it and read C_ℓ and P(k). That sweep is the open item.

### M9 — coarsening the shipped scheme does not degrade, it DIVERGES

`dr_bg_refine` coarsens the daughter *perturbation* grid while leaving the background
untouched (verified: z_eq and the age agree to every printed digit across the sweep), so it
is exactly the experiment a reduced daughter needs. At Γ = 1e7, m = 0.06, P(k) at z = 0:

| refine | daughter points | wall time | max ΔP/P vs 51 pts |
|---|---|---|---|
| 4 | 51 | 663 s | (reference) |
| 10 | 21 | 219 s | 4.3e-2 |
| 20 | 11 | 105 s | **9.3e+10** |
| 40 | 6 | 87 s | **NaN — ETD non-finite state** |

P_cb tracks P(k) to four digits, so this is not the parent's ill-conditioned δ leaking in.
The k-dependence is the tell: the error is 1e4–5e10 for k ≤ 0.017 and falls to **1.3–2.8%**
for k ≥ 0.04. That is an unopposed exponential growth that free streaming (which enters the
hierarchy as k) outruns at high k.

### M10 — the mechanism is a positive eigenvalue at ℓ = 1, and only ℓ = 1

`CLASS_DNCDM_SPECTRUM` dumps the per-multipole collision matrix M_ℓ column by column.
Eigenvalues at a = 9.0e-3:

| daughter points | ℓ=0 | ℓ=1 | ℓ=2 | ℓ=3 |
|---|---|---|---|---|
| 51 | 0 | **+9.07e-5** | 0 | 0 |
| 11 | 0 | **+5.42e-3** | 0 | 0 |

Every other multipole is exactly dissipative (max Re λ = 0, the conserved zeros). ℓ = 1 is
the **momentum** channel, whose exact zero is the momentum conservation identity — and the
coarse deposit breaks it, turning the zero into a growing mode.

Quantitatively: the collision window spans Δτ = 2620 Mpc, so the predicted runaway is
exp(λΔτ) = 1.27 at 51 points and **1.5e6** at 11 points. P ∝ δ², and the measured
√(ΔP/P) is **3.1e5**. Same order, from a completely independent route.

The shipped 201-point grid measures λ = +4.14e-7, i.e. 0.1% growth across the window — it is
positive but harmless, which is exactly what a grid rule tuned to the stability boundary
should produce. 51 points grows 27%, 21 points grows 39×, and 11 points grows 1.5e6.

**The `balanced_gather` hypothesis is falsified at every resolution.** Enabling the log-odds
gather, which makes Λ vanish exactly at equilibrium, changes nothing: 4.350e-2 → 4.362e-2 at
21 points, 9.335e10 → 1.045e11 at 11, and both crash identically at 6. Detailed balance is
not what fails here; the ℓ=1 momentum identity is.

### M11 — the law, and it explains BOTH of the repo's historical grid rules

λ(ℓ=1) was measured at two Γ and four resolutions. Both exponents come out clean:

```
    λ  =  C · Γ^1.02 · N^(−s)          s = 2.09 (coarse) … 3.93 (fine)
```

Γ-scaling at fixed N: λ(1e7)/λ(1e5) = 103.9 at N = 11 and 117.7 at N = 51, i.e. **p = 1.01
and 1.03** — exactly linear in Γ, as the collision rate is.

The grid requirement is then a **stability** condition, λ·Δτ ≲ 1, not an accuracy one:

| local slope s | regime | implied rule N ∝ Γ^(1/s) | the repo's measured rule |
|---|---|---|---|
| 2.09 | coarse (N = 11→21) | Γ^0.48 | **Γ^0.50**, before `lump_limiter` |
| 3.09 | N = 21→51 | Γ^0.32 | |
| 3.93 | fine (N = 51→201, measured) | **Γ^0.25** | **Γ^0.25**, after `lump_limiter` |

Both empirical laws are the same eigenvalue condition read at two points on one
convergence curve. And `lump_limiter` is what moved the operating point between them: the
momentum misplacement the lumped loss books drops from 2.08e-1 to 5.85e-2 at N = 51 when it
is enabled, which slides the requirement from the coarse branch to the fine one.

The suspect term is identified. The kernel books the momentum its **lumped daughter loss**
deliberately misplaces, and that residual tracks the eigenvalue:

| N | booked momentum misplacement | λ(ℓ=1) |
|---|---|---|
| 201 | 5.25e-3 | 4.14e-7 |
| 101 | 1.56e-2 | 6.28e-6 |
| 51 | 5.85e-2 | 9.07e-5 |
| 21 | 1.34e-1 | 1.40e-3 |
| 11 | 6.59e-1 | 5.42e-3 |

roughly λ ∝ (misplacement)². **Confirmed by switching the term off:** at N = 51,
`lump_limiter = no` raises the booked misplacement 5.85e-2 → 2.08e-1 and the eigenvalue
9.07e-5 → **7.19e-4**, an 8× jump in λ for a 3.6× jump in the residual — the square law,
measured. So the lumped daughter loss *owns* the ℓ=1 eigenvalue, and `lump_limiter` reducing
it 8× is precisely why the grid rule moved from Γ^0.50 to Γ^0.25 when that shipped.

The lumping is deliberate and documented — it exists because
the two-bin deposit makes the loss non-Metzler, so exact placement and positivity are
mutually exclusive. What was not previously connected is that its momentum residual is what
makes the ℓ=1 zero positive, and therefore what sets the whole Γ-dependent cost of this
campaign. **The trilemma is positivity vs exact ℓ=1 momentum vs grid cost, and the code
currently pays in grid cost.**

### M12 — the grid axis, measured offline without a single CMB run

λ(ℓ=1) is the objective that sets the grid rule, and it needs only a background state —
no evolution, no k-modes. `grid_eval` (scratchpad) assembles M_ℓ column by column on an
arbitrary daughter grid, reading f̄ at full resolution through the kernel's
`fermion_bg/boson_bg` override so no PSD interpolation contaminates the ranking.
**Validated against CLASS's own `CLASS_DNCDM_SPECTRUM` dumps at N = 11/21/51/101/201 to
four significant figures across four decades of λ.** Candidates then cost seconds each.

⚠ Unit trap that cost a debugging round: `DNCDMSpecies` stores `Gamma_ = Gamma_raw ·
(1e3/_c_)`, so the kernel wants Γ = 33.36 Mpc⁻¹ for an ini `Gamma = 1e7`, not 1e7. The
symptom was λ too large by exactly 2.998e5 = c in km/s.

**Grid SHAPE is already optimal — a negative result.** Six equidistribution monitors
(f̄_φ, its energy weighting, the Bose factor f̄_φ(1+f̄_φ), √f̄_φ, |dln f̄_φ/dlnq|, and both
daughters summed) were tried at four blend strengths. **Every one is worse than the
shipped log-uniform grid**, and each family degrades monotonically as points are
concentrated harder. log q is the right variable: the error is scale-free, so clustering
at the boson peak only coarsens the tails, where f̄_φ falls decades per e-fold. This kills
the moving-front / injection-adapted grid idea outright — and note the mode sits at the
*peak* of f̄_φ (q ≈ 0.3), not at the injection front (q ≈ 9–13), so a front-following mesh
would have been aimed at the wrong feature entirely.

**The CUT-OFFS are the lever, and q_min is the big one.** At N = 51, λ against the shipped
q_min = 1e-2: 1.67× better at 0.03, 3.27× at 0.1, **6.92× at 0.3**. q_max moves the other
way from the reference in arXiv:1706.02123 — *smaller* is better here (1.85× at q_max = 12),
which is in tension with fixing the known dr_q_max energy leak (that wants ≥ 195, and
q_max = 200 costs 1.7× in λ). The leak only bites at late times when the collision is
nearly off, which is presumably why it has survived.

Growth exp(λΔτ) over the window, Γ = 1e7 — the runaway factor:

| N | q_min = 1e-2 (shipped) | q_min = 0.1 | q_min = 0.3 |
|---|---|---|---|
| 21 | 76.5 | 6.52 | 2.9 |
| 31 | 4.33 | 1.72 | **1.33** |
| 41 | 1.85 | 1.23 | 1.11 |
| 51 | **1.33** | 1.09 | 1.04 |

The accuracy cost of raising q_min is small: the fraction of the boson's ENERGY below the
cut is 7e-6 at 0.1, 1.6e-4 at 0.2, 9.4e-4 at 0.3. So **q_min = 0.3 at 31 points matches the
shipped grid at 51** — a 1.65× point reduction for 0.09% of the sector's energy, with no
code change at all (`dr_q_min` is already an ini key).

**Confirmed against P(k)** at 51 points, Γ = 1e7, vs the full 201-point run:

| q_min | growth | ΔP/P |
|---|---|---|
| 1e-2 (shipped) | 1.33 | 1.83e-3 |
| **0.1** | 1.09 | **3.94e-4** |
| 0.3 | 1.04 | 5.61e-3 |

q_min = 0.3 overshoots — growth keeps falling but truncation takes over and z_eq shifts. The
optimum is ~0.1, and at fixed N the error tracks the runaway, so **51 points at q_min = 0.1
is 4.6× more accurate than the shipped cut-off and 7× faster than 201 points (663 s vs
4655 s), for one ini key.** Below 51 points the error climbs steeply (2.3e-3 at 41,
3.4e-3 at 26, 6.8e-3 at 21).

### M13 — how much the eigenvalue is worth, isolated with the Γ lever

λ ∝ Γ, so running identical grids at two Γ separates the runaway from the discretisation
error with no model and no fitting. 21 vs 51 daughter points, q_min = 1e-2:

| Γ | growth | ΔP/P |
|---|---|---|
| 1e5 | ≈ 1 | **1.98e-4** |
| 1e7 | 76.5 | **4.32e-2** |

**218× for the same discretisation.** So the discretisation error at 21 points is ~2e-4, and
essentially all of the Γ = 1e7 error is the unstable mode.

⚠ An earlier attempt to split the two by fitting `err ≈ A(growth−1)^p` at fixed N and
extrapolating across N was **wrong** — it over-predicted the runaway at 21 and 26 points
badly enough to give negative residuals. Two points, extrapolated threefold. The Γ lever is
the controlled version and should be used instead.

**This is the case for the deflation route.** If the ℓ=1 eigenvalue is removed, 21 points
should reach ~2e-4 at Γ = 1e7 — i.e. 201 → 21 daughter points, and from the measured timings
4655 s → 219 s, a **~20× speed-up**, against ~7× from the cut-off alone.

### Two consequences

**For the reduced model.** `F_1(x) = x` *is* the ℓ=1 momentum sum rule — analytic, exact at
every x, independent of resolution. An operator built on the transfer functions has that
zero by construction, so it cannot develop this mode, and its grid requirement should not
scale with Γ at all. First implementation milestone, and it is one number: build the reduced
operator, dump its spectrum, check max Re λ = 0 at ℓ = 1 for any daughter DOF.

**For the existing code, independently of any of this.** The target is now a specific
number: the ℓ=1 momentum residual of the lumped daughter loss. `lump_limiter` already cut it
once and bought a factor Γ^0.25 over Γ^0.50; cutting it again buys more, on a measured square
law (λ ∝ residual², N ∝ Γ^(1/s)). That is a far smaller change than the reduced model and it
should be tried first. The full causal chain, every link measured:

```
    lumped daughter loss  →  misplaces ℓ=1 momentum  →  positive ℓ=1 eigenvalue
    λ ∝ (residual)² ∝ Γ   →  runaway exp(λΔτ)        →  grid rule N ∝ Γ^(1/s)
```

### M5 — the daughter δf(q) has fine q-structure, and it does not matter

An SVD of δf_i(q, ℓ=2) over the window needs >32 q-shapes to reach 50% variance. That
structure is injection-time phase: a daughter at comoving q was created when the
emission band `[(ε−q₁)/2, (ε+q₁)/2]` swept over q, and free-streams with phase
`k(τ − τ_inj(q))`. Because free streaming carries no q (M3), that phase never reaches the
observables — the metric integrates it away, and it is not a reason to keep a q-grid.

But note what M4 adds: the *collision* does not integrate it away, and the collision is
where the daughter grid is actually needed. So M5 says a reduced daughter is viable for
transport and says nothing about whether it is viable for the collision. **Target the
moments, but prove the collision separately.**

### M6 — the reference itself is only converged to a few percent at high Γ

r1 vs r2 differ only in the daughter *perturbation* grid (201 → 101 points), same
background. Scale-relative max over 2 ≤ ℓ ≤ 2500:

| cell | C_ℓ^TT | C_ℓ^TE | P(k), z=0 |
|---|---|---|---|
| G1e5_m0.06 | 1.2e-3 | 2.1e-3 | 2.8e-4 |
| G1e7_m0.06 | 2.7e-3 | 5.1e-3 | 4.5e-4 |
| G1e7_m0.3 | 1.6e-2 | 2.5e-2 | 1.3e-3 |
| G1e9_m0.3 | 4.3e-2 | 6.8e-2 | 1.3e-3 |

The daughter grid requirement is a *discretisation* cost, not a physics one, and at
Γ = 1e9 the shipped grid is not converged in C_ℓ^TT to better than ~4%. This sets the
accuracy target honestly: the reduced model has to match a converged reference, and it
may well be better conditioned than an under-resolved exact run, because it has no
daughter grid to ratchet.

---

## 2. The scheme

### 2.1 Variables

| species | representation | count (l_max = 17) |
|---|---|---|
| ν_H (parent) | unchanged q-resolved hierarchy | 16 × 18 = 288 |
| ν_l | two q-moments per multipole | 2 × 18 = 36 |
| φ | two q-moments per multipole | 2 × 18 = 36 |

**360 variables against the present 7812.** The two daughter moments are

```
    N_{i,ℓ}(τ) = ∫ dq q²   δf_i(q, ℓ)        (number-weighted)
    E_{i,ℓ}(τ) = ∫ dq q³   δf_i(q, ℓ)        (energy-weighted)
```

Two, not one, and not three, for a reason given in §2.5: they are exactly the moments
the collision must conserve, and they span the equilibrium response.

### 2.2 Free streaming — exact

By M3 both moments obey the standard massless hierarchy with no closure error. `E` is
what the metric reads (δρ, θ, σ all come from `E_{ℓ=0,1,2}`); `N` is carried because the
collision needs it.

### 2.3 The decay source — exact, via two Legendre transfer functions

For a parent bin at (q, ε = √(q² + a²m²), x = q/ε) and rest-frame emission cosine u:

```
    q₂ = (ε/2)(1 + x u)        cos θ₁₂ = (u + x)/(1 + x u)
    q₃ = (ε/2)(1 − x u)        cos θ₁₃ = (x − u)/(1 − x u)
```

Define

```
    G_ℓ(x) = ∫₋₁¹ (du/2)          P_ℓ(cos θ₁₂)      number-weighted transfer
    F_ℓ(x) = ∫₋₁¹ (du/2) (1 + x u) P_ℓ(cos θ₁₂)      energy-weighted transfer
```

`F_ℓ` is *exactly* the function `DNCDM_DR_Species::AddCouplingDerivs` already computes —
verified against its closed forms and its upward recursion to 1e-13 for x ≥ 0.5. `G_ℓ`
comes from the same Legendre sweep at no extra cost. Both daughters share the same pair:
the boson's integrand is the fermion's under u → −u, and the map leaves the integral
invariant (checked to 2e-16).

The source into daughter i, multipole ℓ, is then a sum over parent bins:

```
    dN_{i,ℓ}/dτ|src = Σ_b  dq_b q_b² (a²mΓ/ε_b) · w^N_b · G_ℓ(x_b) · δf_H(q_b, ℓ)
    dE_{i,ℓ}/dτ|src = Σ_b  dq_b q_b² (a²mΓ/ε_b) · (ε_b/2) · w^E_b · F_ℓ(x_b) · δf_H(q_b, ℓ)
```

with `w^{N,E}_b` the band-averaged quantum-statistics weight `(1 − f̄_l(q₂))(1 + f̄_φ(q₃))`
folded into the same u-quadrature. No k-dependence, no daughter grid, one 1-D quadrature
of ~16 nodes per parent bin per RHS — and the whole thing can be tabulated on the
background time grid because it depends on a alone.

### 2.4 Conservation is analytic, not discrete

The present kernel needs a discrete transition network with transpose scatter weights to
make number/energy/momentum conservation hold on the grid. The reduced scheme gets it
from three identities, verified to machine precision at every x:

```
    G_0(x) = 1                number:   one daughter of each kind per decay
    F_0(x) = 1                energy:   ⟨q₂⟩ + ⟨q₃⟩ = ε
    F_1(x) = x                momentum: ⟨q₂ cos θ₁₂⟩ + ⟨q₃ cos θ₁₃⟩ = q
```

These hold for the continuum kernel at any resolution, so the reduced operator is
exactly conservative by construction. This is a structural improvement, not just a
simplification: the `Γ^0.5` daughter-grid law was a one-signed discretisation ratchet in
a scheme whose conservation was only as good as its grid.

### 2.5 The inverse decay — the RTA, with rates derived not fitted

Linearising the band factor `Λ = f_l f_φ − f_H + f_H(f_l − f_φ)` groups the terms by
which species carries the perturbation:

```
    δf_H(q₁) × [ −1 + f̄_l(q₂) − f̄_φ(q₃) ]      parent loss   — local in q₁, exact
    δf_l(q₂) × [  f̄_φ(q₃) + f̄_H(q₁) ]          fermion loss  — needs a daughter shape
    δf_φ(q₃) × [  f̄_l(q₂) − f̄_H(q₁) ]          boson loss    — needs a daughter shape
```

The parent's loss is local in its own momentum and carries no angular factor, so it is
exact and multipole-independent — which is precisely why `CollisionDiagonal` is already
ℓ-independent. It needs no approximation at all.

The two daughter-loss terms are the only place the reduced scheme approximates anything.
Projecting the q-dependent rate onto the (N, E) basis gives, per daughter and per
multipole, a 2×2 relaxation matrix plus a cross-term to the other daughter:

```
    d/dτ (N, E)_{i,ℓ} ⊃ − R_i(a) (N, E)_{i,ℓ} + C_{ij}(a) (N, E)_{j,ℓ}
```

with `R_i`, `C_ij` band integrals of background quantities against the two basis shapes.
Everything is a function of a alone.

For ℓ ≥ 2 there is no conservation constraint, so this is a pure relaxation toward
isotropy at rate `R_i(a)` — the RTA. Its rate is not postulated: it is the decay rate
times `F_2(x)` averaged over the parent, and

```
    F_ℓ(x) → c_ℓ x^ℓ    as x → 0     (c_2 = 4/5, c_3 = 4/7, measured)
```

so with the 1/γ time dilation already in the rate, the shear isotropisation rate comes
out as **Γ γ⁻³** — the Hannestad & Raffelt kinematic scaling, produced by the Legendre
machinery rather than assumed. That is the answer to "can the decay-only Legendre trick
be used together with an RTA": they are the same object. The transfer function *is* the
transport rate.

### 2.6 Detailed balance must be exact, and two moments is why

The occupations are extreme — measured peak bare values `f̄_H ≈ 0.3`, `f̄_l ≈ 1.01`
(Pauli-saturated), `f̄_φ ≈ 140` (strongly Bose-enhanced) — so `Λ` is a small residual of
large terms. An approximation that breaks the cancellation injects error proportional to
the *large* terms. This is the same trap the `balanced_gather` work hit, and it is the
single largest risk in the design.

An equilibrium perturbation is a shift of (T, μ_i) preserving η_H = η_l + η_φ, i.e.
`δf_i ∈ span{∂f̄_i/∂μ, ∂f̄_i/∂T}` — a two-dimensional space. Hence:

* **one** daughter moment cannot represent an equilibrium perturbation, so a one-moment
  scheme necessarily breaks detailed balance;
* **two** moments span it exactly, so the reduced operator can be constructed to
  annihilate equilibrium perturbations identically.

The construction: build the basis as the equilibrium response pair, and define the
reduced loss coefficients from the reduced gain by the equilibrium ratio, so that
`L·v_eq = 0` holds by algebra rather than by cancellation. This must be a unit test with
a machine-precision tolerance, not a tolerance-tuned check.

### 2.7 The window — a physically triggered approximation switch

Two switches, in CLASS's existing `ppw->approx` style:

* **collision off** while `max_q [aΓ (am/ε)] / (aH) < ε_on` (before) or once the parent's
  share of the total density falls below `ε_off` (after). Outside the window the daughters
  are ordinary free-streaming radiation and their two-moment hierarchies are exact.
* **parent dropped** once ρ_H/ρ_tot < ε_drop: its perturbations no longer source the
  metric at any level the observables can see, and dropping it removes the last stiff
  eigenvalue. The daughters' hierarchies continue unchanged.

Both thresholds are convergence knobs with a physical meaning, and both must be swept.

### 2.8 A small system is a stiff-solver system

With 16 parent bins and a handful of daughter moments, the collision block at fixed ℓ is a
dense matrix of background quantities a few tens of rows on a side. Three integrators then
become available that are not available at 7812 variables:

* **ndf15** (Thomas Tram's stiff evolver, `evolver=1`). It was re-measured *dead* on the
  full system at Γ = 1e8 — but the reason is system size, not stiffness: its cost is
  dominated by numerical-Jacobian formation and the LU, both super-linear in N. At
  N ≈ 360 that is affordable, and the stiffness it exists to absorb is exactly the −aΓ
  eigenvalue this model is paying for. **Reducing the system is what revives ndf15**, and
  if it works it removes the need for any bespoke exponential machinery. This is the
  cheapest option and should be tried first. Caveat from the injection-adapted-grid work:
  ndf15's numjac sparsity lock needs a static source structure, so the window switch of
  §2.7 must not change the sparsity pattern mid-run — switch coefficients to zero rather
  than removing equations.
* **ETD** on the exact block rather than just its diagonal — the existing `evolver=3`
  already integrates the diagonal analytically; a block of this size can be exponentiated
  whole.
* Plain `rkdp45`, if the window switch alone makes the stiff phase short enough to brute-force.

The point is that the reduction and the integrator question are the same question. Whichever
of the three wins, the in-window step count stops being the thing that scales with Γ.

## 3. Cost

| | present | reduced |
|---|---|---|
| variables / k-mode | 7 812 | 360 in-window, 72 outside |
| collision per RHS | ~16 × 201 × 4 nodes × 18 ℓ, memory-bound (97% of runtime) | ~16 × 16 × 18, cacheable, k-independent |
| steps | see M2 | 3–17× fewer, and Γ-independent if §2.8 lands |

Present cost at Γ = 1e7, m = 0.06, r1: 6 710 s × 32 cores ≈ 60 core-hours per model.
At Γ = 1e9: ≈ 210 core-hours. A ΛCDM+ncdm run is of order a core-minute.

The reduced model should land within a small factor of a standard massive-neutrino run.
The honest claim is **two to three orders of magnitude**, with the residual uncertainty
in the in-window step count; the range is wide because the two ends (§2.8 landing or
not) differ by the Γ scaling itself.

---

## 4. Validation

Ordered so that each step can fail cheaply before the next is attempted.

1. **Identities, as unit tests at machine precision.** `G_0 = 1`, `F_0 = 1`, `F_1 = x`
   over the x range; sector number/energy/momentum sum rules on the assembled operator;
   `L · v_eq = 0` for the equilibrium perturbation. These are the tests that would have
   caught the historical ratchets.
2. **Operator against operator, no cosmology.** Drive the reduced operator and
   `DecayTransitionKernel` from the same background state and the same δf, and compare
   the moment rates directly. This isolates the closure from everything else and is where
   the boson's harder behaviour (M4) will show first.
3. **Trajectory, against the existing dumps.** The `a4PiDot_coll` columns and the
   `CLASS_DNCDM_PERT_PSD` dumps already give a per-time-step reference at Γ = 1e5…1e9 and
   two masses. Reproduce M4's table with derived rather than fitted coefficients — if the
   derived rates do not match the fitted ones, the derivation is wrong and the fit was a
   coincidence.
4. **Observables.** C_ℓ and P(k) against the harness matrix, with the scale-relative
   metric (`collect.py --cl`), and against the *converged* reference where M6 says the
   shipped one is not.
5. **Convergence knobs.** l_max of the daughter hierarchies, ε_on/ε_off/ε_drop, and the
   number of daughter q-moments (2 → 3 → 4 recovers the exact scheme in the limit). The
   moment count is the honest convergence dimension: if 2 is not enough, 3 is available
   and still costs nothing.

Reference runs come from the exact code; there is no external oracle for this model, so
validation is by limits, convergence and internal identity — never by declaring one code
correct.

---

## 5. Risks, and what falsifies the design

* **Detailed balance (§2.6).** The largest risk. Falsified by test 1: if `L·v_eq` cannot
  be made to vanish at machine precision in the two-moment basis, the basis is wrong.
* **The boson.** M4 shows φ is consistently the worse of the two daughters, and M5 shows
  why (it starts empty and is built purely by injection, so its spectrum has the widest
  history). If test 2 says two moments are not enough, φ gets three and ν_l keeps two —
  the scheme does not require them to match.
* **The inverse-decay legs may not be reducible at all.** This is now the headline risk,
  not a detail: M4 says the daughters' collision term is not dominated by the decay
  source, and the occupations (f̄_φ ≈ 140) say why. If two moments cannot carry those legs,
  the honest fallback is a coarse daughter grid with the same analytic source — which
  keeps M1/M2/M3 and the conservation identities, and loses only the largest factor.
* **High Γ.** The collision is ~10× stronger in the ratio sense at Γ = 1e9 than at 1e7,
  and the parent's momentum profile gets richer with Γ (M4). Every check must be run at
  Γ ≥ 1e8, which needs the HPC.
* **The parent's inverse-decay gain.** The one term that reads a daughter shape back,
  and the term the fluid test (M0) failed on. Bounded by test 2.
* **The window thresholds.** Setting them from the background is straightforward; setting
  them from the background *conservatively enough for every (Γ, m)* is not. They must be
  swept, not tuned on one cell.

## 6. Explicitly out of scope

Newtonian gauge and tensors (the exact composite is synchronous scalar only, and the
reduced one inherits that); the high-Γ tight-coupling branch (a real opportunity — the
physics saturates above Γ ≳ H(a_nr) while the cost keeps growing — but a separate design);
and any change to the background.


---

## 7. The deflation: validated design (2026-08-12)

Adopted in preference to the grid work, which is capped at ~7× (M12) while this is worth
~20× (M13). `dr_q_min = 0.1` is parked as a pending test, not adopted.

### 7.1 What must be deflated

Measured at N = 51, Γ = 1e7, over the collision window — **this corrects an earlier claim
of mine that only ℓ=1 misbehaves**, which came from a single-`a` snapshot:

| ℓ | positive eigenvalues |
|---|---|
| 0 | none for a ≲ 1e-2, then **1–2** late in the window (up to +1.88e-4) |
| 1 | exactly **1**, throughout (6.6e-6 → 1.10e-4, peaking at a ≈ 1.3e-2) |
| ≥ 2 | none at any a, to machine zero (checked ℓ = 2,3,4,6,8,12) |

"Up to two at ℓ=0 and one at ℓ=1" — identical to arXiv:1706.02123. So the recipe is theirs:
deflate every eigenvalue with Re λ > 0 at ℓ ∈ {0,1}, to its asymptotic value of zero.

### 7.2 Two-sided, not one-sided — measured

| deflation | kills λ | momentum residual ‖wᵀM‖ | ‖ΔM‖/‖M‖ at N=51 |
|---|---|---|---|
| **two-sided** `M − λ u vᵀ/(vᵀu)` | ✅ machine zero | **unchanged**, 1.6e-8 → 1.9e-8 | 9.5e-4 |
| one-sided (right vector only) | ✅ | **degrades**, 1.6e-8 → 9.4e-7 | 0.29 |

One-sided needs only a power iteration but changes the operator by 30–70% and wrecks
momentum conservation. Two-sided is surgical, so the **left eigenvector is required** — which
is why the matrix must be assembled rather than merely applied.

### 7.3 Finding the eigenpair without LAPACK

CLASSpp has no BLAS/LAPACK/Eigen. Plain power iteration finds the wrong mode (the dominant
eigenvalue by magnitude is the most negative, ≈ −0.34). Shift-invert on (M − σI) is the
tool, and `ludcmp`/`lubksb` in `include/evolver_ndf15.h` supply the dense LU — no new
dependency. The left vector comes from a second LU of (M − σI)ᵀ (~10 ms; a transposed
back-substitution entry point would avoid it but is not needed).

**Measured σ requirement: σ ≈ 3–10 λ.** Validated against dgeev:

| N | λ | σ that converge | deflated λ_max |
|---|---|---|---|
| 11 | 5.4e-3 | 3e-3 … 3e-2 | 2e-16 |
| 21 | 1.4e-3 | 1e-3 … 1e-2 | 7e-18 |
| 51 | 9.1e-5 | 3e-4 … 1e-3 | 3e-16 |

A **descending σ-ladder** with an acceptance test — eigenpair residuals ‖Mu−λu‖, ‖Mᵀv−λv‖
small **and** λ > 0 — is robust and self-terminating. It correctly finds nothing on fine
grids, where λ (4e-7 at N=201) is buried in the cluster of exact conservation zeros and is
harmless anyway (growth 1.001).

### 7.4 Cost, and the one kernel change needed

Assembly is N operator applications per (a, ℓ). `DumpCollisionSpectrum` calls
`PrepareTransitions` per column, which is ~400× more than needed: the only state the
ℓ-recurrence carries between columns is `next_l_`, since the ℓ=0 call re-seeds the Legendre
arrays for every node. **A `ResetMultipoleRecurrence()` entry point (one line, `next_l_ = 0`)
lets `PrepareTransitions` be hoisted out of the column loop.** With that, the precompute is
~30 s per model — a one-off, amortised over ~1500 k-modes, against a 4655 s run.

The correction itself is rank-1 and matrix-free: `dF −= λ u (v·F)/(v·u)`, O(N) per RHS.

### 7.5 Open design decisions

* Where the precompute runs (after the background solve; the species needs a post-background
  hook) and where (λ, u, v) are stored — per background row over the window only, ~10 MB.
* Interpolation in `a`: λ is smooth, but eigenvectors need a sign convention (align each with
  its predecessor) and the near-degeneracy with the zero cluster must be watched.
* Whether the ℓ=0 second mode ever needs deflating in practice, or is always negligible.
* ndf15 interaction: its numjac sparsity lock needs a static source structure, so the
  correction must be present (possibly with λ = 0) at every step rather than switched on.
