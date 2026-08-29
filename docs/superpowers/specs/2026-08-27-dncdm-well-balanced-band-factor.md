# DNCDM: the well-balanced band factor, and the lumped loss it retires

**Written 2026-08-27**, at the end of a session whose brief was to break the Γ ceiling
in the DNCDM background; **§5b added 2026-08-28**, when the perturbation chain rule §6
had deferred turned out to be both fixable and cheap. **Status: implemented and
measured**, on branch `ratchet-b-balanced-gather` (`9d6a4a98`, `balanced_pert` on top).
Everything below is a record of what was built and what it measures, not a plan.

Evidence tags follow `2026-08-17-dncdm-minimal-knob-set.md`: **[M]** measured in this
session, **[P]** earlier campaign, not re-measured here, **[?]** no evidence.

**What this supersedes.** The `balanced_gather` field doc on branch
`dncdm-inverse-qs-archive` says the flag is "NOT SHIPPABLE YET, and the reason is
structural rather than a loose end", citing a 9.0e-2 drift of `CollisionDiagonal` from
finite differences. **That diagnosis is wrong** and §4 says why: the drift was a wrong
linearisation of the lumped deposit, not a structural conflict. The same doc's companion
claim — that retiring the lumped loss and fixing the gather "share a cause and should be
solved together" — is right, and §3 gives the mechanism it was missing.

---

## 1. What the problem was

The daughter grid had to be refined as Γ rose: `bins/decade ~ Γ^0.22` with the harmonic
lump limiter in place, `~Γ^0.5` before it **[P]**. Nothing physical justifies that. The
sector *saturates* above Γ ≳ H(a_nr) — `a_½`, `q*` and `Ecomov` are all Γ-independent
there — so a Γ-dependent grid is a statement about the discretisation, not the model.

The tell is one-sidedness. The error at fixed grid grows as Γ¹ and **every entry is on
the same side of the limit**. That is a ratchet: a per-collision bias accumulating with
the collision count, which is proportional to Γ. Two of them, in fact.

Cost was measured at Γ^0.86 = steps Γ^0.62 × grid Γ^0.22 **[P]**. Only the grid term is
addressed here; the step term is untouched.

## 2. Ratchet B — the band factor reads high because e^−q is convex

```
Λ = f_l f_φ − f_H + f_H(f_l − f_φ)
```

vanishes identically in the continuum when the daughters are FD/BE with
μ_H = μ_l + μ_φ, because **q2\* + q3\* = ε1 exactly** at every node — that identity is
built into the kinematics, not approximated. On the grid f_l(q2\*) and f_φ(q3\*) are
reached by a two-bin **linear-in-q** gather, and e^−q is convex, so both come out high
and Λ keeps a one-signed O(dq²) residual.

That residual is a spurious *inverse decay*. It re-creates a parent at ε1 out of
daughters at q2\*+q3\*, energy-conserving at that instant — but the parent's comoving
energy then grows as a·m before it decays again. Each spurious cycle pumps energy into
the sector, and the cycle count goes as Γ.

**The fix is to interpolate the variable that is linear at equilibrium.** For a fermion
η = ln((1−f)/f), for a boson η = ln((1+f)/f); at equilibrium both are (q−μ)/T. The
two-bin split has Σ w_e = 1 and Σ w_e q_e = q\*, so it reproduces a linear function
exactly, and η_l + η_φ telescopes onto ε1. `Config::balanced_gather`.

The **scatter weights are untouched**. This changes the value of a coefficient, never how
a transition is split between bins, so every discrete number/energy identity is exactly
as before — asserted under the flag in `decay_kernel_test`, at round-off **[M]**. The
separation it relies on is the one the `fermion_bg`/`boson_bg` override already used:
*where* the background occupation is read is a free choice; the deposit is not.

Equilibrium residual, parent leg (which is Λ alone), 64 / 128 daughter bins **[M]**:

| gather | 64 | 128 | |
|---|---|---|---|
| linear-in-f | 8.95e-06 | 2.31e-06 | falling as dq² |
| log-odds | **1.55e-18** | **2.86e-18** | flat, at round-off |

Flatness is the claim, not smallness: a scheme that merely converged faster would still
fall with the grid.

## 3. Ratchet A — and why fixing B is what unlocks it

The lumped daughter loss evaluates Λ's f-linear part at the deposit bin's **own**
occupation rather than at the gathered one, blended by the harmonic limiter θ = 1 − H/A.
It exists for positivity: with the exact two-bin split the loss charged to bin k is
proportional to its neighbour's occupation as well as its own — a negative off-diagonal —
so at an injection front an empty bin is billed for particles it does not hold and the
right-hand side drives it negative. That is a property of the RHS, so no integrator can
repair it; rejection, clamping f at 0 and evolving √f were each tried and each stalled
the solver **[P]**.

It was recorded as an *impossibility*: conserving both number and energy fixes the two
deposit weights, so exact energy placement and positivity are mutually exclusive for a
two-bin deposit. The price is an O(dq) energy misplacement, booked in
`split_energy_residual()` — and it is the **dominant** term in the grid requirement at
high Γ.

**The impossibility is a property of the linear gather, not of the deposit.** The
off-diagonal is the loss term g·G evaluated at f_k = 0. For the linear gather
G|_{f_k=0} = Σ_{e′≠k} w_{e′} f_{e′}, still positive — that is the bill. The log-odds
gather is a weighted **geometric** mean, so **G|_{f_k=0} = 0 identically**. The bill is
not smaller, it is gone, and the RHS at an empty bin reduces to its gain:

```
fermion, at f_l[k]  = 0:   df_k = +  w_e C f_H (1 + G_φ) / meas   ≥ 0
boson,   at f_φ[k]  = 0:   df_k = + 2w_e C f_H (1 − G_l) / meas   ≥ 0
```

for any admissible state (f_H ≥ 0, G_l ≤ 1, G_φ ≥ −1) — exactly the condition the lumped
scheme also needs. The exact split is Metzler again, so the lumping can go.

Measured on front states, with the known-broken corner run as a **control** so the three
positive results cannot be vacuous **[M]**:

| gather | loss | worst empty-bin RHS |
|---|---|---|
| linear | lumped (shipped) | 0 |
| linear | exact | **−1.57e+05** (control: must break) |
| balanced | lumped | 0 |
| balanced | **exact** | **0** |

The two flags are therefore one scheme, and the input layer couples them:
`balanced_gather = yes` turns lumping off unless `lumped_loss` is set explicitly. They
stay separate in `Config` because they are independent mechanisms and the 2×2 is what the
measurement needs.

## 4. The diagonal, which is what actually blocked this

`CollisionDiagonal` is the operator the ETD evolver integrates analytically, on the
background and the perturbations. It read

```
dΛ_k/df_k = g [ θ + (1−θ) dG/df_k + dev · dθ/df_k ]
```

which is the same algebra as the truth **only while the band-factor gather and the
lumping's reference gather coincide**. They do not under `balanced_gather`, and the
rewrite multiplies the chain factor by (1−θ) — so it **drops it entirely at θ = 1**,
which is exactly where the limiter puts a front. That, not a structural conflict, is the
9.0e-2 drift the archive recorded.

Differentiating the deposit `Λ_k = Λ + g θ (f_st[k] − G_lin)` literally:

```
dΛ_k/df_k = g [ dG/df_k + θ (1 − w_k) + (f_k − G_lin) dθ/df_k ]
dG/df_k   = w_k · G(1∓G) / (f_k(1∓f_k))        (upper sign fermion)
```

Three terms, each of which has been got wrong at some point: the chain factor above; the
`(1 − w_k)` — the deviation is taken against the *linear* state-grid gather whatever Λ
used, because that is the identity Σ w_e d_e = 0 the lumped loss rests on; and the dθ
term, whose omission is a 4.1e-2 error **[P]**.

Against central differences **[M]**:

| state | linear | balanced |
|---|---|---|
| smooth random | 2.18e-09 | 7.29e-09 |
| log-uniform over 8 decades | 5.31e-06 | **7.69e-07** |
| same, chain factor dropped | — | 3.53e-01 |

The 8-decade state is new. The chain factor is 1 on a flat stencil, so the pre-existing
smooth-state test could not see it at all. Its finite difference needs a **relative and
coarse** step: df_c is dominated by transitions out of far larger bins, so the quotient is
a cancellation and its round-off branch takes over below h/y ~ 1e-4 (measured
5.3e-6 / 8.3e-5 / 4.2e-4 / 1.8e-3 at h/y = 1e-4 … 1e-7 — growing as h *shrinks*).

**One diagonal, two consumers.** On the background kernel `CollisionDiagonal` is
∂(df_i)/∂(f_i) and the chain factor belongs in it. On a perturbation kernel the daughter
occupations are a fixed table, the state is F, and the F-gather is linear whatever the
band factor was built from — so there it must stay the linear expression or it will not
match the operator ETD subtracts it from. `chain_diag_` separates them, keyed on whether
bg views were supplied, which is a construction-time property.

## 5. What it buys

`Ecomov` against the converged limit, background only, `momenta_bins = 32`,
`dr_q_min`=1e-2, `dr_q_max`=50 **[M]**:

| bins/dec | Γ=1e7 (limit 1.2514) | | Γ=1e8 (1.2331) | | Γ=1e9 (1.2234) | |
|---|---|---|---|---|---|---|
| | lin+lump | bal+exact | lin+lump | bal+exact | lin+lump | bal+exact |
| 14.3 | 1.2726 | **1.2515** | 1.3937 | **1.2325** | 2.9460 | **1.2233** |
| 28.6 | 1.2540 | **1.2514** | 1.2435 | **1.2324** | 1.3209 | **1.2233** |
| 57.2 | 1.2519 | **1.2514** | 1.2336 | **1.2326** | 1.2295 | **1.2234** |

`bal+exact` is converged at 14.3 bins/decade at every Γ — flat across a 4× refinement —
and it reproduces conv2's published Γ=1e7 limit 1.2514 exactly, which is the check that it
converges to the *right* place rather than a convenient one. **The Γ-dependence of the
daughter-grid requirement is gone.**

At the same grid the two schemes cost the same: the log-odds inverse costs an `exp` per
node, the exact split saves the θ and dθ machinery, and the two roughly cancel. Wall
clock at 53 bins on an idle machine, background only **[M]**:

| Γ | lin+lump | | bal+exact | |
|---|---|---|---|---|
| | Ecomov | wall | Ecomov | wall |
| 1e8 | 1.3937 | 3.7 s | 1.2325 | 3.8 s |
| 1e9 | 2.9460 | 19.0 s | 1.2233 | 19.2 s |
| 1e10 | 5.9685 | 77.2 s | 1.2188 | 72.6 s |
| 1e11 | 8.9422 | 355.7 s | 1.2144 | 302.5 s |
| | fitted | **Γ^0.654** | | **Γ^0.627** |

**That fit closes the loop on the cost model.** Cost was measured at Γ^0.86 = steps
Γ^0.62 × grid Γ^0.22 **[P]**, in a different session and a different config. Here, at a
*fixed* grid, both schemes come out at Γ^0.63–0.65 — which is the step term, recovered
independently — and adding the grid law back gives 0.65 + 0.22 = 0.87. So for the shipped
scheme the fixed-grid column above is **not** at fixed accuracy (it is 7× wrong at 1e11
and the grid would have to grow), while for `bal+exact` it is. The cost exponent at fixed
accuracy goes **Γ^0.86 → Γ^0.63**, and the difference is exactly the ratchet.

⚠ **The flatness claim is anchored at Γ ≤ 1e9.** Above it `bal+exact` still drifts
—1.2233 (1e9), 1.2188 (1e10), 1.2144 (1e11), about 0.36% per decade and not obviously
converging. See §6: whether that is under-resolution at 14.3 bins/decade, a residual bias,
or a real Γ-dependence of the saturated limit is the same open question as the 0.8% gap
against conv2 §4's universality prediction, and it is not settled here.

## 5a. The moment converges; the SHAPE is a separate question **[M]**

`Ecomov` is a moment, and the two-bin deposit conserves number and energy *per
transition* whatever the gather or the quadrature does. So §5 is structurally blind to the
daughter PSD's shape — which is what the perturbation hierarchy reads (dlnf/dlnq). Three
measurements, all at Γ=1e9, `bal+exact`:

**Do not measure shape with |Δln f|.** The daughter PSD falls ~87 decades across its
injection front, so a sub-bin shift in the front's position makes ln f differ by *tens*
right there while the distribution is essentially unchanged. That metric reads ~1e-2 and
is flat in resolution — it is measuring the cliff, not the physics. Use the total-variation
distance of the normalised energy distribution w = dq q³ f: bounded, and it is the shape
the source integrates against.

**Grids must NEST or the metric measures its own interpolation.** N_q = 53/106/212/424 do
*not* nest (h = ln(q_max/q_min)/(N−1) gives steps /52, /105, /211, /423). N−1 doubling —
53/105/209/417 — makes every coarse point exactly a fine point.

With both fixed, TV against the N_q=53 answer at Γ=1e8:

| N_q | 105 | 209 | 261 | 313 | 365 | 417 |
|---|---|---|---|---|---|---|
| f_l | 1.7e-4 | 1.6e-4 | 1.6e-4 | 3.1e-4 | 4.4e-4 | 5.4e-4 |
| f_φ | 2.5e-4 | 4.1e-4 | 7.7e-4 | 1.2e-3 | 1.5e-3 | 1.7e-3 |

The shape converges far more slowly than the moment: the increments peak near N≈300 and
then shrink, so 14.3 bins/decade sits about 2e-3 (TV) from the refined answer. That is a
small number and it does not undermine §5, but it is a *different* number from the moment's
1e-4, and a background convergence ladder cannot be used to set a shape-sensitive knob.

**`emission_gauss = 1` is a floor no other knob lifts.** ρ-weighted ⟨|Δln f|⟩ at fixed
daughter grid (dr_N_q = 106), against momenta_bins = 128 + Gauss-3:

| momenta_bins | 32 | 64 | 128 |
|---|---|---|---|
| Gauss-1 (f_φ) | 7.2e-4 | 6.8e-4 | 6.8e-4 |
| Gauss-2 (f_φ) | 5.2e-5 | 4.6e-6 | 1.1e-6 |
| Gauss-3 (f_φ) | 4.9e-5 | 4.1e-6 | — |

With the midpoint rule the parent grid is *irrelevant* to the daughter shape; switch to two
points and momenta_bins becomes the binding constraint, converging at ~O(h³·⁵). Three
points buy almost nothing over two. So the pairing is **Gauss-2 with momenta_bins ≥ 64**,
which is where `2026-08-12`'s q-oscillation work arrived from a different direction.

## 5b. The perturbation chain rule — shipped, and what it costs **[M, 2026-08-28]**

`Config::balanced_pert`. §6's first bullet used to say this was not done and needed
"a limiter nobody has designed"; that bullet is superseded by this section.

**The defect was 1.1e-1, and it was not only an inconsistency.**
`ApplyPerturbationOperator` at l = 0 is a *directional derivative* of
`ComputeBackgroundDerivs`, so it can be pinned against a central difference of the
background right-hand side — two code paths, one of them not differentiated by hand at
all. That test did not exist. It reads **[M]**:

| gather | operator vs. FD Jacobian |
|---|---|
| linear | 8.8e-12 |
| balanced, F gathered linearly | **1.1e-1** |
| balanced + chain rule | **3.3e-11** |

`test_perturbation_diagonal` could not see it: it pins the operator against
`CollisionDiagonal`, and both were written from the same derivation, so a wrong gather in
both agrees with itself. A latent bug came with it — on a SINGLE-GRID kernel `chain_diag_`
was already true, so `CollisionDiagonal` carried the chain factor while the operator did
not, and the two disagreed by **81%** **[M]**. Production never hit it (the per-k kernels
have separate bg grids, so `chain_diag_` was false there), but each half of the old rule
was wrong for one of its two consumers.

**The divergence is not a limiter problem; it is a change of variable.** The product the
operator actually forms is

```
dG/df_e * F_e  =  G(1∓G) (F_e/f_e) / (1∓f_e)  =  G(1∓G) Ψ_e / (1∓f_e)
```

so the chain rule **is** linear interpolation of Ψ = F/f, rescaled by the gathered
background — and Ψ is smooth across the injection front where F and f share the same
cliff. That is the same argument that motivated `balanced_gather`: F = f Ψ inherits f's
convexity, so gathering F linearly reads high for exactly the reason gathering f linearly
did. The exact object is bounded (X_e f_e = G(1∓G)/(1∓f_e) ≤ 1); only floating point
diverges, when F_e carries round-off from a scale f_e has fallen 40 decades below.

**The cap must be a function of the background alone.** Never of F: the operator has to
stay linear in F, or `ReducedCollisionOperator`, which assembles it one column at a time,
is assembling something that does not exist. So no clip on |δη|. What is capped is the
dimensionless X = (balanced sensitivity)/(linear sensitivity), which is 1 for the linear
gather, by a SOFT MINIMUM `X/(1+(X/cap)⁴)^(1/4)`. The obvious harmonic `X/(1+X/cap)` is
wrong for a tight cap: it biases *every* coefficient by X/cap, and at cap = 1e4 it moved
the operator 2.2e-4 off its own FD Jacobian on a state where X was only ~10 and no capping
was wanted at all. The quartic leaves that state at 2.5e-13 **[M]**.

**The cap is the model, not a safety valve.** X reaches **2.6e7** within the first minute
of the fiducial at Γ=1e8, 53 daughter bins — the injection front puts a 1e15 stencil ratio
inside one cell, and a geometric mean is genuinely that sensitive to its smallest member.
Effectively uncapped (1e12) a single-k solve goes from 50 s to **over 12 minutes, and never
finished**. It is also physically empty: G is suppressed by the same factor, so the node
contributes nothing to Λ. The cap converges — 2 k-modes, dr_N_q = 53, against cap = 1e5 **[M]**:

| config | wall | d_dr_l | d_dr_phi | d_cdm / d_tot | phi |
|---|---|---|---|---|---|
| linear F gather (`bal`) | 50 s | 2.02e-02 | 3.24e-02 | 3.0e-05 | 1.2e-06 |
| chain, cap 1e3 | 65 s | 1.9e-05 | 3.6e-04 | 8.2e-07 | 3.3e-08 |
| chain, cap 1e4 (default) | 74 s | 5.6e-06 | 2.9e-04 | 3.8e-07 | 1.6e-08 |
| chain, cap 1e5 | 154 s | — | — | — | — |

`chain_cap` defaults to 1e4 and is read wherever `ChainFactor` is called, which is
`chain_diag_` -- so it binds the **background diagonal** as well, on any kernel with
`balanced_gather` and a single grid (`dr_bg_refine = 1`, the production setting).

⚠ **Corrected 2026-08-28.** This section used to say the cap was read ONLY under
`balanced_pert`, and that on a `balanced_gather`-only kernel `CollisionDiagonal` is the
true derivative of the background RHS which `inv_cap = 0` leaves bit-for-bit uncapped.
That was the shipped behaviour and it was a bug: X diverges as a bin empties, so the
uncapped diagonal reached **6.8e84** and `etd` could not integrate it --
`hpc_ratchet`'s `cmb_G05.0_m0.3_x_bal_q2x` died in 3 s. Capping bends nothing physical,
because `X_e f_e = G(1-+G)/(1-+f_e)` is bounded and a perturbation of an empty bin is
itself empty; measured, capped vs uncapped background at `dr_N_q` = 47/71/101 leaves
rho_ur and H identical and moves the extinct (~1e-32) parent density by <= 6e-6 **[M]**.

**Conservation is untouched**, because only the GATHER changed — the same separation
`balanced_gather` relies on. The l = 0 number identities hold because all three legs are
built from one node source term whatever formed it, and the q-weighted ones depend only on
Σ w_e q_e = q\*. Asserted on a deliberate front state at 0.0 / 1.8e-16 **[M]**.

**What it buys, and what it does not.** The daughter perturbations converge markedly
faster in the daughter grid (nested 53 → 105, so every coarse point is a fine point) **[M]**:

| scheme | d_dr_l | d_dr_phi | d_cdm |
|---|---|---|---|
| linear F gather | 1.35e-02 | 2.08e-02 | 1.59e-03 |
| chain rule | **5.67e-03** | **5.69e-03** | 1.61e-03 |

— 2.4× on the fermion and 3.7× on the boson. But the OBSERVABLES do not move: against the
finest chain run, the linear F gather is off by 3.0e-05 in d_cdm/d_tot and 1.2e-06 in phi,
where the `balanced_gather` background step is 1.8e-02 — 600× larger. And the observables'
own daughter-grid convergence is **unchanged** (1.59e-03 against 1.61e-03), so the chain
rule does not let anyone coarsen `dr_N_q`.

So: a real defect, a real bug, a correct fix at ~1.5× — and on this evidence it does not
pay for itself on observables at Γ=1e8. It is off by default for that reason. Whether it
matters at Γ ≳ 1e10, where the daughters carry far more of the budget, is **[?]**; a
perturbation run there is ~19 h and was not attempted.

⚠ Scope: one configuration (m = 0.3, Γ = 1e8, dr_N_q = 53/105), two k-modes, mTk at z = 0.
Not C_ℓ, not P(k) over a full k-range, not a second Γ.


## 6. What is deliberately not done

* ~~**The perturbation chain rule.**~~ **DONE 2026-08-28 — see §5b.** The text below is
  kept because its diagnosis was half right (the divergence is real) and half wrong (the
  limiter it asks for is not what the fix needed). `ApplyPerturbationOperator` still gathers F linearly,
  so it is the exact Jacobian of the *linear* gather, not of this one. Making it match
  needs dG/df_e = w_e G/f_e, which diverges as a bin empties — a real property of the
  geometric mean, not a coding slip — and needs a limiter on G/f_e that nobody has
  designed or measured **[?]**. The background band factor is where the Γ-amplified error
  was measured to live; a Γ=1e7 single-k run is fine without it (rc = 0, P(k) moves
  0.055% max / 0.047% median against the shipped scheme) **[M]**.
* **`PositivityDecomposition`** is left describing the linear-gather decomposition and
  says so in its doc. Under `balanced_gather` the band factor is not affine in the bin's
  own occupation, so `gain + diag·f` is a tangent line rather than an identity. It has no
  consumer; extend it before giving it one.
* **Bit-identity of the linear path.** Not preserved: ρ / number / pressure move 1.5e-12
  on a background table, round-off the ODE amplifies. Keeping the linear gather lambda
  byte-for-byte and allocating the new scratch only under the flag were both tried and
  neither recovers the md5 — merely putting a branch in that loop reschedules the
  surrounding reductions. Verify the OFF path on the physics, not on a checksum.
* **The Γ=1e9 limit itself.** `bal+exact` converges on 1.2233–1.2234, which sits ~0.8%
  below the 1.233 that conv2 §4's universality argument predicts; that prediction was
  never measured at 1e9 **[?]**. It is either a small real Γ-dependence of the saturated
  limit or a residual bias, and the two are distinguishable by running m = 0.06 (which
  saturates at Γ ≈ 4.6e4) on the fixed scheme. Not done here.
