# The reduced DNCDM operator as a Galerkin congruence — results

**Status:** built and measured, 2026-08-12. Supersedes §2.1 and §2.3–2.5 of
`2026-08-12-dncdm-reduced-perturbation-model-design.md`; §0, §1, §2.2, §2.6, §2.7,
§2.8 and §3–§6 of that note stand.
**Code:** `tools/reduced_collision_operator.{h,cpp}`, `tools/reduced_operator_test.cpp`
(`make test-reduced-operator`, 0.4 s, no cosmology). All three files still present, and
the target is registered in both `CMakeLists.txt` and Makefile `TEST_TARGETS` — verified
2026-08-22.

**The deflation was settled after this note, and against it (added 2026-08-22).** §7 of
the design note is the one section the header above does not rule on, and §5 below
deliberately declines to settle it: it reports the deflation worth ~2x rather than the
~20x M13 projected, and says outright that this is "not a verdict on the method" because
the Γ = 10^11 row had not run. It was settled five days later.
`2026-08-17-dncdm-minimal-knob-set.md` §4 recommends removing the deflation entirely on
the Γ = 10^8–10^11 scan — switching it off costs 5.6e-4 to 4.3e-3, confined to a few low
ℓ and ~100x below cosmic variance there, while saving 42–65% of every reduced cell — and
the code agrees: `deflate_collision` and the whole shift-invert/inverse-iteration
machinery are gone. **Read §5 as the partial evidence it says it is, then go to that
note.**

---

## 0. The one-paragraph version

The reduced operator does not need to be derived. It can be *assembled from the shipped
kernel* as `M~ = P M R`, and doing so settles the structural half of the design outright:
conservation, detailed balance and dissipativity are all inherited rather than re-proved.
The ℓ=0 spectrum is a machine zero at every resolution and every number of daughter
degrees of freedom; the ℓ=1 abscissa is set by the *quadrature* grid's momentum residual,
falls as `N^-3.85`, and reaches `1.4e-3 · aH` while the integrator's state stays at 20
variables. **The design's load-bearing basis, as written in §2.1, is wrong**, and the fix
is one weight function.

---

## 1. The construction, and why it replaces §2.3–2.5

§2.3–2.5 propose deriving the reduced operator: Legendre transfer functions `G_ℓ`/`F_ℓ`
for the decay source, band-averaged quantum-statistics weights, relaxation matrices for
the inverse-decay legs. That is a second implementation of the physics, and it can
disagree with the shipped kernel for reasons that have nothing to do with the closure
under test.

Instead, define a reconstruction `R` from daughter moments to a daughter grid function
and the moment map `P` back, and take

```
    M~_ℓ  =  P · M_ℓ · R
```

with `M_ℓ` the operator `DecayTransitionKernel::ApplyPerturbationOperator` already
applies. Nothing is re-derived. Three properties follow, and they are the whole design:

**Conservation is inherited.** The sector's conserved functionals are moments: number is
`∫dq q² δf`, energy and (at ℓ=1) momentum are `∫dq q³ δf`. Those *are* the retained
coordinates, so `wᵀP = wᵀ` exactly and whatever conservation `M` has, `M~` has — for any
daughter DOF count ≥ 2. No transfer-function identity is needed; `F_1(x) = x` is a
consequence, not an ingredient.

**The quadrature grid decouples from the state dimension.** `M_ℓ` is assembled on a grid
that exists only inside the assembly; the ODE state is the moment vector. The lumped
daughter loss misplaces momentum by O(dq), and that residual is what drives the spurious
positive ℓ=1 eigenvalue (M10/M11, `a8688d02`) — so it can be driven down by refining a
grid the integrator never sees. **This, not an analytic identity, is the mechanism that
breaks the `N ∝ Γ^0.25` rule.**

**Dissipativity is inherited if `P` and `R` are adjoint.** In the entropy inner product
`⟨δf, δf'⟩ = Σ dq q² δf δf' / w`, with `w = f̄(1−f̄)` for a fermion and `f̄(1+f̄)` for a
boson, `P = R†` makes `M~` a congruence, so its Rayleigh quotient at `z` equals `M`'s at
`Rz` and the numerical range is a strict subset:

```
    max Re λ(M~)  ≤  max Re λ(M)      identically, at every background and resolution
```

That inequality — not "the reduced bound is negative" — is the correct statement of the
milestone. The reduction cannot manufacture a growing mode.

### The basis, and the §2.1 correction

`R` reconstructs `δf_i(q) = w_i(q) · Σ_j c_j (q/q_ref)^j`. With the entropy weight and
two moments the span is `{w, w·q}`, which is **exactly** the space of equilibrium
perturbations: `δf = (∂f̄/∂μ)δμ + (∂f̄/∂T)δT` is `w/T` and `w(q−μ)/T²`. So §2.6's
requirement — that the operator annihilate equilibrium perturbations identically, because
`Λ` is a small residual of large terms (`f̄_φ ≈ 140`) — holds by construction rather than
by cancellation.

§2.1 specifies the same two moments but an **unweighted** reconstruction. Measured
below, that basis misses an equilibrium perturbation by a factor of 560–1200. The moments
were right; the reconstruction was not.

---

## 2. What was measured

Kernel units: `Γ = 33.36/Mpc` (an ini `Gamma = 1e7`), `m = M/T = 356.74` (0.06 eV),
`a = 9e-3` — inside the collision window per M1. Production kernel config
(`lump_limiter`, `stratified_quadrature`, `n_gauss = MinGaussNodes`). Parent grid is the
real one, `quadrature_strategy = 3` with 16 bins. Reduced state: **20 variables**
(16 parent + 2×2 daughter moments) against the exact operator's 68–418.

### R1 — the harness is self-consistent

`max |P R − I| = 2.3e-12` over three weights × `n_moments` 1..4.

### R2 — the assembly reproduces the exact operator on its own subspace

Driving both with a `δf` that lies in the reduced span, over ℓ = 0..4:

```
    max relative | M~ z  −  P M R z |  =  6.3e-16
```

So the closure error is *entirely* the projection of out-of-span `δf`, and nothing else.
This is design §4 test 2, run with no free constant — the thing every previous attempt
(a fitted q-weight, the 3×3 fluid RTA) could not do honestly.

### R3 — detailed balance: the §2.1 basis fails, the entropy basis is exact

On a background in exact detailed balance (`μ_H = μ_l + μ_φ`, one temperature,
`balanced_gather = yes`). Relative reconstruction error of an equilibrium perturbation in
the physical entropy norm, and the identity `M~(P v_eq) = P(M v_eq)`:

| basis | nm | ‖RP v−v‖/‖v‖ (ν_l) | (φ) | M~Pv − PMv |
|---|---|---|---|---|
| **entropy** | 2 | **8.7e-16** | **5.6e-16** | **1.2e-13** |
| **entropy** | 3 | 5.4e-15 | 3.1e-15 | 1.5e-12 |
| occupation | 2 | 7.5e-02 | 1.3e-01 | 9.6e+00 |
| occupation | 3 | 3.9e-02 | 1.0e-01 | 4.8e+00 |
| **flat (§2.1)** | 2 | **5.6e+02** | **7.8e+02** | 2.2e+02 |
| **flat (§2.1)** | 3 | 7.6e+02 | 1.2e+03 | 1.3e+02 |

The exact operator's own residual on `v_eq` is 1.2e-5 / 1.3e-7 / 5.0e-4 for the two μ and
the T direction — *not* machine zero, because the two-bin gather leaves the discrete `Λ` a
small residual even under `balanced_gather`. No reduction can beat that floor, and the
entropy basis sits on it.

### R4 — conservation, and it is better than the exact operator's

At `a = 9e-3`, `N_q = 101`, peak `f̄_φ = 104`, `f̄_l = 0.56`, `aH = 2.39e-4/Mpc`. Read
number/energy at ℓ=0 and momentum at ℓ=1 (the others are not identities):

| operator | ℓ=0 number | ℓ=0 energy | ℓ=1 momentum |
|---|---|---|---|
| exact (218 vars) | 6.6e-16 | 1.8e-03 | 1.9e-03 |
| reduced (20 vars) | 3.5e-16 | **2.1e-07** | **9.9e-08** |

Four orders better, and not by magic: the reduced operator only ever acts on `δf` in the
span of `R`, and the badly-conserving directions — a `δf` concentrated in one daughter bin
— are not representable. The reduced *dynamics* therefore conserves better than the exact
discrete dynamics it is derived from.

### R5 — numerical-range containment holds; the flat basis violates it

`max Re λ(M~) ≤ max Re λ(M)` held at every `n_moments` ∈ {2,3,4} and ℓ ∈ {0,1,2} for the
entropy basis, and was violated by 3–4 orders at every row for the flat basis
(`+4.7e+01` against the exact `+2.5e-02`). The test has power.

### R6 — THE MILESTONE

`max Re λ` by shifted QR on the metric-transformed operator (pinned first against a
companion matrix with roots {1,2,3,−4}). **The state is 20 variables in every row; only
the quadrature grid moves.**

| N_q | state | exact ℓ=0 | exact ℓ=1 | **red ℓ=0** | **red ℓ=1** | red ℓ=1 / aH |
|---|---|---|---|---|---|---|
| 26 | 20/68 | +7.70e-04 | +9.02e-04 | +6.3e-17 | +8.65e-04 | 3.6e+00 |
| 51 | 20/118 | +4.08e-03 | +4.66e-03 | +8.3e-17 | +7.13e-05 | 3.0e-01 |
| 101 | 20/218 | +5.95e-03 | +6.40e-03 | −1.6e-17 | +4.98e-06 | 2.1e-02 |
| 201 | 20/418 | +6.76e-03 | +7.04e-03 | +1.3e-17 | **+3.26e-07** | **1.4e-03** |

* **ℓ=0 is a machine zero at every grid.** The number and energy functionals are retained
  coordinates, so their conservation zeros survive the reduction exactly.
* **ℓ=1 falls as `N^-3.85`** — the same exponent the campaign measured for the exact
  operator's spurious eigenvalue (M11's `N^-3.93`) — tracking its own momentum residual,
  which falls 1.7e-5 → 6.4e-9 over the same sweep.
* At `N_q = 201` the reduced abscissa is `1.4e-3` of `aH`: growth of 0.1% over a Hubble
  time.

> **⚠ Retraction.** An earlier version of this note read the ratio of the two columns and
> claimed the reduced abscissa was "21 600× below the exact operator's". **That number is
> void.** It divided by the exact column, which §3 already flagged as anomalous on this
> analytic background — so the ratio measured the artefact, not the reduction. §7 below
> repeats the sweep on real states, where the exact column behaves and the two columns come
> out **comparable**. What the reduction buys is not a smaller eigenvalue at equal grid; it
> is the same eigenvalue at 20 state variables instead of 418.

And against the DOF count at `N_q = 201`:

| nm | state | red ℓ=0 | red ℓ=1 |
|---|---|---|---|
| 2 | 20 | +1.3e-17 | +3.26e-07 |
| 3 | 22 | +3.1e-17 | **−8.3e-09** |
| 4 | 24 | +1.7e-16 | −2.4e-07 |
| 5 | 26 | −1.7e-16 | −2.8e-07 |

`max Re λ ≤ 0` for **any** number of daughter degrees of freedom, which is the load-bearing
half: a reduction whose stability depended on the moment count would merely have moved the
Γ ratchet from the grid onto the moments.

---

## 3. What this does NOT settle

**The accuracy of the closure on a real trajectory.** Everything above is algebra —
conservation, adjointness, congruence, spectra — and holds for arbitrary positive `f̄`,
which is the standard `decay_kernel_test` holds its own identities to. The background here
is analytic (Pauli-saturated `ν_l`, `f̄_φ ≈ 100`, an injection front, 15% out of detailed
balance). It is representative in structure, not a trajectory. **M7 and M8's question — how
large is the error when the daughters' real `δf` is projected onto two moments — is
untouched**, and it is the remaining risk.

**One flag, stated because it is contrary.** The *exact* operator's abscissa RISES with
`N_q` here (7.7e-4 → 6.8e-3), where M11 measured `N^-3.9` on real states. The likely cause
is this background's `q^-6` front, which a finer grid resolves into steeper gradients, plus
the metric floor on unpopulated bins. It does not touch the reduced-operator results, which
are structural, but it says the exact-operator column of R6 should not be quoted, and it is
a second reason to repeat the sweep on states from a CLASS run.

**Cost, on paper only.** The assembly is `size()` = 20 operator applications per background
row, against the deflation precompute's 118–602 — so the reduced operator's precompute is
*cheaper than the deflation it would replace*. Storage is `(l_max+1) × n_rows × 20²`
doubles ≈ 170 MB at `l_max = 17` and 3000 active rows, reducible because most blocks are
ℓ-independent (the parent loss and both daughter self-losses carry no angular factor).

---

## 4. M8, re-run in the entropy inner product — the accuracy question

`CLASS_DNCDM_CLOSURE="a_lo,a_hi[,n_rungs][,path]"` (`species/dncdm_inv_species.cpp`) runs
M8 with the basis weight as a variable, on **real states**: it takes the daughters' actual
`δf` out of `y` mid-run, replaces each by its projection onto `n_moments` basis functions,
keeps the parent exact, applies the same kernel, and emits the absolute ℓ-moment rates for
both. One rung per log-spaced `a` across the window, claimed by whichever k-mode reaches it
first. Driven here by a Γ=1e7, m=0.06 run — M8's own cell — 16 rungs over `a = 1e-3…4e-2`.

**Absolute moments are emitted, not ratios.** M8's pointwise relative error is an artefact
on rates that nearly cancel (the ℓ=0 parent number rate) or vanish outright (the parent,
once extinct): the probe's raw dE_H reaches 8.2 at the last rung purely because the
denominator is going to zero — the same trap that makes `d_ncdm_dncdm1` reach 1e12 in
healthy runs. Normalised below against the window's own RMS of the exact rate.

### 4.1 Like-for-like with M8 (pointwise relative median, ℓ=2)

| basis | DOF | dE_H | dE_l | dE_φ |
|---|---|---|---|---|
| **M8's, quoted** | 1 | 9.8 | 1.8 | 15.5 |
| **M8's, quoted** | 2 | 7.8 | 2.6 | 12.0 |
| **M8's, quoted** | 32 | 0.9 | 0.7 | 0.4 |
| **M8's, quoted** | 134 | 0.48 | 0.19 | 0.06 |
| flat (this probe) | 1 | 12.6 | 5.7 | 4.7 |
| flat (this probe) | 2 | 13.6 | 7.0 | 8.1 |
| **entropy** | 2 | **0.89** | **0.16** | **0.64** |
| **entropy** | 3 | 0.21 | 0.072 | 0.14 |
| **entropy** | 4 | **0.045** | **0.035** | **0.068** |

The `flat` row reproduces M8's magnitudes at the same DOF, which is the control that makes
the comparison mean something: **M8 measured its projection, not the operator.** In the
operator's own inner product, **4 moments beat M8's 134-group projection**, and 2 moments
sit between its 32- and 100-group rows. That is 30–100× in DOF efficiency.

### 4.2 The honest metric (scale-relative to the window RMS, ℓ=2)

Median / p90 / max of `|E' − E|` over the window, in units of the RMS of the exact rate:

| basis | DOF | dE_l med/p90/max | dE_φ med/p90/max | dE_H med/p90/max |
|---|---|---|---|---|
| entropy | 1 | 0.94 / 1.27 / 1.37 | 2.57 / 13.4 / 15.1 | 2.06 / 8.3 / 17.5 |
| **entropy** | **2** | 0.074 / 0.31 / 0.42 | 0.40 / **4.5** / 6.7 | 0.65 / **4.8** / 6.6 |
| **entropy** | **4** | 0.026 / 0.12 / 0.15 | 0.19 / 1.07 / 2.4 | 0.11 / 0.72 / 1.1 |
| entropy | 8 | 0.017 / 0.063 / 0.14 | 0.033 / 0.14 / 0.40 | 0.089 / 0.39 / 0.51 |
| occupation | 2 | 0.19 / 0.65 / 0.72 | 1.18 / 3.8 / 7.7 | 0.70 / 7.2 / 8.5 |
| flat | 2 | 3.0 / 8.1 / 11.3 | 9.4 / 22 / 30 | 10.6 / 24.6 / 33.5 |

**This is the number that decides the moment count, and it is not 2.** The medians at
nm = 2 are good (7–65%), but the *tail* is not: p90 reaches 4.5–4.8× the window RMS. Four
moments bring p90 to ~1×, eight to ~0.4×. The entropy weight is worth 15–25× over the flat
basis at every DOF, but it does not make two moments sufficient.

**Revised reading of M8.** Its verdict — "the two-moment design of §2.1 is dead by two
orders of magnitude" — is overturned as to the two orders; the gap is a factor of a few, and
it closes at 4–8 moments. Even at 8 moments per daughter the state is 2×8×18 = 288 daughter
variables against the shipped 7524, so the cost case survives intact.

**Caveats.** One k (0.0056), one cell (Γ=1e7, m=0.06), 16 rungs, parent kept exact
throughout (as M8 did). The monomial basis is the limit above nm ≈ 8 — the Gram matrix runs
out of conditioning and `SetBackground` throws, which the probe skips rather than papers
over. Higher DOF needs orthogonal polynomials in the weighted measure.

## 5. The deflation matrix, read against a reference

Four of 35 cells have landed (Γ=1e8, m=0.3), and `ref` is not among them. `hpc_cmb`'s
`G1e8_m0.3_r1` is the same ini apart from `dr_N_q` 149 vs 153 and the absent (defaulted)
`deflate_collision` key, so it stands in — a 2.7% grid difference, not the same file.
Harness scale-relative metric:

| variant | pts | defl | TT max | TT med | P(k) max |
|---|---|---|---|---|---|
| c4off | 39 | no | 1.74e-01 | 1.77e-03 | 1.00e-02 |
| c4on | 39 | **yes** | 6.22e-02 | 1.78e-03 | 6.93e-03 |
| c8on | 20 | yes | 4.21e-01 | 3.20e-03 | 2.66e-02 |
| **c2on** | **77** | **yes** | **6.44e-03** | **3.31e-05** | **4.56e-04** |
| r2 | ~75 | no | 1.12e-02 | 5.87e-05 | — |

The controlled comparison is **c2on vs r2** — same resolution, deflation on vs off:
**1.7× on TT max and 1.8× on the median.** At refine = 4 deflation is worth 2.8× on TT max
and *nothing* on the median (1.77 vs 1.78e-3), i.e. it removes the low-ℓ runaway and leaves
the bulk error, which at this Γ is discretisation-limited.

So the deflation is real but is worth **~2×, not the ~20×** M13 projected from the Γ=1e7
laptop measurement. **This does not settle it**: λ ∝ Γ, and the row the README says matters
most — Γ = 1e11, where the Γ^0.25 rule asks for 585 background points — has not run. Read
this as "at Γ=1e8 the eigenvalue is not what binds the grid", not as a verdict on the method.

For the reduced model the reading is favourable and specific: if the bulk error at these
grids is discretisation rather than stability, then a scheme that only removes the
eigenvalue caps out around 2×, while one that also decouples the quadrature grid from the
state does not. That is the distinction §1 draws.

## 6. Steps 1–3 of the previous plan, executed

### 6.1 Orthogonal basis (was step 2)

The monomial Gram is a Hankel moment matrix and its Cholesky fails around `n_moments = 8`.
Replaced by orthogonal polynomials in the measure `dq q² w(q)`, built by the **Stieltjes
three-term recurrence**, which never forms `t^j` at all. `psi_0 = 1` and `psi_1 = t − α₀`
by construction, so `span{psi_0, psi_1} = span{1, t}` and the number/energy functionals stay
exactly representable — their *coordinates* change (energy is `q_ref(m_1 + α₀ m_0)`) and
`Conservation()` carries that. The Gram becomes diagonal, so reconstruction is a scaling.

Conditioning improved across the board: `‖PR − I‖` 2.3e-12 → **1.8e-14**, equilibrium
reconstruction 5.4e-15 → **7.1e-16**, spectra unchanged (same span, as required).

**The DOF axis is not unbounded, and the test now measures the ceiling instead of asserting
past it:** `max Re λ ≤ 0` holds cleanly for `n_moments = 2…9` and breaks at 10 (ℓ=0 drifts
to +2.5e-6). High-order orthogonal polynomials oscillate with a large dynamic range, so the
reconstruction and the projection of the operator's output start to cancel. That is an
accuracy ceiling of the reduction, not of the recurrence, and 2–9 covers what §4 needs.

### 6.2 The closure probe, widened (was step 1)

Rungs are now keyed on **(a, k)**, not `a` alone — with an `a`-only key the first k-mode to
reach a rung claimed it and the entire first sweep came out at one k, leaving the
k-dependence unmeasured. Run at two cells, 12 a-rungs × 8 k-bins, `n_moments` 1…9.

Scale-relative **per k** (so a k-dependent amplitude cannot leak into the comparison),
median / p90 of the ℓ=2 energy-moment rate error:

| cell | basis | nm | dE_l | dE_φ | dE_H |
|---|---|---|---|---|---|
| Γ=1e7, m=0.06 | entropy | 2 | 0.035 / 0.18 | 0.10 / 2.08 | 0.18 / 2.53 |
| Γ=1e7, m=0.06 | entropy | 4 | 0.017 / 0.068 | 0.071 / 0.42 | 0.077 / 0.81 |
| Γ=1e7, m=0.06 | entropy | 8 | 0.004 / 0.038 | 0.013 / 0.083 | 0.043 / 0.44 |
| Γ=1e8, m=0.3 | entropy | 2 | 0.039 / 0.17 | 0.12 / 1.51 | 0.079 / 1.67 |
| Γ=1e8, m=0.3 | entropy | 4 | 0.015 / 0.083 | 0.045 / 0.23 | 0.109 / 0.70 |
| Γ=1e8, m=0.3 | entropy | 8 | 0.003 / 0.035 | 0.008 / 0.066 | 0.026 / 0.21 |
| Γ=1e8, m=0.3 | flat | 8 | 0.57 / 20.3 | 1.21 / 71.5 | 1.97 / 133 |

Three things, all new:

* **The closure error does not grow with Γ.** The Γ=1e8 rows match the Γ=1e7 rows at every
  `nm`. The grid requirement scales as `Γ^0.25`; the closure error does not scale at all.
  That is the asymmetry the whole design is trading on, measured directly.
* **The flat basis does not converge at all** — 0.89 at nm=1 and 0.57 at nm=8, with the p90
  *rising*. M8's groups converged slowly (~1/N); plain q-moments in the wrong measure do not
  converge. The inner product is not a refinement, it is the difference between a convergent
  and a non-convergent scheme.
* **The error peaks at intermediate k** — at Γ=1e7 the worst bin is k ≈ 5.6e-3 (nm=2: median
  1.75, p90 3.58) while both k ≤ 1.1e-3 and k = 0.16 are an order better. That is where free
  streaming and the collision rate cross, and it is where the daughters carry the most
  non-equilibrium q-structure. At nm=8 the worst k is down to 0.063 / 0.30.

**Verdict: 4 moments is adequate, 6–8 is comfortable, 2 is not.**

### 6.3 The spectrum on real states (was step 3) — and the flag is settled

`CLASS_DNCDM_CLOSURE_SPECTRUM=1` adds a spectrum row per rung. Swept over `dr_bg_refine`,
which moves only the perturbation grid — the background is identical in every row, so this
is the same physical state on every grid. Γ=1e7, m=0.06, at `a = 1.125e-2`:

| refine | N_pt | N | exact ℓ=0 | exact ℓ=1 | red ℓ=0 | red ℓ=1 | red ℓ=1 mom |
|---|---|---|---|---|---|---|---|
| 8 | 26 | 68 | 7.9e-17 | 1.342e-03 | −5.0e-18 | 2.347e-03 | 3.2e-05 |
| 4 | 51 | 118 | 1.2e-17 | 1.842e-04 | −8.5e-17 | 2.017e-04 | 3.0e-06 |
| 2 | 101 | 218 | 1.3e-17 | 1.429e-05 | −5.0e-18 | 1.356e-05 | 2.0e-07 |
| 1 | 201 | 418 | 1.6e-17 | 9.020e-07 | 7.2e-17 | 9.194e-07 | 1.4e-08 |

* **The `N^-3.9` law reproduces on real states**: slope 3.57 over 26 → 201, and the
  201-point value 9.0e-7 sits right on M11's measured +4.14e-7. §3's flag is **settled** —
  the rising-abscissa anomaly was the analytic background's `q^-6` front, nothing more.
* **ℓ=0 is a machine zero for the exact operator too**, at every grid. At Γ=1e7 in this
  window the ℓ=0 mode the deflation also targets simply is not there, consistent with
  `a8688d02`'s note that ℓ=0 only misbehaves late in the window.
* **The reduced and exact ℓ=1 abscissae are comparable at equal grid** — this is the
  retraction in §2. The reduced value is bounded by its own momentum residual, as designed,
  and that residual is the same discretisation the exact operator pays.

**What survives, stated exactly.** The reduced operator reaches the ℓ=1 abscissa of a
201-point grid (9.2e-7) while carrying **20 state variables**; the exact operator needs
**418** for the same number, and at 68 state variables it is 1500× worse. That is a ~20×
state reduction at equal stability — the design's actual claim, and it holds.

**But the Γ scaling needs restating too.** The reduced abscissa is set by its quadrature
grid's momentum residual, and that residual still obeys `λ ∝ Γ·N_q^-3.9`. So `N_q` must
still grow as `Γ^0.25` — what changes is *what that costs*. Refining `N_q` now costs only
the per-background-row precompute, which is k-independent and amortises over ~1500 k-modes,
instead of multiplying the ODE state by `l_max` and by every k. **The Γ ratchet does not
disappear; it moves off the integrator and onto a precompute.** That is still the win, but
it is a weaker statement than "the grid requirement stops scaling with Γ", and the earlier
note said the stronger thing.

## 7. The reduced species — started, and where it stands

### 7.1 The basis has to be frozen, and freezing was measured before it was assumed

`psi_j` is built from `w = f̄(1∓f̄)`, which evolves. A basis that tracked `a` would make the
state variables mean something different at every step, and `dm_j/dτ` would acquire an
integral of `(∂psi_j/∂τ)δf` that the free-streaming hierarchy does not contain. So a
shippable reduced species must freeze the basis — which the closure measurements of §6.2 did
*not* test, since they rebuilt it at every rung.

`CLASS_DNCDM_CLOSURE_FREEZE=1` measures the cost. Γ=1e7, ℓ=2, scale-relative per k,
median / p90, basis frozen at the first rung (a = 1e-3) and reused across the window:

| nm | dE_l live | dE_l frozen | dE_φ live | dE_φ frozen |
|---|---|---|---|---|
| 2 | 0.035 / 0.18 | 0.170 / 0.46 | 0.102 / 2.08 | 0.378 / 2.55 |
| 4 | 0.017 / 0.068 | 0.041 / 0.10 | 0.071 / 0.42 | 0.045 / 0.32 |
| 6 | 0.009 / 0.051 | 0.014 / 0.082 | 0.016 / 0.21 | 0.046 / 0.29 |
| 8 | 0.004 / 0.038 | 0.010 / 0.067 | 0.013 / 0.083 | 0.027 / 0.18 |

Freezing costs a factor ~2.5 at nm = 2 and ~1.2–1.5 at nm ≥ 4. **Frozen at 6 moments is as
accurate as live at 4**, which is the direct justification for the shipped default of 6.

### 7.2 What is implemented

`DrPsdSpecies` gained a reduced-moment representation (`SetReducedBasis`, `reduced()`),
branching in `RegisterPerturbationIndices`, `PerturbDerivs`, `ApplyInitialConditions`,
`StressEnergy` and `FillSources`. It compiles and the whole existing suite is green.

It is a representation swap, not new physics, for two exact reasons:

* **Free streaming is unchanged.** `ε = q` makes `qk/ε` exactly `k`, carrying no q, so every
  q-moment obeys the same hierarchy with no closure error (M3). Only the metric driver
  changes, from the pointwise `∂f̄/∂lnq` to its projection `D_j`.
* **The metric only ever reads the energy moment.** With `ε = q` all four stress-energy sums
  carry the weight `q³`, so `δρ`, `θ` and `σ` are the energy moment at ℓ = 0, 1, 2 and
  nothing else. In moment coordinates that is `q_ref(m_1 + α₀ m_0)` — verified against
  `Σ dq q³ δf` to **4.4e-15** for arbitrary (not in-span) `δf`, `reduced_operator_test` q1b.
  That identity is load-bearing and silent when wrong: the hierarchy would still run and
  quietly conserve the wrong quantity.

### 7.3 IT RUNS — and what the observable says

Wired end-to-end: `dr_reduced_moments = N` (0 = off, default), `dr_reduced_weight`
(entropy | occupation | flat), the frozen basis built in `ProcessBackgroundTable`, and the
collision applied **matrix-free** in `ApplyKernelPerturbDerivs` — reconstruct, one ordinary
kernel sweep, project. Γ=1e7, m=0.06, `dr_N_q = 201`, `dr_bg_refine = 4` (51 daughter grid
points), reference = the identical run with the key off:

| n_moments | daughter vars/ℓ | wall | speed-up | max ΔP/P (rtol 3e-6) | max ΔP/P (rtol 3e-8) |
|---|---|---|---|---|---|
| exact | 102 | 609 s | — | — | — |
| 2 | 4 | | | 1.39e-03 | 5.10e-04 |
| 4 | 8 | 161 s | **3.78×** | 1.33e-03 | **4.16e-04** |
| 6 | 12 | 169 s | **3.60×** | 1.39e-03 | **4.37e-04** |
| 8 | 16 | 188 s | 3.24× | 1.28e-03 | 4.18e-04 |

**The rtol column is the point.** At the shipped `tol_perturb_integration = 3e-6` the
discrepancy is 1.4e-3 and is *completely flat* in n_moments (1.39e-3 at 2, 1.28e-3 at 8), in
the freeze point (1.37–1.45e-3 across a factor 20 in `a_freeze`), and in the basis weight
(entropy 1.39e-3, occupation 1.45e-3). Nothing about the reduction moved it. Tightening rtol
100× drops it to **4.4e-4**: most of it was never the reduction at all, but the fact that a
28-variable and a 418-variable system get different step sequences at the same relative
tolerance. That is a comparison artefact, and it would have been reported as a closure error.

**The calibration that makes 4.4e-4 meaningful**: at this same cell the exact scheme's own
daughter-grid convergence scatter — `hpc_cmb` r1 vs r2, 201 → 101 points — is **4.5e-4** in
P(k). So at 6 moments the reduced model sits at the exact scheme's own grid uncertainty,
with 8.5× fewer daughter variables and 3.6× faster.

At rtol 3e-8 an n_moments dependence finally appears — 5.10e-4 at 2, then 4.16 / 4.37 /
4.18e-4 at 4 / 6 / 8 — but it is **small and it saturates immediately**. Read honestly: only
`n_moments = 2` is distinguishable, by ~1e-4; from 4 up the spread (4.16–4.37e-4) is
scatter, not a trend, and sits at the reference's own 4.5e-4 grid uncertainty. **On the
observable, 4 moments already suffice at this cell** — weaker evidence for 6 than §6.2's
offline closure metric suggested. The offline metric is the more sensitive instrument and
should keep setting the default; but this measurement cannot resolve 4 from 8, and saying
otherwise would be reading scatter as signal.

Two further measured facts:

* **The error is transmitted purely through the metric.** ΔP/P for total and for `P_cb`
  agree to every printed digit, and both vanish at low k (2.7e-8 at k = 1e-5) and peak at
  k ≈ 0.11. It is the daughters' free-streaming shear feeding the metric, not their direct
  density contribution.
* **`flat` diverges in a real run** (ΔP/P = 5.9e+04), exactly as §2's numerical-range
  containment predicted — it violates the bound by 3–4 orders and breaks equilibrium. The
  offline dissipativity test has real predictive power, which is worth knowing.

### 7.4 The simplification found while scoping the collision

The collision term is not wired. Scoping it produced a design change worth recording:

**The reduced collision should be applied MATRIX-FREE, not from a precomputed table.**
`dF = P·M·R·z` needs one reconstruction onto the quadrature grid, **one** ordinary kernel
sweep — the same call the exact code already makes — and one projection back. Assembling
`M̃_ℓ` per RHS instead would cost `size()` ≈ 28 kernel applications, i.e. 28× *more* than the
exact code; and precomputing it over the background table costs `(l_max+1) × n_rows ×
size()²` ≈ 340 MB at `l_max = 17`, which needs an interpolation table in `ln a` to be
affordable.

This relocates the claimed saving and §3's cost table should be read accordingly: **the
per-RHS cost is essentially unchanged** (the kernel still sweeps the quadrature grid), and
the win is that the ODE carries `16 + 2×6 = 28` variables per multipole instead of 418 —
fewer steps, a far cheaper Jacobian, and §2.8's stiff-solver options back on the table. The
precomputed-table route remains available as an optimisation once the matrix-free version is
validated, and only then is the collision itself cheaper.

## 8. What to do next, in order

1. **C_ℓ, not just P(k).** Every observable number above is P(k) from an `output = mPk`
   run. C_ℓ^TT is ~6× more sensitive than P(k) on the r1/r2 anchor (M12's caveat), so one
   lensed-C_ℓ run at Γ=1e7 is the next measurement, and it may move the 4-vs-6 verdict.
2. **A second cell.** Everything in §7.3 is Γ=1e7, m=0.06. Γ=1e8/m=0.3 has a converged
   `hpc_defl` reference on disk to compare against, and the closure work (§6.2) says the
   error should not grow with Γ — that prediction is now cheap to test end-to-end.
3. **Fix the comparison protocol, not just the tolerance.** rtol 3e-8 is a diagnostic
   crutch: it makes the two systems' step sequences agree well enough to see the closure
   error, at 4× the cost. A production comparison wants the reduced run at its own sensible
   tolerance, which means the error budget has to separate integration from closure — the
   `dr_bg_refine`-style convergence rung, not a tolerance ladder.
4. **The Jacobian diagonal is a guess.** `ApplyKernelPerturbDiagonal` hands each moment the
   Rayleigh quotient of the per-bin collision diagonal. That is a defensible
   integrating-factor choice and ETD stays convergent for any diagonal, but nobody has
   measured whether a better one buys step size. This is the first place to look for the
   rest of the speed-up.
5. **Then** §2.7's window switch and §2.8's integrator question, which is where a
   28-variable system is supposed to pay off properly — the 3.6× measured here is from the
   ODE dimension alone, with the collision still sweeping the full quadrature grid.

The honest status: **the reduced species runs**, reproduces P(k) to the exact scheme's own
grid uncertainty at 4–6 moments, and is 3.6× faster with 8.5× fewer daughter variables — on
one cell, one observable, and only once the integration tolerance is tightened enough to
stop hiding the signal.
