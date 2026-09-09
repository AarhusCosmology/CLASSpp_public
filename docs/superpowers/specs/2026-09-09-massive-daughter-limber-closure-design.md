# Which continuously-sourced hierarchies want the Limber closure, and which do not

**Status (2026-09-09).** Shipped: the Limber closure in
`WdmDecayProductSpecies::PerturbDerivs` (the `dcdm_wdm` massive daughter),
replacing the Ma & Bertschinger truncation outright, no knob — and the two
closure helpers moved out of `DarkRadiationSpecies` into `species/limber_closure.h`,
where three species now share them. **Deliberately NOT shipped: any change to
`DrPsdSpecies`** (`dr_representation = psd`, both branches). That was measured
here and the mechanism does not apply to it; §4 is the evidence and §4.4 is the
one number that argues the other way, recorded rather than buried.

Resolves issue #422, the follow-up #420 left open. Every number tagged **[M]** was
measured for this note on 2026-09-08/09 on one binary carrying both closures
behind an env switch (built in ONE build directory, so the arms differ by the
closure and nothing else); **[C]** is read off the code; **[P]** is from an
earlier round and not re-measured; **[?]** has no evidence behind it.

## 1 · The question

#420 replaced MB's truncation in `DarkRadiationSpecies` after showing that MB's
premise — one coherent mode released at `tau = 0` — fails for a species injected
continuously, which is a superposition of shells free-streaming at their own
ages. #422 named two hierarchies with the same shape and the same truncation that
were never measured:

* `DrPsdSpecies` — `dr_representation = psd`, both the q-resolved branch and the
  reduced (Galerkin) one.
* the `dcdm_wdm` massive daughter, `WdmDecayProductSpecies`.

The answers are **no** and **yes** respectively, and the reason they differ is the
useful part of this note.

## 2 · The test, and what it reads on species whose answer is already known

Two diagnostics, both from #422, on a single-k run (`k = 0.1/Mpc`, flat,
`tau0 = 1.4e4 Mpc`) with every fluid/streaming approximation disabled so each
hierarchy is carried to its own `l_max`:

**(a) Ladder shape.** `F_l / F_0` at the last output time, and the number of sign
changes across the ladder. A superposition of shells is monotone and
single-signed; a coherent free-streaming field alternates.

**(b) Closure inversion.** From a run at `l_max = 40..60`, ask each closure what
`F_{L+1}` it would have asserted at a *lower* `L`, from that run's own
`F_{L-1}, F_L`, and compare against the truth the run actually carries:

```
F_{L+1}^MB     = (2L+1) (eps/q) cotKgen F_L - F_{L-1}      cotKgen = 1/(k tau), flat
F_{L+1}^Limber = (2L+1) Cot(L)            F_L - F_{L-1}    Cot(L)  = 1/L(l), flat
```

reported as the median over the second half of the run of
`|pred - true| / max(|F_{L-1}|, |F_L|, |F_{L+1}|)` — scale-relative, because these
ladders cross zero and a max-relative metric would read noise.

Calibration on the two species whose answer is known: **[M]**

| species | ladder sign changes | MB residual | Limber residual | correct closure |
|---|---|---|---|---|
| `dcdm_dr` DR daughter (fixed by #420) | **0 / 17** | 1.905 | **0.0003** | Limber |
| `ur` (coherent free-streaming) | 9 / 17 | **0.555** | 1.276 | MB |

MB reads `-1.034` on `dcdm_dr`, reproducing the "rock-steady -1.05, right
magnitude and wrong sign" of #422. The two diagnostics agree with each other and
with #420 on both controls, so they are the instrument for §3 and §4.

**What the inversion cannot do, and why the two sections weigh it differently.**
The residual saturates near 1 wherever the ladder is *evanescent* — below its own
Bessel turning point, `k chi < l`, `|F_{L+1}| << |F_{L-1}|`, so every prediction is
dominated by its `-F_{L-1}` term and scores the scale itself. That is precisely
the `dcdm_wdm` regime (`q/eps ~ 0.003-0.02` puts `k chi ~ 2-40` against
`l_max = 10`), where it reads MB 0.53-3.98 against Limber 0.37-0.65 — "neither is
any good", which is not an answer. So:

* **§3 decides `dcdm_wdm` on the observable `l_max` scan alone**, which shows the
  #420 signature outright: orders of magnitude, plus MB converging to the Limber
  answer when pushed. The inversion is not used there.
* **§4 decides `dr_psd` on the structural test**, because its observable scan does
  **not** show that signature — it shows a factor of 2 (§4.3, §4.4) — and the
  structural test there is unambiguous and is taken at `k chi ~ 1e3`, three orders
  of magnitude above the turning point, where the diagnostic is in its valid
  regime and does separate the two controls by four orders of magnitude.

The two sections use different instruments because the two species sit in
different regimes, not because each was given the instrument that agreed with it.

## 3 · `dcdm_wdm`: the closure applies, and the transport rate comes with it

### 3.1 The one new term

Written so both closures share a shape, the truncation is

```
dF_L/dtau = (qk/eps) F_{L-1} - D F_L
```

and inverting it against the exact hierarchy gives
`F_{L+1} = -F_{L-1} + D F_L (2L+1)/[(L+1)(qk/eps)]`. So the coefficient of `F_L`
is `(2L+1) C` with `C = D eps/[(L+1) q k]`: **`D` carries the transport rate, `C`
does not.** MB takes `C = eps/(q k tau)` — cot_K at the distance a particle of
*this* momentum has already streamed, `chi = (q/eps) tau` — whose `eps/q` cancels
the transport rate and leaves the mass-blind `D = (1+L) k cotKgen` the code had. **[C]**

The Limber argument is a multipole order and not a distance, so it carries no
`eps/q`, and the transport rate **survives**:

```
D_Limber = (1 + L) (qk/eps) Cot(L)        <- what shipped
```

That factor is the whole of the fix, and it is not cosmetic. Running the naive
port instead — `D = (1+L) k Cot(L)`, the massless `DarkRadiationSpecies` line
copied across — is **worth nothing at all**: **[M]**

| closure at `l_max_ncdm = 10` | TT | P(k) |
|---|---|---|
| MB (before) | 1.03e-05 | 1.03e-02 |
| Limber without `q/eps` (control) | 1.13e-05 | **1.16e-02** — worse than MB |
| Limber with `q/eps` (shipped) | **9.64e-08** | **1.64e-04** |

The control is the point: a species this non-relativistic has `eps/q ~ 50`, so a
closure that gets the argument right and the rate wrong lands back where it
started. `species/limber_closure_test.cpp` pins all three coefficients
algebraically.

### 3.2 Truncation independence

`Gamma = 97.78 km/s/Mpc` (10 Gyr), `vkick = 0.02`, `Omega_ini = 0.026`,
`momenta_bins = 32`, `tCl,pCl,lCl,mPk` + lensing, `l_max_scalars = 2500`, rkdp45
at `tol_perturb_integration = 1e-6`. Each closure against **its own**
`l_max_ncdm = 100`; both observables are positive definite, so the metric is the
max relative difference over the whole spectrum. **[M]**

| `l_max_ncdm` | MB TT | MB P(k) | Limber TT | Limber P(k) |
|---|---|---|---|---|
| **10 (default)** | 1.03e-05 | 1.03e-02 | **9.64e-08** | **1.64e-04** |
| 15 | 2.66e-06 | 5.20e-03 | 1.97e-08 | 6.06e-05 |
| 30 | 2.60e-07 | 2.82e-03 | 7.13e-10 | 6.07e-06 |
| 60 | 2.88e-09 | 9.93e-04 | 1.15e-11 | 7.19e-07 |

MB's **1% P(k) error at the default `l_max_ncdm`** sits at `k > 0.1 h/Mpc`,
peaking at `k = 2.7 h/Mpc` — which is the small-scale suppression this model
exists to predict. Below `k = 0.01 h/Mpc` it is 5e-13. The CMB is barely touched
(TT `l < 30` is 5e-9); the error is a free-streaming error and it shows up where
free streaming does. **[M]**

### 3.3 Robustness across the model

`l_max_ncdm = 10` against 100, per arm, varying one axis at a time: **[M]**

| | MB TT | Limber TT | gain | MB P(k) | Limber P(k) | gain |
|---|---|---|---|---|---|---|
| `vkick = 0.005` | 2.26e-07 | 7.95e-09 | 28x | 7.69e-03 | 1.64e-04 | 47x |
| `vkick = 0.02` | 1.03e-05 | 9.64e-08 | 106x | 1.02e-02 | 1.64e-04 | 62x |
| `vkick = 0.05` | 3.76e-05 | 1.85e-07 | 204x | 1.04e-02 | 1.63e-04 | 64x |
| `vkick = 0.1` | 6.97e-05 | 2.89e-07 | 241x | 1.04e-02 | 1.60e-04 | 65x |
| `vkick = 0.3` | 1.69e-04 | 5.26e-07 | 322x | 1.03e-02 | 1.58e-04 | 65x |
| `Gamma = 9.778` | 1.48e-06 | 1.28e-08 | 115x | 1.81e-03 | 1.63e-05 | 111x |
| `Gamma = 977.8` | 1.72e-05 | 2.36e-07 | 73x | 9.18e-03 | 7.09e-04 | 13x |
| `Gamma = 9778` | 3.41e-06 | 9.81e-08 | 35x | 6.71e-03 | 4.92e-05 | 136x |

The gain never falls below 13x anywhere in the range the model is used, and MB's
P(k) error is ~1% almost everywhere in it.

It is also independent of the momentum grid — the daughter's own resolution knob
is not what the closure is standing in for: **[M]**

| `momenta_bins` | MB TT | Limber TT | MB P(k) | Limber P(k) | gain |
|---|---|---|---|---|---|
| 16 | 1.09e-05 | 9.89e-08 | 1.02e-02 | 9.14e-05 | 110x / 112x |
| 32 | 1.03e-05 | 9.64e-08 | 1.02e-02 | 1.64e-04 | 106x / 62x |
| 96 (the shipped default) | 9.89e-06 | 9.72e-08 | 1.04e-02 | 4.49e-04 | 102x / 23x |

**The TT ratio is floored by the integrator, the P(k) one is not.** Repeating the
`momenta_bins = 32` row under ndf15 (`evolver = 1`) instead of rkdp45: **[M]**

| evolver | MB TT | Limber TT | MB P(k) | Limber P(k) | gain |
|---|---|---|---|---|---|
| rkdp45 @ 1e-6 | 1.03e-05 | 9.64e-08 | 1.02e-02 | 1.64e-04 | 106x / 62x |
| ndf15 @ 1e-6 | 9.64e-06 | 1.32e-06 | 1.02e-02 | 1.67e-04 | **7x** / 61x |

**MB's error is the same under both** (that is the number this change removes), and
the P(k) gain is the same under both. Only the *Limber* arm's TT residual moves,
from 9.6e-08 to 1.3e-06 — and the ndf15-vs-rkdp45 spread on that same spectrum is
**9.1e-06**, i.e. above it. So ndf15's 7x is a lower bound set by its own
integration error, not a measurement of the closure; rkdp45 at
`tol_perturb_integration = 1e-6` can resolve the residual and reads 106x. Quote
the P(k) column when the integrator is not pinned. **[M]**

### 3.4 This is validation, not self-consistency

Each closure converging against its own reference proves only that it converges.
Pushing MB up instead settles which limit they share. `vkick = 0.3` (the arm
where the two still differ most at `l_max = 100`), everything against
**Limber@400**: **[M]**

| arm | `l_max_ncdm` | TT | P(k) |
|---|---|---|---|
| MB | 100 | 2.20e-06 | 1.09e-03 |
| MB | 400 | 3.01e-08 | 2.21e-04 |
| MB | 800 | 6.88e-10 | 9.22e-05 |
| MB | **1600** | **9.74e-13** | **7.50e-08** |
| Limber | 100 | 1.55e-10 | 2.83e-07 |
| Limber | 200 | 1.06e-11 | 3.60e-08 |

MB converges **to the Limber answer**. The two closures describe the same
physics; MB needs `l_max_ncdm` of order 1600 to reach what the Limber closure
reaches at 10. (Richardson off the MB 400/800/1600 triple gives `p = 5.42` and a
TT limit 1.7e-11 from Limber@400; the P(k) sequence is not a clean power law —
`p = 1.58`, extrapolating to 4.6e-05 — so the direct 1600 row is the stronger
statement and the extrapolation is quoted only for TT.) **[M]**

## 4 · `dr_psd`: the mechanism does not apply — no change

### 4.1 Both branches read as coherent, not as a superposition

`dncdm1` with `Gamma = 1e5`, `m = 0.3 eV`, `inverse_decays = yes`,
`balanced_gather = yes`, `emission_gauss = 2`, `momenta_bins = 64`,
`dr_N_q = 100`, `l_max_ncdm = 40`: **[M]**

| branch | ladder sign changes | MB residual | Limber residual |
|---|---|---|---|
| q-resolved (`:560`), 8 bins across the grid | 9–37 / 40 | **0.28 – 1.50** | 1.39 – 1.95 |
| reduced/Galerkin (`:534`), 6 moments x 2 daughters | — | **0.026 – 1.02** | 1.72 – 1.97 |

The reduced branch is the sharper reading. Limber's implied `F_{L+1}` is wrong by
up to **82x** there — the `Limber/true` ratio runs 5.0, 34.3, 22.5, 1.0, -0.55,
-0.72 across the six fermion moments and 81.7, 34.6, 22.1, 1.4, -0.53, -0.72
across the six boson ones, so it is the three lowest-order moments that are badly
wrong and the top three that merely have the wrong sign. MB's residual on the same
moments is as low as 0.026. Neither branch shows anything resembling `dcdm_dr`'s
`0 / 17` and 0.0003.

### 4.2 The discriminator is the representation, not the cosmology

Run the **same** cosmology — same `Gamma`, same mass, same everything — in
`dr_representation = integrated`, whose daughter is a `DarkRadiationSpecies`
already carrying the Limber closure since #420: **[M]**

| same cosmology, different representation | sign changes | MB residual | Limber residual |
|---|---|---|---|
| `integrated` (momentum-integrated daughter) | **1 / 40** | 1.401 | **0.332** |
| `psd` (q-resolved daughter) | 9–37 / 40 | **0.28–1.50** | 1.39–1.95 |

So the two representations of one physical sector want **opposite** closures, and
`dr_representation` is what decides it. The mechanism: the emitted comoving
momentum grows with the scale factor (`q2 ~ a M / 2` for a non-relativistic
parent), so a *fixed comoving-q bin* is fed over a narrow window in `a` and is a
near-coherent shell — MB's premise. The momentum-**integrated** hierarchy sums
those windows and only then becomes the broad superposition #420 is about. **[M]**

That also settles the standalone case for free: with no source `DrPsdSpecies` is
a plain free-streaming massless PSD species (`dr_psd_species.h`), i.e. exactly
the coherent field MB is written for. **[C]**

### 4.3 The observable scan, at the campaign's own settings

Moments rung (`dr_reduced_moments = 6`, what `hpc_gal2`/`hpc_prod2` run), same
cosmology, each arm against its own `l_max_ncdm = 100`: **[M]**

| `l_max_ncdm` | MB TT | MB P(k) | Limber TT | Limber P(k) |
|---|---|---|---|---|
| 8 | 1.54e-03 | 1.94e-02 | 7.04e-04 | 9.58e-03 |
| 10 | 1.84e-03 | 1.17e-02 | 6.09e-04 | 5.67e-03 |
| 15 | 7.04e-04 | 3.85e-03 | 2.12e-04 | 2.16e-03 |
| 20 | 1.59e-04 | 2.14e-03 | 1.27e-04 | 8.99e-04 |
| 30 | 5.09e-06 | 8.16e-04 | 2.37e-05 | 4.10e-04 |
| 60 | 2.49e-07 | 9.25e-05 | 2.10e-07 | 4.90e-05 |

Cross-closure at matched `l_max` falls to TT 3.1e-08 / P(k) 1.1e-05 at 100 —
both closures converge to the same answer, as two consistent closures of one
hierarchy must.

The same scan on the **exact (q-resolved) rung** — the reference method, and the
branch whose per-bin coherence §4.1 measured — reads the same picture, so the
conclusion does not rest on the Galerkin reduction (`dr_N_q = 60`,
`momenta_bins = 32`, each arm against its own `l_max_ncdm = 40`): **[M]**

| `l_max_ncdm` | MB TT | MB P(k) | Limber TT | Limber P(k) | gain |
|---|---|---|---|---|---|
| 10 | 1.80e-03 | 1.18e-02 | 5.50e-04 | 5.97e-03 | 3.3x / 2.0x |
| **17 (the campaign's)** | 2.69e-04 | 2.87e-03 | 9.31e-05 | 1.42e-03 | 2.9x / 2.0x |

with the cross-closure spread falling 1.25e-03 -> 1.87e-04 -> 2.21e-06 on TT
across `l_max = 10, 17, 40`. Same factor of 2-3, same shared limit, on both
branches.

**A finding for the campaign that is not about the closure.** `hpc_gal2` and
`hpc_prod2` run `l_max_ncdm = 17` (all 329 of their inis). Measured directly on
the exact rung at that value, the truncation error there is **2.9e-03 on P(k)**
and 2.7e-04 on TT under MB (1.4e-03 / 9.3e-05 under Limber) — so a few parts in
a thousand on P(k) under **either** closure, and the moments rung brackets the
same. That is a resolution setting, not a closure choice, and it is worth knowing
before those numbers go in a figure. **[M]**

### 4.4 The number that argues the other way, and why it did not win

Limber is nonetheless 2–3x better on the observable at low `l_max`, and its
convergence is monotone where MB's is not (MB reads 1.54e-03 at `l_max = 8` and
1.84e-03 at 10). Non-monotone MB convergence is one of #420's signatures, so this
deserves saying out loud rather than omitting.

It did not carry the decision, for three reasons:

1. **Wrong size.** #420 measured 125x on TT and 500x on P(k) for `dcdm_dr`, and §3
   measures 13–322x for `dcdm_wdm`. A factor of 2 is not that signature; it is
   what two consistent closures of an under-resolved ladder differ by.
2. **Wrong direction on the structural test.** The closure that wins the
   observable by 2x is the one whose asserted `F_{L+1}` is off by 22–82x (§4.1).
   Shipping a per-bin-wrong closure because the *sum* converges faster is a
   trade nobody could defend later, and the ladder is not the only consumer.
3. **The 2x is inside a ~1% error either way** (§4.3). The lever that matters at
   `l_max_ncdm = 17` is `l_max_ncdm`, and it is available to anyone who wants it.

If someone later wants the 2x, the thing to measure first is *why* a per-bin
coherent ladder converges faster under a superposition closure — not to flip the
flag. **[?]**

## 5 · What else was looked at and left alone

* **The `dncdm` parent** (`dncdm_species.cpp:642,663`) is re-sourced by inverse
  decays and so has the surface shape, but it is not a Limber candidate: at
  `q/eps ~ 1e-4` it has `k chi ~ 0.2 – 6` against `l_max = 40`, i.e. entirely
  below its own Bessel turning point, where the ladder is evanescent and neither
  closure's premise is in play. Its high-`l` multipoles do run away there: at
  `a = 1` on an extinct parent, `F_40/F_0` reaches **1e58 to 1e73** depending on
  the momentum bin (e.g. `F_0 = 1.9e-71`, `F_40 = 4.4e-07`) and `|F|` peaks near
  `l = 19..36` rather than at `l = 0`. That is upward-recursion instability on a
  species already known not to be trustworthy once its occupation is on
  `kFParentFloor`, and it is unrelated to this note. It is **not** reproduced at
  the default `l_max_ncdm = 10`, which was checked rather than assumed: there
  `|F|` peaks at `l = 0` and `F_lmax/F_0 = -0.074`. Worth its own issue. **[M]**
* **`s_l[l_max]` is missing** from the `F_{L-1}` term of every NCDM-family
  truncation (`ncdm_base_species.cpp:806`, `ncdm_species.cpp:413`,
  `dncdm_species.cpp:642,663`, `dr_psd_species.cpp:534,560`, and the line this
  note changes), while `photons`, `ultra_relativistic`, `interacting_species` and
  `dark_radiation_species` all carry it. The correct truncation has it. It is
  uniform across the family, so it is inherited from upstream rather than a
  CLASS++ regression, and it only bites at `K != 0`. **This change does not touch
  it** — the edited line keeps its `F_{L-1}` coefficient exactly as it was, so the
  diff is the damping factor and nothing else. Separate issue. **[C]**
* **`l_max_dr > 17` on a dncdm composite** was the #421 heap overflow; the
  `l_max_dr = 40` run in §4.2 is exactly that case and it ran clean, which is an
  incidental confirmation of PR #424. **[M]**

## 5.5 · How the change was verified

Both binaries built in **one** build directory, so they differ by the change and
not by the build: **[M]**

| check | result |
|---|---|
| `make test` | 43/43 pass, including `test-limber-closure` with the new massive-daughter assertions |
| LCDM (`tCl,pCl,lCl,mPk` + lensing), master vs branch | `cl_lensed.dat` and `pk.dat` **byte-identical** |
| `dcdm_dr`, master vs branch | **byte-identical** — the helper move out of `DarkRadiationSpecies` into an inline header is a pure refactor and does not move a digit, including under `-ffast-math` |
| `dcdm_wdm`, master vs branch | TT 1.02e-05 (at `l = 2500`), P(k) 1.02e-02 (at `k = 2.70 h/Mpc`) — the intended change, and nothing else moved |
| ASan (`-fsanitize=address`), `dcdm_wdm` at `l_max_ncdm = 10` and `100` | clean, real output produced |

The `dcdm_dr` row is the load-bearing one: it is the species whose code this PR
touches without intending to change its numbers.

## 6 · Rejected alternatives

* **A precision knob to select the closure**, for `dcdm_wdm`. Same answer as
  #420: MB is not a valid option for this species, and a knob would only preserve
  the ability to reproduce a wrong number.
* **Applying the change to `DrPsdSpecies` for symmetry.** §4. The two species
  have the same code shape and opposite physics, which is the whole content of
  this note.
* **Caching the Limber order on `NCDMBaseSpecies::PerturbLayout`.** It would put
  a field on every NCDM species that only one of them uses.
  `WdmDecayProductSpecies::PerturbLayout` derives from it instead, mirroring
  `DarkRadiationSpecies`. It cannot live on the species: `RegisterPerturbationIndices`
  runs once per `perturb_vector`, i.e. concurrently across k-modes. **[C]**
* **Leaving the helpers as `DarkRadiationSpecies` statics** and calling them from
  `wdm_decay_product.cpp`. #422 asked for the move; a second species calling into
  a first species' statics is a dependency neither wants.

## 7 · Limits of what is claimed

* The `dcdm_wdm` evidence is one cosmology varied along three axes — `vkick`,
  `Gamma` and `momenta_bins`, plus an evolver control. It was **not** varied in
  `Omega_ini` or in the background cosmology. **[?]**
* `L` in the curved formula is the **flat** moment ratio, and no curved cosmology
  was run — unchanged from #420, and now shared by one more species. **[?]**
* The `dr_psd` measurement is at `Gamma = 1e5`, `m = 0.3 eV`, on both branches.
  The campaign runs `Gamma = 1e9..1e12`, where the decay is earlier still — which
  pushes further towards coherence, i.e. towards this note's conclusion, but was
  not measured. Its `l_max_ncdm = 17` truncation error was measured only at
  `Gamma = 1e5`, and the daughter grid there (`dr_N_q = 60/100`) is far coarser
  than the campaign's 219. **[?]**
