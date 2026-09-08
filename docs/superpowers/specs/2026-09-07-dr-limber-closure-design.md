# The decay-radiation hierarchy is closed at the wrong Bessel argument

**Status (2026-09-07).** Shipped: the Limber closure described below, in
`DarkRadiationSpecies::PerturbDerivs`, replacing the Ma & Bertschinger truncation
outright (no knob). It covers `dcdm_dr` *and* the default `dncdm_dr`
(`dr_representation = integrated`), whose daughter is a `DarkRadiationSpecies`
instance (`species/dncdm_dr_species.cpp`). Not shipped, and not attempted here:
the same treatment for `DrPsdSpecies` (`dr_representation = psd`) and for the
`dcdm_wdm` massive daughter — see *Where else this applies*. Every number below
tagged **[M]** was measured for this note on 2026-09-07; **[C]** is read off the
code; **[?]** has no evidence behind it.

Resolves issue #418. The issue proposed a different (algebraically equivalent in
flat space) form; *Rejected alternatives* says why this one shipped instead.

## The problem

`species/dark_radiation_species.cpp` closed its `l = l_max` branch with

```cpp
dy[base + l] = k * (s_l[l] * y[base + l - 1] - (1. + l) * cotKgen * y[base + l]);
```

Invert that against the exact hierarchy `dF_l/dtau = k/(2l+1) [l s_l F_{l-1} -
(l+1) s_{l+1} F_{l+1}]` and it is precisely the hyperspherical Bessel recurrence

    s_{l+1} F_{l+1} + s_l F_{l-1} = (2l+1) k cotKgen F_l                    [C]

with `cotKgen = cot_K(chi)/k` evaluated at `chi = tau`
(`source/perturbations_module.cpp:5647`). That is the correct recurrence — for a
**single coherent free-streaming mode released at tau = 0**. So the closure's
entire content is the *argument* it evaluates `cot_K` at, and that is the only
thing wrong with it.

Decay radiation is not a single coherent mode. It is injected continuously at
l = 0,1, so it is a superposition of shells, each free-streaming after its own
release time:

    F_l(tau) = Int_0^tau S(tau') j_l(k(tau - tau')) dtau'

The recurrence holds pointwise in the argument, so it does not survive a
superposition over different arguments `k(tau - tau')`.

The failure is not subtle. At `k tau = 325`, `l = 17`, MB's prefactor
`(2l+1)/(k tau) = 0.108` is negligible, so MB collapses to `F_{l+1} ~ -F_{l-1}`.
That is right for the sign-alternating ladder a genuinely free-streaming species
has, and guaranteed wrong for the monotone, single-signed ladder decay radiation
has. `ur` never hits this because RSA/UFA retires its hierarchy while `k tau <~
2 l_max + 1` still holds; DR has neither approximation and runs to
`k tau = 14208` at `l_max_dr = 17` — three orders of magnitude outside the
closure's validity. **[C]**

## The fix

Evaluate `cot_K` at the distance where the Bessel weight of the superposition
actually sits, rather than at the current distance. The code shape is unchanged;
one factor is substituted.

**The order.** Apply the recurrence *inside* the integral. With `y = k(tau-tau')`,
`w(y) = S(tau - y/k)`:

    F_{l+1} + F_{l-1} = (2l+1) Int w(y) j_l(y) / y dy

so the argument the ladder behaves at is, exactly and with no asymptotics,

    L = (2l+1) F_l / (F_{l+1} + F_{l-1}) = Int w j_l dy / Int w j_l y^-1 dy

`F_l` itself is an incomplete integral, but the closure needs only this ratio,
and for a weight slowly varying across the Bessel oscillation the full-range
Weber–Schafheitlin moments `Int_0^inf y^mu j_l(y) dy = 2^(mu-1) sqrt(pi)
Gamma((l+mu+1)/2) / Gamma((l-mu+2)/2)` give `w` dropping out:

    L(l) = 2(l+1)/l [Gamma((l+1)/2)/Gamma(l/2)]^2
         = l + 1/2 - 3/(8l) + 3/(16 l^2) + O(l^-3)

**This is the Limber approximation.** `j_l(x) ~ sqrt(pi/(2l+1)) delta(l+1/2 - x)`
substituted into that ratio gives `L = l + 1/2` identically, and the closure
coefficient `(2l+1)/L = 2` exactly — plain linear extrapolation of the ladder.
The Gamma-ratio is the *extended*-Limber correction to it. The transfer module
already carries half of this derivation: `IPhiFlat` in `transfer_limber`
(`source/transfer_module.cpp:2675`) is the `mu = 0` moment `Int j_l dy`, the
numerator of `L`. **[C]**

**Curvature.** Taking the turning points from `TransferModule::transfer_limber`
(`:2648-2657`), with `q^2 = k^2 + K` for scalars, all three geometries collapse
to one expression:

| | turning point | `cot_K(chi_L)/k` |
|---|---|---|
| closed | `sin(sqrt(K) chi_L) = sqrt(K) L / q` | `sqrt(q^2 - K L^2)/(k L)` |
| open | `sinh(sqrt(-K) chi_L) = sqrt(-K) L / q` | `sqrt(q^2 - K L^2)/(k L)` |
| flat | `chi_L = L/q` | `1/L` |

    cotK_eff(l) = sqrt(1 - K (L^2 - 1)/k^2) / L

which is the `s_l` expression evaluated at non-integer order `L`, over `L`. It is
exact in flat space, reduces to `1/L` there, and clamps to zero above
`nu = q/sqrt(K)` in a closed universe, where the turning point does not exist and
`s_l` clamps too. `species/dr_closure_test.cpp` checks the closed form against a
literal `asin`/`asinh` evaluation of the turning point in both curved geometries,
to 1e-12. **[M]**

## Three different "effective l" live in this tree — do not harmonise them

| quantity | where | value |
|---|---|---|
| `sqrt(l(l+1))` | `transfer_limber`, closed only | `l + 1/2 - 1/(8l)` |
| `l + 1/2` | `transfer_limber` flat/open, `transfer_limber2` | `l + 1/2` |
| `L(l)` | this closure | `l + 1/2 - 3/(8l)` |

The first two are the *same* quantity (a Bessel turning point) written two ways,
and their disagreement at O(1/l) is a real internal inconsistency in
`transfer_limber` — pre-existing, affects lensing output, **not touched here.**
The third is a *different* quantity: a moment ratio `M_0/M_{-1}`, not a peak
location. Its `-3/(8l)` is not in conflict with the other `-1/(8l)`, and
"fixing" one to match the other would be wrong. **[C]**

## What was measured

Config: `Gamma_dcdm = 100`, `Omega_ini_dcdm = 0.05`, `omega_cdm = 0.09`,
`h = 0.6732`, `N_ur = 2.0308`, one 0.06 eV ncdm, `tCl,pCl,lCl,mPk` + lensing,
`l_max_scalars = 2500`, ndf15, `tol_perturb_integration = 1e-6`. Both closures
built in one build directory; metric is max relative difference over the whole
spectrum (both observables are positive definite, so no zero-crossing artefact).
All runs ASan-clean. **[M]**

**Truncation independence** — each closure against its *own* `l_max_dr = 200`:

| `l_max_dr` | MB TT | MB P(k) | Limber TT | Limber P(k) |
|---|---|---|---|---|
| **17 (default)** | 2.09e-04 | 1.48e-03 | **1.67e-06** | **2.96e-06** |
| 25 | 6.83e-05 | 1.01e-03 | 8.15e-07 | 1.07e-06 |
| 80 | 2.03e-05 | 7.02e-04 | 2.16e-09 | 1.21e-07 |

The two closures disagree at the default by TT 2.08e-04 / P(k) 1.47e-03. That is
the systematic removed. MB converges non-monotonically in P(k) (80 is barely
better than 25); the Limber closure falls monotonically. **[M]**

**This is validation, not self-consistency.** The open caveat in #418 was that
both closures still differ at `l_max_dr = 200`, because the ladder never reaches
its natural edge. Pushing MB up instead settles it — MB converges *to the Limber
answer*: **[M]**

| MB `l_max_dr` | vs Limber-200, TT | vs Limber-200, P(k) |
|---|---|---|
| 17 | 2.09e-04 | 1.48e-03 |
| 200 | 3.65e-06 | 4.22e-04 |
| 400 | 7.94e-07 | 1.89e-04 |
| 800 | 1.26e-07 | 7.47e-05 |
| 1600 | 1.27e-08 | 2.28e-05 |

Richardson-extrapolating the MB sequence off the 400/800/1600 doubling triple
(measured order p = 2.88 for TT, 4.31 but unstable for P(k)) gives a limit that
agrees with the Limber closure to **5.0e-09 in TT**. For P(k) the agreement is
2.0e-5 at the measured p, and ranges 6.8e-6 to 3.2e-5 as p is varied over 1..8 —
i.e. consistent with zero within the extrapolation's own uncertainty. So the two
closures describe the same physics; the Limber one reaches it at `l_max_dr = 17`
where MB needs `l_max_dr` of order 1600. **[M]**

**Robustness in the decay rate.** The concern that a fast decay makes injection
effectively instantaneous — restoring MB's single-coherent-mode premise — was
tested and is false over four decades. `l_max_dr = 17` against `200`, per arm: **[M]**

| `Gamma_dcdm` | MB TT | MB P(k) | Limber TT | Limber P(k) | P(k) gain |
|---|---|---|---|---|---|
| 10 | 3.40e-05 | 1.94e-04 | 1.29e-06 | 4.52e-07 | 430x |
| 100 | 2.09e-04 | 1.48e-03 | 1.67e-06 | 2.96e-06 | 500x |
| 1000 | 1.86e-04 | 3.41e-03 | 1.16e-06 | 6.49e-06 | 525x |
| 10000 | 1.63e-04 | 3.62e-03 | 1.69e-06 | 6.66e-06 | 543x |
| 100000 | 1.13e-04 | 3.79e-03 | 2.32e-06 | 6.99e-06 | 542x |

**Cost** is unaffected: `L` depends only on `l_max`, so it is computed once at
`RegisterPerturbationIndices` and cached on the layout; the RHS gains one
compare and one divide. **[M]**

## Rejected alternatives

* **The expanded form proposed in #418** — reconstruct `F_{l_max+1}` explicitly and
  feed it to the exact recursion. Algebraically identical in flat space, but it
  reintroduces `s_l[l_max+1]`, which is off the end of an array sized
  `max_l_max + 1`, and it discards MB's compact form which has *already* absorbed
  `s_{l+1} F_{l+1}` correctly for curved space. Keeping the compact form makes the
  curvature generalisation one factor instead of an array resize. **[C]**
* **`L = l`** (linear extrapolation with `c = 2 + 1/l`). An earlier env-gated
  experiment measured 57x on TT and 48x on P(k), against 125x/498x for the exact
  moment ratio. The correct `L` is free — it is a compile-time-constant call per
  run — so there is no reason to take the cheaper one. **[P]**
* **A precision knob to select the closure.** MB is not a valid option for this
  species; a knob would only preserve the ability to reproduce a wrong number.

## Limits of what is claimed

* Steps 1–4 are the `k tau -> inf` limit. `x_eff(17)` measures 16.91 against the
  asymptotic 17.48 at `k tau ~ 2l` (3% low), so the asymptote wants
  `k tau >~ 10 l`. A finite-`x` correction is where the injection history `S`
  re-enters. **[P]**
* `L` in the curved formula is the **flat** moment ratio. The genuinely curved
  object would be a hyperspherical moment ratio; taking the flat one is right to
  O(1/l), and the curvature factor itself is only appreciable when
  `K(L^2-1)/k^2` is not small. No curved cosmology was run. **[?]**
* The evidence is `dcdm_dr` only. `dncdm_dr` in the `integrated` representation
  shares the code and is therefore fixed, but was **not** independently measured
  (the attempt is described below). **[?]**

## Where else this applies

`DrPsdSpecies` (`dr_representation = psd`, `species/dr_psd_species.cpp:534,560`)
and the `dcdm_wdm` massive daughter (`species/wdm_decay_product.cpp:586`) have
the same continuously-sourced structure and the same MB truncation. Neither is
changed here and neither has been measured. **[?]**

## A separate bug found on the way

Auditing the dncdm path turned up a **pre-existing heap-buffer-overflow** that is
not related to this closure and reproduces on master with the MB code.
`ppw->max_l_max` is built by name-picking species
(`source/perturbations_module.cpp:1861`): it pulls in `ppr->l_max_dr` only when a
species keyed `"DCDM_DR"` is present. A DNCDM composite also owns a
`DarkRadiationSpecies` daughter running to `l_max_dr`, but is keyed `dncdm1`, so
`s_l` is allocated `max_l_max + 1 = 18` entries while that daughter indexes
`s_l[l+1]` for `l < l_max_dr`, i.e. up to `s_l[l_max_dr]`. Any dncdm run with
`l_max_dr > 17` reads
heap garbage: silently wrong spectra, and nondeterministic "singular matrix" /
"step size too small" aborts at roughly a 50% rate on identical input, which
survive single-threading. ASan names it at
`species/dark_radiation_species.cpp:107` against the allocation at
`perturbations_module.cpp:1883`. **[M]**

This is exactly the failure mode the "no species-type picking in modules" rule
exists to prevent, and the tensor branch immediately below carries a comment
about a previous instance of the same sizing bug. The durable fix is for
`max_l_max` to come from the species themselves — a `MaxMultipole(const
precision*)` virtual that composites forward to their children — rather than from
a hard-coded key list. Filed separately; **not fixed in this change.**
