# One effective multipole for the Limber delta

**Status (2026-09-08).** Shipped: `TransferModule::LimberDeltaPoint(l, q, sgnK, K)`
returning the evaluation distance and the curvature factor at a single
`L = l + 1/2`, `transfer_limber` reduced to one call to it, a
`LimberSourceCrossesEquator` warning for closed models past the equator, and
`LimberGridQMin`, which makes the closed full-Limber grid converge to the flat
one as `K -> 0`. Fixes #423.
Flat cosmologies are **byte-identical** (`cl.dat` and `cl_lensed.dat`, `diff` clean)
**[M]**; curved output moves by `<= 3.7e-4` in `C_l^phiphi` for `|Omega_k| <= 0.1`
**[M]**. Unit test: `source/limber_point_test.cpp` (`test-limber-point`).

## The bug

`transfer_limber` used a different effective multipole in *four* places within
nine lines:

| | expression | effective `L` |
|---|---|---|
| flat turning point | `(l+0.5)/q` | `l + 1/2` |
| closed turning point | `asin(sqrt(l*(l+1))/q*sqrt(K))` | `l + 1/2 - 1/(8l)` |
| open turning point | `asinh((l+0.5)/q*sqrt(-K))` | `l + 1/2` |
| curved amplitude | `pow(1. - K*l*l/q/q, -0.25)` | `l` |

Issue #423 reports the first three. The fourth is the same quantity again, and it
is the *larger* term. Substituting `U = Phi * sin_K(chi)` turns the hyperspherical
radial equation into the flat one with `chi -> sin_K(chi)`,

    U'' + [q^2 - L^2 / sin_K^2(chi)] U = 0,

so the delta sits at `sin_K(chi_L) = L/q` and the WKB weight of a delta at a
turning point carries `cos_K(chi_L)^(-1/2) = (1 - K L^2/q^2)^(-1/4)` -- the *same*
`L` in both. Using `l` there is an `O(1/l)` relative error, against `O(1/l^2)` for
the turning point:

    delta(ln A) = s / (4 L (1-s)),    s = K L^2/q^2 = sin^2_K(chi_L)

and it applies to **open** geometries too, which #423 reports as clean. **[C]**

## Which convention, and why not on accuracy

In units of `1/(8l)` about `l + 1/2`, the candidates are `-1` (classical turning
point `sqrt(l(l+1))`), `0` (Langer, `l + 1/2`), `+1` (`M_1/M_0`, the point at which
the first-order term of the source expansion drops out of the convolution) and `-3`
(`M_0/M_{-1}`, the #420 closure -- a different quantity, see its spec). **[M]**

None of them is optimal. Inverting `S(x_*) M_0/q = Int S j_l dchi` numerically for
the point that would make Limber *exact* gives, in the same units:

```
  source          l=40     l=100    l=400
  chi^-2         -8.88     -9.19    -7.16
  exp(-c/3000)   -7.04     -7.11    -7.14
  exp(-c/1e4)    +1.65     +1.66    +1.89
```

The optimal point is a property of the *source*, not of the Bessel function, and
sits an order of magnitude from every candidate: the leading residual is the `S''`
term, which no relocation removes -- that is what the extended-Limber series
(`transfer_limber2`) is for. All three conventions differ by less than the
approximation's own error, so this is a consistency question. **[M]**

`l + 1/2` wins on three consistency grounds:

1. WKB for a radial equation needs the Langer replacement `l(l+1) -> (l+1/2)^2` to
   reproduce `j_l`; the classical turning point is the wrong one for Bessel
   asymptotics.
2. `IPhiFlat` is the asymptotic expansion of `Int_0^inf j_l(x) dx` -- verified:
   `sqrt(pi/2l)(1 - 1/4l + 1/32l^2)` is that Gamma ratio to `O(l^-2)` **[M]**. It is
   the weight of the prescription `j_l -> sqrt(pi/(2l+1)) delta(x - (l+1/2))`;
   location and weight belong to the same prescription.
3. Standardising the other way would move every *flat* lensing spectrum and diverge
   from `transfer_limber2`, from LoVerde & Afshordi 0809.5112 which the code cites,
   and from `class_public`, whose `transfer_limber` carries the identical code --
   this is inherited, not a CLASS++ regression. **[C]**

## Why fix a sub-per-mille effect: continuity at Omega_k = 0

Open -> flat was continuous; closed -> flat was not, and the offset survived the
limit. Measured on `C_l^phiphi` over `41 <= l <= 2000`, fixed binary vs reference
binary built in the **same** build directory (the flat runs are bit-identical, so
the physical `Omega_k` dependence divides out of this column):

| `Omega_k` | convention difference | at l |
|---|---|---|
| `-1e-3` | 2.216e-04 | 46 |
| `-1e-4` | 2.374e-04 | 46 |
| `-1e-5` | 2.388e-04 | 46 |
| `-1e-6` | 2.389e-04 | 46 |

It plateaus instead of vanishing: a 2.4e-4 step in `C_l^phiphi` near `l = 46`
across `Omega_k = 0`, which any numerical `dC_l/dOmega_k` straddling zero picks up
with a contamination that scales as `1/h`. `l = 46` is where Limber switches on
(`l_switch_limber = 40`). **[M]**

## What was measured

`h=0.6732, omega_b=0.02238, omega_cdm=0.1201, tau_reio=0.0543, n_s=0.966,
A_s=2.1e-9, N_ur=3.044`, `output = tCl,pCl,lCl`, `lensing = yes`,
`l_max_scalars = 2500`, default precision. Metric is `max |dC_l|` over the local
envelope (`max |C|` within `+-40` in l): `TE`, `TPhi` and `EPhi` cross zero, and a
plain max-rel-diff reports 3-4% there for an arbitrarily small absolute shift.

| `Omega_k` | `C_l^phiphi` | at l | lensed `TT` | at l |
|---|---|---|---|---|
| `+0.2` | 3.724e-04 | 151 | 3.243e-05 | 2057 |
| `+0.1` | 2.199e-04 | 149 | 2.279e-05 | 1845 |
| `+0.05` | 1.224e-04 | 148 | 1.698e-05 | 2403 |
| `+0.01` | 2.703e-05 | 156 | 3.935e-06 | 2281 |
| `0` | **0 (byte-identical)** | -- | **0** | -- |
| `-0.01` | 2.132e-05 | 222 | 2.142e-06 | 2225 |
| `-0.05` | 1.514e-04 | 161 | 2.315e-05 | 2094 |
| `-0.1` | 3.729e-04 | 161 | 7.087e-05 | 2478 |
| `-0.2` | 1.528e-03 | 136 | 3.674e-04 | 2500 |
| `-0.3` | 1.465e-01 | 1957 | 5.179e-03 | 2411 |
| `-0.5` | 2.550e-01 | 1779 | 1.657e-02 | 1977 |

The last two rows are not this change going wrong; they are the closed-space Limber
path leaving its domain of validity. See below. **[M]**

## The closed-space validity limit (not fixed here)

`asin` returns only the *near* turning point. A closed universe has a second one at
`pi R - chi_L`, and the source support reaches past the equator once
`chi_max > (pi/2) R`. Measured from the background table for this cosmology:

| `Omega_k` | `R` (Mpc) | equator `(pi/2)R` | `chi_rec` | `chi_max` | past equator |
|---|---|---|---|---|---|
| `-0.05` | 19915 | 31283 | 14064 | 14345 | no |
| `-0.1` | 14082 | 22120 | 14247 | 14528 | no |
| `-0.2` | 9958 | 15641 | 14648 | 14928 | no |
| `-0.3` | 8130 | 12771 | 15103 | 15384 | **yes** |

Past that threshold (`|Omega_k| ~ 0.25` here) two things break together: the
far-side stationary contribution is missing entirely, and the surviving modes are
all near-grazing (`q >= L sqrt(K)` cuts away the modes a flat-like integral would
weight, because `L/chi_rec < L sqrt(K)`), where `s -> 1` and the WKB amplitude
`(1-s)^(-1/4)` diverges. That is why the convention choice, worth 1e-4 in the
physical range, is worth 15-25% there: for the edge mode `nu = l+1` at `l = 1957`,
`(1 - l^2/nu^2)/(1 - L^2/nu^2) = 2.003`, so the two amplitudes differ by
`2.003^(1/4) = 1.19`. **[M]**

Neither version crashes: `Omega_k` from `+0.2` to `-0.5`, exit 0, no NaN or inf in
any output column **[M]**. It degrades silently, so `transfer_init` now warns, gated
on `transfer_verbose > 0` as the other module warnings are. Verified to fire at
`Omega_k = -0.25` and `-0.3` (`tau0 - tau_rec = 1.51e4` against an equator at
`1.277e4`) and to stay silent at `-0.2` and at `+0.3` **[M]**. The compared
distance is `tau0 - tau_rec`, not `tau0`: the CMB lensing source is tabulated from
recombination to today, so that is where its support ends.

**No cheap correct fix exists.** The complete treatment adds the mirror term,
`W [S(chi_L) + (-1)^(nu-l-1) S(pi R - chi_L)]` -- the two Airy regions of a single
eigenfunction of definite parity, `Phi_l^nu(pi - x) = (-1)^(nu-l-1) Phi_l^nu(x)`.
It is a dozen lines, but validating it needs an oracle, and exact integration at
these l is itself unconverged at default precision (measured: turning Limber off
at `Omega_k = -1e-6` *increased* the flat-limit disagreement to 6.9e-2 **[M]**).
Since the regime is excluded by roughly 25 sigma, a warning is the honest
deliverable and the mirror term is left unimplemented. **[C]**

Open geometries have no such limit at any `Omega_k`: `asinh` is monotone and
unbounded and `(1 + |K| L^2/q^2)^(-1/4) < 1` is bounded and smooth. **[C]**

## Two pre-existing precision defects found while checking the flat limit

A nearly-flat closed model does not reproduce the flat one. Both causes are
unrelated to the effective multipole and are **not fixed here**; both are worth
their own issue.

### 1. Curved runs under-sample the radial functions at low l

`hyper_sampling_flat = 8.0` carries the comment "should remain >7.5", but the
curved paths default to `hyper_sampling_curved_low_nu = 7.0` and
`hyper_sampling_curved_high_nu = 3.0`. For `nu < hyper_flat_approximation_nu`
(4000) the curved run evaluates the true hyperspherical table at that sampling,
and those are exactly the modes that dominate low l. Raising both to 8 at
`Omega_k = -1e-6`, against the flat run:

| quantity, `l <= 10` | sampling 3/7 (default) | sampling 8 | sampling 16 |
|---|---|---|---|
| `C_l^EE` | 2.31e-01 | 4.26e-03 | 4.52e-03 |
| `C_l^TT` | 1.57e-02 | 1.59e-03 | 1.55e-03 |
| `C_l^phiphi` | 1.01e-02 | 3.97e-04 | 3.99e-04 |

A **23% error in `C_l^EE` at low l in curved models**, saturating by 8, i.e. the
default is 54x worse than converged and the flat run is the accurate one. **[M]**

**Shipped: `low_nu` 7 -> 8.0, `high_nu` left at 3.0** (owner's call, 2026-09-08).
Separating the two knobs against a `sampling = 20` reference shows they own
different regimes, and only one of them is free:

| | `Omega_k = -0.01`, EE `l<=10` | time | `Omega_k = -1e-6`, EE `l<=10` | time |
|---|---|---|---|---|
| default (3 / 7) | 5.73e-04 | 0.44 s | 2.32e-01 | 0.27 s |
| `high_nu = 8` | 5.73e-04 | 0.67 s | 7.19e-04 | 0.27 s |
| `low_nu = 8` | 2.29e-04 | 0.44 s | 2.32e-01 | 0.27 s |
| both = 8 | 2.29e-04 | 0.67 s | 2.66e-04 | 0.28 s |

`low_nu` carries physical curvature and is free there, so it is raised. `high_nu`
carries the near-flat regime, where it is also free, but costs **+52%** at
`|Omega_k| = 0.01` (best-of-5, closed; +12% open) while buying nothing measurable
there, so it is **kept at 3.0**: a linearly sampled MCMC almost never visits
`|Omega_k| < 1e-4`, and 52% on every closed run is too high a price for a regime
nobody samples. The reasoning is recorded at the default itself, including when to
raise it -- numerical `dC_l/dOmega_k` with a step below ~1e-4, where `high_nu = 8`
is worth 15x/51x on `C_l^phiphi`. **[C]**

What `low_nu` buys, at physical curvature, against the `sampling = 20` reference
(`l <= 10`):

| | TT | EE | `phiphi` |
|---|---|---|---|
| `Omega_k = -0.01` | 6.67e-05 -> 2.86e-05 | 5.73e-04 -> 2.29e-04 | 2.45e-05 -> 9.70e-06 |
| `Omega_k = +0.01` | 1.58e-04 -> 5.91e-05 | 1.50e-03 -> 5.35e-04 | 3.86e-05 -> 1.56e-05 |

2.5-2.8x for free. Flat is untouched -- it reads `hyper_sampling_flat` -- and
stays byte-identical. **[M]**

Flat-limit gap that remains at `l <= 10` with the shipped defaults, against the
flat run:

| `Omega_k` | TT | EE | `phiphi` |
|---|---|---|---|
| `-1e-4` | 1.56e-03 | 3.92e-03 | 6.51e-04 |
| `-1e-5` | 1.56e-03 | 4.23e-03 | 1.48e-03 |
| `-1e-6` | 1.57e-02 | **2.31e-01** | 1.01e-02 |
| `-1e-8` | 1.64e-03 | 2.35e-03 | 2.88e-04 |

The `-1e-6` row is the accepted cost of keeping `high_nu = 3`: the band is worst
when it lands exactly on the low-l power, which happens near `|Omega_k| ~ 1e-6`
and passes above and below it. Everything else is the `nu = 3` anchor of the main
q list and the perturbation `k_min`. **[M]**

### 2. `q_logstep_limber = 1.1` is too coarse, in every geometry

`C_l^phiphi` above `l_switch_limber` is computed on a separate grid built by
`transfer_get_q_limber_list` (`:904`), a pure geometric progression
`q_i = q_min * q_logstep_limber^i` anchored at a geometry-dependent `q_min`:
`k_min` for flat, `sqrt(k_min^2+K)` for open, and **`3 sqrt(K)` for closed**,
which tends to 0 rather than to `k_min` as `K -> 0`. The closed grid therefore
never converges to the flat one, and the two quadratures sample the same
integrand at different phases.

The give-away was that `C_l^phiphi` is **bit-identical** under changes to
`q_linstep` -- it does not use that grid at all -- while `TT` moves by 1e-5. The
flat-vs-curved difference is invariant (to four digits) under `q_linstep`,
`q_logstep_spline`, `k_step_sub/super`, `perturb_sampling_stepsize`, `l_linstep`
and `hyper_sampling_curved_*`, and varies non-monotonically with `Omega_k`
(1.46e-3, 8.4e-4, 1.5e-4, 5.7e-4 at 1e-6, 1e-8, 1e-10, 1e-12) -- the signature of
a quadrature phase, not a formula difference. `q_logstep_limber` is the one knob
that moves it. **[M]**

The absolute error this leaves in a **flat** run's `C_l^phiphi`, against
`q_logstep_limber = 1.01`:

| `q_logstep_limber` | `l > 1500` | `100 < l < 1000` |
|---|---|---|
| 1.1 (default) | 3.14e-03 | 3.26e-03 |
| 1.05 | 7.39e-04 | 7.94e-04 |
| 1.02 | 1.48e-04 | 1.00e-04 |

O(h^2) convergence, and `want_lcmb_full_limber` defaults to **true**, so every
`lCl` run carries a ~0.3% `C_l^phiphi` quadrature error at the default. **[M]**
That absolute error is accepted as a default (owner's call, 2026-09-08); what is
*not* acceptable is that the two geometries disagree about where the nodes go.

**Fixed: `LimberGridQMin`.** All three geometries already used one rule,
`q_min = sqrt(k_min^2 + K)`; the non-convergence was inherited from the closed
perturbation `k_min = sqrt((8-1e-4) K)` (`perturbations_module.cpp:1425`), the
`nu = 3` mode, which vanishes as `K -> 0` while flat truncates at
`k_min_tau0/tau0`. Open needs nothing: its `k_min = sqrt(-K + (k_min_tau0/tau0/ar)^2)`
collapses to the flat anchor exactly. The closed anchor is now clamped by the
flat one -- written as a pure clamp of the existing `3 sqrt(K)`, so it bites only
below the crossover `|Omega_k| ~ 1.3e-4`, and the nodes it drops carry no
integrand (the Limber transfer vanishes for `q < l_switch_limber/tau0`).

`C_l^phiphi`, `l > 1500`, closed vs flat:

| `Omega_k` | before | after |
|---|---|---|
| `-1e-4` | -- | 3.03e-04 |
| `-1e-5` | -- | 2.78e-05 |
| `-1e-6` | 1.459e-03 | 1.39e-05 |
| `-1e-8` | 8.386e-04 | 1.78e-05 |
| `-1e-10` | 1.494e-04 | 1.78e-05 |
| `-1e-12` | 5.682e-04 | 1.79e-05 |

Erratic at the 1e-3 level before, a smooth monotone approach to a stable 1.8e-5
floor after -- ~80x. **Byte-identical** at `Omega_k` = 0, +-0.01, +-0.05, +-0.1,
+-0.2, -0.3, -0.5: no cosmology anyone runs moves at all. **[M]**

What the residual 1.8e-5 is *not*: it is flat in `K`, so it is not this anchor.
`TT` and `EE` at `l > 1500` sit at 4e-05 and 7e-05, untouched by this change --
they ride the main `transfer_get_q_list`, whose closed branch still anchors at
`nu_min = 3`. That one cannot take the same clamp without also moving the integer
`nu` quantisation and the first-point correction in
`spectra_module.cpp:1302`, so it is left alone. **[C]**

At low l the flat limit is still governed by defect 1 above, not by any grid:
at `Omega_k = -1e-6`, `l <= 10` disagreement is 2.31e-01 (EE) and 1.57e-02 (TT)
before and after. Raising `hyper_sampling_curved_*` to 8 is the remaining lever
and has not been taken. **[M]**

## Not touched

- `sqrt(l(l+1))` at `transfer_module.cpp:2526-2531` and `:3028-3032`. There it is
  the turning point used to map `Phi_l^nu` onto a rescaled `j_l` in the flat
  approximation, applied self-consistently on both sides of that mapping; it is a
  sampling/rescaling choice, not an evaluation point. **[C]**
- The `L(l) = M_0/M_{-1}` of the #420 dark-radiation closure. Different quantity;
  see `2026-09-07-dr-limber-closure-design.md`, which says so explicitly.
- The flat branch's arithmetic. `*trsf /= (l + 0.5)` is kept rather than folded
  into the common `amplitude / (tau0_minus_tau_limber * q)` so that the flat path
  cannot pick up a 1-ULP round trip. That is what makes the byte-identical result
  above possible.

## Guarantees the fix relies on

`nu = q/sqrt(K)` is an integer `>= l+1` in a closed universe, so `L/nu < 1` and
both `asin` and the `(1 - K L^2/q^2)^(-1/4)` factor are always finite. Enforced at
`transfer_module.cpp:3473` (the l list is trimmed to `l < nu`) and `:1480`
(transfers zeroed for `l >= nu`). `test_closed_edge_mode_is_representable` pins it
at `nu = l+1`. **[C]**
