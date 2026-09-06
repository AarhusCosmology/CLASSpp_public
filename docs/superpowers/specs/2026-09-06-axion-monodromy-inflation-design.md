# Axion monodromy in the inflation module

**Status, 2026-09-06.** All of section 3 (the `monodromy` potential shape), section 4
(the input plumbing) and section 5 (the two `primordial_inflation_one_k` fixes) shipped
together. **Section 5.2 is weaker than first written**: measured afterwards, the step
estimator of 5.1 is the fix that is load-bearing, and the time re-anchoring of 5.2 is not
required by any case in this study. It is kept because the guard it addresses is wrong on
its own terms, not because anything here needs it; 5.2 says so. Section 8 records what
was measured when asking whether these two call sites should move to the `evolver_*`
interface, and why they did not.

**Review, 2026-09-06 (PR #415).** Two of three automated review findings were acted on:
the monodromy parameter-domain checks now in section 3, and a symbol collision in
`explanatory.ini` (it used `x` for the field variable and then again for the resonance
depth, which is a conformal time $-k\tau$). The third — a claimed missing `return` in
`read_potential_shape` — was **not** a defect: `ThrowFormattedSevere` is `[[noreturn]]`,
and compiling that translation unit with `-Wall -Wextra -Werror=return-type` emits no
return-type diagnostic **[M]**. A dead `return` there would read as a fallback default and
would suppress a real warning if the attribute were ever dropped. Section 6 lists the shortcomings that were *found and not fixed*; none of them
is a blocker for the physics and each says so explicitly. The measured numbers in
section 7 are from `notebooks/axion_monodromy_flauger2009.py`, which regenerates all of
them, and the regression tests are
`TestReviewRegressions.test_axion_monodromy_reproduces_its_analytic_delta_ns` and
`...test_unknown_inflation_potential_is_rejected` in `python/test_class.py`.

Evidence tags: **[M]** measured in this work, **[C]** stated by the code, **[?]** no
evidence, do not act on without checking.

---

## 1. Why

The inflation module (`P_k_ini type = inflation_V` / `inflation_H` / `inflation_V_end`)
has been in CLASS++ since the port and has essentially never been used. It solves the
inflaton background and the Mukhanov–Sasaki equation mode by mode, which is exactly the
machinery a *feature* model needs and exactly what an analytic $P(k)$ parametrisation
cannot give you. Exercising it on a model with a real feature was the point; axion
monodromy inflation (Flauger, McAllister, Pajer, Westphal & Xu,
[arXiv:0907.2916](https://arxiv.org/abs/0907.2916)) is a good choice because the paper
derives a closed-form answer *and* shows where its own numerics stopped working, so
there is something to check against and something to beat.

## 2. The model, and the two conventions it collides with

$$V(\phi) = \mu^3\phi + \Lambda^4\cos(\phi/f)
          = \mu^3\left[\phi + b f\cos(\phi/f)\right],\qquad b \equiv \Lambda^4/(\mu^3 f).$$

Two conventions had to be reconciled, and both are the kind of thing that produces a
plausible wrong answer rather than an error.

* **Units.** The paper uses reduced Planck units; the module uses the *non-reduced*
  Planck mass — `primordial_inflation_derivs_member` writes $aH=\sqrt{(8\pi/3)(\ldots)}$
  and `primordial_inflation_get_epsilon` writes $\epsilon=(1/16\pi)(V'/V)^2$ **[C]**.
  Hence $\phi_{\rm CLASS}=\phi_{\rm paper}/\sqrt{8\pi}$ and
  $V_{\rm CLASS}=V_{\rm paper}/(8\pi)^2$. The potential is form-invariant under that
  rescaling — only $\mu^3$ and $f$ move — which is why the shape below is written once
  and the notebook does the conversion at the boundary.
* **Direction of roll.** `primordial_inflation_check_potential` requires
  $\mathrm{d}V/\mathrm{d}\phi<0$ with $\phi$ increasing **[C]**; the paper's inflaton
  rolls *down* towards zero. The shape therefore works in $x = V_4-\phi$.

## 3. The `monodromy` shape

```
V(phi)  = V0 * ( x^p + b f cos(x/f) ),        x = V4 - phi
```

with `V0` = $\mu^{4-p}$, `V1` = $f$, `V2` = $b$, `V3` = $p$, `V4` = $x_*$. For $p=1$ this
is the paper's eq. (2.2) exactly.

*Why $V_4$ is a parameter at all.* `inflation_V` pins $\phi_{\rm pivot}=0$
(`primordial_inflation_solve_inflation`) **[C]**, so the model needs to be told where the
pivot sits on the monodromy branch. `V4` *is* $\phi_*$, and it also fixes the phase of the
modulation at the pivot, so no separate phase parameter is needed.

*Why the $V_i$ slots and not named keys.* `natural` and `higgs_inflation` already
overload `V_0`…`V_4` with shape-specific meanings and document the correspondence in
`primordial_inflation_potential`. Introducing a second convention for one new shape would
have been the larger change. Rejected: named inputs (`monodromy_f`, …), which read better
in isolation and worse next to the two shapes beside them.

*Why $p$ is additive, not multiplicative.* $V_0 x^p + \Lambda^4\cos(x/f)$ keeps
$b=\Lambda^4/(\mu^{4-p}f)$ unambiguous. The multiplicative form
$\mu^{4-p}x^p[1+\delta\cos(x/f)]$ used in the later literature makes the modulation
amplitude drift with $x$ in a way that is a modelling choice, not a generalisation.
Not implemented: the *drifting* templates of
[arXiv:1412.1814](https://arxiv.org/abs/1412.1814), where the frequency itself runs; that
is a second shape, not a parameter of this one.

*Domain.* Two layers, both `class_test` on the numeric channel — a sampler should reject
the point, not abort.

* **State**, in `primordial_inflation_potential`: `x <= 0` throws. Without it `pow(x, p)`
  returns NaN for non-integer $p$.
* **Parameters**, in `check_monodromy_parameters` at input time (added in review of
  PR #415): $f>0$, $p>0$, $V_0>0$, $x_*>0$. These do not depend on $\phi$, so they belong
  where they are read rather than on every call.

Neither layer can be left to `primordial_inflation_check_potential` downstream. Its `V<=0`
and `dV>=0` tests are *ordered comparisons*, which IEEE-754 says are false against a NaN
and `-ffast-math` — which this project builds with — says nothing about at all. So
`f = 0`, which reaches it as a NaN through `cos(x/f)`, is exactly the input those tests
cannot be relied on to catch. Note that the stock defaults (`V_1 = -1.12e-14`, the
polynomial shape's linear Taylor coefficient) are *not* a valid monodromy, so
`potential = monodromy` with no `V_i` set is a caught error rather than a silent one.

## 4. The input plumbing, which was broken

`input_module.cpp` read the `potential` key and **threw the value away**:

```cpp
pfc->get<std::string>("potential");
/* only polynomial coded so far: no need to interpret the value **/
```

So `potential = natural` silently ran the polynomial shape with `V_1` reinterpreted as the
linear Taylor coefficient — a wrong answer, no warning **[M]**: this is how the first
`monodromy` run failed, with `dV/dphi = 0.00598` where 0.00598 was the $f$ that had been
passed as `V_1`. `natural` was in the `potential_shape` enum but reachable from no input
key at all, since `inflation_V_end`'s `full_potential` accepted only `polynomial` and
`higgs_inflation`.

Both keys now go through one `read_potential_shape()` helper that accepts all four shapes
and rejects anything else with `class_stop_severe` (a structural error — the name of a
potential is not a value a sampler varies). The `PSR_i`/`R_i` reparametrisations, which
are the slow-roll parameters *of the Taylor expansion*, are now refused for non-polynomial
shapes instead of silently overwriting the `V_i` the user set.

## 5. Two integration bugs in `primordial_inflation_one_k`

Both are pre-existing and both are latent for a smooth potential. A modulated one turns
5.1 into a hard failure, which is how they were found: `f = 0.01\,M_p`, `b = 0.08` died
with `generic_integrator: Step size too small: step:2.04e-16, minimum:2.22e-16, in
interval: [7.04071e+06:7.04389e+06]` **[M]**. 5.2 was found while diagnosing 5.1 and is a
correctness fix, not a fix this study needs — see the measurement at the end of 5.2.

### 5.1 The step-size estimator ignored the expansion

```cpp
dtau = pt_stepsize * 2*pi / max( sqrt(|k^2 - z''/z|), k );      // before
dtau = pt_stepsize * 2*pi / max({ sqrt(|k^2 - z''/z|), k, aH }); // after
```

The estimator uses the mode's *instantaneous* effective frequency, which passes through
zero whenever $z''/z$ crosses $k^2$. When it does, `max` falls back to $k$ and the step
jumps — measured 15× in one step, from 202 to 3176 **[M]** — while $a$ is growing
exponentially, so the adaptive integrator is handed several e-folds in one call and gives
up. For a smooth potential the crossing happens once, near horizon exit, and the jump is
survivable. For a modulated one $z''/z$ contains $a^2 V''\sim\pm 2(aH)^2$ and crosses
repeatedly.

`aH` never binds while the mode is sub-horizon (there $aH \le k/\texttt{ratio\_min}$), so
this is a no-op before horizon crossing. It is the same criterion
`primordial_inflation_evolve_background` already applies to the background
(`bg_stepsize * a/a'`).

Rejected: bounding $\mathrm{d}(z''/z)/\mathrm{d}\tau$ directly, which is the quantity that
actually matters but is not available without differentiating the RHS; and adding the
background oscillation frequency $\omega=\dot\phi/f$ explicitly, which is
potential-specific and would have to be plumbed out of the shape.

### 5.2 The minimum-step guard was measured against an arbitrary time origin

`generic_integrator` guards with `class_test(fabs(hnext/x1) <= hmin)` where `x1` is the
**absolute** start of the interval and `hmin` is machine epsilon **[C]**. `one_k` starts
conformal time at 0 and lets it accumulate; after ~1500 steps $\tau\sim10^7$, so the
smallest representable substep is $\sim10^{-9}$ *in absolute time* and a legitimate one is
rejected. The guard therefore tightens without bound the longer a mode is integrated, and
has nothing to do with the physics.

The RHS does not depend explicitly on $\tau$ — the code says so itself — so each interval
is now re-anchored at `[dtau, 2*dtau]`. The guard then means what it says: reject when a
substep falls to $\sim10^{-16}$ *of the step we asked for*. `generic_integrator` was left
alone; it is shared with background, thermodynamics and perturbations, and `rkqs`'s own
`xnew == *x` underflow test plus `_MAXSTP_` are the scale-free protections that actually
work.

Rejected: `x1 = 0`, which makes the guard `inf` and would never fire — and materialising
an infinity under `-ffast-math` is a trap this repo has been bitten by before.

**How much this one is worth: nothing measurable, so far. [M]** With 5.1 in place, letting
`tau` accumulate as before is fine at every point in this study, including the two most
expensive — `f = 2\times10^{-4}` at `ratio_min = 300` (32.7 s anchored / 33.8 s
accumulating) and the same at `ratio_min = 3000`, ten times the requirement and the
longest `tau` reached anywhere here (270 s / 308 s). Conversely 5.2 without 5.1 still dies
with the same "step size too small". So the ordering is unambiguous: **5.1 is the fix,
5.2 is insurance.** It is kept because a guard whose threshold depends on how long the
integration has been running is wrong whether or not something currently trips it, and
re-anchoring costs one line and no time.

### 5.3 A one-line correctness fix that came with them

`dlnPdN` divided the change in curvature by `dtau` *after* `dtau` had been overwritten
with the next step. Harmless while the step is smooth, wrong by 15× at exactly the moment
described in 5.1. It now divides by the step actually taken.

## 6. Found, not fixed

| # | Finding | Evidence | Why not fixed here |
|---|---|---|---|
| 1 | The inflation module refuses to run without tensor modes *and* a tensor $C_\ell$ source, so `output = mPk` alone is impossible. It is structural — `primordial_inflation_spectra` writes `is_non_zero_[index_md_tensors_]` unconditionally. | **[M]** two successive `class_test_severe`s before any physics runs | Guarding three writes is easy but is a separate change; the workaround (`output = tCl`, `modes = s,t`, stop at `compute(level=['primordial'])`) costs nothing. |
| 2 | If `k_pivot` falls outside the perturbation module's $[k_{\min},k_{\max}]$ — e.g. `l_max_scalars = 200` with the default `k_pivot = 0.05` — the run dies deep in `primordial_spectrum_at_k` with `k out of range`, not at input time. | **[M]** reproduced on the stock polynomial potential | Pre-existing, unrelated to this work, and the message does name the range. |
| 3 | `inflation_V` reads `V_0`…`V_4` while `inflation_V_end` reads `Vparam0`…`Vparam4` for the same struct members. | **[C]** | Renaming either is a user-visible input break. |
| 4 | `potential = natural` is now reachable but cannot work under `inflation_V`, because $\phi_{\rm pivot}=0$ and $V'(0)=0$ fails the negative-slope check. It needs `inflation_V_end`. | **[M]** | Correct behaviour: it fails loudly instead of silently running a different potential, which is what it did before. |
| 5 | `primordial_inflation_tol_integration = 1e-3` has *no effect at all* on a smooth potential (identical spectra from 1e-3 to 1e-8) but is the accuracy limit for a modulated one. | **[M]** | It is a precision parameter and the notebook documents the value this model needs. Changing the default would slow every user for a case almost none of them run. |

## 8. Should these two call sites move to the `evolver_*` interface?

Asked because the inflation module is now the **only** direct caller of
`generic_integrator` left: outside `tools/evolver_rkck.cpp` it appears in exactly two
functions, `primordial_inflation_one_k` (132 lines) and
`primordial_inflation_evolve_background` (272 lines) **[C]**. Every other module —
background, thermodynamics, perturbations — picks an `evolver_*` by function pointer from
`ppr->evolver_*` and calls it. The mentions in `perturbations_module.cpp` are stale
doc-comments.

The answer measured here is **not yet**, for three reasons.

### 8.1 The step-size bug is not the integrator's to fix **[M]**

`one_k` was patched to call `evolver_rkdp45` / `evolver_tsit5` per segment instead of
`generic_integrator`, and the `aH` cap of 5.1 was removed. At $f=0.01$, $b=0.08$:

| integrator | without the `aH` cap |
|---|---|
| `generic_integrator` (rkck) | `Step size too small: 2.04e-16` |
| `evolver_rkdp45` | `the monodromy potential is only defined for phi < V_4` |
| `evolver_tsit5` | same |

The modern evolvers fail *differently*, not less. Handed an outer step spanning many
e-folds they march the background off the end of the potential, and no embedded error
estimator rejects that: walking outside the model's domain is not a local truncation
error. The cap belongs in the module whichever integrator sits underneath it.

### 8.2 A per-segment swap is 2.1x slower **[M]**

Same experiment, both fixes in place, wall clock for the whole primordial module:

| | rkck | rkdp45 | tsit5 |
|---|---|---|---|
| $f=2\times10^{-3}$, 800 k/decade | 0.60 s | 1.09 s | 1.14 s |
| $f=5\times10^{-4}$, 1600 k/decade | 2.70 s | 5.65 s | 5.73 s |

The evolvers restart their step controller on entry (`absh` from
`|x_end - x_ini|/10`, or the first sampling gap). Calling one ~1600 times per mode throws
away precisely the step history that makes it good. A per-segment swap buys nothing and
costs a factor two.

### 8.3 What is actually missing is event termination

Both call sites are *integrate until a condition on the state*, and the shared evolver
signature takes a known `x_end` plus a pre-built `x_sampling`:

* `one_k` stops when `k/aH < ratio_max` **and** `|dlnPdN| < tol_curvature`.
* `evolve_background` stops when `aH`, `phi`, `a` or `d2a/dt2` reaches a target, in
  either time direction, with a per-step domain check on `epsilon`, `V` and `V'`.

Notably the backward direction is *not* the obstacle: every evolver except `evolver_rk`
already carries a `tdir` and handles `x_end < x_ini` **[C]**, so
`evolve_background`'s hand-written backward branch would simply disappear.

**The house pattern for this is to precompute the boundary and then call the evolver on
known intervals**, and it is used three times **[C]**:

* thermodynamics root-solves the RECFAST trigger crossings analytically
  (`thermodynamics_module.cpp:3166,3173`) and calls the evolver once per phase;
* perturbations builds `interval_limit[]` in `perturb_find_approximation_switches` and
  calls the evolver once per interval;
* background integrates to a precomputed `lna` grid.

### 8.4 The three options, costed

* **(a) Leave it.** The module is correct, converged (section 7) and takes 0.3 s for 900
  modes. Cost 0. What is given up: `ppr->evolver_*` cannot select an integrator here, so
  a stiff inflationary model has no ndf15 to fall back on.
* **(b) Add an event/stop callback to the shared evolver signature.** Six evolvers
  (header + implementation), four existing call sites passing `nullptr`, and a post-step
  check in each driver — cheap in `evolver_erk_impl.h`, which already has dense output at
  the accept point (line 271) so the crossing can be refined for free, and separately in
  the 1596-line `evolver_ndf15.cpp` and in `evolver_etd.cpp`. Gain: the two functions
  collapse to roughly one evolver call each, ~400 lines to ~150, the outer-step heuristic
  in `one_k` goes away with them, and every module gets event termination. This is the
  right shape, but it is a contract change to the hottest interface in the code and one
  module is thin motivation for it.
* **(c) Follow the house pattern.** Pre-integrate the background alone — three equations,
  cheap — to find the `tau` where `aH = k/ratio_max`, then one evolver call per mode over
  that known range with a sampling grid, checking at the samples that the curvature froze.
  No signature change, and it is what thermodynamics and perturbations already do. Costs a
  duplicated background solve per mode, and the `tau` search is itself an event problem,
  just a much cheaper and monotonic one.

**Recommendation: (a) now, (b) when a second module wants it.** (c) is the fallback if the
inflation module has to be modernised on its own. Nothing here should be done as part of
the axion-monodromy change.

## 7. Verification

Everything below is from `notebooks/axion_monodromy_flauger2009.py` **[M]**, 2026-09-06.

**The headline.** $\delta n_s$ measured from the computed $P_s(k)$, against the paper's
eq. (2.9), at $b=0.08$ over $f = 2\times10^{-4}$ to $0.1\,M_p$ (24 values):
numerical/analytic has **min 0.996, median 1.013, max 1.022**. The paper's own numerics
(its fig. 1) fall off the curve below $f\approx2\times10^{-3}$ and plunge to zero; CLASS++
tracks it a further decade down.

**The estimator, and its controls.** $\delta n_s$ is the amplitude of the projection of
$P_s^{b}/P_s^{b=0}-1$ onto $\cos(\phi_k/f+\beta)$ with $\phi_k$ from eq. (2.5), $\beta$
free, and a quadratic in $\ln k$ to absorb any smooth mismatch. $\phi_*$ is **scanned, not
fixed at 11**: eq. (2.5) is a slow-roll relation good to $O(\epsilon)\approx0.4\%$, which
at $f=10^{-3}$ is ~40 radians of accumulated phase across the fitting window — a fit at
fixed $\phi_*$ returns nearly zero and looks like a physics failure. This cost an hour;
the recovered $\phi_*$ runs 11.1–11.4, not 11.0.

Because the estimator maximises an amplitude over a nuisance parameter it can only bias
*upwards*, so two controls run beside it: the residual after subtracting the best fit
($\le 0.16$ of the amplitude, and $\le0.02$ wherever the scan is used) and the amplitude
the same estimator returns at $1/3$ of the true frequency ($\le 0.09$ of the amplitude).

**Convergence** at the hardest point used, $f=5\times10^{-4}$: doubling
`k_per_decade_primordial`, tripling `primordial_inflation_ratio_min`, taking
`primordial_inflation_tol_integration` to $10^{-7}$, `pt_stepsize` to 0.002, `bg_stepsize`
to 0.001 or `ratio_max` to 0.005 all move $\delta n_s$ by $<0.1\%$. The one exception is
the **default** `tol_integration = 1e-3`, which leaves a 3.7% mode-to-mode scatter in
$P(k)$ (residual 0.037 against an amplitude of 0.045).

**The resonance-depth floor.** A mode resonates with the oscillating background at
$-k\tau = 1/(2f\phi_*)$, and the module starts each mode at
$-k\tau=$ `primordial_inflation_ratio_min`, default 100. Plotted against
`ratio_min` $\times 2f\phi_*$ the three tested $f$ collapse onto one curve with a sharp
threshold at 1: at $f=2\times10^{-4}$ (resonance at $-k\tau=227$) the default returns
**5.7%** of the right answer, `ratio_min = 200` returns 24%, and `ratio_min = 300`
returns 101.9%. This is almost certainly the same mechanism behind the paper's own
small-$f$ breakdown.

**Where the last 1–2% is.** Varying $b$ at fixed $f=0.006$ splits the excess cleanly:
numerical/analytic $= 1.0065 + 0.123\,b$. The slope is the $O(b^2)$ term the
first-order derivation drops (1.0% at $b=0.08$); the intercept is a $b$-independent
0.65%, the slow-roll corrections dropped when the Hankel index was set to $\nu_0=3/2$
($\epsilon = 1/2\phi_*^2 = 0.41\%$ here). Both have the right size and the right scaling,
and the excess vanishing as $b\to0$ is the check that would have caught an error in
either the shape or the notebook.

**Downstream.** A lensed $C_\ell^{TT}$ at $f=0.01$, $b=0.08$ carries the modulation at up
to 12.7% (rms 7.2% over $2\le\ell\le2500$). Doubling the primordial $k$-grid and halving
`perturb_sampling_stepsize` changes that by 0.006%, so the spline of an oscillating
$P(k)$ onto the perturbation grid is not aliasing at this frequency.

**No regression.** The stock polynomial `inflation_V` case returns
`A_s = 2.18996e-09, n_s = 0.955167, alpha_s = -5.43475e-05` before and after, to every
printed digit; `make test` is 34/34.
