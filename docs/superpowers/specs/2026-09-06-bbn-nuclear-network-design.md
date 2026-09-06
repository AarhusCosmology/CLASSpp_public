# A nuclear reaction network for primordial abundances

**Status — 2026-09-06. SHIPPED, with three deviations from the design below.**
The solver is in `tools/bbn_{nuclides,reactions,rates,plasma,weak,solver}.*`,
selected by `YHe = network`, with `YHe = BBN` still the default and unchanged.
Four test executables cover it. Numbers formerly tagged `[P]`/`[?]` in the
acceptance table are now `[M]`; the design text below is left as written except
where a **MEASURED** note corrects it, so a reader can see what was predicted
against what came out.

Headline, measured over the 42 grid points of `bbn/sBBN_2017.dat` spanning
omega_b in [0.00998, 0.03493] and DeltaN in [-3, +7], at the table's own
tau_n = 880.2 s:

| Quantity | Target | Measured | |
|---|---|---|---|
| Y_p vs table, `rho_He4/rho_b` | < 5e-4 abs | mean +9.5e-5, rms 1.2e-4, **max 2.2e-4** | **[M]** |
| Y_p vs table, `4 Y_He4` | — | mean +1.4e-3, **max 1.6e-3** (excluded) | **[M]** |
| D/H at omega_b = 0.02242 | within 2% | **2.498e-5** (published span 2.44–2.58e-5) | **[M]** |
| He3/H | within 5% | **1.037e-5** (published ~1.04e-5) | **[M]** |
| Li7/H | reported only | **4.86e-10** (codes ~4.7–5.5e-10) | **[M]** |
| eta_10 | — | **6.124** at omega_b = 0.02237 | **[M]** |
| Runtime, one solve | assumed ms **[?]** | **19 ms** (0.29 s before section 11) | **[M]** |

**Deviation 1 — the QED plasma correction is NOT implemented.** Section 5 lists
it as included. It was dropped during implementation because no formulation could
be verified to the precision it would need, and shipping an unverifiable ~1e-4
correction is worse than shipping a documented omission. The measured residual
above is inside tolerance without it.

**Deviation 2 — neutrinos use exact instantaneous decoupling, N_eff = 3.0 +
DeltaN_eff, not 3.044 + DeltaN_eff.** Section 5's added paragraph argued for
3.044. That is measured WRONG against the table: with a 3.044 baseline the
worst-case Y_p residual degrades from 2.2e-4 to 1.23e-3, a factor of 5.6 **[M]**.
The physical reason is the one section 8 already gives — incomplete decoupling is
an energy TRANSFER from pairs to neutrinos, not an addition — so adding 0.044 as
extra radiation double-counts it. The table's `DeltaN_eff = N_eff - 3.046` axis
labels the CMB-era N_eff, not the BBN-epoch energy density.

**Deviation 3 — a REACLIB fit ceiling was needed, and is not in the design.**
See the new section 4.2. Without it the network cannot start where it must.

**One design justification did not survive measurement.** Section 3 chose the
coupled system over the two-stage split because only the coupled form carries
nuclear back-heating without adding it by hand. Measured by zeroing that term:
it is worth **2.9e-10 relative in Y_p and 4.2e-9 in D/H** **[M]** -- six orders
below the acceptance targets. The coupled structure is still what shipped, and it
is still the right shape (it matches the reference codes and costs nothing), but
the advantage claimed for it is not measurable and section 3 should be read with
that in mind.

Evidence tags below: **[M]** measured here, **[C]** stated by the code, **[P]**
from the literature, not re-measured, **[?]** no evidence.

## 1. What this replaces, and why

`YHe = BBN` today reads `bbn/sBBN_2017.dat`: a 48x11 grid of Y_p over
(omega_b, Delta N_eff), computed externally with PArthENoPE v1.2 for a neutron
lifetime of 880.2 s, bicubic-splined in `ThermodynamicsModule::thermodynamics_helium_from_bbn()`.
It works, and for LCDM it is accurate. Its limits are structural:

* **It is two-dimensional.** Any parameter that is not omega_b or Delta N_eff —
  the neutron lifetime, a reaction rate, an electron-neutrino chemical potential —
  cannot be varied, because it was frozen when the grid was computed.
* **It gives one number.** D/H, He3/H and Li7/H are the observables that most
  sharply constrain BBN, and the table cannot produce them at all.
* **It hard-fails outside its grid.** Four `class_test` calls reject omega_b or
  Delta N_eff outside the tabulated box, which is a real obstruction for a sampler
  exploring exotic radiation content.
* **Delta N_eff is a lossy summary.** A decaying species that injects energy
  *during* BBN is not describable by a single number, however the number is
  computed.

A solved network removes the first three outright and opens a seam for the fourth.

## 2. Scope

**In:** Y_p, D/H, He3/H and Li7/H from a solved 9-nuclide network, exposed as
CLASS derived parameters and through `classy`; an optional abundance-history
output; opt-in via `YHe = network`, with the table left as the default.

**Out, deliberately** — named here so a reader can tell a dropped section from an
implemented one:

* Thermal (finite-temperature) radiative corrections to the weak rates.
* The O(e^2) QED plasma correction — moved here from section 5 during
  implementation; see Deviation 1 in the status header.
* Full incomplete-neutrino-decoupling with spectral distortions. Carried only as
  `N_eff = 3.044` through an effective neutrino energy density.
* Electron-neutrino chemical potential `xi_nu_e`. (The existing table's own header
  comment flags this as neglected there too.)
* Driving H from the real CLASS background. The `rho_extra` seam exists and is
  exercised, but the prototype fills it from `Delta N_eff` only. See section 8.
* Nuclides beyond Li7/Be7; no CNO.
* Any change to the default. `YHe = BBN` behaves exactly as it does today.

## 3. Structure: one coupled stiff system

Two structures were considered. **Rejected:** a two-stage split solving the
plasma first and splining it into a stiff network solve. It gives a purely
nuclear Jacobian and independently testable halves, but it cannot carry the
back-heating of the plasma by nuclear energy release without adding it back by
hand.

**Chosen:** a single coupled stiff system, as Kawano, PArthENoPE and AlterBBN all
do. Back-heating is then structural rather than a correction, and the system is
directly comparable to the reference codes. The standard objection — that
coupling temperature to every abundance produces a dense Jacobian row and column,
defeating sparse grouping — does not bite at this size: the system is 12
equations, and a dense 12x12 LU costs nothing. The objection is real for a
100-nuclide network; it is not real here.

The one concern that does survive is the error norm, since abundances span twenty
decades. `evolver_ndf15` weights with `wt = max(|y|, 1e-15)` (`tools/evolver_ndf15.cpp:92,235`)
**[C]**. Everything validated here — Y_p ~ 0.25, D/H ~ 4e-5, He3/H ~ 1e-5,
(Li7+Be7)/H ~ 5e-10 — sits far above that floor, so the norm is genuinely relative
for all of them. Only Li6 (~1e-14) approaches it, and Li6 is not validated. This
is why abundances are evolved directly and **not** as `ln Y_i`: log variables move
the floor somewhere less convenient and misbehave at `Y = 0`.

### 3.1 The system

Independent variable `x = ln a`; twelve equations.

```
y = [ T_gamma, h, t, Y_n Y_p Y_d Y_t Y_he3 Y_he4 Y_li6 Y_li7 Y_be7 ]

    h ≡ n_b / T_gamma^3          so  eta = pi^2 h / (2 zeta(3))

dt/dx  = 1/H
dh/dx  = -3 h ( 1 + dlnT/dx )
dT/dx  = [ -3(rho + p) + 3 rho_b - n_b Σ_i m_i (dY_i/dx) ] / (∂rho/∂T)
dY_i/dx = (1/H) · [ network RHS, per unit time ]

H^2 = (8 pi G / 3) ( rho_gamma + rho_e± + rho_nu + rho_b + rho_extra(a) )
```

`rho`, `p` and `∂rho/∂T` run over the thermally coupled sector only: photons,
electron-positron pairs including the O(e^2) QED plasma correction, and baryons
(rest mass plus 3/2 n_b T). Neutrinos are decoupled and excluded from the energy
balance.

The `Σ_i m_i dY_i/dx` term **is** the nuclear back-heating, and it is the reason
this structure was chosen. **MEASURED: zeroing this term changes Y_p by 2.9e-10
and D/H by 4.2e-9 [M].** The reason is arithmetic that should have been done
before choosing: the binding energy released is ~1.7 MeV per baryon over all of
BBN, against a photon energy of ~4e8 MeV per baryon, so the plasma cannot notice
it. The structure is fine; the justification was not. Because `rho_b` is built from real nuclear masses,
binding-energy release requires no separate term: it appears automatically as
Σ m_i Y_i falls.

Weak-interaction energy transfer to the neutrino sector is neglected. It is small
**[P]** and this is an assumption, not a result.

### 3.2 The eta normalization pre-pass

The boundary condition on eta is at the **end** of the integration — the value
implied by omega_b today — while the solve needs `h` at the start. Since
`rho_b/rho_tot ~ 1e-6` during BBN, `T(x)` is independent of `h` far below any
tolerance here, so a one-variable pre-pass integrating `T` alone measures

```
R = (aT)^3_final / (aT)^3_initial      ⇒     h_ini = h_final · R
```

exactly. This is a scalar normalization pass, not a physics split: the main solve
remains fully coupled.

The alternative — assuming the textbook `11/4` — is wrong at the QED-correction
level, worth roughly 0.1% in eta, hence ~0.16% in D/H and ~4e-5 in Y_p **[P]**.
Those sit inside every tolerance in section 7, so this pre-pass buys accuracy that
is not strictly required. It is kept because it costs microseconds and because
comparing `R` against `(11/4)^(1/3)` is then a free unit test of the whole plasma
sector.

## 4. Supplying the network

This was the hard part of the design. The split is:

**Topology in code, forward rates in a data file, reverse rates never in either.**

* `tools/bbn_nuclides.h` — a `constexpr` table: name, A, Z, spin degeneracy
  `g = 2J+1`, mass excess.
* `tools/bbn_reactions.h` — a `constexpr` table of reaction topology: reactant and
  product nuclide indices with multiplicities. **Q-values are computed from the
  mass excesses, never stored**, so a Q-value cannot drift out of agreement with
  the masses used for detailed balance.
* `bbn/rates_<compilation>.dat` — forward rate coefficients only, in JINA
  **REACLIB** 7-parameter form:

  ```
  lambda = Σ_terms exp( a0 + a1/T9 + a2/T9^(1/3) + a3 T9^(1/3)
                           + a4 T9 + a5 T9^(5/3) + a6 ln T9 )
  ```

  REACLIB is the format the nuclear-astrophysics community publishes rates in, so
  a collaborator's compilation drops in without translation.

* **Reverse rates are always derived by detailed balance** from the code-side
  degeneracies, mass numbers and computed Q-values. They are never read from a
  file and there is no syntax for them. A hand-edited or mis-transcribed rate file
  therefore cannot violate thermodynamic consistency; the worst it can do is be
  wrong in one direction, which shows up as a wrong abundance rather than as a
  slow unphysical drift away from equilibrium.

* **Matching is by name, and absence is fatal.** A reaction present in
  `kReactions` with no entry in the file aborts; a file entry naming an unknown
  reaction aborts. There are no defaults and no skipping — a missing rate must
  fail loudly, not quietly produce a plausible number.

Selected by a new `bbn_rates_file` precision parameter, mirroring `sBBN_file`
including its `class_dir` prefixing.

### 4.1 The network

Nine nuclides — n, p, D, T, He3, He4, Li6, Li7, Be7 — and roughly fifteen
reactions, the classic Kawano key set: `p(n,g)d`, `d(p,g)he3`, `d(d,n)he3`,
`d(d,p)t`, `he3(n,p)t`, `t(d,n)he4`, `he3(d,p)he4`, `he3(a,g)be7`, `t(a,g)li7`,
`be7(n,p)li7`, `li7(p,a)he4`, `d(a,g)li6`, `li6(p,a)he3`, `t(p,g)he4`, plus free
neutron decay and the two-directional weak `n <-> p` rates.

This set reaches D/H at roughly the 1% level and Y_p at roughly 1e-4 **[P]**.
Heavier nuclides change nothing that is validated here.

### 4.2 The REACLIB fit ceiling (not in the original design)

REACLIB's seven-parameter form is a FIT, valid over roughly T9 in [0.01, 10], and
its `a4 T9` and `a5 T9^(5/3)` terms are fitting conveniences with no physical
content. Extrapolated to the T9 = 100 start of section 6 they produce `exp(207)`
-- a rate wrong by ninety orders of magnitude, which presents as the integrator
failing rather than as a bad rate **[M]**.

Starting the integration lower is not the fix: at T9 = 10 the weak rates are only
~1.6 times H **[M]**, so the n/p ratio is already leaving equilibrium and the
equilibrium initial condition of section 6 would be wrong.

What ships instead: the fit ARGUMENT is held at T9 = 10 above that, while the
Q-value and thermal factors in the reverse rate always see the true temperature.
This is safe because of what the rate does up there -- above T9 ~ 10 every heavy
nuclide is in nuclear statistical equilibrium, where the abundances are fixed by
detailed balance and the rate magnitude only sets how fast equilibrium is
restored. Verified rather than argued: Y_p and D/H move by less than 1e-5
relative as T9_initial goes 100 -> 200 -> 400 -> 1000 **[M]**.

## 5. Weak rates

Born `n <-> p` rates as two-temperature phase-space integrals over the electron
(T_gamma) and neutrino (T_nu) distributions — the two temperatures differ after
e± annihilation and must not be conflated — normalized so that the free-decay
rate reproduces the input neutron lifetime. On top of Born:

* zero-temperature radiative and Coulomb corrections,
* finite-nucleon-mass corrections,
* ~~the O(e^2) QED plasma correction to the electron-positron energy density and
  pressure~~ **NOT IMPLEMENTED — see Deviation 1. No formulation of it could be
  verified to the precision it would need, and an unverifiable 1e-4 correction is
  worse than a documented omission. The measured residual is inside tolerance
  without it; it is the obvious first thing to add if that ever stops being
  true,**
* incomplete neutrino decoupling carried as `N_eff = 3.044` via an effective
  neutrino energy density.

**MEASURED, and this paragraph is WRONG — see Deviation 2 in the status header.**
What shipped is exact instantaneous decoupling, `N_eff = 3.0 + DeltaN_eff`. A
3.044 baseline degrades the worst-case Y_p residual from 2.2e-4 to 1.23e-3 **[M]**.
The paragraph is kept, rather than deleted, because the reasoning it contains is
the reasoning that turned out to be wrong, and that is worth being able to see.

**The two N_eff reference values must not be conflated.** `sBBN_2017.dat` defines
its axis as `Delta N_eff = N_eff - 3.046`, and `thermodynamics_helium_from_bbn()`
subtracts exactly that **[C]**. This prototype keeps that interface — the
`Delta N_eff` handed to the solver retains the 3.046 reference, so the table path
and the network path are driven by an identical number — while the Standard-Model
baseline inside the plasma is 3.044. Total neutrino content is therefore
`3.044 + Delta N_eff`. The residual 0.002 is deliberate and physical, the
difference between the older and current incomplete-decoupling results, not an
off-by-one; it is worth about 1e-5 in Y_p **[P]**.

Approximate Y_p budget, from the literature, **[P] — not measured here**:

| Ingredient | ~ Delta Y_p |
|---|---|
| Born, normalized to tau_n | baseline |
| zero-T radiative + Coulomb | -4e-4 |
| finite nucleon mass | +1e-4 |
| QED plasma correction | -1e-4 |
| N_eff 3.046 -> 3.044 | -1e-5 |

These signs and magnitudes are the reason the target in section 7 is 5e-4 and not
1e-4: dropping them entirely would miss by more than the budget, and including
them should land inside it. Whether it does is what `bbn_solver_test` measures.

## 6. Initial and final conditions, and reported quantities

Start at T9 = 100 (T ~ 8.6 MeV). `Y_n` from weak equilibrium
`1/(1 + exp(Delta m_np / T))`, exact to well below tolerance there since the weak
rates vastly exceed H. Every heavier nuclide starts at its **NSE (Saha) value**,
computed with the same detailed-balance machinery that produces reverse rates —
reuse, not a second code path, and therefore covered by the same test. Integrate
to T9 = 0.01.

Reported, with conventions made explicit because they are a known source of
silent factor errors:

```
D/H    = Y_d / Y_p
He3/H  = (Y_he3 + Y_t) / Y_p       tritium beta-decays to He3 long after BBN
Li7/H  = (Y_li7 + Y_be7) / Y_p     Be7 electron-captures to Li7 long after BBN
```

### 6.1 The helium convention: an open question, now measured

CLASS uses `fHe = YHe / (_not4_ (1 - YHe))` with `_not4_ = m_He/m_H`
(`source/thermodynamics_module.cpp:3019`) **[C]**, i.e. it treats `YHe` as a mass
fraction. A network produces `n_i/n_b`, and the two candidate conversions

```
(a)  Y_p = 4 · Y_he4                          nucleon-number fraction
(b)  Y_p = Y_he4 m_He4 / Σ_i Y_i m_i          true mass fraction
```

differ through the He4 binding energy. The ratio is `m_He4 / (4 Σ_i Y_i m_i)`,
which for a realistic post-BBN composition is about 0.994 — a 0.5 to 0.8%
relative difference, i.e. **1.2e-3 to 1.9e-3 absolute in Y_p**. Against the 5e-4
tolerance of section 7 that is a factor 2.5 to 4, so `sBBN_2017.dat` does
discriminate between them, though not by a wide margin: the comparison must be
run at the table's own tau_n = 880.2 s for the discrimination to be trustworthy.

**RESOLVED BY MEASUREMENT, and (b) wins [M].** Over the table's whole grid:

| Convention | mean | rms | max abs |
|---|---|---|---|
| (b) `rho_He4 / rho_b` | +9.5e-5 | 1.2e-4 | **2.2e-4** |
| (a) `4 Y_He4` | +1.4e-3 | 1.4e-3 | 1.6e-3 |

A factor of seven apart, so the table discriminates cleanly after all. (b) is
also what the rest of CLASS assumes: recfast forms `fHe = YHe/(_not4_ (1-YHe))`
with `_not4_ = m_He/m_H`, a mass ratio **[C]**. `ThermodynamicsModule` is fed
`Yp_mass`; both are computed and returned, so the comparison stays reproducible.

## 7. Validation

Four test executables. Each must be registered in **both** `CMakeLists.txt` and
the Makefile `TEST_TARGETS` list — CI runs `make test` and only builds targets
named there.

| Test | Checks |
|---|---|
| `bbn_plasma_test` | `(aT)^3` ratio across e± annihilation vs `(11/4)^(1/3)` up to the QED correction; g\*(T) limits; t(T) against the analytic radiation-era result |
| `bbn_rates_test` | detailed-balance round-trip; the network relaxes to Saha/NSE when driven at fixed high T |
| `bbn_weak_test` | Born rate reproduces `1/tau_n` as T -> 0; n/p ratio -> `exp(-Delta m/T)` at high T |
| `bbn_solver_test` | Y_p against `sBBN_2017.dat` over its grid; D/H, He3/H, Li7/H against published values |

**Acceptance — MEASURED 2026-09-06 [M]:**

| Quantity | Target | Measured | Comparison |
|---|---|---|---|
| Y_p | within 5e-4 absolute | **max 2.2e-4** over 42 grid points | `bbn/sBBN_2017.dat` at tau_n = 880.2 s, that file's stated assumption |
| D/H | within 2% | **2.498e-5** | published span 2.44–2.58e-5 **[P]** |
| He3/H | within 5% | **1.037e-5** | published ~1.04e-5 **[P]** |
| Li7/H | reported, not gated | **4.86e-10** | codes only, ~4.7–5.5e-10 **[P]** |

Two further measurements worth keeping:

* **The (aT)^3 entropy transfer converges to (11/4)^(1/3) = 1.4010197 from
  BELOW, at second order in the starting temperature.** Fractional deficit
  -8.09e-5 at T9_initial = 100, -2.02e-5 at 200, -5.05e-6 at 400, -7.86e-7 at
  1000: a factor of four per doubling **[M]**. Below, not above, because at a
  finite start the pairs are not perfectly relativistic and slightly less
  entropy is available to transfer.
* **In a default CLASS run the network sits 1.0e-4 from the table** at
  omega_b = 0.02237 when both use tau_n = 880.2 s. Left at the shipped default
  of 878.4 s (PDG 2022) the gap is 3e-4, essentially all of it the deliberate
  change of neutron lifetime at dY_p/dtau_n ~ 2e-4 per second **[M]**.

Two rules govern how these are read:

* **No code is an oracle.** PArthENoPE, PRIMAT and AlterBBN are independent
  implementations of contested physics. Agreement raises confidence; disagreement
  opens an investigation. None of them is ground truth, and neither is this one.
* **Li7 will disagree with observation by roughly a factor 3.** That is the
  cosmological lithium problem, not a bug in this code. Li7 is therefore compared
  only against other codes, never against data, and the test documents this so a
  future reader does not "fix" it.

**Convergence.** Self-convergence in `tol_bbn_integration`, `bbn_T9_initial` and
`bbn_T9_final` is reported as an order measured from a doubling **triple**. A
single refine pair `|f(2h) - f(h)|` is the *coarse* run's error, not the fine
one's, and is not accepted as an error bar here.

## 8. The rho_extra seam

`H` takes an additive `rho_extra(a)` supplied by the caller. The prototype fills
it from `Delta N_eff` alone, reproducing exactly what the table's interface
offers. The seam exists so that the day a real `rho_dncdm(a)` or a dark-radiation
density should drive BBN, the interface does not have to be reopened.

The reason this is a seam and not a wire is worth recording: **CLASS's background
cannot supply H during BBN.** It contains no electrons or positrons and pins
`T_gamma = T_cmb/a` with a fixed N_eff at all times, so its H at T ~ MeV is wrong
by a large factor. Connecting them means reconciling two inconsistent `T(a)`
relations through comoving entropy — real work, easy to get subtly wrong, and out
of scope here.

## 9. Integration into CLASS

* **Input.** `YHe = BBN` keeps the table, unchanged. `YHe = network` runs the
  solver. New `tau_n` input, default 878.4 s; 880.2 s reproduces the table.
* **Precision.** `bbn_rates_file`, `tol_bbn_integration`, `bbn_T9_initial`,
  `bbn_T9_final`, `bbn_history_file`, and `weak_table_points` on `BbnInput`.
* **Evolver.** **`evolver_bbn` was NOT implemented**, though this section
  originally listed it and described the selection pattern. `bbn_solver.cpp`
  calls `evolver_ndf15` directly.

  Nothing was lost by dropping it. Section 11 measures the residual stiffness
  after the Jacobian diagonal is removed at ~1e10 per unit tau **[M]**, so every
  explicit evolver in the tree would need ~1e10 steps where ndf15 needs a few
  hundred: there is no alternative worth selecting. Adding the knob would also
  mean joining the `evolver_erk_configure` process-wide-state protocol that the
  three other evolving modules each re-implement with the same warning comment,
  for a choice with one correct answer.
* **Exposure.** `ThermodynamicsModule` gains public `DoverH_`, `He3overH_`,
  `Li7overH_` alongside the existing `YHe_`. Public members are what
  `generate_wrapper.py` lifts into `cclassy.pxd` — which is generated, never
  hand-edited. `classy.pyx` `get_current_derived_parameters` gains the matching
  names; a verbose line and a thermodynamics-output header line carry them for
  non-Python users.
* **Abundance history.** Optional `bbn_output`, writing `Y_i(T9)` trajectories
  through the evolver's existing dense-output sampling hook. Cheap, and the
  natural artifact for inspecting freeze-out.

## 10. Cost

A 12-equation stiff solve over roughly four decades in temperature.

**MEASURED: 0.29 s per solve as first written [M]**, not the milliseconds guessed
here — a thousand times the estimate. Section 11 brought it to **19 ms**.

In situ, `background + thermodynamics` alone costs 12.9 ms including process
start; with `YHe = network` the same run costs 36.1 ms **[M]**. So the network
went from roughly 23x the entire background-and-thermodynamics chain to about
1.8x it. It is no longer the dominant cost of a CLASS run, and it is not free
either.

---

## 11. Making it fast (2026-09-06, after the first implementation)

The first working version took 0.29 s. Profiling rather than guessing, per
right-hand-side evaluation **[M]**:

| Component | Time | Share |
|---|---|---|
| weak rates (6 x 256-node quadratures) | 21.7 us | **87%** |
| plasma (128-node quadrature) | 2.0 us | 8% |
| reverse rates (14 reactions) | 1.0 us | 4% |
| forward rates (14 reactions) | 0.24 us | 1% |

and 38 700 right-hand-side calls, spread over **ten** evolver runs.

**The two-stage split of section 3 would not have helped.** It removes the plasma
recomputation and nothing else -- 8%, so at most 1.09x -- because the stiff
network solve still needs the weak rates at every step however the plasma is
obtained. Code-wise it is a wash or slightly worse, since it adds splines to
manage. What makes the split *look* attractive is the thing actually worth
having: a known trajectory to precompute along. The eta normalization pre-pass
already computes exactly that, so the benefit was available without the
restructuring.

Four changes, in descending order of what they bought:

1. **Splined weak rates (`BbnWeakTable`), 11.4x.** Built on the pre-pass
   trajectory, which supplies the T_nu that makes the rates a function of the
   integration coordinate alone. Splined in ln(lambda), because lambda_pn falls
   through seventy decades. At the default 512 nodes this shifts Y_p by 1.8e-10
   and D/H by 2.1e-9 **[M]**; past ~512 the difference stops falling because it
   is then the evolver's own rtol noise rather than interpolation error. The
   reference path (`weak_table_points = 0`) is kept so the comparison is measured
   and not assumed.
2. **Ten evolver runs down to three, ~1.8x.** The eta pre-pass was inside the
   composition-refinement loop although it is eta-independent to 1e-6, so it ran
   five times instead of once. The loop itself converged a 0.18% correction to a
   1e-12 criterion from an all-hydrogen starting guess; it now starts from a
   hydrogen-plus-helium estimate and stops at 1e-8, which pins eta far below
   every other error here.
3. **Detailed-balance constants precomputed.** `C * T^(3/2 surplus) * exp(-Q/T)`
   factorizes, and only the last two factors carry temperature — so the
   degeneracies, masses and unit conversion are formed once in the constructor
   instead of costing a `pow()` per nuclide per reaction per step.
4. **One REACLIB evaluation per reaction.** `Reverse` re-evaluated the forward
   fit that the caller had just computed.

Net: **0.29 s -> 19 ms, 15x [M]**, with the Y_p-versus-table comparison
unchanged to every digit (mean +9.52e-5, rms 1.16e-4, max 2.20e-4).

What was NOT done, and why: dropping the weak quadrature from 256 to 64 nodes
would have given ~2.9x on its own for a Y_p error of ~3e-7 **[M]**, which is
harmless. It is unnecessary now that the integrals are evaluated ~512 times per
solve instead of ~38 700, and the accuracy is worth more than the microseconds.
