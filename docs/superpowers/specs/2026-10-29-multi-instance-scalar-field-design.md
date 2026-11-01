# Multi-instance scalar fields (#346)

**Status (2026-10-29):** shipped as described. Coupled multi-field models, and
closing the budget with a dot-syntax field, were left out on purpose (see
"Not done").

**Revised 2026-10-31, on the rebase onto #493:** Pass 1 of `ConstructSpecies`
now runs every factory, with no scalar-field exception (the "Closure" bullet
below). The verbose `Omega = inf` is fixed. The gauge disagreement under "Not
done" turned out to be under-resolved oscillations, not a gauge limitation.
Later the same day the module-level scalar-field checks for output moved into
the species, and a column title two columns share became an error ("Output").

## Problem

`ScalarFieldSpecies` was single-instance: one field, configured by the legacy
`*_scf` keys and registered as `"ScalarField"`. Models with several independent
fields (multi-field quintessence, an early axion plus late quintessence,
assisted-quintessence-like set-ups) could not be built in one run. #344 left
`scalar_field` out of the single-instance dot-syntax translation to keep this
option open.

## Approach

Each `<N>.type = scalar_field` builds one more field, registered under the key
`N`. The legacy field is unchanged and still lives under `"ScalarField"`.

* **One build body for both spellings.** `ScalarFieldSpecies::InputKeys` names
  the keys a field is read from: `Legacy()` gives the `*_scf` spelling,
  `Instance(N)` gives `N.Omega`, `N.parameters`, `N.tuning_index`,
  `N.shooting_parameter`, `N.attractor_ic`, `N.potential` and, for the axion
  potential, `N.m`, `N.f`, `N.n`, `N.Theta_initial` (the legacy key without its
  `_scf` / `_axion` suffix, as the ncdm dot fields drop `_ncdm`).
  `CreateConfigured` (default potential) and `CreateAxion` take the keys, so the
  two spellings cannot drift apart. Error messages name the key the user
  actually wrote.
* **Per-instance state comes for free.** Background, integration and
  perturbation indices were already per-object members registered through the
  species hooks, so a second object gets its own `index_bg_*`, `index_bi_*`,
  `idx_phi*` and transfer-source slots. Background and perturbation sums
  already loop over the collection. The one ordering rule, scalar fields last
  in the Newtonian-gauge dispatch, already used `dynamic_cast` and so covers
  every field.
* **Shooting.** Each field with `N.Omega > 0` and no `N.shooting_parameter`
  reports the target `{N.Omega, N.shooting_parameter}`. `DoShooting` already
  solves any number of targets jointly and writes each resolved unknown back
  under its own key.
* **Closure.** Only the legacy field can close the budget (`Omega_scf < 0`).
  Dot-syntax fields are always explicit. `ConstructSpecies` used to skip the
  closure species' factory in Pass 1, which would also have skipped the dot
  fields when the legacy field closes. Pass 1 now runs every factory without
  the closure override, and a factory builds its closure instance only when
  handed the override, and then nothing else. Lambda and the fluid already
  behaved this way, since their closure is the absence of their Omega key and
  that reads as 0 in Pass 1. So in Pass 1 the scalar-field factory builds the
  dot fields and skips a negative-`Omega_scf` legacy field; in Pass 2 it builds
  the legacy field alone. The dot fields' Omega is therefore part of the sum the
  closure value is computed from. (The 2026-10-29 version special-cased the
  scalar field's factory name in Pass 1 instead; the generic rule names no
  species.)
* **Output.** Column names carry a tag: `scf` for the legacy field, so its
  output is unchanged, and `N` otherwise (`phi_N`, `(.)rho_N`, `d_N`,
  `t__N`, `delta_N`).
  * **Repeated titles are an error.** A tag is a user's name, so a field named
    `b` would write a second `(.)rho_b`, and one named `tot` a second
    `(.)rho_tot` (and `d_tot`). classy's dicts and the readers of the `.dat`
    files key on the title, so one column silently stood for the other. The
    ncdm dot instances (`(.)rho_<N>`) had the same hole. Every title string
    a species writes into (background, transfer functions, `k_output_values`
    scalars and tensors) is now checked once it is complete
    (`first_duplicate_title`), as a severe error naming the title. The check
    sits where the titles are made, so it costs nothing in a run that writes
    none.
  * **No scalar-field checks in the modules.** The `background_verbose > 2`
    block and the `background.dat` units note used to find the field by a
    `dynamic_cast` in the background and output modules. Now the species does
    both itself. `PrintBackgroundDetails(pvecback_today, rho_crit_today)` is a
    print hook beside `PrintNeffInfo`/`PrintMassInfo`, and composites forward
    it. The units note goes through `BackgroundColumnWriter::Note`, next to
    the `phi`/`V` columns it describes, and a note several fields give is kept
    once. The details block prints each potential's own parameter names
    (`[m [1/Mpc], f, n, Theta_initial]` for the axion, where it said
    `[lambda, alpha, A, B]`), with `%g` instead of `%.3f`, which showed a shot
    axion mass of 1.5e-4/Mpc as 0.000. It drops the Omega_Lambda line, which
    the budget table already prints at this verbosity. Type 3 now gets both,
    since it forwards to its field. The one check left is the ordering
    of scalar fields last in the Newtonian-gauge dispatch
    (`perturbations_module.cpp`). It is physics in the derivative loop, not
    output, and the alternatives are a new species property or a second pass
    over the species per derivative evaluation.
* **Type 3 (`scf_veta`).** The composite still owns the legacy field. The guard
  that stops the standalone factory from building a second one now covers the
  legacy field only, so dot fields can sit beside the composite.
* `ScalarField` and `scf` are reserved as instance names, since they would
  collide with the legacy field's key and column tag. `N.Omega` is required (severe). A negative value is rejected
  as a computation error, because it is a numeric value a sampler can propose
  (STYLE.md §8). `N.Omega = 0` builds no field, as `Omega_scf = 0` does.

## Alternatives rejected

* **Translating the legacy keys into a synthetic dot instance**, the way
  `N_ncdm` becomes `ncdm__1.*`. This would rename the legacy species key and
  its columns, breaking `ScalarField.f_peak`, the Type-3 composite's
  construction and every existing `background.dat` reader. Keeping the legacy
  path as one spelling of the same build body gives the same code sharing with
  no behaviour change.
* **Allowing a dot-syntax closure field** (`N.Omega < 0`). Two-pass
  construction would then have to pick which instance closes the budget and
  defer only that one, for no model that needs it: one closing field is enough,
  and the legacy keys already provide it.
* **Coupled fields.** Out of scope (#346). Each field sees only the metric.

## Verification

All numbers measured on 2026-10-29 **[M]**.

* Unit test `test-scf-multi-instance`: keys, column tags and shooting targets
  per instance; a resolved `N.shooting_parameter` disarms shooting;
  non-attractor initial conditions; Pass 1 and Pass 2 behaviour with a legacy
  closure field; rejections. Through `InputModule`: two fields beside Lambda,
  and a dot field beside a legacy closure field (the budget sums to 1).
* **Legacy is byte-identical.** Old and new binaries, comparing
  `background.dat`, `cl.dat`, `cl_lensed.dat`, `pk.dat` and the full verbose
  log, give identical files for a closure quintessence run (`Omega_scf = -1`)
  and a shot axion run (`Omega_scf = 0.05`, `scf_potential = axion`).
* **A dot field equals the legacy field.** The same axion as `ax.*`: identical
  shooting log. Spectra are bit-identical when the instance sorts where
  `"ScalarField"` does (name `ScalarFieldX`). Under the name `ax`, which sorts
  after the other species and so reorders the perturbation vector, TT differs
  by ≤ 2e-5 for ℓ < 1000 and ≤ 5e-5 up to ℓ = 2500. That is evolver-tolerance
  noise from the reordering.
* **Two fields add up.** In θ = φ/f, the axion equation does not depend on f,
  so two identical fields with f/√2 at half the density each are one field
  with f. With a fixed mass (m = 0.1/Mpc, f = 0.3, Θ = 2), the field is
  Omega ≈ 0.075, oscillates from z ~ 100, and changes TT by O(1) against
  no field. The split agrees with the single field to:

  | gauge | TT | EE | φφ | P(k) |
  |---|---|---|---|---|
  | synchronous | 2.5e-4 | 2.1e-4 | 1.2e-5 | 3.8e-4 |
  | Newtonian | 4.6e-5 | 1.0e-4 | 7.3e-6 | 8.1e-4 |

  (maximum relative difference, ℓ ≤ 2500). With shooting (Omega 0.025 + 0.025
  against 0.05, f = 0.5), the two-target Newton solve converges and TT agrees
  to 1.2e-4.
* `test/dotsyntax_scalar_fields.ini` (axion + exponential quintessence, both
  shot) converges in ten builds, the discovery build included, to Omega = 0.05
  and 0.02.

Re-measured on 2026-10-31 after the rebase onto #493 and the revisions above
**[M]**, against a master binary from the same build directory:

* Legacy runs are byte-identical to master in `background.dat`, `cl.dat`,
  `cl_lensed.dat` and `pk.dat` for LCDM, a closure field, the shot tracker of
  `test/scenarios/gauge_scf.ini`, a shot axion, a fluid closing beside a fixed
  Lambda, the pheno-axion fluid and Type 3. The verbose log differs only in the
  fixed `Omega = inf` lines.
* The same axion under the dot name `ScalarFieldX` is bit-identical to the
  legacy field; under `ax` TT differs by 7e-5. Two shot axions with f/√2 at
  Omega 0.025 each match one shot axion with f to 8e-5 in TT (both gauges).
  Synchronous and Newtonian gauge agree to ≤ 1.2e-4 in TT for the legacy axion,
  the dot axion and the split, and to 2.2e-4 for
  `test/dotsyntax_scalar_fields.ini`.
* After the output changes the legacy files stay byte-identical, the
  `background.dat` header included, except the Type 3 header, which gains the
  units note. The 43 inis under `test/` and `explanatory.ini`, run with
  background, transfer and `k_output_values` output switched on, raise no
  repeated title (`dncdm_dr_bare_omega` fails as it should, on its own
  check).

## Not done / found on the way

* **Under-resolved oscillations, not a gauge limitation.** The 2026-10-29
  version reported that single fields disagree between gauges, 6% in TT. The
  disagreement is real but appears only for the fixed-mass axion (m = 0.1/Mpc):
  12% at ℓ = 2, 1% at ℓ = 10, < 1e-4 from ℓ = 50 **[M]**. The shot axion
  (Omega 0.05, f = 0.5) and the tracker of `gauge_scf.ini` agree to 1e-4 in TT
  at every ℓ **[M]**, as #348 established for the tracker. The field oscillates
  with period Δln a ≈ 2πH/m ≈ 0.014 today, and the default
  `back_integration_stepsize = 7e-3` samples a period twice. With
  `back_integration_stepsize = 1e-3`, `tol_perturb_integration = 1e-7` and
  `perturb_sampling_stepsize = 0.01`, the gauges agree to 8e-5 at ℓ = 2, and both
  move: synchronous by 10% at ℓ = 2, Newtonian by 1-3% at ℓ = 3-20 **[M]**. No
  single knob is enough (the background step alone leaves 0.9% at ℓ = 2).
  A field with m ≫ H0 needs that resolution in either gauge, or the KG-to-fluid
  switch closed as YAGNI in #365. This change does not touch it.
* The `background_verbose > 2` details block printed `Omega = inf`: it read
  `pvecback`, the evolver's workspace, which holds whatever time the evolver
  evaluated last. It now reads today's row of the background table. The values
  this exposes include a known open problem: the explanatory.ini tracker as
  closure field (`Omega_scf = -1`, λ = 10) reaches Omega_scf = 0.030 against a
  wished 0.688, because a closure field is not shot.
