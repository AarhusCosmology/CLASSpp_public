# Shooting: scaled tolerances, reuse of the converged build, and the decaying-NCDM budget

**Status (2026-09-29).** Shipped on branch `fix-dncdm-deg-budget`. One PR (the user's
choice: review and CI overhead). Departures from the design as first written are marked
*Shipped:* where they occur:
- The reserve has its own internal key.
- The seed switches at Γ = H0 alone.
- There is no per-target residual-scale field.
- A fifth defect was found and fixed (E): the integrated representation's two momentum
  grids, now guarded.

Measured results are in *Verification*. Tags: **[M]** measured, **[C]** checked in the
code, **[P]** predicted, **[?]** open.

## Problem

Four defects, found while explaining a 19 lnL "jump" at Γ = 0 in a CLiENT profile of the
A1/NO decaying-neutrino proxy (hpc_client, 2026-09-29).

1. **A `deg`-normalised `ncdm_decay_dr` is left out of the flatness budget.** [C][M]
   `CompositeSpecies::GetOmega0` sums the children. The parent is normalised by `deg` and
   never calls `SetOmega0`, so it reports 0. `DarkRadiationSpecies` and `DrPsdSpecies`
   report 0 by design. The closure (`input_module.cpp`, pass 1) then gives Ω_Λ the whole
   sector, and the universe over-closes by the sector's density today:
   H(0) = H0 √(1 + Ω_sector). This holds for all three representations. For m = 0.05 eV,
   one thermal ν_l and φ empty (the proxy run's configuration):

   | Γ12 = Γ / 10¹² km/s/Mpc | Ω_sector today | H(0)/H0 − 1 |
   |---|---|---|
   | 0 | 1.172e-3 | 5.9e-4 |
   | 1e-10 | 6.9e-4 | 3.5e-4 |
   | ≥ 1e-6 | 3.0e-5 | 1.5e-5 |

   At Γ = 0 this is a 0.04 km/s/Mpc relabelling of H0, worth 3.5 lnL on
   CamSpec + DESI + Pantheon+ at typical points. Correcting H0 by hand brings the proxy at
   Γ = 0 onto the ΛCDM network to −0.15 ± 0.23 over 189 points. [M]

2. **Shooting tolerances are 1e-3 absolute.** [C][M] `DoShooting` calls
   `fzero_Newton(..., 1e-3, 1e-3, ...)`. Until #271 a single target went through a
   bracket hunt and Ridders' method, stopping at 1e-5 × |x| (1e-5 relative on the
   unknown). #271 sent every target through Newton with the old n-D values. class_public
   reads `tol_shooting_deltax = 1e-4` and `tol_shooting_deltaF = 1e-6` from the precision
   structure. On master, `100*theta_s` = 1.04110 / 1.04200 / 1.03900 land at 1.041067 /
   1.041970 / 1.038959: 3–4e-5 off, about 0.1σ of Planck. [M] Any target in Ω units at
   the 1e-3 level "converges" at almost any guess.

3. **The converged evaluation is built twice.** [C][M] `DoShooting` returns a fresh
   `InputModule`, and the `Cosmology` then rebuilds every module the residual already
   built. The verbose log of a `100*theta_s` run prints the converged background and
   thermodynamics twice.

4. **`dr_representation = proxy` has no shooting hooks**, and the `Omega_dncdmdr` guess
   is wrong below Γ ≈ H0. [C][M] With the proxy, `Omega_dncdmdr` and
   `Omega_ini`/`Neff_ini` are read but never shot. At Γ = 0 the constructor's unshot guess
   gives ΔN_eff = 95.8, and the run dies in BBN. `DegGuessFromOmegaToday` below Γ = H0
   divides today's (matter) target by a radiation-scaled density per unit `deg`. For the
   integrated representation at Γ = 0 that guess is ~94× too large, and Newton then
   diverges to N_eff = −5251.

## Design

### A. Every target carries its scales

`DoShooting` solves for x̃ = x / unknown_scale with F̃ = F / residual_scale, and transforms
the seed Jacobian consistently: dx̃/dF̃ = (dx/dF) · residual_scale / unknown_scale. The two
tolerances become precision parameters with class_public's names and defaults:
`tol_shooting_deltax = 1e-4` on the scaled unknowns and `tol_shooting_deltaF = 1e-6` on
the scaled residuals. One tolerance then means the same for every target.

*Shipped:*
- **Unknowns.** Each is scaled by its seed's magnitude (1 for a zero seed), automatically,
  so no species code is involved.
- **Residuals.** No `residual_scale` field. Every existing target (θ_s, `dcdm_dr`,
  `dcdm_wdm`, `scalar_field`, `scalar_tensor`, the decay sector) already returns its
  residual in its natural unit, Ω or 100θ_s, so the scale is 1 by convention. The
  convention is written into `precision.h`. [C]

- **Residual scale** is the unit in which an error of 1 is an order-unity *relative*
  change of what the data measure. For `100*theta_s` that is 1. For a density target or
  closure it is also 1, in units of the critical density. The error enters H² as ΔΩ
  whatever the species' own size, so an absolute ΔΩ is the physical measure.
- **Unknown scale** is the magnitude of the seed. For `h` that is 1. For `deg` or a
  reserve Ω it is the seed itself, so `tol_shooting_deltax` is relative there.

Only one solver path exists (Newton). The seed Jacobian `dxdF` enters it only as the
first Jacobian probe, delx = −dxdF · F(x0). That is why the transformation above is all
the scaling it needs. [C]

### B. The converged evaluation is the production build

`ShootingResidual` builds a whole `Cosmology`, and the target pulls whatever modules it
needs (background; thermodynamics for θ_s; perhaps perturbations for a future target).
The workspace keeps the last *successful* evaluation's `Cosmology`, together with the
exact strings it wrote for the unknowns (`%.17g`). After the solve, `DoShooting` formats
the solution the same way. If the strings match, it returns that `Cosmology`, and
`Cosmology::GetInputModule` adopts it whole: the input module and every module already
built. Anything not yet built is built lazily from the adopted input module, as before.

This is exact, not approximate. Both the evaluation and the production path build from
the same file content: the same keys, the same unknown strings, `is_shooting` set. The
test is on strings, so it is also solver-agnostic. When the solver stops on the step
criterion at a point it never evaluated, the strings differ and the old rebuild happens.

Every getter calls `GetInputModule()` before testing its own pointer, so adoption
happens on the first call, whichever getter it is. (`main/class.cpp` evaluates the
`OutputModule` arguments in unspecified order.)

*Shipped:*
- The residual function frees the previous evaluation before building the next, so at
  most one extra `Cosmology` is alive.
- `ShootingResult` lives at namespace scope, because `generate_wrapper.py` ends a class
  at its first `};`.
- `InputModule::shooting_evaluations_` and `shooting_build_kept_` record what happened.

### C. One decay-sector base; the budget closes in every mode

A `DecayingNcdmSector : CompositeSpecies` holds the parent pointer and implements the
four hooks once, for the integrated, psd and proxy composites. The sector's density is
`Rho(bg_today)`: the composite already sums its children. There are three modes:

| mode | set by | unknown | residual (Ω units, scale 1) |
|---|---|---|---|
| combined | `Omega_dncdmdr` | `deg` | Ω_sector(today) − target |
| initial | `Omega_ini` / `Neff_ini` | reserve `Omega_dncdmdr_reserve` | reserve − Ω_sector(today) |
| **deg (new)** | `deg` | reserve `Omega_dncdmdr_reserve` | reserve − Ω_sector(today) |

*Shipped:* The reserve has its own internal key, written only by the shooter; a user
setting it gets a severe error. Reusing `Omega_dncdmdr` would be read as a combined-mode
target by the parent's constructor in any `deg`-mode shooting build without an explicit
`deg` (the default deg = 1), which would then re-guess `deg`. The initial mode moved to
the same key, and the constructor's `is_shooting` exception for `Omega_dncdmdr` alongside
`Omega_ini` went with it. [C]

The deg and initial modes are the same fixed point: the reserve feeds only Ω_Λ, and Ω_Λ
feeds the sector only through H(a) late on. A change of 1e-3 in the reserve moves the
sector's density today by ~1e-3 × Γt × δt/t, so the map contracts by ≲ 1e-3. [P]
dF/dR ≈ 1 is therefore the seed Jacobian.

**Seeds decide only the cost, never the result**: the residual is always evaluated. A
seed within `tol_shooting_deltaF` costs one build, and B makes that build the production
run. The seeds are:

- **R_stable** (Γ ≤ H0): the sector with no decays, i.e. the parent's quadrature density
  at a = 1 plus the daughters' initial populations redshifted as radiation. Exact at Γ = 0.
- **E_∞ · R_rad** (Γ ≥ H(a_nr)): R_rad is the whole sector's initial energy redshifted
  as radiation. With inverse decays, E_∞ = 1.233: the parent is held in chemical
  equilibrium until it goes non-relativistic, a limit that is independent of Γ and of the
  mass (measured on converged psd backgrounds, m = 0.06 and 0.3, 2026-08-07). Without
  inverse decays, E_∞ = 1: a relativistic decay conserves energy.
- **Between the two**, R_stable. This is also where the background is cheap.

*Shipped:* One switch at Γ = H0: R_stable below it, E_∞ · R_rad above it. The seed only
has to be good where the background is expensive, and that is far above H0. Between H0 and
H(a_nr), E_∞ · R_rad is the closer of the two (S/R_rad = 2.2 at Γ12 = 1e-8), and each
evaluation there costs ≤ 2 s. [M]

What the seeds buy, for the proxy's configuration: [M]

| Γ12 | S / R_rad | background wall time, 1 thread |
|---|---|---|
| 0 – 1e-12 | 47.6 – 47.3 | 0.2 – 1.4 s |
| 1e-10 – 1e-8 | 28.1 – 2.2 | 1.5 – 1.6 s |
| 1e-6 | 1.246 | 1.9 s |
| 1e-4 – 100 | 1.2235 – 1.2254 | 12 – 57 s |

Wherever the background is expensive, the seed 1.233 · R_rad is within ~2e-7 of the
true value, so it is one build. Where the seed is poor, a background takes ≤ 2 s.

Each composite reports its daughters' initial radiation through one narrow virtual.
Integrated: 0 (the daughter fluid starts empty). psd: the `DrPsdSpecies` children's
f_ini populations. Proxy: `DaughterRho` of its seeded coarse PSDs.

### D. The combined-mode guess below H0

Below Γ = H0, `DegGuessFromOmegaToday` divides the target by the parent's stable density
per unit `deg` *today*. The high-Γ branch is unchanged.

*Shipped:* Below H0, `DNCDMSector` first subtracts the daughters' initial populations from
the target, since `deg` does not scale them. At Γ = 0 the guess is then exact: one build for
the integrated representation and the proxy alike. Before, it took three, probing
N_eff = −5249 on the way. [M]

### E. The integrated representation's two grids (found while testing)

With the default quadrature (`qm_auto`), a decaying parent has 9 background and 3
perturbation momentum points. `DNCDMSpecies` evolves ln f on the perturbation grid, but its
background quadrature runs over all 9 bins, and the 6 extra bins keep their initial weights.
The parent's density therefore never falls, and `DNCDM_DR_Species` sources the decay
radiation from it indefinitely. [M]

| Γ [km/s/Mpc] | Ω_DR today, two grids | Ω_DR today, one grid (16 bins) |
|---|---|---|
| 1e4 | 0.041 | 4.1e-5 |
| 1e6 | 0.697 | 1.3e-5 |

The parent itself holds only Ω ≈ 1.2e-3. psd and proxy already refused two grids with a
severe error. The integrated representation now does too, with the same message (the
user's choice over silently forcing one grid). `test/scenarios/dncdm_dr_combined.ini`
reached its target only through this bug. With one grid it needs N_eff ≈ 11, beyond BBN's
helium table, so as a closure smoke test it now fixes `YHe`.

## Behaviour changes

- `100*theta_s` shoots now land within 1e-6 (they were 3–4e-5 off).
- Every `deg`-mode `ncdm_decay_dr` run gets Ω_Λ smaller by the sector's density today.
- The proxy honours `Omega_dncdmdr` and `Omega_ini`/`Neff_ini`.
- Cost:
  - No targets: unchanged.
  - A `100*theta_s` shoot sheds its duplicate build.
  - A high-Γ `deg` run is one build. [P]

## Verification

`test-shooting` (C++) was written first and watched to fail, 12 failures, each for the
reason named in its comment. Every case now passes: [M]

| case | result |
|---|---|
| `100*theta_s` = 1.04200 | 1.04200007 (was 1.04196953) |
| H(0)/H0 − 1, `deg` mode: integrated at Γ = 0 / 1e6, proxy at Γ = 0 | 2e-16 / 4.3e-7 / 1.6e-10 (were 5.8e-4 / 6.6e-6 / 5.9e-4) |
| integrated at Γ = 0 | one evaluation, and it is kept |
| proxy at Γ = 1e6 km/s/Mpc | one evaluation, kept, 1.5e-7 |
| combined mode at Γ = 0, integrated and proxy | on target, one evaluation each (the proxy was ignored before: Ω = 0.110 for 1.17e-3) |
| proxy with `Neff_ini` | budget closed |
| integrated representation with two grids | refused as severe |

The rest of the verification: [M]
- **C++ suite:** 47/47 `ctest` pass.
- **Python, shooting and decay tests** in `test_class.py`: pass. Four inputs needed a
  single grid.
- **Python, other non-scenario tests:** pass, except `test_rs_drag_matches_reference`,
  which fails identically with master as the candidate. It is the local
  fast-math-against-no-fast-math reference pair, not this change.
- **Scenario matrix:** the whole TEST_LEVEL=2 matrix, as CI runs it
  (COMPARE_OUTPUT_REF against a master `classyref`), passes 1512/1512. That includes all
  252 scalar-field scenarios, the only ones that shoot.
- **A/B, the proxy's run configuration** (Γ12 = 0.26): 61.5 s on the branch, 62.4 s on
  master, with one shooting build. Ω_Λ drops by the daughters' 3.0e-5, and the C_l move by
  rms 1.5e-5 (TT) and 2.9e-5 (EE).
