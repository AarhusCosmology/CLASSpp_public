# DNCDM_DR budget closure + uniform `GetOmega0()` semantics — Design

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:writing-plans to turn this into a
> task-by-task implementation plan, then superpowers:subagent-driven-development /
> superpowers:executing-plans to implement it.

**Date:** 2026-06-01
**Branch (to create):** `dncdm-dr-closure` (off `master`; #280 background-module-cleanup is merged)
**Status:** design approved by user, ready to plan

---

## 1. Problem

A `DNCDM_DR` (decaying non-cold dark matter → dark radiation) model does **not** close the
flatness budget: `Ω_total(z=0) ≈ 1 + Ω_dr`, i.e. ~1% off for typical parameters. The decay
radiation produced during integration is never reserved by the closure species (Λ).

Root cause, established by tracing the code:

- Flatness closure is computed once at input time, `input_module.cpp:244`:
  `Ω_closure = 1 − Ω_k − Σ_species GetOmega0()` (Pass-1 sum at line 229), **before** background
  integration.
- `DCDM_DR` closes because `DCDMSpecies::GetOmega0()` returns `Omega0_dcdmdr_` — the **combined**
  dcdm+dr today density (a value known at input time) — so the closure reserves the whole sector,
  and shooting then drives the integrated `(ρ_dcdm+ρ_dr)` to that same number.
- `DNCDM_DR` inherits `CompositeSpecies::GetOmega0()` = Σ children. The DR child returns 0 by
  design (`dark_radiation_species.h`: "decay product, starts at zero"), and in combined/initial
  modes the DNCDM matter child never gets an `Omega0_` either, so the composite reports
  matter-only (or 0). The decay radiation is structurally absent from the reservation.

Two secondary defects surfaced while tracing:

- **`dncdm().GetOmega0()` returns 0** in combined/initial modes (the matter child only sets `deg`,
  never `Omega0_`). This silently drops the decaying species' matter from `fnu`
  (`GetOmega0NcdmTot()` → HMcode/Halofit) **and** from the verbose budget-print neutrino line.
- The verbose budget print dropped the DNCDM decay-DR line entirely. This is **already fixed as an
  uncommitted working-tree edit from this session** in `background_output_budget()`
  (`background_module.cpp` ~1427): guard on `GetNDecayDr() > 0` instead of `count("DCDM_DR")`, since
  `Omega0_dr_` already aggregates DCDM_DR + all DNCDM_DR. This edit travels with this PR (carry it
  onto the `dncdm-dr-closure` branch).

## 2. Principle

`GetOmega0()` means **the species' total energy-density fraction today**, Ω = ρ(a=1)/ρ_crit,0 —
uniformly, composite or not. For a decaying-sector composite that total is the combined matter+dr
("`Omega_dncdmdr`"). This is the single quantity the closure reserves; at shooting convergence it
equals the integrated combined density, so `Λ = 1 − Ω_k − Σ GetOmega0()` makes the universe flat
(to shooting tolerance). It is precisely the existing DCDM pattern, generalised.

**Why no `ClosureOmega0()` and no `fnu` rework is needed:** the only matter-only consumer,
`BackgroundModule::GetOmega0NcdmTot()` (→ `fnu`), loops `GetNcdmSpecies()`, which returns each
composite's `dncdm()` **matter child**, never the composite. So a composite returning "combined"
cannot pollute `fnu`. The matter-part machinery (`MatterRho`/`MatterRhoDelta`) already serves the
*perturbation* matter tally (`delta_rho_m`) and is untouched here.

## 3. Scope

In scope (one PR):

1. Uniform `GetOmega0()` = total-today-fraction; `DNCDM_DR_Species` override returns the combined.
2. Fix `dncdm()` matter child `Omega0_` backfill.
3. Closure: rely on the composite `GetOmega0()` (line 229); delete the dead `omega0_ncdm_tot`
   tally (line 221/231 — declared, accumulated, never read).
4. Shooting: **combined** mode closes via (1) with no shooting change; **initial** mode
   (`Omega_ini`/`omega_ini`/`Neff_ini`) gains the DCDM-style `Omega_dncdmdr` fixed-point shoot.
5. Reject bare `dncdm*.Omega` (today-matter) on a decaying (`ncdm_decay_dr`) species.

Explicitly **out of scope** (YAGNI, per user): an `N_eff_today` / final-DR shooting mode for either
DNCDM or DCDM. The existing `Neff_ini` (initial-abundance ΔN_eff) is untouched.

## 4. Components

### 4.1 `GetOmega0()` semantics

- **`DNCDM_DR_Species::GetOmega0()` override** (new): returns the sector combined today fraction =
  the closure-reserved `Omega_dncdmdr`. Sourced from the value the composite carries (combined mode:
  the user input; initial mode: the shooting-pinned `Omega_dncdmdr` written into the file content
  each iteration). This mirrors `DCDMSpecies::GetOmega0() → Omega0_dcdmdr_`. The composite **surfaces
  the pinned value** rather than summing children, because pre-integration (Pass-1 closure) the
  children can't yet report the emergent DR.
- The matter child (`DNCDMSpecies::GetOmega0()`) keeps the base meaning: its own today matter
  fraction (see 4.2). The DR child keeps returning its own today fraction (0 pre-integration; the
  composite does not depend on it for the reservation).

### 4.2 `dncdm()` matter-child `Omega0_` backfill

Mirror how `Omega0_dcdm_` is set in `background_solve_evolver` (`background_module.cpp` ~925–934):
after integration, set each DNCDM matter child's `Omega0_` from the integrated matter density at
a=1, `ρ_dncdm,matter(today)/H0²`. Use the existing matter-vs-radiation split convention for
`EnergyType::Other` (the background `accumulate` lambda buckets `ρ − 3p` into matter,
`background_module.cpp:375-378`); the matter child's today fraction is that matter part.

Effect: `dncdm().GetOmega0()` reports the real decaying-species matter today, restoring it to
`fnu`/`GetOmega0NcdmTot()` and to the budget-print neutrino line. Behaviour-identical for stable
NCDM (their `Omega0_` is already set at construction).

### 4.3 Closure (input Pass-1)

- No change to the closure formula; it already sums `GetOmega0()` (line 229), which now includes the
  DNCDM_DR combined.
- **Delete** the dead `omega0_ncdm_tot` local (`input_module.cpp:221` decl, `:231` accumulate) — it
  is never read (leftover from the removed `pba->Omega0_ncdm_tot`). Keep the `n_ncdm`/`n_dncdm`
  counters (still used at line 274).

### 4.4 Shooting

Per-flavor taxonomy (mirrors `DCDM_DR_Species`):

| Mode | Input key(s) | abundance knob `deg` | `Omega_dncdmdr` (closure reserve) | shooting unknown |
|------|--------------|----------------------|-----------------------------------|------------------|
| **Combined** | `Omega_dncdmdr` / `omega_dncdmdr` | shot → hit combined target | given | `<flavor>.deg` |
| **Initial** | `Omega_ini` / `omega_ini` / `Neff_ini` | derived from initial (not shot) | fixed-point shot → integrated combined | `<flavor>.Omega_dncdmdr` |

- **Combined mode** is already fully wired (`GetShootingTargets` → `{<flavor>.Omega_dncdmdr,
  <flavor>.deg, Omega_dncdmdr_pending}`; residual `(ρ_dncdm+ρ_dr)/H0² − target`). It begins to
  *close* purely because of 4.1. No shooting-code change.
- **Initial mode** (currently sets `deg` directly via `ApplyDncdmInitialClosure` →
  `SetDeg_from_Omega_ini` and never shoots) gains:
  - `GetShootingTargets()` returns a target with `unknown_param = <flavor>.Omega_dncdmdr`,
    `target_name` marking the fixed-point branch.
  - `ComputeShootingResidual()` fixed-point branch: `−(ρ_dncdm+ρ_dr)/H0² + GetOmega0()` (drives the
    integrated combined to the reserved `Omega_dncdmdr`), parallel to `DCDM_DR` residual line 370.
  - `ComputeShootingGuess()` provides a finite positive guess + Jacobian seed for `Omega_dncdmdr`
    (a decay-factor-scaled estimate from the given initial, parallel to `DCDM_DR` guess lines
    340–351). Exact heuristic is an implementation detail; it only needs to seed Newton.
  - Parser: allow `<flavor>.Omega_dncdmdr` to **co-occur** with `Omega_ini`/`omega_ini`/`Neff_ini`
    **during shooting iterations** (the shooter writes it), exactly as `deg` is already allowed to
    co-occur with `Omega_dncdmdr` today (`dncdm_species.cpp:161`). Outside a shooting build the
    mutual exclusion stays.
  - The composite reads the shot `Omega_dncdmdr` so `GetOmega0()` (4.1) returns it at Pass-1.

### 4.5 Reject bare `.Omega` on a decaying species

In `DNCDMSpecies` construction, if a `ncdm_decay_dr` instance supplies `Omega`/`omega` (today
matter), throw a clear `invalid_argument` directing the user to `Omega_ini`/`Neff_ini` (pin the
sector by initial abundance) or `Omega_dncdmdr` (pin by combined today). Today `.Omega` calls
`SetOmega0` and leaves the budget open — a silent non-flat universe.

## 5. Data flow (one shooting build, initial mode)

```
DoShooting collects {<flavor>.Omega_dncdmdr unknown, guess} from the discovery module
 └─ each Newton iteration: ShootingResidual writes <flavor>.Omega_dncdmdr into fc
     └─ fresh Cosmology → input_read_parameters
         ├─ DNCDMSpecies reads Omega_ini (→ deg) AND Omega_dncdmdr (→ composite reserve)
         ├─ Pass-1 closure: Λ = 1 − Ω_k − Σ GetOmega0()   (reserves the combined)  [line 229]
         └─ background integration → bg_today
     └─ residual = −(ρ_dncdm+ρ_dr)/H0² + GetOmega0()      (fixed-point on Omega_dncdmdr)
 └─ converged: Omega_dncdmdr == integrated combined; closure consistent ⇒ Ω_total = 1 (to tol)
```

## 6. Testing / verification

- **Budget closes:** `dncdm_dr.ini` (rewritten to a supported mode — `Omega_ini` and/or
  `Omega_dncdmdr`) prints `TOTAL Ω ≈ 1` at `background_verbose=2`, in combined and initial modes.
  Add a regression scenario for each mode.
- **Combined ↔ initial agreement:** a combined run and an initial run tuned to the same physical
  sector agree on `Ω_total` and on the decay-DR / matter split (within shooting tol).
- **`fnu` restored:** a DNCDM_DR model with `output=mPk` and HMcode/Halofit has nonzero
  `GetOmega0NcdmTot()` matching the integrated dncdm matter today (previously 0).
- **Bare `.Omega` rejected:** a `ncdm_decay_dr` instance with `.Omega` errors with the guidance
  message.
- **No regressions elsewhere:** full scenario grid (`python/test_class.py -m test_scenario` vs the
  `classyref` build) + the 13-scenario `compare_tol.py` drift harness
  (`/tmp/bgcleanup/driftcheck.sh` pattern) — stable-NCDM, DCDM_DR, idm/idr unchanged. Stable-NCDM
  `fnu` must be bit-stable (the matter-children path is behaviour-identical for them).
- **Shooting convergence / tolerance:** confirm `fzero_Newton` converges for both modes; per the
  user, a small residual flatness violation at the shooting tolerance is acceptable.

## 7. Files (anticipated)

- `species/dncdm_dr_species.{h,cpp}` — `GetOmega0()` override; initial-mode `GetShootingTargets` /
  `ComputeShootingResidual` (fixed-point branch) / `ComputeShootingGuess`.
- `species/dncdm_species.{h,cpp}` — reject bare `.Omega`; relax `Omega_dncdmdr`↔initial co-occurrence
  during shooting; carry the shot `Omega_dncdmdr` for the composite; expose the combined for the
  override.
- `source/background_module.cpp` — backfill `dncdm()` `Omega0_` post-integration (near the existing
  `Omega0_dcdm_`/`Omega0_dr_` block).
- `source/input_module.cpp` — delete the dead `omega0_ncdm_tot` tally.
- `test/scenarios/` — combined + initial DNCDM_DR closing scenarios; update `dncdm_dr.ini`.

No new source files ⇒ no build-system list changes.

## 8. Non-goals / risks

- **Non-goal:** `N_eff_today`/final-DR mode (dropped). `Neff_ini` semantics unchanged.
- **Risk:** the initial-mode guess heuristic must seed Newton well enough to converge across
  Gamma/abundance ranges; reuse `DCDM_DR`'s decay-factor scaling as the template.
- **Risk:** `fnu` for stable NCDM must stay bit-identical — it must, since `GetOmega0NcdmTot()`
  reads matter children whose `Omega0_` are set at construction exactly as before. Verify explicitly.
