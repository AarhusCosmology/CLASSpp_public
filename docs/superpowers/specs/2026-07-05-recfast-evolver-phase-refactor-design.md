# RECFAST evolver path: collapse phase modes onto one state-based RHS

- **Date:** 2026-07-05
- **Branch:** `recfast-evolver-thermodynamics`
- **PR:** #362 (`[codex] Use evolver callbacks for RECFAST thermodynamics`)
- **Status:** design approved, ready for implementation plan

## Goal

Simplify the structure of the RECFAST evolver integration introduced in PR #362,
without changing the physics. The phase-driven evolver approach (identify the two
recombination transition redshifts up front, fill analytic regimes directly, drive
the smooth intervals with the shared evolver + dense output) is kept. What changes
is the *derivative and callback machinery* wrapped around it.

## Background

PR #362 replaced the last RECFAST `generic_integrator` consumer with the shared
`evolver_*` interface. The physics is faithful and the approach is sound. But the
derivative layer grew two orthogonal mode enums —

```cpp
enum class RecfastHydrogenMode { legacy_trigger, frozen, evolved };
enum class RecfastHeliumMode   { legacy_trigger, evolved };
```

— threaded through a core RHS `thermodynamics_recfast_derivs_member(..., hydrogen_mode, helium_mode)`,
plus three static derivs shims, three member adapters, and three output callbacks.

### Key insight (corrected): the regime switches are double-encoded

Every regime switch in master's `thermodynamics_derivs_with_recfast_member` is a pure
function of **state** (`x_H`, `x_He`, `Tmat`), decided *inside* the RHS:

- Hydrogen freeze: `if (x_H > recfast_x_H0_trigger) dy[0] = 0`
- Helium 1.4 corrections: `Heflag` active only for `5e-9 <= x_He <= recfast_x_He0_trigger2`
- Temperature: tight-coupling branch gated by `timeTh < H_frac * timeH`

PR #362 added *redshift* phase boundaries on the outside (root-found
`z_helium_ode_start`, `z_hydrogen_ode_start`) but left these state re-tests on the
inside — so the hydrogen-Saha and helium Saha↔ODE transitions are now **encoded
twice**: once as a redshift boundary in the driver, once as a state threshold in the
RHS. The two copies must be kept in sync.

An earlier version of this spec proposed reverting the RHS to master verbatim,
believing the state tests were inert given the driver's seeding. **That is false for
the helium corrections.** The helium ODE is seeded at `x_He = recfast_x_He0_trigger =
recfast_x_He0_trigger2 = 0.995` *exactly*, and `x_He` drifts a hair above 0.995 during
integration, so master's `x_He > trigger2 → Heflag = 0` gate chatters right at the
seed. That perturbs the trajectory at the ~1e-5 level, which the **default** ndf15
evolver's finite-difference Jacobian amplifies to ~3e-3 in `c_b^2` (explicit rkdp45
only sees ~1e-4). Verified by bisection: this single gate is the entire divergence;
the hydrogen-freeze re-test happens to stay inert only because hydrogen recombination
is monotone at its seed.

## Target design

### 1. Own the regime switches in the phase, not in RHS state (the perturbations pattern)

Follow the approximation-scheme pattern the perturbations module uses: the phase owns
the approximation state, set once at the boundary; the RHS *reads* it instead of
re-deriving it from `y`. This removes the double-encoding and the seed flicker.

- Delete `RecfastHydrogenMode`, `RecfastHeliumMode`, and the fat
  `thermodynamics_recfast_derivs_member(..., hydrogen_mode, helium_mode)`; keep a single
  `thermodynamics_derivs_with_recfast_member(double z, double* y, double* dy, void*)`.
- Near the top of that RHS, derive two flags from the phase and drop the two
  double-encoded state re-tests:
  ```cpp
  const RecfastPhase phase      = ptpaw->recfast_phase;
  const bool hydrogen_frozen    = (phase != RecfastPhase::full);
  const bool helium_corrections = (phase != RecfastPhase::analytic);
  ...
  if (hydrogen_frozen) dy[0] = 0.; else { /* Peebles ... */ }
  ...
  int Heflag = 0;
  if (helium_corrections && (x_He >= 5.e-9)) Heflag = ppr->recfast_Heswitch;
  ```
- The `x_He >= 5e-9` floor and the `Tmat` tight-coupling test stay dynamic state
  tests: neither is double-encoded (no redshift boundary shadows them), the `Tmat`
  switch is a single monotone crossing that is genuinely state-dependent (it depends
  on `x_e(z)`, not `z` alone), and both sit well away from any seed.
- Unrelated tidy in the same function: express the `He_Boltz` overflow guard as
  `exp(std::min(680., preco->Bfact / Tmat))` instead of an if/else. Clamping the
  *argument* is bit-identical to the original `exp(680.)` saturation and is
  `-ffast-math` safe (it never forms the `+inf` that `-ffinite-math-only` assumes away
  — clamping the *result*, `std::min(1e295, exp(...))`, would not be).

This reproduces PR #362's current output **bit-for-bit** (helium/full keep the ungated
Heflag; analytic keeps it off because `x_He >= trigger2` there): thermodynamics A/B
`max|rel| = 0` for both evolvers. It is a pure structural change, not a physics change.

### 2. Carry the phase in the workspace

```cpp
enum class RecfastPhase { analytic, helium, full };

struct thermodynamics_parameters_and_workspace {
  ...
  int recfast_output_index_offset = 0;
  RecfastPhase recfast_phase = RecfastPhase::full;
};
```

The analytic regime does not run through the evolver (the driver calls the RHS
directly for the `cb2`/`dTb` column), but it sets `recfast_phase = analytic` before
that call so the RHS applies the right approximation flags — hence three values.

### 3. One static derivs shim, switching on phase

Replace `thermodynamics_derivs_with_recfast_minus_z`, `..._helium`, and their
member adapters (`..._full_member`, `..._helium_member`) with a single shim. It
handles the `minus_z` sign convention and the 2-var vs 3-var packing:

```cpp
static void thermodynamics_recfast_derivs(double minus_z, double* y, double* dy, void* params) {
  auto ws = static_cast<thermodynamics_parameters_and_workspace*>(params);
  auto* module = ws->thermodynamics_module;
  const double z = -minus_z;

  if (ws->recfast_phase == RecfastPhase::helium) {
    // y = {x_He, Tmat}; reconstruct frozen x_H from Saha
    double y_full[3]  = { module->..._hydrogen_saha_xH(ws->preco, z), y[0], y[1] };
    double dy_full[3];
    module->thermodynamics_derivs_with_recfast_member(z, y_full, dy_full, params);
    dy[0] = -dy_full[1];
    dy[1] = -dy_full[2];
  } else {
    // y = {x_H, x_He, Tmat}
    double dy_dz[_RECFAST_INTEG_SIZE_];
    module->thermodynamics_derivs_with_recfast_member(z, y, dy_dz, params);
    for (int i = 0; i < _RECFAST_INTEG_SIZE_; i++) dy[i] = -dy_dz[i];
  }
}
```

The static `thermodynamics_derivs_with_recfast(z, y, dy, params)` (member forwarder,
non-`minus_z`) is retained for the analytic-phase direct call at the `cb2` step.

### 4. One static output shim, switching on phase

Replace `thermodynamics_recfast_output_helium` / `..._output_full` /
`..._output_none` with a single shim:

```cpp
static void thermodynamics_recfast_output(double minus_z, double y[], double dy[],
                                          int index_x, void* params) {
  auto ws = ...; auto* module = ws->thermodynamics_module; auto* preco = ws->preco;
  const double z = -minus_z;
  const int sample_index = ws->recfast_output_index_offset + index_x;

  if (ws->recfast_phase == RecfastPhase::helium) {
    const double y_full[3] = { module->..._hydrogen_saha_xH(preco, z), y[0], y[1] };
    const double xe = module->..._xe_after_helium_ode(preco, z, y_full);
    module->..._store_row(preco, sample_index, z, xe, y[1], -dy[1]);
  } else {
    const double xe = module->..._xe_after_full_ode(preco, z, y);
    module->..._store_row(preco, sample_index, z, xe, y[2], -dy[2]);
  }
}
```

### 5. Driver changes

- Set `tpaw.recfast_phase = RecfastPhase::helium` before the helium evolver call and
  `RecfastPhase::full` before the full call; pass the single `thermodynamics_recfast_derivs`
  and `thermodynamics_recfast_output` to both.
- Remove the dead zero-sample helium path (`thermodynamics_recfast_output_none`,
  `helium_boundary_sample`, the `helium_output`/`helium_sampling`/`helium_sampling_count`
  selection). Guard with `class_test(helium_sample_count > 0, ...)`: at `recfast_Nz0 = 20000`
  over `z_initial = 1e4` (grid spacing 0.5) the helium phase always contains thousands of
  samples; a zero count would mean a nonsensically small `Nz0` (< ~10).
- The `used_in_output_*` arrays stay (the evolver requires them) but can be plain
  local `int[...]` rather than `std::vector` if convenient. Not load-bearing.

### Function inventory (before → after)

| Symbol | Fate |
|--------|------|
| `RecfastHydrogenMode`, `RecfastHeliumMode` | **delete** |
| `thermodynamics_recfast_derivs_member(..., modes)` | **delete** (body folds back into member RHS) |
| `thermodynamics_derivs_with_recfast_member` | keep; **revert body to master** |
| `thermodynamics_derivs_with_recfast` (static forwarder) | keep (analytic-phase `cb2` call) |
| `thermodynamics_derivs_with_recfast_minus_z`, `..._helium` (static) | **replace** with one `thermodynamics_recfast_derivs` |
| `thermodynamics_derivs_with_recfast_full_member`, `..._helium_member` | **delete** (logic moves into the shim) |
| `thermodynamics_recfast_output_full`, `..._helium` (static) | **replace** with one `thermodynamics_recfast_output` |
| `thermodynamics_recfast_output_none` | **delete** (dead path removed) |
| `thermodynamics_recfast_timescale` | keep |
| Saha helpers, `_store_row`, `_xe_after_helium_ode`, `_xe_after_full_ode` | keep unchanged |

## Intended vs unintended behavior change

- **None.** The refactor deletes the mode machinery and moves the two double-encoded
  regime switches onto the phase, choosing the flag semantics that reproduce PR #362's
  current behavior exactly (helium/full: Heflag on for `x_He >= 5e-9`, no upper gate;
  analytic: off). The `He_Boltz` limiter change is bit-identical. Verified: thermodynamics
  A/B `max|rel| = 0` for both ndf15 and rkdp45, and `compare_evolver.sh` unchanged.
- Note: this deliberately does **not** adopt master's `x_He <= x_He0_trigger2` Heflag
  gate — that gate chatters at the helium seed and is the ~3e-3 ndf15 divergence
  documented above. Keeping PR #362's ungated ODE-phase behavior is both smoother for
  the default evolver and keeps the refactor a true no-op.

## Tolerance study (`tol_thermo_integration`)

PR #362 tightened `tol_thermo_integration` from `1e-2` to `1e-6`. The old value was a
per-bin RK correction tolerance over tiny `Δz = 0.5` bins; the new one is a true adaptive
ODE tolerance spanning whole recombination phases, so the two are not comparable and the
value deserves an empirical justification rather than a guess. This work picks the value
by a convergence study instead of asserting it.

**Method.** With the refactor in place (so the RHS is master's), treat a very tight run as
truth and sweep candidates:

- Reference: `tol_thermo_integration = 1e-9`.
- Candidates: `1e-4, 1e-5, 1e-6, 1e-7, 1e-8`.
- For each, run `explanatory.ini` with thermodynamics output enabled and record, versus the
  reference on a common `z` grid: `max |Δx_e/x_e|` and `max |ΔTb/Tb|` over recombination
  (say `800 < z < 3000`), plus the thermodynamics wall time.
- Spot-check one downstream observable: `max |ΔC_ell^TT / C_ell^TT|` (handling TE-style
  zero crossings is unnecessary here since TT has no sign changes), for the loosest few
  candidates.

**Decision criterion.** Choose the loosest tol that both (a) keeps `x_e`/`Tb` deviation
from the reference `<= 1e-5` and `C_ell^TT` deviation `<= 1e-4`, and (b) keeps `c_b^2` — the
noisiest column, since it carries the numerical `dT_b/dz` derivative and feeds the baryon
perturbation equations — well converged. Runtime is not a tiebreaker here: thermodynamics is
computed once and the wall time is flat from 1e-6 to 1e-8, so where two tols both clear the
observable bar, prefer the one that also tightens `c_b^2`. Run the sweep for the default
evolver (`ndf15`) and confirm the chosen value also converges under `rkdp45`.

**Outcome (measured).** `1e-4` fails outright (unphysical `z_rec <= 500`); `1e-5` is marginal
(`x_e` ~5e-4); `1e-6` meets the `x_e`/`Tb`/TT bars but leaves `c_b^2` at ~6e-3; `1e-7` tightens
`x_e` ~10x and `c_b^2` ~14x (to ~4e-4) at no measurable runtime cost; `1e-8` is diminishing
returns. `C_ell^TT` is insensitive across the whole range (`< 2e-5` vs `1e-9`), and `rkdp45`
agrees with `ndf15` to ~1e-6 at `1e-7`. **Chosen: `1e-7`** — the `c_b^2` margin is the
deciding factor. `precision.h::tol_thermo_integration` set to `1e-7`; sweep table recorded in
the PR description.

## Non-goals

- No change to the phase-driven evolver approach, the 2-var/3-var split, or the analytic
  regime handling (all kept from PR #362).
- Not unifying helium and full into a single 3-variable system (explicitly deferred).

## Verification plan

Two changes ship here and are verified separately so each stays interpretable:

**A. Structural refactor — must be inert.** With `tol_thermo_integration` left at the branch
tip's `1e-6`, the refactor only deletes branches inert at default precision and reorganizes
callbacks, so output must be **unchanged versus the current branch tip** (PR #362 HEAD).

1. `make -j` — clean build.
2. A/B the thermodynamics table against the pre-refactor branch tip (same tol):
   - `./class explanatory.ini` and the recfast evolver `.ini`s used in the PR, with
     thermodynamics output columns enabled;
   - compare `x_e`, `Tb`, `c_b^2`, `dkappa/dtau` HEAD-before vs HEAD-after → expect
     agreement to machine precision (only inert code removed at defaults).
3. `test/scenarios/compare_evolver.sh ./class` → same result as PR #362 (the pre-existing
   `RUNFAIL(ndf15) type3_scf_veta` singular-matrix failure, reproduced on master).
4. Exercise both non-default evolvers (`evolver=ndf15` default and `evolver=rkdp45`) on a
   recfast run; confirm no regression.

**B. Tolerance change — a deliberate, converged output change.** Only after A is confirmed
inert, run the convergence sweep (see *Tolerance study*) and set the chosen value. Verify:

5. The chosen tol's `x_e`/`Tb` deviation from the `1e-9` reference meets the decision
   criterion, and record the sweep table.
6. Sanity check a CMB observable (TT) against master at ~0.1% per repo convention — not a
   gate (PR #362 already validated vs master), confirming the whole path is sound end to end.
