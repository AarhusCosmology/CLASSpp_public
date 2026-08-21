# Upgrading the bundled HyRec to HYREC-2

**Date:** 2026-08-21
**Issues:** #396 (HyRec is less accurate than RECFAST, Δχ² ≈ +25), #369 (HyRec reconstructs H(z) from a CPL pair)
**Branch:** `hyrec2-upgrade`

## Problem

`recombination = hyrec` is currently the *less* accurate option. Measured against Planck 2018
`plik_lite` TTTEEE over 144 cosmologies it costs Δχ² ≈ +25 relative to the RECFAST default, with a
mean TT shift of +0.85%. HyRec against RECFAST should be a ~0.1% effect.

Two defects stack, and reading the code turned up a sharper statement of the first than #396 has.

### A. The bundled code is HyRec Nov2011, compiled in its RECFAST mode

`hyrec/hyrec.h` declares `HYREC_VERSION "Nov2011"`. Worse, `hyrec/history.h:20` says

```c
#define MODEL RECFAST     /* default setting: FULL */
```

so `recombination = hyrec` today runs HyRec's *effective three-level atom with the F = 1.14 fudge
factor* — not its radiative-transfer model. Upstream `class_public` v2.9.4 carries the identical
line, so this is inherited from CLASS v2 rather than introduced here, but the consequence is the
same: the "more accurate" option is a cruder RECFAST (no v1.5 fudge functions) wrapped in HyRec's
helium treatment and its own H(z). That is why it lands ~1% away from our RECFAST through
recombination and into the residual ionization tail, which is what sets the C_ℓ damping.

HYREC-2 defaults to `SWIFT`, a fit to the full radiative-transfer calculation.

### B. HyRec rebuilds H(z) itself, and double-counts the massive neutrino

`thermodynamics_module.cpp` fills `param.omh2` with `GetOmega0NcdmTot()`, which HyRec then
redshifts as matter, while `param.Nnueff = Neff_` already counts the same neutrino as radiation.
With one 0.12 eV neutrino, HyRec's internal H is +0.29% off ours at z = 1100; without a massive
neutrino it is exact. This is #369's mechanism reached through a different term than the CPL one
that issue names.

The same reconstruction is why `FluidSpecies::HyrecCplApproximation` exists: axion-EDE fluids have
no faithful `(w0, wa)` representation and are rejected outright, and scalar-field dark energy is
not a `Fluid` at all, so its density is silently ignored.

## Decisions

| Decision | Choice |
|---|---|
| Where the library lives | Vendored in-tree under `hyrec/`, as today. Not a submodule. |
| Integration architecture | class_public v3's: HyRec supplies **derivatives**, CLASS integrates. |
| Recombination model | `SWIFT` (HYREC-2's default; `FULL` is unreachable — see below). |
| H(z) | CLASS's background, at every step. HyRec never computes H. |
| Energy deposition | Port Galli et al. 2013 (arXiv:1306.0563) properly, for both models. |

### Why vendored, not a submodule

The build already compiles `hyrec/*.c` directly and installs `hyrec/*.dat` into the wheel. A
submodule would complicate sdists and offline builds for a library of nine small files that changes
roughly once a decade. Upstream `nanoomlee/HYREC-2` and `class_public` v3.3.4's
`external/HyRec2020/` are byte-identical except for whitespace and one `malloc` fix inside a
CAMB-only wrapper we do not call; all four `.dat` files match exactly. So we vendor upstream
verbatim and record the commit in `hyrec/PROVENANCE.md`, making the next upgrade a `diff` rather
than an archaeology exercise.

### Why derivatives rather than a one-shot history

HYREC-2 offers both: `rec_build_history(data, model, hubble_array)` builds a whole history with an
optionally injected H(z) array, or the individual derivative routines can be called from a host
integrator. class_public v3 uses the latter, and CLASS++ is already shaped for it.

`thermodynamics_recombination_with_recfast` (as refactored by #362) evolves exactly
`(x_H, x_He, T_mat)` through three phases — `analytic` (both Saha), `helium` (He ODE, H on its Saha
branch), `full` (both ODEs) — driven by the shared evolver. class_public's six-phase `ap_*` scheme
is the same idea at finer granularity, and it swaps *only* the hydrogen and helium derivatives
between RECFAST and HyRec. Everything else — the T_mat equation, the Saha phases, the evolver, the
dense-output sampling — is shared.

The one-shot route would have worked too and been a smaller diff, but it keeps HyRec's own
integrator, its own T_mat, and its own grid, so the two recombination options would agree on less.
The derivative route makes `recfast` and `hyrec` differ in exactly one place: the atomic physics.

## Design

### 1. The vendored library

`hyrec/` keeps its name; its contents are replaced with HYREC-2 verbatim.

- **In:** `energy_injection.{c,h}`, `helium.{c,h}`, `history.{c,h}`, `hydrogen.{c,h}`,
  `hyrectools.{c,h}`; data files `Alpha_inf.dat`, `R_inf.dat`, `two_photon_tables.dat`, and the new
  `fit_swift.dat`.
- **Out:** `hyrec.c` (standalone `main`), `hyrec.h`, `input.dat`, the old `Makefile`.
- **New:** `hyrec/PROVENANCE.md` — upstream URL, commit hash, and the statement that the files are
  unmodified.

No edits to the vendored sources. The two things that would otherwise force edits are handled on
our side:

- **C++ linkage.** The `extern "C"` block lives in our wrapper header, which is the only place that
  includes HYREC-2's headers.
- **`exit(1)`.** HYREC-2 calls it in four places, none of them reachable from the derivative
  architecture: `hyrec_xe`/`hyrec_Tm` (we interpolate nothing — CLASS owns the history),
  `rec_get_cosmoparam` (reads the standalone driver's stdin), the CAMB Fortran wrappers, and
  `interp_Dfnu` (FULL-mode radiative transfer only). A `test-hyrec` case pins that none of them can
  be entered.

`CMakeLists.txt` gains `hyrec/energy_injection.c` in `CLASS_HYREC_FILES`. The wheel's
`hyrec/*.dat` glob picks up `fit_swift.dat` with no change.

**Input surface.** HYREC-2 takes a *directory* and appends its own filenames, so the three
precision strings `Alpha_inf hyrec file`, `R_inf hyrec file` and `two_photon_tables hyrec file`
collapse into one `hyrec_path` (default `/hyrec/`). The parser rejects the three retired names with
a message naming the replacement, rather than ignoring them silently.

### 2. The recombination-model seam

A `RecombinationModel` interface with two members:

```cpp
class RecombinationModel {
 public:
  virtual ~RecombinationModel() = default;
  virtual double dx_H_dz (const RecombinationState& s, const EnergyDeposition& dep) const = 0;
  virtual double dx_He_dz(const RecombinationState& s, const EnergyDeposition& dep) const = 0;
};
```

where `RecombinationState` carries `{z, x_H, x_He, x, n_H, H, T_mat, T_rad}` — everything both
models need, gathered once per derivative evaluation.

`RecfastModel` holds the existing TLA, the v1.4 helium corrections and the v1.5 fudge functions,
moved across unchanged. `HyrecModel` owns the `HYREC_DATA` lifetime (RAII: allocate on
construction, `hyrec_free` in the destructor). The chosen model is constructed once by
`thermodynamics_recombination_*` and lives for the whole recombination integration — the atomic and
SWIFT-fit tables are read from disk exactly once per CLASS run, as today. It forwards:

- `dx_H_dz  = -1/(1+z) · rec_dxHIIdlna(data, model, x, x_H, n_H·1e-6, H, T_mat·kBoltz, T_rad·kBoltz, 0, z)`
- `dx_He_dz = -1/(1+z) · rec_helium_dxHeIIdlna(data, z, 1-x_H, x_He·fHe, H) / fHe`

with HyRec's own guards: `model = PEEBLES` when `T_rad·kBoltz ≤ TR_MIN` or `T_mat/T_rad ≤
T_RATIO_MIN`, and `dx_He_dz = 0` once `x_He·fHe < XHEII_MIN` (1e-6). Note the definition mismatch:
CLASS's `x_He` is HeII over total helium, HyRec's `xHeII` is relative to hydrogen, hence the `fHe`
factors.

`thermodynamics_derivs_with_recfast_member` calls the model rather than branching on
`pth->recombination`, so no algorithm-name test appears in module code.
`thermodynamics_recombination_with_hyrec` is **deleted** — about 280 lines, including the
hand-rolled contiguous rate-table buffer and the `fscanf` loops, all of which HYREC-2 does itself.

**Cosmology fields.** `HyrecModel` fills `REC_COSMOPARAMS` once at construction: `T0`, `obh2`,
`ocbh2`, `okh2`, `YHe`, `fHe`, `Neff`, `nH0` (**cm⁻³**, not m⁻³ as in the old struct), `fsR = meR =
1` (CLASS++ has no varying constants), and `dlna` defensively even though nothing on this path
reads it. `hyrec_allocate(data, z_initial, 0.)` reads the atomic and SWIFT-fit tables; its output
arrays go unused, and the radiation-field tables are only allocated in `FULL` mode.

`ocbh2` feeds SWIFT's fiducial-difference correction and must be **non-free-streaming** matter —
class_public passes `Omega0_nfsm`, i.e. excluding massive neutrinos. Putting `GetOmega0NcdmTot()`
here would repeat defect B in a different slot. CLASS++ has no such aggregate today; it should be
added as a species-level query rather than by extending the `all_species_.count("CDM")` string
branching, which is already a known smell in this module.

**Errors.** HYREC-2 reports through `data->error` and `data->error_message`. Every call site checks
and raises via `class_test` (`runtime_error` → `CosmoComputationError`, so the point is rejected
rather than the run aborted) — a numeric failure, not a structural one, per the #382 severity
convention.

**Why `SWIFT` and not `FULL`.** `FULL` needs the radiative-transfer history indexed on HyRec's own
uniform `lna` grid, which a callback from an adaptive evolver cannot supply — steps get rejected
and retried, and there is no monotone step index to hand it. class_public carries the same
restriction (`class_test(MODEL == FULL, ...)`). `SWIFT` is a fit *to* `FULL`, verified stateless in
`iz`, so this costs accuracy only at the level the HYREC-2 authors accepted. `PEEBLES` remains
reachable as HyRec's own low-temperature fallback.

### 3. Energy deposition

HyRec wants two *separate* rates — `inj_params->ion` and `inj_params->exclya`, each per hydrogen
atom per second — and applies the `(1−C)`-style branching itself inside the multilevel atom. CLASS++
today lumps both into one `chi_ion_H` and splits them with RECFAST's Peebles `C`:

```
energy_rate · chi_ion_H / n · (1/L_H_ion + (1−C)/L_H_alpha)
```

`thermodynamics_energy_injection` stays the single source of the injected rate, including the halo
term and the non-on-the-spot deposition integral. HYREC-2's own `energy_injection.c` is unusable
here: `update_dEdtdV_dep` is stateful in `lna` and our evolver takes rejected steps.

What is missing is the channel split. CLASS++ carries two analytic fits, commented as Slatyer et al.
2013; class_public carries the same two expressions under the name `chi_Galli_analytic` and sets
`chi_lya = 0`. They disagree with the actual Galli et al. 2013 table at low ionization —
`chi_heat → 0` in the fit versus 0.152 in the table, `chi_ion_H` 0.369 versus 0.351 — so the fits'
provenance is not solid enough to build on.

So: a `DepositionModel` supplying `{heat, ion_H, ion_He, lya, lowE}` as functions of `x`.

- `Galli2013` — the 18-row, 6-column table from arXiv:1306.0563, embedded as a `constexpr` array
  (1.3 KB; no new data file, no new path parameter) and splined in `x`. **Default.**
- `GalliAnalytic` — the two fits shipped today, with `lya = 0`, retained so existing results remain
  reproducible.

Both recombination models consume it: RECFAST as `ion_H/E_ion + lya·(1−C)/E_lya`, HyRec as
`ion = ion_H/(n_H·E_ion)` and `exclya = lya/(n_H·E_lya)`. Neither model fudges a channel it does
not have.

This changes `recombination = recfast` results whenever annihilation or decay is switched on. The
shift is measured and reported in PR 1, but it is not a gate: `Galli2013` is the published table
and stays the default on its own merits, with `GalliAnalytic` available for reproducing
pre-change results.

The `annihilation_f_halo > 0 ⇒ recombination = hyrec` requirement is unchanged: it exists because
halo annihilation pushes values outside the RECFAST fits' validity, which this work does not touch.

### 4. What #369 costs: nothing

The derivative is handed `Hz` from `pvecback` at every step, exactly as RECFAST already is.
`rec_HubbleRate` is never called, so:

- Defect B disappears — there is no internal H(z) to double-count the neutrino in.
- #369 closes — no CPL reconstruction at all. `FluidSpecies::HyrecCplApproximation`,
  `AxionEDEFluid`'s override, and the `class_test` that rejects axion-EDE under HyRec are deleted,
  along with their test coverage in `fluid_test.cpp` and `axion_ede_fluid_test.cpp`.
- Scalar-field dark energy, silently ignored by HyRec today, is simply included.

## Sequencing

Two PRs, because the deposition change moves RECFAST results and should not hide inside a HyRec
diff.

**PR 1 — deposition-channel model.** `DepositionModel`, the embedded Galli table, the `chi_type`
input, and the RECFAST derivs rewired to consume it. Touches nothing else. Its effect on the
annihilation and decay scenarios is measured in isolation and reported.

**PR 2 — HYREC-2.** The vendored library, the `RecombinationModel` seam, `HyrecModel`, the deletion
of `thermodynamics_recombination_with_hyrec` and the CPL machinery, the `hyrec_path` input change,
and the tests below. Closes #396 and #369.

## Validation

**Parity against class_public v3.3.4.** Both sides would be HYREC-2/SWIFT with CLASS's H(z) and a
CLASS-integrated T_mat, so `x_e(z)` should agree to a few × 10⁻⁴ across recombination and the
residual tail. This checks *our seam*, not the physics: where the two disagree we investigate
rather than assume class_public is right.

**The issue's own metric** (`benchmarks/precision/s21_hyrec.py`, on the campaign branch). CLASS++'s
HyRec against CLASS++'s RECFAST should land near where class_public's does:

| | today | target |
|---|---|---|
| mean ΔTT | +0.880% | ≈ +0.04% |
| max ΔTT | 2.21% | ≈ 0.10% |
| S vs own RECFAST | 22.8 | ≈ 1.5 |
| Δχ² vs `plik_lite`, 144 cosmologies | +25 | ≲ 1 |

**Ionization history by z**, against class_public HyRec, same cosmology. Today's divergence is
+0.45% at z = 1098, −1.14% at z = 800, −1.31% at z = 200; all three should drop below ~0.05%.

**H(z).** Assert in a unit test that `HyrecModel` receives the background's `H` — with one 0.12 eV
neutrino the old path was +0.29% off at z = 1100, and the new path has no independent H to compare.

**New scenarios** in `test/scenarios/`: `recombination = hyrec` with a massive neutrino (the
combination #396 notes is untested — `grep -rl "hyrec" test/` currently returns nothing), and
`recombination = hyrec` with `Omega_scf` (previously ignored) or an axion-EDE fluid (previously
rejected).

**Unit target `test-hyrec`**, registered in **both** `CMakeLists.txt` and the Makefile
`TEST_TARGETS` — CI only builds targets named in the latter, and forgetting it is a recurring trap
here. Covers: the vendored tables load from `hyrec_path`; `rec_dxHIIdlna` returns finite values
across the SWIFT range; the `PEEBLES` fallback triggers below `TR_MIN`; the helium cutoff triggers
below `XHEII_MIN`; a forced HYREC-2 error surfaces as an exception rather than `exit()`.

**Build hygiene.** The project compiles with `-ffast-math`; confirm HYREC-2 has no non-finite test
that folds away under it, and that the SWIFT interpolation is unaffected.

**Regression.** `recombination = recfast` is untouched by PR 2, so `COMPARE_OUTPUT_REF` should stay
valid; verify rather than assume. Measure thermodynamics wall time before and after — HyRec's rate
interpolation now runs at every evolver step instead of once per history.

## Out of scope

- Porting class_public v3's `injection` module (PBH accretion, `f_eff` files, `chi` from external
  files). Only the deposition-channel split is needed here.
- Varying fundamental constants (`fsR`, `meR`). HYREC-2 supports them; CLASS++ has no such input,
  so both are pinned to 1.
- Making the recombination model runtime-selectable beyond `recfast` / `hyrec`.
- #399 (halofit `nonlinear_min_k_max`) and #398 (tensor TCA), which are unrelated.
