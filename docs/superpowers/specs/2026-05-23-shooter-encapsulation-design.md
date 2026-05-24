# Shooter encapsulation — design (PR D of #268 follow-ups)

## Goal

Move the species-specific "shooting" (root-finding for input parameters that can't be set directly)
out of `InputModule`'s enum-driven dispatch into per-species hooks, and lift the shooting
orchestration out of the `InputModule` constructor into a static `InputModule::DoShooting` invoked
lazily by `Cosmology`. After this PR:

- `InputModule`'s constructor is pure: read precisions + parameters, then `ConstructSpecies()`. **No
  shooting.** A species whose *target-form* input is set but whose direct *unknown* is absent builds
  itself from its own guess and records the target.
- `enum target_names`, `kTargetNamestrings_`, `kUnknownNamestrings_`,
  `input_auxillary_target_conditions`, and the per-target `switch` statements in
  `input_try_unknown_parameters` / `input_get_guess` are **deleted**. Each species owns its targets.
- The one cosmological target (`100*theta_s`, unknown `h`) stays module-level.

This is **not** behaviour-neutral for dncdm shooting (see "Latent bug" below); the dcdm / scf /
theta_s paths are behaviour-preserving and test-covered.

## Background: the current shooting machinery

`InputModule(fc)` runs `input_init()` → (a) read precisions, (b) **shoot**, (c)
`input_read_parameters()`; then `ConstructSpecies()`. Shooting therefore happens *before species
exist*, which is why it is enum/string-driven over the raw `FileContent`.

A single `enum target_names` (12 values) drives five concerns spread across `InputModule`:

| # | Concern | Today |
|---|---|---|
| 1 | detect requested targets in the fc | scan loop over `kTargetNamestrings_` (`input_init`) |
| 2 | which fc param each target varies | `kUnknownNamestrings_[index_target]` |
| 3 | is a detected target active? | `input_auxillary_target_conditions` switch |
| 4 | residual = computed − target | `input_try_unknown_parameters` switch (builds a fresh `Cosmology`, reads its species' `Rho`) |
| 5 | initial guess + Jacobian seed | `input_get_guess` switch (builds a fresh `InputModule`, uses species methods) |

The residual is evaluated by building a **whole fresh cosmology each iteration** (`new InputModule(fc)`
with the unknown set → `Cosmology` → compute background) and reading off the computed quantity. The
current `switch (target_name) { case Omega_dcdmdr: static_cast<DCDM_DR_Species&>(...) }` downcasts to
species types in module code — the anti-pattern this refactor removes.

Species owning targets: `DCDM_DR` ↔ {`Omega_dcdmdr`, `omega_dcdmdr`, `Omega_ini_dcdm`,
`omega_ini_dcdm`}; `ScalarField` ↔ {`Omega_scf`}. `theta_s` is cosmological (varies `h`, residual from
thermo `rs/ra`) and is owned by no species.

**The 6 dncdm enum targets are dead** (`Omega_dncdmdr`, `omega_dncdmdr`, `deg_ncdm_decay_dr`,
`Omega_ini_dncdm`, `Neff_ini_dncdm`, `omega_ini_dncdm`): the #260 dot-syntax migration moved decaying-
ncdm input to per-instance dot keys and **hard-rejects** every flat key
(`DNCDMSpecies` `RejectLegacyDecayDrKeys`, called from its `CreateAll`). The per-flavor amount is set
directly — `nu_i.deg`, `nu_i.Omega_ini`/`omega_ini`/`Neff_ini` → `Omega_ini_pending_` →
`SetDeg_from_Omega_ini` (an analytic closure at `a_ini`, **not** root-finding). So no dncdm root-find
shooting remains; the flat enum targets are unreachable and are simply deleted with the enum (Task 8).

## New: per-flavor dncdm today-density shooting

Since the flat dncdm targets are dead (above) and no per-flavor *today-density* fit exists yet, PR D
**adds** one: new per-instance dot inputs `nu_i.Omega_dncdmdr` / `nu_i.omega_dncdmdr` (= `/h²`) specify
the target today-density of that flavor's `dncdm + dr`, and the shooter varies that flavor's
`nu_i.deg` to hit it. For a decaying flavor the actual today-density (after decay) is not analytic in
`deg`, so this genuinely needs root-finding (unlike `nu_i.Omega_ini`, which is the analytic
`SetDeg_from_Omega_ini` closure). It is encapsulated in `DNCDM_DR`'s shooter hooks (mirroring DCDM_DR),
mutually exclusive with the other amount inputs (`deg`/`Omega`/`Omega_ini`/`omega_ini`/`Neff_ini`). No
presence issue — the flavor is already constructed from `nu_i.type = ncdm_decay_dr`. Residual/guess are
ported from the (corrected, A_s-removed) dncdm math in `input_try_unknown_parameters` /
`input_get_guess`. This was never exercised before, so it is verified by construction + convergence
check + author eyeball (no master reference).

(Aside on the old code: `kUnknownNamestrings_` even had a stray `"A_s"` — leftover from the retired
`sigma8 → A_s` shooting target — that misaligned every dncdm row, e.g. `Omega_dncdmdr → "A_s"`. Moot
now that those flat targets are deleted.)

## Architecture

### 0. Presence follows construction (enables instance-based detection)

Instance-based shooter detection requires the species to be *constructed* before it can report a
target. But "ini-form" targets (`Omega_ini_dcdm`/`omega_ini_dcdm`, `Omega_ini_dncdm`/`Neff_ini_dncdm`/
`omega_ini_dncdm`) give the *initial* density and leave the *today* density absent, while presence is
gated on the today density (`has_dcdm = Omega0_dcdmdr != 0`, and `CreateAll` builds only when
`has_dcdm`). The old enum-scan detected these before species existed; the instance-based design can't.

**Decision (per the project direction that `has_*` flags should eventually disappear — a species is
present iff constructed, and constructed iff *any* input says so):** a shooting composite is built when
**any** of its trigger inputs is present — the today-density target *or* an ini-form target. Concretely
`ConstructSpecies`'s presence determination sets `has_dcdm`/`has_dr` (resp. the dncdm presence) true
when an ini-form target is present in the fc, not only when the today density is set; the species then
self-constructs, guessing whichever density is the unknown (decision 1). `has_*` flags are still
written (they gate downstream physics) but the *build decision* now follows the union of inputs. This
is the minimal step toward presence=construction; fully removing `has_*` is a later PR.

### 1. `InputModule` constructor becomes pure (species guess on construction)

`input_init()` drops the shooting block (the `unknown_parameters_size` detection loop +
`FixUnknownParameters`). It becomes: read precisions → `input_read_parameters()`. Then
`ConstructSpecies()` runs as today.

Each shooting-capable species, in its `CreateAll`, reads **both** its direct unknown (e.g.
`Omega_ini_dcdm`) and its target-form input (e.g. `Omega_dcdmdr`):

- direct unknown present → build from it (today's behaviour); still record the target value if its
  target-form input is also present (needed for the residual during iterations).
- direct unknown absent but target-form present → call `ComputeShootingGuess` for its own guess,
  build from the guess, and record the `ShootingTarget` (`{target_name, unknown_param, target_value}`)
  as *needing shooting* (reported by `GetShootingTargets`).

The user supplies either the target or the direct unknown, never both, so "target-form present" ⟺
"this species needs shooting". A freshly-built module is therefore always valid (possibly carrying
guesses). The `file_content_.is_shooting` flag and its recursion-guard role go away.

### 2. `InputModule::DoShooting` (static) — orchestration

```cpp
static std::unique_ptr<InputModule> DoShooting(std::unique_ptr<InputModule> im);
```

1. **Collect targets:** query `im->all_species_` via `GetShootingTargets()` (non-empty only for
   species that guessed) + check the fc for the lone cosmological target `100*theta_s`. Seed the
   guesses/`dxdy` via each species' `ComputeShootingGuess` (and the module-level `theta_s` guess).
2. **None →** `return im;` (move-through; the common, zero-overhead case).
3. **Some →** run the existing Newton / ridder root finder over the flattened unknown vector
   (deterministic order = `theta_s` first, then `all_species_` lex
   order, each species filling its slots in order). Each iteration:
   - write the trial unknown value(s) into the fc keyed by `unknown_param`,
   - build a fresh `InputModule` and its `BackgroundModule` (+ `ThermodynamicsModule` only if
     `theta_s` is active) **directly** — no `Cosmology` (so iteration builds never re-enter the lazy
     shooting path; an iteration module has its unknown set, so it reports nothing to shoot anyway),
   - gather residuals in the same fixed order: `theta_s` from thermo (module-level);
     each species via `ComputeShootingResidual`.
   On convergence, write the solved unknown(s) into the fc, build **one** final `InputModule` from the
   resolved fc, and return it. The returned module is always fully resolved — no guessed state escapes.

### 3. Per-species hooks (the physics; no module downcasts)

Every target is **scalar** — the multi-value (comma-list) dncdm input was removed in favour of
dot-syntax, so each decaying-ncdm flavor is its own `DNCDM_DR` instance with its own per-instance keys
and shoots one scalar. The current `target_sizes` / comma-list plumbing is therefore dropped.

New `species/shooting_target.h`:

```cpp
struct ShootingResidualContext {
  const background* pba;
  const BackgroundModule* bgm;   // computed: gives the today background row + H0
};

struct ShootingTarget {
  std::string target_name;     // input vocabulary, e.g. "Omega_dcdmdr"  (was the enum + kTargetNamestrings_)
  std::string unknown_param;   // fc param this target varies, e.g. "Omega_ini_dcdm"  (was kUnknownNamestrings_)
  double target_value;         // the requested value
};
```

On `BaseSpecies` (default no-ops) — three single-purpose hooks:

```cpp
// WHAT this species is fitting. Non-empty iff it guessed (target set, unknown absent):
// drives target detection, the unknown vector, and the lazy-shooting trigger.
virtual std::vector<ShootingTarget> GetShootingTargets() const { return {}; }

// Initial guess + Jacobian seed for each unknown, same order as GetShootingTargets.
// Context-limited to SpeciesBuildContext (what CreateAll sees: pfc/pba/ppr/ncdm_settings) +
// the species' own constructed state. Also called by CreateAll to build from the guess.
virtual void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                  std::vector<double>& guess, std::vector<double>& dxdy) const {}

// computed-today minus requested target, one entry per target, same order. Valid whenever the
// target-form input is set (so it works during iterations where DoShooting has set the unknown).
virtual void ComputeShootingResidual(const ShootingResidualContext& ctx,
                                     std::vector<double>& residuals_out) const {}
```

Overridden by `DCDM_DR_Species`, `ScalarFieldSpecies`, `DNCDM_DR_Species`. The residual locates its
own data from `ctx.bgm` (e.g. `bgm->all_species_.at("DCDM_DR")`; each `DNCDM_DR` instance reads its
own flavor) — the downcast lives in **species code that knows its own type**, satisfying the
"no species-picking in modules" rule. A dedicated `ComputeShootingGuess` hook (rather than inlining
the guess in `CreateAll`) keeps the guess formulas — today scattered through `input_get_guess` —
visible and isolated; `CreateAll` calls it to obtain the value it builds from when the unknown is
absent.

`theta_s` is **not** a species hook: its guess (the `3.54·θ²−5.455·θ+2.548` fit, varying `h`) and
residual (`100·rs/ra − target`, from `ThermodynamicsModule`) stay in `InputModule` as the one
module-level target that `DoShooting` always handles.

### 4. `Cosmology` shoots lazily

A single `input_module_ptr_` member holds the raw (possibly-guessed) module at construction and is
replaced **in place** by the resolved module the first time anything is used:

```cpp
InputModulePtr& GetInputModule() {
  if (!shot_) { input_module_ptr_ = InputModule::DoShooting(std::move(input_module_ptr_)); shot_ = true; }
  return input_module_ptr_;
}
```

Every other getter (`GetBackgroundModule`, …) currently reads `input_module_ptr_` **directly**; each
is changed to read `GetInputModule()` so shooting fires before any module is built. Both `Cosmology`
constructors are kept (`Cosmology(unique_ptr<InputModule>)` is still used by `input_prepare_pk_eq` in
`nonlinear_module`); both just store into `input_module_ptr_`.

## Components touched

| File | Change |
|---|---|
| `species/shooting_target.h` | **new** — `ShootingTarget` + `ShootingResidualContext` |
| `species/base_species.h` | add `GetShootingTargets` / `ComputeShootingGuess` / `ComputeShootingResidual` no-op virtuals |
| `species/dcdm_dr_species.{h,cpp}` | the three hooks + guess-driven `CreateAll` (Omega_dcdmdr / omega_dcdmdr / Omega_ini_dcdm / omega_ini_dcdm) |
| `species/scalar_field.{h,cpp}` | same, for `Omega_scf` (unknown `scf_shooting_parameter`) |
| `species/dncdm_dr_species.{h,cpp}` | same, one scalar target per flavor — with the **corrected** unknown mapping |
| `source/input_module.{h,cpp}` | delete `enum target_names`, `kTargetNamestrings_`, `kUnknownNamestrings_`, `input_auxillary_target_conditions`, the two per-target switches, `FixUnknownParameters`'s enum plumbing; remove shooting from the ctor; add static `DoShooting` + module-level `theta_s` guess/residual; keep `input_find_root` / `fzero_Newton` / ridder helpers |
| `source/cosmology.{h,cpp}` | lazy `GetInputModule` (DoShooting + `shot_`); route other getters through it |

## Verification

- `make class` clean; `make classy`.
- **theta_s** (most common shoot), **Omega_dcdmdr** (`test_dcdm_dr_matches_reference`), **Omega_scf**
  (`test_class.py:161`) — strict ~0.1% behaviour-preserving (correct mappings, test-covered). Add an
  explicit `100*theta_s` reference test if none exists.
- Full 84-scenario grid for no crashes.
- **dncdm shooting (Omega_dncdmdr / deg_ncdm_decay_dr / Omega_ini_dncdm / Neff_ini_dncdm / …):**
  correct-by-construction — the refactor fixes the `kUnknownNamestrings_` misalignment, so results
  change vs the (broken, untested) master path. Verify each converges to a sensible value; author
  eyeballs. Add at least one dncdm-shooting smoke ini.

## Open implementation details (resolve in the plan)

- **Exact dncdm per-instance keys.** Each `DNCDM_DR` flavor shoots one scalar via its own dot-syntax
  keys (the comma-list vector form is gone). Confirm the precise target/unknown key names a flavor
  reads (e.g. `<name>.Omega_dncdmdr` → `<name>.deg`, and the deg/Neff/Omega_ini variants) when
  porting the guess/residual math.
- **`ComputeShootingGuess` for dncdm needs a constructed species** (`GetIni`, `ComputeMomenta`,
  `Gamma`, `M`, `GetDeg`). Confirm the build order inside `DNCDM_DR::CreateAll`: build the dncdm child
  from its direct inputs, call `ComputeShootingGuess` for the unknown (deg/Omega), then finalise — and
  that `SpeciesBuildContext` (bgm = nullptr at input time) carries everything the guess needs.
- **`theta_s` + species combined shoot.** Confirm the multi-unknown ordering (theta_s `h` plus species
  unknowns) matches today so combined shoots (e.g. `100*theta_s` with `Omega_dcdmdr`) still converge.
- Whether `file_content_.is_shooting` can be deleted entirely once shooting leaves the ctor.
