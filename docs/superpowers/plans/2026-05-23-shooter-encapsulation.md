# Shooter encapsulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move species-specific shooting (root-finding for inputs that can't be set directly) out of `InputModule`'s enum/switch dispatch into per-species hooks, and lift shooting orchestration out of the `InputModule` constructor into a static `InputModule::DoShooting` invoked lazily by `Cosmology`.

**Architecture:** Each shooting-capable species owns three hooks — `GetShootingTargets` (what it fits), `ComputeShootingGuess` (initial guess, also used by `CreateAll` to self-construct when the unknown is absent), `ComputeShootingResidual` (computed-minus-target). The `InputModule` ctor becomes pure (read + construct, guessing where needed); `Cosmology::GetInputModule` runs `DoShooting` once on first use. The one cosmological target (`100*theta_s`, unknown `h`) stays module-level. The `enum target_names` and its 5 string tables/switches are deleted.

**Tech Stack:** C++17, CLASSpp Makefile + classy, pytest (`python/test_class.py`).

**Spec:** `docs/superpowers/specs/2026-05-23-shooter-encapsulation-design.md`. Read its "Latent bug" and "Open implementation details" sections — they capture why dncdm results change and what to confirm while porting.

**Verification model:** `theta_s` / `Omega_dcdmdr` / `Omega_scf` are correct + test-covered → strict ~0.1% behaviour-preserving. dncdm shooting is latently misaligned + untested on master → the refactor fixes it by construction; verify convergence to sensible values and have the author eyeball (do NOT match master). Full 84-scenario grid for no crashes.

**Incremental safety:** the new machinery is built dormant alongside the old (Tasks 1–5), wired as a no-op while the old path still resolves first (Task 6), then activated by removing ctor-shooting (Task 7), then the dead code is deleted (Task 8). The build + tests stay green at every commit.

---

## File structure

| File | Responsibility after this PR |
|---|---|
| `species/shooting_target.h` | **new** — `ShootingTarget` descriptor + `ShootingResidualContext` |
| `species/base_species.h` | `GetShootingTargets` / `ComputeShootingGuess` / `ComputeShootingResidual` no-op virtuals |
| `species/dcdm_dr_species.{h,cpp}` | the 3 hooks + guess-driven `CreateAll` (targets: `Omega_dcdmdr`/`omega_dcdmdr`/`Omega_ini_dcdm`/`omega_ini_dcdm`) |
| `species/scalar_field.{h,cpp}` | the 3 hooks + guess-driven construction (target `Omega_scf`, unknown `scf_shooting_parameter`) |
| `species/dncdm_dr_species.{h,cpp}` | the 3 hooks + guess-driven construction, one scalar target per flavor (corrected unknown mapping) |
| `source/input_module.{h,cpp}` | add static `DoShooting` + module-level `theta_s`; delete the enum, tables, `input_auxillary_target_conditions`, the two switches, `FixUnknownParameters`, `input_try_unknown_parameters`, `input_get_guess`; remove shooting from the ctor; keep `input_find_root`/`fzero_Newton`/ridder |
| `source/cosmology.{h,cpp}` | lazy `GetInputModule` (DoShooting + `shot_`); other getters route through it |
| `species/all_species.h` | add `species/shooting_target.h` include if needed (it's pulled via base_species.h) |
| `Makefile`, `setup.py`, `CLASS.xcodeproj/project.pbxproj` | no new `.cpp` (shooting_target.h is header-only) — **no build-system change needed** |

---

## Task 1: `ShootingTarget` descriptor + `BaseSpecies` hooks (dormant)

**Files:** Create `species/shooting_target.h`; Modify `species/base_species.h`

- [ ] **Step 1: Create `species/shooting_target.h`**

```cpp
#pragma once
#include <string>

struct background;

/** Context for ComputeShootingResidual: a freshly-COMPUTED cosmology's "today" state. */
struct ShootingResidualContext {
  const background* pba       = nullptr;
  const double* bg_today      = nullptr;  // last row of BackgroundModule::background_table_
};

/** One scalar fittable target a species owns. (All shooting targets are scalar: the
 *  multi-value dncdm input was removed in favour of dot-syntax, one flavor per instance.) */
struct ShootingTarget {
  std::string target_name;    // input vocabulary, e.g. "Omega_dcdmdr"
  std::string unknown_param;  // fc param this target varies, e.g. "Omega_ini_dcdm"
  double target_value = 0.;   // requested value
};
```

- [ ] **Step 2: Add the three no-op hooks to `BaseSpecies`**

Edit `species/base_species.h`. Add `#include "shooting_target.h"` near the top, and forward-declare `struct SpeciesBuildContext;` (defined in `species_build_context.h`). Add these virtuals in the public section (next to the other Stage hooks, e.g. just after the `CopyPerturbationsAcrossSwitch` declaration):

```cpp
  // ── Shooter / root-finding ────────────────────────────────────────────────
  // Reported after construction; non-empty iff this species guessed its unknown
  // (target-form input set, direct unknown absent). Drives target detection, the
  // unknown vector, and the lazy-shooting trigger.
  virtual std::vector<ShootingTarget> GetShootingTargets() const {
    return {};
  }
  // Initial guess + Jacobian seed for each unknown, same order as GetShootingTargets.
  // Also called by CreateAll to obtain the value to build from when the unknown is absent.
  virtual void ComputeShootingGuess(const SpeciesBuildContext& /*ctx*/,
                                    std::vector<double>& /*guess*/,
                                    std::vector<double>& /*dxdy*/) const {}
  // computed-today minus requested target, appended one entry per target (same order).
  virtual void ComputeShootingResidual(const ShootingResidualContext& /*ctx*/,
                                       std::vector<double>& /*residuals_out*/) const {}
```

Ensure `<vector>` and `<string>` are included (base_species.h already uses `std::vector`).

- [ ] **Step 3: Build**

Run: `make class -j 2>&1 | tail -5`
Expected: builds clean (hooks are unused no-ops).

- [ ] **Step 4: Commit**

```bash
git add species/shooting_target.h species/base_species.h
git commit -m "species: add ShootingTarget descriptor + shooter hooks (no-op defaults)"
```

---

## Task 2: `DCDM_DR_Species` shooter hooks + guess-driven `CreateAll` (dormant)

**Files:** Modify `species/dcdm_dr_species.{h,cpp}`

The DCDM_DR targets and unknowns (from `kUnknownNamestrings_`, all correct for dcdm):
`Omega_dcdmdr`→`Omega_ini_dcdm`, `omega_dcdmdr`→`Omega_ini_dcdm`, `Omega_ini_dcdm`→`Omega_dcdmdr`,
`omega_ini_dcdm`→`omega_dcdmdr`.

- [ ] **Step 1: Read the current code to port**

```bash
sed -n '280,287p' species/dcdm_dr_species.cpp          # DCDM_DR_Species::CreateAll
sed -n '3550,3602p' source/input_module.cpp            # residual cases: Omega_dcdmdr/omega_dcdmdr/Omega_ini_dcdm/omega_ini_dcdm
sed -n '3692,3757p' source/input_module.cpp            # guess cases for the same
sed -n '1133,1160p' source/input_module.cpp            # how Omega0_dcdmdr / Omega_ini_dcdm are read into pba
```

- [ ] **Step 2: Declare the hooks in `species/dcdm_dr_species.h`**

Add to the public section (after `WriteOutputColumns`):

```cpp
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess, std::vector<double>& dxdy) const override;
  void ComputeShootingResidual(const ShootingResidualContext& ctx,
                               std::vector<double>& residuals_out) const override;
```

Add a private member to remember the active target (set in `CreateAll`):

```cpp
 private:
  // Set when a target-form input was given. unknown_param empty => no shooting target.
  ShootingTarget shooting_target_{};
  bool needs_shooting_ = false;   // true iff the direct unknown was absent (we guessed)
```

(Include `shooting_target.h` and `species_build_context.h` as needed.)

- [ ] **Step 3: Implement `ComputeShootingResidual` in `species/dcdm_dr_species.cpp`**

Port the residual arithmetic from `input_module.cpp:3550-3602`, reading `this`'s children and the
context's `bg_today` instead of `bam->all_species_.at("DCDM_DR")`:

```cpp
void DCDM_DR_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                              std::vector<double>& residuals_out) const {
  if (shooting_target_.unknown_param.empty())
    return;  // not fitting anything
  // Use ctx.pba consistently for all background params (it == the freshly-computed cosmology).
  const background& ba  = *ctx.pba;
  const double* bg      = ctx.bg_today;
  const double H0       = ba.H0;
  const double rho_dcdm = dcdm_->Rho(bg);
  const double rho_dr   = (ba.has_dr == _TRUE_) ? dr_sp_->Rho(bg) : 0.;
  const std::string& t  = shooting_target_.target_name;
  double residual;
  if (t == "Omega_dcdmdr") {
    residual = (rho_dcdm + rho_dr) / (H0 * H0) - shooting_target_.target_value;
  } else if (t == "omega_dcdmdr") {
    residual = (rho_dcdm + rho_dr) / (H0 * H0) - shooting_target_.target_value / ba.h / ba.h;
  } else {  // Omega_ini_dcdm / omega_ini_dcdm: target is the ini, we vary the today density
    residual = -(rho_dcdm + rho_dr) / (H0 * H0) + ba.Omega0_dcdmdr;
  }
  residuals_out.push_back(residual);
}
```

(Confirm `dcdm_`/`dr_sp_` are reachable; they are non-owning members. `pba_` is a member.)

- [ ] **Step 4: Implement `ComputeShootingGuess`**

Port the guess arithmetic from `input_module.cpp:3692-3757` (the `Omega_dcdmdr`/`omega_dcdmdr` and
`Omega_ini_dcdm`/`omega_ini_dcdm` cases), reading `ctx.pba` and `shooting_target_.target_value`:

```cpp
void DCDM_DR_Species::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                           std::vector<double>& guess,
                                           std::vector<double>& dxdy) const {
  const background& ba = *ctx.pba;
  const std::string& t = shooting_target_.target_name;
  const double tv      = shooting_target_.target_value;
  if (t == "Omega_dcdmdr" || t == "omega_dcdmdr") {
    double Omega_M = ba.Omega0_cdm + ba.Omega0_idm_dr + ba.Omega0_dcdmdr + ba.Omega0_b;
    double gamma = ba.Gamma_dcdm / ba.H0, a_decay = 1.0;
    if (gamma > 1) a_decay = pow(1 + (gamma * gamma - 1.) / Omega_M, -1. / 3.);
    double tgt = (t == "omega_dcdmdr") ? tv / ba.h / ba.h : tv;
    guess.push_back(tgt / a_decay);
    dxdy.push_back((tgt / a_decay) / tv);
  } else {  // Omega_ini_dcdm / omega_ini_dcdm
    double Omega0_dcdmdr = (t == "omega_ini_dcdm") ? tv / (ba.h * ba.h) : tv;
    double Omega_M = ba.Omega0_cdm + ba.Omega0_idm_dr + Omega0_dcdmdr + ba.Omega0_b;
    double gamma = ba.Gamma_dcdm / ba.H0, a_decay = 1.0;
    if (gamma > 1) a_decay = pow(1 + (gamma * gamma - 1.) / Omega_M, -1. / 3.);
    guess.push_back(tv * a_decay);
    double d = a_decay;
    if (gamma > 100) d *= gamma / 100;
    dxdy.push_back(d);
  }
}
```

- [ ] **Step 5: Implement `GetShootingTargets`**

```cpp
std::vector<ShootingTarget> DCDM_DR_Species::GetShootingTargets() const {
  if (needs_shooting_)
    return {shooting_target_};
  return {};
}
```

- [ ] **Step 6: Guess-driven construction in `DCDM_DR_Species::CreateAll`**

In `CreateAll` (currently `species/dcdm_dr_species.cpp:280-286`), after building the composite, detect
the target/unknown situation from `ctx.pfc` and, if the unknown is absent, set the unknown in `pba`
from the species' own guess so the rest of construction proceeds normally. Use `parser_read_double`
with a `flag` to detect presence. The active target is whichever of the four target-form keys is set;
its unknown is the partner per the mapping above. Record `shooting_target_` + `needs_shooting_`.

Pseudocode for the factory body (fill exact reads against Step 1's `:1133-1160`):

```cpp
auto composite = std::make_unique<DCDM_DR_Species>(ctx.pba, ctx.bgm);
// Determine which target-form key is set and its partner unknown:
//   Omega_dcdmdr / omega_dcdmdr  -> unknown Omega_ini_dcdm
//   Omega_ini_dcdm / omega_ini_dcdm -> unknown Omega_dcdmdr
// If the target key is set AND its unknown key is absent in ctx.pfc:
//   composite->shooting_target_ = {target_name, unknown_param, target_value};
//   composite->needs_shooting_ = true;
//   std::vector<double> g, d; composite->ComputeShootingGuess(ctx, g, d);
//   ctx.pfc->set(unknown_param, std::to_string(g[0]));   // build proceeds from the guess
//   (also write the guessed value into the pba field the species reads, if read already)
```

**Write-back by unknown (not always `Omega_ini_dcdm`):** the guess `g[0]` must be written into the
`pba` field that corresponds to the *unknown*, plus the fc key:
- unknown `Omega_ini_dcdm` → `pba->Omega_ini_dcdm = g[0]`
- unknown `Omega_dcdmdr`   → `pba->Omega0_dcdmdr  = g[0]`
- unknown `omega_dcdmdr`   → `pba->Omega0_dcdmdr  = g[0] / (h*h)`  (omega→Omega, matching how `input_read_parameters` interprets the fc key)

and `ctx.pfc->set(unknown_param, std::to_string(g[0]))`. (Writing `g[0]` to `Omega_ini_dcdm` in the
inverse case would clobber the user's target.) The correct unknown per target (from the original
`kUnknownNamestrings_`): `Omega_dcdmdr`→`Omega_ini_dcdm`, `omega_dcdmdr`→`Omega_ini_dcdm`,
`Omega_ini_dcdm`→`Omega_dcdmdr`, `omega_ini_dcdm`→`omega_dcdmdr`. `pba` is non-const underneath; mutate
via `const_cast<background*>(ctx.pba)` (mirrors `input_module.cpp:1152`).

**Ini-form presence (per spec decision 0):** for the inverse targets (`Omega_ini_dcdm`/
`omega_ini_dcdm`) the today density `Omega0_dcdmdr` is 0, so `has_dcdm`/`has_dr` are false and
`CreateAll` wouldn't build DCDM_DR at all. In `source/input_module.cpp` `ConstructSpecies` presence
determination (`:241-243`), set `has_dcdm`/`has_dr` `_TRUE_` when an ini-form dcdm target key
(`Omega_ini_dcdm` or `omega_ini_dcdm`) is present in the fc, not only when `Omega0_dcdmdr != 0`. Then
`DCDM_DR::CreateAll` runs, reads the ini-form target value from `ctx.pfc`, sets `pba->Omega_ini_dcdm`
(the user's value), guesses `Omega0_dcdmdr` (the unknown), and reports the target. The discovery
build's budget uses the guessed `Omega0_dcdmdr` (fine — it's refined by shooting).

**Dormant-safety:** the common targets (`Omega_dcdmdr`/`omega_dcdmdr`) already set `Omega0_dcdmdr != 0`,
so the ini-form presence branch is inert there; and while old ctor-shooting is still active (until
Task 7) the unknown is resolved into the fc before `ConstructSpecies`, so `needs_shooting_` stays false
and nothing changes at runtime.

- [ ] **Step 7: Build (dormant — old ctor-shooting still resolves first, so this never triggers yet)**

Run: `make class -j 2>&1 | tail -5`  → clean.
Run: `cd python && python -m pytest -q test_class.py -k "dcdm_dr_matches_reference" && cd ..`
Expected: PASS (behaviour unchanged — old shooting still active; new hooks unused).

- [ ] **Step 8: Commit**

```bash
git add species/dcdm_dr_species.h species/dcdm_dr_species.cpp
git commit -m "dcdm_dr: shooter hooks (GetShootingTargets/ComputeShootingGuess/ComputeShootingResidual) + guess-driven CreateAll"
```

---

## Task 3: `ScalarFieldSpecies` shooter hooks + guess-driven construction (dormant)

**Files:** Modify `species/scalar_field.{h,cpp}`

Target `Omega_scf`, unknown `scf_shooting_parameter` (writes `pba->scf_parameters[scf_tuning_index]`).

- [ ] **Step 1: Read the current code to port**

```bash
grep -n "CreateAll" species/scalar_field.cpp species/scalar_field.h
sed -n '3579,3587p' source/input_module.cpp     # residual: Omega_scf
sed -n '3716,3734p' source/input_module.cpp     # guess: Omega_scf
sed -n '1238,1290p' source/input_module.cpp     # how Omega0_scf / scf_parameters / scf_tuning_index are read
```

- [ ] **Step 2: Declare the three hooks** in `species/scalar_field.h` (mirror Task 2 Step 2) and add
the `shooting_target_` / `needs_shooting_` private members.

- [ ] **Step 3: Implement `ComputeShootingResidual`** — port `input_module.cpp:3579-3587`:

```cpp
void ScalarFieldSpecies::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                 std::vector<double>& residuals_out) const {
  if (shooting_target_.unknown_param.empty())
    return;
  // Omega_scf is filled to close the budget; compare to pba->Omega0_scf (as today).
  residuals_out.push_back(Rho(ctx.bg_today) / (ctx.pba->H0 * ctx.pba->H0) - ctx.pba->Omega0_scf);
}
```

- [ ] **Step 4: Implement `ComputeShootingGuess`** — port `input_module.cpp:3716-3734`:

```cpp
void ScalarFieldSpecies::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                              std::vector<double>& guess,
                                              std::vector<double>& dxdy) const {
  const background& ba = *ctx.pba;
  if (ba.scf_tuning_index == 0) {
    guess.push_back(sqrt(3.0 / ba.Omega0_scf));
    dxdy.push_back(-0.5 * sqrt(3.0) * pow(ba.Omega0_scf, -1.5));
  } else {
    guess.push_back(ba.scf_parameters[ba.scf_tuning_index]);
    dxdy.push_back(1.);
  }
}
```

- [ ] **Step 5: Implement `GetShootingTargets`** (mirror Task 2 Step 5).

- [ ] **Step 6: Guess-driven construction in `ScalarFieldSpecies::CreateAll`** — if `Omega_scf` is set
and `scf_shooting_parameter` is absent in `ctx.pfc`, record `shooting_target_ = {"Omega_scf",
"scf_shooting_parameter", Omega_scf_value}`, set `needs_shooting_`, compute the guess, and write it
into `pba->scf_parameters[scf_tuning_index]` (and the fc key) so construction proceeds.

- [ ] **Step 7: Build + verify dormant**

Run: `make class -j 2>&1 | tail -5` → clean.
Run: `cd python && python -m pytest -q test_class.py -k "scf or scalar" && cd ..` → PASS.

- [ ] **Step 8: Commit**

```bash
git add species/scalar_field.h species/scalar_field.cpp
git commit -m "scalar_field: shooter hooks + guess-driven construction (Omega_scf)"
```

---

## Task 4: NEW per-flavor dncdm today-density shooting (`nu_i.Omega_dncdmdr` → shoot `nu_i.deg`)

**Files:** Modify `species/dncdm_species.{h,cpp}` (dot-input + deg-guess + pending state),
`species/dncdm_dr_species.{h,cpp}` (the three shooter hooks).

The old flat dncdm targets are dead (rejected). This task **adds** a per-flavor today-density fit:
`nu_i.Omega_dncdmdr` (or `nu_i.omega_dncdmdr` = `/h²`) = target today-density of that flavor's
`dncdm + dr`; unknown varied = that flavor's `nu_i.deg`. The flavor is already constructed from
`nu_i.type = ncdm_decay_dr` (no presence issue). Unlike the COMMON-active state in Tasks 2-3, this
path is **not** dormant once activated — but until Task 7 the old ctor-shooting is still active and a
`nu_i.Omega_dncdmdr` (with no flat key) is simply read by the species (no old-enum target matches it),
so verify it at least builds; full activation is Task 7.

- [ ] **Step 1: Read the code to port + the dncdm input plumbing**
```bash
sed -n '3604,3651p' source/input_module.cpp     # residual: Omega_dncdmdr case (sums rho_dr+rho_dncdm)
sed -n '3759,3800p' source/input_module.cpp     # guess: Omega_dncdmdr case (ComputeMomenta at a_ini + decay corr)
sed -n '117,172p' species/dncdm_species.cpp      # per-instance dot-input reading (deg/Omega/Omega_ini/...)
sed -n '46,72p'   species/dncdm_species.cpp      # ApplyDncdmInitialClosure (SetDeg_from_Omega_ini pattern)
grep -n "Omega_ini_pending_\|read_double\|instance_name\|SetDegAndFactor\|GetIni\|ComputeMomenta" species/dncdm_species.{h,cpp}
```
Confirm how a per-instance dot key maps to the raw fc key (`<instance>.deg`), since that raw key is the
`unknown_param` `DoShooting` will `fc.set`. (The instance reader reads `"deg"` from `<instance>.deg`.)

- [ ] **Step 2 (DNCDMSpecies — pending state + accessor):** in `species/dncdm_species.h` add
`std::optional<double> Omega_dncdmdr_pending_;` (today-density target, Omega units) next to
`Omega_ini_pending_`, with accessor `const std::optional<double>& Omega_dncdmdr_pending() const`.

- [ ] **Step 3 (DNCDMSpecies — deg-guess method):** add
`std::pair<double,double> DegGuessFromOmegaToday(const SpeciesBuildContext& ctx, double Omega_target) const;`
porting `input_get_guess:3759-3800` for THIS flavor: `a_ini = GetIni(pr.a_ini_over_a_today_default*a_today, a_today, pr.tol_ncdm_initial_w)`, `z_ini=1/a_ini-1`, `ComputeMomenta(z_ini, nullptr, &rho_actual, ...)`, `rho_deg1 = rho_actual/GetDeg()` (GetDeg here is the *current* deg — guard >0), `Omega_deg1 = rho_deg1*a_ini^4/H0²`, then the `Gamma()/H0 > 1` fully-decayed branch (with `a_nr=3.15/M()`, `k_rad`, `erfc`, etc.) to get `Omega_ini`; return `{Omega_ini/Omega_deg1, (Omega_ini/Omega_deg1)/Omega_target}` (guess, dxdy). Substitute `ba→*ctx.pba`, `pr→*ctx.ppr`, `dncdm_sp→this`.

- [ ] **Step 4 (DNCDMSpecies — read the new dot target + build from guessed deg):** in the per-instance
read (`:140-168`), additionally `input.read_double("Omega_dncdmdr", ...)` and `"omega_dncdmdr"` (→ `/h²`).
Extend the "at most one of {deg, Omega_ini, omega_ini, Neff_ini}" exclusivity so that **`Omega_dncdmdr`
counts as an amount-spec too** — BUT `deg` may legitimately co-occur with `Omega_dncdmdr` during shooting
(DoShooting sets `nu_i.deg` alongside the `Omega_dncdmdr` target). Resolution: **if `deg` is present, use
it (`SetDegAndFactor`) and do NOT stash a pending target** (this is the iteration/resolved case); **else
if `Omega_dncdmdr`/`omega_dncdmdr` present, stash `Omega_dncdmdr_pending_ = Omega_target` and seed the
build with a guessed deg**: `auto [g, dxdy] = DegGuessFromOmegaToday(ctx, Omega_target); SetDegAndFactor(g);`.
Error only on genuinely conflicting *target* combos (e.g. `Omega_dncdmdr` together with `Omega_ini`/
`omega_ini`/`Neff_ini`/`Omega`).

- [ ] **Step 5 (DNCDM_DR — declare the 3 hooks)** in `species/dncdm_dr_species.h` (mirror Task 2's
declarations; no extra members needed — the state lives on the dncdm child).

- [ ] **Step 6 (DNCDM_DR — GetShootingTargets):**
```cpp
std::vector<ShootingTarget> DNCDM_DR_Species::GetShootingTargets() const {
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return {};
  return {{name() + ".Omega_dncdmdr", name() + ".deg", *dncdm_->Omega_dncdmdr_pending()}};
}
```
(`name()` is the flavor instance key, e.g. `"nu1"`; `unknown_param` is the raw `<instance>.deg` fc key.)

- [ ] **Step 7 (DNCDM_DR — ComputeShootingResidual):** today-density of this flavor's dncdm+dr:
```cpp
void DNCDM_DR_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                               std::vector<double>& residuals_out) const {
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return;
  const double* bg = ctx.bg_today;
  const double H0  = ctx.pba->H0;
  residuals_out.push_back((dncdm_->Rho(bg) + dr_sp_->Rho(bg)) / (H0 * H0)
                          - *dncdm_->Omega_dncdmdr_pending());
}
```

- [ ] **Step 8 (DNCDM_DR — ComputeShootingGuess):** delegate to the child:
```cpp
void DNCDM_DR_Species::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                            std::vector<double>& guess, std::vector<double>& dxdy) const {
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return;
  auto [g, d] = dncdm_->DegGuessFromOmegaToday(ctx, *dncdm_->Omega_dncdmdr_pending());
  guess.push_back(g);
  dxdy.push_back(d);
}
```

- [ ] **Step 9: Build + verify** — `make class -j 2>&1 | tail -8` clean;
`cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py; cd ..`
→ 84 pass (existing scenarios don't use the new key). Smoke that a `nu1.Omega_dncdmdr` input at least
parses+builds (full convergence is checked in Task 7 once shooting is active):
`./class` on a tiny ini with `nu1.type=ncdm_decay_dr, nu1.m=1, nu1.Gamma=1e3, nu1.Omega_dncdmdr=0.01, output=` — it should build (it goes through the OLD shooter which has no matching target, so the species just builds from the guessed deg). Surface to author.

- [ ] **Step 10: Commit**
```bash
git add species/dncdm_species.h species/dncdm_species.cpp species/dncdm_dr_species.h species/dncdm_dr_species.cpp
git commit -m "dncdm: add per-flavor today-density shooting (nu.Omega_dncdmdr -> shoot nu.deg) via DNCDM_DR hooks"
```

---

## Task 5: `InputModule::DoShooting` (static) + module-level `theta_s` (dormant — not yet wired)

**Files:** Modify `source/input_module.{h,cpp}`

- [ ] **Step 1: Declare in `source/input_module.h`** (public, so `Cosmology` can call it):

```cpp
  static std::unique_ptr<InputModule> DoShooting(std::unique_ptr<InputModule> input_module);
```

- [ ] **Step 2: Implement `DoShooting`** in `source/input_module.cpp`. It collects targets, and if any,
runs the existing Newton root finder. Structure (reuse the residual-eval shape of the existing
`input_try_unknown_parameters:3508-3540` for building+computing per iteration, but read residuals via
species hooks + module-level theta_s):

```cpp
std::unique_ptr<InputModule> InputModule::DoShooting(std::unique_ptr<InputModule> im) {
  // 1. Collect targets (deterministic order: theta_s first, then all_species_ lex order).
  struct Slot { std::string unknown_param; std::string target_name; double target_value; bool is_theta_s; };
  std::vector<Slot> slots;
  std::vector<double> xguess, dxdy;
  // theta_s (cosmological): present iff "100*theta_s" is in the fc.
  // (read it via parser_read_double on im->file_content_; push its guess from the theta_s fit.)
  // species:
  for (auto& [key, sp] : im->all_species_) {
    for (const ShootingTarget& t : sp->GetShootingTargets()) {
      slots.push_back({t.unknown_param, t.target_name, t.target_value, false});
      std::vector<double> g, d;
      // rebuild a SpeciesBuildContext from im to call the guess (or have the species cache it)
      sp->ComputeShootingGuess(ctx, g, d);
      xguess.push_back(g[0]); dxdy.push_back(d[0]);
    }
  }
  if (slots.empty()) return im;   // nothing to shoot — the common path

  // 2. Run fzero_Newton over a residual functor that, per evaluation:
  //    - sets each slot's unknown_param in the fc to the trial value,
  //    - builds a fresh InputModule + its BackgroundModule (+ ThermodynamicsModule iff a theta_s slot),
  //    - fills residuals: theta_s = 100*thm->rs_rec_/thm->ra_rec_ - target;
  //      species via ComputeShootingResidual(ctx{pba, bg_today}, out) iterating all_species_ in lex order.
  //  (Mirror input_try_unknown_parameters' fc.set + build + read; mirror input_find_root's Newton call.)

  // 3. On convergence, write solved unknowns into the fc, build one final InputModule, return it.
  // (set is_shooting handling is no longer needed — see Task 8.)
}
```

**Port references:** the per-iteration "set fc + build cosmology + read computed quantity" is
`input_try_unknown_parameters:3494-3540`; the theta_s residual is the `case theta_s` at
`:3545-3548` (`100. * thm->rs_rec_ / thm->ra_rec_ - target`); the theta_s guess + `h` update is
`input_get_guess:3683-3690`; the Newton invocation is `input_find_root:3859-` /
`fzero_Newton` (`tools/evolver_ndf15.cpp:1145`, signature `func, x_inout, dxdF, x_size, tolx, tolF,
param, fevals, errmsg`). Build the iteration's `BackgroundModule` directly
(`new BackgroundModule(iteration_im_ptr)`) and `ThermodynamicsModule` only if a theta_s slot is
present — **no `Cosmology`** (keeps iteration builds off the lazy-shooting path).

`bg_today` = `bgm->background_table_.data() + (bgm->bt_size_ - 1) * bgm->bg_size_`.

- [ ] **Step 3: Build (unused — nothing calls DoShooting yet)**

Run: `make class -j 2>&1 | tail -5` → clean.

- [ ] **Step 4: Commit**

```bash
git add source/input_module.h source/input_module.cpp
git commit -m "input_module: add static DoShooting (species hooks + module-level theta_s); not yet wired"
```

---

## Task 6: Wire `Cosmology` to shoot lazily — as a no-op (old ctor-shooting still active)

**Files:** Modify `source/cosmology.{h,cpp}`

At this point the InputModule ctor STILL shoots, so a Cosmology's stored module is already resolved →
`DoShooting` finds no targets → returns it unchanged. This step is therefore behaviour-neutral and
proves the wiring before we flip the switch in Task 7.

- [ ] **Step 1: `source/cosmology.h`** — add `bool shot_ = false;` and keep both ctors storing into
`input_module_ptr_`. Change `GetInputModule()` to:

```cpp
InputModulePtr& GetInputModule() {
  if (!shot_) {
    input_module_ptr_ = InputModule::DoShooting(std::move(input_module_ptr_));
    shot_ = true;
  }
  return input_module_ptr_;
}
```

(Move the `GetInputModule` body to `cosmology.cpp` if it can't be inline due to `InputModule`
completeness; the header includes `input_module.h`, so inline is fine.)

- [ ] **Step 2: Route every other getter through `GetInputModule()`** in `source/cosmology.cpp`:
replace each direct `input_module_ptr_` use in `GetBackgroundModule` / `GetThermodynamicsModule` /
`GetPerturbationsModule` / `GetPrimordialModule` / `GetNonlinearModule` / `GetTransferModule` /
`GetSpectraModule` / `GetLensingModule` with `GetInputModule()`.

- [ ] **Step 3: Build + full verification (must be unchanged)**

```bash
make class -j 2>&1 | tail -5 && make classy 2>&1 | tail -3
cd python && python -m pytest -q test_class.py -k "dcdm_dr_matches_reference or scf or theta" && \
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py; cd ..
```
Expected: all PASS, identical to before (DoShooting is a no-op while ctor-shooting resolves first).

- [ ] **Step 4: Commit**

```bash
git add source/cosmology.h source/cosmology.cpp
git commit -m "cosmology: lazy shooting via GetInputModule (no-op while ctor-shooting active)"
```

---

## Task 7: Remove shooting from the `InputModule` ctor (activate the new path)

**Files:** Modify `source/input_module.cpp`

- [ ] **Step 1: Delete the in-ctor shooting** in `input_init()` — the target-detection loop
(`:545-574`) and the `if ((unknown_parameters_size > 0) && !file_content_.is_shooting) FixUnknownParameters(...)`
block (`:582-586`). Leave `input_read_precisions()` and `input_read_parameters()` in place. Now a
target-form-only run reaches `ConstructSpecies` with the unknown absent → species self-construct from
their guess and report targets.

- [ ] **Step 2: Verify `input_read_parameters` tolerates an absent unknown.** With ctor-shooting gone,
the outer `input_read_parameters` runs once with e.g. `Omega_ini_dcdm` absent (it stays default; the
species guesses in `CreateAll`). Confirm no `class_test` requires the unknown to be present. Build:

Run: `make class -j 2>&1 | tail -10`
Expected: clean. If a `class_test` fires for an absent unknown, relax it (the species now supplies the
guess) — note the exact test in the commit.

- [ ] **Step 3: Strict behaviour-preserving tests (theta_s / dcdm / scf)**

```bash
make classy 2>&1 | tail -3
cd python && python -m pytest -v test_class.py -k "dcdm_dr_matches_reference or scf or theta or rs_drag" ; cd ..
```
Expected: PASS. These exercise the now-active DoShooting path for the correct, test-covered targets.
If theta_s isn't covered, add the smoke test from Task 9 Step 1 first and run it here.

- [ ] **Step 4: No-crash grid + dcdm/scf shooting smoke**

```bash
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py; cd ..
./class explanatory.ini > /tmp/exp.log 2>&1 && echo OK || tail -20 /tmp/exp.log
```
Expected: 84 pass; explanatory.ini (which sets `100*theta_s`) runs clean.

- [ ] **Step 5: dncdm shooting — correct-by-construction (author eyeballs)**

Run a dncdm-shooting ini (target form), confirm convergence + sensible `Omega_dncdmdr` today:

```bash
cat > /tmp/dncdm_shoot.ini <<'EOF'
h=0.6736
omega_b=0.02237
omega_cdm=0.12
output=
background_verbose=2
input_verbose=2
dncdm1.type=ncdm_decay_dr
dncdm1.m=1.0
dncdm1.Gamma=1e3
Omega_dncdmdr=0.01
EOF
./class /tmp/dncdm_shoot.ini 2>&1 | tail -20
```
(Adjust keys to the confirmed dot-syntax from Task 4 Step 1.) Expected: shooting iterations converge;
`Omega_dncdmdr (computed) ≈ 0.01`. Surface to the author to eyeball (this path was broken on master).

- [ ] **Step 6: Commit**

```bash
git add source/input_module.cpp
git commit -m "input_module: remove shooting from the ctor; Cosmology::DoShooting is now the live path"
```

---

## Task 8: Delete the dead enum-dispatch machinery

**Files:** Modify `source/input_module.{h,cpp}`

- [ ] **Step 1: Delete from `source/input_module.h`** — `enum target_names` (`:25-38`),
`#define _NUM_TARGETS_`, the `fzerofun_workspace` struct (`:62-70`) if now unused, the
`kTargetNamestrings_`/`kUnknownNamestrings_` declarations, and the declarations of
`FixUnknownParameters`, `input_try_unknown_parameters`, `input_fzerofun_1d`, `input_get_guess`,
`input_auxillary_target_conditions`, and `shooting_workspace_`. **Keep** `input_find_root`,
`class_fzero_ridder` (if still used by `DoShooting`/`input_find_root`).

- [ ] **Step 2: Delete from `source/input_module.cpp`** — the definitions of `kTargetNamestrings_`/
`kUnknownNamestrings_` (`:49-73`), `FixUnknownParameters`, `input_try_unknown_parameters`,
`input_get_guess` (old), `input_fzerofun_1d`, `input_auxillary_target_conditions`, and the
`shooting_workspace_` initialisation in the ctor (`:194`). Fold whatever `input_find_root` /
`class_fzero_ridder` logic `DoShooting` still needs into `DoShooting` (or keep them as private static
helpers it calls).

- [ ] **Step 3: Remove `is_shooting`** from `FileContent` if no longer referenced.

```bash
grep -rn "is_shooting\|target_names\|fzerofun_workspace\|kTargetNamestrings_\|kUnknownNamestrings_\|input_try_unknown_parameters\|input_auxillary_target_conditions\|FixUnknownParameters" source/ include/ | grep -v "DoShooting"
```
Expected: no matches (or only inside `DoShooting`).

- [ ] **Step 4: Build + full verification**

```bash
make class -j 2>&1 | tail -10 && make classy 2>&1 | tail -3
cd python && python -m pytest -q test_class.py -k "dcdm_dr_matches_reference or scf or theta or rs_drag" && \
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py; cd ..
```
Expected: clean build (no new warnings), all PASS.

- [ ] **Step 5: Commit**

```bash
git add source/input_module.h source/input_module.cpp
git commit -m "input_module: delete dead enum-based shooting dispatch (target_names, tables, switches)"
```

---

## Task 9: Verification scenarios + clang-format + final check

**Files:** Modify `python/test_class.py` (add a theta_s reference test if absent); run clang-format.

- [ ] **Step 1: Add a `100*theta_s` reference regression** to `python/test_class.py` (mirror
`test_dcdm_dr_matches_reference` at `:742-763`), comparing candidate vs `classyref` for a scenario
with `'100*theta_s': 1.041783` (+ default cosmology, `'output': 'tCl'`), if no such test exists:

```python
def test_theta_s_shooting_matches_reference(self):
    scenario = {'100*theta_s': 1.041783, 'output': 'tCl', 'l_max_scalars': 200}
    candidate, reference = self._compute_candidate_and_reference(scenario)
    try:
        status = self.compare_output(reference, "Reference", candidate, "Candidate",
                                     COMPARE_CL_RELATIVE_ERROR, COMPARE_PK_RELATIVE_ERROR)
        self.assertTrue(status, "Reference comparison failed for theta_s shooting")
    finally:
        reference.struct_cleanup(); reference.empty()
        candidate.struct_cleanup(); candidate.empty()
```

- [ ] **Step 2: Run the shooting-relevant tests**

```bash
cd python && python -m pytest -v test_class.py -k "theta_s_shooting or dcdm_dr_matches_reference or scf" ; cd ..
```
Expected: PASS.

- [ ] **Step 3: clang-format the modified C++ files**

```bash
git diff --name-only master...HEAD -- '*.cpp' '*.h' | xargs clang-format -i
git diff --name-only master...HEAD -- '*.cpp' '*.h' | xargs clang-format --dry-run --Werror 2>&1 | head
```
Expected: format check clean.

- [ ] **Step 4: Full final verification**

```bash
make class -j 2>&1 | grep -iE "warning|error" ; echo "build done"
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py; cd ..
```
Expected: clean build, 84 pass.

- [ ] **Step 5: Commit**

```bash
git add python/test_class.py $(git diff --name-only -- '*.cpp' '*.h')
git commit -m "test: theta_s shooting reference regression; clang-format shooter refactor"
```

---

## Self-Review

**Spec coverage:**
- InputModule ctor becomes pure → Task 7. ✓
- Static `DoShooting`, lazy via Cosmology getter → Tasks 5, 6. ✓
- Per-species `GetShootingTargets`/`ComputeShootingGuess`/`ComputeShootingResidual` + guess-driven `CreateAll` → Tasks 1–4. ✓
- `theta_s` module-level → Task 5 (DoShooting). ✓
- Delete enum/tables/switches/`input_auxillary_target_conditions`/`FixUnknownParameters` → Task 8. ✓
- Scalar targets (no comma-list) → descriptor in Task 1; per-flavor dncdm in Task 4. ✓
- `kUnknownNamestrings_` `A_s` bug corrected by construction → Task 4 (corrected mapping) + verification model. ✓
- Lazy Cosmology single member + `shot_`, getters routed → Task 6. ✓
- Verification (theta_s/dcdm/scf strict; dncdm correct-by-construction; 84-grid) → Tasks 6–9. ✓

**Placeholder scan:** Novel code (descriptor, hooks, DoShooting skeleton, Cosmology, residual/guess for dcdm+scf) is written out. The dncdm guess body and the per-species `CreateAll` edits are specified as ports from exact source line ranges with explicit substitutions + a confirm-the-keys grep step (the residual/guess math is real in-repo code, not vague) — matching the PR-B plan convention. The two genuinely deferred specifics (exact dncdm dot keys; DNCDM_DR build order) are isolated to Task 4 with the exact greps to resolve them.

**Type consistency:** `ShootingTarget {target_name, unknown_param, target_value}`, `ShootingResidualContext {pba, bg_today}`, and the three hook signatures are used identically across Tasks 1–6. `GetShootingTargets`/`ComputeShootingGuess`/`ComputeShootingResidual` names match throughout.

**Risks flagged inline:** `input_read_parameters` tolerating an absent unknown (Task 7 Step 2); DNCDM_DR build order for the guess (Task 4 Step 6); exact dncdm dot keys (Task 4 Step 1).
