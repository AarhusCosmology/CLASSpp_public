# DNCDM_DR Closure + Uniform `GetOmega0()` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `DNCDM_DR` models close the flatness budget (`Ω_total(z=0)=1`) by giving the composite the same combined-density reservation `DCDM_DR` already has, fix the dncdm matter-child `Omega0_` omission, add the initial-mode fixed-point shoot, and reject the ambiguous bare `.Omega` key.

**Architecture:** `GetOmega0()` becomes uniformly "the species' total density fraction today." `DNCDM_DR_Species` overrides it to return the combined `Omega_dncdmdr`, so the Pass-1 closure (`input_module.cpp:244`) reserves the whole sector and shooting drives the integrated `(ρ_dncdm+ρ_dr)` to match — mirroring `DCDMSpecies::GetOmega0()→Omega0_dcdmdr_`. The dncdm matter child's `Omega0_` is backfilled post-integration (like `Omega0_dcdm_`). The initial mode gains the DCDM-style `Omega_dncdmdr` fixed-point shoot.

**Tech Stack:** C++17; build `make class -j`; behavioural tests are scenario `.ini` runs of `./class` whose verbose budget print is asserted by a small Python harness. Full regression via `python/test_class.py` + `test/scenarios/compare_tol.py`.

**Design spec:** `docs/superpowers/specs/2026-06-01-dncdm-dr-closure-design.md`

**Branch:** `dncdm-dr-closure` (already created off `master`; the budget-print `GetNDecayDr()` fix and the spec are already committed there).

---

## File structure

- `test/scenarios/assert_budget.py` (create) — runs `./class` on an ini with `background_verbose=2`, parses the verbose budget print, asserts `TOTAL Ω≈1` and optional named-line checks. The single behavioural test driver.
- `test/scenarios/dncdm_dr.ini` (modify) — currently uses the to-be-rejected `dncdm1.Omega`; repoint to a supported mode.
- `test/scenarios/dncdm_dr_combined.ini`, `dncdm_dr_initial.ini` (create) — closing scenarios for the two modes.
- `species/dncdm_dr_species.{h,cpp}` (modify) — `GetOmega0()` override; initial-mode shooting hooks.
- `species/dncdm_species.{h,cpp}` (modify) — reject bare `.Omega`; relax `Omega_dncdmdr`↔initial co-occurrence during shooting; expose the combined for the composite; flag initial mode for shooting.
- `source/background_module.cpp` (modify) — backfill dncdm matter-child `Omega0_` post-integration.
- `source/input_module.cpp` (modify) — delete the dead `omega0_ncdm_tot` tally.

---

## Task 1: Budget-closure assertion harness

**Files:**
- Create: `test/scenarios/assert_budget.py`

- [ ] **Step 1: Write the harness**

```python
#!/usr/bin/env python3
"""Run ./class on a scenario with background_verbose=2 and assert the printed
budget equation. Usage:
    assert_budget.py <ini> [--tol 0.02] [--positive-line "Neutrino Species"] [--expect-error SUBSTR]
- default: assert the "TOTAL ... Omega = X" line has |X-1| <= tol.
- --positive-line NAME: also assert a "-> NAME ... Omega = Y" line exists with Y > 0.
- --expect-error SUBSTR: instead expect ./class to FAIL with SUBSTR in its output.
Exits 0 on success, 1 on failure (prints reason + captured output tail).
"""
import argparse, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLASS = os.path.join(ROOT, "class")

def run_class(ini_text):
    with tempfile.TemporaryDirectory() as d:
        ini = os.path.join(d, "run.ini")
        with open(ini, "w") as f:
            f.write(ini_text)
        p = subprocess.run([CLASS, ini], capture_output=True, text=True)
        return p.returncode, p.stdout + p.stderr

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ini")
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--positive-line", default=None)
    ap.add_argument("--expect-error", default=None)
    a = ap.parse_args()

    with open(a.ini) as f:
        base = [ln for ln in f.read().splitlines()
                if not re.match(r"\s*(background_verbose|root)\s*=", ln)]
    with tempfile.TemporaryDirectory() as d:
        ini_text = "\n".join(base) + f"\nbackground_verbose = 2\nroot = {d}/out_\n"
        rc, out = run_class(ini_text)

    if a.expect_error is not None:
        if rc == 0 or a.expect_error not in out:
            print(f"FAIL: expected error containing {a.expect_error!r}; rc={rc}\n{out[-1500:]}")
            return 1
        print(f"OK: rejected with {a.expect_error!r}")
        return 0

    if rc != 0:
        print(f"FAIL: ./class exited {rc}\n{out[-1500:]}")
        return 1
    m = re.search(r"TOTAL\s+Omega = ([-0-9.eE+]+)", out)
    if not m:
        print(f"FAIL: no TOTAL line in budget output\n{out[-1500:]}")
        return 1
    total = float(m.group(1))
    if abs(total - 1.0) > a.tol:
        print(f"FAIL: TOTAL Omega = {total}, |1-total| > tol {a.tol}\n{out[-1500:]}")
        return 1
    if a.positive_line:
        pm = re.search(rf"->\s+{re.escape(a.positive_line)}\s+Omega = ([-0-9.eE+]+)", out)
        if not pm or float(pm.group(1)) <= 0.0:
            print(f"FAIL: line {a.positive_line!r} missing or <= 0\n{out[-1500:]}")
            return 1
    print(f"OK: TOTAL Omega = {total}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Make it executable and smoke-test against an existing closing model**

Run:
```bash
chmod +x test/scenarios/assert_budget.py
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/gauge_dcdm.ini --tol 0.02
```
Expected: `OK: TOTAL Omega = 1` (DCDM_DR already closes — confirms the harness reads the print correctly).

- [ ] **Step 3: Commit**

```bash
git add test/scenarios/assert_budget.py
git commit -m "test: add budget-closure assertion harness for ./class scenarios"
```

---

## Task 2: Combined-mode closure — `DNCDM_DR_Species::GetOmega0()` override

The composite must report the sector's combined today fraction so the Pass-1 closure (`input_module.cpp:229`) reserves it. The combined value lives on the DNCDM child as `Omega_dncdmdr_pending_` (set whenever `Omega_dncdmdr` is in the file content — true in combined mode, and after Task 5 in initial-shooting iterations).

**Files:**
- Create: `test/scenarios/dncdm_dr_combined.ini`
- Modify: `species/dncdm_dr_species.h` (add override decl), `species/dncdm_dr_species.cpp` (add override)
- Test: `test/scenarios/assert_budget.py`

- [ ] **Step 1: Write the failing test (combined-mode scenario)**

Create `test/scenarios/dncdm_dr_combined.ini`:
```ini
# DNCDM_DR closure smoke (combined-today normalization). Background-only.
T_cmb = 2.7255
omega_b = 0.022032
omega_cdm = 0.10
input_verbose = 1
dncdm1.type = ncdm_decay_dr
dncdm1.m = 1.0
dncdm1.T = 0.71611
dncdm1.Gamma = 1e2
dncdm1.Omega_dncdmdr = 0.12
```

- [ ] **Step 2: Run to verify it fails (budget open)**

Run:
```bash
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_combined.ini --tol 0.02
```
Expected: `FAIL: TOTAL Omega = ~1.0x` (composite still reports matter-only; DR not reserved).

- [ ] **Step 3: Add the `GetOmega0()` override declaration**

In `species/dncdm_dr_species.h`, in the `public:` accessor region (near the `dr()`/`dncdm()` accessors, ~line 69), add:
```cpp
  /** Total density fraction of the whole decaying sector today (matter + decay
   *  radiation) = the closure-reserved Omega_dncdmdr. Mirrors DCDMSpecies::GetOmega0().
   *  Surfaces the pinned/shot combined rather than summing children, because at
   *  Pass-1 closure time the emergent DR is not yet integrated. */
  double GetOmega0() const override;
```

- [ ] **Step 4: Implement the override**

In `species/dncdm_dr_species.cpp`, add (near the other `DNCDM_DR_Species::` member definitions):
```cpp
double DNCDM_DR_Species::GetOmega0() const {
  // Combined mode: Omega_dncdmdr is the user input; initial-shooting iterations:
  // the shooter writes <flavor>.Omega_dncdmdr, read into Omega_dncdmdr_pending_.
  if (dncdm_->Omega_dncdmdr_pending().has_value())
    return *dncdm_->Omega_dncdmdr_pending();
  // Pre-shooting discovery build with no pinned value yet: fall back to the
  // child-sum (matter child + DR child(0)). Replaced by the pinned value once
  // DoShooting writes the unknown.
  return CompositeSpecies::GetOmega0();
}
```

- [ ] **Step 5: Build and verify the combined budget now closes**

Run:
```bash
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_combined.ini --tol 0.02
```
Expected: `OK: TOTAL Omega = ~1.0`.

- [ ] **Step 6: Commit**

```bash
git add species/dncdm_dr_species.h species/dncdm_dr_species.cpp test/scenarios/dncdm_dr_combined.ini
git commit -m "DNCDM_DR: GetOmega0() returns combined sector density so closure reserves the DR"
```

---

## Task 3: Backfill the dncdm matter-child `Omega0_` post-integration

In combined/initial modes the matter child only sets `deg`, never `Omega0_`, so `dncdm().GetOmega0()` returns 0 — dropping the decaying species from `fnu` and the budget-print neutrino line. Set it from the integrated density today, mirroring `Omega0_dcdm_` (`background_module.cpp:925-934`).

**Files:**
- Modify: `source/background_module.cpp` (the post-integration block ~925-934)
- Test: `test/scenarios/assert_budget.py` (the `--positive-line` check)

- [ ] **Step 1: Confirm the failure — neutrino line is 0 in combined mode**

Run:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_combined.ini --tol 0.02 --positive-line "Neutrino Species"
```
Expected: `FAIL: line 'Neutrino Species' missing or <= 0` (matter child `Omega0_` unset).

- [ ] **Step 2: Add the backfill**

In `source/background_module.cpp`, immediately after the existing DNCDM_DR `Omega0_dr_` accumulation loop (the loop at ~931-934 that reads `dncdm_dr->dr().bi_rho_index()`), set each composite's matter-child `Omega0_` from its integrated density today. Use the species' own `Rho()` at the final (today) background row:
```cpp
  // Backfill each DNCDM matter child's today density fraction (parallels Omega0_dcdm_).
  // In combined/initial modes the child only carries `deg`; without this its GetOmega0()
  // stays 0, dropping it from fnu (GetOmega0NcdmTot) and the budget-print neutrino line.
  {
    const double* bg_today = background_table_.data() + (bt_size_ - 1) * bg_size_;
    for (auto& [key, sp] : all_species_) {
      if (auto* dncdm_dr = dynamic_cast<DNCDM_DR_Species*>(sp.get())) {
        const double omega0_matter = dncdm_dr->dncdm().Rho(bg_today) / pba->H0 / pba->H0;
        dncdm_dr->dncdm().SetOmega0(omega0_matter, pba->h);
      }
    }
  }
```
(`SetOmega0(double, double)` already exists on `NCDMBaseSpecies`, `ncdm_base_species.h:194`. `dncdm()` returns a mutable `DNCDMSpecies&`, `dncdm_dr_species.h:58`.)

- [ ] **Step 3: Build and verify the neutrino line is now positive**

Run:
```bash
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_combined.ini --tol 0.02 --positive-line "Neutrino Species"
```
Expected: `OK: TOTAL Omega = ~1.0` (TOTAL still closes; neutrino line now > 0).

- [ ] **Step 4: Verify stable-NCDM is unaffected (Omega0_ already set at construction)**

Run:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/ncdm_dncdm_idmdr_combined.ini --tol 0.02 --positive-line "Neutrino Species"
```
Expected: `OK` (the combined multi-species scenario still closes with a positive neutrino bucket).

- [ ] **Step 5: Commit**

```bash
git add source/background_module.cpp
git commit -m "background: backfill DNCDM matter-child Omega0_ from integrated today density"
```

---

## Task 4: Delete the dead `omega0_ncdm_tot` tally

`omega0_ncdm_tot` (`input_module.cpp:221` decl, `:231` accumulate) is never read. Remove it; keep the `n_ncdm`/`n_dncdm` counters (used at line 274).

**Files:**
- Modify: `source/input_module.cpp:218-238`

- [ ] **Step 1: Remove the declaration and accumulation**

In `source/input_module.cpp`, delete the `double omega0_ncdm_tot = 0.0;` declaration (line 221) and the `omega0_ncdm_tot += e.species->GetOmega0();` line (231). Leave `n_ncdm`, `n_dncdm`, and the surrounding `is_ncdm_family` logic intact. The block becomes:
```cpp
  double omega0_sum      = 0.0;
  int n_ncdm             = 0;
  int n_dncdm            = 0;
  for (const auto& entry : kAllSpeciesFactories) {
    if (entry.name == closure_name)
      continue;
    auto produced             = entry.create_all(ctx);
    const bool is_ncdm_family = (entry.name == "NCDM" || entry.name == "DNCDM_DR" ||
                                 entry.name == "NCDMInt");
    for (auto& e : produced) {
      omega0_sum += e.species->GetOmega0();
      if (is_ncdm_family) {
        n_ncdm += 1;
        if (entry.name == "DNCDM_DR")
          n_dncdm += 1;
      }
      all_species_.insert(e.key, std::move(e.species));
    }
  }
```

- [ ] **Step 2: Build and confirm no behavioural change**

Run:
```bash
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_combined.ini --tol 0.02
```
Expected: `OK` (unchanged — the removed variable fed nothing).

- [ ] **Step 3: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: drop dead omega0_ncdm_tot tally (never read since pba->Omega0_ncdm_tot removed)"
```

---

## Task 5: Initial-mode `Omega_dncdmdr` fixed-point shoot

Give the initial-abundance modes (`Omega_ini`/`omega_ini`/`Neff_ini`) the DCDM-style fixed-point shoot so closure reserves the true combined. The given initial seeds `deg` (already done by `ApplyDncdmInitialClosure`); we additionally shoot `<flavor>.Omega_dncdmdr` so that `GetOmega0()` reserves the integrated combined. This mirrors `DCDM_DR_Species`' "else" residual branch (`dcdm_dr_species.cpp:340-351,369-370`).

**Files:**
- Create: `test/scenarios/dncdm_dr_initial.ini`
- Modify: `species/dncdm_species.h`/`.cpp` (track initial mode; relax co-occurrence; read shot `Omega_dncdmdr`), `species/dncdm_dr_species.cpp` (`GetShootingTargets`/`ComputeShootingGuess`/`ComputeShootingResidual` initial branch)
- Test: `test/scenarios/assert_budget.py`

- [ ] **Step 1: Write the failing test (initial-mode scenario)**

Create `test/scenarios/dncdm_dr_initial.ini`:
```ini
# DNCDM_DR closure smoke (initial-abundance normalization). Background-only.
T_cmb = 2.7255
omega_b = 0.022032
omega_cdm = 0.10
input_verbose = 1
dncdm1.type = ncdm_decay_dr
dncdm1.m = 1.0
dncdm1.T = 0.71611
dncdm1.Gamma = 1e2
dncdm1.Omega_ini = 0.12
```

- [ ] **Step 2: Run to verify it fails (no shoot → budget open)**

Run:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_initial.ini --tol 0.02
```
Expected: `FAIL: TOTAL Omega = ~1.0x` (initial mode sets `deg` directly and never shoots).

- [ ] **Step 3: Expose the initial-mode flag + a shot-`Omega_dncdmdr` reader on the child**

In `species/dncdm_species.h`, add a public predicate (near `Omega_ini_pending()`/`Neff_ini_pending()`, ~line 38):
```cpp
  /** True iff this flavor is normalized by initial abundance (Omega_ini/omega_ini/Neff_ini),
   *  i.e. the mode that needs the Omega_dncdmdr fixed-point shoot for closure. */
  bool InitialAbundanceMode() const {
    return Omega_ini_pending_.has_value() || Neff_ini_pending_.has_value();
  }
```

In `species/dncdm_species.cpp`, in the constructor parsing block, relax the mutual exclusion so the shooter-written `Omega_dncdmdr` may co-occur with the initial keys **inside a shooting build**, and stash it. Replace the conflict check (`dncdm_species.cpp:162-166`) with:
```cpp
  // Omega_dncdmdr conflicts with Omega_ini/omega_ini/Neff_ini/Omega/omega — EXCEPT inside a
  // shooting build, where DoShooting writes <flavor>.Omega_dncdmdr as the initial-mode
  // fixed-point unknown alongside the user's Omega_ini/Neff_ini.
  const bool in_shooting = pfc && pfc->is_shooting;
  if (has_Omega_dncdmdr && (has_Omega_ini || has_omega_ini || has_Neff_ini || has_Omega) &&
      !in_shooting) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': Omega_dncdmdr conflicts with Omega_ini/omega_ini/Neff_ini/Omega/omega");
  }
```
(`pfc` is the `FileContent*` available in the constructor — confirm its name where `read_double` is sourced; it is the parser handle used for `input.read_double(...)`.)

Then, where the modes are resolved (`dncdm_species.cpp:177-208`), ensure that when **both** an initial key and `Omega_dncdmdr` are present (shooting iteration), the initial key still drives `deg` (via the existing `Omega_ini_pending_`/`Neff_ini_pending_` path) **and** `Omega_dncdmdr_pending_` is set from the shot value so `GetOmega0()` reserves it. Concretely, before the `if (has_deg) … else if (has_Omega_dncdmdr) … else {…}` chain, add:
```cpp
  // Initial-mode fixed-point shoot: keep the initial-abundance closure (sets deg) and also
  // record the shot Omega_dncdmdr so the composite reserves the right combined.
  if (has_Omega_dncdmdr && (has_Omega_ini || has_omega_ini || has_Neff_ini)) {
    Omega_dncdmdr_pending_ = Omega_dncdmdr_local;
    if (has_Omega_ini)        Omega_ini_pending_  = Omega_ini_local;
    else if (has_omega_ini)   Omega_ini_pending_  = omega_ini_local / settings.h / settings.h;
    else                      Neff_ini_pending_   = Neff_ini_local;
    dq_ = ComputeDq();
    return;  // skip the single-mode chain below; ApplyDncdmInitialClosure sets deg from the initial
  }
```

- [ ] **Step 4: Emit the initial-mode shooting target**

In `species/dncdm_dr_species.cpp`, extend `GetShootingTargets()` to also fire in initial mode (unknown = `Omega_dncdmdr`):
```cpp
std::vector<ShootingTarget> DNCDM_DR_Species::GetShootingTargets() const {
  if (dncdm_->Omega_dncdmdr_pending().has_value() && !dncdm_->InitialAbundanceMode()) {
    // Combined mode: shoot deg to hit the given combined.
    return {{name() + ".Omega_dncdmdr", name() + ".deg", *dncdm_->Omega_dncdmdr_pending()}};
  }
  if (dncdm_->InitialAbundanceMode()) {
    // Initial mode: fixed-point shoot on Omega_dncdmdr (closure reserve).
    // target_value is only a Newton seed; the residual is a fixed point (see below).
    double seed = dncdm_->Omega_dncdmdr_pending().value_or(0.12);
    return {{name() + ".Omega_dncdmdr_fixedpoint", name() + ".Omega_dncdmdr", seed}};
  }
  return {};
}
```

- [ ] **Step 5: Provide the initial-mode guess**

In `species/dncdm_dr_species.cpp`, extend `ComputeShootingGuess()` so initial mode seeds the `Omega_dncdmdr` unknown. The fixed-point residual is ~linear in the unknown (the unknown sets Λ only, weakly coupled to the integrated combined through H), so seed value ≈ the initial and Jacobian seed 1:
```cpp
void DNCDM_DR_Species::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                            std::vector<double>& guess,
                                            std::vector<double>& dxdy) const {
  if (dncdm_->InitialAbundanceMode()) {
    // Seed Omega_dncdmdr with the given initial abundance; Newton refines (d residual/d x ~= 1).
    double seed = 0.12;
    if (dncdm_->Omega_ini_pending().has_value())
      seed = *dncdm_->Omega_ini_pending();
    guess.push_back(seed);
    dxdy.push_back(1.0);
    return;
  }
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return;
  auto [g, d] = dncdm_->DegGuessFromOmegaToday(ctx, *dncdm_->Omega_dncdmdr_pending());
  guess.push_back(g);
  dxdy.push_back(d);
}
```

- [ ] **Step 6: Add the fixed-point residual branch**

In `species/dncdm_dr_species.cpp`, extend `ComputeShootingResidual()` to handle the fixed-point target:
```cpp
double DNCDM_DR_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                 const ShootingTarget& target) const {
  const double* bg = ctx.bg_today;
  const double H0  = ctx.pba->H0;
  const double combined = (dncdm_->Rho(bg) + dr_sp_->Rho(bg)) / (H0 * H0);
  if (target.target_name == name() + ".Omega_dncdmdr_fixedpoint") {
    // Initial mode: drive the reserved Omega_dncdmdr (= GetOmega0()) to the integrated combined.
    return -combined + GetOmega0();
  }
  // Combined mode: drive integrated combined to the given target.
  return combined - target.target_value;
}
```

- [ ] **Step 7: Build and verify the initial budget now closes**

Run:
```bash
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_initial.ini --tol 0.02 --positive-line "Neutrino Species"
```
Expected: `OK: TOTAL Omega = ~1.0`.

- [ ] **Step 8: Verify combined mode still closes (no regression)**

Run:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_combined.ini --tol 0.02 --positive-line "Neutrino Species"
```
Expected: `OK`.

- [ ] **Step 9: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp species/dncdm_dr_species.cpp test/scenarios/dncdm_dr_initial.ini
git commit -m "DNCDM_DR: initial-mode Omega_dncdmdr fixed-point shoot so closure reserves the combined"
```

---

## Task 6: Reject bare `.Omega` (today-matter) on a decaying species

`dncdm*.Omega` calls `SetOmega0` (today-matter) and leaves the budget open. Reject it with guidance.

**Files:**
- Modify: `species/dncdm_species.cpp` (the `has_Omega` handling, ~135-137)
- Create: implicit test via `assert_budget.py --expect-error`

- [ ] **Step 1: Write the failing test (a `.Omega` scenario must error)**

Create `test/scenarios/dncdm_dr_bare_omega.ini`:
```ini
T_cmb = 2.7255
omega_b = 0.022032
omega_cdm = 0.10
dncdm1.type = ncdm_decay_dr
dncdm1.m = 1.0
dncdm1.T = 0.71611
dncdm1.Gamma = 1e2
dncdm1.Omega = 0.001
```
Run:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_bare_omega.ini --expect-error "Omega_ini"
```
Expected: `FAIL` (today it runs and produces a non-closing budget, not an error).

- [ ] **Step 2: Reject `.Omega` for decaying species**

In `species/dncdm_species.cpp`, replace the `if (has_Omega) { SetOmega0(...); }` block (lines 135-137) with a hard rejection (a `ncdm_decay_dr` instance is always decaying):
```cpp
  if (has_Omega || has_omega) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': a decaying species (ncdm_decay_dr) cannot be normalized by today-matter 'Omega'/'omega' "
        "(the decay radiation would be left out of the flatness budget). Use 'Omega_ini'/'omega_ini'/"
        "'Neff_ini' to pin the initial abundance, or 'Omega_dncdmdr'/'omega_dncdmdr' to pin the "
        "combined matter+radiation density today.");
  }
```
Then remove the now-dead `has_Omega` references in the downstream conflict/co-occurrence checks (the `|| has_Omega` clauses), since `has_Omega` can no longer be true past this point.

- [ ] **Step 3: Build and verify the rejection**

Run:
```bash
make class -j 2>&1 | tail -1
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_bare_omega.ini --expect-error "Omega_ini"
```
Expected: `OK: rejected with 'Omega_ini'`.

- [ ] **Step 4: Repoint the legacy scenario off the rejected key**

Edit `test/scenarios/dncdm_dr.ini`: change `dncdm1.Omega = 0.001` to `dncdm1.Omega_ini = 0.001`. Then:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr.ini --tol 0.02
```
Expected: `OK`.

- [ ] **Step 5: Commit**

```bash
git add species/dncdm_species.cpp test/scenarios/dncdm_dr_bare_omega.ini test/scenarios/dncdm_dr.ini
git commit -m "DNCDM: reject bare .Omega (today-matter) on decaying species; repoint legacy scenario"
```

---

## Task 7: Full regression + stable-species invariance

**Files:** none (verification only)

- [ ] **Step 1: Budget closes across all three DNCDM scenarios**

Run:
```bash
for s in dncdm_dr dncdm_dr_combined dncdm_dr_initial; do
  python3 test/scenarios/assert_budget.py test/scenarios/$s.ini --tol 0.02 --positive-line "Neutrino Species" || echo "FAILED: $s"
done
python3 test/scenarios/assert_budget.py test/scenarios/dncdm_dr_bare_omega.ini --expect-error "Omega_ini" || echo "FAILED: bare_omega"
```
Expected: `OK` for all four (three closing + one rejected).

- [ ] **Step 2: DCDM_DR and combined multi-species still close**

Run:
```bash
python3 test/scenarios/assert_budget.py test/scenarios/gauge_dcdm.ini --tol 0.02
python3 test/scenarios/assert_budget.py test/scenarios/ncdm_dncdm_idmdr_combined.ini --tol 0.02 --positive-line "Neutrino Species"
```
Expected: `OK` for both (no regression to the DCDM path or the mixed sector).

- [ ] **Step 3: classy builds and a DNCDM_DR + mPk model runs (fnu path exercised)**

Run:
```bash
make classy 2>&1 | tail -3
python3 - <<'PY'
from classy import Class
c = Class()
c.set({'output':'mPk','P_k_max_1/Mpc':1.0,'non linear':'halofit',
       'omega_b':0.022032,'omega_cdm':0.10,
       'dncdm1.type':'ncdm_decay_dr','dncdm1.m':1.0,'dncdm1.T':0.71611,
       'dncdm1.Gamma':1e2,'dncdm1.Omega_ini':0.12})
c.compute()
print("OK: pk(0.1,0) =", c.pk(0.1, 0.0))
c.struct_cleanup(); c.empty()
PY
```
Expected: `OK: pk(...)` prints without error (the matter-child `Omega0_` backfill feeds `fnu`/halofit without NaN/crash).

- [ ] **Step 4: Stable-NCDM regression — bit-stable vs master**

Build a baseline `class` from `origin/master` into a temp path, run a stable-NCDM scenario through both, and diff with `compare_tol.py` (the matter-children `fnu` path must be behaviour-identical for stable NCDM):
```bash
git stash -u 2>/dev/null; git worktree add /tmp/dncdm-base origin/master 2>/dev/null || true
(cd /tmp/dncdm-base && make class -j >/dev/null 2>&1)
mkdir -p /tmp/dncdm-cmp/{base,new}
/tmp/dncdm-base/class <(sed 's#^root.*#root = /tmp/dncdm-cmp/base/o_#' test/scenarios/ncdm_dncdm_idmdr_combined.ini) >/dev/null 2>&1 || true
./class <(sed 's#^root.*#root = /tmp/dncdm-cmp/new/o_#' test/scenarios/ncdm_dncdm_idmdr_combined.ini) >/dev/null 2>&1 || true
python3 test/scenarios/compare_tol.py /tmp/dncdm-cmp/base /tmp/dncdm-cmp/new
git worktree remove /tmp/dncdm-base --force 2>/dev/null; git stash pop 2>/dev/null || true
```
Expected: `compare_tol.py` reports OK / negligible drift for the stable-NCDM observables (the combined scenario mixes stable + decaying; the stable parts must not move). If process substitution `<(...)` is awkward for `./class`, write the edited inis to temp files instead.

- [ ] **Step 5: Canonical scenario grid vs classyref**

Run the standard suite (decaying species are background-only / absent from classyref, so this guards the rest of the grid):
```bash
cd python && python test_class.py -m test_scenario 2>&1 | tail -20
```
Expected: the grid passes as before (no new failures vs the pre-existing baseline).

- [ ] **Step 6: Final commit (if any verification fixups were needed)**

```bash
git add -A && git commit -m "test: DNCDM_DR closure regression scenarios green; stable-species invariance verified" || echo "nothing to commit"
```

---

## Self-review notes (author)

- **Spec coverage:** §4.1 → Task 2; §4.2 → Task 3; §4.3 → Task 4; §4.4 combined → Task 2, initial → Task 5; §4.5 → Task 6; §6 testing → Tasks 1 & 7. Budget-print fix is already committed on the branch (spec §1).
- **Naming consistency:** `GetOmega0` (override, Task 2), `Omega_dncdmdr_pending()` / `InitialAbundanceMode()` / `Omega_ini_pending()` (Task 5), `SetOmega0` / `dncdm()` (Task 3), `target_name == name()+".Omega_dncdmdr_fixedpoint"` used identically in Tasks 4-6 of §Task 5.
- **Implementation note for the executor:** the exact parser-handle name in `dncdm_species.cpp` (called `pfc`/`input` above) and the precise location of the mode-resolution chain must be read from the file; mirror `DCDM_DR_Species` (`dcdm_dr_species.cpp:314-371,411-489`) for any ambiguity in the shooting wiring.
