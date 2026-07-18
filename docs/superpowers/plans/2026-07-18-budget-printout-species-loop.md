# Budget printout: species-owned loop — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-type string/`static_cast` branch chain in
`BackgroundModule::background_output_budget()` with a loop over `all_species_`,
fixing #377 (dcdm_wdm missing → TOTAL≠1; NCDM-family mislabeled "Neutrino
Species") and the latent Type3 gap — without adding any virtual to `BaseSpecies`.

**Architecture:** Two POD types (`BudgetBucket`, `BudgetLine`) and a free
classifier `BudgetBucketOf()` in `base_species.{h,cpp}`; a purpose-specific
`CompositeSpecies::AppendBudgetLines()` method so a composite owns its sector's
breakdown (children never exposed); a rewritten `background_output_budget()` that
loops, dispatches composite vs plain via `dynamic_cast`, reads each species' today
`Rho()/rho_crit` from the last background-table row, and prints grouped by bucket
with a local display-name map.

**Tech Stack:** C++17, CMake (primary build) + Makefile shim, `assert()`-based
unit tests run under `ctest`.

## Global Constraints

- **No new virtuals on `BaseSpecies`.** Copied verbatim from the design.
- **No general `children()` accessor.** Composites expose only the
  purpose-specific `AppendBudgetLines`.
- **CLASSpp is C++-only** — plain C++ headers, no `extern "C"`, no
  `#ifdef __cplusplus`.
- **Never `git add -A`** — stage explicit paths (in-source CMake/Xcode artifacts
  get swept in otherwise).
- **New test executables must be registered in BOTH `CMakeLists.txt` and the
  Makefile `TEST_TARGETS`** — CI runs `make test` and only builds targets named
  in `TEST_TARGETS`.
- **No bit-identical requirement.** The budget is stdout-only (no output-file
  golden). Verify with `TOTAL ≈ 1` and correct per-species labels.
- Primary build: `make class` (wraps `cmake --build build/cmake --target class`).
  The `class` binary is written to the repo root; run inis as `./class <file>`.

---

### Task 1: Budget seam — types, classifier, composite method + unit test

**Files:**
- Modify: `species/base_species.h` (add types + free-function declaration after
  the `class BaseSpecies { … };` close, currently line 693)
- Modify: `species/base_species.cpp` (define `BudgetBucketOf`)
- Modify: `species/composite_species.h` (add inline `AppendBudgetLines`)
- Create: `species/species_budget_test.cpp`
- Modify: `CMakeLists.txt` (register `test-species-budget`)
- Modify: `Makefile` (add `test-species-budget` to `TEST_TARGETS`)

**Interfaces:**
- Produces:
  - `enum class BudgetBucket { Radiation, NonRelativistic, Ncdm, Other };`
    (namespace scope, in `base_species.h`)
  - `struct BudgetLine { std::string label; double omega; BudgetBucket bucket; };`
  - `BudgetBucket BudgetBucketOf(const BaseSpecies& s);`
  - `void CompositeSpecies::AppendBudgetLines(const double* pvecback_today, double rho_crit, std::vector<BudgetLine>& out) const;`

- [ ] **Step 1: Write the failing test**

Create `species/species_budget_test.cpp`:

```cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "base_species.h"
#include "composite_species.h"

// Minimal concrete species: configurable name + energy_type, preset Rho.
// Implements only the 7 pure virtuals of BaseSpecies; Rho ignores pvecback.
class BudgetFake : public BaseSpecies {
 public:
  BudgetFake(std::string name, EnergyType et, double rho)
      : BaseSpecies(std::move(name), et), rho_(rho) {}
  void RegisterBackgroundIndices(int&) override {}
  void ComputeBackground(double, const double*, double*) override {}
  double Rho(const double*) const override { return rho_; }
  double P(const double*) const override { return 0.; }
  void PerturbDerivs(const PerturbLayout&, double, const double*, double*,
                     const perturb_parameters_and_workspace&) const override {}
  StressEnergyContribution StressEnergy(const PerturbLayout&, const perturb_vector*,
                                        const double*, const double*,
                                        const perturb_workspace*) const override {
    return {};
  }
  double GetOmega0() const override { return 0.; }

 private:
  double rho_;
};

// CompositeSpecies::children_ is protected; this exposes AddChild for the test.
class TestComposite : public CompositeSpecies {
 public:
  TestComposite() : CompositeSpecies("test_composite", EnergyType::Other) {}
  void AddChild(std::unique_ptr<BaseSpecies> c) { children_.push_back(std::move(c)); }
};

int main() {
  // 1) BudgetBucketOf maps energy_type for non-NCDM species.
  BudgetFake rad("g", BaseSpecies::EnergyType::Radiation, 1.0);
  BudgetFake mat("cdm", BaseSpecies::EnergyType::Matter, 1.0);
  BudgetFake de("lambda", BaseSpecies::EnergyType::DarkEnergy, 1.0);
  BudgetFake oth("scf", BaseSpecies::EnergyType::Other, 1.0);
  assert(BudgetBucketOf(rad) == BudgetBucket::Radiation);
  assert(BudgetBucketOf(mat) == BudgetBucket::NonRelativistic);
  assert(BudgetBucketOf(de) == BudgetBucket::Other);
  assert(BudgetBucketOf(oth) == BudgetBucket::Other);

  // 2) AppendBudgetLines: one line per child, omega = Rho/rho_crit,
  //    label = child name(), bucket = BudgetBucketOf(child). (Rho ignores the
  //    null pvecback here.)
  TestComposite comp;
  comp.AddChild(std::make_unique<BudgetFake>("DCDM", BaseSpecies::EnergyType::Matter, 2.0));
  comp.AddChild(std::make_unique<BudgetFake>("DR", BaseSpecies::EnergyType::Radiation, 3.0));
  std::vector<BudgetLine> lines;
  comp.AppendBudgetLines(/*pvecback_today=*/nullptr, /*rho_crit=*/10.0, lines);
  assert(lines.size() == 2);
  assert(lines[0].label == "DCDM");
  assert(lines[0].bucket == BudgetBucket::NonRelativistic);
  assert(std::fabs(lines[0].omega - 0.2) < 1e-12);
  assert(lines[1].label == "DR");
  assert(lines[1].bucket == BudgetBucket::Radiation);
  assert(std::fabs(lines[1].omega - 0.3) < 1e-12);

  printf("species_budget_test: all assertions passed\n");
  return 0;
}
```

Register it. In `CMakeLists.txt`, after the line
`add_executable(test-ncdm-family species/ncdm_family_test.cpp)` (currently 235),
add:

```cmake
  add_executable(test-species-budget species/species_budget_test.cpp)
```

and add `test-species-budget` to the `foreach(_t IN ITEMS …)` list on the next
line (append it before the closing `)`).

In `Makefile`, change the end of the `TEST_TARGETS` block so the last two lines read:

```make
	test-dncdm-switch-copy \
	test-ncdm-family \
	test-species-budget
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
make test-species-budget
```
Expected: FAIL at compile — `BudgetBucket`, `BudgetLine`, `BudgetBucketOf`, and
`CompositeSpecies::AppendBudgetLines` are undeclared.

- [ ] **Step 3: Add the POD types + classifier declaration**

In `species/base_species.h`, immediately after the closing `};` of
`class BaseSpecies` (line 693, the one preceded by `friend class SpeciesCollection;`),
add:

```cpp

// ── Background budget printout (background_verbose > 1) ───────────────────────
// Coarse presentation buckets for background_output_budget(). Distinct from
// EnergyType because NCDM-family species are EnergyType::Other yet warrant their
// own section, and composites (also EnergyType::Other) split across buckets by
// child.
enum class BudgetBucket { Radiation, NonRelativistic, Ncdm, Other };

// One printed budget line: a (sub-)species' today density fraction and its bucket.
struct BudgetLine {
  std::string label;
  double omega;
  BudgetBucket bucket;
};

// Classify a species into a budget bucket. NCDM-family -> Ncdm (the same family
// test GetNcdmSpecies uses); otherwise mapped from energy_type(). Defined in
// base_species.cpp.
BudgetBucket BudgetBucketOf(const BaseSpecies& s);
```

- [ ] **Step 4: Define the classifier**

Replace the entire contents of `species/base_species.cpp` (currently just
`#include "base_species.h"`) with:

```cpp
#include "base_species.h"

#include "ncdm_base_species.h"

BudgetBucket BudgetBucketOf(const BaseSpecies& s) {
  // NCDM-family species carry EnergyType::Other (they tally into both rho_r and
  // rho_m), so route them by family, not energy_type. Same dynamic_cast idiom as
  // GetNcdmSpecies.
  if (dynamic_cast<const NCDMBaseSpecies*>(&s))
    return BudgetBucket::Ncdm;
  switch (s.energy_type()) {
    case BaseSpecies::EnergyType::Radiation:
      return BudgetBucket::Radiation;
    case BaseSpecies::EnergyType::Matter:
      return BudgetBucket::NonRelativistic;
    case BaseSpecies::EnergyType::DarkEnergy:
    case BaseSpecies::EnergyType::Other:
      return BudgetBucket::Other;
  }
  return BudgetBucket::Other;  // unreachable; silences -Wreturn-type
}
```

- [ ] **Step 5: Add the composite method**

In `species/composite_species.h`, inside `class CompositeSpecies`, in the
`// ── Omega0 (closure) ─` region (right after the closing `}` of the
`GetOmega0()` override, currently near line 93), add:

```cpp
  /** Append one budget line per child (label = child name(), omega = today
   *  Rho/rho_crit, bucket = BudgetBucketOf(child)) for background_output_budget().
   *  The composite owns its sector breakdown; children_ is never exposed. */
  void AppendBudgetLines(const double* pvecback_today,
                         double rho_crit,
                         std::vector<BudgetLine>& out) const {
    for (const auto& c : children_)
      out.push_back({c->name(), c->Rho(pvecback_today) / rho_crit, BudgetBucketOf(*c)});
  }
```

(`base_species.h` — which declares `BudgetLine`/`BudgetBucketOf` — and `<vector>`
are already included by `composite_species.h`.)

- [ ] **Step 6: Run test to verify it passes**

Run:
```bash
make test-species-budget && ctest --test-dir build/cmake -R test-species-budget --output-on-failure
```
Expected: PASS — `species_budget_test: all assertions passed`, ctest reports
`1/1 Test #…: test-species-budget … Passed`.

- [ ] **Step 7: Commit**

```bash
git add species/base_species.h species/base_species.cpp species/composite_species.h \
        species/species_budget_test.cpp CMakeLists.txt Makefile
git commit -m "Budget seam: BudgetBucket/BudgetLine + BudgetBucketOf + CompositeSpecies::AppendBudgetLines (#377)"
```

---

### Task 2: Rewrite `background_output_budget()` to loop over species

**Files:**
- Modify: `source/background_module.cpp` (add includes near top; replace the whole
  `background_output_budget()` function, currently ~lines 1123–1225)

**Interfaces:**
- Consumes (from Task 1): `BudgetBucket`, `BudgetLine`, `BudgetBucketOf`,
  `CompositeSpecies::AppendBudgetLines`.

- [ ] **Step 1: Add includes**

In `source/background_module.cpp`, in the `#include` block near the top: after
`#include <algorithm>` (line 8) add:

```cpp
#include <map>
#include <vector>
```

and in the species-header group (after `#include "../species/bisection.h"` is not
right — put it with the species headers, e.g. after
`#include "../species/type3_species.h"`) add:

```cpp
#include "../species/composite_species.h"
```

- [ ] **Step 2: Replace the function**

Replace the entire `background_output_budget()` function — from
`void BackgroundModule::background_output_budget() {` through its matching closing
`}` — with:

```cpp
void BackgroundModule::background_output_budget() {
  if (pba->background_verbose <= 1)
    return;

  // Today densities live in the last background-table row; Rho()/rho_crit is a
  // species' Omega today and closes the budget to 1 by construction.
  const double* bg_today = background_table_.data() + (bt_size_ - 1) * bg_size_;
  const double rho_crit  = bg_today[index_bg_rho_crit_];

  // Display-only: terse internal species names -> readable labels. Any name not
  // in the map prints verbatim via name(). This map never drives logic.
  static const std::map<std::string, const char*> kPretty = {
      {"Baryons", "Baryons"},
      {"Photons", "Photons"},
      {"CDM", "Cold Dark Matter"},
      {"UR", "Ultra-relativistic relics"},
      {"Lambda", "Cosmological Constant"},
      {"Fluid", "Dark Energy Fluid"},
      {"ScalarField", "Scalar Field"},
      {"DCDM", "Decaying Cold Dark Matter"},
      {"DR", "Dark Radiation"},
      {"IDR", "Interacting Dark Radiation"},
      {"IDM_DR", "Interacting Dark Matter (DR)"},
      {"IDR_DRMD", "Dark Radiation (DRMD)"},
      {"IDM_DRMD", "Interacting Dark Matter (DRMD)"},
  };
  auto pretty = [&](const std::string& name) -> const char* {
    auto it = kPretty.find(name);
    return it != kPretty.end() ? it->second : name.c_str();
  };

  // Gather one line per (sub-)species. Composites own their sector breakdown;
  // plain species contribute a single line. Curvature is not a species.
  std::vector<BudgetLine> lines;
  for (const auto& [key, sp] : all_species_) {
    if (const auto* comp = dynamic_cast<const CompositeSpecies*>(sp.get()))
      comp->AppendBudgetLines(bg_today, rho_crit, lines);
    else
      lines.push_back({sp->name(), sp->Rho(bg_today) / rho_crit, BudgetBucketOf(*sp)});
  }
  if (pba->sgnK != 0)
    lines.push_back({"Spatial Curvature", pba->Omega0_k, BudgetBucket::Other});

  struct Section {
    BudgetBucket bucket;
    const char* header;
    const char* total_label;
  };
  static const Section kSections[] = {
      {BudgetBucket::Radiation, " ---> Relativistic Species ", "Radiation"},
      {BudgetBucket::NonRelativistic, " ---> Nonrelativistic Species ", "Non-relativistic"},
      {BudgetBucket::Ncdm, " ---> Non-Cold Dark Matter (NCDM) ", "Non-Cold Dark Matter"},
      {BudgetBucket::Other, " ---> Other Content ", "Other Content"},
  };

  printf(" ---------------------------- Budget equation ----------------------- \n");

  double totals[4] = {0., 0., 0., 0.};
  for (const auto& section : kSections) {
    bool any = false;
    for (const auto& l : lines)
      if (l.bucket == section.bucket) {
        any = true;
        break;
      }
    if (!any)
      continue;
    printf("%s\n", section.header);
    for (const auto& l : lines) {
      if (l.bucket != section.bucket)
        continue;
      printf("-> %-30s Omega = %-15g , omega = %-15g\n",
             pretty(l.label), l.omega, l.omega * pba->h * pba->h);
      totals[static_cast<int>(section.bucket)] += l.omega;
    }
  }

  printf(" ---> Total budgets \n");
  double grand = 0.;
  for (const auto& section : kSections) {
    const double t = totals[static_cast<int>(section.bucket)];
    grand += t;
    if (t != 0.)
      printf(" %-32s Omega = %-15g , omega = %-15g \n",
             section.total_label, t, t * pba->h * pba->h);
  }
  printf(" TOTAL                            Omega = %-15g , omega = %-15g \n",
         grand, grand * pba->h * pba->h);
  printf(" -------------------------------------------------------------------- \n");
}
```

- [ ] **Step 3: Build**

Run:
```bash
make class
```
Expected: builds to `./class` with no errors and no new warnings on
`background_module.cpp`.

- [ ] **Step 4: Verify defect 1 — dcdm_wdm closes the budget**

Run:
```bash
S="/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/e76b6160-fb6e-4f69-9b74-a4adc2719137/scratchpad"
cp test/dotsyntax_dcdm_wdm.ini "$S/bw.ini"
grep -q '^background_verbose' "$S/bw.ini" \
  && sed -i '' 's/^background_verbose.*/background_verbose = 2/' "$S/bw.ini" \
  || printf '\nbackground_verbose = 2\n' >> "$S/bw.ini"
./class "$S/bw.ini" | sed -n '/Budget equation/,/-----------------------------------/p'
```
Expected: a `Non-Cold Dark Matter (NCDM)` section listing the dcdm_wdm daughter,
a `Nonrelativistic` line for the `DCDM` parent, and `TOTAL … Omega = ~1.0`
(the pre-fix value was `0.974`).

- [ ] **Step 5: Verify defect 2 — NCDM-family self-labels**

Run:
```bash
S="/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/e76b6160-fb6e-4f69-9b74-a4adc2719137/scratchpad"
for f in dotsyntax_ncdm dotsyntax_axion dotsyntax_dncdm; do
  cp "test/$f.ini" "$S/$f.ini"
  grep -q '^background_verbose' "$S/$f.ini" \
    && sed -i '' 's/^background_verbose.*/background_verbose = 2/' "$S/$f.ini" \
    || printf '\nbackground_verbose = 2\n' >> "$S/$f.ini"
  echo "===== $f ====="; ./class "$S/$f.ini" | sed -n '/Budget equation/,/-----------------------------------/p'
done
```
Expected: **no** line says `Neutrino Species`; each NCDM-family species prints its
own name; the `dncdm` case shows the dncdm child under NCDM and its dark radiation
under Radiation; each `TOTAL ≈ 1.0`.

- [ ] **Step 6: Verify no regression — plain ΛCDM still closes and reads well**

Run:
```bash
S="/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/e76b6160-fb6e-4f69-9b74-a4adc2719137/scratchpad"
cp explanatory.ini "$S/exp.ini"
grep -q '^background_verbose' "$S/exp.ini" \
  && sed -i '' 's/^background_verbose.*/background_verbose = 2/' "$S/exp.ini" \
  || printf '\nbackground_verbose = 2\n' >> "$S/exp.ini"
./class "$S/exp.ini" | sed -n '/Budget equation/,/-----------------------------------/p'
```
Expected: `Photons`, `Ultra-relativistic relics` under Relativistic; `Baryons`,
`Cold Dark Matter` under Nonrelativistic; `Cosmological Constant` under Other;
`TOTAL … Omega = ~1.0`.

- [ ] **Step 7: Commit**

```bash
git add source/background_module.cpp
git commit -m "background_output_budget: loop over species instead of per-type branches (fixes #377)"
```

---

### Task 3: Remove the now-dead `PrintOmegaInfo`

**Files:**
- Modify: `species/ncdm_base_species.h` (remove declaration, line 137)
- Modify: `species/ncdm_base_species.cpp` (remove definition, lines ~636–638)

**Interfaces:** none produced; removes the `"Neutrino Species"` hardcode whose
only caller (the old budget loop) was deleted in Task 2.

- [ ] **Step 1: Confirm it is dead**

Run:
```bash
grep -rn "PrintOmegaInfo" source/ species/ python/ include/
```
Expected: exactly two hits — the declaration in `ncdm_base_species.h` and the
definition in `ncdm_base_species.cpp`. (If any other caller appears, stop and
report; do not remove.)

- [ ] **Step 2: Remove the declaration**

In `species/ncdm_base_species.h`, delete the line:

```cpp
  void PrintOmegaInfo() const;
```

- [ ] **Step 3: Remove the definition**

In `species/ncdm_base_species.cpp`, delete the function:

```cpp
void NCDMBaseSpecies::PrintOmegaInfo() const {
  printf("-> %-26s Omega = %-15g , omega = %-15g\n", "Neutrino Species", Omega0_, omega0_);
}
```

- [ ] **Step 4: Rebuild and re-verify**

Run:
```bash
make class && grep -rn "PrintOmegaInfo" source/ species/ python/ include/ || echo "PrintOmegaInfo fully removed"
S="/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/e76b6160-fb6e-4f69-9b74-a4adc2719137/scratchpad"
./class "$S/dotsyntax_ncdm.ini" | sed -n '/Budget equation/,/-----------------------------------/p'
```
Expected: builds cleanly; `PrintOmegaInfo fully removed`; NCDM budget still prints
each species' own name with `TOTAL ≈ 1.0`.

- [ ] **Step 5: Run the full unit-test suite**

Run:
```bash
make test
```
Expected: all tests pass, including `test-species-budget`.

- [ ] **Step 6: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp
git commit -m "Remove dead NCDMBaseSpecies::PrintOmegaInfo (superseded by the budget loop, #377)"
```

---

## Self-Review

**Spec coverage:**
- Design §1 (POD types + `BudgetBucketOf`) → Task 1 Steps 3–4. ✓
- Design §2 (`CompositeSpecies::AppendBudgetLines`) → Task 1 Step 5. ✓
- Design §3 (rewrite `background_output_budget`, loop + pretty map + grouped
  print + curvature line + NCDM rename) → Task 2. ✓
- Design §4 (remove dead `PrintOmegaInfo`) → Task 3. ✓
- Defect 1 (dcdm_wdm/TOTAL) → Task 2 Step 4. ✓
- Defect 2 (labels) → Task 2 Step 5 + Task 3. ✓
- Type3 latent gap → covered by the generic composite loop (Task 2); not
  separately asserted because `test/` has no Type3 ini (acceptable per design
  non-goals).
- Include pruning: intentionally **not** done — several composite headers in
  `background_module.cpp` are still referenced by `GetOmega0Species` and the
  `verbose>2` detail block, so leaving all includes is the safe, focused choice.

**Placeholder scan:** none — every code step contains complete code.

**Type consistency:** `BudgetBucket`, `BudgetLine`, `BudgetBucketOf`, and
`AppendBudgetLines(const double*, double, std::vector<BudgetLine>&)` are used
identically in Task 1 (definition + test) and Task 2 (consumption). Enum order
`{Radiation, NonRelativistic, Ncdm, Other}` matches the `totals[4]` indexing and
the `kSections` order in Task 2.
