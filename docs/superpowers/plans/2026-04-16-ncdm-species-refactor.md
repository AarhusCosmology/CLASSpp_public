# NCDM Species Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Dissolve `NonColdDarkMatter` and `DarkRadiation` tool classes into self-contained species classes with `CreateAll` factories, making each species fully responsible for reading its own parameters from `FileContent`.

**Architecture:** Introduce `NCDMBaseSpecies` (inherits `BaseSpecies`) owning all single-species quadrature/distribution-function data; `NCDMSpecies` and `DNCDMSpecies` inherit from it. `DarkRadiationSpecies` absorbs `DarkRadiation` directly. Each concrete species class gets a `CreateAll` (or `Create`) static factory that reads its own parameters from `FileContent*`. `InputModule` lines 1022–1193 collapse to three `CreateAll` calls.

**Tech Stack:** C++17, CLASS++ build system (`make`), Python/classy integration tests (`cd python && python -m pytest test_class.py -v`)

**Spec:** `docs/superpowers/specs/2026-04-16-ncdm-species-refactor-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| CREATE | `species/ncdm_base_species.h` | `NCDMBaseSpecies` declaration + `NcdmSettings` struct |
| CREATE | `species/ncdm_base_species.cpp` | All single-species quadrature, dist. function, momenta |
| MODIFY | `species/ncdm_species.h` | Inherit `NCDMBaseSpecies`; add `CreateAll` |
| MODIFY | `species/ncdm_species.cpp` | Remove `ncdm_`/`ncdm_id_`; call base methods |
| MODIFY | `species/dncdm_species.h` | Inherit `NCDMBaseSpecies`; add decay data; `CreateAll` |
| MODIFY | `species/dncdm_species.cpp` | Absorb `DecayDRProperties`; `CreateAll` reads input |
| MODIFY | `species/dncdm_dr_species.h` | Add `CreateAll` factory |
| MODIFY | `species/dncdm_dr_species.cpp` | `CreateAll` wraps `DNCDMSpecies::CreateAll` |
| MODIFY | `species/dark_radiation_species.h` | Absorb `DarkRadiation` data; add `Create` |
| MODIFY | `species/dark_radiation_species.cpp` | Single-channel quadrature + integration |
| MODIFY | `source/input_module.h` | Remove `ncdm_`, `dr_` members; update `target_names` |
| MODIFY | `source/input_module.cpp` | Replace lines 1022–1193 with `CreateAll` calls |
| MODIFY | `source/base_module.h` | Remove `ncdm_`, `dr_` members and includes |
| MODIFY | `source/background_module.cpp` | Replace 5 `ncdm_->*` calls with species iteration |
| MODIFY | `source/perturbations_module.cpp` | Replace `ncdm_->q_size_ncdm_[n]` with `sp->q_size()` |
| MODIFY | `source/nonlinear_module.cpp` | Replace `ncdm_->GetMassInElectronvolt(id)` with `sp->GetMassInElectronvolt()` |
| MODIFY | `Makefile` | Remove tool `.opp` targets; add `ncdm_base_species` |
| DELETE | `tools/non_cold_dark_matter.h/.cpp` | Gone |
| DELETE | `tools/dark_radiation.h/.cpp` | Gone |

---

## Task 1: Create `NCDMBaseSpecies` header

**Files:**
- Create: `species/ncdm_base_species.h`

- [ ] **Step 1: Write the header**

```cpp
// species/ncdm_base_species.h
#pragma once
#include <memory>
#include <tuple>
#include <vector>

#include "../species/base_species.h"
#include "background.h"
#include "parser.h"

class BackgroundModule;

struct NcdmSettings {
  double h;
  double T_cmb;
  double tol_ncdm;
  double tol_ncdm_bg;
  double tol_M_ncdm;
};

/**
 * Abstract base for all NCDM flavors. Owns per-species quadrature,
 * distribution function, and thermodynamic parameters. Subclasses
 * implement the perturbation interface.
 */
class NCDMBaseSpecies : public BaseSpecies {
 public:
  // ── Public accessors ──────────────────────────────────────────────────────
  double GetOmega0() const;
  double GetNeff(double z) const;
  double GetMassInElectronvolt() const { return m_in_eV_; }
  double GetDeg() const { return deg_; }
  double GetIni(double a, double a_today, double tol_ncdm_initial_w) const;
  double GetRescalingFactor(const double* pvecback_begin) const;
  std::tuple<double, double> GetRescaledParameters(double a, const double* lnf_array) const;

  // Called from NCDMSpecies/DNCDMSpecies background methods
  int ComputeMomenta(double z,
                     double* n,
                     double* rho,
                     double* p,
                     double* drho_dM,
                     double* pseudo_p) const;

  int q_size() const { return static_cast<int>(q_.size()); }
  int q_size_bg() const { return static_cast<int>(q_bg_.size()); }

  void PrintNeffInfo() const;
  void PrintMassInfo() const;
  void PrintOmegaInfo() const;

  // Background overrides (shared by NCDMSpecies and DNCDMSpecies)
  void SetBackgroundModule(const BackgroundModule* bgm) override { bgm_ = bgm; }

 protected:
  NCDMBaseSpecies(std::string name,
                  EnergyType energy_type,
                  FileContent* pfc,
                  int species_index,
                  const NcdmSettings& settings);

  // ── Quadrature — perturbation sampling ───────────────────────────────────
  std::vector<double> q_;
  std::vector<double> w_;
  std::vector<double> dlnf0_dlnq_;

  // ── Quadrature — background sampling ─────────────────────────────────────
  std::vector<double> q_bg_;
  std::vector<double> w_bg_;

  double factor_ = 0.;

  // ── Species thermodynamic parameters ─────────────────────────────────────
  double m_in_eV_ = 0.;
  double M_       = 0.;   // dimensionless mass: m / (k_B * T_ncdm)
  double deg_     = 1.;
  double T_       = 0.71611;   // T_ncdm / T_cmb
  double ksi_     = 0.;        // mu / T_ncdm

  const BackgroundModule* bgm_ = nullptr;
  double T_cmb_    = 0.;
  double h_        = 0.;

 private:
  void ReadParameters(FileContent* pfc, int species_index, const NcdmSettings& settings);
  void InitQuadrature(const NcdmSettings& settings);
  void InitDistribution(FileContent* pfc, int species_index);

  int ComputeMomentaMass(double M,
                         double z,
                         double* n,
                         double* rho,
                         double* p,
                         double* drho_dM,
                         double* pseudo_p) const;

  static int DistributionFunction(void* params, double q, double* f0);
  static int TestFunction(void* params, double q, double* test);

  int quadrature_strategy_ = 3;
  int input_q_size_        = -1;
  double qmax_             = 15.;
  std::vector<double> psd_parameters_;
  bool got_file_           = false;
  std::string psd_file_;

  mutable ErrorMsg error_message_;
};
```

- [ ] **Step 2: Create stub `.cpp` that compiles**

```cpp
// species/ncdm_base_species.cpp
#include "ncdm_base_species.h"

NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 int species_index,
                                 const NcdmSettings& settings)
    : BaseSpecies(std::move(name), energy_type),
      T_cmb_(settings.T_cmb),
      h_(settings.h) {
  ReadParameters(pfc, species_index, settings);
}

void NCDMBaseSpecies::ReadParameters(FileContent* /*pfc*/, int /*species_index*/,
                                     const NcdmSettings& /*settings*/) {
  // TODO: implement in Task 2
}

void NCDMBaseSpecies::InitQuadrature(const NcdmSettings& /*settings*/) {}
void NCDMBaseSpecies::InitDistribution(FileContent* /*pfc*/, int /*species_index*/) {}

int NCDMBaseSpecies::ComputeMomenta(double, double*, double*, double*, double*, double*) const {
  return _FAILURE_;
}
double NCDMBaseSpecies::GetOmega0() const { return 0.; }
double NCDMBaseSpecies::GetNeff(double) const { return 0.; }
double NCDMBaseSpecies::GetIni(double a, double, double) const { return a; }
double NCDMBaseSpecies::GetRescalingFactor(const double*) const { return 1.; }
std::tuple<double,double> NCDMBaseSpecies::GetRescaledParameters(double, const double*) const { return {0.,0.}; }
void NCDMBaseSpecies::PrintNeffInfo() const {}
void NCDMBaseSpecies::PrintMassInfo() const {}
void NCDMBaseSpecies::PrintOmegaInfo() const {}
int NCDMBaseSpecies::ComputeMomentaMass(double,double,double*,double*,double*,double*,double*) const { return _FAILURE_; }
int NCDMBaseSpecies::DistributionFunction(void*, double, double*) { return _FAILURE_; }
int NCDMBaseSpecies::TestFunction(void*, double, double*) { return _FAILURE_; }
```

- [ ] **Step 3: Add to `Makefile`**

Find the block that lists `species/*.cpp` objects and add `ncdm_base_species`:

```makefile
# In the OFILES_CPP block (search for "ncdm_species.opp"):
OFILES_CPP += $(WRKDIR)/ncdm_base_species.opp
```

- [ ] **Step 4: Compile to confirm the stub builds**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: compiles without errors (stubs return `_FAILURE_` but link fine).

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp Makefile
git commit -m "feat: add NCDMBaseSpecies skeleton"
```

---

## Task 2: Implement `NCDMBaseSpecies` — parameter reading, quadrature, momenta

**Files:**
- Modify: `species/ncdm_base_species.h` (minor additions)
- Modify: `species/ncdm_base_species.cpp` (full implementation)

The implementation is a single-species extract of `NonColdDarkMatter`. The source of truth is `tools/non_cold_dark_matter.cpp` — read it alongside this task.

- [ ] **Step 1: Implement `ReadParameters`**

Copy `NonColdDarkMatter::background_ncdm_init` from `tools/non_cold_dark_matter.cpp` (the section that reads `m_ncdm`, `T_ncdm`, `deg_ncdm`, `ksi_ncdm`, `ncdm_psd_parameters`, file-based PSD, quadrature settings for a single species at `species_index`). Remove all loops over `N_ncdm_` — there is only one species per instance. The per-species indexed reads like `m_ncdm__1` become reads for index `species_index+1`.

Key reads to perform (using `parser_read_list_of_doubles`):
- `m_ncdm` → `m_in_eV_`
- `T_ncdm` → `T_`
- `deg_ncdm` → `deg_`
- `ksi_ncdm` → `ksi_`
- `quadrature_strategy_ncdm` → `quadrature_strategy_`
- `N_momentum_bins_ncdm` → `input_q_size_`
- `maximum_q_ncdm` → `qmax_`
- `ncdm_psd_parameters` → `psd_parameters_`
- `ncdm_psd_files` → `psd_file_` / `got_file_`

After reading, call:
```cpp
InitQuadrature(settings);
InitDistribution(pfc, species_index);
// Set M_ from m_in_eV_ and T_ncdm:
M_ = m_in_eV_ / (_k_B_ / _eV_) / (T_ * T_cmb_);
// Set factor_:
factor_ = deg_ * 4. * _PI_ / pow(2. * _PI_ * _hbar_, 3);
```

- [ ] **Step 2: Implement `InitQuadrature`**

Copy `background_ncdm_init` quadrature section from `non_cold_dark_matter.cpp`. The key difference: call `get_qsampling_trapezoidal` or the relevant quadrature routine once, storing results into `q_`, `w_`, `q_bg_`, `w_bg_`, `dlnf0_dlnq_` (all single vectors, not vector-of-vectors).

- [ ] **Step 3: Implement `ComputeMomentaMass` and `ComputeMomenta`**

Copy `background_ncdm_momenta_mass` from `non_cold_dark_matter.cpp`. Rename to `ComputeMomentaMass`. Remove the `n_ncdm` index — use `q_bg_`, `w_bg_`, `M_` directly:

```cpp
int NCDMBaseSpecies::ComputeMomentaMass(double M, double z, double* n, double* rho,
                                        double* p, double* drho_dM, double* pseudo_p) const {
  // Same body as background_ncdm_momenta_mass but using q_bg_, w_bg_ directly
  // (no indexing by n_ncdm)
  ...
}

int NCDMBaseSpecies::ComputeMomenta(double z, double* n, double* rho,
                                    double* p, double* drho_dM, double* pseudo_p) const {
  return ComputeMomentaMass(M_, z, n, rho, p, drho_dM, pseudo_p);
}
```

- [ ] **Step 4: Implement `GetOmega0`, `GetNeff`, `GetIni`**

```cpp
double NCDMBaseSpecies::GetOmega0() const {
  double n, rho, p, drho_dM, pseudo_p;
  ComputeMomenta(0., &n, &rho, &p, &drho_dM, &pseudo_p);
  return rho / (3. * h_ * h_);  // same formula as NonColdDarkMatter::GetOmega0 for one species
}

double NCDMBaseSpecies::GetNeff(double z) const {
  double n, rho, p, drho_dM, pseudo_p;
  ComputeMomenta(0., &n, &rho, &p, &drho_dM, &pseudo_p);
  // Same Neff formula from non_cold_dark_matter.cpp GetNeff, for one species
  ...
}

double NCDMBaseSpecies::GetIni(double a, double a_today, double tol) const {
  // Single-species version of NonColdDarkMatter::GetIni:
  // shrink a until w < tol, or return a unchanged
  ...
}
```

Copy the implementations from `non_cold_dark_matter.cpp` (methods `GetOmega0`, `GetNeff`, `GetIni`), removing the loop over `N_ncdm_`.

- [ ] **Step 5: Implement `GetRescalingFactor`, `GetRescaledParameters`**

Copy directly from `non_cold_dark_matter.cpp` replacing `n_ncdm` indexing with direct member access.

- [ ] **Step 6: Implement `SetOmega0`, `SetDegAndFactor`, `SetDeg_from_Omega_ini` as private helpers**

These are called from `CreateAll` factories (Tasks 3 and 4). Declare them `protected` in the header:

```cpp
// In ncdm_base_species.h (protected section):
void SetOmega0(double Omega0, double h);
void SetDegAndFactor(double deg);
void SetDeg_from_Omega_ini(double z_ini, double H0, double Omega_ini);
```

Copy implementations from `non_cold_dark_matter.cpp`, stripping the `n_ncdm` index.

- [ ] **Step 7: Implement `PrintNeffInfo`, `PrintMassInfo`, `PrintOmegaInfo`**

Per-species print helpers. Copy from `non_cold_dark_matter.cpp` `PrintNeffInfo` etc., removing the loop. Each prints one line for its own species.

- [ ] **Step 8: Implement `DistributionFunction` and `TestFunction` static callbacks**

These are the `background_ncdm_distribution` and `background_ncdm_test_function` statics. Copy from `non_cold_dark_matter.cpp` replacing the `background_parameters_for_distributions` struct usage. Add a private nested struct `DistributionParams` to `NCDMBaseSpecies`:

```cpp
// Private in ncdm_base_species.h:
struct DistributionParams {
  const NCDMBaseSpecies* sp;
  // interpolation state for file-based PSD:
  int tablesize = 0;
  std::vector<double> q, f0, d2f0;
  int last_index = 0;
};
```

- [ ] **Step 9: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: clean compile. The existing `NCDMSpecies`/`DNCDMSpecies` still use `NonColdDarkMatter` — no regression yet.

- [ ] **Step 10: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp
git commit -m "feat: implement NCDMBaseSpecies — quadrature, distribution function, momenta"
```

---

## Task 3: Migrate `NCDMSpecies` to inherit from `NCDMBaseSpecies`

**Files:**
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp`

- [ ] **Step 1: Rewrite `ncdm_species.h`**

```cpp
// species/ncdm_species.h
#pragma once
#include <memory>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

class NCDMSpecies : public NCDMBaseSpecies {
 public:
  NCDMSpecies(FileContent* pfc, int species_index, const NcdmSettings& settings,
              const background* pba, const BackgroundModule* bgm);

  static std::vector<std::unique_ptr<NCDMSpecies>>
      CreateAll(FileContent* pfc, const NcdmSettings& settings,
                const background* pba, const BackgroundModule* bgm);

  bool IsFreestreaming() const override { return true; }

  // ── Background ──────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override { return pvecback[index_bg_rho_]; }
  double P(const double* pvecback) const override { return pvecback[index_bg_p_]; }
  double DpDloga(const double* pvecback) const override {
    return pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_];
  }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr,
                                   int& index_pt, const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  double Delta(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double Theta(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double DeltaP(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double RhoPlusPShear(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void WriteOutputColumns(PerturbColumnWriter&, const PerturbationsModule&,
                          enum file_format, TransferColumnSection) const override;
  void PrintVariables(PerturbColumnWriter&, double, const double*,
                      const PerturbationsModule&, const perturb_workspace*) const override;

  // Accessors used by PerturbationsModule for tensor modes and indexing:
  int bg_number_index() const { return index_bg_number_; }
  int bg_pseudo_p_index() const { return index_bg_pseudo_p_; }

 private:
  const background* pba_;

  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;
  int index_pt_psi0_     = -1;
};
```

Note: `ncdm_id_` is gone — it was only needed to index into `NonColdDarkMatter`. Order in `all_species_` is now determined by insertion order from `CreateAll`.

- [ ] **Step 2: Update `ncdm_species.cpp` — constructor**

```cpp
NCDMSpecies::NCDMSpecies(FileContent* pfc, int species_index,
                         const NcdmSettings& settings,
                         const background* pba, const BackgroundModule* bgm)
    : NCDMBaseSpecies("NCDM_" + std::to_string(species_index),
                      EnergyType::Other, pfc, species_index, settings),
      pba_(pba) {
  bgm_ = bgm;
}
```

- [ ] **Step 3: Implement `NCDMSpecies::CreateAll`**

```cpp
std::vector<std::unique_ptr<NCDMSpecies>>
NCDMSpecies::CreateAll(FileContent* pfc, const NcdmSettings& settings,
                       const background* pba, const BackgroundModule* bgm) {
  std::vector<std::unique_ptr<NCDMSpecies>> result;

  // Read the type list to determine which indices are standard NCDM
  // (vs decay_dr which belong to DNCDMSpecies)
  // Read N_ncdm_standard from pfc — copy the type-discrimination logic
  // from NonColdDarkMatter::background_ncdm_init:
  int N_ncdm = 0;
  int flag1;
  char errmsg[2048];
  parser_read_int(pfc, "N_ncdm", &N_ncdm, &flag1, errmsg);

  // Read type list (standard/decay_dr per entry):
  std::vector<std::string> type_list(N_ncdm, "standard");
  // (Read "type_ncdm" comma-separated list if present, same as in
  // NonColdDarkMatter::background_ncdm_init)

  int standard_index = 0;
  for (int n = 0; n < N_ncdm; ++n) {
    if (type_list[n] != "decay_dr") {
      result.push_back(std::make_unique<NCDMSpecies>(
          pfc, n, settings, pba, bgm));
      standard_index++;
    }
  }
  return result;
}
```

> **Note:** The exact type-discrimination parsing (reading `type_ncdm` comma-separated list) is in `NonColdDarkMatter::background_ncdm_init`. Copy that logic directly into `CreateAll`.

- [ ] **Step 4: Update all methods in `ncdm_species.cpp` that used `ncdm_->*`**

Replace every `ncdm_->background_ncdm_momenta(ncdm_id_, ...)` with `ComputeMomenta(...)`.
Replace `ncdm_->q_size_ncdm_[ncdm_id_]` with `q_size()`.
Replace `ncdm_->w_ncdm_[ncdm_id_]` with `w_`.
Replace `ncdm_->q_ncdm_[ncdm_id_]` with `q_`.
Replace `ncdm_->dlnf0_dlnq_ncdm_[ncdm_id_]` with `dlnf0_dlnq_`.
Replace `ncdm_->factor_ncdm_[ncdm_id_]` with `factor_`.
Replace `ncdm_->M_ncdm_[ncdm_id_]` with `M_`.
Replace `ncdm_->GetRescalingFactor(ncdm_id_, ...)` with `GetRescalingFactor(...)`.
Replace `ncdm_->GetRescaledParameters(ncdm_id_, ...)` with `GetRescaledParameters(...)`.

- [ ] **Step 5: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: clean compile. `DNCDMSpecies` still uses `NonColdDarkMatter`; that's fine.

- [ ] **Step 6: Run a smoke test with standard NCDM**

Create `test_ncdm_standard.ini`:
```ini
output = tCl,mPk
N_ncdm = 1
m_ncdm = 0.06
deg_ncdm = 3.0
T_ncdm = 0.71611
P_k_max_1/Mpc = 1.0
```

Run: `./class test_ncdm_standard.ini`
Expected: runs to completion without error.

- [ ] **Step 7: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp
git commit -m "refactor: NCDMSpecies inherits NCDMBaseSpecies, adds CreateAll factory"
```

---

## Task 4: Migrate `DNCDMSpecies` to inherit from `NCDMBaseSpecies`

**Files:**
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`

- [ ] **Step 1: Rewrite `dncdm_species.h`**

```cpp
// species/dncdm_species.h
#pragma once
#include <memory>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;
class BackgroundColumnWriter;

class DNCDMSpecies : public NCDMBaseSpecies {
 public:
  DNCDMSpecies(FileContent* pfc, int ncdm_index, int dncdm_index,
               const NcdmSettings& settings,
               const background* pba, const BackgroundModule* bgm);

  static std::vector<std::unique_ptr<DNCDMSpecies>>
      CreateAll(FileContent* pfc, const NcdmSettings& settings,
                const background* pba, const BackgroundModule* bgm);

  // ── Background ──────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override { return pvecback[index_bg_rho_]; }
  double P(const double* pvecback) const override { return pvecback[index_bg_p_]; }
  double DpDloga(const double* pvecback) const override {
    return pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_];
  }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr,
                                   int& index_pt, const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  double Delta(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double Theta(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double DeltaP(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double RhoPlusPShear(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;

  bool IsFreestreaming() const override { return true; }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // Accessors for DNCDM_DR_Species coupling:
  int bg_number_index() const { return index_bg_number_; }
  int bg_pseudo_p_index() const { return index_bg_pseudo_p_; }
  int bg_lnf_index() const { return index_bg_lnf_decay_dr1_; }
  int bg_dlnfdlnq_index() const { return index_bg_dlnfdlnq_decay_; }
  int bg_dlnfdlnq_sep_index() const { return index_bg_dlnfdlnq_sep_; }
  int bi_lnf_index() const { return index_bi_lnf_decay_dr1_; }
  int bi_dlnfdlnq_sep_index() const { return index_bi_dlnfdlnq_separate_decay_; }

  double Gamma() const { return Gamma_; }
  const std::vector<double>& dq() const { return dq_; }

 private:
  const background* pba_;

  // ── Decay parameters (absorbed from DecayDRProperties) ───────────────────
  double Gamma_              = 0.;
  std::vector<double> dq_;
  int q_offset_              = 0;
  int quadrature_strategy_decay_ = 3;

  // ── Background indices ────────────────────────────────────────────────────
  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;

  int index_bi_lnf_decay_dr1_           = -1;
  int index_bi_dlnfdlnq_separate_decay_ = -1;

  int index_bg_lnf_decay_dr1_  = -1;
  int index_bg_dlnfdlnq_decay_ = -1;
  int index_bg_dlnfdlnq_sep_   = -1;

  int index_pt_psi0_ = -1;
};
```

- [ ] **Step 2: Implement `DNCDMSpecies::CreateAll`**

This is where the logic from `input_module.cpp` lines 1080–1190 moves. The factory reads DNCDM-specific parameters for each decay_dr species:

```cpp
std::vector<std::unique_ptr<DNCDMSpecies>>
DNCDMSpecies::CreateAll(FileContent* pfc, const NcdmSettings& settings,
                        const background* pba, const BackgroundModule* bgm) {
  std::vector<std::unique_ptr<DNCDMSpecies>> result;
  char errmsg[2048];

  // 1. Discover which ncdm indices are decay_dr (same type-list logic as NCDMSpecies::CreateAll)
  int N_ncdm = 0, flag1;
  parser_read_int(pfc, "N_ncdm", &N_ncdm, &flag1, errmsg);
  // Read type list ... collect dncdm_indices = {n : type[n] == "decay_dr"}

  if (dncdm_indices.empty()) return result;
  const int N_dncdm = dncdm_indices.size();

  // 2. Construct instances (they read their base NCDM params via NCDMBaseSpecies ctor)
  for (int i = 0; i < N_dncdm; ++i) {
    result.push_back(std::make_unique<DNCDMSpecies>(
        pfc, dncdm_indices[i], i, settings, pba, bgm));
  }

  // 3. Read DNCDM-specific Omega/deg/Gamma parameters (was lines 1080-1190 in input_module.cpp)
  //    Copy the Omega_dncdmdr / deg_ncdm_decay_dr / Omega_ini_dncdm /
  //    omega_ini_dncdm / Neff_ini_dncdm / Gamma_ncdm_decay_dr reading logic here,
  //    calling result[i]->SetOmega0(...), result[i]->SetDegAndFactor(...),
  //    result[i]->SetDeg_from_Omega_ini(...), result[i]->Gamma_ = ...
  //    All class_test validation moves here too.

  return result;
}
```

> **Copy:** The parameter-reading block (lines 1080–1190 of `source/input_module.cpp`) moves verbatim into this factory, replacing `ncdm_->` calls with `result[i]->` calls.

- [ ] **Step 3: Implement constructor and update all methods**

Constructor calls `NCDMBaseSpecies(...)` for base params, then reads `Gamma_ncdm_decay_dr[i]` and sets `Gamma_` and `dq_`.

Replace all `ncdm_->*` references in `dncdm_species.cpp` with direct member access (same mapping as Task 3 Step 4).

- [ ] **Step 4: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: clean compile.

- [ ] **Step 5: Run DNCDM smoke test**

Create `test_dncdm.ini`:
```ini
output = tCl
N_ncdm = 1
N_ncdm_decay_dr = 1
m_ncdm_decay_dr = 0.1
T_ncdm_decay_dr = 0.71611
deg_ncdm_decay_dr = 1.0
Gamma_ncdm_decay_dr = 1e3
Omega_dncdmdr = 0.12
```

Run: `./class test_dncdm.ini`
Expected: runs to completion without error.

- [ ] **Step 6: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp
git commit -m "refactor: DNCDMSpecies inherits NCDMBaseSpecies, CreateAll absorbs input_module DNCDM logic"
```

---

## Task 5: Add `DNCDM_DR_Species::CreateAll` factory

**Files:**
- Modify: `species/dncdm_dr_species.h`
- Modify: `species/dncdm_dr_species.cpp`

- [ ] **Step 1: Add declaration to `dncdm_dr_species.h`**

```cpp
// In class DNCDM_DR_Species:
static std::vector<std::unique_ptr<DNCDM_DR_Species>>
    CreateAll(FileContent* pfc, const NcdmSettings& settings,
              const background* pba, const BackgroundModule* bgm);
```

- [ ] **Step 2: Implement `DNCDM_DR_Species::CreateAll`**

```cpp
std::vector<std::unique_ptr<DNCDM_DR_Species>>
DNCDM_DR_Species::CreateAll(FileContent* pfc, const NcdmSettings& settings,
                             const background* pba, const BackgroundModule* bgm) {
  auto dncdm_vec = DNCDMSpecies::CreateAll(pfc, settings, pba, bgm);
  std::vector<std::unique_ptr<DNCDM_DR_Species>> result;
  result.reserve(dncdm_vec.size());
  for (auto& dncdm : dncdm_vec) {
    result.push_back(std::make_unique<DNCDM_DR_Species>(std::move(dncdm), pba, bgm));
  }
  return result;
}
```

> **Note:** This requires updating `DNCDM_DR_Species`'s constructor to accept `std::unique_ptr<DNCDMSpecies>` instead of `(int, shared_ptr<NonColdDarkMatter>, ...)`. Update constructor signature and implementation accordingly.

- [ ] **Step 3: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```

- [ ] **Step 4: Commit**

```bash
git add species/dncdm_dr_species.h species/dncdm_dr_species.cpp
git commit -m "feat: add DNCDM_DR_Species::CreateAll factory"
```

---

## Task 6: Migrate `DarkRadiationSpecies` to absorb `DarkRadiation`

**Files:**
- Modify: `species/dark_radiation_species.h`
- Modify: `species/dark_radiation_species.cpp`

- [ ] **Step 1: Rewrite `dark_radiation_species.h`**

Replace the `shared_ptr<DarkRadiation> dr_` member with direct data members absorbed from `DarkRadiation`:

```cpp
// species/dark_radiation_species.h
#pragma once
#include <memory>
#include <vector>

#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"
#include "parser.h"

class BackgroundModule;
class DCDMSpecies;
enum class DRType { fermion, boson };

class DarkRadiationSpecies : public BaseSpecies {
 public:
  // Returns nullptr if dcdm == nullptr or DR not requested
  static std::unique_ptr<DarkRadiationSpecies>
      Create(FileContent* pfc, const DCDMSpecies* dcdm,
             double T_cmb, double h);

  bool IsFreestreaming() const override { return true; }

  void SetBackgroundModule(const BackgroundModule* bgm) override { bgm_ = bgm; }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override { return pvecback[index_bg_rho_]; }
  double P(const double* pvecback) const override { return pvecback[index_bg_rho_] / 3.; }
  double DpDloga(const double* pvecback) const override {
    return -4. / 3. * pvecback[index_bg_rho_];
  }

  void RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr,
                                   int& index_pt, const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  double Delta(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double Theta(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double DeltaP(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  double RhoPlusPShear(const perturb_vector*, const double*, const double*, const perturb_workspace*) const override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  void WriteOutputColumns(PerturbColumnWriter&, const PerturbationsModule&,
                          enum file_format, TransferColumnSection) const override;
  void PrintVariables(PerturbColumnWriter&, double, const double*,
                      const PerturbationsModule&, const perturb_workspace*) const override;

  int bg_rho_dr_species_index() const { return index_bg_rho_dr_species_; }
  int bi_rho_dr_species_index() const { return index_bi_rho_dr_species_; }

  // Used by BackgroundModule for DR energy output
  void IntegrateDistribution(double z, double* number, double* rho, double* p) const;

 private:
  DarkRadiationSpecies(double deg, DRType dr_type, const DCDMSpecies* dcdm,
                       const background* pba, const BackgroundModule* bgm,
                       std::vector<double> q, std::vector<double> dq,
                       std::vector<double> w, double factor);

  // ── Absorbed from DarkRadiation ───────────────────────────────────────────
  std::vector<double> q_;
  std::vector<double> dq_;
  std::vector<double> w_;
  double deg_    = 1.;
  double factor_ = 0.;
  int N_q_       = 5;
  DRType dr_type_ = DRType::fermion;

  const DCDMSpecies* dcdm_ = nullptr;
  const background* pba_   = nullptr;
  const BackgroundModule* bgm_ = nullptr;

  int index_bg_rho_dr_species_ = -1;
  int index_bi_rho_dr_species_ = -1;
  int index_pt_F0_dr_species_  = -1;

  static int InitialDistribution(void* params, double q, double* f0);

  mutable ErrorMsg error_message_;
};
```

- [ ] **Step 2: Implement `DarkRadiationSpecies::Create`**

Copy `DarkRadiation::Init` from `tools/dark_radiation.cpp` into this factory, removing the multi-channel (`N_species_`) logic. The factory reads `ncdm_dr_squeezing`, quadrature settings, constructs `q_`, `dq_`, `w_`, and returns the species:

```cpp
std::unique_ptr<DarkRadiationSpecies>
DarkRadiationSpecies::Create(FileContent* pfc, const DCDMSpecies* dcdm,
                              double T_cmb, double h) {
  if (dcdm == nullptr) return nullptr;

  // Read DR parameters (from dark_radiation.cpp Init method):
  int N_q = 5, quadrature_strategy = 3;
  double qmax = 15., deg = 1.;
  // ... (parser_read_int / parser_read_double calls from dark_radiation.cpp)

  // Set up quadrature (copy from DarkRadiation::Init):
  std::vector<double> q, dq, w;
  // get_qsampling_trapezoidal(...) or equivalent

  double factor = deg * 4. * _PI_ / pow(2. * _PI_ * _hbar_, 3);

  return std::unique_ptr<DarkRadiationSpecies>(
      new DarkRadiationSpecies(deg, DRType::fermion, dcdm,
                               /*pba=*/nullptr, /*bgm=*/nullptr,
                               std::move(q), std::move(dq), std::move(w), factor));
}
```

- [ ] **Step 3: Implement `IntegrateDistribution` and remaining methods**

Copy from `dark_radiation.cpp` `IntegrateDistribution`, stripping the `index_dr` parameter (single-channel only).

Update all method bodies in `dark_radiation_species.cpp` to use direct member access (`q_`, `w_`, `factor_`, etc.) instead of `dr_->*`.

- [ ] **Step 4: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```

- [ ] **Step 5: Run DCDM+DR smoke test**

Create `test_dcdm_dr.ini`:
```ini
output = tCl
Omega_dcdmdr = 0.25
Omega_ini_dcdm = 0.30
Gamma_dcdm = 1.0
```

Run: `./class test_dcdm_dr.ini`
Expected: runs to completion without error.

- [ ] **Step 6: Commit**

```bash
git add species/dark_radiation_species.h species/dark_radiation_species.cpp
git commit -m "refactor: DarkRadiationSpecies absorbs DarkRadiation, adds Create factory"
```

---

## Task 7: Update `InputModule` — replace parameter-reading block and `ConstructSpecies`

**Files:**
- Modify: `source/input_module.h`
- Modify: `source/input_module.cpp`

- [ ] **Step 1: Update `input_module.h` — remove `target_names` entries that move into factories**

Remove `Omega_dncdmdr`, `omega_dncdmdr`, `deg_ncdm_decay_dr`, `Omega_ini_dncdm`, `Neff_ini_dncdm`, `omega_ini_dncdm` from the `target_names` enum and update `_NUM_TARGETS_`. These parameters are now handled inside `DNCDMSpecies::CreateAll`.

Remove `std::shared_ptr<NonColdDarkMatter> ncdm_;` and `std::shared_ptr<DarkRadiation> dr_;` from `InputModule` (we'll add them back temporarily as `nullptr` if anything still needs them — clean up in Task 9).

Add includes for the new factories:
```cpp
#include "../species/ncdm_species.h"
#include "../species/dncdm_dr_species.h"
#include "../species/dark_radiation_species.h"
```

- [ ] **Step 2: Replace `input_module.cpp` lines 1022–1193 with `CreateAll` calls**

Delete lines 1022–1193 (the DCDM DR assembly and all NCDM/DNCDM parameter reading). Replace with:

```cpp
  // ── Standard NCDM ────────────────────────────────────────────────────────
  NcdmSettings ncdm_settings;
  ncdm_settings.h           = pba->h;
  ncdm_settings.T_cmb       = pba->T_cmb;
  ncdm_settings.tol_ncdm    = ppr->tol_ncdm;
  ncdm_settings.tol_ncdm_bg = ppr->tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = ppr->tol_M_ncdm;

  auto ncdm_species_list = NCDMSpecies::CreateAll(pfc, ncdm_settings, pba, nullptr);
  auto dncdm_dr_list     = DNCDM_DR_Species::CreateAll(pfc, ncdm_settings, pba, nullptr);

  // DCDM dark radiation (DCDMSpecies is inserted before this point, look it up):
  const DCDMSpecies* dcdm_sp = nullptr;
  if (all_species_.count("DCDM_DR")) {
    // DCDMSpecies is a sub-component of DCDM_DR_Species;
    // reach it via a downcast or provide an accessor on DCDM_DR_Species
  }
  auto dr_sp = DarkRadiationSpecies::Create(pfc, dcdm_sp, pba->T_cmb, pba->h);
```

> **Note:** The DCDM/DR species are constructed before NCDM in `ConstructSpecies`. Since `DarkRadiationSpecies::Create` needs a `DCDMSpecies*`, you need to access it from `all_species_["DCDM_DR"]` via a cast or accessor. Add `dcdm()` accessor to `DCDM_DR_Species` if not present.

- [ ] **Step 3: Update `ConstructSpecies()`**

Replace the `if (pba->has_ncdm == _TRUE_ && ncdm_ != nullptr)` block with:

```cpp
  // Insert standard NCDM species
  for (auto& sp : ncdm_species_list) {
    std::string key = sp->name();
    all_species_[key] = std::move(sp);
  }
  // Insert DNCDM+DR composite species
  for (auto& sp : dncdm_dr_list) {
    std::string key = sp->name();
    all_species_[key] = std::move(sp);
  }
  // Insert DCDM dark radiation
  if (dr_sp) {
    all_species_["DR"] = std::move(dr_sp);
  }
```

Set `pba->N_ncdm` from the species count:
```cpp
  pba->N_ncdm = static_cast<int>(ncdm_species_list.size() + dncdm_dr_list.size());
```

- [ ] **Step 4: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: clean compile (there may be linker warnings about unused `ncdm_` / `dr_` refs in `BaseModule` — fix by temporarily keeping them as `nullptr` members).

- [ ] **Step 5: Run full scenario suite**

```bash
./class explanatory.ini
./class test_ncdm_standard.ini
./class test_dncdm.ini
./class test_dcdm_dr.ini
```
Expected: all complete without error.

- [ ] **Step 6: Commit**

```bash
git add source/input_module.h source/input_module.cpp
git commit -m "refactor: InputModule delegates NCDM/DNCDM/DR construction to CreateAll factories"
```

---

## Task 8: Replace `ncdm_->*` calls in module files

**Files:**
- Modify: `source/background_module.cpp`
- Modify: `source/perturbations_module.cpp`
- Modify: `source/nonlinear_module.cpp`

This task removes the 7 remaining usages of `ncdm_` from non-input modules.

- [ ] **Step 1: Add `GetNcdmSpecies` helper to `BackgroundModule`**

Add a private helper that returns all NCDMBaseSpecies from `all_species_`:

```cpp
// In background_module.cpp (or as a free function):
static std::vector<NCDMBaseSpecies*> GetNcdmSpecies(
    const std::map<std::string, std::unique_ptr<BaseSpecies>>& all_species) {
  std::vector<NCDMBaseSpecies*> result;
  for (auto& [key, sp] : all_species) {
    if (auto* ncdm = dynamic_cast<NCDMBaseSpecies*>(sp.get())) {
      result.push_back(ncdm);
    }
    // Also check DNCDM_DR_Species composites (access dncdm() sub-species):
    if (auto* dncdm_dr = dynamic_cast<DNCDM_DR_Species*>(sp.get())) {
      result.push_back(&dncdm_dr->dncdm());
    }
  }
  return result;
}
```

- [ ] **Step 2: Replace `background_module.cpp` line 586 — Neff**

```cpp
// Before:
Neff += ncdm_->GetNeff(0.);
ncdm_->PrintNeffInfo();

// After:
for (auto* sp : GetNcdmSpecies(all_species_)) {
  Neff += sp->GetNeff(0.);
  sp->PrintNeffInfo();
}
```

- [ ] **Step 3: Replace `background_module.cpp` line 624 — mass info**

```cpp
// Before:
ncdm_->PrintMassInfo();

// After:
for (auto* sp : GetNcdmSpecies(all_species_)) {
  sp->PrintMassInfo();
}
```

- [ ] **Step 4: Replace `background_module.cpp` line 1441 — GetIni**

```cpp
// Before:
a = ncdm_->GetIni(a, pba->a_today, ppr->tol_ncdm_initial_w);

// After:
for (auto* sp : GetNcdmSpecies(all_species_)) {
  a = sp->GetIni(a, pba->a_today, ppr->tol_ncdm_initial_w);
}
// Note: GetIni returns the smaller of a and the required start — keep looping
// to find the minimum across all species. Change logic to:
for (auto* sp : GetNcdmSpecies(all_species_)) {
  double a_sp = sp->GetIni(ppr->a_ini_over_a_today_default * pba->a_today,
                           pba->a_today, ppr->tol_ncdm_initial_w);
  a = std::min(a, a_sp);
}
```

- [ ] **Step 5: Replace `background_module.cpp` lines 1978–1979 — Omega budget**

```cpp
// Before:
ncdm_->PrintOmegaInfo();
budget_neutrino += ncdm_->GetOmega0();

// After:
for (auto* sp : GetNcdmSpecies(all_species_)) {
  sp->PrintOmegaInfo();
  budget_neutrino += sp->GetOmega0();
}
```

- [ ] **Step 6: Replace `perturbations_module.cpp` line 3452 — tensor q_size**

```cpp
// Before:
ppv->q_size_ncdm[n] = ncdm_->q_size_ncdm_[n];

// After:
ppv->q_size_ncdm[n] = ncdm_sp->q_size();
// (ncdm_sp is the NCDMBaseSpecies* already in scope from the ncdm_vec loop)
```

- [ ] **Step 7: Replace `nonlinear_module.cpp` line 1013 — GetMassInElectronvolt**

```cpp
// Before:
double m_ncdm_in_electronvolt = ncdm_->GetMassInElectronvolt(ncdm_sp->ncdm_id());

// After:
double m_ncdm_in_electronvolt = ncdm_sp->GetMassInElectronvolt();
// (ncdm_sp is already a NCDMBaseSpecies* or NCDMSpecies* in scope)
```

- [ ] **Step 8: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: clean compile.

- [ ] **Step 9: Run full scenario suite**

```bash
./class explanatory.ini
./class test_ncdm_standard.ini
./class test_dncdm.ini
./class test_dcdm_dr.ini
```

- [ ] **Step 10: Commit**

```bash
git add source/background_module.cpp source/perturbations_module.cpp source/nonlinear_module.cpp
git commit -m "refactor: replace ncdm_->* module calls with NCDMBaseSpecies species iteration"
```

---

## Task 9: Remove `ncdm_` and `dr_` from `InputModule` and `BaseModule`

**Files:**
- Modify: `source/input_module.h`
- Modify: `source/base_module.h`

- [ ] **Step 1: Remove from `base_module.h`**

```cpp
// Delete these lines:
const std::shared_ptr<NonColdDarkMatter> ncdm_;
const std::shared_ptr<DarkRadiation> dr_;

// And from the constructor initializer list:
: ncdm_(input_module->ncdm_), dr_(input_module->dr_), ...
// becomes:
: all_species_(input_module->all_species_), ...
```

Remove `#include "non_cold_dark_matter.h"` and `#include "dark_radiation.h"` from `base_module.h` and `input_module.h`.

- [ ] **Step 2: Remove from `input_module.h`**

```cpp
// Delete:
std::shared_ptr<NonColdDarkMatter> ncdm_;
std::shared_ptr<DarkRadiation> dr_;
```

- [ ] **Step 3: Compile — expect errors where `ncdm_` and `dr_` are still referenced**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4) 2>&1 | grep error
```

Fix any remaining references (should be zero if Tasks 3–8 are complete).

- [ ] **Step 4: Run full scenario suite**

```bash
./class explanatory.ini
./class test_ncdm_standard.ini
./class test_dncdm.ini
./class test_dcdm_dr.ini
```

- [ ] **Step 5: Commit**

```bash
git add source/input_module.h source/base_module.h
git commit -m "refactor: remove ncdm_ and dr_ shared_ptrs from InputModule and BaseModule"
```

---

## Task 10: Delete tool classes and clean up

**Files:**
- Delete: `tools/non_cold_dark_matter.h`, `tools/non_cold_dark_matter.cpp`
- Delete: `tools/dark_radiation.h`, `tools/dark_radiation.cpp`
- Modify: `Makefile`
- Modify: any remaining `#include` of these files

- [ ] **Step 1: Check for remaining includes**

```bash
grep -rn "non_cold_dark_matter\|dark_radiation\.h" source/ species/ include/ tools/ \
  --include="*.h" --include="*.cpp" | grep -v "^Binary"
```
Expected: zero results.

- [ ] **Step 2: Delete the tool files**

```bash
rm tools/non_cold_dark_matter.h tools/non_cold_dark_matter.cpp
rm tools/dark_radiation.h tools/dark_radiation.cpp
```

- [ ] **Step 3: Update `Makefile` — remove deleted targets**

Find and delete the lines that add `non_cold_dark_matter.opp` and `dark_radiation.opp` to `OFILES_CPP`.

- [ ] **Step 4: Compile**

```bash
make -j$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
```
Expected: clean compile.

- [ ] **Step 5: Run Python integration test suite**

```bash
cd python && python -m pytest test_class.py -v 2>&1 | tail -20
cd ..
```
Expected: all tests pass. Pay attention to NCDM-related test cases.

- [ ] **Step 6: Commit**

```bash
git add -u  # stages deletions
git add Makefile
git commit -m "chore: delete NonColdDarkMatter and DarkRadiation tool classes"
```

---

## Self-Review Checklist

**Spec coverage:**

| Spec section | Covered by |
|---|---|
| NCDMBaseSpecies data members | Tasks 1–2 |
| NCDMBaseSpecies protected constructor | Tasks 1–2 |
| NCDMSpecies inherits NCDMBaseSpecies | Task 3 |
| NCDMSpecies::CreateAll | Task 3 |
| DNCDMSpecies inherits NCDMBaseSpecies | Task 4 |
| DNCDMSpecies::CreateAll absorbs input_module lines 1080–1190 | Task 4 |
| DNCDM_DR_Species::CreateAll | Task 5 |
| DarkRadiationSpecies absorbs DarkRadiation | Task 6 |
| DarkRadiationSpecies::Create(pfc, dcdm*) | Task 6 |
| InputModule lines 1022–1193 replaced | Task 7 |
| ncdm_, dr_ removed from InputModule/BaseModule | Task 9 |
| Tool files deleted | Task 10 |
| NcdmSettings moved to ncdm_base_species.h | Task 1 |
| DecayDRProperties absorbed into DNCDMSpecies | Task 4 |
| NCDMType enum deleted | Task 10 |
| SourceType enum deleted | Task 6 |
| NoNcdmRequested / NoDrRequested deleted | Task 10 |
| background_parameters_for_distributions privatised | Task 2 |

**Type consistency notes:**
- `NCDMBaseSpecies::ComputeMomenta` signature is used identically in Tasks 2, 3, 4.
- `NCDMBaseSpecies::q_size()` used in Task 3, 4, 8 — consistent.
- `DarkRadiationSpecies::Create(pfc, dcdm*, T_cmb, h)` signature consistent across Tasks 6 and 7.
- `DNCDM_DR_Species::CreateAll` calls `DNCDMSpecies::CreateAll` internally — consistent with Task 5 calling Task 4's factory.
