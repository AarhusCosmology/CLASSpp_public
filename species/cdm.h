#pragma once

#include "base_species.h"
#include "species_build_context.h"

struct background;

/**
 * Cold Dark Matter (CDM) species.
 *
 * Background:  rho_cdm = Omega0_cdm * H0^2 / a_rel^3,  p = 0
 * Perturbations: pressureless fluid, no anisotropic stress.
 *   Newtonian gauge: evolves delta_cdm and theta_cdm.
 *   Synchronous gauge: evolves delta_cdm only (theta = 0 by gauge choice).
 */
class CDMSpecies : public BaseSpecies {
 public:
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;
    int idx_theta = -1;  // newtonian gauge only
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  explicit CDMSpecies(const background& pba, double omega0_cdm);

  double GetOmega0() const override {
    return Omega0_cdm_;
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

 private:
  double Omega0_cdm_;
  double H0_;
  int index_bg_rho_cdm_ = -1;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;  // registered only when gauge != synchronous
};
