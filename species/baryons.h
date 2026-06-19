#pragma once
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"
#include "species_build_context.h"

/** Baryons: rho ~ a^{-3}. Two perturbation variables: delta_b, theta_b. */
class BaryonsSpecies : public BaseSpecies {
 public:
  explicit BaryonsSpecies(const background& pba)
      : BaseSpecies("Baryons", EnergyType::Matter), pba_(pba) {}

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;  ///< index_pt_delta_b
    int idx_theta = -1;  ///< index_pt_theta_b
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  double GetOmega0() const override {
    return pba_.Omega0_b;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_b * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* /*pvecback*/) const override {
    return 0.;
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_b", 0.);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_b", pvecback[index_bg_rho_]);
  }

  // ── Perturbation index registration ────────────────────────────────────────

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  // ── PerturbDerivs ──────────────────────────────────────────────────────────

  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  // ── Source filling and initial conditions ──────────────────────────────────

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

  // ── Stress-energy observables ───────────────────────────────────────────────

  double Delta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  double Theta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  /** Baryon pressure perturbation: rho_b * (delta_p_b / rho_b) from pre-computed context. */
  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& /*layout*/,
                       const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

 private:
  const background& pba_;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
