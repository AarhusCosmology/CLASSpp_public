#pragma once
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"
#include "species_build_context.h"

/** Photons: rho ~ a^{-4}. Boltzmann hierarchy with l_max_g multipoles. */
class PhotonsSpecies : public BaseSpecies {
 public:
  explicit PhotonsSpecies(const background& pba)
      : BaseSpecies("Photons", EnergyType::Radiation), pba_(pba) {}

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    // Scalar / vector / tensor — same struct holds whichever mode this pv is for.
    int idx_delta = -1, idx_theta = -1, idx_shear = -1, idx_l3 = -1;
    int idx_pol0 = -1, idx_pol1 = -1, idx_pol2 = -1, idx_pol3 = -1;
    int l_max     = -1;
    int l_max_pol = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  double GetOmega0() const override {
    return pba_.Omega0_g;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_g * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_rho_] / 3.;
  }
  double PPrime(double a,
                double H,
                const double* /*pvecback_B*/,
                const double* pvecback) const override {
    return a * H * (-4. / 3. * pvecback[index_bg_rho_]);
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", 0.);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", pvecback[index_bg_rho_]);
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

  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& layout,
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

  void PerturbVectorDerivs(const BaseSpecies::PerturbLayout& layout,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) override;

  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& layout,
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

  // ── MarkUsedInSources ──────────────────────────────────────────────────────

  /** Photon temperature l>=3 and polarization l=1,3+ are not needed in source
   *  evaluation when both rsa and tca are off. Called only in scalar mode
   *  (the dispatch loop in the module is guarded by if (_scalars_)). */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& layout,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;

  /** In tensor mode, photon temperature l=0,2,4 and pol l=0,2,4 are needed;
   *  mark higher multipoles as unused when both rsa and tca are off. */
  void MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& layout,
                               const perturb_workspace* ppw,
                               int* used_in_sources) const override;

  // ── Switch-copy hook ────────────────────────────────────────────────────────

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

  double DeltaRho(const BaseSpecies::PerturbLayout& layout,
                  const perturb_vector* pv,
                  const double* y,
                  const double* pvecback,
                  const perturb_workspace* ppw) const override;

  double RhoPlusPTheta(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  /** δp_g = δρ_g / 3 = ρ_g * δ_g / 3. Returns 0 when RSA is active. */
  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override {
    thm_ = thm;
  }
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
  }

 private:
  const background& pba_;
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_             = nullptr;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
