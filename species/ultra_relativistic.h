#pragma once

#include "base_species.h"
#include "species_build_context.h"

struct background;

/**
 * Ultra-relativistic neutrinos / relics (UR) species.
 *
 * Background:  rho_ur = Omega0_ur * H0^2 / a_rel^4,
 *              p_ur   = rho_ur / 3,
 *              dp/dloga = -4/3 * rho_ur.
 * Perturbations: full Boltzmann hierarchy in l (delta, theta, shear, l3..l_max).
 *   Supports UFA (fluid approximation) and RSA (radiation streaming).
 *   Non-standard ceff2_ur and cvis2_ur terms from ppt->three_ceff2_ur etc.
 */
class UltraRelativisticSpecies : public BaseSpecies {
 public:
  UltraRelativisticSpecies(const background& pba, double omega0_ur);

  double GetOmega0() const override {
    return Omega0_ur_;
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  bool IsFreestreaming() const override {
    return true;
  }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────

  /** Per-species layout holding perturbation slot indices for UR. */
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;  ///< index_pt_delta_ur
    int idx_theta = -1;  ///< index_pt_theta_ur
    int idx_shear = -1;  ///< index_pt_shear_ur
    int idx_l3    = -1;  ///< index_pt_l3_ur (only when ufa_off and l_max>=3)
    int l_max     = -1;  ///< 2 when ufa_on, else ppr->l_max_ur
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  /** Registers delta_ur, theta_ur, shear_ur, and (when UFA is off) l3..l_max_ur. */
  // Layout-based override (primary path)
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  // Layout-based PerturbDerivs (primary path)
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void PerturbDerivs(double /*tau*/,
                     const double* /*y*/,
                     double* /*dy*/,
                     const perturb_parameters_and_workspace& /*ppaw*/) override {}

  // Layout-based Delta/Theta/DeltaP/RhoPlusPShear
  double Delta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Delta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double Theta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* /*pv*/,
                const double* /*y*/,
                const double* /*pvecback*/,
                const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  // Layout-based FillSources and ApplyInitialConditions
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) override;
  void FillSources(const double* /*y*/,
                   const double* /*dy*/,
                   PerturbSourceContext& /*ctx*/) override {}

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  void ApplyInitialConditions(double* /*y*/, const PerturbIcContext& /*ctx*/) override {}

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  /** Scalar mode: UR l>=3 multipoles not needed in sources when both rsa and ufa are off. */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& layout,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;

  /** Copy UR perturbations across an approximation switch (UFA on/off transition). */
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  // ── Tensor mode ────────────────────────────────────────────────────────────
  //
  // The tensor "relativistic neutrino" hierarchy is NOT a UR-species concern: it
  // is a single hierarchy sourced by the gravitational waves and shared by ur and
  // massless-approximated ncdm. perturb_vector owns it (pv->tensor_ur_layout) and
  // the perturbations module registers/evolves it directly, so URSpecies has no
  // tensor index-registration or tensor-derivs overrides.

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      enum file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

 private:
  double Omega0_ur_;
  double H0_;
};
