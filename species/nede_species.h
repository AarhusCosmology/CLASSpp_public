#pragma once

#include <string_view>

#include "base_species.h"
#include "species_build_context.h"

struct background;

/**
 * Cold New Early Dark Energy (Niedermann & Sloth, arXiv:1910.10739 and 2006.06686; the
 * NEDE model of the H0 Olympics, arXiv:2107.10291 sec. 2.4.5), following TriggerCLASS.
 *
 * A vacuum energy decays in a first-order phase transition once the expansion rate falls
 * to H_* = H_over_m_NEDE m_NEDE, where m_NEDE is the mass of a subdominant trigger field:
 *
 *   a <  a_*:  rho = f_NEDE H_*^2,                        w = -1
 *   a >= a_*:  rho = f_NEDE H_*^2 (a_* / a)^{3(1+w_NEDE)},   w = c_s^2 = w_NEDE
 *
 * f_NEDE is NEDE's share of the total density at the transition, where rho_tot = H_*^2.
 *
 * NEDE has no perturbations before the transition and a fluid's delta and theta after it,
 * so it owns an approximation with these two regimes. Before the transition the
 * trigger's delta_phi is evolved instead. The transition happens on a surface of constant
 * trigger field, and the junction conditions of 2006.06686 seed the fluid from it:
 *
 *   delta = -3 (1+w) aH delta_phi/phi',   theta = k^2 delta_phi/phi',   shear = 0.
 *
 * The trigger's potential m^2 phi^2/2 makes its equations linear, so delta_phi/phi' does
 * not depend on the field's amplitude: the trigger is a clock rather than a component. It
 * is left out of the energy budget and frozen after the transition. Synchronous gauge only.
 */
class NEDESpecies : public BaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "nede";

  NEDESpecies(const background& pba,
              double f_NEDE,
              double m_NEDE,
              double w_NEDE,
              double H_star,
              double a_star);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  double GetOmega0() const override;
  std::optional<double> GetParam(const std::string& name) const override;

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  enum Regime { kBeforeTransition = 0, kAfterTransition = 1 };

  int ApproximationCount() const override {
    return 1;
  }
  int ApproximationRegimeAt(int which, double k, const double* pvecback) const override;

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta_phi       = -1;  ///< trigger, before the transition
    int idx_delta_phi_prime = -1;
    int idx_delta           = -1;  ///< NEDE fluid, after it
    int idx_theta           = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

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
                     const perturb_parameters_and_workspace& ppaw) const override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

 private:
  double RhoAt(double a) const;
  double WAt(double a) const {
    return (a < a_star_) ? -1. : w_NEDE_;
  }

  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;

  double m_NEDE_;    ///< trigger mass [1/Mpc]
  double w_NEDE_;    ///< equation of state and sound speed after the transition
  double rho_star_;  ///< rho before the transition, f_NEDE H_*^2
  double a_star_;    ///< scale factor of the transition

  int index_bg_phi_prime_ = -1;  ///< trigger phi'
  int index_bi_phi_       = -1;
  int index_bi_phi_prime_ = -1;
};
