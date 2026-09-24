#pragma once

#include <string_view>

#include "base_species.h"
#include "modified_gravity.h"
#include "species_build_context.h"

struct background;

/**
 * A scalar field non-minimally coupled to gravity, in the Jordan frame:
 *
 *   S = ∫ d^4x √-g [ f(φ) R/2 − (∂φ)²/2 − V(φ) ] + S_m ,   M_pl = 1,
 *   f(φ) = N_st + ξ_st φ² ,   V(φ) = 3 H0² Omega_Lambda_st + λ_st M_pl² φ⁴ / 4 ,
 *
 * with matter minimally coupled. That covers Early Modified Gravity (H0 Olympics,
 * arXiv:2107.10291 sec. 2.4.6; Braglia et al. 2011.12934) at N = 1; the massless
 * non-minimally coupled scalar at λ = 0; induced gravity, i.e. Brans–Dicke with
 * ω = 1/(4ξ), at N = 0; and Rock 'n' Roll early dark energy at ξ = 0. φ is in units
 * of the reduced Planck mass, and λ_st is the dimensionless coupling of the action.
 *
 * The coupling modifies gravity, so besides being a species this implements
 * ModifiedGravity:
 *   - Background: CloseFriedmann solves the Jordan-frame Friedmann equations for H and
 *     H', then reports an effective density and pressure, (.)rho_st and (.)p_st, that
 *     make GR's expressions exact. It is dark energy for the module's energy buckets;
 *     instead it scales the radiation and matter buckets by 1/f, so N_eff at BBN and
 *     z_eq see G/f. The potential's constant, Omega_Lambda_st, is shot for so that the
 *     effective density vanishes today, i.e. H(z = 0) = H0, with Lambda closing the
 *     GR budget as usual.
 *   - Perturbations: the scalar Einstein equations in synchronous gauge, and δφ with
 *     the Klein–Gordon equation in trace form. The species' stress-energy is the
 *     canonical scalar's; the coupling enters only through the Einstein equations.
 *
 * Synchronous gauge, scalar modes, adiabatic initial conditions and a flat universe
 * only; the factory refuses anything else.
 *
 * Design: docs/superpowers/specs/2026-09-24-scalar-tensor-gravity-design.md
 */
class ScalarTensorSpecies : public BaseSpecies, public ModifiedGravity {
 public:
  static constexpr std::string_view kTypeName = "scalar_tensor";

  ScalarTensorSpecies(const background& pba,
                      double xi,
                      double N,
                      double lambda,
                      double V_Lambda,
                      double phi_ini,
                      double phi_prime_ini);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  /** The coupling f(φ) = F/M_pl² and the potential V(φ), in CLASS units (1/Mpc²). */
  double F(double phi) const {
    return N_ + xi_ * phi * phi;
  }
  double V(double phi) const;

  // ── Budget, shooting, derived parameters ────────────────────────────────────
  /** The effective density today is shot to zero; see the class comment. */
  double GetOmega0() const override {
    return 0.;
  }
  /** Refused: there is no density before the background is solved, because this
   *  species changes H itself. */
  double BackgroundDensityOverH0Sq(double a, double H0) const override;
  std::optional<double> GetParam(const std::string& name) const override;
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

  // ── Background ──────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  /** The canonical scalar's p' only: the effective pressure's derivative would need
   *  H'', and the two consumers of p_tot' (N-body gauge, PPF) are refused. */
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void ProcessBackgroundTable(const double* background_table,
                              int n_rows,
                              int row_stride,
                              const double* z_table) override;

  /** omega / pi for the field's oscillation in conformal time, from the period of
   *  the quartic oscillator at the amplitude its energy implies; 0 for lambda = 0.
   *  See BaseSpecies::SourceSamplingRate. */
  double SourceSamplingRate(const double* pvecback) const override;

  double PhiPrimePrime(const double* pvecback) const {
    return pvecback[index_bg_phi_prime_prime_];
  }
  int bi_phi_index() const {
    return index_bi_phi_;
  }
  int bi_phi_prime_index() const {
    return index_bi_phi_prime_;
  }

  // ── ModifiedGravity ─────────────────────────────────────────────────────────
  void CloseFriedmann(double a,
                      double K,
                      const double* pvecback_B,
                      double* pvecback,
                      double& rho_tot,
                      double& p_tot,
                      double& rho_r,
                      double& rho_m) const override;
  double HPrime(double k, const double* y, const perturb_workspace* ppw) const override;
  double EtaPrime(double k, const double* y, const perturb_workspace* ppw) const override;
  double HPrimePrime(double k, const double* y, const perturb_workspace* ppw) const override;
  double AlphaPrime(double k, const double* y, const perturb_workspace* ppw) const override;

  // ── Perturbations ───────────────────────────────────────────────────────────
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_phi       = -1;  // δφ
    int idx_phi_prime = -1;  // δφ'
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
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      const BaseSpecies::PerturbLayout* base,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

 private:
  double dV(double phi) const;
  double ddV(double phi) const;

  /** The coupling and its perturbation at one (k, tau): delta f = f_phi delta phi and
   *  delta f' = f_phi delta phi' + f_phiphi phi' delta phi. */
  struct Coupling {
    double a = 0., aH = 0., f = 1., f_prime = 0.;
    double dphi = 0., dphi_prime = 0., df = 0., df_prime = 0.;
  };
  Coupling CouplingAt(const double* y, const perturb_workspace* ppw) const;
  /** delta phi'' from the perturbed Klein-Gordon equation in trace form. */
  double DeltaPhiPrimePrime(double k, const double* y, const perturb_workspace* ppw) const;

  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;

  double xi_            = 0.;
  double N_             = 1.;
  double lambda_class_  = 0.;  // λ M_pl², the φ⁴ coefficient ×4 in CLASS units
  double V_Lambda_      = 0.;
  double phi_ini_       = 0.;
  double phi_prime_ini_ = 0.;

  bool needs_shooting_ = false;

  // Derived today, for a massless field (Solar System tests).
  double G_eff_today_       = 1.;
  double gamma_PPN_minus_1_ = 0.;

  int index_bg_phi_             = -1;
  int index_bg_phi_prime_       = -1;
  int index_bg_phi_prime_prime_ = -1;
  int index_bg_trace_           = -1;  // T of every species plus the canonical scalar
  int index_bg_f_               = -1;
  int index_bg_p_               = -1;
  int index_bi_phi_             = -1;
  int index_bi_phi_prime_       = -1;
};
