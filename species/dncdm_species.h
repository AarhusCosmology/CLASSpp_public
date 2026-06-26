#pragma once
#include <cmath>
#include <memory>
#include <optional>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "../species/species_build_context.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Decaying Non-Cold Dark Matter (DNCDM).
 * Inherits NCDMBaseSpecies; owns per-species quadrature, distribution function,
 * Gamma decay rate, and dq volume elements absorbed from DecayDRProperties.
 */
class DNCDMSpecies : public NCDMBaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "ncdm_decay_dr";

  // Reads all DNCDM-specific parameters from the dot-syntax instance
  // identified by instance_name (e.g. "dncdm1").
  DNCDMSpecies(FileContent* pfc,
               const std::string& instance_name,
               const NcdmSettings& settings,
               const background* pba,
               const BackgroundModule* bgm);

  // Accessors for deferred closure (used by CreateAll after construction)
  const std::optional<double>& Omega_ini_pending() const {
    return Omega_ini_pending_;
  }
  const std::optional<double>& Neff_ini_pending() const {
    return Neff_ini_pending_;
  }
  const std::optional<double>& Omega_dncdmdr_pending() const {
    return Omega_dncdmdr_pending_;
  }

  /** True iff this flavor is normalized by initial abundance (Omega_ini/omega_ini/Neff_ini) —
   *  the mode that needs the Omega_dncdmdr fixed-point shoot for closure (vs combined mode,
   *  which shoots deg). */
  bool InitialAbundanceMode() const {
    return Omega_ini_pending_.has_value() || Neff_ini_pending_.has_value();
  }

  /** Backfill this species' today density fraction Omega0_ from its integrated density
   *  at a=1. In combined/initial modes the matter child is normalized via `deg`, never
   *  SetOmega0, so GetOmega0() stays 0 — which drops it from fnu (GetOmega0NcdmTot) and
   *  the budget-print neutrino line. BackgroundModule calls this post-integration.
   *  (SetOmega0 is protected; this public wrapper lets the module trigger the backfill.) */
  void BackfillOmega0FromToday(const double* pvecback_today, double H0, double h) {
    SetOmega0(Rho(pvecback_today) / (H0 * H0), h);
  }

  // Compute a (deg_guess, dxdy) pair for Newton shooting that varies deg to hit
  // a target today-density Omega_target = (rho_dncdm + rho_dr) / H0^2 at z=0.
  // Ported from input_module.cpp:3759-3800 (single-flavor, no loop).
  std::pair<double, double> DegGuessFromOmegaToday(const SpeciesBuildContext& ctx,
                                                   double Omega_target) const;

  struct Named {
    std::string key;
    std::unique_ptr<DNCDMSpecies> species;
  };

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Background ──────────────────────────────────────────────────────────
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
  double PPrime(double a,
                double H,
                const double* /*pvecback_B*/,
                const double* pvecback) const override {
    return a * H * (pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_]);
  }

  // ── Perturbations ────────────────────────────────────────────────────────

  // Layout-based scalar register (writes both layout and legacy pv arrays).
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  // Scalar PerturbDerivs (called by DNCDM_DR_Species composite with my.dncdm).
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  /** Fused override: one RescaledPerturbations call for delta/theta/shear;
   *  DeltaP uses its own independent loop (different weights), matching the
   *  individual methods' per-term expressions and operand order exactly. */
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  bool IsFreestreaming() const override {
    return true;
  }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Accessors for DNCDM_DR_Species coupling ───────────────────────────────
  int bg_number_index() const {
    return index_bg_number_;
  }
  int bg_pseudo_p_index() const {
    return index_bg_pseudo_p_;
  }
  int bg_lnf_index() const {
    return index_bg_lnf_decay_dr1_;
  }
  int bg_dlnfdlnq_index() const {
    return index_bg_dlnfdlnq_decay_;
  }
  int bg_dlnfdlnq_sep_index() const {
    return index_bg_dlnfdlnq_sep_;
  }
  int bi_lnf_index() const {
    return index_bi_lnf_decay_dr1_;
  }
  int bi_dlnfdlnq_sep_index() const {
    return index_bi_dlnfdlnq_separate_decay_;
  }

  double Gamma() const {
    return Gamma_;
  }
  const std::vector<double>& dq() const {
    return dq_;
  }
  double GetMass() const {
    return M_;
  }
  const std::vector<double>& GetQ() const {
    return q_;
  }

  // Override GetRescaledParameters to use dq_ (not standard w_ weights)
  std::tuple<double, double> GetRescaledParameters(double a,
                                                   const double* lnf_array) const override;

  /**
   * Returns rescaled (delta, theta, shear) for this decaying NCDM flavor.
   * Rescaling subtracts a common lnN from every lnf to prevent exp(lnf)
   * underflow near the precision floor; lnN cancels in the delta/theta/shear
   * ratios, so this is mathematically equivalent to the unrescaled form but
   * numerically stable.
   *
   * @param layout this species' per-pv NCDM layout (provides index_per_q).
   * @param a Scale factor at the current integration time.
   * @param k Fourier wavenumber of the perturbation mode.
   * @param ppw Perturbation workspace holding the current phase-space state.
   */
  std::tuple<double, double, double> RescaledPerturbations(
      const NCDMBaseSpecies::PerturbLayout& layout,
      double a,
      double k,
      const perturb_workspace* ppw) const;

 protected:
  double GetDlnf0Dlnq(int iq, const double* pvecback) const override {
    return pvecback[index_bg_dlnfdlnq_decay_ + iq];
  }
  double GetW0ForGwSource(int iq, const double* pvecback) const override {
    return dq_[iq] * std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]);
  }

 private:
  const background* pba_;

  // Deferred closure stash (instance-name constructor only).
  // Set when the user specifies Omega_ini/omega_ini or Neff_ini; CreateAll
  // applies SetDeg_from_Omega_ini once a_ini is available.
  std::optional<double> Omega_ini_pending_;
  std::optional<double> Neff_ini_pending_;
  // Today-density target: shoot deg so that (rho_dncdm+rho_dr)/H0^2 == this value at z=0.
  std::optional<double> Omega_dncdmdr_pending_;

  // Absorbed from DecayDRProperties
  double Gamma_ = 0.;
  std::vector<double> dq_;

  // Background indices
  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;

  int index_bi_lnf_decay_dr1_           = -1;
  int index_bi_dlnfdlnq_separate_decay_ = -1;

  int index_bg_lnf_decay_dr1_  = -1;
  int index_bg_dlnfdlnq_decay_ = -1;
  int index_bg_dlnfdlnq_sep_   = -1;

  // Perturbation indices
  int index_pt_psi0_ = -1;
};
