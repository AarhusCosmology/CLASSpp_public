#pragma once
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "../species/species_build_context.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;
class BackgroundColumnWriter;

/**
 * Decaying Non-Cold Dark Matter (DNCDM).
 * Inherits NCDMBaseSpecies; owns per-species quadrature, distribution function,
 * Gamma decay rate, and dq volume elements absorbed from DecayDRProperties.
 */
class DNCDMSpecies : public NCDMBaseSpecies {
 public:
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

  struct Named {
    std::string key;
    std::unique_ptr<DNCDMSpecies> species;
  };

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Background ──────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double DpDloga(const double* pvecback) const override {
    return pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_];
  }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  double Delta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  bool IsFreestreaming() const override {
    return true;
  }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Accessors for DNCDM_DR_Species coupling ───────────────────────────────
  int ncdm_id() const {
    return ncdm_id_;
  }
  void SetNcdmId(int id) override {
    ncdm_id_ = id;
  }
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
   */
  std::tuple<double, double, double> RescaledPerturbations(double a,
                                                           double k,
                                                           const perturb_workspace* ppw) const;

 private:
  int ncdm_id_ = -1;  // perturbation-array slot index; assigned by CreateAll via SetNcdmId
  const background* pba_;

  // Deferred closure stash (instance-name constructor only).
  // Set when the user specifies Omega_ini/omega_ini or Neff_ini; CreateAll
  // applies SetDeg_from_Omega_ini once a_ini is available.
  std::optional<double> Omega_ini_pending_;
  std::optional<double> Neff_ini_pending_;

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
