#pragma once
#include <memory>
#include <tuple>
#include <vector>

#include "../species/ncdm_base_species.h"
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
  // ncdm_index:  overall NCDM index (position in N_ncdm list)
  // dncdm_index: 0-based index among decay_dr species only
  DNCDMSpecies(FileContent* pfc,
               int ncdm_index,
               int dncdm_index,
               const NcdmSettings& settings,
               const background* pba,
               const BackgroundModule* bgm);

  static std::vector<std::unique_ptr<DNCDMSpecies>> CreateAll(FileContent* pfc,
                                                              const NcdmSettings& settings,
                                                              const background* pba,
                                                              const BackgroundModule* bgm);

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

 private:
  int ncdm_id_;  // overall NCDM index
  const background* pba_;

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
