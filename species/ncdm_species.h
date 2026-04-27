#pragma once
#include <memory>
#include <string>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "../species/species_build_context.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

class NCDMSpecies : public NCDMBaseSpecies {
 public:
  NCDMSpecies(FileContent* pfc,
              int species_index,
              const NcdmSettings& settings,
              const background* pba,
              const BackgroundModule* bgm,
              const std::string& suffix = "_standard");

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  bool IsFreestreaming() const override {
    return true;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
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
  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  double Delta(const perturb_vector*,
               const double*,
               const double*,
               const perturb_workspace*) const override;
  double Theta(const perturb_vector*,
               const double*,
               const double*,
               const perturb_workspace*) const override;
  double DeltaP(const perturb_vector*,
                const double*,
                const double*,
                const perturb_workspace*) const override;
  double RhoPlusPShear(const perturb_vector*,
                       const double*,
                       const double*,
                       const perturb_workspace*) const override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void WriteOutputColumns(PerturbColumnWriter&,
                          const PerturbationsModule&,
                          enum file_format,
                          TransferColumnSection) const override;
  void PrintVariables(PerturbColumnWriter&,
                      double,
                      const double*,
                      const PerturbationsModule&,
                      const perturb_workspace*) const override;

  int ncdm_id() const {
    return ncdm_id_;
  }
  int bg_number_index() const {
    return index_bg_number_;
  }
  int bg_pseudo_p_index() const {
    return index_bg_pseudo_p_;
  }

 protected:
  int ncdm_id_;  // species index (0-based), used for pv->index_ncdm_ etc.
  const background* pba_;

  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;
  int index_pt_psi0_     = -1;
};
