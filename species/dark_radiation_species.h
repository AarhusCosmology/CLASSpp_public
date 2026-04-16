#pragma once
#include <memory>

#include "../species/base_species.h"
#include "../tools/dark_radiation.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;
class DCDMSpecies;

/**
 * Dark Radiation (from DCDM decay and/or decaying NCDM).
 * rho_dr stored per decay-channel in the ODE integration vector; total also stored.
 */
class DarkRadiationSpecies : public BaseSpecies {
 public:
  DarkRadiationSpecies(std::shared_ptr<DarkRadiation> dr,
                       const background* pba,
                       const BackgroundModule* bgm,
                       const DCDMSpecies* dcdm = nullptr)
      : BaseSpecies("DR", EnergyType::Radiation), dr_(std::move(dr)), pba_(pba), bgm_(bgm),
        dcdm_(dcdm) {}

  bool IsFreestreaming() const override {
    return true;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_rho_] / 3.;
  }
  double DpDloga(const double* pvecback) const override {
    return -4. / 3. * pvecback[index_bg_rho_];
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
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

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

  int bg_rho_dr_species_index() const {
    return index_bg_rho_dr_species_;
  }
  int bi_rho_dr_species_index() const {
    return index_bi_rho_dr_species_;
  }

 private:
  std::shared_ptr<DarkRadiation> dr_;
  const background* pba_;
  const BackgroundModule* bgm_;
  const DCDMSpecies* dcdm_ = nullptr;  // optional: set when created inside DCDM_DR_Species

  // Background indices (per-channel, then total)
  int index_bg_rho_dr_species_ = -1;  // first of N_decay_dr contiguous slots

  // Integration indices
  int index_bi_rho_dr_species_ = -1;  // first of N_decay_dr ODE slots

  // Perturbation indices
  int index_pt_F0_dr_species_ = -1;  // per-species multipoles (N_decay_dr*(l_max_dr+1))
};
