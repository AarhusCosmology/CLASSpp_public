#ifndef NCDM_INTERACTING_SPECIES_H
#define NCDM_INTERACTING_SPECIES_H

#include "ncdm_species.h"

class NCDMInteractingSpecies : public NCDMSpecies {
 public:
  // Constructor inherits directly from NCDMSpecies to reuse base setup
  NCDMInteractingSpecies(FileContent* pfc,
                         int species_index,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm);

  // Factory method to read N_ncdm_interacting and create instances
  static std::vector<std::unique_ptr<NCDMInteractingSpecies>> CreateAll(
      FileContent* pfc,
      const NcdmSettings& settings,
      const background* pba,
      const BackgroundModule* bgm);

  // Override PerturbDerivs to append collision terms
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  double GetGeff() const {
    return G_eff_;
  }

 private:
  double G_eff_ = 0.0;
};

#endif  // NCDM_INTERACTING_SPECIES_H