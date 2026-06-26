#ifndef NCDM_INTERACTING_SPECIES_H
#define NCDM_INTERACTING_SPECIES_H

#include <string_view>

#include "ncdm_species.h"
#include "species/species_build_context.h"

class NCDMInteractingSpecies : public NCDMSpecies {
 public:
  static constexpr std::string_view kTypeName = "ncdm_self_interacting";

  // New input path: parameters read from PFC under <instance_name>.<field>
  NCDMInteractingSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm);

  // Factory method to read N_ncdm_interacting and create instances
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  double FitIntegralOfl(double z) const;

  // Layout-based PerturbDerivs: runs NCDMSpecies hierarchy + collision terms.
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;

  double GetGeff() const {
    return G_eff_;
  }

 private:
  double G_eff_              = 0.0;
  bool use_alpha_correction_ = _FALSE_;
};

#endif  // NCDM_INTERACTING_SPECIES_H