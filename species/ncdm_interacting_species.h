#ifndef NCDM_INTERACTING_SPECIES_H
#define NCDM_INTERACTING_SPECIES_H

#include <string_view>

#include "ncdm_species.h"
#include "species/species_build_context.h"

class NCDMInteractingSpecies : public NCDMSpecies {
 public:
  /** Deliberately NOT inheriting the parent's opt-in: this sector adds stiffness
   *  and has not been benchmarked against ndf15 with an explicit evolver.
   *  See BaseSpecies::SupportsExplicitPerturbationEvolver(). */
  bool SupportsExplicitPerturbationEvolver() const override {
    return false;
  }
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
  bool use_alpha_correction_ = false;
};

#endif  // NCDM_INTERACTING_SPECIES_H