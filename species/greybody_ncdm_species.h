#pragma once
#include <optional>
#include <string>
#include <vector>

#include "../species/greybody_moments.h"
#include "../species/ncdm_species.h"
#include "../species/species_build_context.h"

/** NCDM with a grey-body phase-space distribution:
 *  f0(q) = 2/(2π)³ (q/q0)^(α-1) / (exp(α x q) + 1).
 *  Parameterized directly by (alpha, q0, x) or by inverting velocity moments
 *  (r, M2, M3). All perturbation/background behavior is inherited from NCDMSpecies. */
class GreyBodyNCDMSpecies : public NCDMSpecies {
 public:
  GreyBodyNCDMSpecies(FileContent* pfc,
                      const std::string& instance_name,
                      const NcdmSettings& settings,
                      const background* pba,
                      const BackgroundModule* bgm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // Per-species derived parameters surfaced via BackgroundModule::GetSpeciesParam.
  std::optional<double> GetParam(const std::string& name) const override;

 protected:
  double EvaluatePsdAnalytic(double q) const override;
  quadrature_method DefaultQuadratureStrategy() const override;
  void FillQuadratureParams(GBQuadParams& p) const override;

 private:
  // Computes the species' dimensionless moments from its own quadrature.
  void GreyBodyMoments(double& M2, double& M3, double& M4) const;

  greybody::GreyBodyParams gb_ = greybody::GreyBodyParams::FromDirect(1., 1., 1.);
  bool gb_ready_               = false;
};
