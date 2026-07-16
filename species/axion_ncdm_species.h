#pragma once
#include <string>
#include <string_view>
#include <vector>

#include "../species/ncdm_species.h"
#include "../species/species_build_context.h"

/** NCDM with a Bose-Einstein phase-space distribution — a thermal axion (or
 *  any thermal boson) that decoupled from the SM plasma while relativistic:
 *    f0(q) = 1/(2π)³ · 1/(e^q − 1)   per internal dof.
 *  NOTE: deg = 1 means ONE bosonic dof (a real scalar) — unlike ncdm_standard,
 *  where deg = 1 means a full neutrino family (2 fermionic dof).
 *  Temperature is given directly (T = T_a/T_cmb) or via the entropy dof at
 *  decoupling: T = (43/11 / gstar_dec)^(1/3). See arXiv:1307.0615.
 *  All background/perturbation behavior is inherited from NCDMSpecies. */
class AxionNCDMSpecies : public NCDMSpecies {
 public:
  static constexpr std::string_view kTypeName = "ncdm_axion";

  AxionNCDMSpecies(FileContent* pfc,
                   const std::string& instance_name,
                   const NcdmSettings& settings,
                   const background* pba,
                   const BackgroundModule* bgm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

 protected:
  double EvaluatePsdAnalytic(double q) const override;
};
