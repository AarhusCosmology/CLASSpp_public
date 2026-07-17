#include "axion_ncdm_species.h"

#include <cmath>
#include <stdexcept>

#include "species/ncdm_family.h"
#include "species/species_input.h"

AxionNCDMSpecies::AxionNCDMSpecies(FileContent* pfc,
                                   const std::string& instance_name,
                                   const NcdmSettings& settings,
                                   const background* pba,
                                   const BackgroundModule* bgm)
    : NCDMSpecies(pfc, instance_name, settings, pba, bgm, NCDMBaseSpecies::DeferInit{}) {
  SpeciesInput input(pfc, instance_name);

  // Exactly one of T (= T_axion/T_cmb) and gstar_dec (entropy dof at
  // decoupling) fixes the axion temperature. No default: silently inheriting
  // the base's neutrino value 0.71611 would be a wrong-physics trap.
  auto T_opt     = input.get<double>("T");
  auto gstar_opt = input.get<double>("gstar_dec");
  if (T_opt && gstar_opt) {
    throw std::invalid_argument("axion species '" + instance_name +
                                "': give either T or gstar_dec, not both");
  }
  if (!T_opt && !gstar_opt) {
    throw std::invalid_argument("axion species '" + instance_name +
                                "': one of T (= T_axion/T_cmb) or gstar_dec is required");
  }
  if (T_opt && *T_opt <= 0.) {
    throw std::invalid_argument("axion species '" + instance_name +
                                "': T (= T_axion/T_cmb) must be positive");
  }
  if (gstar_opt) {
    constexpr double kGstarSToday = 43. / 11.;  // 2 + (7/8)*6*(4/11)
    if (*gstar_opt <= kGstarSToday) {
      throw std::invalid_argument("axion species '" + instance_name +
                                  "': gstar_dec must exceed today's g*S = 43/11; the axion "
                                  "cannot be hotter than the photons it decoupled from");
    }
    T_ = std::cbrt(kGstarSToday / *gstar_opt);
  }

  if (ksi_ != 0.) {
    throw std::invalid_argument("axion species '" + instance_name +
                                "': chemical potential ksi is not supported "
                                "(Bose-Einstein with mu > 0 is ill-defined)");
  }
  if (input.get_or("use_psd_file", 0) != 0) {
    throw std::invalid_argument("axion species '" + instance_name +
                                "': use_psd_file is incompatible — the Bose-Einstein PSD is "
                                "built in; use ncdm_standard for file-based PSDs");
  }

  // Object fully configured: build quadrature with the Bose-Einstein PSD
  // (dispatched via the virtual EvaluatePsdAnalytic) and resolve mass/Omega.
  BuildQuadratureAndMass(settings);
  ResolveMassOmegaClosure(settings);
}

double AxionNCDMSpecies::EvaluatePsdAnalytic(double q) const {
  // Bose-Einstein, one bosonic dof per deg unit. expm1 keeps small-q accuracy;
  // every caller evaluates at q > 0 (quadrature nodes are interior and the
  // dlnf0/dlnq stencil is bounded by dq <= (0.5-eps)*q at the first node), so
  // the integrable q -> 0 divergence is never materialized.
  return 1.0 / std::pow(2. * _PI_, 3) / std::expm1(q);
}

std::vector<Named> AxionNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  return CreateAllNcdmInstances<AxionNCDMSpecies>(ctx);
}
