#include "greybody_ncdm_species.h"

#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "species/species_input.h"

namespace {
constexpr double kGreyBodyLeadingFactor = 2.0;  // the leading 2 in f0
}

GreyBodyNCDMSpecies::GreyBodyNCDMSpecies(FileContent* pfc,
                                         const std::string& instance_name,
                                         const NcdmSettings& settings,
                                         const background* pba,
                                         const BackgroundModule* bgm)
    : NCDMSpecies(pfc, instance_name, settings, pba, bgm, NCDMBaseSpecies::DeferInit{}) {
  SpeciesInput input(pfc, instance_name);

  std::string mode;
  if (!input.read_string("parameterization", mode)) {
    throw std::invalid_argument("grey-body species '" + instance_name +
                                "': required field 'parameterization' is missing "
                                "(set it to 'direct' or 'moments')");
  }
  if (mode != "direct" && mode != "moments") {
    throw std::invalid_argument("grey-body species '" + instance_name +
                                "': parameterization must be 'direct' or 'moments'");
  }

  if (mode == "direct") {
    double alpha = 0., q0 = 0., x = 0.;
    bool ok  = input.read_double("alpha", alpha);
    ok      &= input.read_double("q0", q0);
    ok      &= input.read_double("x", x);
    if (!ok) {
      throw std::invalid_argument("grey-body species '" + instance_name +
                                  "': direct mode requires alpha, q0, x");
    }
    gb_ = greybody::GreyBodyParams::FromDirect(alpha, x, q0);
  }
  else {
    double r = 0., M2 = 0., M3 = 0.;
    bool ok  = input.read_double("r", r);
    ok      &= input.read_double("M2", M2);
    ok      &= input.read_double("M3", M3);
    if (!ok) {
      throw std::invalid_argument("grey-body species '" + instance_name +
                                  "': moments mode requires r, M2, M3");
    }
    gb_ = greybody::GreyBodyParams::FromMoments(r, M2, M3);

    // Optional mass-via-M2 convenience: m = m_M2 / M2.
    double m_M2 = 0.;
    if (input.read_double("m_M2", m_M2)) {
      m_in_eV_ = m_M2 / M2;
    }
    if (m_in_eV_ != 0.0 && m_in_eV_ < 0.01) {
      throw std::invalid_argument(
          "grey-body species '" + instance_name +
          "': inferred mass < 0.01 eV; moment inversion is unreliable in the low-mass limit");
    }
  }
  gb_ready_ = true;

  // Object is now fully configured: build quadrature with the grey-body PSD
  // (dispatched via the virtual EvaluatePsdAnalytic) and resolve mass/Omega.
  BuildQuadratureAndMass(settings);
  ResolveMassOmegaClosure(settings);
}

double GreyBodyNCDMSpecies::EvaluatePsdAnalytic(double q) const {
  if (!gb_ready_) {
    // Invariant: BuildQuadratureAndMass (the only caller path into here) is
    // deferred until gb_ is configured, so this is unreachable. Fail loudly
    // rather than silently returning the wrong (Fermi-Dirac) PSD.
    throw std::logic_error("GreyBodyNCDMSpecies::EvaluatePsdAnalytic called before gb_ready_");
  }
  // Log-space evaluation for numerical stability at large parameters. The
  // denominator log(exp(t)+1) is computed as the stable softplus t + log1p(e^-t)
  // (t = alpha*x*q >= 0 here, since all quadrature nodes q are positive) so the
  // exponential never overflows for large t.
  double t    = gb_.x_times_alpha() * q;
  double logf = (gb_.alpha() - 1.) * std::log(q) - gb_.alpham1_logq0() -
                (t + std::log1p(std::exp(-t)));
  return kGreyBodyLeadingFactor / std::pow(2. * _PI_, 3) * std::exp(logf);
}

quadrature_method GreyBodyNCDMSpecies::DefaultQuadratureStrategy() const {
  return qm_GB_Laguerre;
}

void GreyBodyNCDMSpecies::FillQuadratureParams(GBQuadParams& p) const {
  p.active        = true;
  p.alpha         = gb_.alpha();
  p.x_times_alpha = gb_.x_times_alpha();
}

void GreyBodyNCDMSpecies::GreyBodyMoments(double& M2, double& M3, double& M4) const {
  // Report the moments of the *actual* sampled distribution (q_bg_/w_bg_, whose
  // weights already fold in f0), not the analytic gb_.moments(): this reflects
  // what the background integration used, including quadrature resolution.
  const std::vector<double>& qv = q_bg_;
  const std::vector<double>& wv = w_bg_;
  M2 = M3 = M4 = 0.;
  for (size_t i = 0; i < qv.size(); ++i) {
    double q2  = qv[i] * qv[i];
    M2        += q2 * wv[i];
    M3        += q2 * qv[i] * wv[i];
    M4        += q2 * q2 * wv[i];
  }
  double pref  = 0.5 * std::pow(2. * _PI_, 3);
  M2          *= pref;
  M3          *= pref;
  M4          *= pref;
}

std::optional<double> GreyBodyNCDMSpecies::GetParam(const std::string& name) const {
  if (name == "alpha")
    return gb_.alpha();
  if (name == "q0")
    return gb_.q0();
  if (name == "x")
    return gb_.x();
  if (name == "M_2" || name == "M_3" || name == "M_4") {
    double M2, M3, M4;
    GreyBodyMoments(M2, M3, M4);
    if (name == "M_2")
      return M2;
    if (name == "M_3")
      return M3;
    return M4;
  }
  return BaseSpecies::GetParam(name);
}

std::vector<Named> GreyBodyNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  const auto instances = ctx.pfc->instances_with("type", "ncdm_greybody");
  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);
    auto sp =
        std::make_unique<GreyBodyNCDMSpecies>(ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
