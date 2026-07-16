#include "axion_ede_fluid.h"

#include <cmath>

#include "background.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

namespace {
// ln(1 + e^t), overflow-safe for large positive t.
double LnOnePlusExp(double t) {
  return (t > 0.) ? t + std::log1p(std::exp(-t)) : std::log1p(std::exp(t));
}
}  // namespace

AxionEDEFluid::AxionEDEFluid(const background& pba,
                             double omega0_fld,
                             double a_c,
                             double n_axion,
                             double nu,
                             double w_i,
                             double w_f,
                             double theta_i)
    : FluidSpecies(pba,
                   omega0_fld,
                   PhenoAxion,
                   /*w0_fld=*/w_i,
                   /*wa_fld=*/0.,
                   /*cs2_fld=*/1.,
                   /*Omega_EDE=*/0.),
      a_c_(a_c), n_axion_(n_axion), nu_(nu), w_i_(w_i), w_f_(w_f), theta_i_(theta_i), pba_ref_(pba),
      omega_ac_(OmegaAcFromOmegaZero(omega0_fld, a_c, nu, w_i, w_f)) {}

double AxionEDEFluid::WFinal(double n) {
  return (n - 1.) / (n + 1.);
}

double AxionEDEFluid::Integral3OnePlusWOverA(
    double a, double a_c, double nu, double w_i, double w_f) {
  // int_a^1 3(1+w)/a' da' = 3(1+w_i) ln(1/a) + nu [ ln(1+a_c^-r) - ln(1+(a/a_c)^r) ],
  // r = 3(w_f - w_i)/nu. (Derived fresh; agrees with AxiCLASS input.c:4154 for w_i = -1.)
  const double r = 3. * (w_f - w_i) / nu;
  return 3. * (1. + w_i) * std::log(1. / a) +
         nu * (LnOnePlusExp(-r * std::log(a_c)) - LnOnePlusExp(r * std::log(a / a_c)));
}

double AxionEDEFluid::OmegaZeroFromOmegaAc(
    double omega_ac, double a_c, double nu, double w_i, double w_f) {
  return omega_ac / std::exp(Integral3OnePlusWOverA(a_c, a_c, nu, w_i, w_f));
}

double AxionEDEFluid::OmegaAcFromOmegaZero(
    double omega0, double a_c, double nu, double w_i, double w_f) {
  return omega0 * std::exp(Integral3OnePlusWOverA(a_c, a_c, nu, w_i, w_f));
}

void AxionEDEFluid::ComputeWFld(double a,
                                double* w_fld,
                                double* dw_over_da_fld,
                                double* integral_fld) const {
  const double dw_tot = w_f_ - w_i_;
  const double r      = 3. * dw_tot / nu_;

  // a = 0 is a real call (background_checks probes w(a->0)); the sigmoid limit
  // is exact there: w = w_i, dw/da = 0.
  if (a <= 0.) {
    *w_fld          = w_i_;
    *dw_over_da_fld = 0.;
    *integral_fld   = 0.;
    return;
  }
  const double t = r * std::log(a_c_ / a);  // large t covers underflow near a = 0
  if (t > 700.) {
    *w_fld          = w_i_;
    *dw_over_da_fld = 0.;
    *integral_fld   = Integral3OnePlusWOverA(a, a_c_, nu_, w_i_, w_f_);
    return;
  }
  const double x  = std::exp(t);  // (a_c/a)^r
  *w_fld          = w_i_ + dw_tot / (1. + x);
  *dw_over_da_fld = dw_tot * (r / a) * x / ((1. + x) * (1. + x));
  *integral_fld   = Integral3OnePlusWOverA(a, a_c_, nu_, w_i_, w_f_);
}

double AxionEDEFluid::Cs2(double k2, double a) const {
  // GDM effective sound speed, Poulin et al. 1806.10608 (AxiCLASS perturbations.c:7802).
  const double wbar = omega_axion_ * std::pow(a, -3. * (n_axion_ - 1.) / (n_axion_ + 1.));
  const double t    = 2. * a * a * wbar * wbar;
  return ((n_axion_ - 1.) * t + k2) / ((n_axion_ + 1.) * t + k2);
}

void AxionEDEFluid::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                           double* y,
                                           const PerturbIcContext& ctx) {
  // EDE-adapted adiabatic ICs (AxiCLASS perturbations.c:5784-5795): regular as
  // w -> -1 because delta, theta carry explicit (1+w) suppression... delta does;
  // theta's (1+w) lives in rho_plus_p_theta downstream.
  const auto& layout = static_cast<const FluidSpecies::PerturbLayout&>(base);
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return;

  double w_fld, dw_over_da, integral;
  ComputeWFld(ctx.a, &w_fld, &dw_over_da, &integral);
  const double cs2   = Cs2(ctx.k * ctx.k, ctx.a);
  const double denom = 32. + 6. * cs2 + 12. * w_f_;

  y[layout.idx_delta] = 0.5 * ctx.ktau_two * (1. + w_fld) * (-4. + 3. * cs2) / denom *
                        ctx.ppr->curvature_ini * ctx.s2_squared;
  y[layout.idx_theta] = -0.5 * ctx.k * ctx.ktau_three * cs2 / denom * ctx.ppr->curvature_ini *
                        ctx.s2_squared;
}

void AxionEDEFluid::SetBackgroundModule(const BackgroundModule* bgm) {
  FluidSpecies::SetBackgroundModule(bgm);
  DeriveAxionScales();
}

void AxionEDEFluid::DeriveAxionScales() {
  // Port of AxiCLASS background.c 982-1023 (Eqs. 27, 28, 30 of 1806.10608).
  // E(a_c)^2 = sum over the other species scaled to a_c by EnergyType
  // (Radiation a^-4, Matter a^-3, DarkEnergy flat) + our own Omega_fld_ac.
  // Composites and EnergyType::Other are skipped, matching the approximate
  // nature of AxiCLASS's E(a_c) (this only calibrates omega_axion, not rho_fld).
  double omega_r = 0., omega_m = 0., omega_de = 0.;
  for (const auto& sp : bgm_->all_species_) {
    if (sp.get() == static_cast<const BaseSpecies*>(this))
      continue;
    switch (sp->energy_type()) {
      case EnergyType::Radiation:
        omega_r += sp->GetOmega0();
        break;
      case EnergyType::Matter:
        omega_m += sp->GetOmega0();
        break;
      case EnergyType::DarkEnergy:
        omega_de += sp->GetOmega0();
        break;
      case EnergyType::Other:
        break;
    }
  }

  const double n = n_axion_;
  if (n > 50.) {  // AxiCLASS guard: the Gamma-function factor degenerates
    m_fld_       = 0.;
    alpha_fld_   = 0.;
    omega_axion_ = 0.;
    return;
  }

  const double a_eq = omega_r / omega_m;
  const double p    = (a_c_ < a_eq) ? 0.5 : 2. / 3.;
  const double Eac  = std::sqrt(omega_r * std::pow(a_c_, -4.) + omega_m * std::pow(a_c_, -3.) +
                                omega_de + omega_ac_);
  const double xc   = p / Eac;
  const double f    = 7. / 8.;
  const double cos_initial = std::cos(theta_i_);
  const double sin_initial = std::sin(theta_i_);

  m_fld_           = std::pow(1. - cos_initial, (1. - n) / 2.) *
                     std::sqrt((1. - f) * (6. * p + 2.) * theta_i_ / (n * sin_initial)) / xc;
  alpha_fld_       = std::sqrt(6. * omega_ac_) / m_fld_ / std::pow(1. - cos_initial, n / 2.);
  const double Gac = std::sqrt(_PI_) * std::tgamma((n + 1.) / (2. * n)) /
                     std::tgamma(1. + 1. / (2. * n)) * std::pow(2., -(n * n + 1.) / (2. * n)) *
                     std::pow(3., 0.5 * (1. / n - 1.)) * std::pow(a_c_, 3. - 6. / (1. + n)) *
                     std::pow(std::pow(a_c_, 6. * n / (1. + n)) + 1., 0.5 * (1. / n - 1.));
  omega_axion_     = pba_ref_.H0 * m_fld_ * std::pow(1. - cos_initial, 0.5 * (n - 1.)) * Gac;
}
