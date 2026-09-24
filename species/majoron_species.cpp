#include "majoron_species.h"

#include <cmath>

#include "background.h"
#include "errors.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"

namespace {

/* K_2(x) = K_0(x) + (2/x) K_1(x), from Abramowitz & Stegun 9.8.1-9.8.8 (libc++ has no
   std::cyl_bessel_k): relative error below 1.3e-7 against scipy for 1e-8 < x < 630. */
double BesselK2(double x) {
  if (x <= 2.) {
    const double t = x * x / (3.75 * 3.75);
    const double I0 =
        1. +
        t * (3.5156229 +
             t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
    const double I1 =
        x *
        (0.5 + t * (0.87890594 +
                    t * (0.51498869 +
                         t * (0.15084934 + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
    const double y = x * x / 4.;
    const double K0 =
        -log(x / 2.) * I0 +
        (-0.57721566 +
         y * (0.42278420 +
              y * (0.23069756 +
                   y * (0.03488590 + y * (0.00262698 + y * (0.00010750 + y * 0.0000074))))));
    const double K1 =
        log(x / 2.) * I1 +
        (1. +
         y * (0.15443144 +
              y * (-0.67278579 +
                   y * (-0.18156897 + y * (-0.01919402 + y * (-0.00110404 + y * -0.00004686)))))) /
            x;
    return K0 + 2. / x * K1;
  }
  const double u = 2. / x;
  const double K0 =
      1.25331414 +
      u * (-0.07832358 +
           u * (0.02189568 +
                u * (-0.01062446 + u * (0.00587872 + u * (-0.00251540 + u * 0.00053208)))));
  const double K1 =
      1.25331414 +
      u * (0.23498619 +
           u * (-0.03655620 +
                u * (0.01504268 + u * (-0.00780353 + u * (0.00325614 + u * -0.00068245)))));
  return exp(-x) / sqrt(x) * (K0 + 2. / x * K1);
}

}  // namespace

MajoronSpecies::MajoronSpecies(const background& pba, double m_phi, double Gamma_eff)
    : UltraRelativisticSpecies(pba, 0., "Majoron", "nuphi"),
      Gamma_phi_(Gamma_eff * 4e-14 * 4e-14 * (m_phi / 0.1) * m_phi / (16. * _PI_)),
      x0_(pba.T_cmb * _k_B_ / (_eV_ * m_phi)), H0_(pba.H0),
      Omega_nu_per_F_(pba.Omega0_g * 15. / (_PI_ * _PI_)) {
  const double s  = sqrt(Gamma_eff);
  const double e1 = exp(2.79816 * s + 1.);
  const double e2 = exp(1.41748 * s + 1.);
  alpha1_         = 0.473770 - 0.0066356 * e1 / (1. + 3.25434 * s * e1);
  ten_alpha2_     = pow(10., -1.7329 - 0.459096711 * e2 / (1. + 8.92332 * s * e2));
  alpha3_         = pow(10.,
                        -1.23587 * pow(Gamma_eff, -0.2151) - 1.87406e-6 * pow(Gamma_eff, -1.52905) -
                            1.20918 * pow(Gamma_eff, -0.2182));
  alpha4_         = -0.011777 + 0.15448 * s / (1. + 0.738962 * s) +
                    0.04227 * exp(-pow((4. * s - 1.) / 16.1314, 2));
}

void MajoronSpecies::Fit(double x, double* F, double* dF, double* d2F) const {
  const double E   = exp(-(x - 1.) / alpha4_);
  const double D   = 1. + alpha3_ * (1. + x) * E;
  const double dD  = alpha3_ * E * (1. - (1. + x) / alpha4_);
  const double d2D = alpha3_ * E * ((1. + x) / alpha4_ - 2.) / alpha4_;
  *F               = alpha1_ - ten_alpha2_ / D;
  *dF              = ten_alpha2_ * dD / (D * D);
  *d2F             = ten_alpha2_ * (d2D - 2. * dD * dD / D) / (D * D);
}

double MajoronSpecies::GetOmega0() const {
  double F, dF, d2F;
  Fit(x0_, &F, &dF, &d2F);
  return Omega_nu_per_F_ * F;
}

double MajoronSpecies::GetRadiationOmega0() const {
  /* Before the majorons, x -> infinity. */
  return Omega_nu_per_F_ * (alpha1_ - ten_alpha2_);
}

double MajoronSpecies::BackgroundDensityOverH0Sq(double a, double /*H0*/) const {
  double F, dF, d2F;
  Fit(x0_ / a, &F, &dF, &d2F);
  return Omega_nu_per_F_ * F / (a * a * a * a);
}

// ── Background ─────────────────────────────────────────────────────────────────

void MajoronSpecies::RegisterBackgroundIndices(int& index_bg) {
  UltraRelativisticSpecies::RegisterBackgroundIndices(index_bg);
  class_define_index(index_bg_p_, true, index_bg, 1);
}

void MajoronSpecies::ComputeBackground(double a, const double* /*pvecback_B*/, double* pvecback) {
  /* rho = rho_gamma F(x) / (pi^2/15); continuity then gives p = (rho/3) (1 + x F'/F). */
  const double x = x0_ / a;
  double F, dF, d2F;
  Fit(x, &F, &dF, &d2F);
  const double rho        = Omega_nu_per_F_ * H0_ * H0_ / (a * a * a * a) * F;
  pvecback[index_bg_rho_] = rho;
  pvecback[index_bg_p_]   = rho / 3. * (1. + x * dF / F);
}

double MajoronSpecies::P(const double* pvecback) const {
  return pvecback[index_bg_p_];
}

double MajoronSpecies::PPrime(double a,
                              double H,
                              const double* /*pvecback_B*/,
                              const double* pvecback) const {
  /* With G = x F'/F and x = x0/a: d ln rho/d ln a = -4 - G and dG/d ln a = -x G', so
     dp/d ln a = (rho/3) [-(4 + G)(1 + G) - x G']. */
  const double x = x0_ / a;
  double F, dF, d2F;
  Fit(x, &F, &dF, &d2F);
  const double G     = x * dF / F;
  const double dG_dx = dF / F + x * d2F / F - x * (dF / F) * (dF / F);
  const double rho   = pvecback[index_bg_rho_];
  return a * H * rho / 3. * (-(4. + G) * (1. + G) - x * dG_dx);
}

// ── Perturbations ──────────────────────────────────────────────────────────────

double MajoronSpecies::CollisionRate(double a) const {
  /* x = m_phi/T_nu, with T_nu = (4/11)^{1/3} T_gamma as Escudero & Witte take it. */
  const double x     = a / (pow(4. / 11., 1. / 3.) * x0_);
  const double Gamma = Gamma_phi_ * x * x * x * BesselK2(x);  // eV
  return a * Gamma * _eV_ * 2. * _PI_ / _h_P_ / _c_ * _Mpc_over_m_;
}

void MajoronSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                   double tau,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) const {
  UltraRelativisticSpecies::PerturbDerivs(base, tau, y, dy, ppaw);
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_delta < 0)  // radiation streaming approximation
    return;
  /* dF_l/dtau -= a Gamma F_l for l >= 2, with F_2 = 2 shear; under the ur fluid
     approximation l_max is 2 and only the shear is left. */
  const double rate     = CollisionRate(ppaw.ppw->scalar_ctx.a);
  dy[layout.idx_shear] -= rate * y[layout.idx_shear];
  for (int l = 3; l <= layout.l_max; ++l)
    dy[layout.idx_delta + l] -= rate * y[layout.idx_delta + l];
}

// ── Factory ────────────────────────────────────────────────────────────────────

std::vector<Named> MajoronSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  const auto m_phi     = ctx.pfc->get<double>("m_majoron");
  const auto Gamma_eff = ctx.pfc->get<double>("Gamma_eff_majoron");
  if (!m_phi && !Gamma_eff)
    return result;
  class_test_severe(!m_phi || !Gamma_eff,
                    "the majoron needs both 'm_majoron' (in eV) and 'Gamma_eff_majoron'");
  class_test(*m_phi <= 0., "m_majoron = %g eV must be positive", *m_phi);
  class_test(*Gamma_eff <= 0., "Gamma_eff_majoron = %g must be positive", *Gamma_eff);
  /* The fluid carries the three neutrinos, so the default N_ur = 3.046 would count them
     twice: N_ur is only the extra free-streaming radiation, and must be given. */
  class_test_severe(!ctx.pfc->get<double>("N_ur") && !ctx.pfc->get<double>("N_eff") &&
                        !ctx.pfc->get<double>("Omega_ur") && !ctx.pfc->get<double>("omega_ur"),
                    "with the majoron, which carries the three neutrinos, set 'N_ur' to the "
                    "extra free-streaming radiation (0 for none)");

  result.push_back({"Majoron", std::make_unique<MajoronSpecies>(*ctx.pba, *m_phi, *Gamma_eff)});
  return result;
}
