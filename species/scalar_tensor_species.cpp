#include "scalar_tensor_species.h"

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>

#include "background.h"
#include "background_module.h"
#include "constants.h"
#include "errors.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"

namespace {

/* The reduced Planck mass sqrt(hbar c / 8 pi G) as an inverse length, in 1/Mpc: the
   factor that turns the action's dimensionless lambda into CLASS units. */
double ReducedPlanckMassInverseMpc() {
  const double hbar = _h_P_ / (2. * _PI_);
  return sqrt(pow(_c_, 3) / (8. * _PI_ * _G_ * hbar)) * _Mpc_over_m_;
}

}  // namespace

ScalarTensorSpecies::ScalarTensorSpecies(const background& pba,
                                         double xi,
                                         double N,
                                         double lambda,
                                         double V_Lambda,
                                         double phi_ini,
                                         double phi_prime_ini)
    : BaseSpecies("ScalarTensor", EnergyType::DarkEnergy), pba_(pba), xi_(xi), N_(N),
      lambda_class_(lambda * pow(ReducedPlanckMassInverseMpc(), 2)), V_Lambda_(V_Lambda),
      phi_ini_(phi_ini), phi_prime_ini_(phi_prime_ini) {}

double ScalarTensorSpecies::V(double phi) const {
  return V_Lambda_ + 0.25 * lambda_class_ * pow(phi, 4);
}
double ScalarTensorSpecies::dV(double phi) const {
  return lambda_class_ * pow(phi, 3);
}
double ScalarTensorSpecies::ddV(double phi) const {
  return 3. * lambda_class_ * phi * phi;
}

// ── Budget, shooting, derived parameters ──────────────────────────────────────

double ScalarTensorSpecies::BackgroundDensityOverH0Sq(double /*a*/, double /*H0*/) const {
  class_stop_severe(
      "the scalar-tensor species has no density before the background is solved, since it "
      "changes H itself; a species that locates an event from the densities in advance (NEDE) "
      "cannot be combined with it");
}

std::optional<double> ScalarTensorSpecies::GetParam(const std::string& name) const {
  if (name == "G_eff_today")
    return G_eff_today_;
  if (name == "gamma_PPN_minus_1")
    return gamma_PPN_minus_1_;
  return std::nullopt;
}

std::vector<ShootingTarget> ScalarTensorSpecies::GetShootingTargets() const {
  if (needs_shooting_)
    return {{"Omega_scalar_tensor", "Omega_Lambda_st", 0.}};
  return {};
}

void ScalarTensorSpecies::ComputeShootingGuess(const SpeciesBuildContext& /*ctx*/,
                                               std::vector<double>& guess,
                                               std::vector<double>& dxdy) const {
  /* The unknown is the potential's constant in units of 3 H0^2, where Newton's step
     tolerance (1e-3) means something: in 1/Mpc^2 it declared convergence after one
     step. The effective density responds to it as 1/f; f at phi_ini stands in for f
     today. */
  guess.push_back(0.);
  dxdy.push_back(F(phi_ini_));
}

double ScalarTensorSpecies::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                    const ShootingTarget& target) const {
  return Rho(ctx.bg_today) / (ctx.pba->H0 * ctx.pba->H0) - target.target_value;
}

// ── Background ────────────────────────────────────────────────────────────────

void ScalarTensorSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_, true, index_bg, 1);
  class_define_index(index_bg_p_, true, index_bg, 1);
  class_define_index(index_bg_phi_, true, index_bg, 1);
  class_define_index(index_bg_phi_prime_, true, index_bg, 1);
  class_define_index(index_bg_phi_prime_prime_, true, index_bg, 1);
  class_define_index(index_bg_trace_, true, index_bg, 1);
  class_define_index(index_bg_f_, true, index_bg, 1);
}

void ScalarTensorSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_phi_, true, index_bi, 1);
  class_define_index(index_bi_phi_prime_, true, index_bi, 1);
}

void ScalarTensorSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  ctx.pvecback_integration[index_bi_phi_]       = phi_ini_;
  ctx.pvecback_integration[index_bi_phi_prime_] = phi_prime_ini_;
}

void ScalarTensorSpecies::ComputeBackground(double /*a*/,
                                            const double* pvecback_B,
                                            double* pvecback) {
  const double phi              = pvecback_B[index_bi_phi_];
  pvecback[index_bg_phi_]       = phi;
  pvecback[index_bg_phi_prime_] = pvecback_B[index_bi_phi_prime_];
  pvecback[index_bg_f_]         = F(phi);
  /* The effective density and pressure need H, and H needs every other species:
     CloseFriedmann writes them. Until then this species adds nothing to the totals. */
  pvecback[index_bg_rho_] = 0.;
  pvecback[index_bg_p_]   = 0.;
}

void ScalarTensorSpecies::CloseFriedmann(double a,
                                         double K,
                                         const double* pvecback_B,
                                         double* pvecback,
                                         double& rho_tot,
                                         double& p_tot,
                                         double& rho_r,
                                         double& rho_m) const {
  /* Jordan frame, M_pl = 1, densities in CLASS units (rho = rho_phys / 3), dots d/dt:
       (00)  f (H^2 + K/a^2) + f_dot H = rho_bare
       (ij)  -f (2 H_dot + 3 H^2 + K/a^2) = 3 p_bare + f_ddot + 2 H f_dot
     where the bare totals are every other species plus the canonical scalar. */
  const double phi       = pvecback_B[index_bi_phi_];
  const double phi_prime = pvecback_B[index_bi_phi_prime_];
  const double f         = F(phi);
  const double f_phi     = 2. * xi_ * phi;
  const double f_phiphi  = 2. * xi_;
  class_test(f <= 0.,
             "the effective Planck mass vanished: f(phi) = %g at phi = %g, a = %g",
             f,
             phi,
             a);

  const double phi_dot  = phi_prime / a;
  const double kinetic  = 0.5 * phi_dot * phi_dot;
  const double rho_bare = rho_tot + (kinetic + V(phi)) / 3.;
  const double p_bare   = p_tot + (kinetic - V(phi)) / 3.;
  const double f_dot    = f_phi * phi_dot;

  /* (00) is a quadratic in H; its positive root, in the form without cancellation.
     Without one, GR's sqrt(rho_tot) below would silently return |H|. */
  const double c     = rho_bare - f * K / (a * a);
  const double disc2 = f_dot * f_dot + 4. * f * c;
  class_test(disc2 < 0. || (f_dot >= 0. && c <= 0.),
             "the Jordan-frame Friedmann equation has no expanding solution at a = %g "
             "(bare density %g, f = %g, f_dot = %g)",
             a,
             rho_bare,
             f,
             f_dot);
  const double disc = sqrt(disc2);
  const double H    = (f_dot >= 0.) ? 2. * c / (f_dot + disc) : (disc - f_dot) / (2. * f);

  /* Klein-Gordon, box phi = V_phi - f_phi R/2, with R eliminated by the trace of the
     field equations, R = (3 box f - T)/f:  C box phi = S. box phi = -(phi_ddot + 3 H
     phi_dot) on the background, (d phi)^2 = -phi_dot^2, T is the total trace. */
  const double trace    = 3. * (-rho_bare + 3. * p_bare);
  const double C        = 1. + 1.5 * f_phi * f_phi / f;
  const double S        = dV(phi) - f_phi / (2. * f) * (3. * f_phiphi * (-2. * kinetic) - trace);
  const double phi_ddot = -S / C - 3. * H * phi_dot;
  const double f_ddot   = f_phi * phi_ddot + f_phiphi * phi_dot * phi_dot;
  const double H_dot = -(3. * p_bare + f_ddot + 2. * H * f_dot + f * (3. * H * H + K / (a * a))) /
                       (2. * f);

  pvecback[index_bg_phi_prime_prime_] = a * a * (phi_ddot + H * phi_dot);  // a' = a^2 H
  pvecback[index_bg_trace_]           = trace;

  /* The effective fluid: whatever GR's H^2 = rho_tot - K/a^2 and
     H_dot = -3/2 (rho_tot + p_tot) + K/a^2 need beyond the other species. */
  const double rho_eff     = H * H + K / (a * a) - rho_tot;
  const double p_eff       = -(2. * H_dot + 3. * H * H + K / (a * a)) / 3. - p_tot;
  pvecback[index_bg_rho_]  = rho_eff;
  pvecback[index_bg_p_]    = p_eff;
  rho_tot                 += rho_eff;
  p_tot                   += p_eff;
  /* Radiation and matter drive the expansion with G/f. */
  rho_r /= f;
  rho_m /= f;
}

void ScalarTensorSpecies::BackgroundDerivs(double /*tau*/,
                                           const double* y,
                                           double* dy,
                                           const double* pvecback) {
  dy[index_bi_phi_]       = y[index_bi_phi_prime_];
  dy[index_bi_phi_prime_] = pvecback[index_bg_phi_prime_prime_];
}

double ScalarTensorSpecies::PPrime(double a,
                                   double H,
                                   const double* pvecback_B,
                                   const double* pvecback) const {
  /* p_phi = (phi'^2/(2a^2) - V)/3, and (1/a^2)' = -2 aH / a^2. */
  const double phi       = pvecback_B[index_bi_phi_];
  const double phi_prime = pvecback_B[index_bi_phi_prime_];
  const double a_H       = a * H;
  return (phi_prime * (pvecback[index_bg_phi_prime_prime_] - a_H * phi_prime) / (a * a) -
          dV(phi) * phi_prime) /
         3.;
}

double ScalarTensorSpecies::SourceSamplingRate(const double* pvecback) const {
  /* phi'' = -lambda phi^3 (per unit cosmic time) oscillates with period 7.4163/(sqrt(lambda) A)
     at amplitude A, and the energy lambda A^4/4 = phi_dot^2/2 + lambda phi^4/4 gives A.
     Twice the frequency: twenty source samples per period. */
  if (lambda_class_ <= 0.)
    return 0.;
  const double a       = pvecback[bgm_->index_bg_a_];
  const double phi     = pvecback[index_bg_phi_];
  const double phi_dot = pvecback[index_bg_phi_prime_] / a;
  const double energy  = 0.5 * phi_dot * phi_dot + 0.25 * lambda_class_ * pow(phi, 4);
  const double A       = pow(4. * energy / lambda_class_, 0.25);
  return 2. * a * sqrt(lambda_class_) * A / 7.4163;
}

void ScalarTensorSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_st", 0.);
  w.Add("(.)p_st", 0.);
  w.Add("phi_st", 0.);
  w.Add("phi'_st", 0.);
  w.Add("f_st", 0.);
}

void ScalarTensorSpecies::WriteBackgroundData(const double* pvecback,
                                              BackgroundColumnWriter& w) const {
  w.Add("(.)rho_st", pvecback[index_bg_rho_]);
  w.Add("(.)p_st", pvecback[index_bg_p_]);
  w.Add("phi_st", pvecback[index_bg_phi_]);
  w.Add("phi'_st", pvecback[index_bg_phi_prime_]);
  w.Add("f_st", pvecback[index_bg_f_]);
}

void ScalarTensorSpecies::ProcessBackgroundTable(const double* background_table,
                                                 int n_rows,
                                                 int row_stride,
                                                 const double* /*z_table*/) {
  /* Today's Newton constant in a Cavendish experiment, and gamma_PPN - 1, for a
     massless field (Boisseau et al., gr-qc/0001066), with G the constant of the action. */
  const double phi   = background_table[(n_rows - 1) * row_stride + index_bg_phi_];
  const double f     = F(phi);
  const double f_phi = 2. * xi_ * phi;
  const double f2    = f_phi * f_phi;
  G_eff_today_       = (2. * f + 4. * f2) / (f * (2. * f + 3. * f2));
  gamma_PPN_minus_1_ = -f2 / (f + 2. * f2);
}

// ── Perturbations ─────────────────────────────────────────────────────────────
//
// Synchronous gauge (Ma & Bertschinger), M_pl = 1, perturbations of the CLASS-unit
// totals delta_rho, (rho+p)theta, delta_p, (rho+p)sigma, which include the canonical
// scalar. Perturbing G_mu^nu = [T_mu^nu + nabla^nu nabla_mu f - delta_mu^nu box f]/f,
// the four Einstein equations gain the terms in delta f = f_phi delta phi below, and
// reduce to CLASS's GR lines at f = 1, delta f = 0. Derivation in the design spec.

ScalarTensorSpecies::Coupling ScalarTensorSpecies::CouplingAt(const double* y,
                                                              const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(
      *ppw->pv->species_layouts[collection_index_]);
  const double* pvecback = ppw->pvecback.data();
  const double phi       = pvecback[index_bg_phi_];
  const double phi_prime = pvecback[index_bg_phi_prime_];
  const double f_phi     = 2. * xi_ * phi;
  const double f_phiphi  = 2. * xi_;

  Coupling c;
  c.a          = pvecback[bgm_->index_bg_a_];
  c.aH         = c.a * pvecback[bgm_->index_bg_H_];
  c.f          = F(phi);
  c.f_prime    = f_phi * phi_prime;
  c.dphi       = y[layout.idx_phi];
  c.dphi_prime = y[layout.idx_phi_prime];
  c.df         = f_phi * c.dphi;
  c.df_prime   = f_phi * c.dphi_prime + f_phiphi * phi_prime * c.dphi;
  return c;
}

double ScalarTensorSpecies::DeltaPhiPrimePrime(double k,
                                               const double* y,
                                               const perturb_workspace* ppw) const {
  /* Perturbed Klein-Gordon in trace form, C box phi = S (see CloseFriedmann), with
       box phi -> -[delta phi'' + 2 aH delta phi' + k^2 delta phi + h' phi'/2] / a^2,
       (d phi)^2 = -phi'^2/a^2 -> -2 phi' delta phi'/a^2, delta T = 3 (-delta rho + 3 delta p).
     It needs only h' and the totals, both ready before the species' derivatives run. */
  const Coupling c       = CouplingAt(y, ppw);
  const double* pvecback = ppw->pvecback.data();
  const double phi       = pvecback[index_bg_phi_];
  const double phi_prime = pvecback[index_bg_phi_prime_];
  const double a2        = c.a * c.a;
  const double f_phi     = 2. * xi_ * phi;
  const double f_phiphi  = 2. * xi_;
  const double f         = c.f;

  const double C       = 1. + 1.5 * f_phi * f_phi / f;
  const double C_phi   = 3. * f_phi * f_phiphi / f - 1.5 * pow(f_phi, 3) / (f * f);
  const double X2      = -phi_prime * phi_prime / a2;
  const double dX2     = -2. * phi_prime * c.dphi_prime / a2;
  const double T       = pvecback[index_bg_trace_];
  const double dT      = 3. * (-ppw->delta_rho + 3. * ppw->delta_p);
  const double S       = dV(phi) - f_phi / (2. * f) * (3. * f_phiphi * X2 - T);
  const double dS      = ddV(phi) * c.dphi -
                         (f_phiphi / (2. * f) - f_phi * f_phi / (2. * f * f)) * c.dphi *
                             (3. * f_phiphi * X2 - T) -
                         f_phi / (2. * f) * (3. * f_phiphi * dX2 - dT);
  const double h_prime = ppw->pvecmetric[ppw->index_mt_h_prime];

  return -2. * c.aH * c.dphi_prime - k * k * c.dphi - 0.5 * h_prime * phi_prime -
         a2 / C * (dS - C_phi * c.dphi * S / C);
}

double ScalarTensorSpecies::HPrime(double k, const double* y, const perturb_workspace* ppw) const {
  /* 00: k^2 eta - aH h'/2 = [-1.5 a^2 delta rho + (3 aH df' + k^2 df + 3 aH^2 df + f' h'/2)/2]/f */
  const Coupling c = CouplingAt(y, ppw);
  const double k2  = k * k;
  const double eta = y[ppw->pv->index_pt_eta];
  return (k2 * eta + (1.5 * c.a * c.a * ppw->delta_rho -
                      0.5 * (3. * c.aH * c.df_prime + k2 * c.df + 3. * c.aH * c.aH * c.df)) /
                         c.f) /
         (0.5 * c.aH + 0.25 * c.f_prime / c.f);
}

double ScalarTensorSpecies::EtaPrime(double k,
                                     const double* y,
                                     const perturb_workspace* ppw) const {
  /* 0i: k^2 eta' = [1.5 a^2 (rho+p) theta + k^2 (df' - aH df)/2]/f */
  const Coupling c = CouplingAt(y, ppw);
  const double k2  = k * k;
  return (1.5 * c.a * c.a * ppw->rho_plus_p_theta + 0.5 * k2 * (c.df_prime - c.aH * c.df)) /
         (c.f * k2);
}

double ScalarTensorSpecies::HPrimePrime(double k,
                                        const double* y,
                                        const perturb_workspace* ppw) const {
  /* ii: h'' = -2 aH h' + 2 k^2 eta
               - [9 a^2 delta p + 3 df'' + 3 aH df' + 2 k^2 df + f' h' + 3 (2 aH' + aH^2) df]/f,
     with (aH)' = a^2 H^2 + a H'. */
  const Coupling c       = CouplingAt(y, ppw);
  const double* pvecback = ppw->pvecback.data();
  const double k2        = k * k;
  const double eta       = y[ppw->pv->index_pt_eta];
  const double phi       = pvecback[index_bg_phi_];
  const double phi_prime = pvecback[index_bg_phi_prime_];
  const double f_phi     = 2. * xi_ * phi;
  const double f_phiphi  = 2. * xi_;
  const double aH_prime  = c.aH * c.aH + c.a * pvecback[bgm_->index_bg_H_prime_];
  const double h_prime   = ppw->pvecmetric[ppw->index_mt_h_prime];
  const double ddf       = f_phi * DeltaPhiPrimePrime(k, y, ppw) +
                           2. * f_phiphi * phi_prime * c.dphi_prime +
                           f_phiphi * pvecback[index_bg_phi_prime_prime_] * c.dphi;
  return -2. * c.aH * h_prime + 2. * k2 * eta -
         (9. * c.a * c.a * ppw->delta_p + 3. * ddf + 3. * c.aH * c.df_prime + 2. * k2 * c.df +
          c.f_prime * h_prime + 3. * (2. * aH_prime + c.aH * c.aH) * c.df) /
             c.f;
}

double ScalarTensorSpecies::AlphaPrime(double k,
                                       const double* y,
                                       const perturb_workspace* ppw) const {
  /* Traceless ij: alpha' = -2 aH alpha + eta - [4.5 (a^2/k^2)(rho+p) sigma + df + f' alpha]/f.
     In Newtonian gauge this is the scalar's gravitational slip, psi - phi = -df/f. */
  const Coupling c   = CouplingAt(y, ppw);
  const double eta   = y[ppw->pv->index_pt_eta];
  const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
  const double a2_k2 = c.a * c.a / (k * k);
  return -2. * c.aH * alpha + eta -
         (4.5 * a2_k2 * ppw->rho_plus_p_shear + c.df + c.f_prime * alpha) / c.f;
}

void ScalarTensorSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                      perturb_vector* /*pv*/,
                                                      const precision* /*ppr*/,
                                                      int& index_pt,
                                                      const perturb_workspace* /*ppw*/,
                                                      int /*gauge*/) {
  auto& layout         = static_cast<PerturbLayout&>(base);
  layout.idx_phi       = index_pt++;
  layout.idx_phi_prime = index_pt++;
}

void ScalarTensorSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                        double /*tau*/,
                                        const double* y,
                                        double* dy,
                                        const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout           = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw = ppaw.ppw;
  dy[layout.idx_phi]           = y[layout.idx_phi_prime];
  dy[layout.idx_phi_prime]     = DeltaPhiPrimePrime(ppw->scalar_ctx.k, y, ppw);
}

BaseSpecies::StressEnergyContribution ScalarTensorSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  /* The canonical scalar, in CLASS units; the coupling acts through the Einstein
     equations, not here. */
  const auto& layout      = static_cast<const PerturbLayout&>(base);
  const double phi        = pvecback[index_bg_phi_];
  const double phi_prime  = pvecback[index_bg_phi_prime_];
  const double a2         = ppw->scalar_ctx.a2;
  const double dphi       = y[layout.idx_phi];
  const double dphi_prime = y[layout.idx_phi_prime];
  const double kinetic    = 0.5 * phi_prime * phi_prime / a2;

  StressEnergyContribution se;
  se.rho              = (kinetic + V(phi)) / 3.;
  se.p                = (kinetic - V(phi)) / 3.;
  se.delta_rho        = (phi_prime * dphi_prime / a2 + dV(phi) * dphi) / 3.;
  se.rho_plus_p_theta = ppw->scalar_ctx.k2 * phi_prime * dphi / (3. * a2);
  se.delta_p          = (phi_prime * dphi_prime / a2 - dV(phi) * dphi) / 3.;
  return se;
}

void ScalarTensorSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                 double* y,
                                                 const PerturbIcContext& /*ctx*/) {
  /* Adiabatic: the field is frozen by Hubble friction, so delta phi = -phi' delta tau
     vanishes at leading order. */
  const auto& layout      = static_cast<const PerturbLayout&>(base);
  y[layout.idx_phi]       = 0.;
  y[layout.idx_phi_prime] = 0.;
}

void ScalarTensorSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                        const BaseSpecies::PerturbLayout& new_base,
                                                        const double* old_y,
                                                        double* new_y,
                                                        const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_layout          = static_cast<const PerturbLayout&>(old_base);
  const auto& new_layout          = static_cast<const PerturbLayout&>(new_base);
  new_y[new_layout.idx_phi]       = old_y[old_layout.idx_phi];
  new_y[new_layout.idx_phi_prime] = old_y[old_layout.idx_phi_prime];
}

void ScalarTensorSpecies::PrintVariables(PerturbColumnWriter& w,
                                         const BaseSpecies::PerturbLayout* base,
                                         double /*tau*/,
                                         const double* y,
                                         const PerturbationsModule& /*mod*/,
                                         const perturb_workspace* /*ppw*/) const {
  double dphi = 0., dphi_prime = 0.;
  if (!w.IsTitleMode()) {
    const auto& layout = static_cast<const PerturbLayout&>(*base);
    dphi               = y[layout.idx_phi];
    dphi_prime         = y[layout.idx_phi_prime];
  }
  w.Add("delta_phi_st", dphi, true);
  w.Add("delta_phi_prime_st", dphi_prime, true);
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> ScalarTensorSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  const auto xi = ctx.pfc->get<double>("xi_st");
  if (!xi)
    return result;

  const auto phi_ini = ctx.pfc->get<double>("phi_ini_st");
  class_test_severe(!phi_ini,
                    "the scalar-tensor species (xi_st) needs 'phi_ini_st', the field at the "
                    "start of the integration in units of the reduced Planck mass");
  const double N             = ctx.pfc->get_or("N_st", 1.);
  const double lambda        = ctx.pfc->get_or("lambda_st", 0.);
  const double phi_prime_ini = ctx.pfc->get_or("phi_prime_ini_st", 0.);
  const auto Omega_Lambda    = ctx.pfc->get<double>("Omega_Lambda_st");

  /* Implemented for synchronous gauge, scalar modes, adiabatic initial conditions and a
     flat universe. Anything else would run GR's equations or unchecked ones. (N-body-
     gauge sources are refused by the perturbations module, which knows every route.) */
  auto key_has = [&](const char* key, std::initializer_list<const char*> needles) {
    const auto value = ctx.pfc->get<std::string>(key);
    if (!value)
      return false;
    for (const char* needle : needles)
      if (value->find(needle) != std::string::npos)
        return true;
    return false;
  };
  class_test_severe(key_has("gauge", {"new", "New"}),
                    "the scalar-tensor species is implemented in synchronous gauge only");
  class_test_severe(key_has("modes", {"t", "T", "v", "V"}),
                    "the scalar-tensor species is implemented for scalar modes only (its "
                    "tensor modes would need the friction f'/f)");
  class_test_severe(key_has("ic", {"bi", "BI", "cdi", "CDI", "nid", "NID", "niv", "NIV"}),
                    "the scalar-tensor species is implemented for adiabatic initial "
                    "conditions only");
  /* Severe although Omega_k is numeric: this is a structural limit, and as a computation
     error a sampler varying Omega_k would reject every point instead of reporting it. */
  class_test_severe(ctx.pba->Omega0_k != 0.,
                    "the scalar-tensor species is implemented for a flat universe only");

  auto species = std::make_unique<ScalarTensorSpecies>(*ctx.pba,
                                                       *xi,
                                                       N,
                                                       lambda,
                                                       3. * ctx.pba->H0 * ctx.pba->H0 *
                                                           Omega_Lambda.value_or(0.),
                                                       *phi_ini,
                                                       phi_prime_ini);
  /* Shoot the potential's constant unless given: DoShooting writes the resolved value
     back as Omega_Lambda_st for the final build. */
  species->needs_shooting_ = !Omega_Lambda.has_value();
  result.push_back({"ScalarTensor", std::move(species)});
  return result;
}
