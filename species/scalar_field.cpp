#include "scalar_field.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

#include "background_module.h"
#include "errors.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "species_input.h"

ScalarFieldSpecies::ScalarFieldSpecies(const background& pba,
                                       double omega0_scf,
                                       std::vector<double> scf_parameters,
                                       int scf_tuning_index,
                                       bool attractor_ic_scf,
                                       double phi_ini_scf,
                                       double phi_prime_ini_scf,
                                       ScalarFieldPotential potential,
                                       double beta,
                                       std::string instance_name)
    : BaseSpecies(instance_name.empty() ? "ScalarField" : instance_name, EnergyType::Other),
      pba_(pba), Omega0_scf_(omega0_scf), scf_parameters_(std::move(scf_parameters)),
      scf_tuning_index_(scf_tuning_index), attractor_ic_scf_(attractor_ic_scf),
      phi_ini_scf_(phi_ini_scf), phi_prime_ini_scf_(phi_prime_ini_scf),
      potential_(std::move(potential)), beta_(beta),
      column_tag_(instance_name.empty() ? "scf" : std::move(instance_name)) {}

ScalarFieldSpecies::InputKeys ScalarFieldSpecies::InputKeys::Legacy() {
  return {/*instance_name=*/"",
          "Omega_scf",
          "scf_potential",
          "scf_parameters",
          "scf_tuning_index",
          "scf_shooting_parameter",
          "attractor_ic_scf",
          "m_axion",
          "f_axion",
          "n_axion",
          "Theta_initial_scf"};
}

ScalarFieldSpecies::InputKeys ScalarFieldSpecies::InputKeys::Instance(const std::string& name) {
  // The dot fields drop the legacy key's _scf / _axion suffix, as the ncdm fields do.
  const std::string p = name + ".";
  return {name,
          p + "Omega",
          p + "potential",
          p + "parameters",
          p + "tuning_index",
          p + "shooting_parameter",
          p + "attractor_ic",
          p + "m",
          p + "f",
          p + "n",
          p + "Theta_initial"};
}

double ScalarFieldSpecies::V_scf(double phi) const {
  return potential_.V(phi, scf_parameters_);
}
double ScalarFieldSpecies::dV_scf(double phi) const {
  return potential_.dV(phi, scf_parameters_);
}
double ScalarFieldSpecies::ddV_scf(double phi) const {
  return potential_.ddV(phi, scf_parameters_);
}

void ScalarFieldSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  /** - Fix initial value of \f$ \phi, \phi' \f$
   * set directly in the radiation attractor => fixes the units in terms of rho_ur
   *
   * TODO:
   * - There seems to be some small oscillation when it starts.
   * - Check equations and signs. Sign of phi_prime?
   * - is rho_ur all there is early on?
   */
  double* pvecback_integration = ctx.pvecback_integration;
  const double rho_rad         = ctx.rho_rad;
  double scf_lambda            = scf_parameters()[0];
  if (attractor_ic_scf()) {
    const double attractor_denom = 3 * pow(scf_lambda, 2) - 12;
    if (attractor_denom < 0) {
      /** - --> If there is no attractor solution for scf_lambda, assign some value. Otherwise would give a nan.*/
      pvecback_integration[bi_phi_index()] = 1. / scf_lambda;  //seems to the work
      if (pba_.background_verbose > 0) {
        printf(" No attractor IC for lambda = %.3e ! \n ", scf_lambda);
      }
    }
    else {
      pvecback_integration[bi_phi_index()] = -1. / scf_lambda *
                                             log(rho_rad * 4. / attractor_denom) * phi_ini_scf();
    }
    pvecback_integration[bi_phi_prime_index()] = 2 * ctx.a_ini *
                                                 sqrt(V_scf(pvecback_integration[bi_phi_index()])) *
                                                 phi_prime_ini_scf();
  }
  else {
    printf("Not using attractor initial conditions\n");
    /** - --> If no attractor initial conditions are assigned, gets the provided ones. */
    pvecback_integration[bi_phi_index()]       = phi_ini_scf();
    pvecback_integration[bi_phi_prime_index()] = phi_prime_ini_scf();
  }
  class_test(!isfinite(pvecback_integration[bi_phi_index()]) ||
                 !isfinite(pvecback_integration[bi_phi_prime_index()]),
             "initial phi = %e phi_prime = %e -> check initial conditions",
             pvecback_integration[bi_phi_index()],
             pvecback_integration[bi_phi_prime_index()]);
}

void ScalarFieldSpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  const double phi                  = pvecback_B[index_bi_phi_scf_];
  const double phi_prime            = pvecback_B[index_bi_phi_prime_scf_];
  pvecback[index_bg_phi_scf_]       = phi;
  pvecback[index_bg_phi_prime_scf_] = phi_prime;
  pvecback[index_bg_V_scf_]         = V_scf(phi);
  pvecback[index_bg_dV_scf_]        = dV_scf(phi);
  pvecback[index_bg_ddV_scf_]       = ddV_scf(phi);
  pvecback[index_bg_rho_] =
      ((1. - 2. * beta_) * phi_prime * phi_prime / (2. * a * a) + V_scf(phi)) / 3.;
  pvecback[index_bg_p_] = ((1. - 2. * beta_) * phi_prime * phi_prime / (2. * a * a) - V_scf(phi)) /
                          3.;
  pvecback[index_bg_p_prime_scf_] = 0.;
}

void ScalarFieldSpecies::BackgroundDerivs(double /*tau*/,
                                          const double* y,
                                          double* dy,
                                          const double* pvecback) {
  const double a         = pvecback[bgm_->index_bg_a_];
  const double H         = pvecback[bgm_->index_bg_H_];
  const double phi       = y[index_bi_phi_scf_];
  const double phi_prime = y[index_bi_phi_prime_scf_];
  /** phi'' + 2*a*H*phi' + a^2*dV = 0 */
  dy[index_bi_phi_scf_]       = phi_prime;
  dy[index_bi_phi_prime_scf_] = -a * (2. * H * phi_prime + a * dV_scf(phi) / (1. - 2. * beta_));
}

double ScalarFieldSpecies::PPrime(double a,
                                  double H,
                                  const double* pvecback_B,
                                  const double* /*pvecback*/) const {
  const double phi       = pvecback_B[index_bi_phi_scf_];
  const double phi_prime = pvecback_B[index_bi_phi_prime_scf_];
  return phi_prime * (-(1. - 2. * beta_) * phi_prime * H / a - 2. / 3. * dV_scf(phi));
}

void ScalarFieldSpecies::FinalizeBackground(double a,
                                            double H,
                                            const double* pvecback_B,
                                            double* pvecback) {
  pvecback[index_bg_p_prime_scf_] = PPrime(a, H, pvecback_B, pvecback);
}

void ScalarFieldSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Note(
      "The units of phi, tau in the derivatives and the potential V are the following:\n"
      "--> phi is given in units of the reduced Planck mass m_Pl = (8 pi G)^(-1/2)\n"
      "--> tau in the derivative of V(phi) is given in units of Mpc.\n"
      "--> the potential V(phi) is given in units of m_Pl^2/Mpc^2.");
  const std::string& t = column_tag_;
  w.Add("(.)rho_" + t, 0.);
  w.Add("(.)p_" + t, 0.);
  w.Add("(.)p_prime_" + t, 0.);
  w.Add("phi_" + t, 0.);
  w.Add("phi'_" + t, 0.);
  w.Add("V_" + t, 0.);
  w.Add("V'_" + t, 0.);
  w.Add("V''_" + t, 0.);
}

void ScalarFieldSpecies::WriteBackgroundData(const double* pvecback,
                                             BackgroundColumnWriter& w) const {
  const std::string& t = column_tag_;
  w.Add("(.)rho_" + t, Rho(pvecback));
  w.Add("(.)p_" + t, P(pvecback));
  w.Add("(.)p_prime_" + t, pvecback[index_bg_p_prime_scf_]);
  w.Add("phi_" + t, pvecback[index_bg_phi_scf_]);
  w.Add("phi'_" + t, pvecback[index_bg_phi_prime_scf_]);
  w.Add("V_" + t, pvecback[index_bg_V_scf_]);
  w.Add("V'_" + t, pvecback[index_bg_dV_scf_]);
  w.Add("V''_" + t, pvecback[index_bg_ddV_scf_]);
}

void ScalarFieldSpecies::PrintBackgroundDetails(const double* pvecback_today,
                                                double rho_crit_today) const {
  if (name() == "ScalarField") {
    printf("    Scalar field details:\n");
  }
  else {
    printf("    Scalar field '%s' details:\n", name().c_str());
  }
  printf("     -> Omega_%s = %g, wished %g\n",
         column_tag_.c_str(),
         Rho(pvecback_today) / rho_crit_today,
         GetOmega0());
  printf("     -> parameters: [%s] = \n", potential_.parameter_names.c_str());
  printf("                    [");
  for (size_t i = 0; i < scf_parameters_.size(); i++) {
    printf("%s%g", i == 0 ? "" : ", ", scf_parameters_[i]);
  }
  printf("]\n");
}

void ScalarFieldSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                       const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void ScalarFieldSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                     perturb_vector* /*pv*/,
                                                     const precision* /*ppr*/,
                                                     int& index_pt,
                                                     const perturb_workspace* /*ppw*/,
                                                     int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_phi = index_pt;
  ++index_pt;

  layout.idx_phi_prime = index_pt;
  ++index_pt;
}

void ScalarFieldSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                       double /*tau*/,
                                       const double* y,
                                       double* dy,
                                       const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  /* Use BackgroundModule's canonical pvecback slots for background quantities */
  const double phi_prime_bg = ppw->pvecback[index_bg_phi_prime_scf_];
  const double dV_bg        = ppw->pvecback[index_bg_dV_scf_];
  const double ddV_bg       = ppw->pvecback[index_bg_ddV_scf_];

  const double k2                = ctx.k2;
  const double a2                = ctx.a2;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double metric_continuity = ctx.metric_continuity;

  if (ctx.gauge == static_cast<int>(possible_gauges::synchronous)) {
    dy[layout.idx_phi]       = y[layout.idx_phi_prime];
    dy[layout.idx_phi_prime] = -2. * a_prime_over_a * y[layout.idx_phi_prime] -
                               metric_continuity * phi_prime_bg -
                               (k2 / (1. - 2. * beta_) + a2 * ddV_bg) * y[layout.idx_phi];
  }
  else {
    /* In Newtonian gauge, evolve q = delta_phi' - phi'_bg (psi + 3 phi)
       instead of delta_phi'. The background KG equation removes the otherwise
       required psi' term exactly:

         delta_phi' = q + phi'_bg (psi + 3 phi)
         q' = -2 H q - (k^2 + a^2 V_,phiphi) delta_phi
              + a^2 V_,phi (3 phi - psi).

       Thus free-streaming anisotropic stress enters only through the algebraic
       psi constraint; no shear time derivative is needed. */
    const double phi = y[ppw->pv->index_pt_phi];
    const double psi = ppw->pvecmetric[ppw->index_mt_psi];
    const double q   = y[layout.idx_phi_prime];

    dy[layout.idx_phi]       = q + phi_prime_bg * (psi + 3. * phi);
    dy[layout.idx_phi_prime] = -2. * a_prime_over_a * q - (k2 + a2 * ddV_bg) * y[layout.idx_phi] +
                               a2 * dV_bg * (3. * phi - psi);
  }
}

void ScalarFieldSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                     const double* y,
                                     const double* /*dy*/,
                                     PerturbSourceContext& ctx) const {
  const auto& layout         = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const double* pvecback     = ppw->pvecback.data();

  // Scalar field sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_) {
    return;
  }

  const double phi_prime_bg = pvecback[index_bg_phi_prime_scf_];
  const double dV_bg        = pvecback[index_bg_dV_scf_];
  const double rho_scf      = Rho(pvecback);
  const double p_scf        = P(pvecback);
  const perturbs* ppt       = p_mod->GetPerturbs();
  const double a2           = ctx.a2;
  const double k2 = ctx.k *
                    ctx.k;  // PerturbSourceContext has no k2 field (unlike PerturbScalarContext)
  double delta_phi_prime = y[layout.idx_phi_prime];
  if (ppt->gauge == possible_gauges::newtonian) {
    const double phi  = y[ppw->pv->index_pt_phi];
    const double psi  = ppw->pvecmetric[ppw->index_mt_psi];
    delta_phi_prime  += phi_prime_bg * (psi + 3. * phi);
  }

  // ── delta_scf ─────────────────────────────────────────────────────────────
  if (index_tp_delta_ >= 0) {
    double delta_rho_scf;
    if (ppt->gauge == possible_gauges::synchronous) {
      delta_rho_scf = 1. / 3. *
                          (1. / a2 * phi_prime_bg * delta_phi_prime + dV_bg * y[layout.idx_phi]) +
                      3. * ctx.a_prime_over_a * (1. + p_scf / rho_scf) *
                          ctx.theta_over_k2;  // N-body gauge correction
    }
    else {
      delta_rho_scf = 1. / 3. *
                          (1. / a2 * phi_prime_bg * delta_phi_prime + dV_bg * y[layout.idx_phi] -
                           1. / a2 * phi_prime_bg * phi_prime_bg *
                               ppw->pvecmetric[ppw->index_mt_psi]) +
                      3. * ctx.a_prime_over_a * (1. + p_scf / rho_scf) *
                          ctx.theta_over_k2;  // N-body gauge correction
    }
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_delta_,
                          ctx.index_tau,
                          ctx.index_k,
                          delta_rho_scf / rho_scf);
  }

  // ── theta_scf ─────────────────────────────────────────────────────────────
  if (index_tp_theta_ >= 0) {
    const double rho_plus_p_theta_scf = 1. / 3. * k2 / a2 * phi_prime_bg * y[layout.idx_phi];

    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_theta_,
                          ctx.index_tau,
                          ctx.index_k,
                          rho_plus_p_theta_scf / (rho_scf + p_scf) +
                              ctx.theta_shift);  // N-body gauge correction
  }
}

void ScalarFieldSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                double* y,
                                                const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_) {
    return;
  }
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_phi >= 0) {
    y[layout.idx_phi] = 0.;
  }
  if (layout.idx_phi_prime >= 0) {
    y[layout.idx_phi_prime] = 0.;
  }
}

void ScalarFieldSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                       double* y,
                                                       const PerturbIcContext& ctx) {
  const auto& l          = static_cast<const PerturbLayout&>(base);
  const double* pvecback = ctx.ppw->pvecback.data();
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double dV        = pvecback[index_bg_dV_scf_];
  if (l.idx_phi >= 0) {
    y[l.idx_phi] += ctx.alpha * phi_prime;
  }
  if (l.idx_phi_prime >= 0) {
    /* Transform the synchronous delta_phi' into the Newtonian q state. Using
       alpha' = psi - H alpha cancels psi explicitly and avoids requiring the
       initial total shear here. */
    const double phi    = y[ctx.ppw->pv->index_pt_phi];
    y[l.idx_phi_prime] += (-3. * ctx.a_prime_over_a * phi_prime - ctx.a * ctx.a * dV) * ctx.alpha -
                          3. * phi_prime * phi;
  }
}

void ScalarFieldSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                            const PerturbationsModule& mod,
                                            file_format fmt,
                                            BaseSpecies::TransferColumnSection section) const {
  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers) {
      w.Add("d_" + column_tag_, index_tp_delta_, index_tp_delta_ >= 0);
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers) {
      w.Add("t__" + column_tag_, index_tp_theta_, index_tp_theta_ >= 0);
    }
  }
}

void ScalarFieldSpecies::PrintVariables(PerturbColumnWriter& w,
                                        const BaseSpecies::PerturbLayout* base,
                                        double /*tau*/,
                                        const double* y,
                                        const PerturbationsModule& mod,
                                        const perturb_workspace* ppw) const {
  double delta_scf = 0.;
  double theta_scf = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const double* pvecmetric = ppw->pvecmetric.data();
    const double k           = ppw->scalar_ctx.k;
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];
    const double a2          = ppw->scalar_ctx.a2;
    const perturbs* ppt      = mod.GetPerturbs();

    // beta_ is the Type-3 momentum-transfer coupling and is 0 for a plain scalar
    // field, so this is the standard expression unless Type3Species built us.
    const double phi_prime_bg = pvecback[index_bg_phi_prime_scf_];
    const double dV_bg        = pvecback[index_bg_dV_scf_];
    const double rho_scf      = Rho(pvecback);
    const double p_scf        = P(pvecback);

    const auto& layout = static_cast<const PerturbLayout&>(*base);

    double delta_phi_prime = y[layout.idx_phi_prime];
    if (ppt->gauge == possible_gauges::newtonian) {
      const double phi  = y[pv->index_pt_phi];
      delta_phi_prime  += phi_prime_bg * (pvecmetric[ppw->index_mt_psi] + 3. * phi);
    }

    double delta_rho_scf = 0.;
    if (ppt->gauge == possible_gauges::synchronous) {
      delta_rho_scf = 1. / 3. *
                      ((1. - 2. * beta_) / a2 * phi_prime_bg * delta_phi_prime +
                       dV_bg * y[layout.idx_phi]);
    }
    else {
      delta_rho_scf = 1. / 3. *
                      ((1. - 2. * beta_) / a2 * phi_prime_bg * delta_phi_prime +
                       dV_bg * y[layout.idx_phi] -
                       1. / a2 * phi_prime_bg * phi_prime_bg * pvecmetric[ppw->index_mt_psi]);
    }

    const double rho_plus_p_theta_scf = 1. / 3. * ppw->scalar_ctx.k2 / a2 * phi_prime_bg *
                                        y[layout.idx_phi];

    delta_scf = delta_rho_scf / rho_scf;
    theta_scf = rho_plus_p_theta_scf / (rho_scf + p_scf);

    if (ppt->gauge == possible_gauges::synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double w_scf  = p_scf / rho_scf;
      delta_scf          += alpha * (-3. * H * (1. + w_scf));
      theta_scf          += k * k * alpha;
    }
  }

  w.Add("delta_" + column_tag_, delta_scf, true);
  w.Add("theta_" + column_tag_, theta_scf, true);
}

BaseSpecies::StressEnergyContribution ScalarFieldSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* pv,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;

  // rho, p: inline from Rho/P (both read pvecback slots)
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_p_];

  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double dV        = pvecback[index_bg_dV_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;
  const bool newtonian   = ppw->scalar_ctx.gauge == static_cast<int>(possible_gauges::newtonian);

  double psi             = 0.;
  double delta_phi_prime = y[layout.idx_phi_prime];
  if (newtonian) {
    const double phi  = y[pv->index_pt_phi];
    psi               = phi - 4.5 * (a2 / k2) * ppw->rho_plus_p_shear;
    delta_phi_prime  += phi_prime * (psi + 3. * phi);
  }

  // delta_rho: body of DeltaRho
  double delta_rho = (1. / 3.) * ((1. - 2. * beta_) / a2 * phi_prime * delta_phi_prime +
                                  dV * y[layout.idx_phi]);
  if (newtonian) {
    delta_rho -= (1. / 3.) * (1. / a2) * phi_prime * phi_prime * psi;
  }
  se.delta_rho = delta_rho;

  // rho_plus_p_theta: body of RhoPlusPTheta
  se.rho_plus_p_theta = (1. / 3.) * k2 / a2 * phi_prime * y[layout.idx_phi];

  // delta_p: body of DeltaP
  double delta_p = (1. / 3.) *
                   ((1. - 2. * beta_) / a2 * phi_prime * delta_phi_prime - dV * y[layout.idx_phi]);
  if (newtonian) {
    delta_p -= (1. / 3.) * (1. / a2) * phi_prime * phi_prime * psi;
  }
  se.delta_p = delta_p;

  // rho_plus_p_shear: always 0 for scalar field
  se.rho_plus_p_shear = 0.;

  return se;
}

void ScalarFieldSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                       const BaseSpecies::PerturbLayout& new_base,
                                                       const double* old_y,
                                                       double* new_y,
                                                       const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_phi < 0 || new_l.idx_phi < 0) {
    return;
  }
  new_y[new_l.idx_phi]       = old_y[old_l.idx_phi];
  new_y[new_l.idx_phi_prime] = old_y[old_l.idx_phi_prime];
}

// ── Shooting hooks ────────────────────────────────────────────────────────────

std::vector<ShootingTarget> ScalarFieldSpecies::GetShootingTargets() const {
  if (needs_shooting_) {
    return {shooting_target_};
  }
  return {};
}

void ScalarFieldSpecies::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                              std::vector<double>& guess,
                                              std::vector<double>& dxdy) const {
  // The guess is potential-specific, so it lives in the potential bundle. The
  // default bundle reproduces the historical lambda/Omega guess; the composite's
  // 1EXP bundle returns the frozen-field V0 guess (~3 H0^2 Omega).
  const auto [xguess, dxdF] =
      potential_.shooting_guess(Omega0_scf_, ctx.pba->H0, scf_parameters_, scf_tuning_index_);
  guess.push_back(xguess);
  dxdy.push_back(dxdF);
}

double ScalarFieldSpecies::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                   const ShootingTarget& target) const {
  // Omega_scf is filled to close the budget; the requested value is target.target_value.
  const background& ba = *ctx.pba;
  return Rho(ctx.bg_today) / (ba.H0 * ba.H0) - target.target_value;
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> ScalarFieldSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // ── Dot-syntax instances (#346): "<N>.type = scalar_field", one independent field
  // per instance, registered under N. They never close the budget, so they are built
  // in Pass 1 only (ConstructSpecies runs this factory in Pass 1 even when the legacy
  // field is the closure species); Pass 2, which carries the closure override, builds
  // the legacy field alone.
  if (!ctx.omega0_closure_override.has_value()) {
    for (const auto& name : ctx.pfc->instances_with("type", std::string(kTypeName))) {
      SpeciesInput input(ctx.pfc, name);
      (void) input.get<std::string>("type");  // consume "<N>.type"
      // "ScalarField" is the legacy field's key, "scf" its column tag.
      class_test_severe(name == "ScalarField" || name == "scf",
                        "'%s' is reserved for the scalar field configured by the legacy *_scf "
                        "keys; give the dot-syntax instance another name",
                        name.c_str());
      const InputKeys keys    = InputKeys::Instance(name);
      const double omega0_scf = input.require<double>("Omega");
      class_test(omega0_scf < 0.,
                 "%s = %g: a dot-syntax scalar field cannot close the budget; give it "
                 "Omega > 0, or close with the legacy Omega_scf < 0",
                 keys.Omega.c_str(),
                 omega0_scf);
      if (omega0_scf == 0.) {
        continue;  // as for the legacy field, Omega = 0 builds no field
      }
      for (auto& e : CreateConfigured(ctx, keys, omega0_scf, /*shooting_allowed=*/true)) {
        result.push_back(std::move(e));
      }
    }
  }

  // ── The legacy field (*_scf keys) ─────────────────────────────────────────
  // When the Type-3 (scf_veta) coupling is active, the Type3 composite owns the
  // legacy scalar field; the standalone factory must not also build one.
  if (auto veta = ctx.pfc->get<double>("scf_veta"); veta.has_value() && *veta != 0.) {
    return result;  // the Type3 composite owns the legacy scalar field
  }

  // Decide our Omega0_scf and whether shooting may apply.
  // Three cases:
  //   (i) Closure-Scf: Pass 2 of ConstructSpecies hands us the closure value via
  //       omega0_closure_override. Use it; no shooting.
  //  (ii) Non-closure, user gave Omega_scf > 0 in the input file → use it; shooting
  //       applies when scf_shooting_parameter is absent.
  // (iii) Otherwise (no override, Omega_scf absent or negative or zero) → no Scf species.
  //       (A negative Omega_scf in the input is the user's "scf is closure" signal,
  //       which is handled separately by ConstructSpecies via the override path.)

  double omega_scf_val   = 0.;
  auto omega_scf_opt     = ctx.pfc->get<double>("Omega_scf");
  bool omega_scf_present = omega_scf_opt.has_value();
  if (omega_scf_present) {
    omega_scf_val = *omega_scf_opt;
  }

  double omega0_scf     = 0.;
  bool shooting_allowed = false;

  if (ctx.omega0_closure_override.has_value()) {
    omega0_scf = *ctx.omega0_closure_override;
  }
  else if (omega_scf_present && omega_scf_val > 0.) {
    omega0_scf       = omega_scf_val;
    shooting_allowed = true;
  }
  else {
    return result;
  }

  if (omega0_scf == 0.) {
    return result;
  }

  for (auto& e : CreateConfigured(ctx, InputKeys::Legacy(), omega0_scf, shooting_allowed)) {
    result.push_back(std::move(e));
  }
  return result;
}

std::vector<Named> ScalarFieldSpecies::CreateConfigured(const SpeciesBuildContext& ctx,
                                                        const InputKeys& keys,
                                                        double omega0_scf,
                                                        bool shooting_allowed) {
  std::vector<Named> result;

  // ── potential dispatch ────────────────────────────────────────────────────
  if (auto pot_opt = ctx.pfc->get<std::string>(keys.potential)) {
    if (*pot_opt == "axion") {
      return CreateAxion(ctx, keys, omega0_scf);
    }
    class_stop_severe("unknown %s '%s': supported values are 'axion' (absent = historical default)",
                      keys.potential.c_str(),
                      pot_opt->c_str());
  }

  // ── parameters (variable-length list) ────────────────────────────────────
  std::vector<double> scf_parameters;
  try {
    if (auto scf_params_opt = ctx.pfc->get<std::vector<double>>(keys.parameters)) {
      scf_parameters = *scf_params_opt;
    }
  }
  catch (const std::exception& e) {
    class_stop_severe("%s parse error: %s", keys.parameters.c_str(), e.what());
  }

  // ── tuning_index + bounds check ──────────────────────────────────────────
  int scf_tuning_index = 0;
  scf_tuning_index     = ctx.pfc->get_or(keys.tuning_index, scf_tuning_index);
  class_test_severe(scf_tuning_index >= static_cast<int>(scf_parameters.size()),
                    "Tuning index %s = %d is larger than the number of entries "
                    "%zu in %s. Check your .ini file.",
                    keys.tuning_index.c_str(),
                    scf_tuning_index,
                    scf_parameters.size(),
                    keys.parameters.c_str());

  // shooting_parameter (when present) OVERWRITES the slot at tuning_index.
  // This is the resolved value DoShooting writes back into the fc for the final build.
  if (auto scf_shooting_opt = ctx.pfc->get<double>(keys.shooting_parameter)) {
    scf_parameters[scf_tuning_index] = *scf_shooting_opt;
  }

  // Preserve the legacy SCF parser warning: lambda<3 in the exponential-quintessence
  // potential won't track unless overwritten by the shooter / tuning function.
  if (!scf_parameters.empty()) {
    const double scf_lambda = scf_parameters[0];
    if (std::fabs(scf_lambda) < 3. && ctx.pba->background_verbose > 1) {
      printf(
          "lambda = %e <3 won't be tracking (for exp quint) unless overwritten "
          "by tuning function\n",
          scf_lambda);
    }
  }

  // ── attractor_ic (y/n string) ────────────────────────────────────────────
  const bool attractor_ic_scf = ctx.pfc->get_yes(keys.attractor_ic, true);
  double phi_ini_scf          = 1.;
  double phi_prime_ini_scf    = 1.;
  if (!attractor_ic_scf) {
    class_test_severe(scf_parameters.size() < 2,
                      "Since you are not using attractor initial conditions, you must specify "
                      "phi and its derivative phi' as the last two entries in %s. "
                      "See explanatory.ini for more details.",
                      keys.parameters.c_str());
    phi_ini_scf       = scf_parameters[scf_parameters.size() - 2];
    phi_prime_ini_scf = scf_parameters[scf_parameters.size() - 1];
  }

  auto species = std::make_unique<ScalarFieldSpecies>(*ctx.pba,
                                                      omega0_scf,
                                                      std::move(scf_parameters),
                                                      scf_tuning_index,
                                                      attractor_ic_scf,
                                                      phi_ini_scf,
                                                      phi_prime_ini_scf,
                                                      DefaultScalarFieldPotential(),
                                                      /*beta=*/0.,
                                                      keys.instance_name);

  // Detect guess-driven construction: Omega user-set to >0 but shooting_parameter
  // absent. The nonzero check matches the historic shooting condition (a zero target
  // meant "no shooting").
  if (shooting_allowed) {
    bool shooting_param_present = ctx.pfc->get<double>(keys.shooting_parameter).has_value();

    if (!shooting_param_present) {
      species->shooting_target_ = {keys.Omega, keys.shooting_parameter, omega0_scf};
      species->needs_shooting_  = true;

      std::vector<double> g;
      std::vector<double> d;
      species->ComputeShootingGuess(ctx, g, d);

      // Seed the guess into this discovery-build species' scf_parameters_[tuning_index]
      // so its ComputeBackground / shooter residual see a valid potential. Do NOT write it
      // into the file content: the fc is user-facing state checked for unread parameters,
      // and DoShooting writes the *resolved* shooting_parameter back into the fc (which
      // the final build then reads). A guess here would never be read.
      species->scf_parameters_[species->scf_tuning_index_] = g[0];
    }
  }

  std::string key = species->name();
  result.push_back({std::move(key), std::move(species)});
  return result;
}

std::vector<Named> ScalarFieldSpecies::CreateAxion(const SpeciesBuildContext& ctx,
                                                   const InputKeys& keys,
                                                   double omega0_scf) {
  std::vector<Named> result;

  const auto f_opt     = ctx.pfc->get<double>(keys.f_axion);
  const auto theta_opt = ctx.pfc->get<double>(keys.Theta_initial);
  const double n       = ctx.pfc->get_or(keys.n_axion, 1.);
  const char* pot      = keys.potential.c_str();
  class_test_severe(!f_opt,
                    "%s = axion requires '%s' (decay constant, reduced-Planck units)",
                    pot,
                    keys.f_axion.c_str());
  class_test(*f_opt <= 0.,
             "%s = axion requires %s > 0 (decay constant, reduced-Planck units)",
             pot,
             keys.f_axion.c_str());
  class_test_severe(!theta_opt, "%s = axion requires '%s'", pot, keys.Theta_initial.c_str());
  class_test(*theta_opt <= 0. || *theta_opt >= _PI_,
             "%s = axion requires %s in (0, pi): phi_ini = Theta*f",
             pot,
             keys.Theta_initial.c_str());
  class_test(n < 1., "%s must be >= 1 (V = m^2 f^2 (1-cos(phi/f))^n)", keys.n_axion.c_str());
  const double f     = *f_opt;
  const double theta = *theta_opt;

  // A (1-cos)^n potential has no exponential attractor: frozen ICs are forced.
  class_test_severe(ctx.pfc->get_yes(keys.attractor_ic),
                    "%s = yes is incompatible with %s = axion "
                    "(frozen ICs phi = %s * %s, phi' = 0 are set automatically)",
                    keys.attractor_ic.c_str(),
                    pot,
                    keys.Theta_initial.c_str(),
                    keys.f_axion.c_str());
  class_test_severe(ctx.pfc->get<std::vector<double>>(keys.parameters).has_value(),
                    "%s is not used with %s = axion; give %s, %s, %s, %s instead",
                    keys.parameters.c_str(),
                    pot,
                    keys.m_axion.c_str(),
                    keys.f_axion.c_str(),
                    keys.n_axion.c_str(),
                    keys.Theta_initial.c_str());
  ctx.pfc->get<int>(keys.tuning_index);  // consume; the axion branch forces 0

  std::vector<double> params = {ctx.pfc->get_or(keys.m_axion, 0.), f, n, theta};

  // Resolved shooting value (written back by DoShooting) wins over any seed.
  bool shooting_param_present = false;
  if (auto sp = ctx.pfc->get<double>(keys.shooting_parameter)) {
    params[0]              = *sp;
    shooting_param_present = true;
  }

  auto species = std::make_unique<ScalarFieldSpecies>(*ctx.pba,
                                                      omega0_scf,
                                                      std::move(params),
                                                      /*scf_tuning_index=*/0,
                                                      /*attractor_ic_scf=*/false,
                                                      /*phi_ini_scf=*/theta * f,
                                                      /*phi_prime_ini_scf=*/0.,
                                                      AxionScalarFieldPotential(),
                                                      /*beta=*/0.,
                                                      keys.instance_name);
  if (!shooting_param_present) {
    species->shooting_target_ = {keys.Omega, keys.shooting_parameter, omega0_scf};
    species->needs_shooting_  = true;
    if (species->scf_parameters_[0] <= 0.) {
      std::vector<double> g;
      std::vector<double> d;
      species->ComputeShootingGuess(ctx, g, d);
      species->scf_parameters_[0] = g[0];
    }
  }
  std::string key = species->name();
  result.push_back({std::move(key), std::move(species)});
  return result;
}

std::vector<Named> ScalarFieldSpecies::CreateAllForComposite(const SpeciesBuildContext& ctx,
                                                             double beta) {
  // Builds the scalar-field child the Type3 composite owns: an injectable pure-1EXP
  // potential V = V0 * exp(-lambda*phi) with non-attractor FROZEN initial conditions
  // (phi_ini = 1e-4, phi'_ini = 0), shooting V0 to hit Omega_scf. Frozen ICs keep the
  // field negligible during radiation domination (so the Omega_r radiation-start check
  // is not tripped), and the (1-2beta) kinetic term in the background then has no early
  // effect. Distinct from the standalone CreateAll, which uses the exponential-
  // quintessence default potential + attractor ICs and must stay byte-identical.
  std::vector<Named> result;

  // Resolve Omega0_scf (closure override, or user Omega_scf > 0 ⇒ shooting V0).
  double omega_scf_val   = 0.;
  auto omega_scf_opt     = ctx.pfc->get<double>("Omega_scf");
  bool omega_scf_present = omega_scf_opt.has_value();
  if (omega_scf_present) {
    omega_scf_val = *omega_scf_opt;
  }

  double omega0_scf     = 0.;
  bool shooting_allowed = false;
  if (ctx.omega0_closure_override.has_value()) {
    omega0_scf = *ctx.omega0_closure_override;
  }
  else if (omega_scf_present && omega_scf_val > 0.) {
    omega0_scf       = omega_scf_val;
    shooting_allowed = true;
  }
  else {
    return result;
  }
  if (omega0_scf == 0.) {
    return result;
  }

  // scf_parameters for a Type-3 run = [V0, lambda]: tuning index 0 shoots V0, entry 1
  // is the 1EXP slope lambda. The V0 entry is a placeholder; the shooting guess
  // (3 H0^2 Omega) provides the real starting value, and DoShooting writes the resolved
  // V0 back into scf_shooting_parameter for the final build.
  std::vector<double> scf_parameters;
  try {
    if (auto scf_params_opt = ctx.pfc->get<std::vector<double>>("scf_parameters")) {
      scf_parameters = *scf_params_opt;
    }
  }
  catch (const std::exception& e) {
    class_stop_severe("scf_parameters parse error: %s", e.what());
  }
  class_test_severe(scf_parameters.size() < 2,
                    "Type-3 (scf_veta) coupling requires scf_parameters = 'V0, lambda' (two "
                    "entries): tuning index 0 shoots V0, entry 1 is the 1EXP slope lambda.");
  const int scf_tuning_index = 0;                         // shoot V0
  ctx.pfc->get_or("scf_tuning_index", scf_tuning_index);  // consume (composite forces 0)

  // The composite controls the scalar-field ICs (always non-attractor frozen), so it
  // ignores attractor_ic_scf's value, but must still consume it to satisfy the
  // unread-parameter check (the standalone CreateAll, which would read it, is skipped
  // when scf_veta is active).
  ctx.pfc->get<std::string>("attractor_ic_scf");

  if (auto scf_shooting_opt = ctx.pfc->get<double>("scf_shooting_parameter")) {
    scf_parameters[scf_tuning_index] = *scf_shooting_opt;  // resolved V0
  }

  // Injectable pure-1EXP potential, params = [V0, lambda].
  ScalarFieldPotential exp1exp;
  exp1exp.V = [](double phi, const std::vector<double>& p) { return p[0] * std::exp(-p[1] * phi); };
  exp1exp.dV = [](double phi, const std::vector<double>& p) {
    return -p[1] * p[0] * std::exp(-p[1] * phi);
  };
  exp1exp.ddV = [](double phi, const std::vector<double>& p) {
    return p[1] * p[1] * p[0] * std::exp(-p[1] * phi);
  };
  exp1exp.shooting_guess =
      [](double omega, double H0, const std::vector<double>&, int) -> std::pair<double, double> {
    // Frozen field: rho_phi,today ~ V0 and rho/H0^2 = Omega, so V0 ~ 3 H0^2 Omega.
    const double v0 = 3.0 * H0 * H0 * omega;
    return {v0, 3.0 * H0 * H0};
  };

  auto field = std::make_unique<ScalarFieldSpecies>(*ctx.pba,
                                                    omega0_scf,
                                                    std::move(scf_parameters),
                                                    scf_tuning_index,
                                                    /*attractor_ic_scf=*/false,
                                                    /*phi_ini_scf=*/1e-4,
                                                    /*phi_prime_ini_scf=*/0.,
                                                    std::move(exp1exp),
                                                    beta);

  if (shooting_allowed) {
    bool shooting_param_present = ctx.pfc->get<double>("scf_shooting_parameter").has_value();
    if (!shooting_param_present) {
      field->shooting_target_ = {"Omega_scf", "scf_shooting_parameter", omega_scf_val};
      field->needs_shooting_  = true;

      std::vector<double> g;
      std::vector<double> d;
      field->ComputeShootingGuess(ctx, g, d);                   // 1EXP guess: V0 ~ 3 H0^2 Omega
      field->scf_parameters_[field->scf_tuning_index_] = g[0];  // seed V0 guess
    }
  }

  result.push_back({"ScalarField", std::move(field)});
  return result;
}
