#include "type3_species.h"

#include <stdexcept>
#include <string>

#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

Type3Species::Type3Species(const background& pba,
                           double omega0_cdm,
                           std::unique_ptr<ScalarFieldSpecies> scf)
    : CompositeSpecies("CDM_SCF_Momentum", EnergyType::Other) {
  auto cdm = std::make_unique<CDMSpecies>(pba, omega0_cdm, /*coupled=*/true);
  cdm_     = cdm.get();
  scf_     = scf.get();
  children_.push_back(std::move(cdm));
  children_.push_back(std::move(scf));
}

void Type3Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  cdm_->WriteBackgroundColumnTitles(w);
  scf_->WriteBackgroundColumnTitles(w);
}

void Type3Species::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  cdm_->WriteBackgroundData(pvecback, w);
  scf_->WriteBackgroundData(pvecback, w);
}

void Type3Species::WriteOutputColumns(PerturbColumnWriter& writer,
                                      const PerturbationsModule& mod,
                                      file_format fmt,
                                      TransferColumnSection section) const {
  cdm_->WriteOutputColumns(writer, mod, fmt, section);
  scf_->WriteOutputColumns(writer, mod, fmt, section);
}

// PrintVariables is composite-owned (NOT forwarded): it takes no layout, so each
// child would re-derive its layout from pv->species_layouts[collection_index_],
// but the children are not separately registered in the collection and have no
// valid collection_index_. Mirror idm_dr_idr_species.cpp: read this composite's
// own nested layout and replicate each child's print logic against it.
void Type3Species::PrintVariables(PerturbColumnWriter& w,
                                  double /*tau*/,
                                  const double* y,
                                  const PerturbationsModule& mod,
                                  const perturb_workspace* ppw) const {
  double delta_cdm = 0., theta_cdm = 0.;
  double delta_scf = 0., theta_scf = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const double* pvecmetric = ppw->pvecmetric.data();
    const double k           = ppw->scalar_ctx.k;
    const double a2          = ppw->scalar_ctx.a2;
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const perturbs* ppt      = mod.GetPerturbs();

    const auto& my_lay = static_cast<const CompositeSpecies::PerturbLayout&>(
        *pv->species_layouts[collection_index_]);
    const auto& cdm_lay = cdm_layout(my_lay);
    const auto& scf_lay = scf_layout(my_lay);

    // ── CDM (replicates CDMSpecies::PrintVariables; coupled child => theta dynamical) ──
    delta_cdm = y[cdm_lay.idx_delta];
    theta_cdm = (cdm_lay.idx_theta >= 0) ? y[cdm_lay.idx_theta] : 0.;

    // ── Scalar field (replicates ScalarFieldSpecies::PrintVariables) ──────────
    const double phi_prime_bg = pvecback[scf_->index_bg_phi_prime_scf()];
    const double dV_bg        = pvecback[scf_->index_bg_dV_scf()];
    const double rho_scf      = scf_->Rho(pvecback);
    const double p_scf        = scf_->P(pvecback);

    double delta_phi_prime = y[scf_lay.idx_phi_prime];
    if (ppt->gauge == possible_gauges::newtonian) {
      const double phi  = y[pv->index_pt_phi];
      delta_phi_prime  += phi_prime_bg * (pvecmetric[ppw->index_mt_psi] + 3. * phi);
    }

    double delta_rho_scf = 0.;
    if (ppt->gauge == possible_gauges::synchronous) {
      delta_rho_scf = 1. / 3. *
                      ((1. - 2. * scf_->beta()) / a2 * phi_prime_bg * delta_phi_prime +
                       dV_bg * y[scf_lay.idx_phi]);
    }
    else {
      delta_rho_scf = 1. / 3. *
                      ((1. - 2. * scf_->beta()) / a2 * phi_prime_bg * delta_phi_prime +
                       dV_bg * y[scf_lay.idx_phi] -
                       1. / a2 * phi_prime_bg * phi_prime_bg * pvecmetric[ppw->index_mt_psi]);
    }

    const double rho_plus_p_theta_scf = 1. / 3. * ppw->scalar_ctx.k2 / a2 * phi_prime_bg *
                                        y[scf_lay.idx_phi];

    delta_scf = delta_rho_scf / rho_scf;
    theta_scf = rho_plus_p_theta_scf / (rho_scf + p_scf);

    // ── Synchronous => Newtonian gauge correction (each child's own formula) ──
    if (ppt->gauge == possible_gauges::synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      delta_cdm          -= 3. * H * a * alpha;
      theta_cdm          += k * k * alpha;
      const double w_scf  = p_scf / rho_scf;
      delta_scf          += alpha * (-3. * H * (1. + w_scf));
      theta_scf          += k * k * alpha;
    }
  }

  w.Add("delta_cdm", delta_cdm, true);
  w.Add("theta_cdm", theta_cdm, true);
  w.Add("delta_scf", delta_scf, true);
  w.Add("theta_scf", theta_scf, true);
}

// ── Coupling free functions (arXiv:1604.04222, synchronous gauge) ────────────
// Zbar = -phi_prime_bg/a;  D = 3*rho_cdm*a2 - 2*beta*Zbar^2. Every term carries an
// explicit beta factor => both vanish identically at beta = 0.
double Type3CouplingDeltaPhiPrime(double beta, double phi_prime_bg, double theta_cdm) {
  return 2. * (beta / (1. - 2. * beta)) * phi_prime_bg * theta_cdm;
}

double Type3CouplingDeltaThetaCdm(double beta,
                                  double k2,
                                  double a_prime_over_a,
                                  double rho_cdm,
                                  double a2,
                                  double Zbar,
                                  double dV,
                                  double phi,
                                  double phi_prime,
                                  double phi_prime_bg,
                                  double theta_cdm) {
  const double one_m_2b = 1. - 2. * beta;
  const double D        = 3. * rho_cdm * a2 - 2. * beta * Zbar * Zbar;
  const double term1    = -k2 * 2. * beta * (one_m_2b * Zbar * Zbar * phi_prime + dV * phi) /
                          (one_m_2b * D);
  /* phi_prime_bg is factored into the parenthesis analytically (k2*phi/phi_prime_bg *
     phi_prime_bg = k2*phi) so the coupling stays finite when the background field starts
     at rest (phi_prime_bg = 0, e.g. thawing/non-attractor ICs). Dividing by phi_prime_bg
     here gave 0/0 = NaN, which poisoned the whole theta_cdm Jacobian row and tripped the
     "singular matrix" abort in the evolver. */
  const double term2 = 4. * beta * dV * (k2 * phi - 2. * beta * theta_cdm * phi_prime_bg) /
                       (one_m_2b * D);
  const double term3 = -(6. * beta * a_prime_over_a * Zbar * Zbar + 4. * beta * dV * phi_prime_bg) *
                       theta_cdm / D;
  return term1 + term2 + term3;
}

void Type3Species::AddCouplingDerivs(double /*tau*/,
                                     const double* y,
                                     double* dy,
                                     const perturb_parameters_and_workspace& ppaw) const {
  const double beta = scf_->beta();
  if (beta == 0.)
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv.get();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const auto& my                  = static_cast<const CompositeSpecies::PerturbLayout&>(
      *pv->species_layouts[collection_index_]);
  const auto& cdm_lay = cdm_layout(my);
  const auto& scf_lay = scf_layout(my);

  const double* pvecback      = ppw->pvecback.data();
  const double a              = ctx.a;
  const double a2             = ctx.a2;
  const double a_prime_over_a = ctx.a_prime_over_a;
  const double k2             = ctx.k2;
  const double phi_prime_bg   = pvecback[scf_->index_bg_phi_prime_scf()];
  const double dV             = pvecback[scf_->index_bg_dV_scf()];
  const double rho_cdm        = cdm_->Rho(pvecback);
  const double Zbar           = -phi_prime_bg / a;

  const double phi       = y[scf_lay.idx_phi];
  const double phi_prime = y[scf_lay.idx_phi_prime];
  const double theta_cdm = y[cdm_lay.idx_theta];

  dy[scf_lay.idx_phi_prime] += Type3CouplingDeltaPhiPrime(beta, phi_prime_bg, theta_cdm);
  dy[cdm_lay.idx_theta]     += Type3CouplingDeltaThetaCdm(beta,
                                                          k2,
                                                          a_prime_over_a,
                                                          rho_cdm,
                                                          a2,
                                                          Zbar,
                                                          dV,
                                                          phi,
                                                          phi_prime,
                                                          phi_prime_bg,
                                                          theta_cdm);
}

// The coupling's -(2 beta/3) Zbar^2 theta_cdm contribution to (rho+p)theta.
// Shared by StressEnergy (ICs) and DelegateTally (hot path). 0 at beta = 0.
double Type3Species::CrossTermRhoPlusPTheta(const BaseSpecies::PerturbLayout& base,
                                            const double* y,
                                            const double* pvecback,
                                            const perturb_workspace* ppw) const {
  const double beta = scf_->beta();
  if (beta == 0.)
    return 0.;
  const auto& my            = static_cast<const CompositeSpecies::PerturbLayout&>(base);
  const double a            = ppw->scalar_ctx.a;
  const double phi_prime_bg = pvecback[scf_->index_bg_phi_prime_scf()];
  const double Zbar         = -phi_prime_bg / a;
  const double theta_cdm    = y[cdm_layout(my).idx_theta];
  return -(2. * beta / 3.) * Zbar * Zbar * theta_cdm;
}

BaseSpecies::StressEnergyContribution Type3Species::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* pv,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  StressEnergyContribution se  = CompositeSpecies::StressEnergy(base, pv, y, pvecback, ppw);
  se.rho_plus_p_theta         += CrossTermRhoPlusPTheta(base, y, pvecback, ppw);
  return se;
}

void Type3Species::DelegateTally(const BaseSpecies::PerturbLayout& base,
                                 const perturb_vector* pv,
                                 const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw,
                                 StressEnergyContribution& total,
                                 StressEnergyContribution& total_cold,
                                 StressEnergyContribution& total_warm) const {
  // Per-child sum into total/cold/warm (the matter-tally classification).
  CompositeSpecies::DelegateTally(base, pv, y, pvecback, ppw, total, total_cold, total_warm);
  // The coupling cross-term is part of the total Einstein (rho+p)theta source but
  // belongs to neither matter bucket (it is the scalar field's theta_cdm-induced
  // momentum, and the scalar field is not cold/warm matter), so add to `total` only.
  total.rho_plus_p_theta += CrossTermRhoPlusPTheta(base, y, pvecback, ppw);
}

std::vector<Named> Type3Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // Coupling active iff scf_veta is present and nonzero.
  auto veta_opt = ctx.pfc->get<double>("scf_veta");
  if (!veta_opt.has_value() || *veta_opt == 0.)
    return result;
  const double beta = *veta_opt;

  if (beta >= 0.5)
    throw std::invalid_argument(
        "scf_veta (beta) must be < 1/2: beta >= 1/2 has a ghost / strong-coupling pathology.");

  // Synchronous gauge only. Mirror the parser's Newtonian resolution
  // (input_module.cpp): it resolves Newtonian on any of "newtonian"/"Newtonian"/"new",
  // so a bare "new"/"newt" must be rejected here too — a find("newton") check lets it
  // slip past and the coupling would silently run in Newtonian gauge.
  if (auto gauge = ctx.pfc->get<std::string>("gauge")) {
    if (gauge->find("newtonian") != std::string::npos ||
        gauge->find("Newtonian") != std::string::npos || gauge->find("new") != std::string::npos)
      throw std::invalid_argument(
          "Type-3 (scf_veta) coupling is implemented in synchronous gauge only.");
  }

  // CDM Omega from the resolved coupled-species budget (same source CDMSpecies uses).
  const double omega0_cdm = (ctx.omega_budget && ctx.omega_budget->cdm.has_value())
                                ? *ctx.omega_budget->cdm
                                : 0.;

  // Build the scalar-field child via the shared scalar-field input path, then
  // hand its single instance to the composite. ScalarFieldSpecies::CreateAllForComposite
  // returns at most one "ScalarField" entry (with shooting state set); it must be
  // present here.
  std::vector<Named> scf_built = ScalarFieldSpecies::CreateAllForComposite(ctx, beta);
  if (scf_built.empty())
    throw std::invalid_argument(
        "scf_veta is set but no scalar field was configured (set Omega_scf and scf_parameters).");
  auto scf_ptr = std::unique_ptr<ScalarFieldSpecies>(
      static_cast<ScalarFieldSpecies*>(scf_built.front().species.release()));

  result.push_back({"CDM_SCF_Momentum",
                    std::make_unique<Type3Species>(*ctx.pba, omega0_cdm, std::move(scf_ptr))});
  return result;
}
