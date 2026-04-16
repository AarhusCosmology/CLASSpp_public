#include "dcdm_dr_species.h"

#include <cmath>

#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations_module.h"

DCDM_DR_Species::DCDM_DR_Species(std::shared_ptr<DarkRadiation> dr,
                                 const background* pba,
                                 const BackgroundModule* bgm)
    : CompositeSpecies("DCDM_DR", BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm) {
  auto dcdm  = std::make_unique<DCDMSpecies>(*pba);
  auto dr_sp = std::make_unique<DarkRadiationSpecies>(dr, pba, bgm, dcdm.get());
  dcdm_      = dcdm.get();
  dr_sp_     = dr_sp.get();
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(dr_sp));
}

void DCDM_DR_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

void DCDM_DR_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm", 0.);
  w.Add("(.)rho_dr", 0.);
  for (int j = 0; j < pba_->N_decay_dr; ++j) {
    char tmp[40];
    snprintf(tmp, 40, "(.)rho_dr[%d]", j);
    w.Add(tmp, 0.);
  }
}

void DCDM_DR_Species::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm", dcdm().Rho(pvecback));
  w.Add("(.)rho_dr", dr().Rho(pvecback));
  for (int j = 0; j < pba_->N_decay_dr; ++j) {
    char tmp[40];
    snprintf(tmp, 40, "(.)rho_dr[%d]", j);
    w.Add(tmp, pvecback[dr().bg_rho_dr_species_index() + j]);
  }
}

void DCDM_DR_Species::SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {
  // Initialize children first (DCDM)
  CompositeSpecies::SetBackgroundInitialConditions(a_rel, pvecback_integration);

  // Then add the DR initial condition from DCDM decay
  if (pba_->has_dcdm == _TRUE_) {
    const double Omega_rad    = pba_->Omega0_g + pba_->Omega0_ur;
    const double rho_dcdm_ini = pvecback_integration[dcdm_->bi_rho_index()];
    double f                  = 1. / 3. * std::pow(a_rel, 6) * rho_dcdm_ini * pba_->Gamma_dcdm /
                                std::pow(pba_->H0, 3) / std::sqrt(Omega_rad);
    pvecback_integration[dr_sp_->bi_rho_dr_species_index()] = f * std::pow(pba_->H0, 2) /
                                                              std::pow(a_rel, 4);
  }
}

void DCDM_DR_Species::BackgroundDerivs(double tau,
                                       const double* y,
                                       double* dy,
                                       const double* pvecback) {
  // Children handle their own dilution terms
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);

  // DCDM->DR decay source for first DR channel
  const double a                         = pvecback[bgm_->index_bg_a_];
  dy[dr_sp_->bi_rho_dr_species_index()] += a * pba_->Gamma_dcdm * dcdm_->Rho(pvecback);
}

void DCDM_DR_Species::FillSources(const double* y,
                                  const double* /*dy*/,
                                  PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const perturb_vector* pv   = ppw->pv;
  const double* pvecback     = ppw->pvecback;

  const double a_prime_over_a = ctx.a_prime_over_a;
  const double a2_rel         = ctx.a2_rel;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_dcdm ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_dcdm_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_dcdm_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[pv->index_pt_delta_dcdm] +
                              (3. * a_prime_over_a + ctx.a_rel * pba_->Gamma_dcdm) *
                                  ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_dcdm ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_dcdm_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_dcdm_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[pv->index_pt_theta_dcdm] + ctx.theta_shift);  // N-body gauge correction
  }

  // ── delta_dr ───────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_dr_ == _TRUE_) {
    const double r_dr = (a2_rel / pba_->H0) * (a2_rel / pba_->H0) *
                        pvecback[bgm_->index_bg_rho_dr_];
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_dr_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[pv->index_pt_F0_dr_sum] / r_dr +
                              4. * a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_dr ───────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_dr_ == _TRUE_) {
    const double r_dr = (a2_rel / pba_->H0) * (a2_rel / pba_->H0) *
                        pvecback[bgm_->index_bg_rho_dr_];
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_dr_,
                          ctx.index_tau,
                          ctx.index_k,
                          3. / 4. * ctx.k * y[pv->index_pt_F0_dr_sum + 1] / r_dr +
                              ctx.theta_shift);  // N-body gauge correction
  }
}

void DCDM_DR_Species::AddCouplingDerivs(double /*tau*/,
                                        const double* y,
                                        double* dy,
                                        const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  if (pv->index_pt_delta_dcdm < 0 || pv->index_pt_F0_dr_species < 0)
    return;

  const double* pvecback = ppw->pvecback;
  const double a         = ctx.a;
  const double k         = ctx.k;

  const int base         = pv->index_pt_F0_dr_species;  // first DR channel, index_dr=0
  const double rprime_dr = pba_->Gamma_dcdm * dcdm_->Rho(pvecback) * std::pow(a, 5) /
                           (pba_->H0 * pba_->H0);

  const double delta_dcdm = y[pv->index_pt_delta_dcdm];
  const double theta_dcdm = y[pv->index_pt_theta_dcdm];

  // Add DCDM source to DR l=0 and l=1
  const double dl0  = rprime_dr * (delta_dcdm + ctx.metric_euler / (k * k));
  const double dl1  = rprime_dr / k * theta_dcdm;
  dy[base + 0]     += dl0;
  dy[base + 1]     += dl1;

  // Keep sum slots consistent
  dy[pv->index_pt_F0_dr_sum + 0] += dl0;
  dy[pv->index_pt_F0_dr_sum + 1] += dl1;
}
