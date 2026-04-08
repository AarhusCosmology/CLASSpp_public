#include "ncdm_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include <cmath>

// ── Background ─────────────────────────────────────────────────────────────

void NCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_   = index_bg++;
  index_bg_rho_      = index_bg++; // base class protected
  index_bg_p_        = index_bg++; // base class protected
  index_bg_pseudo_p_ = index_bg++;

  if (pba_->has_ncdm_decay_dr == _TRUE_ && ncdm_) {
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
      index_bg_lnf_decay_dr1_  = index_bg; index_bg += ncdm_->q_size_ncdm_[ncdm_id_];
      index_bg_dlnfdlnq_decay_ = index_bg; index_bg += ncdm_->q_size_ncdm_[ncdm_id_];
      index_bg_dlnfdlnq_sep_   = index_bg; index_bg += ncdm_->q_size_ncdm_[ncdm_id_];
    }
  }
}

void NCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  if (pba_->has_ncdm_decay_dr == _TRUE_ && ncdm_) {
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
      index_bi_lnf_decay_dr1_           = index_bi; index_bi += ncdm_->q_size_ncdm_[ncdm_id_];
      index_bi_dlnfdlnq_separate_decay_ = index_bi; index_bi += ncdm_->q_size_ncdm_[ncdm_id_];
    }
  }
}

void NCDMSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  if (!ncdm_) return;
  double z = 1. / a_rel - 1.;

  // For decaying NCDM: update distribution function weights from ODE vars
  if (pba_->has_ncdm_decay_dr == _TRUE_ && ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
    const auto& dncdm_props = ncdm_->decay_dr_map_.at(ncdm_id_);
    for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id_]; ++i) {
      double f_from_lnf = std::exp(
          pvecback_B[index_bi_lnf_decay_dr1_ + i]);
      ncdm_->SetBackgroundWeight(ncdm_id_, i, f_from_lnf * dncdm_props.dq[i]);
      pvecback[index_bg_lnf_decay_dr1_ + i] =
          pvecback_B[index_bi_lnf_decay_dr1_ + i];
      pvecback[index_bg_dlnfdlnq_sep_ + i] =
          pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
    }
  }

  double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
  ncdm_->background_ncdm_momenta(ncdm_id_, z, &number_ncdm, &rho_ncdm, &p_ncdm,
                                  nullptr, &pseudo_p_ncdm);
  pvecback[index_bg_number_]   = number_ncdm;
  pvecback[index_bg_rho_]      = rho_ncdm;
  pvecback[index_bg_p_]        = p_ncdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_ncdm;
}

void NCDMSpecies::BackgroundDerivs(double /*tau*/, const double* /*y*/, double* dy,
                                    const double* pvecback) {
  if (!pba_->has_ncdm_decay_dr || !ncdm_ || ncdm_->ncdm_types_[ncdm_id_] != NCDMType::decay_dr) return;
  const double a = pvecback[bgm_->index_bg_a_];
  const auto& dncdm_props = ncdm_->decay_dr_map_.at(ncdm_id_);
  const double M_ncdm = ncdm_->M_ncdm_[ncdm_id_];
  const double Gamma  = dncdm_props.Gamma;
  for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id_]; ++i) {
    const double q       = ncdm_->q_ncdm_[ncdm_id_][i];
    const double epsilon = std::sqrt(q * q + a * a * M_ncdm * M_ncdm);
    dy[index_bi_lnf_decay_dr1_ + i] =
        -a * a * M_ncdm * Gamma / epsilon;
    dy[index_bi_dlnfdlnq_separate_decay_ + i] =
        a * a * M_ncdm * Gamma * q * q / std::pow(epsilon, 3);
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void NCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr, int& index_pt,
                                               const perturb_workspace* ppw,
                                               int /*gauge*/) {
  if (!pba_->has_ncdm || !ncdm_) return;

  // Note: we still set pv->index_pt_psi0_ncdm1 to the start of the first NCDM species.
  // This is for backward compatibility with code that hasn't been refactored yet.
  if (ncdm_id_ == 0) {
    pv->index_pt_psi0_ncdm1 = index_pt;
  }
  index_pt_psi0_ = index_pt;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);

  if (fa_on) {
    pv->l_max_ncdm[ncdm_id_] = 2;
    pv->q_size_ncdm[ncdm_id_] = 1;
  } else {
    pv->l_max_ncdm[ncdm_id_] = ppr->l_max_ncdm;
    pv->q_size_ncdm[ncdm_id_] = ncdm_->q_size_ncdm_[ncdm_id_];
  }
  pv->index_ncdm_[ncdm_id_].clear();
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq)
    pv->index_ncdm_[ncdm_id_].push_back(index_pt + iq * (pv->l_max_ncdm[ncdm_id_] + 1));
  index_pt += (pv->l_max_ncdm[ncdm_id_] + 1) * pv->q_size_ncdm[ncdm_id_];
}

void NCDMSpecies::PerturbDerivs(double tau, const double* y, double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  if (!pba_->has_ncdm || !ncdm_) return;

  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector*    pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l   = ppw->s_l;
  const double k      = ctx.k;
  const double a2     = ctx.a2;
  const double a_prime_over_a = ctx.a_prime_over_a;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;
  const double metric_shear      = ctx.metric_shear;
  const double metric_ufa_class  = ctx.metric_ufa_class;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);

  if (fa_on) {
    // Fluid approximation (ncdmfa_on)
    const double rho_ncdm  = pvecback[index_bg_rho_];
    const double p_ncdm    = pvecback[index_bg_p_];
    const double pseudo_p  = pvecback[index_bg_pseudo_p_];
    const double w_ncdm    = p_ncdm / rho_ncdm;
    const double ca2_ncdm  = w_ncdm / 3. / (1. + w_ncdm)
                              * (5. - pseudo_p / p_ncdm);
    const double ceff2     = ca2_ncdm;
    const double cvis2     = 3. * w_ncdm * ca2_ncdm;
    const double pseudo_p_over_p = pseudo_p / p_ncdm;
    const int idx = pv->index_ncdm_.at(ncdm_id_)[0];

    dy[idx] = -(1. + w_ncdm) * (y[idx + 1] + metric_continuity)
              - 3. * a_prime_over_a * (ceff2 - w_ncdm) * y[idx];

    dy[idx + 1] = -a_prime_over_a * (1. - 3. * ca2_ncdm) * y[idx + 1]
                  + ceff2 / (1. + w_ncdm) * k * k * y[idx]
                  - k * k * y[idx + 2]
                  + metric_euler;

    // CLASS fluid approximation for shear
    dy[idx + 2] = -3. * (a_prime_over_a * (2./3. - ca2_ncdm - pseudo_p_over_p / 3.)
                          + 1. / tau) * y[idx + 2]
                  + 8. / 3. * cvis2 / (1. + w_ncdm) * s_l[2]
                    * (y[idx + 1] + metric_ufa_class);
  } else {
    // Exact Boltzmann hierarchy per momentum bin
    const double M_ncdm = ncdm_->M_ncdm_[ncdm_id_];
    for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
      const double q = ncdm_->q_ncdm_[ncdm_id_][iq];
      double dlnf0_dlnq;
      if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::standard) {
        dlnf0_dlnq = ncdm_->dlnf0_dlnq_ncdm_[ncdm_id_][iq];
      } else {
        // decay_dr: use time-dependent dlnfdlnq from pvecback
        dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_sep_ + iq];
      }

      const double epsilon = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
      const double qk_div_epsilon = k * q / epsilon;
      const int idx = pv->index_ncdm_.at(ncdm_id_)[iq];
      const int lmax = pv->l_max_ncdm[ncdm_id_];

      // l=0 (density)
      dy[idx] = -qk_div_epsilon * y[idx + 1]
                + metric_continuity * dlnf0_dlnq / 3.;
      // l=1 (velocity)
      dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2])
                    - epsilon * metric_euler / (3. * q * k) * dlnf0_dlnq;
      // l=2 (shear)
      dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1]
                                             - 3. * s_l[3] * y[idx + 3])
                    - s_l[2] * metric_shear * 2. / 15. * dlnf0_dlnq;
      // l=3..lmax-1
      for (int l = 3; l < lmax; ++l)
        dy[idx + l] = qk_div_epsilon / (2. * l + 1.)
                      * (l * s_l[l] * y[idx + l - 1]
                         - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
      // l=lmax (truncation)
      dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1]
                       - (1. + lmax) * k * cotKgen * y[idx + lmax];
    }
  }
}

// ── Integrated observables ──────────────────────────────────────────────────

double NCDMSpecies::Delta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty()) return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);
  if (fa_on) return y[pv->index_ncdm_.at(ncdm_id_)[0]];

  double rho_delta_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0;
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::standard) {
      w0 = ncdm_->w_ncdm_[ncdm_id_][iq];
    } else {
      w0 = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    }
    const double q = ncdm_->q_ncdm_[ncdm_id_][iq];
    const double epsilon = std::sqrt(q * q + std::pow(ncdm_->M_ncdm_[ncdm_id_] * pvecback[bgm_->index_bg_a_], 2));
    rho_delta_ncdm += q * q * epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / pvecback[bgm_->index_bg_a_], 4);
  return rho_delta_ncdm * factor / pvecback[index_bg_rho_];
}

double NCDMSpecies::Theta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty()) return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);
  if (fa_on) return y[pv->index_ncdm_.at(ncdm_id_)[0] + 1];

  double rho_plus_p_theta_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0;
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::standard) {
      w0 = ncdm_->w_ncdm_[ncdm_id_][iq];
    } else {
      w0 = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    }
    const double q = ncdm_->q_ncdm_[ncdm_id_][iq];
    rho_plus_p_theta_ncdm += q * q * q * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq] + 1];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / pvecback[bgm_->index_bg_a_], 4);
  const double k = ppw->scalar_ctx.k;
  return rho_plus_p_theta_ncdm * k * factor / (pvecback[index_bg_rho_] + pvecback[index_bg_p_]);
}

double NCDMSpecies::DeltaP(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty()) return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);
  if (fa_on) {
    // cg2_ncdm logic for fluid approx
    double rho_bg = pvecback[index_bg_rho_];
    double p_bg = pvecback[index_bg_p_];
    double pseudo_p_bg = pvecback[index_bg_pseudo_p_];
    double w_ncdm, cg2_ncdm;
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::standard) {
      w_ncdm = p_bg / rho_bg;
      cg2_ncdm = w_ncdm * (1.0 - 1.0 / (3.0 + 3.0 * w_ncdm) * (3.0 * w_ncdm - 2.0 + pseudo_p_bg / p_bg));
    } else {
      double pseudo_p_over_p;
      std::tie(w_ncdm, pseudo_p_over_p) = ncdm_->GetRescaledParameters(ncdm_id_, pvecback[bgm_->index_bg_a_], pvecback + index_bg_lnf_decay_dr1_);
      cg2_ncdm = w_ncdm * (1.0 - 1.0 / (3.0 + 3.0 * w_ncdm) * (3.0 * w_ncdm - 2.0 + pseudo_p_over_p));
    }
    return cg2_ncdm * rho_bg * y[pv->index_ncdm_.at(ncdm_id_)[0]];
  }

  double delta_p_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0;
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::standard) {
      w0 = ncdm_->w_ncdm_[ncdm_id_][iq];
    } else {
      w0 = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    }
    const double q = ncdm_->q_ncdm_[ncdm_id_][iq];
    const double epsilon = std::sqrt(q * q + std::pow(ncdm_->M_ncdm_[ncdm_id_] * pvecback[bgm_->index_bg_a_], 2));
    delta_p_ncdm += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / pvecback[bgm_->index_bg_a_], 4);
  return delta_p_ncdm * factor / 3.;
}

double NCDMSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty()) return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);
  if (fa_on) return (pvecback[index_bg_rho_] + pvecback[index_bg_p_]) * y[pv->index_ncdm_.at(ncdm_id_)[0] + 2];

  double rho_plus_p_shear_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0;
    if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::standard) {
      w0 = ncdm_->w_ncdm_[ncdm_id_][iq];
    } else {
      w0 = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    }
    const double q = ncdm_->q_ncdm_[ncdm_id_][iq];
    const double epsilon = std::sqrt(q * q + std::pow(ncdm_->M_ncdm_[ncdm_id_] * pvecback[bgm_->index_bg_a_], 2));
    rho_plus_p_shear_ncdm += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq] + 2];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / pvecback[bgm_->index_bg_a_], 4);
  return 2.0 / 3.0 * factor * rho_plus_p_shear_ncdm;
}
