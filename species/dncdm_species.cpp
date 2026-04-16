#include "dncdm_species.h"

#include <cmath>
#include <vector>

#include "arrays.h"
#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations_module.h"

// ── Background ─────────────────────────────────────────────────────────────

void DNCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_   = index_bg++;
  index_bg_rho_      = index_bg++;  // base class protected
  index_bg_p_        = index_bg++;  // base class protected
  index_bg_pseudo_p_ = index_bg++;

  if (ncdm_ && ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
    index_bg_lnf_decay_dr1_   = index_bg;
    index_bg                 += ncdm_->q_size_ncdm_[ncdm_id_];
    index_bg_dlnfdlnq_decay_  = index_bg;
    index_bg                 += ncdm_->q_size_ncdm_[ncdm_id_];
    index_bg_dlnfdlnq_sep_    = index_bg;
    index_bg                 += ncdm_->q_size_ncdm_[ncdm_id_];
  }
}

void DNCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  if (ncdm_ && ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
    index_bi_lnf_decay_dr1_            = index_bi;
    index_bi                          += ncdm_->q_size_ncdm_[ncdm_id_];
    index_bi_dlnfdlnq_separate_decay_  = index_bi;
    index_bi                          += ncdm_->q_size_ncdm_[ncdm_id_];
  }
}

void DNCDMSpecies::SetBackgroundInitialConditions(double /*a_rel*/, double* pvecback_integration) {
  if (!ncdm_ || ncdm_->ncdm_types_[ncdm_id_] != NCDMType::decay_dr)
    return;
  const auto& dncdm_props = ncdm_->decay_dr_map_.at(ncdm_id_);
  for (int index_q = 0; index_q < ncdm_->q_size_ncdm_[ncdm_id_]; index_q++) {
    double q  = ncdm_->q_ncdm_[ncdm_id_][index_q];
    double f0 = ncdm_->w_ncdm_[ncdm_id_][index_q] / dncdm_props.dq[index_q];

    pvecback_integration[index_bi_lnf_decay_dr1_ + index_q] = std::log(f0);
    pvecback_integration[index_bi_dlnfdlnq_separate_decay_ + index_q] =
        -q * std::exp(q) / (std::exp(q) + 1);  // Fermi-Dirac
  }
}

void DNCDMSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  if (!ncdm_)
    return;
  double z = 1. / a_rel - 1.;

  // For decaying NCDM: update distribution function weights from ODE vars
  if (ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
    const int q_size        = ncdm_->q_size_ncdm_[ncdm_id_];
    const auto& dncdm_props = ncdm_->decay_dr_map_.at(ncdm_id_);

    std::vector<double> lnf_dlnf_array(2 * q_size);
    std::vector<double> ddlnf_array(q_size);
    std::vector<double> lnq(q_size);

    if (pba_->background_method == bgevo_evolver) {
      bool has_problem = false;
      for (int i = 0; i < q_size; i++) {
        if (pvecback_B[index_bi_lnf_decay_dr1_ + i] <= -460) {
          has_problem = true;
          break;
        }
        lnf_dlnf_array[i] = pvecback_B[index_bi_lnf_decay_dr1_ + i];
        lnq[i]            = std::log(ncdm_->q_ncdm_[ncdm_id_][i]);
      }
      if (has_problem) {
        // Use Fermi-Dirac approximation when f becomes static
        for (int i = 0; i < q_size; i++) {
          double q                               = ncdm_->q_ncdm_[ncdm_id_][i];
          double lnf_fd                          = -std::log(1. + std::exp(q));
          double dlnfdlnq_fd                     = -q * std::exp(q) / (1. + std::exp(q));
          pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_fd;
          pvecback[index_bg_dlnfdlnq_decay_ + i] = dlnfdlnq_fd;
          ncdm_->SetBackgroundWeight(ncdm_id_, i, std::exp(lnf_fd) * dncdm_props.dq[i]);
        }
      }
      else {
        // Find df/dq by first splining and then calculating the derivative
        class_call(array_spline_table_lines(lnq.data(),
                                            q_size,
                                            lnf_dlnf_array.data(),
                                            1,
                                            ddlnf_array.data(),
                                            _SPLINE_EST_DERIV_,
                                            bgm_->error_message_),
                   bgm_->error_message_,
                   bgm_->error_message_);

        class_call(array_derive_spline(lnq.data(),
                                       q_size,
                                       lnf_dlnf_array.data(),
                                       ddlnf_array.data(),
                                       1,
                                       0,
                                       q_size,
                                       bgm_->error_message_),
                   bgm_->error_message_,
                   bgm_->error_message_);

        for (int i = 0; i < q_size; i++) {
          pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf_array[i];
          pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf_array[q_size + i];
          pvecback[index_bg_dlnfdlnq_sep_ + i] = pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
          ncdm_->SetBackgroundWeight(ncdm_id_, i, std::exp(lnf_dlnf_array[i]) * dncdm_props.dq[i]);
        }
      }
    }
    else {
      for (int i = 0; i < q_size; i++) {
        lnf_dlnf_array[i] = pvecback_B[index_bi_lnf_decay_dr1_ + i];
        lnq[i]            = std::log(ncdm_->q_ncdm_[ncdm_id_][i]);
      }

      class_call(array_spline_table_lines(lnq.data(),
                                          q_size,
                                          lnf_dlnf_array.data(),
                                          1,
                                          ddlnf_array.data(),
                                          _SPLINE_EST_DERIV_,
                                          bgm_->error_message_),
                 bgm_->error_message_,
                 bgm_->error_message_);

      class_call(array_derive_spline(lnq.data(),
                                     q_size,
                                     lnf_dlnf_array.data(),
                                     ddlnf_array.data(),
                                     1,
                                     0,
                                     q_size,
                                     bgm_->error_message_),
                 bgm_->error_message_,
                 bgm_->error_message_);

      for (int i = 0; i < q_size; i++) {
        pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf_array[i];
        pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf_array[q_size + i];
        pvecback[index_bg_dlnfdlnq_sep_ + i]   = pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
        ncdm_->SetBackgroundWeight(ncdm_id_, i, std::exp(lnf_dlnf_array[i]) * dncdm_props.dq[i]);
      }
    }
  }

  double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
  ncdm_->background_ncdm_momenta(ncdm_id_,
                                 z,
                                 &number_ncdm,
                                 &rho_ncdm,
                                 &p_ncdm,
                                 nullptr,
                                 &pseudo_p_ncdm);
  pvecback[index_bg_number_]   = number_ncdm;
  pvecback[index_bg_rho_]      = rho_ncdm;
  pvecback[index_bg_p_]        = p_ncdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_ncdm;
}

void DNCDMSpecies::BackgroundDerivs(double /*tau*/,
                                    const double* /*y*/,
                                    double* dy,
                                    const double* pvecback) {
  if (!ncdm_ || ncdm_->ncdm_types_[ncdm_id_] != NCDMType::decay_dr)
    return;
  const double a          = pvecback[bgm_->index_bg_a_];
  const auto& dncdm_props = ncdm_->decay_dr_map_.at(ncdm_id_);
  const double M_ncdm     = ncdm_->M_ncdm_[ncdm_id_];
  const double Gamma      = dncdm_props.Gamma;
  for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id_]; ++i) {
    const double q                            = ncdm_->q_ncdm_[ncdm_id_][i];
    const double epsilon                      = std::sqrt(q * q + a * a * M_ncdm * M_ncdm);
    dy[index_bi_lnf_decay_dr1_ + i]           = -a * a * M_ncdm * Gamma / epsilon;
    dy[index_bi_dlnfdlnq_separate_decay_ + i] = a * a * M_ncdm * Gamma * q * q /
                                                std::pow(epsilon, 3);
  }
}

// ── Background column output ───────────────────────────────────────────────

void DNCDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  if (ncdm_ && ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
    for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id_]; i++) {
      snprintf(tmp, 40, "lnf_dncdm[%d][%d]", ncdm_id_, i);
      w.Add(tmp, 0.);
      snprintf(tmp, 40, "dlnfdlnq_dncdm[%d][%d]", ncdm_id_, i);
      w.Add(tmp, 0.);
      snprintf(tmp, 40, "dlnfdlnq_separate_dncdm[%d][%d]", ncdm_id_, i);
      w.Add(tmp, 0.);
    }
  }
}

void DNCDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, pvecback[bg_number_index()]);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, Rho(pvecback));
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, P(pvecback));
  if (ncdm_ && ncdm_->ncdm_types_[ncdm_id_] == NCDMType::decay_dr) {
    const int bg_lnf_idx          = bg_lnf_index();
    const int bg_dlnfdlnq_idx     = bg_dlnfdlnq_index();
    const int bg_dlnfdlnq_sep_idx = bg_dlnfdlnq_sep_index();
    for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id_]; i++) {
      snprintf(tmp, 40, "lnf_dncdm[%d][%d]", ncdm_id_, i);
      w.Add(tmp, pvecback[bg_lnf_idx + i]);
      snprintf(tmp, 40, "dlnfdlnq_dncdm[%d][%d]", ncdm_id_, i);
      w.Add(tmp, pvecback[bg_dlnfdlnq_idx + i]);
      snprintf(tmp, 40, "dlnfdlnq_separate_dncdm[%d][%d]", ncdm_id_, i);
      w.Add(tmp, pvecback[bg_dlnfdlnq_sep_idx + i]);
    }
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DNCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                               const precision* ppr,
                                               int& index_pt,
                                               const perturb_workspace* ppw,
                                               int /*gauge*/) {
  if (!ncdm_ || ncdm_->ncdm_types_[ncdm_id_] != NCDMType::decay_dr)
    return;

  index_pt_psi0_ = index_pt;

  pv->l_max_ncdm[ncdm_id_]  = ppr->l_max_ncdm;
  pv->q_size_ncdm[ncdm_id_] = ncdm_->q_size_ncdm_[ncdm_id_];

  pv->index_ncdm_[ncdm_id_].clear();
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq)
    pv->index_ncdm_[ncdm_id_].push_back(index_pt + iq * (pv->l_max_ncdm[ncdm_id_] + 1));
  index_pt += (pv->l_max_ncdm[ncdm_id_] + 1) * pv->q_size_ncdm[ncdm_id_];
}

void DNCDMSpecies::PerturbDerivs(double tau,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  if (!ncdm_ || ncdm_->ncdm_types_[ncdm_id_] != NCDMType::decay_dr)
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback;

  // Exact Boltzmann hierarchy per momentum bin (fa_on not supported for decaying)
  const double M_ncdm = ncdm_->M_ncdm_[ncdm_id_];
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    const double q    = ncdm_->q_ncdm_[ncdm_id_][iq];
    double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + iq];

    const double epsilon        = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
    const double qk_div_epsilon = k * q / epsilon;
    const int idx               = pv->index_ncdm_.at(ncdm_id_)[iq];
    const int lmax              = pv->l_max_ncdm[ncdm_id_];

    // l=0 (density)
    dy[idx] = -qk_div_epsilon * y[idx + 1] + metric_continuity * dlnf0_dlnq / 3.;
    // l=1 (velocity)
    dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) -
                  epsilon * metric_euler / (3. * q * k) * dlnf0_dlnq;
    // l=2 (shear)
    dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                  s_l[2] * metric_shear * 2. / 15. * dlnf0_dlnq;
    // l=3..lmax-1
    for (int l = 3; l < lmax; ++l)
      dy[idx + l] = qk_div_epsilon / (2. * l + 1.) *
                    (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
    // l=lmax (truncation)
    dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] - (1. + lmax) * k * cotKgen * y[idx + lmax];
  }
}

void DNCDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  if (!ncdm_ || ncdm_->ncdm_types_[ncdm_id_] != NCDMType::decay_dr ||
      pv->index_ncdm_.at(ncdm_id_).empty()) {
    return;
  }

  const double* pvecback = ctx.ppw->pvecback;
  for (int index_q = 0; index_q < pv->q_size_ncdm[ncdm_id_]; ++index_q) {
    const int idx           = pv->index_ncdm_[ncdm_id_][index_q];
    const double q          = ncdm_->q_ncdm_[ncdm_id_][index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * ncdm_->M_ncdm_[ncdm_id_] *
                                                    ncdm_->M_ncdm_[ncdm_id_]);
    const double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + index_q];
    const int lmax          = pv->l_max_ncdm[ncdm_id_];

    y[idx + 0] = -0.25 * ctx.delta_ur * dlnf0_dlnq;
    if (lmax >= 1)
      y[idx + 1] = -epsilon / 3. / q / ctx.k * ctx.theta_ur * dlnf0_dlnq;
    if (lmax >= 2)
      y[idx + 2] = -0.5 * ctx.shear_ur * dlnf0_dlnq;
    if (lmax >= 3)
      y[idx + 3] = -0.25 * ctx.l3_ur * dlnf0_dlnq;
  }
}

// ── Integrated observables ──────────────────────────────────────────────────

double DNCDMSpecies::Delta(const perturb_vector* pv,
                           const double* y,
                           const double* pvecback,
                           const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const double a        = ppw->scalar_ctx.a;
  double rho_delta_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0             = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) *
                            ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    const double q        = ncdm_->q_ncdm_[ncdm_id_][iq];
    const double epsilon  = std::sqrt(q * q + std::pow(ncdm_->M_ncdm_[ncdm_id_] * a, 2));
    rho_delta_ncdm       += q * q * epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / a, 4);
  return rho_delta_ncdm * factor / pvecback[index_bg_rho_];
}

double DNCDMSpecies::Theta(const perturb_vector* pv,
                           const double* y,
                           const double* pvecback,
                           const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const double a               = ppw->scalar_ctx.a;
  double rho_plus_p_theta_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0              = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) *
                             ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    const double q         = ncdm_->q_ncdm_[ncdm_id_][iq];
    rho_plus_p_theta_ncdm += q * q * q * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq] + 1];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / a, 4);
  const double k      = ppw->scalar_ctx.k;
  return rho_plus_p_theta_ncdm * k * factor / (pvecback[index_bg_rho_] + pvecback[index_bg_p_]);
}

double DNCDMSpecies::DeltaP(const perturb_vector* pv,
                            const double* y,
                            const double* pvecback,
                            const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const double a      = ppw->scalar_ctx.a;
  double delta_p_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0             = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) *
                            ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    const double q        = ncdm_->q_ncdm_[ncdm_id_][iq];
    const double epsilon  = std::sqrt(q * q + std::pow(ncdm_->M_ncdm_[ncdm_id_] * a, 2));
    delta_p_ncdm         += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / a, 4);
  return delta_p_ncdm * factor / 3.;
}

double DNCDMSpecies::RhoPlusPShear(const perturb_vector* pv,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const double a               = ppw->scalar_ctx.a;
  double rho_plus_p_shear_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0              = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) *
                             ncdm_->decay_dr_map_.at(ncdm_id_).dq[iq];
    const double q         = ncdm_->q_ncdm_[ncdm_id_][iq];
    const double epsilon   = std::sqrt(q * q + std::pow(ncdm_->M_ncdm_[ncdm_id_] * a, 2));
    rho_plus_p_shear_ncdm += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq] + 2];
  }
  const double factor = ncdm_->factor_ncdm_[ncdm_id_] * std::pow(pba_->a_today / a, 4);
  return 2.0 / 3.0 * factor * rho_plus_p_shear_ncdm;
}
