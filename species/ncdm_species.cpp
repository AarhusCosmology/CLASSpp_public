#include "ncdm_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include <cmath>

// ── Background ─────────────────────────────────────────────────────────────

void NCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_   = index_bg; index_bg += pba_->N_ncdm;
  index_bg_rho_      = index_bg; index_bg += pba_->N_ncdm;  // base class protected
  index_bg_p_        = index_bg; index_bg += pba_->N_ncdm;  // base class protected
  index_bg_pseudo_p_ = index_bg; index_bg += pba_->N_ncdm;

  if (pba_->has_ncdm_decay_dr == _TRUE_ && ncdm_) {
    index_bg_lnf_decay_dr1_  = index_bg; index_bg += ncdm_->q_total_size_dncdm_;
    index_bg_dlnfdlnq_decay_ = index_bg; index_bg += ncdm_->q_total_size_dncdm_;
    index_bg_dlnfdlnq_sep_   = index_bg; index_bg += ncdm_->q_total_size_dncdm_;
  }
}

void NCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  if (pba_->has_ncdm_decay_dr == _TRUE_ && ncdm_) {
    index_bi_lnf_decay_dr1_           = index_bi; index_bi += ncdm_->q_total_size_dncdm_;
    index_bi_dlnfdlnq_separate_decay_ = index_bi; index_bi += ncdm_->q_total_size_dncdm_;
  }
}

void NCDMSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  if (!ncdm_) return;
  double z = 1. / a_rel - 1.;

  // For decaying NCDM: update distribution function weights from ODE vars
  if (pba_->has_ncdm_decay_dr == _TRUE_) {
    for (const auto& [ncdm_id, dncdm_props] : ncdm_->decay_dr_map_) {
      for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id]; ++i) {
        double f_from_lnf = std::exp(
            pvecback_B[index_bi_lnf_decay_dr1_ + dncdm_props.q_offset + i]);
        ncdm_->SetBackgroundWeight(ncdm_id, i, f_from_lnf * dncdm_props.dq[i]);
        pvecback[index_bg_lnf_decay_dr1_ + dncdm_props.q_offset + i] =
            pvecback_B[index_bi_lnf_decay_dr1_ + dncdm_props.q_offset + i];
        pvecback[index_bg_dlnfdlnq_sep_ + dncdm_props.q_offset + i] =
            pvecback_B[index_bi_dlnfdlnq_separate_decay_ + dncdm_props.q_offset + i];
      }
    }
  }

  for (int n = 0; n < pba_->N_ncdm; ++n) {
    double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
    ncdm_->background_ncdm_momenta(n, z, &number_ncdm, &rho_ncdm, &p_ncdm,
                                    nullptr, &pseudo_p_ncdm);
    pvecback[index_bg_number_   + n] = number_ncdm;
    pvecback[index_bg_rho_      + n] = rho_ncdm;
    pvecback[index_bg_p_        + n] = p_ncdm;
    pvecback[index_bg_pseudo_p_ + n] = pseudo_p_ncdm;
  }
}

void NCDMSpecies::BackgroundDerivs(double /*tau*/, const double* /*y*/, double* dy,
                                    const double* pvecback) {
  if (!pba_->has_ncdm_decay_dr || !ncdm_) return;
  const double a = pvecback[bgm_->index_bg_a_];
  for (const auto& [ncdm_id, dncdm_props] : ncdm_->decay_dr_map_) {
    const double M_ncdm = ncdm_->M_ncdm_[ncdm_id];
    const double Gamma  = dncdm_props.Gamma;
    for (int i = 0; i < ncdm_->q_size_ncdm_[ncdm_id]; ++i) {
      const double q       = ncdm_->q_ncdm_[ncdm_id][i];
      const double epsilon = std::sqrt(q * q + a * a * M_ncdm * M_ncdm);
      dy[index_bi_lnf_decay_dr1_ + dncdm_props.q_offset + i] =
          -a * a * M_ncdm * Gamma / epsilon;
      dy[index_bi_dlnfdlnq_separate_decay_ + dncdm_props.q_offset + i] =
          a * a * M_ncdm * Gamma * q * q / std::pow(epsilon, 3);
    }
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void NCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                               const perturb_workspace* ppw,
                                               int /*gauge*/) {
  if (!pba_->has_ncdm || !ncdm_) return;

  pv->index_pt_psi0_ncdm1 = index_pt;
  index_pt_psi0_ = index_pt;
  pv->N_ncdm = pba_->N_ncdm;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);

  for (int n = 0; n < pba_->N_ncdm; ++n) {
    if (fa_on) {
      pv->l_max_ncdm[n] = 2;
      pv->q_size_ncdm[n] = 1;
    } else {
      pv->q_size_ncdm[n] = ncdm_->q_size_ncdm_[n];
    }
    pv->index_ncdm_[n].clear();
    for (int iq = 0; iq < pv->q_size_ncdm[n]; ++iq)
      pv->index_ncdm_[n].push_back(index_pt + iq * (pv->l_max_ncdm[n] + 1));
    index_pt += (pv->l_max_ncdm[n] + 1) * pv->q_size_ncdm[n];
  }
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

  auto* bgm = ppaw.perturbations_module->GetBackgroundModule().get();
  const double* pvecback = ppw->pvecback;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on);

  if (fa_on) {
    // Fluid approximation (ncdmfa_on)
    for (int n = 0; n < pv->N_ncdm; ++n) {
      const double rho_ncdm  = pvecback[bgm->index_bg_rho_ncdm1_ + n];
      const double p_ncdm    = pvecback[bgm->index_bg_p_ncdm1_   + n];
      const double pseudo_p  = pvecback[bgm->index_bg_pseudo_p_ncdm1_ + n];
      const double w_ncdm    = p_ncdm / rho_ncdm;
      const double ca2_ncdm  = w_ncdm / 3. / (1. + w_ncdm)
                                * (5. - pseudo_p / p_ncdm);
      const double ceff2     = ca2_ncdm;
      const double cvis2     = 3. * w_ncdm * ca2_ncdm;
      const double pseudo_p_over_p = pseudo_p / p_ncdm;
      const int idx = pv->index_ncdm_.at(n)[0];

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
    }
  } else {
    // Exact Boltzmann hierarchy per momentum bin
    for (int n = 0; n < pv->N_ncdm; ++n) {
      const double M_ncdm = ncdm_->M_ncdm_[n];
      for (int iq = 0; iq < pv->q_size_ncdm[n]; ++iq) {
        const double q = ncdm_->q_ncdm_[n][iq];
        double dlnf0_dlnq;
        if (ncdm_->ncdm_types_[n] == NCDMType::standard) {
          dlnf0_dlnq = ncdm_->dlnf0_dlnq_ncdm_[n][iq];
        } else {
          // decay_dr: use time-dependent dlnfdlnq from pvecback
          dlnf0_dlnq = pvecback[bgm->index_bg_dlnfdlnq_ncdm_decay_dr1_
                                 + ncdm_->decay_dr_map_.at(n).q_offset + iq];
        }

        const double epsilon = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
        const double qk_div_epsilon = k * q / epsilon;
        const int idx = pv->index_ncdm_.at(n)[iq];
        const int lmax = pv->l_max_ncdm[n];

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
}

// ── Integrated observables ──────────────────────────────────────────────────

double NCDMSpecies::Delta(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const {
  // Detailed delta requires summing over q-bins with rho_ncdm weights, which requires
  // the full perturb_workspace context. Return 0 here; the perturbations module computes
  // the NCDM source terms directly using pv->index_ncdm_.
  return 0.;
}

double NCDMSpecies::Theta(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const {
  return 0.;
}

double NCDMSpecies::DeltaP(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const {
  return 0.;
}

double NCDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const {
  return 0.;
}
