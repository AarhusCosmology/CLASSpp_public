#include "dncdm_dr_species.h"

#include <cmath>
#include <vector>

#include "background_module.h"
#include "perturbations_module.h"

DNCDM_DR_Species::DNCDM_DR_Species(int ncdm_id,
                                   std::shared_ptr<NonColdDarkMatter> ncdm,
                                   const background* pba,
                                   const BackgroundModule* bgm)
    : CompositeSpecies("DNCDM_DR_" + std::to_string(ncdm_id), BaseSpecies::EnergyType::Other),
      ncdm_id_(ncdm_id), pba_(pba), bgm_(bgm) {
  auto dncdm = std::make_unique<DNCDMSpecies>(ncdm_id, ncdm, pba, bgm);
  auto dr_sp = std::make_unique<DNCDM_DecayRadiationSpecies>(ncdm_id, pba, bgm);
  dncdm_     = dncdm.get();
  dr_sp_     = dr_sp.get();
  children_.push_back(std::move(dncdm));
  children_.push_back(std::move(dr_sp));
}

void DNCDM_DR_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

void DNCDM_DR_Species::SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {
  // Initialize children first (DNCDM)
  CompositeSpecies::SetBackgroundInitialConditions(a_rel, pvecback_integration);
}

void DNCDM_DR_Species::BackgroundDerivs(double tau,
                                        const double* y,
                                        double* dy,
                                        const double* pvecback) {
  // Children handle their own dilution terms (and DNCDM's distribution function decay)
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);

  // DNCDM->DR decay source
  const double a          = pvecback[bgm_->index_bg_a_];
  const auto& dncdm_props = dncdm_->ncdm().decay_dr_map_.at(ncdm_id_);
  const double M_ncdm     = dncdm_->ncdm().M_ncdm_[ncdm_id_];
  const double Gamma      = dncdm_props.Gamma;

  dy[dr_sp_->bi_rho_index()] += a * Gamma * M_ncdm * pvecback[dncdm_->bg_number_index()];
}

void DNCDM_DR_Species::AddCouplingDerivs(double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const precision* ppr            = ppaw.perturbations_module->GetPrecision();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double* pvecback = ppw->pvecback;
  const double a         = ctx.a;
  const double k         = ctx.k;

  const auto& dncdm_props = dncdm_->ncdm().decay_dr_map_.at(ncdm_id_);
  const double M_ncdm     = dncdm_->ncdm().M_ncdm_[ncdm_id_];
  const double Gamma      = dncdm_props.Gamma;
  const int q_size        = dncdm_->ncdm().q_size_ncdm_[ncdm_id_];

  // rprime_dr = a^5 * Gamma * M_ncdm * n_dncdm / H0^2
  const double rprime_dr = std::pow(a, 5) / (pba_->H0 * pba_->H0) * M_ncdm * Gamma *
                           pvecback[dncdm_->bg_number_index()];

  // Ported logic from legacy manual block in perturbations_module.cpp
  auto ComputeFl = [&](int index_q, int lmax, std::vector<double>& output) {
    double q       = dncdm_->ncdm().q_ncdm_[ncdm_id_][index_q];
    double epsilon = sqrt(q * q + a * a * M_ncdm * M_ncdm);
    double x       = q / epsilon;

    if (x < 0.9999) {
      int km = 42 + lmax;
      if (x > 0.9)
        km *= int(-1.0 - 1.8 * log(1. / x - 1.0));
      double Fp2 = 0.;
      double Fp1 = 1.;
      for (int l = km; l >= 0; --l) {
        double Fp = ((2 * l + 3) * Fp1 / x - l * Fp2) / (l + 3.);
        if ((Fp > 1e200) || (l == 0)) {
          Fp1 /= Fp;
          for (int ll = l + 1; ll <= lmax; ++ll)
            output[ll * q_size + index_q] /= Fp;
          Fp = 1.0;
        }
        if (l <= lmax)
          output[l * q_size + index_q] = Fp;
        Fp2 = Fp1;
        Fp1 = Fp;
      }
    }
    else {
      output[0 * q_size + index_q] = 1.;
      if (lmax > 0)
        output[1 * q_size + index_q] = x;
      if (lmax > 1)
        output[2 * q_size + index_q] =
            (x * (5. * x * x - 3.) + 3. * pow(x * x - 1., 2.) * atanh(x)) / (2. * x * x * x);
      for (int l = 3; l <= lmax; ++l) {
        double Fm2                   = output[(l - 2) * q_size + index_q];
        double Fm1                   = output[(l - 1) * q_size + index_q];
        output[l * q_size + index_q] = ((2. * l - 1.) * Fm1 / x - (l + 1.) * Fm2) / (l - 2.);
      }
    }
  };

  std::vector<double> FL(q_size * (pv->l_max_dr + 1));
  for (int index_q = 0; index_q < q_size; ++index_q) {
    ComputeFl(index_q, pv->l_max_dr, FL);
  }

  auto compute_collision_integral = [&](int l) {
    double integral_num   = 0.;
    double integral_denom = 0.;

    if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_off) {
      bool must_rescale = false;
      for (int index_q = 0; index_q < q_size; ++index_q) {
        double dq = dncdm_props.dq[index_q];
        double w0 = dq * exp(pvecback[dncdm_->bg_lnf_index() + index_q]);
        double q  = dncdm_->ncdm().q_ncdm_[ncdm_id_][index_q];

        if (w0 == 0.) {
          must_rescale = true;
          break;
        }

        int psi_ind     = pv->index_ncdm_.at(ncdm_id_)[index_q] + l;
        integral_num   += w0 * q * q * y[psi_ind] * FL[l * q_size + index_q];
        integral_denom += w0 * q * q;
      }
      if (must_rescale) {
        integral_num   = 0.;
        integral_denom = 0.;
        double lnN = dncdm_->ncdm().GetRescalingFactor(ncdm_id_, pvecback + dncdm_->bg_lnf_index());
        for (int index_q = 0; index_q < q_size; ++index_q) {
          double dq       = dncdm_props.dq[index_q];
          double lnf      = pvecback[dncdm_->bg_lnf_index() + index_q];
          double q        = dncdm_->ncdm().q_ncdm_[ncdm_id_][index_q];
          int psi_ind     = pv->index_ncdm_.at(ncdm_id_)[index_q] + l;
          integral_num   += dq * q * q * exp(lnN + lnf) * y[psi_ind] * FL[l * q_size + index_q];
          integral_denom += dq * q * q * exp(lnN + lnf);
        }
      }
      return rprime_dr * integral_num / integral_denom;
    }
    else {
      if (l == 0)
        return rprime_dr * y[pv->index_ncdm_.at(ncdm_id_)[0]];
      else if (l == 1)
        return rprime_dr * y[pv->index_ncdm_.at(ncdm_id_)[0] + 1] / k;
      else
        return 0.;
    }
  };

  const int base = dr_sp_->pt_F0_index();
  for (int l = 0; l <= pv->l_max_dr; ++l) {
    double collision_term = 0.;
    if ((l <= ppr->l_max_dr_col) && (l < 800)) {
      collision_term = compute_collision_integral(l);
    }
    dy[base + l]                   += collision_term;
    dy[pv->index_pt_F0_dr_sum + l] += collision_term;
  }
}