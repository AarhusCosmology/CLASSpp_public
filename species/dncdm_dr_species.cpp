#include "dncdm_dr_species.h"

#include <cmath>
#include <vector>

#include "background_module.h"
#include "perturbations_module.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"

std::vector<Named> DNCDM_DR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  auto dncdm_vec = DNCDMSpecies::CreateAll(ctx);
  std::vector<Named> result;
  result.reserve(dncdm_vec.size());
  for (auto& e : dncdm_vec) {
    result.push_back(
        {e.key, std::make_unique<DNCDM_DR_Species>(std::move(e.species), ctx.pba, ctx.bgm)});
  }
  return result;
}

DNCDM_DR_Species::DNCDM_DR_Species(std::unique_ptr<DNCDMSpecies> dncdm_arg,
                                   const background* pba,
                                   const BackgroundModule* bgm)
    : CompositeSpecies(dncdm_arg->name(), BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm) {
  auto dr_sp = std::make_unique<DarkRadiationSpecies>(dncdm_arg->name() + "_DR", pba, bgm);
  dncdm_     = dncdm_arg.get();
  dr_sp_     = dr_sp.get();
  children_.push_back(std::move(dncdm_arg));
  children_.push_back(std::move(dr_sp));
}

void DNCDM_DR_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

void DNCDM_DR_Species::SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {
  CompositeSpecies::SetBackgroundInitialConditions(a_rel, pvecback_integration);
}

void DNCDM_DR_Species::BackgroundDerivs(double tau,
                                        const double* y,
                                        double* dy,
                                        const double* pvecback) {
  // Children handle their own dilution terms (and DNCDM's distribution function decay)
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);

  // DNCDM->DR decay source
  const double a      = pvecback[bgm_->index_bg_a_];
  const double M_ncdm = dncdm_->GetMass();
  const double Gamma  = dncdm_->Gamma();

  dy[dr_sp_->bi_rho_index()] += a * Gamma * M_ncdm * pvecback[dncdm_->bg_number_index()];
}

// ── Perturbation layout-based overrides ───────────────────────────────────────

void DNCDM_DR_Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                   perturb_vector* pv,
                                                   const precision* ppr,
                                                   int& index_pt,
                                                   const perturb_workspace* ppw,
                                                   int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  /* DR child first, DNCDM child second — contiguous per composite.
     NOTE: differs from the legacy split layout (all-DR-of-all-composites
     first, all-DNCDM-of-all-composites second). See design spec section
     B.3 for the y-vector reorder rationale. */
  dr_sp_->RegisterPerturbationIndices(my.dr, pv, ppr, index_pt, ppw, gauge);
  dncdm_->RegisterPerturbationIndices(my.dncdm, pv, ppr, index_pt, ppw, gauge);
}

void DNCDM_DR_Species::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                     double tau,
                                     const double* y,
                                     double* dy,
                                     const perturb_parameters_and_workspace& ppaw) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dncdm_->PerturbDerivs(my.dncdm, tau, y, dy, ppaw);
  dr_sp_->PerturbDerivs(my.dr, tau, y, dy, ppaw);
  AddCouplingDerivs(my, y, dy, ppaw);
}

void DNCDM_DR_Species::PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                                           double tau,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  // DNCDMSpecies inherits NCDMBaseSpecies::PerturbTensorDerivs.
  // GetDlnf0Dlnq is overridden by DNCDMSpecies to use pvecback.
  dncdm_->PerturbTensorDerivs(my.dncdm, tau, y, dy, ppaw);
}

void DNCDM_DR_Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                              double* y,
                                              const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dncdm_->ApplyInitialConditions(my.dncdm, y, ctx);
  dr_sp_->ApplyInitialConditions(my.dr, y, ctx);
}

void DNCDM_DR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                     double* y,
                                                     const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dncdm_->PerturbSynchronousToNewtonian(my.dncdm, y, ctx);  // NCDMBaseSpecies per-q transform
  const double* pvecback  = ctx.ppw->pvecback;
  const double rho_dr     = dr_sp_->Rho(pvecback);
  const double decay_corr = (rho_dr > 0.) ? ctx.a * dncdm_->Gamma() * dncdm_->Rho(pvecback) / rho_dr
                                          : 0.;
  dr_sp_->PerturbNewtonianReseed(my.dr, y, ctx, decay_corr);
}

double DNCDM_DR_Species::Delta(const BaseSpecies::PerturbLayout& base,
                               const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rho_d   = dncdm_->Rho(pvecback);
  const double rho_dr  = dr_sp_->Rho(pvecback);
  const double rho_tot = rho_d + rho_dr;
  if (rho_tot <= 0.)
    return 0.;
  return (rho_d * dncdm_->Delta(my.dncdm, pv, y, pvecback, ppw) +
          rho_dr * dr_sp_->Delta(my.dr, pv, y, pvecback, ppw)) /
         rho_tot;
}

double DNCDM_DR_Species::Theta(const BaseSpecies::PerturbLayout& base,
                               const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rpp_d   = dncdm_->Rho(pvecback) + dncdm_->P(pvecback);
  const double rpp_dr  = dr_sp_->Rho(pvecback) + dr_sp_->P(pvecback);
  const double rpp_tot = rpp_d + rpp_dr;
  if (rpp_tot <= 0.)
    return 0.;
  return (rpp_d * dncdm_->Theta(my.dncdm, pv, y, pvecback, ppw) +
          rpp_dr * dr_sp_->Theta(my.dr, pv, y, pvecback, ppw)) /
         rpp_tot;
}

double DNCDM_DR_Species::DeltaP(const BaseSpecies::PerturbLayout& base,
                                const perturb_vector* pv,
                                const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  return dncdm_->DeltaP(my.dncdm, pv, y, pvecback, ppw) +
         dr_sp_->DeltaP(my.dr, pv, y, pvecback, ppw);
}

double DNCDM_DR_Species::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                       const perturb_vector* pv,
                                       const double* y,
                                       const double* pvecback,
                                       const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  return dncdm_->RhoPlusPShear(my.dncdm, pv, y, pvecback, ppw) +
         dr_sp_->RhoPlusPShear(my.dr, pv, y, pvecback, ppw);
}

// ── Coupling derivs ───────────────────────────────────────────────────────────

void DNCDM_DR_Species::AddCouplingDerivs(const PerturbLayout& my,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const precision* ppr            = ppaw.perturbations_module->GetPrecision();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double* pvecback = ppw->pvecback;
  const double a         = ctx.a;
  const double k         = ctx.k;

  const double M_ncdm = dncdm_->GetMass();
  const double Gamma  = dncdm_->Gamma();
  const int q_size    = dncdm_->q_size();

  // rprime_dr = a^5 * Gamma * M_ncdm * n_dncdm / H0^2
  const double rprime_dr = std::pow(a, 5) / (pba_->H0 * pba_->H0) * M_ncdm * Gamma *
                           pvecback[dncdm_->bg_number_index()];

  const int lmax = my.dr.l_max;

  auto ComputeFl = [&](int index_q, std::vector<double>& output) {
    double q       = dncdm_->GetQ()[index_q];
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

  std::vector<double> FL(q_size * (lmax + 1));
  for (int index_q = 0; index_q < q_size; ++index_q) {
    ComputeFl(index_q, FL);
  }

  auto compute_collision_integral = [&](int l) {
    double integral_num   = 0.;
    double integral_denom = 0.;

    if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_off) {
      bool must_rescale = false;
      for (int index_q = 0; index_q < q_size; ++index_q) {
        double dq = dncdm_->dq()[index_q];
        double w0 = dq * exp(pvecback[dncdm_->bg_lnf_index() + index_q]);
        double q  = dncdm_->GetQ()[index_q];

        if (w0 == 0.) {
          must_rescale = true;
          break;
        }

        int psi_ind     = my.dncdm.index_per_q[index_q] + l;
        integral_num   += w0 * q * q * y[psi_ind] * FL[l * q_size + index_q];
        integral_denom += w0 * q * q;
      }
      if (must_rescale) {
        integral_num   = 0.;
        integral_denom = 0.;
        double lnN     = dncdm_->GetRescalingFactor(pvecback + dncdm_->bg_lnf_index());
        for (int index_q = 0; index_q < q_size; ++index_q) {
          double dq       = dncdm_->dq()[index_q];
          double lnf      = pvecback[dncdm_->bg_lnf_index() + index_q];
          double q        = dncdm_->GetQ()[index_q];
          int psi_ind     = my.dncdm.index_per_q[index_q] + l;
          integral_num   += dq * q * q * exp(lnN + lnf) * y[psi_ind] * FL[l * q_size + index_q];
          integral_denom += dq * q * q * exp(lnN + lnf);
        }
      }
      return rprime_dr * integral_num / integral_denom;
    }
    else {
      if (l == 0)
        return rprime_dr * y[my.dncdm.index_per_q[0]];
      else if (l == 1)
        return rprime_dr * y[my.dncdm.index_per_q[0] + 1] / k;
      else
        return 0.;
    }
  };

  const int base = my.dr.idx_F0;
  for (int l = 0; l <= lmax; ++l) {
    double collision_term = 0.;
    if ((l <= ppr->l_max_dr_col) && (l < 800)) {
      collision_term = compute_collision_integral(l);
    }
    dy[base + l] += collision_term;
  }
}

// ── Matter tally (warm DNCDM only; DR is radiation → 0) ─────────────────────────

double DNCDM_DR_Species::MatterRhoDelta(const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return dncdm_->IsMatterSpecies()
             ? dncdm_->Rho(pvecback) * dncdm_->Delta(my.dncdm, pv, y, pvecback, ppw)
             : 0.;
}

double DNCDM_DR_Species::MatterRhoPlusPTheta(const perturb_vector* pv,
                                             const double* y,
                                             const double* pvecback,
                                             const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return dncdm_->IsMatterSpecies() ? (dncdm_->Rho(pvecback) + dncdm_->P(pvecback)) *
                                         dncdm_->Theta(my.dncdm, pv, y, pvecback, ppw)
                                   : 0.;
}

// ── Output ────────────────────────────────────────────────────────────────────

void DNCDM_DR_Species::FillSources(const double* y,
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

  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);

  // ── delta_dr (this channel's slot) ───────────────────────────────────────────
  // r_dr == 0 when the channel carries no DR (e.g. Gamma == 0); write 0 then,
  // matching DarkRadiationSpecies::Delta's rho_dr <= 0 convention (avoids 0/0).
  if (p_mod->has_source_delta_dr_ == _TRUE_) {
    const double r_dr = (a2_rel / pba_->H0) * (a2_rel / pba_->H0) * dr_sp_->Rho(pvecback);
    const double src  = (r_dr > 0.)
                            ? y[my.dr.idx_F0] / r_dr +
                                  4. * a_prime_over_a * ctx.theta_over_k2  // N-body gauge corr.
                            : 0.;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_dr_ + dr_sp_->source_slot(),
                          ctx.index_tau,
                          ctx.index_k,
                          src);
  }

  // ── theta_dr (this channel's slot) ───────────────────────────────────────────
  if (p_mod->has_source_theta_dr_ == _TRUE_) {
    const double r_dr = (a2_rel / pba_->H0) * (a2_rel / pba_->H0) * dr_sp_->Rho(pvecback);
    const double src  = (r_dr > 0.) ? 3. / 4. * ctx.k * y[my.dr.idx_F0 + 1] / r_dr +
                                          ctx.theta_shift  // N-body gauge correction
                                    : 0.;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_dr_ + dr_sp_->source_slot(),
                          ctx.index_tau,
                          ctx.index_k,
                          src);
  }
}

void DNCDM_DR_Species::WriteOutputColumns(PerturbColumnWriter& w,
                                          const PerturbationsModule& mod,
                                          enum file_format fmt,
                                          BaseSpecies::TransferColumnSection section) const {
  if (fmt != class_format)
    return;
  const perturbs* ppt = mod.GetPerturbs();
  const int slot      = dr_sp_->source_slot();
  if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
    w.Add("d_" + dr_sp_->name(), mod.index_tp_delta_dr_ + slot, mod.has_source_delta_dr_ == _TRUE_);
  if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
    w.Add("t_" + dr_sp_->name(), mod.index_tp_theta_dr_ + slot, mod.has_source_theta_dr_ == _TRUE_);
}

void DNCDM_DR_Species::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                     const BaseSpecies::PerturbLayout& new_base,
                                                     const double* old_y,
                                                     double* new_y,
                                                     const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  dncdm_->CopyPerturbationsAcrossSwitch(old_l.dncdm, new_l.dncdm, old_y, new_y, ctx);
  dr_sp_->CopyPerturbationsAcrossSwitch(old_l.dr, new_l.dr, old_y, new_y, ctx);
}

// ── Shooter hooks ─────────────────────────────────────────────────────────────

std::vector<ShootingTarget> DNCDM_DR_Species::GetShootingTargets() const {
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return {};
  // target_name: the input key the user specified (per-flavor dot key, e.g. "dncdm1.Omega_dncdmdr")
  // unknown_param: the fc key DoShooting will vary   (e.g. "dncdm1.deg")
  return {{name() + ".Omega_dncdmdr", name() + ".deg", *dncdm_->Omega_dncdmdr_pending()}};
}

void DNCDM_DR_Species::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                            std::vector<double>& guess,
                                            std::vector<double>& dxdy) const {
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return;
  auto [g, d] = dncdm_->DegGuessFromOmegaToday(ctx, *dncdm_->Omega_dncdmdr_pending());
  guess.push_back(g);
  dxdy.push_back(d);
}

double DNCDM_DR_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                 const ShootingTarget& target) const {
  // Today-density of this flavor's dncdm+dr minus the requested target (Omega units).
  const double* bg = ctx.bg_today;
  const double H0  = ctx.pba->H0;
  return (dncdm_->Rho(bg) + dr_sp_->Rho(bg)) / (H0 * H0) - target.target_value;
}
