#include "dncdm_dr_species.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "background_module.h"
#include "errors.h"
#include "perturbations_module.h"
#include "species/dncdm_inv_species.h"
#include "species/dncdm_proxy_species.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"
#include "species/species_input.h"

std::vector<Named> DNCDM_DR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  auto dncdm_vec = DNCDMSpecies::CreateAll(ctx);
  std::vector<Named> result;
  result.reserve(dncdm_vec.size());
  for (auto& e : dncdm_vec) {
    // Per-instance dispatch on the DECAY-RADIATION REPRESENTATION, which is what
    // actually selects the implementation:
    //
    //   integrated (default) -> DNCDM_DR_Species: the daughters are one
    //                           momentum-INTEGRATED ultra-relativistic hierarchy.
    //   psd                  -> DNCDMInvSpecies: the daughters carry a resolved
    //                           per-q PSD, driven by the transition kernel.
    //
    // Inverse decays REQUIRE psd -- the inverse rate is f_l(q2) f_phi(q3) at the
    // transition's own momenta, which an integrated hierarchy simply does not
    // carry -- so inverse_decays=yes implies it. The reverse is not true, and that
    // asymmetry is the point of having a separate key: `inverse_decays = no` with
    // `dr_representation = psd` builds the composite with the kernel's inverse term
    // switched OFF, which is
    //   (a) the paper's decay-only ("dec") rung, and
    //   (b) the closest thing to an independent check this sector can have --
    //       the same pure-decay physics through two unrelated discretisations
    //       (resolved-PSD kernel vs integrated fluid), so agreement tests the
    //       composite's initial conditions, hierarchy, stress-energy assembly and
    //       metric coupling against code that does none of it the same way.
    SpeciesInput in(ctx.pfc, e.key);
    const bool inv        = in.get_flag("inverse_decays", false);
    const auto rep_opt    = in.get<std::string>("dr_representation");
    const std::string rep = rep_opt.value_or(inv ? "psd" : "integrated");
    // Structural: decidable from which keys are set, not from any numeric value a
    // sampler could vary -> severe.
    class_test_severe(rep != "psd" && rep != "integrated" && rep != "proxy",
                      "species '%s': dr_representation (='%s') must be 'psd', 'integrated' or "
                      "'proxy'",
                      e.key.c_str(),
                      rep.c_str());
    class_test_severe(inv && rep == "integrated",
                      "species '%s': inverse_decays = yes requires dr_representation = psd "
                      "(the inverse rate is f_l(q2)*f_phi(q3) at the transition momenta, which "
                      "a momentum-integrated daughter hierarchy does not carry). Drop "
                      "'dr_representation' to get it implicitly.",
                      e.key.c_str());
    if (rep == "proxy") {
      // proxy -> DNCDMProxySpecies: the daughters keep a COARSE background PSD
      // driven by the same transition kernel, but their perturbations collapse to
      // one integrated hierarchy each and the linearised collision operator is
      // replaced by a relaxation-time closure. Three or four orders of magnitude
      // cheaper than `psd`, which is what makes high-Gamma parameter estimation
      // possible at all; see dncdm_proxy_species.h for what is and is not kept.
      result.push_back(DNCDMProxySpecies::Create(std::move(e.species), ctx));
    }
    else if (rep == "psd") {
      result.push_back(DNCDMInvSpecies::Create(std::move(e.species), ctx));
    }
    else {
      // The parent's ln f is evolved on its perturbation grid and published on the background
      // one (DNCDMSpecies::ComputeBackground): with two grids the extra background bins never
      // decay and keep sourcing decay radiation (Omega_DR = 0.70 today at Gamma = 1e6 km/s/Mpc
      // from a 1.2e-3 parent). psd and proxy refuse this the same way.
      class_test_severe(e.species->q_size() != e.species->q_size_bg(),
                        "species '%s': a decaying NCDM needs a single parent momentum grid "
                        "(q_size=%d != q_size_bg=%d); set quadrature_strategy and momenta_bins "
                        "(or momenta_bins_bg == momenta_bins, or tol_ncdm == tol_ncdm_bg)",
                        e.key.c_str(),
                        e.species->q_size(),
                        e.species->q_size_bg());
      auto composite = std::make_unique<DNCDM_DR_Species>(std::move(e.species), ctx.pba, ctx.bgm);
      // The decay radiation's own l_max (#137); the parent read <name>.l_max itself.
      composite->dr().SetLMaxInput(in.get_l_max("l_max_dr"));
      result.push_back({e.key, std::move(composite)});
    }
  }
  return result;
}

DNCDM_DR_Species::DNCDM_DR_Species(std::unique_ptr<DNCDMSpecies> dncdm_arg,
                                   const background* pba,
                                   const BackgroundModule* bgm)
    : DNCDMSector(dncdm_arg.get()), pba_(pba), bgm_(bgm) {
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

void DNCDM_DR_Species::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  CompositeSpecies::SetBackgroundInitialConditions(ctx);
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
  auto& my = static_cast<CompositeSpecies::PerturbLayout&>(base);
  /* DR child first, DNCDM child second — contiguous per composite.
     NOTE: differs from the legacy split layout (all-DR-of-all-composites
     first, all-DNCDM-of-all-composites second) and from children_ order.
     See design spec section B.3 for the y-vector reorder rationale. */
  dr_sp_->RegisterPerturbationIndices(*my.child_layouts[kDr], pv, ppr, index_pt, ppw, gauge);
  dncdm_->RegisterPerturbationIndices(*my.child_layouts[kDncdm], pv, ppr, index_pt, ppw, gauge);
}

void DNCDM_DR_Species::PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                                           double tau,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) const {
  // DNCDMSpecies inherits NCDMBaseSpecies::PerturbTensorDerivs.
  // GetDlnf0Dlnq is overridden by DNCDMSpecies to use pvecback.
  dncdm_->PerturbTensorDerivs(dncdm_layout(base), tau, y, dy, ppaw);
}

void DNCDM_DR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                     double* y,
                                                     const PerturbIcContext& ctx) {
  dncdm_->PerturbSynchronousToNewtonian(dncdm_layout(base),
                                        y,
                                        ctx);  // NCDMBaseSpecies per-q transform
  const double* pvecback = ctx.ppw->pvecback.data();
  const double M_ncdm    = dncdm_->GetMass();
  const double Gamma     = dncdm_->Gamma();

  // Parent: the per-q transform above is the stable-NCDM one, in which f-bar(q) is
  // constant in conformal time. The decaying parent's is not. The convention here
  // is that of ApplyFluidLikeNewtonianShift, delta_N = delta_S + alpha rho'/rho
  // (radiation: -4 aH alpha), whose per-q form is Psi_N = Psi_S + alpha (ln f-bar)'
  // at fixed q. With (ln f-bar)' = -a^2 M Gamma / eps (DNCDMSpecies::BackgroundDerivs):
  //   Psi_0 += -alpha a^2 M Gamma / eps
  // (the per-q analogue of DCDM's -alpha a Gamma inside RhoPrimeOverRho). This
  // species runs only in Psi-space: the F-space representation belongs to the
  // inverse-decay composite, which is synchronous-gauge only.
  const auto& dncdm_lay = dncdm_layout(base);
  for (int iq = 0; iq < dncdm_lay.q_size; ++iq) {
    const double q                = dncdm_->GetQ()[iq];
    const double epsilon          = std::sqrt(q * q + ctx.a * ctx.a * M_ncdm * M_ncdm);
    y[dncdm_lay.index_per_q[iq]] -= ctx.alpha * ctx.a * ctx.a * M_ncdm * Gamma / epsilon;
  }

  // Daughters: delta_N = delta_S + alpha rho_dr'/rho_dr (same convention), with
  //   rho_dr' = -4 aH rho_dr + a Gamma M n_parent      (BackgroundDerivs above),
  // so PerturbNewtonianReseed's (-4 aH + decay_corr) alpha needs
  // decay_corr = +a Gamma M n / rho_dr. It is the parent's rest-mass
  // density M n that is injected, NOT its energy density rho: the two differ by
  // <eps>/(aM), which is huge while the parent is relativistic -- i.e. exactly when
  // the initial conditions are set (rho_dr itself grows from zero there, so the
  // correction is O(aH) with M n and orders of magnitude larger with rho).
  const double rho_dr = dr_sp_->Rho(pvecback);
  const double decay_corr =
      (rho_dr > 0.) ? ctx.a * Gamma * M_ncdm * pvecback[dncdm_->bg_number_index()] / rho_dr : 0.;
  dr_sp_->PerturbNewtonianReseed(dr_layout(base), y, ctx, decay_corr);
}

// ── Coupling derivs ───────────────────────────────────────────────────────────
// Registration order aside, PerturbDerivs, ICs, StressEnergy and the
// approximation-switch copy all use the generic CompositeSpecies child loops.

void DNCDM_DR_Species::AddCouplingDerivs(double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) const {
  const perturb_workspace* ppw    = ppaw.ppw;
  const precision* ppr            = ppaw.perturbations_module->GetPrecision();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const auto& dncdm_lay = dncdm_layout(*ppw->pv->species_layouts[collection_index_]);
  const auto& dr_lay    = dr_layout(*ppw->pv->species_layouts[collection_index_]);

  const double* pvecback = ppw->pvecback.data();
  const double a         = ctx.a;
  const double k         = ctx.k;

  const double M_ncdm = dncdm_->GetMass();
  const double Gamma  = dncdm_->Gamma();
  const int q_size    = dncdm_->q_size();

  // rprime_dr = a^5 * Gamma * M_ncdm * n_dncdm / H0^2
  const double rprime_dr = std::pow(a, 5) / (pba_->H0 * pba_->H0) * M_ncdm * Gamma *
                           pvecback[dncdm_->bg_number_index()];

  const int lmax = dr_lay.l_max;

  auto ComputeFl = [&](int index_q, std::vector<double>& output) {
    double q       = dncdm_->GetQ()[index_q];
    double epsilon = sqrt(q * q + a * a * M_ncdm * M_ncdm);
    double x       = q / epsilon;

    if (x < 0.9999) {
      int km = 42 + lmax;
      if (x > 0.9) {
        km *= static_cast<int>(-1.0 - 1.8 * log(1. / x - 1.0));
      }
      double Fp2 = 0.;
      double Fp1 = 1.;
      for (int l = km; l >= 0; --l) {
        double Fp = ((2 * l + 3) * Fp1 / x - l * Fp2) / (l + 3.);
        if ((Fp > 1e200) || (l == 0)) {
          Fp1 /= Fp;
          for (int ll = l + 1; ll <= lmax; ++ll) {
            output[ll * q_size + index_q] /= Fp;
          }
          Fp = 1.0;
        }
        if (l <= lmax) {
          output[l * q_size + index_q] = Fp;
        }
        Fp2 = Fp1;
        Fp1 = Fp;
      }
    }
    else {
      output[0 * q_size + index_q] = 1.;
      if (lmax > 0) {
        output[1 * q_size + index_q] = x;
      }
      if (lmax > 1) {
        output[2 * q_size + index_q] =
            (x * (5. * x * x - 3.) + 3. * pow(x * x - 1., 2.) * atanh(x)) / (2. * x * x * x);
      }
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

  // The parent is ALWAYS a full per-q hierarchy (DNCDMSpecies registers no fluid
  // variables and never collapses), so the source is always the momentum integral.
  // This used to branch on the global ncdmfa flag and, once it switched on
  // (k tau > ncdm_fluid_trigger_tau_over_tau_k), read Psi_0 / Psi_1 of the
  // LOWEST-q bin as if they were a fluid delta / theta and drop l >= 2 -- wrong in
  // both gauges, and in Newtonian gauge fatal: that bin's Psi_1 ~ (eps/q) psi is
  // huge for a non-relativistic parent and drove the DR velocity to a runaway
  // (P(k) off by 1e16, #136).
  auto compute_collision_integral = [&](int l) {
    double integral_num   = 0.;
    double integral_denom = 0.;
    bool must_rescale     = false;
    for (int index_q = 0; index_q < q_size; ++index_q) {
      double dq = dncdm_->dq()[index_q];
      double w0 = dq * exp(pvecback[dncdm_->bg_lnf_index() + index_q]);
      double q  = dncdm_->GetQ()[index_q];

      if (w0 == 0.) {
        must_rescale = true;
        break;
      }

      int psi_ind     = dncdm_lay.index_per_q[index_q] + l;
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
        int psi_ind     = dncdm_lay.index_per_q[index_q] + l;
        integral_num   += dq * q * q * exp(lnN + lnf) * y[psi_ind] * FL[l * q_size + index_q];
        integral_denom += dq * q * q * exp(lnN + lnf);
      }
    }
    return rprime_dr * integral_num / integral_denom;
  };

  // The collision at l reads the parent's l-th multipole: it stops at the smaller
  // hierarchy, each per species (#137), or earlier at an explicit l_max_dr_col (which
  // perturb_init has checked against both).
  const int l_col = std::min(
      {ppr->l_max_dr_col >= 0 ? ppr->l_max_dr_col : lmax, lmax, dncdm_lay.l_max});
  const int base = dr_lay.idx_F0;
  for (int l = 0; l <= lmax; ++l) {
    double collision_term = 0.;
    if ((l <= l_col) && (l < 800)) {
      collision_term = compute_collision_integral(l);
    }
    dy[base + l] += collision_term;
  }

  // Newtonian gauge: the injection rate a Gamma M n runs on the parent's proper time,
  // so the isotropic source carries a factor (1 + psi), exactly as DCDM_DR_Species
  // adds rprime_dr * metric_euler / k^2 to its DR monopole. metric_euler = k^2 psi in
  // Newtonian gauge and exactly 0 in synchronous gauge; the branch keeps the
  // synchronous-gauge arithmetic untouched.
  if (ppaw.perturbations_module->GetPerturbs()->gauge == possible_gauges::newtonian) {
    dy[base + 0] += rprime_dr * ctx.metric_euler / (k * k);
  }
}

// ── Output ────────────────────────────────────────────────────────────────────
