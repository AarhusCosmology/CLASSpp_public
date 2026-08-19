#include "dncdm_dr_species.h"

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
      result.push_back(
          {e.key, std::make_unique<DNCDM_DR_Species>(std::move(e.species), ctx.pba, ctx.bgm)});
    }
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
  const double* pvecback  = ctx.ppw->pvecback.data();
  const double rho_dr     = dr_sp_->Rho(pvecback);
  const double decay_corr = (rho_dr > 0.) ? ctx.a * dncdm_->Gamma() * dncdm_->Rho(pvecback) / rho_dr
                                          : 0.;
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
    }
    else {
      if (l == 0)
        return rprime_dr * y[dncdm_lay.index_per_q[0]];
      else if (l == 1)
        return rprime_dr * y[dncdm_lay.index_per_q[0] + 1] / k;
      else
        return 0.;
    }
  };

  const int base = dr_lay.idx_F0;
  for (int l = 0; l <= lmax; ++l) {
    double collision_term = 0.;
    if ((l <= ppr->l_max_dr_col) && (l < 800)) {
      collision_term = compute_collision_integral(l);
    }
    dy[base + l] += collision_term;
  }
}

// ── Output ────────────────────────────────────────────────────────────────────

// ── Closure ───────────────────────────────────────────────────────────────────

double DNCDM_DR_Species::GetOmega0() const {
  // Combined mode: Omega_dncdmdr is the user input; initial-shooting iterations:
  // the shooter writes <flavor>.Omega_dncdmdr, read into Omega_dncdmdr_pending_.
  if (dncdm_->Omega_dncdmdr_pending().has_value())
    return *dncdm_->Omega_dncdmdr_pending();
  // Pre-shooting discovery build with no pinned value yet: fall back to the
  // child-sum (matter child + DR child(0)). Replaced by the pinned value once
  // DoShooting writes the unknown.
  return CompositeSpecies::GetOmega0();
}

// ── Shooter hooks ─────────────────────────────────────────────────────────────

std::vector<ShootingTarget> DNCDM_DR_Species::GetShootingTargets() const {
  if (dncdm_->InitialAbundanceMode()) {
    // Initial mode: the given initial abundance drives deg; shoot Omega_dncdmdr itself (the
    // closure reserve) as a fixed point, driven to the integrated combined density.
    // target_value is only a Newton seed (the residual is a fixed point, not a difference to it).
    const double seed = dncdm_->Omega_dncdmdr_pending().value_or(
        dncdm_->Omega_ini_pending().value_or(0.1));
    return {{name() + ".Omega_dncdmdr_fixedpoint", name() + ".Omega_dncdmdr", seed}};
  }
  if (dncdm_->Omega_dncdmdr_pending().has_value()) {
    // Combined mode: shoot deg to hit the user-specified combined today-density.
    // target_name = the input key; unknown_param = the fc key DoShooting varies (.deg).
    return {{name() + ".Omega_dncdmdr", name() + ".deg", *dncdm_->Omega_dncdmdr_pending()}};
  }
  return {};
}

void DNCDM_DR_Species::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                            std::vector<double>& guess,
                                            std::vector<double>& dxdy) const {
  if (dncdm_->InitialAbundanceMode()) {
    // Seed Omega_dncdmdr with the given initial abundance. The fixed-point residual is ~linear
    // in this unknown (it sets only Lambda, weakly coupled to the integrated combined via H),
    // so Newton converges in ~1-2 steps; d(unknown)/d(residual) ~= 1.
    // For a Neff_ini-only spec (no Omega_ini) convert to Omega scale via the photon energy
    // density, matching ApplyDncdmInitialClosure(), rather than falling back to the 0.1 default.
    double seed = 0.1;
    if (dncdm_->Omega_ini_pending().has_value()) {
      seed = *dncdm_->Omega_ini_pending();
    }
    else if (dncdm_->Neff_ini_pending().has_value() && ctx.pba) {
      seed = *dncdm_->Neff_ini_pending() * 7. / 8. * std::pow(4. / 11., 4. / 3.) *
             ctx.pba->Omega0_g;
    }
    guess.push_back(seed);
    dxdy.push_back(1.0);
    return;
  }
  if (!dncdm_->Omega_dncdmdr_pending().has_value())
    return;
  auto [g, d] = dncdm_->DegGuessFromOmegaToday(ctx, *dncdm_->Omega_dncdmdr_pending());
  guess.push_back(g);
  dxdy.push_back(d);
}

double DNCDM_DR_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                 const ShootingTarget& target) const {
  const double* bg      = ctx.bg_today;
  const double H0       = ctx.pba->H0;
  const double combined = (dncdm_->Rho(bg) + dr_sp_->Rho(bg)) / (H0 * H0);
  if (target.target_name == name() + ".Omega_dncdmdr_fixedpoint") {
    // Initial mode: drive the reserved Omega_dncdmdr (= GetOmega0()) to the integrated combined.
    return -combined + GetOmega0();
  }
  // Combined mode: drive the integrated combined to the requested target (Omega units).
  return combined - target.target_value;
}
