#include "dncdm_proxy_species.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

#include "background_module.h"
#include "perturbations_module.h"
#include "species/species_input.h"

namespace {

/** Euler–Mascheroni, for the small-x branch of 𝓕. */
constexpr double kEulerGamma = 0.5772156649015329;

/** Where 𝓕 switches from its series to its tail: the point where the two branches
 *  cross, so the switch is continuous. The paper says "x > 2", where they differ by
 *  a factor 5.7. */
constexpr double kCurlyFSeam = 1.9497103;

/** Coefficient of 𝓕's tail. */
constexpr double kCurlyFTail = 0.409;

/** Floor on a daughter occupation, mirroring DrPsdSpecies::kFFloor: the PSDs are
 *  sourced from zero and the kernel's (1±f) coefficients must see a positive
 *  number, not a signed round-off residue. */
constexpr double kFFloor = 1e-100;

/** Margin on the closed-form kinematic bound when defaulting the daughter q_max
 *  (the same 5% the exact composite uses: Gauss-Legendre nodes approach the band
 *  edge without reaching it). */
constexpr double kQMaxMargin = 1.05;

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// The transport rate
// ─────────────────────────────────────────────────────────────────────────────

double DNCDMProxySpecies::CurlyF(double x) {
  if (x <= 0.) {
    return 0.;
  }
  if (x < kCurlyFSeam) {
    // The clamp only keeps log() finite; 𝓕 ~ -ln x here, which is physical.
    const double xs = (x < 1e-300) ? 1e-300 : x;
    const double x2 = xs * xs;
    return -0.5 + 2. * xs - x2 - x2 * xs / 9. + x2 * x2 / 96. - x2 * x2 * xs / 900. +
           x2 * x2 * x2 / 8640. + (kEulerGamma + std::log(xs)) * (0.5 * x2 - 1.);
  }
  if (x > 40.) {  // e^-40 = 4e-18: dead, and this TU is built -ffinite-math-only
    return 0.;
  }
  // The asymptotic tail. 0.409 is the minimax fit to eq. (14) over 2 < x < 17, which
  // reproduces the "<10%" the paper claims for this branch; the 2.4 it prints is a
  // typo (538% there, and it would make the paper's own two branches disagree by 5.7x
  // where they meet).
  return kCurlyFTail * std::exp(-x) / (x * std::sqrt(x));
}

double DNCDMProxySpecies::TransportRate(int l, const double* pvecback) const {
  if (l < 2) {
    return 0.;  // the l <= 1 collision integrals vanish by energy/momentum conservation
  }
  return TransportRate(l, pvecback[bgm_->index_bg_a_], pvecback);
}

double DNCDMProxySpecies::TransportRate(int l, double a, const double* pvecback) const {
  if (l < 2) {
    return 0.;
  }
  const double X       = a * parent_->GetMass();  // = a m_nuH / T_0
  const double rho_sec = pvecback[index_bg_rho_sec_];
  if (rho_sec <= 0. || X <= 0.) {
    return 0.;
  }
  // Both terms carry the parent's energy fraction: the transport rate is a rate at
  // which the SECTOR loses anisotropic stress, and no parent means no transitions.
  // It is also what keeps this function arithmetically safe -- below the floor the
  // parent's density has underflowed to ~1e-60 of the sector and every ratio built
  // on it is noise, so the term is cut rather than evaluated.
  const double frac_H = parent_->Rho(pvecback) / rho_sec;
  if (!(frac_H > kParentFloor)) {
    return 0.;
  }
  if (C3_ == 0. && C5_ == 0.) {
    return 0.;
  }

  // gamma = <eps>/X, the parent's number-weighted mean Lorentz factor. <eps> =
  // a*rho/n exactly (rho and n share one normalisation and differ only by the eps
  // weight), so no extra quadrature. This is the same variable the reference runs'
  // rate curves are plotted against, which is what makes C3/C5 calibratable straight
  // off those figures.
  const double number = pvecback[parent_->bg_number_index()];
  if (!(number > 0.)) {
    return 0.;
  }
  double gamma = (a * parent_->Rho(pvecback) / number) / X;
  if (gamma < 1.) {
    gamma = 1.;  // gamma >= 1 identically; this only clamps quadrature round-off
  }
  const double g2     = gamma * gamma;
  const double g3i    = 1. / (g2 * gamma);
  const double eps_ne = pvecback[index_bg_eps_ne_];

  double rate3 = 0., rate5 = 0.;
  if (rta_form_ == RtaForm::kCOPW) {
    // Chen, Oldengott, Pierobon & Wong, arXiv:2203.09075 eq. (13), in full:
    //
    //   Gamma_T,l = alpha_l a Gamma^0 (1/12) (rho_H/rho_sec) X^5 F(X)
    //
    // ONE term. The paper derives the gamma^-5 rate from first principles and argues
    // explicitly against the gamma^-3 rate of earlier work, which it calls "not an
    // ab initio result" but a heuristic random-walk argument that its appendix shows
    // to be incomplete. So this branch carries no gamma^-3 piece: adding one would
    // make it a third form rather than the paper's, which defeats the point of being
    // able to select it. dr_rta_C3, dr_rta_n3 and dr_rta_vshut are all ignored here.
    //
    // The (1/12) IS the paper's amplitude -- eq. (13) has no free normalisation -- so
    // dr_rta_C5 defaults to 1 under this form rather than to the fitted 2.5228 (see
    // Create). Leaving the fitted value in place would silently return 2.52x the
    // published rate to anyone selecting `copw` to compare against the paper.
    //
    // rho_sec is the paper's own rho_nuphi, which it defines as the TOTAL
    // (nu_H, nu_l, phi) energy density, not the phi density alone. The phase-space
    // factor Phi(m_nul/m_nuH) is absent because it is identically 1 for the massless
    // daughters this composite carries.
    //
    // This is the one form derived rather than calibrated, which is why it is worth
    // being able to run -- but it does NOT reproduce the measured gamma-dependence
    // inside the converged window: over gamma = 3 -> 6 its X^5 F(X) falls by 10.0x
    // where the runs fall by 16.9x, because X ~ 1 there and F's own decline eats two
    // of the five powers. The gamma^-5 asymptote it is named for needs X << 1, which
    // is exactly where the runs stop converging. So it is the reference, not default.
    const double X2 = X * X;
    rate5           = C5_ * (1. / 12.) * X2 * X2 * X * CurlyF(X);
  }
  else {
    // The gamma^-3 amplitude carries the background's own departure from detailed
    // balance, and that is the whole content of the gamma^-3 / gamma^-5 dispute:
    // Hannestad-Raffelt's kinematic rate survives only out of balance, and Barenboim
    // et al. show the gamma^-3 term cancels identically AT balance. eps_ne falls from
    // 0.72 to 1e-4 across Gamma = 1e5 -> 1e10 in this code, and that fall is what
    // turns the measured effective exponent from 3 into 5.
    //
    // The fitted two-power form, Gamma_T/Gamma = C3 eps_ne^n3 gamma^-3 + C5 gamma^-5,
    // good to a few per cent everywhere with C5 constant across four decades of Gamma
    // and both masses. This is the form the shipped defaults are calibrated for.
    rate3 = C3_ * std::pow(eps_ne > 1e-12 ? eps_ne : 1e-12, n3_) * g3i;
    rate5 = C5_ * g3i / g2;

    // Non-relativistic shut-off. A pure power law in gamma has none, and that is not
    // a detail: the measured Gamma_T/H TURNS OVER below gamma ~ 2.5, and for
    // m = 0.3 eV recombination happens at gamma ~ 2.1 -- inside the turnover. Left
    // out, the rate runs away exactly where the CMB is being formed (measured 2x too
    // high at recombination and 17x too high by a = 4e-3 at Gamma = 1e10).
    //
    // The shut-off variable is the parent's squared velocity, v^2 = 1 - 1/gamma^2. It
    // needs no new scale and it is the right physics: inverse decay is kinematically
    // open only because the parent is boosted, so the rate must vanish with the
    // boost. kCOPW needs no such term -- F(X) shuts the rate off on its own. The
    // exponent is fitted, not derived.
    const double v2  = 1. - 1. / g2;
    const double S   = std::pow(v2 > 1e-12 ? v2 : 1e-12, vshut_);
    rate3           *= S;
    rate5           *= S;
  }

  const double beta3 = l * (l + 1.) / 6.;
  return a * parent_->Gamma() * frac_H * (beta3 * rate3 + AlphaL(l) * rate5);
}

// ─────────────────────────────────────────────────────────────────────────────
// Construction
// ─────────────────────────────────────────────────────────────────────────────

DNCDMProxySpecies::DNCDMProxySpecies(std::unique_ptr<DNCDMSpecies> parent,
                                     std::unique_ptr<DarkRadiationSpecies> fermion,
                                     std::unique_ptr<DarkRadiationSpecies> boson,
                                     std::vector<double> q_d,
                                     std::vector<double> dq_d,
                                     DecayTransitionKernel::Config cfg,
                                     const background* pba,
                                     const BackgroundModule* bgm)
    : CompositeSpecies(parent->name(), BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm),
      q_d_(std::move(q_d)), dq_d_(std::move(dq_d)) {
  parent_  = parent.get();
  fermion_ = fermion.get();
  boson_   = boson.get();
  parent_->SetCollisionOwned(true);  // before any index registration (design §3 table)

  // ONE grid for BOTH daughters. They are emitted in pairs at q2 + q3 = eps1, so a
  // shared grid is what lets the two-bin deposit be exact in energy as well as in
  // number -- the same reason the exact composite gives them identical grids.
  const GridView parent_view{parent_->GetQ().data(), parent_->dq().data(), parent_->q_size()};
  const GridView d_view{q_d_.data(), dq_d_.data(), static_cast<int>(q_d_.size())};
  kernel_ = std::make_unique<DecayTransitionKernel>(parent_view,
                                                    d_view,
                                                    d_view,
                                                    Statistics::Fermion,
                                                    Statistics::Boson,
                                                    cfg);

  // Child order is the child_layouts contract: kParent, kFermion, kBoson (#358).
  children_.push_back(std::move(parent));
  children_.push_back(std::move(fermion));
  children_.push_back(std::move(boson));

  const int Nd = static_cast<int>(q_d_.size());
  fH_bare_.assign(parent_->q_size(), 0.);
  df_H_.assign(parent_->q_size(), 0.);
  df_l_.assign(Nd, 0.);
  df_phi_.assign(Nd, 0.);
  diag_H_.assign(parent_->q_size(), 0.);
  diag_l_.assign(Nd, 0.);
  diag_phi_.assign(Nd, 0.);
}

void DNCDMProxySpecies::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

Named DNCDMProxySpecies::Create(std::unique_ptr<DNCDMSpecies> parent,
                                const SpeciesBuildContext& ctx) {
  const std::string name = parent->name();
  SpeciesInput in(ctx.pfc, name);

  // Same structural guards as the exact composite, for the same reasons: the RTA
  // is written in the synchronous gauge, has no tensor counterpart, and the fluid
  // approximation would throw away the very multipoles it damps.
  if (auto gauge = ctx.pfc->get<std::string>("gauge")) {
    std::string g = *gauge;
    std::transform(g.begin(), g.end(), g.begin(), [](unsigned char c) { return std::tolower(c); });
    class_test_severe(g.rfind("new", 0) == 0,
                      "species '%s': dr_representation = proxy supports the synchronous gauge only",
                      name.c_str());
  }
  if (auto modes = ctx.pfc->get<std::string>("modes")) {
    const bool has_tensor = modes->find('t') != std::string::npos ||
                            modes->find('T') != std::string::npos;
    class_test_severe(has_tensor,
                      "species '%s': dr_representation = proxy does not support tensor modes",
                      name.c_str());
  }
  class_test_severe(in.get<std::string>("fluid_approximation").has_value(),
                    "species '%s': dr_representation = proxy needs the full Boltzmann hierarchy "
                    "for the parent; remove '%s.fluid_approximation'",
                    name.c_str(),
                    name.c_str());
  class_test_severe(parent->q_size() != parent->q_size_bg(),
                    "species '%s': dr_representation = proxy requires a single parent momentum "
                    "grid (q_size=%d != q_size_bg=%d)",
                    name.c_str(),
                    parent->q_size(),
                    parent->q_size_bg());

  // The same two kernel-boundary preconditions DNCDMInvSpecies::Create enforces, for
  // the same reason: this representation converts the parent's stored occupation to
  // bare with the SAME scalar KappaStoredToBare(), and that reading assumes one
  // Fermi-Dirac shape. It matters more here, not less -- `inverse_decays` and
  // `quantum_statistics` both default to TRUE under proxy, so the conversion is on
  // unless the user turns it off.
  class_test_severe(parent->UsesPsdFile(),
                    "species '%s': dr_representation = proxy requires an analytic (non-file) "
                    "parent PSD; remove '%s.use_psd_file' (the kernel boundary's bare-occupation "
                    "reading assumes a single Fermi-Dirac shape, #385)",
                    name.c_str(),
                    name.c_str());
  class_test(parent->GetKsi() != 0.,
             "species '%s': dr_representation = proxy requires the parent's chemical potential "
             "ksi = 0 (FD(q-ksi)+FD(q+ksi) has no single bare-occupation reading at the kernel "
             "boundary, #385)",
             name.c_str());

  // kLMaxScratch sizes the fixed per-multipole stack arrays in the RTA closure. Reject
  // at construction rather than let AddCouplingDerivs bail out: bailing there drops the
  // coupling silently, leaving free-streaming daughters against an interacting
  // background -- a changed model with no diagnostic. Structural (a precision setting,
  // not a sampled value), hence severe.
  class_test_severe(ctx.ppr->l_max_ncdm >= kLMaxScratch || ctx.ppr->l_max_dr >= kLMaxScratch,
                    "species '%s': dr_representation = proxy supports multipoles up to "
                    "l_max = %d; got l_max_ncdm=%d, l_max_dr=%d",
                    name.c_str(),
                    kLMaxScratch - 1,
                    ctx.ppr->l_max_ncdm,
                    ctx.ppr->l_max_dr);

  // ── the shared daughter grid ───────────────────────────────────────────────
  // Log-trapezoid with half-weight endpoints, the same construction DrPsdSpecies
  // uses, so a proxy cell and an exact cell at the same dr_N_q integrate the same
  // way and their difference is the representation and nothing else.
  //
  // The DEFAULT is far coarser than the exact scheme's, because here the grid only
  // has to integrate f-bar's moments -- it no longer sets the perturbation state
  // vector, which is what forced dr_N_q ~ Gamma^0.25 there. Measured ladder at
  // Gamma = 1e7, m = 0.3 (sector comoving energy at a = 1, against the hpc_cmb
  // reference cell's 1.2536 at its own dr_N_q = 101):
  //
  //     dr_N_q     32       48       96      192      384
  //     E_comov  1.3738   1.2800   1.2540   1.2513   1.2509
  //     wall[s]    1.6      2.1      3.8      8.1     18.6
  //
  // 96 reproduces the reference to 3e-4 and costs under four seconds, so that is
  // the default. The requirement grows slowly with Gamma (the emission band
  // narrows); re-run the ladder rather than assuming it if a cell matters.
  //
  // ⚠ dr_q_min IS THE BIGGER LEVER, and it is not the one you would guess. The
  // default here is 5e-2, not the exact scheme's 1e-2, because the proxy only needs
  // the daughters' MOMENTS: q^2 f is negligible over the bottom decade of a log
  // grid, but the emission band still sweeps through all of those fine bins on
  // every node, so they are pure cost. Measured, Gamma = 1e9 background, E_comov
  // against a converged reference:
  //
  //     dr_q_min   1e-3     1e-2     5e-2     1e-1
  //     core-s     62.2     28.9     13.5      8.2
  //     error     +0.41%   +0.17%   +0.09%   +0.04%
  //
  // i.e. raising it is cheaper AND more accurate in both directions tested, and
  // 1e-3 is twice the cost for 2.4x the error. It holds at high Gamma: at
  // Gamma = 1e11 the whole span 3e-2 -> 1e-1 moves E_comov by 0.002% for 2.5x the
  // cost. 5e-2 rather than 1e-1 only as margin for masses/Gammas not swept.
  //
  // Two more, measured on the same ladder, that belong to the PARENT and so are the
  // caller's to set: momenta_bins = 12 beat 16 on both cost and accuracy
  // (21.4 s/+0.13% vs 28.5 s/+0.17%), and quadrature_strategy = 3 is the most
  // accurate of the four available (2 is 1.8x cheaper at +0.44%).
  //
  // Together, Gamma = 1e11 background: dr_N_q 224 / dr_q_min 1e-1 / momenta_bins 12
  // is 28.5 core-s against 431 core-s for dr_N_q 448 / 1e-2 / 16, at the same
  // E_comov to 0.05%.
  const double q_star_max = DecayTransitionKernel::MaxDaughterMomentum(parent->GetQ().back(),
                                                                       parent->GetMass());
  const double q_min      = in.get_or("dr_q_min", 5e-2);
  const double q_max      = in.get_or("dr_q_max", kQMaxMargin * q_star_max);
  const int N_q           = in.get_or("dr_N_q", 96);
  class_test_severe(N_q < 8 || q_min <= 0. || q_max <= q_min,
                    "species '%s': bad daughter grid (dr_N_q=%d, dr_q_min=%g, dr_q_max=%g)",
                    name.c_str(),
                    N_q,
                    q_min,
                    q_max);

  // dr_q_sampling: `log` (default) or `linear`. Same two families DrPsdSpecies
  // offers, and for the same reason -- but the trade is DIFFERENT here. The exact
  // scheme needs the daughter PSD as a resolved function, so it wants log's dynamic
  // range; the proxy only needs its MOMENTS and the pointwise f at the transition
  // momenta, and q^2 f is negligible over the bottom decade of a log grid. Which
  // wins is an empirical question -- sweep it, do not assume.
  const std::string sampling = in.get<std::string>("dr_q_sampling").value_or("log");
  class_test_severe(sampling != "log" && sampling != "linear",
                    "species '%s': dr_q_sampling (='%s') must be 'log' or 'linear'",
                    name.c_str(),
                    sampling.c_str());

  std::vector<double> q_d(N_q), dq_d(N_q);
  if (sampling == "linear") {
    // Uniform trapezoid over [0, q_max] with q = 0 an implicit zero-weight point
    // (q^2 f vanishes there), so the FIRST node takes full weight and only the last
    // is halved. Byte-identical to quadrature.c's qm_trapz with qmin == 0, which is
    // the grid the published arXiv:2011.01502 background figures use.
    // dr_q_min is IGNORED in this mode, by that convention.
    const double h = q_max / N_q;
    for (int i = 0; i < N_q; ++i) {
      q_d[i]  = (i + 1) * h;
      dq_d[i] = (i == N_q - 1) ? 0.5 * h : h;
    }
  }
  else {
    const double h = std::log(q_max / q_min) / (N_q - 1);
    for (int i = 0; i < N_q; ++i) {
      q_d[i] = q_min * std::exp(i * h);
    }
    for (int i = 0; i < N_q; ++i) {
      const double w = (i == 0 || i == N_q - 1) ? 0.5 : 1.0;  // trapezoid endpoints
      dq_d[i]        = w * q_d[i] * h;
    }
  }

  DecayTransitionKernel::Config cfg;
  // Inverse decays are not optional here: the proxy exists to stand in for the
  // ν_H ↔ ν_l φ sector, and a decay-only run is served exactly (and more cheaply) by
  // `dr_representation = integrated`. The `psd` representation keeps the switch,
  // because there it buys a decay-only cross-check against that integrated scheme.
  class_test_severe(in.get<std::string>("inverse_decays").has_value() &&
                        !in.get_flag("inverse_decays", true),
                    "species '%s': dr_representation = proxy is a stand-in for the "
                    "inverse-decay sector and cannot run with inverse_decays = no; use "
                    "dr_representation = integrated for decay-only",
                    name.c_str());
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = in.get_flag("quantum_statistics", true);
  // Stiffness cap, ON here and off in the exact scheme.
  //
  // The per-bin collision rate is K/eps -> a*Gamma once the parent is
  // non-relativistic, and ~90% of a run happens after that -- with the parent's
  // energy fraction below 1e-40, so nothing it does reaches an observable while its
  // eigenvalue still sets the step. Uncapped, a Gamma = 1e7 cell takes 54k rkdp45
  // steps.
  //
  // The cap is the kernel's own device (one global scale on K, so its exact
  // conservation identities survive) and it is CONTINUOUS in a, which matters: an
  // abrupt "switch the collision off below a threshold" makes the RHS discontinuous
  // in the state and an adaptive stepper chatters at the switch instead of stepping
  // over it (measured: rkdp45 stopped converging entirely).
  //
  // The exact scheme leaves it off because capping the reference method is an
  // approximation of unquantified accuracy. The proxy is a cheap stand-in by
  // construction, so it defaults to a value that has been scanned -- see
  // kDefaultRateCap.
  cfg.max_rate = in.get_or("dr_rate_cap", kDefaultRateCap);

  auto fermion = std::make_unique<DarkRadiationSpecies>("dr_" + name + "_l", ctx.pba, ctx.bgm);
  auto boson   = std::make_unique<DarkRadiationSpecies>("dr_" + name + "_phi", ctx.pba, ctx.bgm);

  auto composite        = std::make_unique<DNCDMProxySpecies>(std::move(parent),
                                                              std::move(fermion),
                                                              std::move(boson),
                                                              std::move(q_d),
                                                              std::move(dq_d),
                                                              cfg,
                                                              ctx.pba,
                                                              ctx.bgm);
  composite->f_ini_l_   = in.get_or("dr_f_ini_l", 1.0);
  composite->f_ini_phi_ = in.get_or("dr_f_ini_phi", 0.0);

  // Parsed BEFORE the coefficients below, which default differently per form.
  // Structural (names a closed set of forms, not a value a sampler varies) -> severe.
  // Rejected rather than defaulted: the two forms differ by more than the effects
  // being measured with them, so a silently-ignored typo would be indistinguishable
  // from a physics result.
  const std::string rta_form = in.get<std::string>("dr_rta_form").value_or("powers");
  class_test_severe(rta_form != "powers" && rta_form != "copw",
                    "species '%s': dr_rta_form (='%s') must be 'powers' (fitted, default) "
                    "or 'copw' (analytic, arXiv:2203.09075)",
                    name.c_str(),
                    rta_form.c_str());
  const bool copw      = (rta_form == "copw");
  composite->rta_form_ = copw ? RtaForm::kCOPW : RtaForm::kPowers;

  // Under `copw` these default to the PAPER, not to the fit: eq. (13) is a single
  // gamma^-5 term carrying its own (1/12) amplitude and no free normalisation, and
  // the paper argues against the gamma^-3 rate entirely. So C5 = 1 and C3 = 0 there.
  // Carrying the fitted 2.5228 over would hand 2.52x the published rate to anyone
  // selecting `copw` precisely in order to compare against the paper -- and both
  // remain settable, so nothing is lost by defaulting them honestly.
  //
  // Defaults from a global least-squares fit of
  //   Gamma_T^phys/Gamma = (rho_H/rho_sec) (1-1/gamma^2)^q [C3 eps^n3 gamma^-3 + C5 gamma^-5]
  // to measured Gamma_T/H curves from the moments rung (72 points, 6 cells,
  // Gamma = 1e5..1e10 at m = 0.3, restricted to their own convergence gate
  // gamma <= 20, with
  // eps_ne, rho_H/rho_sec and gamma all read off this code's background). The fit
  // reproduces the measured rate with a 16/84% spread of 0.74-1.46 -- about +-40%.
  //
  // THIS IS A CALIBRATION, NOT A DERIVATION. The gamma^-3 and gamma^-5 powers and
  // the eps dependence come from the physics (Hannestad-Raffelt; Barenboim et al.
  // arXiv:2011.01502 sec.5); C3, C5, n3 and q are fitted numbers, and a single power
  // fits that data as well as the two-term form does, over a window only 0.7 decades
  // wide in gamma. Re-fit them if the measurement
  // improves; do not read them as physical constants.
  //
  // Dropping the shut-off costs a factor ~2 at recombination and ~17 by a = 4e-3 at
  // Gamma = 1e10; dropping the eps dependence (n3 = 0) widens the spread to
  // 0.48-1.72. Both terms earn their place.
  composite->C3_    = in.get_or("dr_rta_C3", copw ? 0.0 : 0.2802);
  composite->C5_    = in.get_or("dr_rta_C5", copw ? 1.0 : 2.5228);
  composite->n3_    = in.get_or("dr_rta_n3", 0.8714);
  composite->vshut_ = in.get_or("dr_rta_vshut", 1.7995);
  return Named{name, std::move(composite)};
}

// ─────────────────────────────────────────────────────────────────────────────
// Background
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMProxySpecies::RegisterBackgroundIndices(int& index_bg) {
  CompositeSpecies::RegisterBackgroundIndices(index_bg);
  const int Nd      = static_cast<int>(q_d_.size());
  index_bg_rate_    = index_bg++;
  index_bg_eps_ne_  = index_bg++;
  index_bg_rho_sec_ = index_bg++;
  index_bg_nu_H_    = index_bg++;
  index_bg_nu_l_    = index_bg++;
  index_bg_nu_phi_  = index_bg++;
  index_bg_n_l_     = index_bg++;
  index_bg_n_phi_   = index_bg++;
  // The per-bin daughter PSDs are NOT published. Every background column is splined
  // and interpolated at EVERY perturbation RHS evaluation, for every k-mode, and the
  // RTA never reads them: it works entirely off the sector's integrated moments. The
  // exact scheme pays that cost and has to -- its perturbations need f-bar per bin.
}

void DNCDMProxySpecies::RegisterIntegrationIndices(int& index_bi) {
  CompositeSpecies::RegisterIntegrationIndices(index_bi);
  const int Nd     = static_cast<int>(q_d_.size());
  index_bi_f_l_    = index_bi;
  index_bi        += Nd;
  index_bi_f_phi_  = index_bi;
  index_bi        += Nd;
}

void DNCDMProxySpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  CompositeSpecies::SetBackgroundInitialConditions(ctx);
  const int Nd = static_cast<int>(q_d_.size());
  for (int i = 0; i < Nd; ++i) {
    const double q = q_d_[i];
    // f_ini_l = 1 seeds nu_l with a FULL Fermi-Dirac at the parent's temperature
    // (the published A = [0,1,1] initial state: phi empty, nu_l thermal). Seeded at
    // the true occupation, never at a large subtracted offset -- f-bar feeds the
    // kernel's (1±f) coefficients directly.
    const double fl = (f_ini_l_ > 0.) ? f_ini_l_ / (std::exp(q) + 1.) : 0.;
    const double fp = (f_ini_phi_ > 0.) ? f_ini_phi_ / std::expm1(q) : 0.;
    ctx.pvecback_integration[index_bi_f_l_ + i]   = fl + kFFloor;
    ctx.pvecback_integration[index_bi_f_phi_ + i] = fp + kFFloor;
  }
  // The daughters' DarkRadiationSpecies rho slots must agree with the PSDs they are
  // driven from, or the sector's energy is double-booked at the first step.
  ctx.pvecback_integration[fermion_->bi_rho_index()] =
      DaughterRho(&ctx.pvecback_integration[index_bi_f_l_], ctx.a_ini, false);
  ctx.pvecback_integration[boson_->bi_rho_index()] =
      DaughterRho(&ctx.pvecback_integration[index_bi_f_phi_], ctx.a_ini, true);
}

double DNCDMProxySpecies::DaughterRho(const double* f, double a, bool boson) const {
  const double fac = BareFactor() * (boson ? 0.5 : 1.0) / (a * a * a * a);
  double sum       = 0.;
  const int Nd     = static_cast<int>(q_d_.size());
  for (int i = 0; i < Nd; ++i) {
    const double q  = q_d_[i];
    sum            += dq_d_[i] * q * q * q * std::max(f[i], 0.);
  }
  return fac * sum;
}

double DNCDMProxySpecies::DaughterNumber(const double* f, bool boson) const {
  double sum   = 0.;
  const int Nd = static_cast<int>(q_d_.size());
  for (int i = 0; i < Nd; ++i) {
    sum += dq_d_[i] * q_d_[i] * q_d_[i] * std::max(f[i], 0.);
  }
  return boson ? 0.5 * sum : sum;
}

void DNCDMProxySpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  CompositeSpecies::ComputeBackground(a, pvecback_B, pvecback);
  ComputeDerivedBackground(a, pvecback_B, pvecback);
}

void DNCDMProxySpecies::ComputeDerivedBackground(double a,
                                                 const double* pvecback_B,
                                                 double* pvecback) {
  pvecback[index_bg_rho_sec_] = parent_->Rho(pvecback) + fermion_->Rho(pvecback) +
                                boson_->Rho(pvecback);

  // Gross decay rate and the departure-from-balance order parameter.
  //
  //   R_gross = sum_q dq q^2 (K/eps) f_H          (pure loss: every parent that
  //                                                decays, ignoring regeneration)
  //   eps_ne  = |dN_H/dtau|_net / R_gross         (1 = free decay, 0 = balance)
  //
  // eps_ne is the perturbation side's only handle on how far the sector is from
  // detailed balance, and it is what multiplies the gamma^-3 term in Gamma_T. It
  // costs one extra quadrature over the parent grid, once per background row.
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const double K     = a * a * m * Gamma;
  const int Np       = parent_->q_size();
  const auto& q      = parent_->GetQ();
  const auto& dq     = parent_->dq();
  const int lnf_idx  = parent_->bg_lnf_index();
  const double kappa = KappaStoredToBare();

  double gross = 0.;
  for (int i = 0; i < Np; ++i) {
    const double eps  = std::sqrt(q[i] * q[i] + a * a * m * m);
    const double fH   = kappa * std::exp(pvecback[lnf_idx + i]);
    gross            += dq[i] * q[i] * q[i] * (K / eps) * fH;
  }
  pvecback[index_bg_rate_] = gross;

  // eps_ne needs the NET parent rate, which needs a kernel pass -- and this
  // function runs on every trial evaluation the evolver makes, four per ETDRK4
  // step. But eps_ne is only ever read back through the interpolated background
  // TABLE, so it only has to be right on stored rows. Computing it on trial rows
  // was ~4 wasted kernel passes per step; gating on StoringBackgroundTable() is
  // exactly the use that flag documents. Trial rows carry the last stored value,
  // which nothing reads.
  double net = eps_ne_cached_;
  if (gross > 0. && (bgm_ == nullptr || bgm_->StoringBackgroundTable())) {
    net = 0.;
    for (int i = 0; i < Np; ++i) {
      fH_bare_[i] = kappa * std::exp(pvecback[lnf_idx + i]);
    }
    kernel_->ComputeBackgroundDerivs(a,
                                     m,
                                     Gamma,
                                     fH_bare_.data(),
                                     &pvecback_B[index_bi_f_l_],
                                     &pvecback_B[index_bi_f_phi_],
                                     df_H_.data(),
                                     df_l_.data(),
                                     df_phi_.data());
    for (int i = 0; i < Np; ++i) {
      net += dq[i] * q[i] * q[i] * df_H_[i];
    }
    net = std::fabs(net) / gross;
    if (net > 1.) {
      net = 1.;
    }
    eps_ne_cached_ = net;
  }
  pvecback[index_bg_eps_ne_] = net;

  // Exchange rates nu_i = (gross transitions) / (comoving number of species i):
  // the fraction of species i the collision reprocesses per unit conformal time.
  // For the parent that is exactly the boosted decay rate a*Gamma*<m/E>; for the
  // daughters it goes to ZERO once the parent is extinct, which is what stops two
  // free-streaming daughters from staying artificially locked for the rest of the
  // run. Computed here, once per background row, rather than in the threaded
  // k-loop.
  const double N_H   = pvecback[parent_->bg_number_index()] / BareFactorNumber(a);
  const double N_l   = DaughterNumber(&pvecback_B[index_bi_f_l_], false);
  const double N_phi = DaughterNumber(&pvecback_B[index_bi_f_phi_], true);
  // Published as physical number densities: WriteBackgroundData needs them and only
  // ever sees pvecback, and the per-bin PSD columns it used to sum are optional.
  pvecback[index_bg_n_l_]    = N_l * BareFactorNumber(a);
  pvecback[index_bg_n_phi_]  = N_phi * BareFactorNumber(a);
  pvecback[index_bg_nu_H_]   = (N_H > 0.) ? gross / N_H : 0.;
  pvecback[index_bg_nu_l_]   = (N_l > 0.) ? gross / N_l : 0.;
  pvecback[index_bg_nu_phi_] = (N_phi > 0.) ? gross / N_phi : 0.;
}

double DNCDMProxySpecies::BareFactorNumber(double a) const {
  // n = BareFactor()/a^3 * sum dq q^2 f, so dividing the published number column by
  // this returns the bare grid moment the gross rate is expressed in.
  return BareFactor() / (a * a * a);
}

void DNCDMProxySpecies::ApplyKernelBackgroundDerivs(double a,
                                                    const double* y,
                                                    double* dy,
                                                    double a_prime_over_a) {
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const int Np       = parent_->q_size();
  const int Nd       = static_cast<int>(q_d_.size());

  // Kernel boundary (#385): the band factor and the (1±f) coefficients need BARE
  // per-dof occupation. The parent's integration slot is stored occupation times
  // DNCDMSpecies::kFScale; the daughters here already store bare occupation.
  const double kappa = KappaStoredToBare() / DNCDMSpecies::kFScale;
  for (int i = 0; i < Np; ++i) {
    fH_bare_[i] = kappa * y[parent_->bi_f_parent_index() + i];
  }
  kernel_->ComputeBackgroundDerivs(a,
                                   m,
                                   Gamma,
                                   fH_bare_.data(),
                                   &y[index_bi_f_l_],
                                   &y[index_bi_f_phi_],
                                   df_H_.data(),
                                   df_l_.data(),
                                   df_phi_.data(),
                                   a_prime_over_a);

  // Assign (=): in collision-owned mode the parent's own background RHS is a no-op
  // and a comoving PSD has no dilution term, so the kernel source IS the full RHS.
  for (int i = 0; i < Np; ++i) {
    dy[parent_->bi_f_parent_index() + i] = df_H_[i] / kappa;
  }
  for (int i = 0; i < Nd; ++i) {
    dy[index_bi_f_l_ + i]   = df_l_[i];
    dy[index_bi_f_phi_ + i] = df_phi_[i];
  }

  // The daughters' DR energy slots, driven from the SAME transitions, so the
  // sector's energy budget closes to integrator accuracy rather than by
  // construction-plus-hope. DarkRadiationSpecies::BackgroundDerivs has already
  // written the -4(a'/a)rho dilution; this is the source on top of it.
  const double fac = BareFactor() / (a * a * a * a);
  double e_l = 0., e_phi = 0.;
  for (int i = 0; i < Nd; ++i) {
    const double q3  = q_d_[i] * q_d_[i] * q_d_[i] * dq_d_[i];
    e_l             += q3 * df_l_[i];
    e_phi           += q3 * df_phi_[i];
  }
  dy[fermion_->bi_rho_index()] += fac * e_l;
  dy[boson_->bi_rho_index()]   += fac * 0.5 * e_phi;

  // Fingerprint the state this kernel pass was prepared at, so the diagonal can
  // reuse the cached transitions instead of paying for a second pass (see
  // BackgroundDerivsDiagonal). Deterministic: same values, same summation order.
  prep_a_  = a;
  double c = 0.;
  for (int i = 0; i < Np; ++i) {
    c += fH_bare_[i];
  }
  for (int i = 0; i < Nd; ++i) {
    c += y[index_bi_f_l_ + i] + y[index_bi_f_phi_ + i];
  }
  prep_sum_ = c;
}

void DNCDMProxySpecies::BackgroundDerivs(double tau,
                                         const double* y,
                                         double* dy,
                                         const double* pvecback) {
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);
  const double a_bg = pvecback[bgm_->index_bg_a_];
  pvecback_row_     = pvecback;  // read by the extinction cut; valid for this call only
  ApplyKernelBackgroundDerivs(a_bg, y, dy, a_bg * pvecback[bgm_->index_bg_H_]);
  pvecback_row_ = nullptr;
}

void DNCDMProxySpecies::BackgroundDerivsDiagonal(double tau,
                                                 const double* y,
                                                 double* diag,
                                                 const double* pvecback) {
  CompositeSpecies::BackgroundDerivsDiagonal(tau, y, diag, pvecback);

  // The collision's per-bin self-coupling. This is the whole stiffness of the
  // background system: K/eps -> a*Gamma once the parent is non-relativistic, and
  // ~90% of a run happens after that. Reporting it lets the exponential evolver
  // integrate it EXACTLY instead of resolving it with the step size, which is the
  // entire reason that evolver exists.
  //
  // Re-prepares the transitions rather than reading the kernel's cache, for the
  // same reason the exact composite does: a diagonal taken from a DIFFERENT state
  // than the RHS it is paired with is the failure mode that stays stable while
  // converging to the wrong attractor.
  const double a_bg  = pvecback[bgm_->index_bg_a_];
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const int Np       = parent_->q_size();
  const int Nd       = static_cast<int>(q_d_.size());
  const double kappa = KappaStoredToBare() / DNCDMSpecies::kFScale;
  for (int i = 0; i < Np; ++i) {
    fH_bare_[i] = kappa * y[parent_->bi_f_parent_index() + i];
  }
  // Reuse the transitions the RHS just prepared, when the state is provably the
  // same one. The exact composite re-prepares unconditionally, and its reason is
  // sound -- a diagonal taken from a DIFFERENT state than the RHS it is paired with
  // is the failure mode that stays stable while converging to the wrong attractor.
  // The fix is not to trust the evolver's calling order but to CHECK it: a exactly
  // equal and the state checksum exactly equal means the kernel's cache is the same
  // object the RHS produced, so re-preparing would recompute it bit for bit.
  // Worth ~1.25x, because the diagonal is one pass in five per ETDRK4 step.
  double c = 0.;
  for (int i = 0; i < Np; ++i) {
    c += fH_bare_[i];
  }
  for (int i = 0; i < Nd; ++i) {
    c += y[index_bi_f_l_ + i] + y[index_bi_f_phi_ + i];
  }
  if (a_bg != prep_a_ || c != prep_sum_) {
    kernel_->ComputeBackgroundDerivs(a_bg,
                                     m,
                                     Gamma,
                                     fH_bare_.data(),
                                     &y[index_bi_f_l_],
                                     &y[index_bi_f_phi_],
                                     df_H_.data(),
                                     df_l_.data(),
                                     df_phi_.data(),
                                     a_bg * pvecback[bgm_->index_bg_H_]);
    prep_a_   = a_bg;
    prep_sum_ = c;
  }
  kernel_->CollisionDiagonal(diag_H_.data(), diag_l_.data(), diag_phi_.data());

  // d(dy_i)/d(y_i) is INVARIANT under the parent's stored/bare rescaling: the write
  // back divides the rate by the same kappa the read multiplies the occupation by.
  for (int i = 0; i < Np; ++i) {
    diag[parent_->bi_f_parent_index() + i] += diag_H_[i];
  }
  for (int i = 0; i < Nd; ++i) {
    diag[index_bi_f_l_ + i]   += diag_l_[i];
    diag[index_bi_f_phi_ + i] += diag_phi_[i];
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbation Jacobian diagonal
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMProxySpecies::PerturbDerivsDiagonal(const BaseSpecies::PerturbLayout& base,
                                              double tau,
                                              const double* y,
                                              double* diag,
                                              const perturb_parameters_and_workspace& ppaw) const {
  CompositeSpecies::PerturbDerivsDiagonal(base, tau, y, diag, ppaw);

  // The RTA is diagonal by construction: every state variable carries a term
  // -(nu_i + Gamma_T,l) y, and the only off-diagonal piece is the target contrast,
  // which is a source. Reporting that diagonal lets the exponential evolver
  // integrate the relaxation EXACTLY instead of resolving it with the step size,
  // which is what the etd scheme exists for.
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* pvecback          = ppw->pvecback.data();

  const NCDMBaseSpecies::PerturbLayout& p_lay      = parent_layout(base);
  const DarkRadiationSpecies::PerturbLayout& f_lay = dr_layout(base, kFermion);
  const DarkRadiationSpecies::PerturbLayout& b_lay = dr_layout(base, kBoson);
  const int Np                                     = p_lay.q_size;
  const int l_max = std::min({p_lay.l_max, f_lay.l_max, b_lay.l_max});
  if (Np <= 0 || l_max < 0 || l_max >= kLMaxScratch) {
    return;
  }

  const double a         = ctx.a;
  const double rate_cap  = kRateCapFactor * std::max(ctx.k, a * pvecback[bgm_->index_bg_H_]);
  auto capped            = [rate_cap](double v) { return (v > rate_cap) ? rate_cap : v; };
  const double rho_sec   = pvecback[index_bg_rho_sec_];
  const bool parent_live = (rho_sec > 0.) && (parent_->Rho(pvecback) > kParentFloor * rho_sec);

  const double nu_H   = parent_live ? capped(pvecback[index_bg_nu_H_]) : 0.;
  const double nu_l   = capped(pvecback[index_bg_nu_l_]);
  const double nu_phi = capped(pvecback[index_bg_nu_phi_]);

  for (int l = 0; l <= l_max; ++l) {
    const double GT         = capped(TransportRate(l, pvecback));
    diag[f_lay.idx_F0 + l] += -(nu_l + GT);
    diag[b_lay.idx_F0 + l] += -(nu_phi + GT);
    if (nu_H > 0. || GT > 0.) {
      for (int iq = 0; iq < Np; ++iq) {
        diag[p_lay.index_per_q[iq] + l] += -(nu_H + GT);
      }
    }
  }
}

void DNCDMProxySpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  // Column layout deliberately mirrors DNCDMInvSpecies: parent, then the two
  // daughters' number/rho/p. A proxy run and an exact run are then readable by the
  // same analysis script.
  parent_->WriteBackgroundColumnTitles(w);
  const std::string& nm = name();
  for (int leg = 0; leg < 2; ++leg) {
    const std::string dn = "dr_" + nm + (leg == 0 ? "_l" : "_phi");
    w.Add("(.)number_" + dn, 0.);
    w.Add("(.)rho_" + dn, 0.);
    w.Add("(.)p_" + dn, 0.);
  }
  w.Add("proxy_rate_" + nm, 0.);
  w.Add("proxy_eps_ne_" + nm, 0.);
  w.Add("proxy_gamma_" + nm, 0.);
}

void DNCDMProxySpecies::WriteBackgroundData(const double* pvecback,
                                            BackgroundColumnWriter& w) const {
  parent_->WriteBackgroundData(pvecback, w);
  const std::string& nm = name();
  for (int leg = 0; leg < 2; ++leg) {
    const bool boson     = (leg == 1);
    const std::string dn = "dr_" + nm + (boson ? "_phi" : "_l");
    const double rho     = boson ? boson_->Rho(pvecback) : fermion_->Rho(pvecback);
    w.Add("(.)number_" + dn, pvecback[boson ? index_bg_n_phi_ : index_bg_n_l_]);
    w.Add("(.)rho_" + dn, rho);
    w.Add("(.)p_" + dn, rho / 3.);  // massless
  }
  const double a     = pvecback[bgm_->index_bg_a_];
  const double X     = a * parent_->GetMass();
  const double n     = pvecback[parent_->bg_number_index()];
  const double rho   = parent_->Rho(pvecback);
  const double gamma = (n > 0. && X > 0.) ? (a * rho / n) / X : 0.;
  w.Add("proxy_rate_" + nm, pvecback[index_bg_rate_]);
  w.Add("proxy_eps_ne_" + nm, pvecback[index_bg_eps_ne_]);
  w.Add("proxy_gamma_" + nm, gamma);
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbations: the RTA closure
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMProxySpecies::AddCouplingDerivs(double /*tau*/,
                                          const double* y,
                                          double* dy,
                                          const perturb_parameters_and_workspace& ppaw) const {
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* pvecback          = ppw->pvecback.data();
  const double a                  = ctx.a;

  const BaseSpecies::PerturbLayout& my             = *ppw->pv->species_layouts[collection_index_];
  const NCDMBaseSpecies::PerturbLayout& p_lay      = parent_layout(my);
  const DarkRadiationSpecies::PerturbLayout& f_lay = dr_layout(my, kFermion);
  const DarkRadiationSpecies::PerturbLayout& b_lay = dr_layout(my, kBoson);

  const int Np    = p_lay.q_size;
  const int l_max = std::min({p_lay.l_max, f_lay.l_max, b_lay.l_max});
  if (Np <= 0 || l_max < 0 || l_max >= kLMaxScratch) {
    return;
  }

  // ── rates ─────────────────────────────────────────────────────────────────
  // The exchange rates are precomputed per background row (ComputeDerivedBackground);
  // all that is left here is the k-dependent ceiling, which every rate shares.
  const double rate_cap = kRateCapFactor * std::max(ctx.k, a * pvecback[bgm_->index_bg_H_]);
  auto capped           = [rate_cap](double v) { return (v > rate_cap) ? rate_cap : v; };

  const double rho_H   = parent_->Rho(pvecback);
  const double rho_sec = pvecback[index_bg_rho_sec_];
  // Once the parent is a 1e-60 fraction of the sector, nothing it does reaches an
  // observable, but its exchange rate is still a*Gamma -- the single stiffest
  // eigenvalue in the whole system, and one the explicit evolver would have to
  // resolve step by step for the ~90% of the run that happens after extinction.
  // Dropping the parent leg outright there is not an approximation with an error
  // budget; it is arithmetic on a weight that has already underflowed the metric.
  const bool parent_live = (rho_sec > 0.) && (rho_H > kParentFloor * rho_sec);

  const double nu_H   = parent_live ? capped(pvecback[index_bg_nu_H_]) : 0.;
  const double nu_l   = capped(pvecback[index_bg_nu_l_]);
  const double nu_phi = capped(pvecback[index_bg_nu_phi_]);
  if (nu_H <= 0. && nu_l <= 0. && nu_phi <= 0.) {
    return;  // collision inactive; the target contrast is undefined and unneeded
  }

  const double a_sq_over_H0 = (a * a) / pba_->H0;
  const double r_l          = fermion_->Rho(pvecback) * a_sq_over_H0 * a_sq_over_H0;
  const double r_phi        = boson_->Rho(pvecback) * a_sq_over_H0 * a_sq_over_H0;

  // ── the parent's moments ──────────────────────────────────────────────────
  // A_l = the state's contribution to the l-th stress-energy moment; B_l = what
  // that moment would be at unit contrast, i.e. for the separable profile
  // F_l(q) = -(1/4) df-bar/dlnq. Their ratio is the parent's contrast.
  //
  // The parent carries the UNNORMALISED F = delta-f in collision-owned mode (#386),
  // so A_l has no f-bar of its own and B_l carries df-bar/dlnq = f-bar dlnf/dlnq.
  //
  // q outer / l inner, with (q/eps)^l built by one multiply per step: the paper's
  // weight costs a pow() per (q, l) if written directly, which at 16 bins x 18
  // multipoles is ~300 pow calls on every RHS evaluation of every k-mode.
  double A_H[kLMaxScratch], B_H[kLMaxScratch];
  for (int l = 0; l <= l_max; ++l) {
    A_H[l] = 0.;
    B_H[l] = 0.;
  }
  const auto& q      = parent_->GetQ();
  const auto& dq     = parent_->dq();
  const double aM    = a * parent_->GetMass();
  const int lnf_idx  = parent_->bg_lnf_index();
  const int dlnf_idx = parent_->bg_dlnfdlnq_index();
  const double Cnorm = parent_->factor() / (pba_->H0 * pba_->H0);

  if (nu_H > 0.) {
    for (int iq = 0; iq < Np; ++iq) {
      const double qq    = q[iq];
      const double eps   = std::sqrt(qq * qq + aM * aM);
      const double dfdln = std::exp(pvecback[lnf_idx + iq]) * pvecback[dlnf_idx + iq];
      const double ratio = qq / eps;
      const int base     = p_lay.index_per_q[iq];
      double w           = dq[iq] * qq * qq * eps;  // (q/eps)^0
      for (int l = 0; l <= l_max; ++l) {
        A_H[l] += w * y[base + l];
        B_H[l] += w * (-0.25) * dfdln;
        w      *= ratio;
      }
    }
    for (int l = 0; l <= l_max; ++l) {
      A_H[l] *= Cnorm;
      B_H[l] *= Cnorm;
    }
  }

  // ── per-multipole relaxation ──────────────────────────────────────────────
  double target[kLMaxScratch], GT[kLMaxScratch];
  for (int l = 0; l <= l_max; ++l) {
    const double A_l   = y[f_lay.idx_F0 + l];
    const double A_phi = y[b_lay.idx_F0 + l];
    // The nu*B-weighted mean. This choice is what makes the l = 0 and l = 1
    // exchange EXACTLY conservative in delta-rho and (rho+P)theta -- not optional,
    // because the Einstein equations see the sector's total and a relaxation that
    // leaked into it would source the metric out of nothing.
    const double den = nu_H * B_H[l] + nu_l * r_l + nu_phi * r_phi;
    target[l]        = (den > 0.) ? (nu_H * A_H[l] + nu_l * A_l + nu_phi * A_phi) / den : 0.;
    // The damping is capped for the same reason the exchange is: alpha_l grows as
    // l^4 (8041 at l = 17), so the top of the hierarchy would otherwise carry an
    // eigenvalue thousands of times the collision rate. A damping already many
    // times faster than k has wiped the moment out; making it faster still changes
    // nothing but the step size.
    GT[l] = capped(TransportRate(l, pvecback));
  }

  for (int l = 0; l <= l_max; ++l) {
    // Daughters: A = F_l and B = r_i by CLASS's own DR convention (its adiabatic
    // seed is F_l = contrast * r_dr for every l), so the contrast equation maps
    // straight onto the state.
    const double F_l      = y[f_lay.idx_F0 + l];
    const double F_phi    = y[b_lay.idx_F0 + l];
    dy[f_lay.idx_F0 + l] += -nu_l * (F_l - target[l] * r_l) - GT[l] * F_l;
    dy[b_lay.idx_F0 + l] += -nu_phi * (F_phi - target[l] * r_phi) - GT[l] * F_phi;
  }

  if (nu_H > 0.) {
    // Parent: relax the per-q profile toward the separable target shape. Integrated
    // against the same weight this returns -nu_H(A_l - B_l*target) - GT*A_l, i.e.
    // exactly the contrast equation above -- which is why the conservation identity
    // holds bin by bin and not only on average.
    for (int iq = 0; iq < Np; ++iq) {
      const double dfdln = std::exp(pvecback[lnf_idx + iq]) * pvecback[dlnf_idx + iq];
      const int base     = p_lay.index_per_q[iq];
      for (int l = 0; l <= l_max; ++l) {
        const double F_tgt  = -0.25 * dfdln * target[l];
        dy[base + l]       += -nu_H * (y[base + l] - F_tgt) - GT[l] * y[base + l];
      }
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbation output (children write their own transfer columns)
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMProxySpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                           const PerturbationsModule& mod,
                                           file_format fmt,
                                           BaseSpecies::TransferColumnSection section) const {
  parent_->WriteOutputColumns(w, mod, fmt, section);
  fermion_->WriteOutputColumns(w, mod, fmt, section);
  boson_->WriteOutputColumns(w, mod, fmt, section);
}

// ─────────────────────────────────────────────────────────────────────────────
// k_output time series
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMProxySpecies::PrintVariables(PerturbColumnWriter& w,
                                       double /*tau*/,
                                       const double* y,
                                       const PerturbationsModule& mod,
                                       const perturb_workspace* ppw) const {
  // a⁴Π_νφ ≡ a⁴·Σ_{i∈{H,l,φ}} (ρ̄_i+p̄_i)σ_i, the sector's total anisotropic
  // stress and the quantity the ℓ ≥ 2 closure approximates. Each child already
  // forms (ρ+p)σ in its own StressEnergy — the parent from its per-q moments,
  // the daughters from their integrated F_2 — so this only sums the three and
  // multiplies by a⁴. Column names match DNCDMInvSpecies exactly, which is the
  // point: it lets one script read a proxy run against an exact one.
  double pi_H = 0., pi_l = 0., pi_phi = 0.;
  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const double a4          = a * a * a * a;
    const auto& my           = static_cast<const CompositeSpecies::PerturbLayout&>(
        *pv->species_layouts[collection_index_]);
    pi_H = a4 *
           parent_->StressEnergy(*my.child_layouts[kParent], pv, y, pvecback, ppw).rho_plus_p_shear;
    pi_l =
        a4 *
        fermion_->StressEnergy(*my.child_layouts[kFermion], pv, y, pvecback, ppw).rho_plus_p_shear;
    pi_phi = a4 *
             boson_->StressEnergy(*my.child_layouts[kBoson], pv, y, pvecback, ppw).rho_plus_p_shear;
  }
  w.Add("a4Pi_" + name(), pi_H + pi_l + pi_phi, true);
  w.Add("a4Pi_" + name() + "_H", pi_H, true);
  w.Add("a4Pi_" + name() + "_l", pi_l, true);
  w.Add("a4Pi_" + name() + "_phi", pi_phi, true);
}
