#include "wdm_decay_product.h"

#include <algorithm>
#include <cmath>

#include "background_module.h"
#include "dcdm.h"
#include "errors.h"
#include "perturbations_module.h"
#include "precision.h"
#include "species/species_input.h"

// ─────────────────────────────────────────────────────────────────────────────
// Fiducial cosmic time t(a) and its inverse (placement-only; NOT physics).
// H(a) = H0 a^-2 sqrt(Or + Om a)  =>  exact radiation+matter cosmic time
//   t(a) = (2 / (3 H0 Om^2)) [ sqrt(Or + Om a)(Om a - 2 Or) + 2 Or^{3/2} ].
// (Lambda is negligible at the epochs the daughter grid cares about.)
// ─────────────────────────────────────────────────────────────────────────────

double WdmDecayProductSpecies::FiducialCosmicTime(double a) {
  const double H0 = kFidH * 1.0e5 / _c_;  // Mpc^-1
  const double Om = kFidOmegaM, Or = kFidOmegaR;
  const double s = std::sqrt(Or + Om * a);
  return (2.0 / (3.0 * H0 * Om * Om)) * (s * (Om * a - 2.0 * Or) + 2.0 * Or * std::sqrt(Or));
}

double WdmDecayProductSpecies::FiducialScaleFactorAtTime(double t_mpc) {
  const double H0 = kFidH * 1.0e5 / _c_;
  const double Om = kFidOmegaM, Or = kFidOmegaR;
  double a_lo = 1e-12, a_hi = 1.0;
  while (FiducialCosmicTime(a_hi) < t_mpc && a_hi < 1e3)
    a_hi *= 2.0;
  // Matter-limit seed; refine with bracketed Newton (dt/da = a / (H0 sqrt(Or+Om a)))
  // and bisection fallback — monotone, so this always converges.
  double a = std::pow(1.5 * t_mpc * H0 * std::sqrt(Om), 2.0 / 3.0);
  if (!(a > a_lo && a < a_hi))
    a = 0.5 * (a_lo + a_hi);
  for (int it = 0; it < 80; ++it) {
    const double f = FiducialCosmicTime(a) - t_mpc;
    if (f > 0.)
      a_hi = a;
    else
      a_lo = a;
    const double dtda = a / (H0 * std::sqrt(Or + Om * a));
    double a_next     = a - f / dtda;
    if (!(a_next > a_lo && a_next < a_hi))
      a_next = 0.5 * (a_lo + a_hi);
    if (std::fabs(a_next - a) <= 1e-14 * a) {
      a = a_next;
      break;
    }
    a = a_next;
  }
  return a;
}

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

WdmDecayProductSpecies::WdmDecayProductSpecies(FileContent* pfc,
                                               const std::string& instance_name,
                                               const NcdmSettings& settings,
                                               const background* pba,
                                               const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings, DeferInit{}),
      pba_(pba) {
  bgm_ = bgm;
  SpeciesInput input(pfc, instance_name);

  // NCDM-family keys with no meaning here: the distribution is built from the
  // decay, not from a thermal PSD. Reject loudly instead of silently ignoring.
  for (const char* key : {"m",
                          "T",
                          "deg",
                          "ksi",
                          "Omega",
                          "omega",
                          "max_q",
                          "quadrature_strategy",
                          "momenta_bins_bg",
                          "use_psd_file",
                          "psd_filename"}) {
    class_test_severe(input.get<std::string>(key).has_value(),
                      "species '%s': parameter '%s' is not supported by dcdm_wdm",
                      instance_name.c_str(),
                      key);
  }

  // Internal thermodynamic convention: T_wdm = T_cmb, deg = 1, no chemical potential.
  T_       = 1.0;
  ksi_     = 0.;
  m_in_eV_ = 0.;         // the daughter mass in eV is not physical (kQKick is a convention)
  SetDegAndFactor(1.0);  // recomputes factor_ with T_ = 1

  // Decay rate: exactly one of Gamma [km/s/Mpc], log10Gamma, lifetime [yr], log10lifetime.
  auto Gamma_value         = input.get<double>("Gamma");
  auto log10Gamma_value    = input.get<double>("log10Gamma");
  auto lifetime_value      = input.get<double>("lifetime");
  auto log10lifetime_value = input.get<double>("log10lifetime");
  int n_provided = (int) Gamma_value.has_value() + (int) log10Gamma_value.has_value() +
                   (int) lifetime_value.has_value() + (int) log10lifetime_value.has_value();
  class_test_severe(n_provided != 1,
                    "species '%s': specify exactly one of Gamma, log10Gamma, lifetime, "
                    "log10lifetime",
                    instance_name.c_str());
  double Gamma_raw = 0.;
  if (Gamma_value) {
    Gamma_raw = *Gamma_value;
  }
  else if (log10Gamma_value) {
    Gamma_raw = std::pow(10., *log10Gamma_value);
  }
  else if (lifetime_value) {
    Gamma_raw = 1. / *lifetime_value / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  else {  // log10lifetime
    Gamma_raw = 1. / std::pow(10., *log10lifetime_value) / (365 * 24 * 60 * 60) * _Mpc_over_m_ *
                1e-3;
  }
  Gamma_ = Gamma_raw * (1.e3 / _c_);

  // Kinematics: exactly one of epsilon (mass retention) or vkick (v/c).
  auto eps_opt   = input.get<double>("epsilon");
  auto vkick_opt = input.get<double>("vkick");
  class_test_severe(eps_opt.has_value() == vkick_opt.has_value(),
                    "species '%s': specify exactly one of epsilon, vkick",
                    instance_name.c_str());
  if (eps_opt) {
    eps_ = *eps_opt;
    class_test(eps_ < 0. || eps_ > 1. - 1e-9,
               "species '%s': epsilon must be in [0, 1 - 1e-9]; for smaller kicks input vkick "
               "directly",
               instance_name.c_str());
    vkick_ = std::sqrt((1. - eps_) * (1. + eps_));
  }
  else {
    vkick_ = *vkick_opt;
    class_test(vkick_ < 1e-8 || vkick_ > 1.,
               "species '%s': vkick must be in [1e-8, 1] (in units of c); for vkick -> 0 the "
               "model is plain CDM",
               instance_name.c_str());
    eps_ = std::sqrt((1. - vkick_) * (1. + vkick_));
  }
  // p/m_d = v/sqrt(1-v^2) = v/eps  =>  M = m_d/T0 = kQKick * eps / v.
  M_ = kQKick * eps_ / vkick_;

  // ── Injection-adapted momentum grid (equal injected decays per bin) ─────────
  // Decay is Poisson in cosmic time, so cumulative decayed fraction F = 1-e^{-Gamma t}
  // is the natural coordinate: uniform-in-F bins carry equal injected decays (= equal
  // injected energy). q = a'(F) * kQKick via the fiducial t(a) map => grid = f(Gamma).
  n_bins_      = input.get_or("momenta_bins", 96);
  q_edge_tol_  = input.get_or("q_edge_tol", 1e-3);
  l_max_input_ = input.get_or("l_max", -1);
  if (auto qmr = input.get<double>("q_min_ratio"))
    q_min_ratio_floor_ = *qmr;
  class_test(n_bins_ < 8 || n_bins_ > 4096,
             "species '%s': momenta_bins out of range",
             instance_name.c_str());
  class_test(!(q_edge_tol_ > 1e-12 && q_edge_tol_ < 1e-2),
             "species '%s': q_edge_tol must be in (1e-12, 1e-2)",
             instance_name.c_str());
  if (q_min_ratio_floor_)
    class_test(!(*q_min_ratio_floor_ > 0. && *q_min_ratio_floor_ < 0.5),
               "species '%s': q_min_ratio must be in (0, 0.5)",
               instance_name.c_str());

  BuildInjectionAdaptedGrid();

  // Normalization: Omega_ini/omega_ini (parent initial abundance, dcdm convention)
  // XOR Omega_dcdmwdm/omega_dcdmwdm (combined sector today). During a shooting
  // build DoShooting legitimately writes <instance>.Omega_dcdmwdm alongside the
  // user's initial key.
  auto Omega_ini_opt  = input.get<double>("Omega_ini");
  auto omega_ini_opt  = input.get<double>("omega_ini");
  auto Omega_comb_opt = input.get<double>("Omega_dcdmwdm");
  auto omega_comb_opt = input.get<double>("omega_dcdmwdm");
  class_test_severe(Omega_ini_opt && omega_ini_opt,
                    "species '%s': specify only one of Omega_ini, omega_ini",
                    instance_name.c_str());
  class_test_severe(Omega_comb_opt && omega_comb_opt,
                    "species '%s': specify only one of Omega_dcdmwdm, omega_dcdmwdm",
                    instance_name.c_str());
  const bool has_initial  = Omega_ini_opt.has_value() || omega_ini_opt.has_value();
  const bool has_combined = Omega_comb_opt.has_value() || omega_comb_opt.has_value();
  const bool in_shooting  = (pfc != nullptr) && pfc->is_shooting;
  class_test_severe(has_initial && has_combined && !in_shooting,
                    "species '%s': Omega_ini conflicts with Omega_dcdmwdm — choose one",
                    instance_name.c_str());
  class_test_severe(!has_initial && !has_combined,
                    "species '%s': specify Omega_ini/omega_ini or Omega_dcdmwdm/omega_dcdmwdm",
                    instance_name.c_str());
  if (Omega_ini_opt)
    Omega_ini_pending_ = *Omega_ini_opt;
  else if (omega_ini_opt)
    Omega_ini_pending_ = *omega_ini_opt / settings.h / settings.h;
  if (Omega_comb_opt)
    Omega_combined_pending_ = *Omega_comb_opt;
  else if (omega_comb_opt)
    Omega_combined_pending_ = *omega_comb_opt / settings.h / settings.h;
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildInjectionAdaptedGrid / SigmaAt / GaussWeights
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::BuildInjectionAdaptedGrid() {
  const int N = n_bins_;
  // Today's decayed fraction in the fiducial background caps the upper edge.
  const double g_at_1 = Gamma_ * FiducialCosmicTime(1.0);
  const double F_at_1 = -std::expm1(-g_at_1);
  double F_lo         = q_edge_tol_;
  double F_hi         = std::min(1.0 - q_edge_tol_, F_at_1);
  if (q_min_ratio_floor_) {  // hard floor on q_lo == floor on a' raises F_lo
    const double g_floor = Gamma_ * FiducialCosmicTime(*q_min_ratio_floor_);
    F_lo                 = std::max(F_lo, -std::expm1(-g_floor));
  }

  F_lo_ = F_lo;
  F_hi_ = F_hi;

  q_.resize(N);
  u_.resize(N);
  dq_.resize(N);
  u_edge_.resize(N + 1);

  if (!(F_hi > F_lo + q_edge_tol_)) {
    // Negligible-decay fallback: too little decays by today to define a span.
    // Daughter is essentially empty; a minimal valid uniform-ln q grid keeps the
    // hierarchy/kernel well-formed. Observationally irrelevant.
    F_lo_             = 0.;
    F_hi_             = 0.;  // signals SigmaAt to use the uniform fallback spacing
    const double q_lo = kQKick / 30.0;
    const double du   = std::log(kQKick / q_lo) / N;
    for (int j = 0; j <= N; ++j)
      u_edge_[j] = std::log(q_lo) + j * du;
    for (int i = 0; i < N; ++i)
      u_[i] = 0.5 * (u_edge_[i] + u_edge_[i + 1]);  // uniform-grid cell centre
  }
  else {
    const double dF = (F_hi - F_lo) / N;
    auto u_at_F     = [&](double F) {
      const double g = -std::log1p(-F);  // -ln(1 - F) = Gamma * t
      return std::log(FiducialScaleFactorAtTime(g / Gamma_) * kQKick);
    };
    for (int j = 0; j <= N; ++j)
      u_edge_[j] = u_at_F(F_lo + j * dF);
    // Grid point at the F-quantile midpoint: q_[i] = q(F_mid[i]) (spec §4.2), NOT the
    // u-midpoint of the edges — F is nonlinear in u, so the two differ at O((dF)^2).
    for (int i = 0; i < N; ++i)
      u_[i] = u_at_F(F_lo + (i + 0.5) * dF);
  }

  for (int i = 0; i < N; ++i) {
    q_[i]  = std::exp(u_[i]);
    dq_[i] = std::exp(u_edge_[i + 1]) - std::exp(u_edge_[i]);
  }

  q_bg_ = q_;
  w_.assign(N, 0.);
  w_bg_.assign(N, 0.);
  dlnf0_dlnq_.assign(N, 0.);
  scratch_w_.assign(N, 0.);
}

double WdmDecayProductSpecies::SigmaAt(double u_cut) const {
  // Fallback (negligible-decay) grid: uniform spacing, constant width.
  if (!(F_hi_ > F_lo_)) {
    const double du = u_edge_[1] - u_edge_[0];
    return std::min(kSigmaMax, std::max(kSigmaMin, du));
  }
  // Analytic local quantile-bin width (spec §4.3 v3):
  //   dF/du = Γ e^{−Γ t_fid(a)} / H_fid(a),  a = e^{u_cut}/kQKick
  //   Δu_loc = [(F_hi−F_lo)/N] · H_fid(a) e^{+Γ t_fid(a)} / Γ
  // Γ t_fid can exceed 700 (fast decays, late times): compare in log form and
  // return the cap BEFORE exponentiating — never materialize inf (-ffast-math).
  const double a     = std::exp(u_cut) / kQKick;
  const double H0    = kFidH * 1.0e5 / _c_;
  const double H     = H0 * std::sqrt(kFidOmegaM / (a * a * a) + kFidOmegaR / (a * a * a * a));
  const double gt    = Gamma_ * FiducialCosmicTime(a);
  const double ln_du = std::log((F_hi_ - F_lo_) / n_bins_ * H / Gamma_) + gt;
  if (ln_du >= std::log(kSigmaMax))
    return kSigmaMax;
  return std::max(kSigmaMin, std::exp(ln_du));
}

void WdmDecayProductSpecies::GaussWeights(double u_cut, double sigma, double* w) const {
  const int N        = n_bins_;
  const double inv_s = 1.0 / (sigma * std::sqrt(2.0));
  // Normal CDF at a cell edge; erf saturates by |x| ~ 6, clamp at 8 so no
  // huge-argument libm call happens under -ffast-math.
  auto cdf = [&](double u_e) {
    const double x = (u_e - u_cut) * inv_s;
    if (x >= 8.)
      return 1.0;
    if (x <= -8.)
      return 0.0;
    return 0.5 * (1.0 + std::erf(x));
  };
  double c_lo = cdf(u_edge_[0]);
  w[0]        = c_lo;  // below-grid tail clamped into the bottom bin
  for (int i = 0; i < N; ++i) {
    const double c_hi = cdf(u_edge_[i + 1]);
    if (i > 0)
      w[i] = 0.;
    w[i] += c_hi - c_lo;
    c_lo  = c_hi;
  }
  w[N - 1] += 1.0 - c_lo;  // above-grid tail clamped into the top bin
  // Σ w = 1 identically for every u_cut (telescoping) — the exact energy
  // sum-rule needs no renormalization.
}

// ─────────────────────────────────────────────────────────────────────────────
// FillInjection
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::FillInjection(double a,
                                           double rho_dcdm,
                                           double* J,
                                           double* dJdlnq) const {
  const int N        = n_bins_;
  const double u_cut = std::log(a * kQKick);

  auto zero_all = [&]() {
    for (int i = 0; i < N; ++i) {
      J[i] = 0.;
      if (dJdlnq != nullptr)
        dJdlnq[i] = 0.;
    }
  };
  if (rho_dcdm <= 0. || Gamma_ == 0.) {
    zero_all();
    return;
  }

  // Conservative cell-integrated Gaussian (erf) energy deposit (spec §4.3 v3):
  // smooth in time (both evolvers; tabulated injection columns stay
  // representable), Σ w = 1 identically (exact sum rule, edges clamped).
  const double a2    = a * a;
  const double A     = a * Gamma_ * rho_dcdm * a2 * a2 / factor_;
  const double sigma = SigmaAt(u_cut);
  GaussWeights(u_cut, sigma, scratch_w_.data());

  // Static-sparsity floor (see kSparsityFloor doc): w -> (w + f)/(1 + N f)
  // keeps Sum w = 1 exactly while making every row structurally nonzero.
  const double fnorm = 1.0 / (1.0 + N * kSparsityFloor);
  for (int i = 0; i < N; ++i)
    scratch_w_[i] = (scratch_w_[i] + kSparsityFloor) * fnorm;

  for (int i = 0; i < N; ++i) {
    const double w       = scratch_w_[i];
    const double epsilon = std::sqrt(q_[i] * q_[i] + a2 * M_ * M_);
    J[i]                 = A * w / (dq_[i] * q_[i] * q_[i] * epsilon);
  }

  if (dJdlnq == nullptr)
    return;
  // g-source: d/dlnq of the deposited profile J(u) = A·φ_σ(u−u_cut)/(q³ε). The
  // 1/(q³ε) measure factor makes this a TWO-term derivative, not just the
  // edge-pdf difference:
  //   ∂J/∂ln q = [A·φ′_σ/(q³ε)]  −  (3 + q²/ε²)·J
  // First term: cell-averaged d/dlnq of the raw kernel, via the normal pdf
  // differenced at the cell edges (no edge special cases; telescopes to ~0
  // under the Σ dq q² ε (·) moment for interior cutoffs). Second term: the
  // measure's own log-derivative (d ln(q³ε)/d ln q = 3 + q²/ε²) pulled back
  // through the quotient rule; this is what supplies the correct ρ-weighted
  // moment Σ dq q² ε dJdlnq ≈ −Σ dq q² ε (3 + q²/ε²) J ≈ −3A (non-relativistic,
  // via integration by parts) — the g-channel's fluid-limit continuity source.
  // Omitting it leaves the g-channel's metric-driven growth dead (measured: a
  // flat −7.6% low-k P_m deficit; see the spec's §4.3 dJ/dlnq block).
  const double inv_s = 1.0 / sigma;
  const double norm  = inv_s / std::sqrt(2.0 * M_PI);
  auto pdf           = [&](double u_e) {
    const double x  = (u_e - u_cut) * inv_s;
    const double x2 = 0.5 * x * x;
    return (x2 >= 60.) ? 0.0 : norm * std::exp(-x2);  // master-style exp clamp
  };
  double p_lo = pdf(u_edge_[0]);
  for (int i = 0; i < N; ++i) {
    const double p_hi         = pdf(u_edge_[i + 1]);
    const double diff         = p_hi - p_lo;
    p_lo                      = p_hi;
    const double epsilon      = std::sqrt(q_[i] * q_[i] + a2 * M_ * M_);
    const double q2_over_eps2 = q_[i] * q_[i] / (epsilon * epsilon);
    // + kSparsityFloor: static-sparsity seed (sign irrelevant, structure only)
    dJdlnq[i] = A * (diff + kSparsityFloor) / (dq_[i] * q_[i] * q_[i] * epsilon) -
                (3. + q2_over_eps2) * J[i];
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Info prints
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::PrintNeffInfo() const {
  printf(
      " -> dcdm_wdm daughter '%s': %d momentum bins, vkick = %g c (epsilon = %g), "
      "Gamma = %g Mpc^-1\n",
      name().c_str(),
      n_bins_,
      vkick_,
      eps_,
      Gamma_);
}

void WdmDecayProductSpecies::PrintMassInfo() const {
  // The daughter mass in eV is not physical here (internal kQKick convention);
  // report the physical kinematics instead.
  printf(" -> dcdm_wdm daughter '%s': m_d/m_parent = %g, vkick = %g c\n",
         name().c_str(),
         eps_ / 2.,
         vkick_);
}

// ─────────────────────────────────────────────────────────────────────────────
// Background (Task 2) — placeholders so the file links.
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_    = index_bg++;
  index_bg_rho_       = index_bg++;
  index_bg_p_         = index_bg++;
  index_bg_pseudo_p_  = index_bg++;
  index_bg_f_         = index_bg;
  index_bg           += n_bins_;
  index_bg_dfdlnq_    = index_bg;
  index_bg           += n_bins_;
  index_bg_inj_       = index_bg;  // filled by DCDM_WDM_Species (needs rho_dcdm)
  index_bg           += n_bins_;
}

void WdmDecayProductSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_f_       = index_bi;
  index_bi         += n_bins_;
  index_bi_dfdlnq_  = index_bi;
  index_bi         += n_bins_;
}

void WdmDecayProductSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  for (int i = 0; i < n_bins_; ++i) {
    // f starts at kFSeed, not exactly 0 (see the kFSeed doc comment: ndf15's
    // error-weight floor makes exactly-zero components force micro-steps at
    // each bin's injection onset). Subtracted again in ComputeBackground.
    ctx.pvecback_integration[index_bi_f_ + i]      = kFSeed;
    ctx.pvecback_integration[index_bi_dfdlnq_ + i] = 0.;
  }
}

void WdmDecayProductSpecies::ComputeBackground(double a,
                                               const double* pvecback_B,
                                               double* pvecback) {
  const double z = 1. / a - 1.;
  for (int i = 0; i < n_bins_; ++i) {
    // Subtract the kFSeed IC offset, then clamp tiny negative excursions from
    // the ODE (f is a density).
    const double f                 = std::max(pvecback_B[index_bi_f_ + i] - kFSeed, 0.);
    w_bg_[i]                       = f * dq_[i];
    pvecback[index_bg_f_ + i]      = f;
    pvecback[index_bg_dfdlnq_ + i] = pvecback_B[index_bi_dfdlnq_ + i];
  }
  double n_wdm, rho_wdm, p_wdm, pseudo_p_wdm;
  ComputeMomenta(z, &n_wdm, &rho_wdm, &p_wdm, nullptr, &pseudo_p_wdm);
  pvecback[index_bg_number_]   = n_wdm;
  pvecback[index_bg_rho_]      = rho_wdm;
  pvecback[index_bg_p_]        = p_wdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_wdm;
  // index_bg_inj_ columns are written by DCDM_WDM_Species::ComputeBackground
  // after the child loop (they need rho_dcdm).
}

double WdmDecayProductSpecies::PPrime(double a,
                                      double H,
                                      const double* /*pvecback_B*/,
                                      const double* pvecback) const {
  // Standard NCDM redshift part + injection contribution to dP/dtau:
  //   dP_inj = factor/(3 a^4) * sum dq q^4/eps * J
  double dp_inj = 0.;
  for (int i = 0; i < n_bins_; ++i) {
    const double q        = q_[i];
    const double q2       = q * q;
    const double epsilon  = std::sqrt(q2 + a * a * M_ * M_);
    dp_inj               += dq_[i] * q2 * q2 / epsilon * pvecback[index_bg_inj_ + i];
  }
  const double a2 = a * a;
  return a * H * (pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_]) +
         factor_ / (3. * a2 * a2) * dp_inj;
}

double WdmDecayProductSpecies::RhoDotOverRho(const double* pvecback, double a_prime_over_a) const {
  const double rho = Rho(pvecback);
  if (rho <= 0.)
    return 0.;
  const double a   = pvecback[bgm_->index_bg_a_];
  const double inj = (parent_ != nullptr) ? a * Gamma_ * parent_->Rho(pvecback) : 0.;
  return (-3. * a_prime_over_a * (rho + P(pvecback)) + inj) / rho;
}

// ── Background output columns ────────────────────────────────────────────────

void WdmDecayProductSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  const std::string& nm = name();
  w.Add("(.)rho_wdm_" + nm, 0.);
  w.Add("(.)p_wdm_" + nm, 0.);
  for (int i = 0; i < n_bins_; ++i)
    w.Add("f_" + nm + "[" + std::to_string(i) + "]", 0.);
}

void WdmDecayProductSpecies::WriteBackgroundData(const double* pvecback,
                                                 BackgroundColumnWriter& w) const {
  const std::string& nm = name();
  w.Add("(.)rho_wdm_" + nm, Rho(pvecback));
  w.Add("(.)p_wdm_" + nm, P(pvecback));
  for (int i = 0; i < n_bins_; ++i)
    w.Add("f_" + nm + "[" + std::to_string(i) + "]", pvecback[index_bg_f_ + i]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbations
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                         perturb_vector* /*pv*/,
                                                         const precision* ppr,
                                                         int& index_pt,
                                                         const perturb_workspace* /*ppw*/,
                                                         int /*gauge*/) {
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);

  // No fluid approximation: always the full Boltzmann hierarchy (DNCDM precedent).
  layout.l_max  = (l_max_input_ > 0) ? l_max_input_ : ppr->l_max_ncdm;
  layout.q_size = q_size();

  layout.index_per_q.clear();
  layout.index_per_q.reserve(layout.q_size);
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  index_pt += layout.total_size();
}

void WdmDecayProductSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                           double /*tau*/,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback.data();
  const double M_wdm     = M_;
  const int lmax         = layout.l_max;

  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q = q_[iq];
    // Unnormalized slope g = df0/dlnq: multiplies exactly the terms that carry
    // dlnf0/dlnq in the normalized (Psi) hierarchy.
    const double g = pvecback[index_bg_dfdlnq_ + iq];

    const double epsilon        = std::sqrt(q * q + a2 * M_wdm * M_wdm);
    const double qk_div_epsilon = k * q / epsilon;
    const int idx               = layout.index_per_q[iq];

    dy[idx]     = -qk_div_epsilon * y[idx + 1] + metric_continuity * g / 3.;
    dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) -
                  epsilon * metric_euler / (3. * q * k) * g;
    dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                  s_l[2] * metric_shear * 2. / 15. * g;
    for (int l = 3; l < lmax; ++l)
      dy[idx + l] = qk_div_epsilon / (2. * l + 1.) *
                    (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
    dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] - (1. + lmax) * k * cotKgen * y[idx + lmax];
  }
}

void WdmDecayProductSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;  // isocurvature: daughter starts empty and unperturbed

  // Bins already populated at tau_ini (very short lifetimes) carry the parent's
  // adiabatic density contrast; superhorizon injection shares delta_dcdm.
  const double* pvecback     = ctx.ppw->pvecback.data();
  const double delta_dcdm_ic = 0.75 * ctx.delta_g_ic;
  for (int iq = 0; iq < layout.q_size; ++iq)
    y[layout.index_per_q[iq]] = pvecback[index_bg_f_ + iq] * delta_dcdm_ic;
  // psi_{l>=1} start at zero.
}

BaseSpecies::StressEnergyContribution WdmDecayProductSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_p_];
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return se;

  const double k  = ppw->scalar_ctx.k;
  const double a2 = ppw->scalar_ctx.a2;

  double s_drho = 0., s_theta = 0., s_dp = 0., s_shear = 0.;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q        = q_[iq];
    const double q2       = q * q;
    const double epsilon  = std::sqrt(q2 + a2 * M_ * M_);
    const int idx         = layout.index_per_q[iq];
    const double dqv      = dq_[iq];
    s_drho               += dqv * q2 * epsilon * y[idx];
    s_theta              += dqv * q2 * q * y[idx + 1];
    s_dp                 += dqv * q2 * q2 / epsilon * y[idx];
    s_shear              += dqv * q2 * q2 / epsilon * y[idx + 2];
  }
  const double fac    = factor_ / (a2 * a2);
  se.delta_rho        = fac * s_drho;
  se.rho_plus_p_theta = fac * k * s_theta;
  se.delta_p          = fac * s_dp / 3.;
  se.rho_plus_p_shear = fac * (2. / 3.) * s_shear;
  return se;
}

void WdmDecayProductSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                         const double* y,
                                         const double* /*dy*/,
                                         PerturbSourceContext& ctx) const {
  if (ctx.index_md != ctx.p_mod->index_md_scalars_)
    return;
  if (index_tp_delta_ < 0 && index_tp_theta_ < 0)
    return;

  const auto& layout     = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  const double* pvecback = ctx.ppw->pvecback.data();
  const double rho       = Rho(pvecback);
  const double p         = P(pvecback);
  const double a2        = ctx.a2;
  const double k         = ctx.k;

  double s_drho = 0., s_theta = 0.;
  if (layout.q_size > 0 && !layout.index_per_q.empty()) {
    for (int iq = 0; iq < layout.q_size; ++iq) {
      const double q        = q_[iq];
      const double q2       = q * q;
      const double epsilon  = std::sqrt(q2 + a2 * M_ * M_);
      const int idx         = layout.index_per_q[iq];
      s_drho               += dq_[iq] * q2 * epsilon * y[idx];
      s_theta              += dq_[iq] * q2 * q * y[idx + 1];
    }
  }
  const double fac = factor_ / (a2 * a2);

  if (index_tp_delta_ >= 0) {
    // delta with N-body gauge correction: delta_Nb = delta - (rhodot/rho) theta_tot/k^2
    const double delta = (rho > 0.) ? fac * s_drho / rho : 0.;
    const double src = (rho > 0.)
                           ? delta - RhoDotOverRho(pvecback, ctx.a_prime_over_a) * ctx.theta_over_k2
                           : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md,
                              ctx.index_ic,
                              index_tp_delta_,
                              ctx.index_tau,
                              ctx.index_k,
                              src);
  }
  if (index_tp_theta_ >= 0) {
    const double rho_plus_p = rho + p;
    const double theta      = (rho_plus_p > 0.) ? fac * k * s_theta / rho_plus_p : 0.;
    const double src        = (rho_plus_p > 0.) ? theta + ctx.theta_shift : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md,
                              ctx.index_ic,
                              index_tp_theta_,
                              ctx.index_tau,
                              ctx.index_k,
                              src);
  }
}

void WdmDecayProductSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                                const PerturbationsModule& mod,
                                                file_format fmt,
                                                BaseSpecies::TransferColumnSection section) const {
  if (fmt != file_format::class_format)
    return;
  const perturbs* ppt = mod.GetPerturbs();
  if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
    w.Add("d_wdm_" + name(), index_tp_delta_, index_tp_delta_ >= 0);
  if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
    w.Add("t_wdm_" + name(), index_tp_theta_, index_tp_theta_ >= 0);
}
