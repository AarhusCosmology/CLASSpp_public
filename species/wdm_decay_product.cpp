#include "wdm_decay_product.h"

#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "dcdm.h"
#include "perturbations_module.h"
#include "species/species_input.h"

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
    if (input.get<std::string>(key).has_value()) {
      throw std::invalid_argument("species '" + instance_name + "': parameter '" +
                                  std::string(key) + "' is not supported by dcdm_wdm");
    }
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
  if (n_provided != 1) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': specify exactly one of Gamma, log10Gamma, lifetime, log10lifetime");
  }
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
  if (eps_opt.has_value() == vkick_opt.has_value()) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of epsilon, vkick");
  }
  if (eps_opt) {
    eps_ = *eps_opt;
    if (eps_ < 0. || eps_ > 1. - 1e-9) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': epsilon must be in [0, 1 - 1e-9]; for smaller kicks "
                                  "input vkick directly");
    }
    vkick_ = std::sqrt((1. - eps_) * (1. + eps_));
  }
  else {
    vkick_ = *vkick_opt;
    if (vkick_ < 1e-8 || vkick_ > 1.) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': vkick must be in [1e-8, 1] (in units of c); for "
                                  "vkick -> 0 the model is plain CDM");
    }
    eps_ = std::sqrt((1. - vkick_) * (1. + vkick_));
  }
  // p/m_d = v/sqrt(1-v^2) = v/eps  =>  M = m_d/T0 = kQKick * eps / v.
  M_ = kQKick * eps_ / vkick_;

  // Momentum grid: uniform in u = ln q on [ln q_lo, ln q_kick + 3 sigma], so the
  // kernel keeps a 3 sigma margin above q_cut(a=1) = kQKick for any bin count:
  //   du (n_bins - 3 kernel_width) = ln(kQKick / q_lo).
  n_bins_       = input.get_or("momenta_bins", 96);
  q_min_ratio_  = input.get_or("q_min_ratio", 1e-4);
  kernel_width_ = input.get_or("kernel_width", 1.0);
  l_max_input_  = input.get_or("l_max", -1);
  if (n_bins_ < 8 || n_bins_ > 4096) {
    throw std::invalid_argument("species '" + instance_name + "': momenta_bins out of range");
  }
  if (!(q_min_ratio_ > 0. && q_min_ratio_ < 0.5)) {
    throw std::invalid_argument("species '" + instance_name + "': q_min_ratio must be in (0, 0.5)");
  }
  if (!(kernel_width_ >= 0.2 && kernel_width_ <= 5.)) {
    throw std::invalid_argument("species '" + instance_name +
                                "': kernel_width must be in [0.2, 5] bins");
  }
  if (n_bins_ <= 3. * kernel_width_ + 4.) {
    throw std::invalid_argument("species '" + instance_name +
                                "': momenta_bins must exceed 3*kernel_width + 4 "
                                "(the grid must resolve the injection kernel)");
  }

  const double q_lo = kQKick * q_min_ratio_;
  du_               = std::log(kQKick / q_lo) / (n_bins_ - 3. * kernel_width_);
  q_.resize(n_bins_);
  u_.resize(n_bins_);
  dq_.resize(n_bins_);
  for (int i = 0; i < n_bins_; ++i) {
    u_[i]  = std::log(q_lo) + (i + 0.5) * du_;
    q_[i]  = std::exp(u_[i]);
    dq_[i] = q_[i] * du_;
  }
  q_bg_ = q_;
  w_.assign(n_bins_, 0.);
  w_bg_.assign(n_bins_, 0.);        // refreshed from f each ComputeBackground
  dlnf0_dlnq_.assign(n_bins_, 0.);  // unused (tensor slots disabled)
  scratch_G_.assign(n_bins_, 0.);

  // Normalization: Omega_ini/omega_ini (parent initial abundance, dcdm convention)
  // XOR Omega_dcdmwdm/omega_dcdmwdm (combined sector today). During a shooting
  // build DoShooting legitimately writes <instance>.Omega_dcdmwdm alongside the
  // user's initial key.
  auto Omega_ini_opt  = input.get<double>("Omega_ini");
  auto omega_ini_opt  = input.get<double>("omega_ini");
  auto Omega_comb_opt = input.get<double>("Omega_dcdmwdm");
  auto omega_comb_opt = input.get<double>("omega_dcdmwdm");
  if (Omega_ini_opt && omega_ini_opt) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify only one of Omega_ini, omega_ini");
  }
  if (Omega_comb_opt && omega_comb_opt) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify only one of Omega_dcdmwdm, omega_dcdmwdm");
  }
  const bool has_initial  = Omega_ini_opt.has_value() || omega_ini_opt.has_value();
  const bool has_combined = Omega_comb_opt.has_value() || omega_comb_opt.has_value();
  const bool in_shooting  = (pfc != nullptr) && pfc->is_shooting;
  if (has_initial && has_combined && !in_shooting) {
    throw std::invalid_argument("species '" + instance_name +
                                "': Omega_ini conflicts with Omega_dcdmwdm — choose one");
  }
  if (!has_initial && !has_combined) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify Omega_ini/omega_ini or Omega_dcdmwdm/omega_dcdmwdm");
  }
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
// FillInjection
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::FillInjection(double a,
                                           double rho_dcdm,
                                           double* J,
                                           double* dJdlnq) const {
  const int N        = n_bins_;
  const double q_cut = a * kQKick;
  const double u_cut = std::log(q_cut);
  const double sigma = kernel_width_ * du_;

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

  // Smooth onset gate: injection ramps on as a Gaussian over the first few
  // kernel widths of the grid instead of switching on discontinuously when
  // q_cut crosses the grid edge (a hard switch-on makes the stiff evolver
  // grind to a halt at rtol ~ 1e-6). Decays suppressed or dropped by the gate
  // are the pre-grid ones, negligible by construction (~1e-6 for defaults).
  // The clamp value 60 keeps every residual on/off step at exp(-60) ~ 9e-27 —
  // orders of magnitude below the evolver's error-control threshold — while
  // staying far away from exp-underflow/denormal territory, which is not
  // IEEE-safe under -ffast-math (auto-vectorized exp, flush-to-zero).
  double gate       = 1.;
  const double x_on = (u_cut - (u_[0] + 3. * sigma)) / sigma;
  if (x_on < 0.) {
    const double x2_on = 0.5 * x_on * x_on;
    if (x2_on >= 60.) {
      zero_all();
      return;
    }
    gate = std::exp(-x2_on);
  }

  double D = 0.;
  for (int i = 0; i < N; ++i) {
    const double x  = (u_[i] - u_cut) / sigma;
    const double x2 = 0.5 * x * x;
    // Clamp: exp of very large negative arguments must never be evaluated —
    // under -ffast-math the auto-vectorized exp returns garbage outside its
    // valid range, and denormal results get flushed to zero mid-expression.
    // exp(-60) ~ 9e-27 is far below any relevant weight.
    const double G        = (x2 >= 60.) ? 0. : std::exp(-x2);
    scratch_G_[i]         = G;
    const double epsilon  = std::sqrt(q_[i] * q_[i] + a * a * M_ * M_);
    D                    += dq_[i] * q_[i] * q_[i] * epsilon * G;
  }
  if (D <= 0.) {
    zero_all();
    return;
  }
  const double a2 = a * a;
  const double A  = gate * a * Gamma_ * rho_dcdm * a2 * a2 / factor_;
  for (int i = 0; i < N; ++i) {
    J[i] = A * scratch_G_[i] / D;
    if (dJdlnq != nullptr)
      dJdlnq[i] = J[i] * (u_cut - u_[i]) / (sigma * sigma);
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
