#include "ncdm_species.h"

#include <cmath>

#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations_module.h"
#include "species/species_input.h"

// ── Constructor ─────────────────────────────────────────────────────────────

NCDMSpecies::NCDMSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings), pba_(pba) {
  bgm_ = bgm;

  // Standard-NCDM closure: NCDMBaseSpecies::ReadParametersByInstance only
  // computes M_ from m_in_eV_. Standard NCDM additionally needs the closure
  // tying Omega0_ to m / factor / deg / M_:
  //   - if m != 0 and Omega0 == 0: derive Omega0_ from rho_ncdm
  //   - if m != 0 and Omega0 != 0: rescale factor_ and deg_ by fnu_factor
  //   - if m == 0:                 derive M_ via MFromOmega and back-compute m_in_eV_
  // DNCDMSpecies skips this and instead uses a deferred Omega_ini closure
  // applied by DNCDMSpecies::CreateAll (see ApplyDncdmInitialClosure).
  const double H0 = settings.h * 1.e5 / _c_;
  if (m_in_eV_ != 0.0) {
    double rho_ncdm = 0.;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
    if (GetOmega0() == 0.0) {
      SetOmega0(rho_ncdm / H0 / H0, settings.h);
    }
    else {
      const double fnu_factor  = H0 * H0 * GetOmega0() / rho_ncdm;
      factor_                 *= fnu_factor;
      deg_                    *= fnu_factor;
    }
  }
  else {
    M_       = MFromOmega(H0, GetOmega0(), settings.tol_M_ncdm);
    m_in_eV_ = _k_B_ / _eV_ * T_ * M_ * T_cmb_;
  }

  // ncdm_id_ stays at its in-class default (-1) until SetNcdmId is called.
}

// ── CreateAll factory ───────────────────────────────────────────────────────

namespace {

struct LegacyFieldMap {
  const char* legacy;  // legacy CSV key, e.g. "m_ncdm"
  const char* dot;     // per-instance dot field, e.g. "m"
};

// Note: use_ncdm_psd_files MUST precede ncdm_psd_filenames so the PSD-filename
// branch in TransmuteLegacyStandardNcdmToDot can read the just-written
// use_psd_file dot keys for each instance.
constexpr LegacyFieldMap kLegacyStandardFields[] = {
    {"m_ncdm", "m"},
    {"T_ncdm", "T"},
    {"deg_ncdm", "deg"},
    {"Omega_ncdm", "Omega"},
    {"omega_ncdm", "omega"},
    {"ksi_ncdm", "ksi"},
    {"ncdm_psd_parameters", "psd_parameters"},
    {"use_ncdm_psd_files", "use_psd_file"},
    {"ncdm_psd_filenames", "psd_filename"},
    {"quadrature_strategy_ncdm_standard", "quadrature_strategy"},
    {"N_momentum_bins_ncdm_standard", "momenta_bins"},
    {"maximum_q_ncdm_standard", "max_q"},
};

// Returns true if any key of the form "<prefix>.*" exists in pfc.
bool any_dot_key_with_prefix(const FileContent& pfc, const std::string& prefix) {
  bool found = false;
  pfc.for_each([&](const std::string& name, const std::string&, bool) {
    if (name.rfind(prefix + ".", 0) == 0) {
      found = true;
    }
  });
  return found;
}

void TransmuteLegacyStandardNcdmToDot(FileContent& pfc) {
  // 1. Resolve N from N_ncdm / N_ncdm_standard.
  int N1 = 0, N2 = 0;
  bool has1 = pfc.read_int("N_ncdm_standard", N1);
  bool has2 = pfc.read_int("N_ncdm", N2);

  if (has1 && has2) {
    throw std::invalid_argument(
        "In input file, you can only enter one of N_ncdm_standard and N_ncdm, choose one");
  }
  int N = has1 ? N1 : (has2 ? N2 : 0);
  if (N <= 0) {
    return;  // nothing to transmute
  }

  // 2. Synthesised prefix collision check.
  for (int i = 1; i <= N; ++i) {
    const std::string prefix = "ncdm__" + std::to_string(i);
    if (any_dot_key_with_prefix(pfc, prefix)) {
      throw std::invalid_argument(
          "legacy NCDM transmutation collides with existing dot-instance '" + prefix +
          "': rename your dot instance or remove the legacy N_ncdm keys");
    }
  }

  // 3. For each legacy CSV field, broadcast or distribute per instance.
  for (const auto& fm : kLegacyStandardFields) {
    std::string raw_value;
    if (!pfc.read_string(fm.legacy, raw_value)) {
      continue;
    }

    std::vector<std::string> values = FileContent::split_csv(raw_value);
    if (values.empty()) {
      continue;
    }

    // psd_filename special case: legacy lists one filename per instance whose
    // use_ncdm_psd_files flag is _TRUE_, in order. Skip the broadcast/distribute
    // path and assign filenames only to instances that opted in. Counts must
    // match exactly: under- or over-supplying filenames is an error.
    if (std::string(fm.dot) == "psd_filename") {
      std::vector<int> use_flags(N, 0);
      int enabled = 0;
      for (int i = 1; i <= N; ++i) {
        int flag = 0;
        pfc.read_int("ncdm__" + std::to_string(i) + ".use_psd_file", flag);
        use_flags[i - 1] = flag;
        if (flag != 0) {
          ++enabled;
        }
      }
      if (static_cast<int>(values.size()) != enabled) {
        throw std::invalid_argument("ncdm_psd_filenames has " + std::to_string(values.size()) +
                                    " entries but " + std::to_string(enabled) +
                                    " instances have use_ncdm_psd_files=1");
      }
      size_t fn_idx = 0;
      for (int i = 0; i < N; ++i) {
        if (use_flags[i] != 0) {
          pfc.set("ncdm__" + std::to_string(i + 1) + ".psd_filename", values[fn_idx]);
          ++fn_idx;
        }
      }
      continue;  // skip the generic broadcast/distribute path
    }

    if (values.size() == 1u) {
      // broadcast scalar to all N instances
      for (int i = 1; i <= N; ++i) {
        pfc.set("ncdm__" + std::to_string(i) + "." + fm.dot, values[0]);
      }
    }
    else if (static_cast<int>(values.size()) == N) {
      for (int i = 0; i < N; ++i) {
        pfc.set("ncdm__" + std::to_string(i + 1) + "." + fm.dot, values[i]);
      }
    }
    else {
      throw std::invalid_argument(std::string(fm.legacy) + " has " + std::to_string(values.size()) +
                                  " entries but N_ncdm_standard = " + std::to_string(N) +
                                  "; expected 1 (broadcast) or " + std::to_string(N));
    }
  }

  // 4. Stamp the type field on each synthesised instance.
  for (int i = 1; i <= N; ++i) {
    pfc.set("ncdm__" + std::to_string(i) + ".type", "ncdm_standard");
  }
}

}  // namespace

std::vector<Named> NCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  TransmuteLegacyStandardNcdmToDot(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_standard");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    // Mark <instance>.type as consumed (instances_with does not mark anything as read)
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);

    auto sp = std::make_unique<NCDMSpecies>(ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    result.push_back({name, std::move(sp)});
  }
  return result;
}

// ── Background ─────────────────────────────────────────────────────────────

void NCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_   = index_bg++;
  index_bg_rho_      = index_bg++;  // base class protected
  index_bg_p_        = index_bg++;  // base class protected
  index_bg_pseudo_p_ = index_bg++;
}

void NCDMSpecies::RegisterIntegrationIndices(int& /*index_bi*/) {
  // Stable NCDM has no integration variables
}

void NCDMSpecies::ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) {
  double z = 1. / a_rel - 1.;

  double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
  ComputeMomenta(z, &number_ncdm, &rho_ncdm, &p_ncdm, nullptr, &pseudo_p_ncdm);
  pvecback[index_bg_number_]   = number_ncdm;
  pvecback[index_bg_rho_]      = rho_ncdm;
  pvecback[index_bg_p_]        = p_ncdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_ncdm;
}

void NCDMSpecies::BackgroundDerivs(double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const double* /*pvecback*/) {
  // Stable NCDM has no background derivatives
}

void NCDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
}

void NCDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, pvecback[bg_number_index()]);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, Rho(pvecback));
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, P(pvecback));
}

void NCDMSpecies::FillSources(const double* /*y*/,
                              const double* /*dy*/,
                              PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;

  // NCDM sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const int n              = ncdm_id_;
  const double* pvecback   = ppw->pvecback;
  const perturb_vector* pv = ppw->pv;
  const double* y          = ppw->pv->y;

  // delta_ncdm[n]: density perturbation
  if (p_mod->has_source_delta_ncdm_ == _TRUE_) {
    const double w = pvecback[index_bg_p_] / pvecback[index_bg_rho_];
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_ncdm1_ + n,
                          ctx.index_tau,
                          ctx.index_k,
                          Delta(pv, y, pvecback, ppw) +
                              3. * ctx.a_prime_over_a * (1. + w) *
                                  ctx.theta_over_k2);  // N-body gauge correction
  }

  // theta_ncdm[n]: velocity perturbation
  if (p_mod->has_source_theta_ncdm_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_ncdm1_ + n,
                          ctx.index_tau,
                          ctx.index_k,
                          Theta(pv, y, pvecback, ppw) + ctx.theta_shift);
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void NCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                              const precision* ppr,
                                              int& index_pt,
                                              const perturb_workspace* ppw,
                                              int /*gauge*/) {
  if (!pba_->has_ncdm)
    return;

  if (ncdm_id_ == 0) {
    pv->index_pt_psi0_ncdm1 = index_pt;
  }
  index_pt_psi0_ = index_pt;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);

  if (fa_on) {
    pv->l_max_ncdm[ncdm_id_]  = 2;
    pv->q_size_ncdm[ncdm_id_] = 1;
  }
  else {
    pv->l_max_ncdm[ncdm_id_]  = ppr->l_max_ncdm;
    pv->q_size_ncdm[ncdm_id_] = q_size();
  }
  pv->index_ncdm_[ncdm_id_].clear();
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq)
    pv->index_ncdm_[ncdm_id_].push_back(index_pt + iq * (pv->l_max_ncdm[ncdm_id_] + 1));
  index_pt += (pv->l_max_ncdm[ncdm_id_] + 1) * pv->q_size_ncdm[ncdm_id_];
}

void NCDMSpecies::PerturbDerivs(double tau,
                                const double* y,
                                double* dy,
                                const perturb_parameters_and_workspace& ppaw) {
  if (!pba_->has_ncdm)
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double a_prime_over_a     = ctx.a_prime_over_a;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double metric_ufa_class   = ctx.metric_ufa_class;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);

  if (fa_on) {
    // Fluid approximation (ncdmfa_on)
    const double rho_ncdm        = pvecback[index_bg_rho_];
    const double p_ncdm          = pvecback[index_bg_p_];
    const double pseudo_p        = pvecback[index_bg_pseudo_p_];
    const double w_ncdm          = p_ncdm / rho_ncdm;
    const double ca2_ncdm        = w_ncdm / 3. / (1. + w_ncdm) * (5. - pseudo_p / p_ncdm);
    const double ceff2           = ca2_ncdm;
    const double cvis2           = 3. * w_ncdm * ca2_ncdm;
    const double pseudo_p_over_p = pseudo_p / p_ncdm;
    const int idx                = pv->index_ncdm_.at(ncdm_id_)[0];

    dy[idx] = -(1. + w_ncdm) * (y[idx + 1] + metric_continuity) -
              3. * a_prime_over_a * (ceff2 - w_ncdm) * y[idx];

    dy[idx + 1] = -a_prime_over_a * (1. - 3. * ca2_ncdm) * y[idx + 1] +
                  ceff2 / (1. + w_ncdm) * k * k * y[idx] - k * k * y[idx + 2] + metric_euler;

    // CLASS fluid approximation for shear
    dy[idx + 2] = -3. * (a_prime_over_a * (2. / 3. - ca2_ncdm - pseudo_p_over_p / 3.) + 1. / tau) *
                      y[idx + 2] +
                  8. / 3. * cvis2 / (1. + w_ncdm) * s_l[2] * (y[idx + 1] + metric_ufa_class);
  }
  else {
    // Exact Boltzmann hierarchy per momentum bin
    const double M_ncdm = M_;
    for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
      const double q          = q_[iq];
      const double dlnf0_dlnq = dlnf0_dlnq_[iq];

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
      dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] -
                       (1. + lmax) * k * cotKgen * y[idx + lmax];
    }
  }
}

void NCDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  if (!pba_->has_ncdm || pv->index_ncdm_.at(ncdm_id_).empty())
    return;
  for (int index_q = 0; index_q < pv->q_size_ncdm[ncdm_id_]; ++index_q) {
    const int idx           = pv->index_ncdm_[ncdm_id_][index_q];
    const double q          = q_[index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);
    const double dlnf0_dlnq = dlnf0_dlnq_[index_q];
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

// ── Output columns ─────────────────────────────────────────────────────────

void NCDMSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                     const PerturbationsModule& mod,
                                     enum file_format fmt,
                                     BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (pba->has_ncdm != _TRUE_)
    return;

  const int n = ncdm_id_;
  char tmp[40];

  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_) {
      snprintf(tmp, sizeof(tmp), "d_ncdm[%d]", n);
      w.Add(tmp, mod.index_tp_delta_ncdm1_ + n, _TRUE_);
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_) {
      snprintf(tmp, sizeof(tmp), "t_ncdm[%d]", n);
      w.Add(tmp, mod.index_tp_theta_ncdm1_ + n, _TRUE_);
    }
  }
  else if (fmt == camb_format) {
    // camb_format: single aggregate column emitted only for n==0
    if (section != TransferColumnSection::velocity && n == 0)
      w.Add("-T_ncdm/k2", mod.index_tp_delta_ncdm1_, _TRUE_);
  }
}

void NCDMSpecies::PrintVariables(PerturbColumnWriter& w,
                                 double /*tau*/,
                                 const double* y,
                                 const PerturbationsModule& mod,
                                 const perturb_workspace* ppw) const {
  const background* pba = mod.GetBackground();
  if (pba->has_ncdm != _TRUE_)
    return;

  const int n       = ncdm_id_;
  double delta_ncdm = 0., theta_ncdm = 0., shear_ncdm = 0., cs2_ncdm = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;
    const double k           = ppw->scalar_ctx.k;
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];

    const double rho_ncdm_bg = Rho(pvecback);
    const double p_ncdm_bg   = P(pvecback);
    const double rho_plus_p  = rho_ncdm_bg + p_ncdm_bg;
    const double w_ncdm      = (rho_ncdm_bg > 0.) ? p_ncdm_bg / rho_ncdm_bg : 0.;

    delta_ncdm = Delta(pv, y, pvecback, ppw);
    theta_ncdm = Theta(pv, y, pvecback, ppw);
    shear_ncdm = (rho_plus_p > 0.) ? RhoPlusPShear(pv, y, pvecback, ppw) / rho_plus_p : 0.;
    const double delta_p_ncdm   = DeltaP(pv, y, pvecback, ppw);
    const double delta_rho_ncdm = rho_ncdm_bg * delta_ncdm;
    constexpr double eps        = 1e-300;
    cs2_ncdm = (std::abs(delta_rho_ncdm) > eps) ? delta_p_ncdm / delta_rho_ncdm : 0.;

    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      const double alpha_corr = pvecmetric[ppw->index_mt_alpha];
      // Store pre-correction delta for cs2 update
      const double delta_ncdm_syn  = delta_ncdm;
      delta_ncdm                  -= 3. * a * H * (1. + w_ncdm) * alpha_corr;
      theta_ncdm                  += k * k * alpha_corr;
      // Update cs2_ncdm with gauge correction
      if (p_ncdm_bg > eps) {
        const double pseudo_p_ncdm         = pvecback[index_bg_pseudo_p_];
        const double p_prime_over_rho_ncdm = -a * H * w_ncdm * (5. - pseudo_p_ncdm / p_ncdm_bg);
        const double delta_ncdm_corrected  = delta_ncdm_syn -
                                             3. * (1. + w_ncdm) * a * H * alpha_corr;
        if (std::abs(delta_ncdm_corrected) > eps)
          cs2_ncdm = (cs2_ncdm * delta_ncdm_syn + p_prime_over_rho_ncdm * alpha_corr) /
                     delta_ncdm_corrected;
      }
    }
  }

  char tmp[40];
  snprintf(tmp, sizeof(tmp), "delta_ncdm[%d]", n);
  w.Add(tmp, delta_ncdm, true);
  snprintf(tmp, sizeof(tmp), "theta_ncdm[%d]", n);
  w.Add(tmp, theta_ncdm, true);
  snprintf(tmp, sizeof(tmp), "shear_ncdm[%d]", n);
  w.Add(tmp, shear_ncdm, true);
  snprintf(tmp, sizeof(tmp), "cs2_ncdm[%d]", n);
  w.Add(tmp, cs2_ncdm, true);
}

// ── Integrated observables ──────────────────────────────────────────────────

double NCDMSpecies::Delta(const perturb_vector* pv,
                          const double* y,
                          const double* pvecback,
                          const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);
  if (fa_on)
    return y[pv->index_ncdm_.at(ncdm_id_)[0]];

  const double a        = ppw->scalar_ctx.a;
  double rho_delta_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    const double w0       = w_[iq];
    const double q        = q_[iq];
    const double epsilon  = std::sqrt(q * q + std::pow(M_ * a, 2));
    rho_delta_ncdm       += q * q * epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double factor = factor_ * std::pow(pba_->a_today / a, 4);
  return rho_delta_ncdm * factor / pvecback[index_bg_rho_];
}

double NCDMSpecies::Theta(const perturb_vector* pv,
                          const double* y,
                          const double* pvecback,
                          const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);
  if (fa_on)
    return y[pv->index_ncdm_.at(ncdm_id_)[0] + 1];

  const double a               = ppw->scalar_ctx.a;
  double rho_plus_p_theta_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    const double w0        = w_[iq];
    const double q         = q_[iq];
    rho_plus_p_theta_ncdm += q * q * q * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq] + 1];
  }
  const double factor = factor_ * std::pow(pba_->a_today / a, 4);
  const double k      = ppw->scalar_ctx.k;
  return rho_plus_p_theta_ncdm * k * factor / (pvecback[index_bg_rho_] + pvecback[index_bg_p_]);
}

double NCDMSpecies::DeltaP(const perturb_vector* pv,
                           const double* y,
                           const double* pvecback,
                           const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);
  if (fa_on) {
    // cg2_ncdm logic for fluid approx
    double rho_bg      = pvecback[index_bg_rho_];
    double p_bg        = pvecback[index_bg_p_];
    double pseudo_p_bg = pvecback[index_bg_pseudo_p_];
    double w_ncdm      = p_bg / rho_bg;
    double cg2_ncdm    = w_ncdm * (1.0 - 1.0 / (3.0 + 3.0 * w_ncdm) *
                                             (3.0 * w_ncdm - 2.0 + pseudo_p_bg / p_bg));
    return cg2_ncdm * rho_bg * y[pv->index_ncdm_.at(ncdm_id_)[0]];
  }

  const double a      = ppw->scalar_ctx.a;
  double delta_p_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    const double w0       = w_[iq];
    const double q        = q_[iq];
    const double epsilon  = std::sqrt(q * q + std::pow(M_ * a, 2));
    delta_p_ncdm         += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double factor = factor_ * std::pow(pba_->a_today / a, 4);
  return delta_p_ncdm * factor / 3.;
}

double NCDMSpecies::RhoPlusPShear(const perturb_vector* pv,
                                  const double* y,
                                  const double* pvecback,
                                  const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);
  if (fa_on)
    return (pvecback[index_bg_rho_] + pvecback[index_bg_p_]) *
           y[pv->index_ncdm_.at(ncdm_id_)[0] + 2];

  const double a               = ppw->scalar_ctx.a;
  double rho_plus_p_shear_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    const double w0        = w_[iq];
    const double q         = q_[iq];
    const double epsilon   = std::sqrt(q * q + std::pow(M_ * a, 2));
    rho_plus_p_shear_ncdm += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq] + 2];
  }
  const double factor = factor_ * std::pow(pba_->a_today / a, 4);
  return 2.0 / 3.0 * factor * rho_plus_p_shear_ncdm;
}
