#include "ncdm_species.h"

#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "perturbations_module.h"
#include "species/ncdm_family.h"
#include "species/species_input.h"

// ── Constructors ────────────────────────────────────────────────────────────

NCDMSpecies::NCDMSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings), pba_(pba) {
  bgm_ = bgm;
  ResolveMassOmegaClosure(settings);
}

NCDMSpecies::NCDMSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm,
                         NCDMBaseSpecies::DeferInit defer)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings, defer),
      pba_(pba) {
  bgm_ = bgm;
  // No closure here: the subclass calls BuildQuadratureAndMass + ResolveMassOmegaClosure
  // after configuring its PSD.
}

void NCDMSpecies::ResolveMassOmegaClosure(const NcdmSettings& settings) {
  // At entry, NCDMBaseSpecies::BuildQuadratureAndMass has already set M_ from
  // m_in_eV_; this step resolves Omega0_ (or rescales deg_/factor_) so the
  // background density matches what was requested.
  // Standard-NCDM closure tying Omega0_ to m / factor / deg / M_:
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
    {"use_ncdm_psd_files", "use_psd_file"},
    {"ncdm_psd_filenames", "psd_filename"},
    {"quadrature_strategy_ncdm_standard", "quadrature_strategy"},
    {"N_momentum_bins_ncdm_standard", "momenta_bins"},
    {"N_momentum_bins_bg_ncdm_standard", "momenta_bins_bg"},
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
  auto n1_opt = pfc.get<int>("N_ncdm_standard");
  auto n2_opt = pfc.get<int>("N_ncdm");
  bool has1   = n1_opt.has_value();
  bool has2   = n2_opt.has_value();
  if (has1)
    N1 = *n1_opt;
  if (has2)
    N2 = *n2_opt;

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
    auto raw_opt = pfc.get<std::string>(fm.legacy);
    if (!raw_opt) {
      continue;
    }
    const std::string& raw_value = *raw_opt;

    std::vector<std::string> values = FileContent::split_csv(raw_value);
    if (values.empty()) {
      continue;
    }

    // psd_filename special case: legacy lists one filename per instance whose
    // use_ncdm_psd_files flag is true, in order. Skip the broadcast/distribute
    // path and assign filenames only to instances that opted in. Counts must
    // match exactly: under- or over-supplying filenames is an error.
    if (std::string(fm.dot) == "psd_filename") {
      std::vector<int> use_flags(N, 0);
      int enabled = 0;
      for (int i = 1; i <= N; ++i) {
        int flag         = 0;
        flag             = pfc.get_or("ncdm__" + std::to_string(i) + ".use_psd_file", flag);
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
  return CreateAllNcdmInstances<NCDMSpecies>(ctx);
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

void NCDMSpecies::ComputeBackground(double a, const double* /*pvecback_B*/, double* pvecback) {
  double z = 1. / a - 1.;

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
  w.Add("(.)number_" + name(), 0.);
  w.Add("(.)rho_" + name(), 0.);
  w.Add("(.)p_" + name(), 0.);
}

void NCDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)number_" + name(), pvecback[bg_number_index()]);
  w.Add("(.)rho_" + name(), Rho(pvecback));
  w.Add("(.)p_" + name(), P(pvecback));
}

void NCDMSpecies::FillSources(const BaseSpecies::PerturbLayout& layout,
                              const double* /*y*/,
                              const double* /*dy*/,
                              PerturbSourceContext& ctx) const {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;

  // NCDM sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const double* pvecback   = ppw->pvecback.data();
  const perturb_vector* pv = ppw->pv.get();
  const double* y          = ppw->pv->y.data();

  // delta_ncdm[n]: density perturbation
  if (index_tp_delta_ >= 0 || index_tp_theta_ >= 0) {
    const auto se = StressEnergy(layout, pv, y, pvecback, ppw);
    if (index_tp_delta_ >= 0) {
      const double w = pvecback[index_bg_p_] / pvecback[index_bg_rho_];
      p_mod->SetSourceValue(ctx.index_md,
                            ctx.index_ic,
                            index_tp_delta_,
                            ctx.index_tau,
                            ctx.index_k,
                            se.delta_rho / pvecback[index_bg_rho_] +
                                3. * ctx.a_prime_over_a * (1. + w) *
                                    ctx.theta_over_k2);  // N-body gauge correction
    }

    // theta_ncdm[n]: velocity perturbation
    if (index_tp_theta_ >= 0) {
      p_mod->SetSourceValue(ctx.index_md,
                            ctx.index_ic,
                            index_tp_theta_,
                            ctx.index_tau,
                            ctx.index_k,
                            se.rho_plus_p_theta /
                                    (pvecback[index_bg_rho_] + pvecback[index_bg_p_]) +
                                ctx.theta_shift);
    }
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void NCDMSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                              perturb_vector* pv,
                                              const precision* ppr,
                                              int& index_pt,
                                              const perturb_workspace* ppw,
                                              int /*gauge*/) {
  // This species only exists when NCDM is present, so the legacy pba->has_ncdm guard is gone.
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);

  index_pt_psi0_ = index_pt;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);
  if (fa_on) {
    layout.l_max  = 2;
    layout.q_size = 1;
  }
  else {
    layout.l_max  = ppr->l_max_ncdm;
    layout.q_size = q_size();
  }
  layout.index_per_q.clear();
  layout.index_per_q.reserve(layout.q_size);
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  index_pt += layout.total_size();
}

void NCDMSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                double tau,
                                const double* y,
                                double* dy,
                                const perturb_parameters_and_workspace& ppaw) const {
  // This species only exists when NCDM is present, so the legacy pba->has_ncdm guard is gone.
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double a_prime_over_a     = ctx.a_prime_over_a;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double metric_ufa_class   = ctx.metric_ufa_class;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback.data();

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
    const int idx                = layout.index_per_q[0];

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
    // Exact Boltzmann hierarchy
    const double M_ncdm = M_;
    const int lmax      = layout.l_max;
    for (int iq = 0; iq < layout.q_size; ++iq) {
      const double q          = q_[iq];
      const double dlnf0_dlnq = dlnf0_dlnq_[iq];

      const double epsilon        = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
      const double qk_div_epsilon = k * q / epsilon;
      const int idx               = layout.index_per_q[iq];

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

void NCDMSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                         double* y,
                                         const PerturbIcContext& ctx) {
  // This species only exists when NCDM is present, so the legacy pba->has_ncdm guard is gone.
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const int lmax = layout.l_max;
  for (int index_q = 0; index_q < layout.q_size; ++index_q) {
    const int idx           = layout.index_per_q[index_q];
    const double q          = q_[index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);
    const double dlnf0_dlnq = dlnf0_dlnq_[index_q];

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
                                     file_format fmt,
                                     BaseSpecies::TransferColumnSection section) const {
  const std::string& nm = name();

  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
      w.Add("d_" + nm, index_tp_delta_, index_tp_delta_ >= 0);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
      w.Add("t_" + nm, index_tp_theta_, index_tp_theta_ >= 0);
  }
  else if (fmt == file_format::camb_format) {
    // file_format::camb_format: per-flavor column (intentional change from aggregate, #309)
    if (section != TransferColumnSection::velocity)
      w.Add("-T_" + nm + "/k2", index_tp_delta_, index_tp_delta_ >= 0);
  }
}

void NCDMSpecies::PrintVariables(PerturbColumnWriter& w,
                                 double /*tau*/,
                                 const double* y,
                                 const PerturbationsModule& mod,
                                 const perturb_workspace* ppw) const {
  double delta_ncdm = 0., theta_ncdm = 0., shear_ncdm = 0., cs2_ncdm = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const double* pvecmetric = ppw->pvecmetric.data();
    const double k           = ppw->scalar_ctx.k;
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];

    const double rho_ncdm_bg = Rho(pvecback);
    const double p_ncdm_bg   = P(pvecback);
    const double rho_plus_p  = rho_ncdm_bg + p_ncdm_bg;
    const double w_ncdm      = (rho_ncdm_bg > 0.) ? p_ncdm_bg / rho_ncdm_bg : 0.;

    const auto& layout          = *pv->species_layouts[collection_index_];
    const auto se               = StressEnergy(layout, pv, y, pvecback, ppw);
    delta_ncdm                  = (rho_ncdm_bg > 0.) ? se.delta_rho / rho_ncdm_bg : 0.;
    theta_ncdm                  = (rho_plus_p > 0.) ? se.rho_plus_p_theta / rho_plus_p : 0.;
    shear_ncdm                  = (rho_plus_p > 0.) ? se.rho_plus_p_shear / rho_plus_p : 0.;
    const double delta_p_ncdm   = se.delta_p;
    const double delta_rho_ncdm = rho_ncdm_bg * delta_ncdm;
    constexpr double eps        = 1e-300;
    cs2_ncdm = (std::abs(delta_rho_ncdm) > eps) ? delta_p_ncdm / delta_rho_ncdm : 0.;

    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == possible_gauges::synchronous) {
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

  const std::string& nm = name();
  w.Add("delta_" + nm, delta_ncdm, true);
  w.Add("theta_" + nm, theta_ncdm, true);
  w.Add("shear_" + nm, shear_ncdm, true);
  w.Add("cs2_" + nm, cs2_ncdm, true);
}

// ── Tensor output column titles ────────────────────────────────────────────

void NCDMSpecies::WriteTensorOutputColumnTitles(std::string& tensor_titles) const {
  const std::string& nm = name();
  std::string title     = "delta_" + nm;
  class_store_columntitle(tensor_titles, title.c_str(), true);
  title = "theta_" + nm;
  class_store_columntitle(tensor_titles, title.c_str(), true);
  title = "shear_" + nm;
  class_store_columntitle(tensor_titles, title.c_str(), true);
}

// Fused stress-energy: eliminates four virtual dispatches with one override.
// The FA branch handles all four quantities with a single ncdmfa check.
// The non-FA branch runs four sequential loops whose per-term expressions and
// accumulation order match the four individual methods exactly, so the result
// is bit-for-bit identical to the base-default delegation.
BaseSpecies::StressEnergyContribution NCDMSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_p_];
  if (layout.index_per_q.empty())
    return se;

  if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on) {
    const double rho_bg      = pvecback[index_bg_rho_];
    const double p_bg        = pvecback[index_bg_p_];
    const double pseudo_p_bg = pvecback[index_bg_pseudo_p_];
    const int idx0           = layout.index_per_q[0];
    const double w_ncdm      = p_bg / rho_bg;
    const double cg2_ncdm    = w_ncdm * (1.0 - 1.0 / (3.0 + 3.0 * w_ncdm) *
                                                   (3.0 * w_ncdm - 2.0 + pseudo_p_bg / p_bg));
    se.delta_rho             = rho_bg * y[idx0];
    se.rho_plus_p_theta      = (rho_bg + p_bg) * y[idx0 + 1];
    se.delta_p               = cg2_ncdm * rho_bg * y[idx0];
    se.rho_plus_p_shear      = (rho_bg + p_bg) * y[idx0 + 2];
    return se;
  }

  // Single fused pass over the q grid: epsilon and the (a0/a)^4 power are
  // computed once per bin instead of once per quantity. The per-bin sums are
  // ~ULP-equivalent (not bit-identical) to the four individual methods —
  // benign vectorization-reduction drift, validated to <<0.1% on Cl/P(k).
  const double a      = ppw->scalar_ctx.a;
  const double ma2    = M_ * M_ * a * a;
  const double factor = factor_ * std::pow(1. / a, 4);

  double s_drho = 0., s_theta = 0., s_dp = 0., s_shear = 0.;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q        = q_[iq];
    const double w0       = w_[iq];
    const double q2       = q * q;
    const double epsilon  = std::sqrt(q2 + ma2);
    const double q4_eps   = q2 * q2 / epsilon;
    const int idx         = layout.index_per_q[iq];
    s_drho               += q2 * epsilon * w0 * y[idx];
    s_theta              += q2 * q * w0 * y[idx + 1];
    s_dp                 += q4_eps * w0 * y[idx];
    s_shear              += q4_eps * w0 * y[idx + 2];
  }
  se.delta_rho        = s_drho * factor;
  se.rho_plus_p_theta = s_theta * ppw->scalar_ctx.k * factor;
  se.delta_p          = s_dp * factor / 3.;
  se.rho_plus_p_shear = 2. / 3. * factor * s_shear;
  return se;
}

// ── CopyPerturbationsAcrossSwitch (FA-collapse) ─────────────────────────────

void NCDMSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                const BaseSpecies::PerturbLayout& new_base,
                                                const double* old_y,
                                                double* new_y,
                                                const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const NCDMBaseSpecies::PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const NCDMBaseSpecies::PerturbLayout&>(new_base);
  if (old_l.q_size < 0 || new_l.q_size < 0)
    return;

  const bool fa_was_on = (old_l.q_size == 1 && old_l.l_max == 2);
  const bool fa_is_on  = (new_l.q_size == 1 && new_l.l_max == 2);

  if (fa_was_on == fa_is_on) {
    // Shape unchanged (no-switch or unchanged approximation): base slot copy.
    NCDMBaseSpecies::CopyPerturbationsAcrossSwitch(old_base, new_base, old_y, new_y, ctx);
    return;
  }

  if (!fa_was_on && fa_is_on) {
    // Collapse full hierarchy → (delta, theta, shear) for fluid approximation.
    const double a            = ctx.a;
    const double k            = ctx.k;
    const double rho          = Rho(ctx.pvecback);
    const double p            = P(ctx.pvecback);
    const double rho_plus_p   = rho + p;
    const double M_local      = M_;
    const double factor_local = factor_ * std::pow(1. / a, 4);

    const int idx_new = new_l.index_per_q[0];
    for (int l = 0; l <= 2; ++l)
      new_y[idx_new + l] = 0.0;

    double delta = 0., theta = 0., shear = 0.;
    for (int iq = 0; iq < old_l.q_size; ++iq) {
      const int idx_old    = old_l.index_per_q[iq];
      const double w0      = w_[iq];
      const double q       = q_[iq];
      const double epsilon = std::sqrt(q * q + a * a * M_local * M_local);

      delta += w0 * std::pow(q, 2) * epsilon * old_y[idx_old];
      theta += w0 * std::pow(q, 3) * old_y[idx_old + 1];
      shear += w0 * std::pow(q, 4) / epsilon * old_y[idx_old + 2];
    }
    delta *= factor_local / rho;
    theta *= k * factor_local / rho_plus_p;
    shear *= 2. / 3. * factor_local / rho_plus_p;

    new_y[idx_new]     = delta;
    new_y[idx_new + 1] = theta;
    new_y[idx_new + 2] = shear;
    return;
  }

  // fa_on → fa_off: unreachable in CLASS once FA is on.
  // Defensive fallback: leave new_y zero-initialized.
}
