#include "photons_optimal.h"

#include "baryons.h"
#include "output_table.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

// ── Perturbation index registration ──────────────────────────────────────────

void PhotonsOptimalSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                        perturb_vector* /*pv*/,
                                                        const precision* ppr,
                                                        int& index_pt,
                                                        const perturb_workspace* ppw,
                                                        int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  /* Initialize all photon indices to sentinel -1 so that any early return
     (RSA, TCA) leaves them in a well-defined state that guards can check. */
  layout.idx_delta = layout.idx_theta = layout.idx_shear = layout.idx_l3 = -1;
  layout.idx_pol0 = layout.idx_pol1 = layout.idx_pol2 = layout.idx_pol3 = -1;

  /* Copy l_max values from precision into layout for later use in derivs/switch. */
  layout.l_max     = ppr->l_max_g;
  layout.l_max_pol = ppr->l_max_pol_g;

  /* RSA active: photons are handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_on)) {
    return;
  }

  layout.idx_delta = index_pt++;
  layout.idx_theta = index_pt++;

  /* TCA active: only delta and theta; shear/polarization are analytic */
  if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
    return;
  }

  /* Full hierarchy: shear (l=2), l3..l_max_g, polarization */
  layout.idx_shear = index_pt++;

  /* l=3 .. l_max_g → (l_max_g - 2) slots */
  if (layout.l_max >= 3) {
    layout.idx_l3  = index_pt;
    index_pt      += layout.l_max - 2;
  }

  layout.idx_pol0 = index_pt++;
  layout.idx_pol1 = index_pt++;
  layout.idx_pol2 = index_pt++;

  /* pol l=3 .. l_max_pol_g → (l_max_pol_g - 2) slots */
  if (layout.l_max_pol >= 3) {
    layout.idx_pol3  = index_pt;
    index_pt        += layout.l_max_pol - 2;
  }
}

void PhotonsOptimalSpecies::RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                              perturb_vector* /*pv*/,
                                                              const precision* ppr,
                                                              int& index_pt,
                                                              const perturb_workspace* ppw,
                                                              int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = layout.idx_theta = layout.idx_shear = layout.idx_l3 = -1;
  layout.idx_pol0 = layout.idx_pol1 = layout.idx_pol2 = layout.idx_pol3 = -1;

  /* Photons read their own l_max values from ppr — no caller pre-setup. */
  layout.l_max     = ppr->l_max_g_ten;
  layout.l_max_pol = ppr->l_max_pol_g_ten;

  /* In vector mode there is no TCA shortcut: either full hierarchy or nothing. */
  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_on)) {
    return;
  }
  if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
    return;
  }

  layout.idx_delta  = index_pt++;
  layout.idx_theta  = index_pt++;
  layout.idx_shear  = index_pt++;
  layout.idx_l3     = index_pt;
  index_pt         += layout.l_max - 2;

  layout.idx_pol0  = index_pt++;
  layout.idx_pol1  = index_pt++;
  layout.idx_pol2  = index_pt++;
  layout.idx_pol3  = index_pt;
  index_pt        += layout.l_max_pol - 2;
}

void PhotonsOptimalSpecies::RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                              perturb_vector* /*pv*/,
                                                              const precision* ppr,
                                                              int& index_pt,
                                                              const perturb_workspace* ppw,
                                                              int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = layout.idx_theta = layout.idx_shear = layout.idx_l3 = -1;
  layout.idx_pol0 = layout.idx_pol1 = layout.idx_pol2 = layout.idx_pol3 = -1;

  /* Photons read their own l_max values from ppr — no caller pre-setup. */
  layout.l_max     = ppr->l_max_g_ten;
  layout.l_max_pol = ppr->l_max_pol_g_ten;

  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_on)) {
    return;
  }
  if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
    return;
  }

  layout.idx_delta  = index_pt++;
  layout.idx_theta  = index_pt++;
  layout.idx_shear  = index_pt++;
  layout.idx_l3     = index_pt;
  index_pt         += layout.l_max - 2;

  layout.idx_pol0  = index_pt++;
  layout.idx_pol1  = index_pt++;
  layout.idx_pol2  = index_pt++;
  layout.idx_pol3  = index_pt;
  index_pt        += layout.l_max_pol - 2;
}

// ── PerturbDerivs ─────────────────────────────────────────────────────────────

void PhotonsOptimalSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                          double /*tau*/,
                                          const double* y,
                                          double* dy,
                                          const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);

  /* RSA active: photons are handled analytically */
  if (layout.idx_delta < 0) {
    return;
  }
  /* delta_g' and theta_g'; under tight coupling nothing else is evolved */
  if (DensityAndVelocityDerivs(layout, y, dy, ppaw)) {
    return;
  }

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ppaw.k;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;
  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];
  const double theta_g = y[layout.idx_theta];

  /* Full hierarchy (TCA off) */
  const double P0 = ScalarPolarisationSource(layout, y, ppw);

  dy[layout.idx_shear] = 0.5 * (8. / 15. * (theta_g + metric_shear) -
                                3. / 5. * k * s_l[3] / s_l[2] * y[layout.idx_l3] -
                                dkappa * (2. * y[layout.idx_shear] - 4. / 5. / s_l[2] * P0));

  /* l = 3 */
  {
    const int l       = 3;
    dy[layout.idx_l3] = k / (2. * l + 1.) *
                            (l * s_l[l] * 2. * s_l[2] * y[layout.idx_shear] -
                             (l + 1.) * s_l[l + 1] * y[layout.idx_l3 + 1]) -
                        dkappa * y[layout.idx_l3];
  }

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < layout.l_max; l++) {
    dy[layout.idx_delta + l] = k / (2. * l + 1.) *
                                   (l * s_l[l] * y[layout.idx_delta + l - 1] -
                                    (l + 1.) * s_l[l + 1] * y[layout.idx_delta + l + 1]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* l = l_max_g (free-streaming truncation) */
  {
    const int l              = layout.l_max;
    dy[layout.idx_delta + l] = k * (s_l[l] * y[layout.idx_delta + l - 1] -
                                    (1. + l) * cotKgen * y[layout.idx_delta + l]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* Polarization l = 0 */
  dy[layout.idx_pol0] = -k * y[layout.idx_pol0 + 1] - dkappa * (y[layout.idx_pol0] - 4. * P0);

  /* Polarization l = 1 */
  dy[layout.idx_pol1] = k / 3. * (y[layout.idx_pol1 - 1] - 2. * s_l[2] * y[layout.idx_pol1 + 1]) -
                        dkappa * y[layout.idx_pol1];

  /* Polarization l = 2 */
  dy[layout.idx_pol2] =
      k / 5. * (2. * s_l[2] * y[layout.idx_pol2 - 1] - 3. * s_l[3] * y[layout.idx_pol2 + 1]) -
      dkappa * (y[layout.idx_pol2] - 4. / 5. * P0);

  /* Polarization l = 3 .. l_max_pol_g - 1 */
  for (int l = 3; l < layout.l_max_pol; l++) {
    dy[layout.idx_pol0 + l] = k / (2. * l + 1.) *
                                  (l * s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (l + 1.) * s_l[l + 1] * y[layout.idx_pol0 + l + 1]) -
                              dkappa * y[layout.idx_pol0 + l];
  }

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l             = layout.l_max_pol;
    dy[layout.idx_pol0 + l] = k * (s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (1. + l) * cotKgen * y[layout.idx_pol0 + l]) -
                              dkappa * y[layout.idx_pol0 + l];
  }
}

void PhotonsOptimalSpecies::PerturbVectorDerivs(
    const BaseSpecies::PerturbLayout& base,
    double /*tau*/,
    const double* y,
    double* dy,
    const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout           = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw = ppaw.ppw;
  const double* s_l            = ppw->s_l.data();

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (layout.idx_delta < 0) {
    return;
  }

  const auto gauge     = ppaw.perturbations_module->GetPerturbs()->gauge;
  const double k       = ppaw.k;
  const double k2      = k * k;
  const double cotKgen = ppw->cotKgen;
  const double ssqrt3  = sqrt(1. - 2. * pba_.K / k2);

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[layout.idx_delta];
  const double theta_g = y[layout.idx_theta];
  const double shear_g = y[layout.idx_shear];

  /* P^{(1)}: polarization combination (Eq. B.23 in arXiv:1305.3261) */
  const double P1 = -_SQRT6_ / 40. *
                    (4. / (3. * k) * theta_g + y[layout.idx_delta + 3] + 2. * y[layout.idx_pol0] +
                     10. / 7. * y[layout.idx_pol2] - 4. / 7. * y[layout.idx_pol0 + 4]);

  /* Vector-mode baryon theta_b lives in the Baryons layout (registered in the
     module's _vectors_ branch via `b_vec_lay.idx_theta`). */
  const auto& species   = ppaw.perturbations_module->all_species_;
  const size_t b_phv_i  = species.index_of("Baryons");
  const auto& b_phv_lay = static_cast<const BaryonsSpecies::PerturbLayout&>(
      *ppw->pv->species_layouts[b_phv_i]);
  const double theta_b_vec = y[b_phv_lay.idx_theta];

  if (gauge == possible_gauges::synchronous) {
    dy[layout.idx_delta] = -4. / 3. * theta_g - dkappa * (delta_g + 2. * _SQRT2_ * theta_b_vec);

    dy[layout.idx_theta] = k2 * (delta_g / 4. - s_l[2] * shear_g) -
                           dkappa * (theta_g + 4.0 / _SQRT6_ * P1) +
                           4.0 / (3.0 * _SQRT2_) * ssqrt3 * y[ppw->pv->index_pt_hv_prime];
  }
  else { /* newtonian */
    dy[layout.idx_delta] = -4. / 3. * theta_g - dkappa * (delta_g + 2. * _SQRT2_ * theta_b_vec) -
                           2. * _SQRT2_ * ppw->pvecmetric[ppw->index_mt_V_prime];

    dy[layout.idx_theta] = k2 * (delta_g / 4. - s_l[2] * shear_g) -
                           dkappa * (theta_g + 4.0 / _SQRT6_ * P1);
  }

  /* photon shear (F_2/2) */
  dy[layout.idx_shear] = 4. / 15. * s_l[2] * theta_g -
                         3. / 10. * k * s_l[3] * y[layout.idx_shear + 1] - dkappa * shear_g;

  /* l = 3 */
  dy[layout.idx_l3] = k / 7. * (6. * s_l[3] * shear_g - 4. * s_l[4] * y[layout.idx_l3 + 1]) -
                      dkappa * y[layout.idx_l3];

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < layout.l_max; l++) {
    dy[layout.idx_delta + l] = k / (2. * l + 1.) *
                                   (l * s_l[l] * y[layout.idx_delta + l - 1] -
                                    (l + 1.) * s_l[l + 1] * y[layout.idx_delta + l + 1]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* l = l_max_g (truncation) */
  {
    const int l              = layout.l_max;
    dy[layout.idx_delta + l] = k * (s_l[l] * y[layout.idx_delta + l - 1] -
                                    (1. + l) * cotKgen * y[layout.idx_delta + l]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* Polarization l = 0 (G_0) */
  dy[layout.idx_pol0] = -k * y[layout.idx_pol0 + 1] - dkappa * (y[layout.idx_pol0] - _SQRT6_ * P1);

  /* Polarization l = 1 .. l_max_pol_g - 1 */
  for (int l = 1; l < layout.l_max_pol; l++) {
    dy[layout.idx_pol0 + l] = k / (2. * l + 1.) *
                                  (l * s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (l + 1.) * s_l[l + 1] * y[layout.idx_pol0 + l + 1]) -
                              dkappa * y[layout.idx_pol0 + l];
  }

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l             = layout.l_max_pol;
    dy[layout.idx_pol0 + l] = k * (s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (1. + l) * cotKgen * y[layout.idx_pol0 + l]) -
                              dkappa * y[layout.idx_pol0 + l];
  }
}

void PhotonsOptimalSpecies::PerturbTensorDerivs(
    const BaseSpecies::PerturbLayout& base,
    double /*tau*/,
    const double* y,
    double* dy,
    const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout           = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw = ppaw.ppw;
  const double* s_l            = ppw->s_l.data();

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (layout.idx_delta < 0) {
    return;
  }

  const double k       = ppaw.k;
  const double k2      = k * k;
  const double cotKgen = ppw->cotKgen;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[layout.idx_delta];
  const double theta_g = y[layout.idx_theta];
  const double shear_g = y[layout.idx_shear];

  /* P^{(2)}: polarization combination */
  const double P2 = TensorPolarisationSource(layout, y);

  dy[layout.idx_delta] = -4. / 3. * theta_g - dkappa * (delta_g + _SQRT6_ * P2) +
                         _SQRT6_ * y[ppw->pv->index_pt_gwdot]; /* gravitational wave source */

  dy[layout.idx_theta] = k2 * (delta_g / 4. - s_l[2] * shear_g) - dkappa * theta_g;

  dy[layout.idx_shear] = 4. / 15. * s_l[2] * theta_g -
                         3. / 10. * k * s_l[3] * y[layout.idx_shear + 1] - dkappa * shear_g;

  /* l = 3 */
  dy[layout.idx_l3] = k / 7. * (6. * s_l[3] * shear_g - 4. * s_l[4] * y[layout.idx_l3 + 1]) -
                      dkappa * y[layout.idx_l3];

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < layout.l_max; l++) {
    dy[layout.idx_delta + l] = k / (2. * l + 1.) *
                                   (l * s_l[l] * y[layout.idx_delta + l - 1] -
                                    (l + 1.) * s_l[l + 1] * y[layout.idx_delta + l + 1]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* l = l_max_g (truncation) */
  {
    const int l              = layout.l_max;
    dy[layout.idx_delta + l] = k * (s_l[l] * y[layout.idx_delta + l - 1] -
                                    (1. + l) * cotKgen * y[layout.idx_delta + l]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* Polarization l = 0 (G_0) */
  dy[layout.idx_pol0] = -k * y[layout.idx_pol0 + 1] - dkappa * (y[layout.idx_pol0] - _SQRT6_ * P2);

  /* Polarization l = 1 .. l_max_pol_g - 1 */
  for (int l = 1; l < layout.l_max_pol; l++) {
    dy[layout.idx_pol0 + l] = k / (2. * l + 1.) *
                                  (l * s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (l + 1.) * s_l[l + 1] * y[layout.idx_pol0 + l + 1]) -
                              dkappa * y[layout.idx_pol0 + l];
  }

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l             = layout.l_max_pol;
    dy[layout.idx_pol0 + l] = k * (s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (1. + l) * cotKgen * y[layout.idx_pol0 + l]) -
                              dkappa * y[layout.idx_pol0 + l];
  }
}

// ── Switch-copy hook ──────────────────────────────────────────────────────────

void PhotonsOptimalSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y,
    double* new_y,
    const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);

  if (old_l.idx_delta < 0 || new_l.idx_delta < 0) {
    return;
  }

  // delta, theta, shear, l3, pol0..pol3, all multipoles up to l_max.
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];

  if (old_l.idx_shear >= 0 && new_l.idx_shear >= 0) {
    new_y[new_l.idx_shear] = old_y[old_l.idx_shear];
    if (old_l.idx_l3 >= 0 && new_l.idx_l3 >= 0) {
      new_y[new_l.idx_l3] = old_y[old_l.idx_l3];
      for (int l = 4; l <= new_l.l_max; ++l) {
        new_y[new_l.idx_delta + l] = old_y[old_l.idx_delta + l];
      }
    }
    if (new_l.idx_pol0 >= 0 && old_l.idx_pol0 >= 0) {
      new_y[new_l.idx_pol0] = old_y[old_l.idx_pol0];
      new_y[new_l.idx_pol1] = old_y[old_l.idx_pol1];
      new_y[new_l.idx_pol2] = old_y[old_l.idx_pol2];
      if (old_l.idx_pol3 >= 0 && new_l.idx_pol3 >= 0) {
        new_y[new_l.idx_pol3] = old_y[old_l.idx_pol3];
        for (int l = 4; l <= new_l.l_max_pol; ++l) {
          new_y[new_l.idx_pol0 + l] = old_y[old_l.idx_pol0 + l];
        }
      }
    }
  }
}

// ── MarkUsedInSources ─────────────────────────────────────────────────────────

void PhotonsOptimalSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                              const perturb_workspace* ppw,
                                              int* used_in_sources) const {
  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  /* Photon temperature l>=3 and pol l=1,3+ are not needed in source evaluation
     when both rsa and tca are off. */
  if (ppw->approx[ppw->index_ap_rsa] != static_cast<int>(rsa_off)) {
    return;
  }
  if (ppw->approx[ppw->index_ap_tca] != static_cast<int>(tca_off)) {
    return;
  }

  for (int idx = g_lay.idx_l3; idx <= g_lay.idx_delta + g_lay.l_max; ++idx) {
    used_in_sources[idx] = 0;
  }

  used_in_sources[g_lay.idx_pol1] = 0;
  for (int idx = g_lay.idx_pol3; idx <= g_lay.idx_pol0 + g_lay.l_max_pol; ++idx) {
    used_in_sources[idx] = 0;
  }
}

// ── MarkTensorUsedInSources ────────────────────────────────────────────────────

void PhotonsOptimalSpecies::MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& base,
                                                    const perturb_workspace* ppw,
                                                    int* used_in_sources) const {
  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  /* In tensor mode, we only need temperature l=0,2,4 and pol l=0,2,4
     when both rsa and tca are off. */
  if (ppw->approx[ppw->index_ap_rsa] != static_cast<int>(rsa_off)) {
    return;
  }
  if (ppw->approx[ppw->index_ap_tca] != static_cast<int>(tca_off)) {
    return;
  }

  used_in_sources[g_lay.idx_theta] = 0;
  used_in_sources[g_lay.idx_l3]    = 0;
  for (int idx = g_lay.idx_delta + 5; idx <= g_lay.idx_delta + g_lay.l_max; ++idx) {
    used_in_sources[idx] = 0;
  }

  used_in_sources[g_lay.idx_pol1] = 0;
  used_in_sources[g_lay.idx_pol3] = 0;
  for (int idx = g_lay.idx_pol0 + 5; idx <= g_lay.idx_pol0 + g_lay.l_max_pol; ++idx) {
    used_in_sources[idx] = 0;
  }
}

// ── Hooks the perturbations module calls ──────────────────────────────────────

void PhotonsOptimalSpecies::StartScalarHierarchy(const PhotonsSpecies::PerturbLayout& base,
                                                 double* y,
                                                 double k,
                                                 double dkappa,
                                                 const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);

  /* tight-coupling approximation for shear_g (previously
     computed in perturb_derivs: perturb_derivs is always
     called at the end of generic_evolver, in order to update
     all quantities in ppw to the time at which the
     approximation is switched off) */
  y[layout.idx_shear] = ppw->tca_shear_g;

  /* The tight-coupling solution of the l = 0 .. 3 equations, with F_2 = 2 s_2 shear_g:
     G_0 = 4 P and G_2 = 4/5 P at first order, so that P = (G_0 + G_2 + F_2)/8 = 5 s_2 shear_g/8,
     the source perturb_sources uses while tight coupling is on; G_1, G_3 and F_3 at second
     order. Every term carries s_2 through F_2 (as in class_public since v3.0); without it P
     jumps at the switch in curved space. */
  const double s2_shear = ppw->s_l[2] * y[layout.idx_shear];
  /* second-order tight-coupling approximation for l=3 */
  y[layout.idx_l3] = 6. / 7. * k / dkappa * ppw->s_l[3] * s2_shear;
  /* first-order tight-coupling approximation for polarization, l=0 */
  y[layout.idx_pol0] = 2.5 * s2_shear;
  /* second-order tight-coupling approximation for polarization, l=1 */
  y[layout.idx_pol1] = k / dkappa * (5. - 2. * ppw->s_l[2]) / 6. * s2_shear;
  /* first-order tight-coupling approximation for polarization, l=2 */
  y[layout.idx_pol2] = 0.5 * s2_shear;
  /* second-order tight-coupling approximation for polarization, l=3 */
  y[layout.idx_pol3] = k / dkappa * 3. * ppw->s_l[3] / 14. * s2_shear;
}

double PhotonsOptimalSpecies::ScalarPolarisationSource(const PhotonsSpecies::PerturbLayout& base,
                                                       const double* y,
                                                       const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* P^(0) = (G_0 + G_2 + F_2) / 8, with F_2 = 2 s_2 shear_g */
  return (y[layout.idx_pol0] + y[layout.idx_pol2] + 2. * ppw->s_l[2] * y[layout.idx_shear]) / 8.;
}

void PhotonsOptimalSpecies::StartTensorHierarchy(const PhotonsSpecies::PerturbLayout& base,
                                                 double* y,
                                                 const TensorTcaClosure& tca) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* Multipoles above l = 0 are higher order in 1/kappa' and stay at zero. */
  y[layout.idx_delta] = tca.F0;
  y[layout.idx_pol0]  = tca.G0;
}

double PhotonsOptimalSpecies::TensorPolarisationSource(const PhotonsSpecies::PerturbLayout& base,
                                                       const double* y) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  return TensorP2(y[layout.idx_delta],
                  y[layout.idx_shear],
                  y[layout.idx_delta + 4],
                  y[layout.idx_pol0],
                  y[layout.idx_pol2],
                  y[layout.idx_pol0 + 4]);
}

double PhotonsOptimalSpecies::TensorGwSource(const PhotonsSpecies::PerturbLayout& base,
                                             const double* y,
                                             double a2,
                                             double rho_g) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* Pi = -4 sqrt(6) [F_0/15 + 2/21 F_2 + 1/35 F_4] (1305.3261), with F_2 = 2 shear */
  return -_SQRT6_ * 4 * a2 * rho_g *
         (1. / 15. * y[layout.idx_delta] + 4. / 21. * y[layout.idx_shear] +
          1. / 35. * y[layout.idx_l3 + 1]);
}

void PhotonsOptimalSpecies::AddTensorPerturbationOutput(OutputTable& row,
                                                        const PhotonsSpecies::PerturbLayout& base,
                                                        const double* y,
                                                        const TensorTcaClosure* tca) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  double delta_g     = 0.;
  double shear_g     = 0.;
  double l4_g        = 0.;
  double pol0_g      = 0.;
  double pol2_g      = 0.;
  double pol4_g      = 0.;
  if (tca != nullptr) {
    delta_g = tca->F0;
    pol0_g  = tca->G0;
  }
  else if (layout.idx_delta >= 0) {
    delta_g = y[layout.idx_delta];
    shear_g = y[layout.idx_shear];
    l4_g    = y[layout.idx_delta + 4];
    pol0_g  = y[layout.idx_pol0];
    pol2_g  = y[layout.idx_pol2];
    pol4_g  = y[layout.idx_pol0 + 4];
  }
  row.Set("delta_g", delta_g);
  row.Set("shear_g", shear_g);
  row.Set("l4_g", l4_g);
  row.Set("pol0_g", pol0_g);
  row.Set("pol2_g", pol2_g);
  row.Set("pol4_g", pol4_g);
}

void PhotonsOptimalSpecies::AddPolarisationOutput(OutputTable& row,
                                                  const PhotonsSpecies::PerturbLayout& base,
                                                  const double* y,
                                                  double k,
                                                  double dkappa,
                                                  const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  double pol0_g      = 0.;
  double pol1_g      = 0.;
  double pol2_g      = 0.;
  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_off)) {
    if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
      /* the state StartScalarHierarchy hands over */
      const double s2_shear = ppw->s_l[2] * ppw->tca_shear_g;
      pol0_g                = 2.5 * s2_shear;
      pol1_g                = k / dkappa * (5. - 2. * ppw->s_l[2]) / 6. * s2_shear;
      pol2_g                = 0.5 * s2_shear;
    }
    else {
      pol0_g = y[layout.idx_pol0];
      pol1_g = y[layout.idx_pol1];
      pol2_g = y[layout.idx_pol2];
    }
  }
  row.Set("pol0_g", pol0_g);
  row.Set("pol1_g", pol1_g);
  row.Set("pol2_g", pol2_g);
}
