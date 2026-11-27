#include "photons_tam.h"

#include <cmath>

#include "errors.h"
#include "output_table.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

/* The total angular momentum hierarchy, ported from Cyril Pitrou's implementation in
   class_public v3.4.0. Equation numbers: (33) and (34) of astro-ph/9709066 for the
   hierarchy, (5), (9) and (10) of 2005.12119 for the coefficients and the truncation. */

// ── Coefficients ──────────────────────────────────────────────────────────────

double PhotonsTamSpecies::Kappa(int s, int m, int l, double q2, double K) {
  if (l <= m || l <= s) {
    return 0.;
  }
  const double l2 = static_cast<double>(l) * l;
  return std::sqrt((l2 - m * m) * (l2 - s * s) / l2 * std::max(q2 - K * l2, 0.));
}

void PhotonsTamSpecies::FillCoefficients(Coefficients& c, double k, double K, int m, int l_max) {
  const double q2 = k * k + K * (1. + m);
  c.q             = std::sqrt(std::max(q2, 0.));
  c.kappa0.assign(l_max + 1, 0.);
  c.kappa2.assign(l_max + 1, 0.);
  for (int l = 0; l <= l_max; ++l) {
    c.kappa0[l] = Kappa(0, m, l, q2, K);
    c.kappa2[l] = Kappa(2, m, l, q2, K);
  }
}

const PhotonsTamSpecies::Coefficients& PhotonsTamSpecies::CoefficientsFor(
    const perturb_workspace* ppw, double k, int m, int l_max) const {
  auto* scratch = static_cast<Scratch*>(ppw->species_scratch[collection_index_].get());
  if (scratch->k != k || scratch->m != m ||
      static_cast<int>(scratch->c.kappa2.size()) < l_max + 1) {
    FillCoefficients(scratch->c, k, pba_.K, m, l_max);
    scratch->k = k;
    scratch->m = m;
  }
  return scratch->c;
}

// ── Perturbation index registration ──────────────────────────────────────────

void PhotonsTamSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                    perturb_vector* /*pv*/,
                                                    const precision* ppr,
                                                    int& index_pt,
                                                    const perturb_workspace* ppw,
                                                    int /*gauge*/) {
  auto& layout     = static_cast<PerturbLayout&>(base);
  layout           = PerturbLayout{};
  layout.l_max     = ppr->l_max_g;
  layout.l_max_pol = ppr->l_max_pol_g;

  /* RSA active: photons are handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_on)) {
    return;
  }
  layout.idx_delta = index_pt++;
  layout.idx_theta = index_pt++;

  /* TCA active: only delta and theta */
  if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
    return;
  }
  layout.idx_shear  = index_pt++;
  layout.idx_l3     = index_pt; /* Theta_3 .. Theta_l_max */
  index_pt         += layout.l_max - 2;
  layout.idx_e2     = index_pt; /* E_2 .. E_l_max_pol */
  index_pt         += layout.l_max_pol - 1;
}

void PhotonsTamSpecies::RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& /*base*/,
                                                          perturb_vector* /*pv*/,
                                                          const precision* /*ppr*/,
                                                          int& /*index_pt*/,
                                                          const perturb_workspace* /*ppw*/,
                                                          int /*gauge*/) {
  class_stop_severe("vector modes are not implemented with hierarchy = tam");
}

void PhotonsTamSpecies::RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                          perturb_vector* /*pv*/,
                                                          const precision* ppr,
                                                          int& index_pt,
                                                          const perturb_workspace* ppw,
                                                          int /*gauge*/) {
  auto& layout     = static_cast<PerturbLayout&>(base);
  layout           = PerturbLayout{};
  layout.l_max     = ppr->l_max_g_ten + kTensorExtraMultipoles;
  layout.l_max_pol = ppr->l_max_pol_g_ten + kTensorExtraMultipoles;

  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_on)) {
    return;
  }
  if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
    return;
  }
  /* m = 2 has no l = 0, 1 */
  layout.idx_t2  = index_pt; /* Theta_2 .. Theta_l_max */
  index_pt      += layout.l_max - 1;
  layout.idx_e2  = index_pt; /* E_2 .. E_l_max_pol */
  index_pt      += layout.l_max_pol - 1;
  layout.idx_b2  = index_pt; /* B_2 .. B_l_max_pol */
  index_pt      += layout.l_max_pol - 1;
}

// ── Scalar modes ──────────────────────────────────────────────────────────────

void PhotonsTamSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
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
  const double cotKgen            = ctx.cotKgen;
  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];
  const double theta_g = y[layout.idx_theta];
  const double* kappa2 = CoefficientsFor(ppw, k, 0, layout.l_max_pol).kappa2.data();
  const double P0      = ScalarPolarisationSource(layout, y, ppw);
  const double* T      = y + layout.idx_delta; /* T[l] = Theta_l for l >= 3 */
  double* dT           = dy + layout.idx_delta;
  const double* E      = y + layout.idx_e2 - 2; /* E[l] = E_l, l >= 2 */
  double* dE           = dy + layout.idx_e2 - 2;
  const int l_max      = layout.l_max;
  const int l_max_pol  = layout.l_max_pol;

  /* Temperature: (33) of astro-ph/9709066 with m = 0, Theta_1 = theta_g/k,
     Theta_2 = 5/2 s_2 shear_g, and q {}_0kappa_l^0 = k l s_l. */

  /* shear_g (l = 2), with Theta_3 = 7/4 F_3 */
  dy[layout.idx_shear] = 4. / 15. * (theta_g + ctx.metric_shear) -
                         6. / 35. * k * s_l[3] / s_l[2] * T[3] -
                         dkappa * (y[layout.idx_shear] - 2. / 5. / s_l[2] * P0);

  /* l = 3 */
  dT[3] = k * (3. / 2. * s_l[3] * s_l[2] * y[layout.idx_shear] - 4. / 9. * s_l[4] * T[4]) -
          dkappa * T[3];

  /* l = 4 .. l_max - 1 */
  for (int l = 4; l < l_max; l++) {
    dT[l] = k * (l / (2. * l - 1.) * s_l[l] * T[l - 1] -
                 (l + 1.) / (2. * l + 3.) * s_l[l + 1] * T[l + 1]) -
            dkappa * T[l];
  }

  /* l = l_max: (9) of 2005.12119 */
  {
    const int l = l_max;
    dT[l] = k * ((2. * l + 1.) / (2. * l - 1.) * s_l[l] * T[l - 1] - (1. + l) * cotKgen * T[l]) -
            dkappa * T[l];
  }

  /* Polarisation: (34) of astro-ph/9709066 with m = 0 (no B) */
  dE[2] = -kappa2[3] / 7. * E[3] - dkappa * (E[2] + _SQRT6_ * P0);

  for (int l = 3; l < l_max_pol; l++) {
    dE[l] = kappa2[l] / (2. * l - 1.) * E[l - 1] - kappa2[l + 1] / (2. * l + 3.) * E[l + 1] -
            dkappa * E[l];
  }

  /* l = l_max_pol: (10) of 2005.12119 */
  {
    const int l = l_max_pol;
    dE[l]       = kappa2[l] * (2. * l + 1.) / (2. * l - 1.) / (l - 2.) * E[l - 1] -
                  (l + 3.) * k * cotKgen * E[l] - dkappa * E[l];
  }
}

double PhotonsTamSpecies::ScalarPolarisationSource(const PhotonsSpecies::PerturbLayout& base,
                                                   const double* y,
                                                   const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* P^(0) = (Theta_2 - sqrt(6) E_2)/10, with Theta_2 = 5/2 s_2 shear_g */
  return (5. / 2. * ppw->s_l[2] * y[layout.idx_shear] - _SQRT6_ * y[layout.idx_e2]) / 10.;
}

void PhotonsTamSpecies::StartScalarHierarchy(const PhotonsSpecies::PerturbLayout& base,
                                             double* y,
                                             double k,
                                             double dkappa,
                                             const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);

  /* The tight-coupling solution of the hierarchy's own equations, the same state that
     PhotonsOptimalSpecies hands over: E_2 = -sqrt(6) P with P = 5 s_2 shear_g / 8 at first
     order; Theta_3 and E_3 at second order, from the l = 3 equations. (class_public sets
     Theta_3 and E_2; E_3 is added here because the optimal basis hands over its counterparts,
     G_1 and G_3.) */
  y[layout.idx_shear]   = ppw->tca_shear_g;
  const double s2_shear = ppw->s_l[2] * y[layout.idx_shear];

  y[layout.idx_l3]     = 3. / 2. * k / dkappa * ppw->s_l[3] * s2_shear;
  y[layout.idx_e2]     = -5. * _SQRT6_ / 8. * s2_shear;
  const double q2      = k * k + pba_.K;
  y[layout.idx_e2 + 1] = Kappa(2, 0, 3, q2, pba_.K) / (5. * dkappa) * y[layout.idx_e2];
}

void PhotonsTamSpecies::AddPolarisationOutput(OutputTable& row,
                                              const PhotonsSpecies::PerturbLayout& base,
                                              const double* y,
                                              double /*k*/,
                                              double /*dkappa*/,
                                              const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  double E2_g        = 0.;
  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_off)) {
    if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
      /* the state StartScalarHierarchy hands over */
      E2_g = -5. * _SQRT6_ / 8. * ppw->s_l[2] * ppw->tca_shear_g;
    }
    else {
      E2_g = y[layout.idx_e2];
    }
  }
  row.Set("E2_g", E2_g);
}

void PhotonsTamSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                          const perturb_workspace* ppw,
                                          int* used_in_sources) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (ppw->approx[ppw->index_ap_rsa] != static_cast<int>(rsa_off) ||
      ppw->approx[ppw->index_ap_tca] != static_cast<int>(tca_off)) {
    return;
  }
  for (int idx = layout.idx_l3; idx <= layout.idx_delta + layout.l_max; ++idx) {
    used_in_sources[idx] = 0;
  }
  for (int idx = layout.idx_e2 + 1; idx <= layout.idx_e2 + layout.l_max_pol - 2; ++idx) {
    used_in_sources[idx] = 0;
  }
}

// ── Vector modes ──────────────────────────────────────────────────────────────

void PhotonsTamSpecies::PerturbVectorDerivs(
    const BaseSpecies::PerturbLayout& /*base*/,
    double /*tau*/,
    const double* /*y*/,
    double* /*dy*/,
    const perturb_parameters_and_workspace& /*ppaw*/) const {
  class_stop_severe("vector modes are not implemented with hierarchy = tam");
}

// ── Tensor modes ──────────────────────────────────────────────────────────────

void PhotonsTamSpecies::PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                                            double /*tau*/,
                                            const double* y,
                                            double* dy,
                                            const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (layout.idx_t2 < 0) {
    return;
  }

  const perturb_workspace* ppw = ppaw.ppw;
  const double k               = ppaw.k;
  const double cotKgen         = ppw->cotKgen;
  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];
  const int l_max       = layout.l_max;
  const int l_max_pol   = layout.l_max_pol;
  const Coefficients& c = CoefficientsFor(ppw, k, 2, std::max(l_max, l_max_pol));
  const double q        = c.q;
  const double* kappa0  = c.kappa0.data();
  const double* kappa2  = c.kappa2.data();
  const double P2       = TensorPolarisationSource(layout, y);
  const double* T       = y + layout.idx_t2 - 2; /* T[l] = Theta_l, l >= 2 */
  double* dT            = dy + layout.idx_t2 - 2;
  const double* E       = y + layout.idx_e2 - 2;
  double* dE            = dy + layout.idx_e2 - 2;
  const double* B       = y + layout.idx_b2 - 2;
  double* dB            = dy + layout.idx_b2 - 2;

  /* Temperature: (33) of astro-ph/9709066 with m = 2, sourced by the gravitational waves */
  dT[2] = -kappa0[3] / 7. * T[3] - dkappa * (T[2] - P2) - y[ppw->pv->index_pt_gwdot];

  for (int l = 3; l < l_max; l++) {
    dT[l] = kappa0[l] / (2. * l - 1.) * T[l - 1] - kappa0[l + 1] / (2. * l + 3.) * T[l + 1] -
            dkappa * T[l];
  }

  /* l = l_max: (9) of 2005.12119 */
  {
    const int l = l_max;
    dT[l]       = kappa0[l] * (2. * l + 1.) / (2. * l - 1.) / (l - 2.) * T[l - 1] -
                  (l + 3.) * k * cotKgen * T[l] - dkappa * T[l];
  }

  /* Polarisation: (34) of astro-ph/9709066 with m = 2; E and B mix through q */
  dE[2] = -q * 4. / 6. * B[2] - kappa2[3] / 7. * E[3] - dkappa * (E[2] + _SQRT6_ * P2);
  dB[2] = q * 4. / 6. * E[2] - kappa2[3] / 7. * B[3] - dkappa * B[2];

  for (int l = 3; l < l_max_pol; l++) {
    const double mix = q * 4. / l / (l + 1.);
    dE[l]            = kappa2[l] / (2. * l - 1.) * E[l - 1] - mix * B[l] -
                       kappa2[l + 1] / (2. * l + 3.) * E[l + 1] - dkappa * E[l];
    dB[l]            = kappa2[l] / (2. * l - 1.) * B[l - 1] + mix * E[l] -
                       kappa2[l + 1] / (2. * l + 3.) * B[l + 1] - dkappa * B[l];
  }

  /* l = l_max_pol: (10) of 2005.12119, whose E-B term, m q / l, differs from the hierarchy's */
  {
    const int l       = l_max_pol;
    const double down = kappa2[l] * (2. * l + 1.) / (2. * l - 1.) / (l - 2.);
    dE[l] = down * E[l - 1] + q * 2. / l * B[l] - (l + 3.) * k * cotKgen * E[l] - dkappa * E[l];
    dB[l] = down * B[l - 1] - q * 2. / l * E[l] - (l + 3.) * k * cotKgen * B[l] - dkappa * B[l];
  }
}

double PhotonsTamSpecies::TensorPolarisationSource(const PhotonsSpecies::PerturbLayout& base,
                                                   const double* y) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* P^(2) = (Theta_2 - sqrt(6) E_2)/10, (37) of astro-ph/9709066 */
  return (y[layout.idx_t2] - _SQRT6_ * y[layout.idx_e2]) / 10.;
}

double PhotonsTamSpecies::TensorGwSource(const PhotonsSpecies::PerturbLayout& base,
                                         const double* y,
                                         double a2,
                                         double rho_g) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* Pi = 8/5 Theta_2, (40) of astro-ph/9709066 */
  return 8. / 5. * a2 * rho_g * y[layout.idx_t2];
}

void PhotonsTamSpecies::StartTensorHierarchy(const PhotonsSpecies::PerturbLayout& base,
                                             double* y,
                                             const TensorTcaClosure& tca) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* The closure's state (F_0, G_0, every other optimal multipole zero) in this basis:
     Theta_2 = -F_0/sqrt(6), E_2 = -G_0, which keeps P^(2); B_2 and the higher multipoles
     are higher order in 1/kappa' and stay at zero. At leading order, Theta_2 = -4/3 h'/kappa'
     and E_2 = sqrt(2/3) h'/kappa' (C. Pitrou). */
  y[layout.idx_t2] = -tca.F0 / _SQRT6_;
  y[layout.idx_e2] = -tca.G0;
}

void PhotonsTamSpecies::AddTensorPerturbationOutput(OutputTable& row,
                                                    const PhotonsSpecies::PerturbLayout& base,
                                                    const double* y,
                                                    const TensorTcaClosure* tca) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  double Theta2_g    = 0.;
  double E2_g        = 0.;
  double B2_g        = 0.;
  if (tca != nullptr) {
    Theta2_g = -tca->F0 / _SQRT6_;
    E2_g     = -tca->G0;
  }
  else if (layout.idx_t2 >= 0) {
    Theta2_g = y[layout.idx_t2];
    E2_g     = y[layout.idx_e2];
    B2_g     = y[layout.idx_b2];
  }
  row.Set("Theta2_g", Theta2_g);
  row.Set("E2_g", E2_g);
  row.Set("B2_g", B2_g);
}

void PhotonsTamSpecies::MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& base,
                                                const perturb_workspace* ppw,
                                                int* used_in_sources) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (ppw->approx[ppw->index_ap_rsa] != static_cast<int>(rsa_off) ||
      ppw->approx[ppw->index_ap_tca] != static_cast<int>(tca_off)) {
    return;
  }
  for (int idx = layout.idx_t2 + 1; idx <= layout.idx_t2 + layout.l_max - 2; ++idx) {
    used_in_sources[idx] = 0;
  }
  for (int idx = layout.idx_e2 + 1; idx <= layout.idx_e2 + layout.l_max_pol - 2; ++idx) {
    used_in_sources[idx] = 0;
  }
  for (int idx = layout.idx_b2; idx <= layout.idx_b2 + layout.l_max_pol - 2; ++idx) {
    used_in_sources[idx] = 0;
  }
}

// ── Switch-copy hook ──────────────────────────────────────────────────────────

void PhotonsTamSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                      const BaseSpecies::PerturbLayout& new_base,
                                                      const double* old_y,
                                                      double* new_y,
                                                      const PerturbSwitchContext& /*ctx*/) const {
  const auto& o = static_cast<const PerturbLayout&>(old_base);
  const auto& n = static_cast<const PerturbLayout&>(new_base);

  /* Tensor modes: Theta_2.., E_2.., B_2.. when both vectors evolve them */
  if (o.idx_t2 >= 0 && n.idx_t2 >= 0) {
    for (int l = 2; l <= n.l_max; ++l) {
      new_y[n.idx_t2 + l - 2] = old_y[o.idx_t2 + l - 2];
    }
    for (int l = 2; l <= n.l_max_pol; ++l) {
      new_y[n.idx_e2 + l - 2] = old_y[o.idx_e2 + l - 2];
      new_y[n.idx_b2 + l - 2] = old_y[o.idx_b2 + l - 2];
    }
    return;
  }

  /* Scalar modes */
  if (o.idx_delta < 0 || n.idx_delta < 0) {
    return;
  }
  new_y[n.idx_delta] = old_y[o.idx_delta];
  new_y[n.idx_theta] = old_y[o.idx_theta];
  if (o.idx_shear >= 0 && n.idx_shear >= 0) {
    new_y[n.idx_shear] = old_y[o.idx_shear];
    for (int l = 3; l <= n.l_max; ++l) {
      new_y[n.idx_delta + l] = old_y[o.idx_delta + l];
    }
    for (int l = 2; l <= n.l_max_pol; ++l) {
      new_y[n.idx_e2 + l - 2] = old_y[o.idx_e2 + l - 2];
    }
  }
}
