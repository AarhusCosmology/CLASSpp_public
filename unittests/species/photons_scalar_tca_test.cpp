/* The scalar photon hierarchy that tight coupling hands over, checked against the
   hierarchy's own equations rather than against transcribed numbers.

   While tight coupling is on, perturb_sources uses P^(0) = 5 s_2 shear_g / 8. In the
   optimal basis, with F_2 = 2 s_2 shear_g, the O(1) balances of the G_0 and G_2
   equations as kappa' -> infinity are

     G_0 = 4 P,   G_2 = 4/5 P,   with   P = (G_0 + G_2 + F_2) / 8,

   which give P = 5 s_2 shear_g / 8. The l = 1 and l = 3 equations then fix the
   next order:

     G_1 = k/(3 kappa') (G_0 - 2 s_2 G_2),
     G_3 = k/(7 kappa') 3 s_3 G_2,
     F_3 = k/(7 kappa') 3 s_3 F_2.

   In flat space s_l = 1 and a hand-over that drops the s_2 factors still passes;
   with curvature it makes P jump at the switch.

   The TAM basis (PhotonsTamSpecies) hands over the same state in its own variables,
   with P = (Theta_2 - sqrt(6) E_2)/10 and Theta_2 = 5/2 s_2 shear_g: its E_2 and l = 3
   equations give E_2 = -sqrt(6) P, Theta_3 = (k/kappa') 3/2 s_3 s_2 shear_g and
   E_3 = {}_2kappa_3^0/(5 kappa') E_2, with {}_2kappa_3^0 = sqrt(5 (k^2 - 8K)). */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "background.h"
#include "constants.h"
#include "perturbations.h"
#include "photons_optimal.h"
#include "photons_tam.h"

namespace {

constexpr int kLmax      = 12;
constexpr int kLmaxPol   = 10;
constexpr double kK      = 0.05;  /* Mpc^-1 */
constexpr double kDkappa = 2.4e3; /* Mpc^-1, so that k/kappa' is small but not negligible */
constexpr double kShear  = 3.7e-4;

bool Close(double a, double b) {
  return std::fabs(a - b) <= 1e-12 * std::max({std::fabs(a), std::fabs(b), 1e-30});
}

/* K/k^2 = 0.1: s_2 = sqrt(0.7), far enough from 1 that a missing factor shows. */
void FillCurvedWorkspace(perturb_workspace& ppw, double K_over_k2) {
  ppw.s_l.assign(kLmax + 2, 1.);
  for (int l = 0; l < kLmax + 2; ++l) {
    ppw.s_l[l] = std::sqrt(std::max(1. - K_over_k2 * (l * l - 1.), 0.));
  }
  ppw.tca_shear_g = kShear;
}

PhotonsOptimalSpecies::PerturbLayout OptimalLayout() {
  PhotonsOptimalSpecies::PerturbLayout lay;
  int i          = 0;
  lay.idx_delta  = i++;
  lay.idx_theta  = i++;
  lay.idx_shear  = i++;
  lay.idx_l3     = i;
  i             += kLmax - 2;
  lay.idx_pol0   = i++;
  lay.idx_pol1   = i++;
  lay.idx_pol2   = i++;
  lay.idx_pol3   = i;
  lay.l_max      = kLmax;
  lay.l_max_pol  = kLmaxPol;
  return lay;
}

void test_optimal_hand_over(double K_over_k2) {
  background pba{};
  const PhotonsOptimalSpecies g(pba);
  perturb_workspace ppw;
  FillCurvedWorkspace(ppw, K_over_k2);
  const auto lay = OptimalLayout();
  std::vector<double> y(lay.idx_pol0 + kLmaxPol + 1, 0.);

  g.StartScalarHierarchy(lay, y.data(), kK, kDkappa, &ppw);

  const double s2 = ppw.s_l[2];
  const double s3 = ppw.s_l[3];
  const double F2 = 2. * s2 * y[lay.idx_shear];
  const double P  = g.ScalarPolarisationSource(lay, y.data(), &ppw);
  const double e  = kK / kDkappa;

  std::printf("K/k^2 = %.2f: P = %.6e, tight coupling's 5 s_2 shear/8 = %.6e\n",
              K_over_k2,
              P,
              5. * s2 * kShear / 8.);
  assert(Close(y[lay.idx_shear], kShear));
  assert(Close(P, 5. * s2 * kShear / 8.));
  assert(Close(y[lay.idx_pol0], 4. * P));
  assert(Close(y[lay.idx_pol2], 4. / 5. * P));
  assert(Close(y[lay.idx_pol1], e / 3. * (y[lay.idx_pol0] - 2. * s2 * y[lay.idx_pol2])));
  assert(Close(y[lay.idx_pol3], e / 7. * 3. * s3 * y[lay.idx_pol2]));
  assert(Close(y[lay.idx_l3], e / 7. * 3. * s3 * F2));
}

void test_tam_hand_over(double K_over_k2) {
  background pba{};
  pba.K = K_over_k2 * kK * kK;
  const PhotonsTamSpecies g(pba);
  perturb_workspace ppw;
  FillCurvedWorkspace(ppw, K_over_k2);

  PhotonsTamSpecies::PerturbLayout lay;
  lay.idx_delta = 0;
  lay.idx_theta = 1;
  lay.idx_shear = 2;
  lay.idx_l3    = 3;
  lay.idx_e2    = lay.idx_delta + kLmax + 1;
  lay.l_max     = kLmax;
  lay.l_max_pol = kLmaxPol;
  std::vector<double> y(lay.idx_e2 + kLmaxPol - 1, 0.);

  g.StartScalarHierarchy(lay, y.data(), kK, kDkappa, &ppw);

  const double s2      = ppw.s_l[2];
  const double s3      = ppw.s_l[3];
  const double P       = g.ScalarPolarisationSource(lay, y.data(), &ppw);
  const double e       = kK / kDkappa;
  const double kappa23 = std::sqrt(5. * (kK * kK - 8. * pba.K));

  std::printf("TAM, K/k^2 = %.2f: P = %.6e, tight coupling's 5 s_2 shear/8 = %.6e\n",
              K_over_k2,
              P,
              5. * s2 * kShear / 8.);
  assert(Close(y[lay.idx_shear], kShear));
  assert(Close(P, 5. * s2 * kShear / 8.));
  assert(Close(y[lay.idx_e2], -_SQRT6_ * P));
  assert(Close(y[lay.idx_l3], e * 3. / 2. * s3 * s2 * kShear));
  assert(Close(y[lay.idx_e2 + 1], kappa23 / (5. * kDkappa) * y[lay.idx_e2]));
}

}  // namespace

int main() {
  test_optimal_hand_over(0.);
  test_optimal_hand_over(0.1);
  test_optimal_hand_over(-0.1);
  test_tam_hand_over(0.);
  test_tam_hand_over(0.1);
  test_tam_hand_over(-0.1);
  std::printf("photons scalar tight-coupling hand-over tests passed\n");
  return 0;
}
