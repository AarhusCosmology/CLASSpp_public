// hyperspherical_Hermite_evaluate (locations found once, node derivatives from tables) must give
// the cubic Hermite interpolants hyperspherical_Hermite_interpolation<4, ...> computes on the fly:
// the same polynomials, so they may differ in rounding only. Checked for every multipole, for
// Phi, dPhi and d2Phi, at arguments inside the table, on nodes, at its ends and outside it, and
// for a suffix of the located arguments (how the transfer module uses them).
// docs/superpowers/specs/2026-10-27-transfer-hermite-kernel-design.md
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "hyperspherical.h"

int main() {
  std::vector<int> l;
  for (int ll = 2; ll <= 10; ll++) {
    l.push_back(ll);
  }
  while (l.back() < 1200) {
    l.push_back(std::min(l.back() + 30, 1200));
  }
  const double x_max            = 2000.;
  const HyperInterpStruct table = HyperInterpStruct::OnFixedSpacing(0,
                                                                    1.,
                                                                    static_cast<int>(l.size()),
                                                                    l.data(),
                                                                    1e-5,
                                                                    x_max,
                                                                    8.,
                                                                    l.back() + 1,
                                                                    1e-8);

  /* increasing arguments, as the transfer module passes them: below the table, dense and
     sparse stretches, two nodes, the last node, and beyond the table */
  std::vector<double> x{0., 1e-6};
  for (int j = 0; j < 4000; j++) {
    x.push_back(0.5 + j * 0.37);
  }
  x.push_back(table.x[1234]);
  x.push_back(table.x[1235]);
  std::sort(x.begin(), x.end());
  x.push_back(table.x.back());
  x.push_back(table.x.back() + 1.);
  const int n = static_cast<int>(x.size());

  HermiteLocations where;
  where.resize(n);
  hyperspherical_Hermite_locate(&table, n, x.data(), where);

  /* Not 1e-15: near x_min the radial equation gives the third derivative at a node as a
     difference of terms ~1e4 times larger (l(l+1) y cot/sin^2 with sin^2 ~ 1e-10), so d2Phi
     there carries ~1e-12 of rounding in either code, and they round differently. An index or
     formula slip is an O(1) difference. */
  constexpr double kTolerance = 1e-10;
  double worst                = 0.;
  std::vector<double> phi(n);
  std::vector<double> dphi(n);
  std::vector<double> d2phi(n);
  std::vector<double> phi0(n);
  std::vector<double> dphi0(n);
  std::vector<double> d2phi0(n);
  const int first = 1500;
  for (int k = 0; k < table.l_size; k++) {
    hyperspherical_Hermite_evaluate<true, true, true>(&table,
                                                      k,
                                                      where,
                                                      0,
                                                      n,
                                                      phi.data(),
                                                      dphi.data(),
                                                      d2phi.data());
    hyperspherical_Hermite_interpolation<4, true, true, true>(&table,
                                                              n,
                                                              k,
                                                              x.data(),
                                                              phi0.data(),
                                                              dphi0.data(),
                                                              d2phi0.data());
    for (const auto& [mine, reference] :
         {std::pair{&phi, &phi0}, std::pair{&dphi, &dphi0}, std::pair{&d2phi, &d2phi0}}) {
      double scale = 0.;
      for (double v : *reference) {
        scale = std::max(scale, std::fabs(v));
      }
      for (int j = 0; j < n; j++) {
        worst = std::max(worst, std::fabs((*mine)[j] - (*reference)[j]) / scale);
        if (std::fabs((*mine)[j] - (*reference)[j]) > kTolerance * scale) {
          std::printf("l = %d, x = %.6g: %.17g vs %.17g\n",
                      l[k],
                      x[j],
                      (*mine)[j],
                      (*reference)[j]);
          assert(false);
        }
      }
    }
    assert(phi[0] == 0. && phi[n - 1] == 0. && d2phi[n - 1] == 0.);

    /* a suffix of the located arguments, and the Phi-only variant */
    std::vector<double> tail(n - first);
    hyperspherical_Hermite_evaluate<true, false, false>(&table,
                                                        k,
                                                        where,
                                                        first,
                                                        n - first,
                                                        tail.data(),
                                                        nullptr,
                                                        nullptr);
    for (int j = 0; j < n - first; j++) {
      assert(tail[j] == phi[first + j]);
    }
  }

  std::printf("Hermite evaluate tests passed (largest difference %.1e of the largest value)\n",
              worst);
  return 0;
}
