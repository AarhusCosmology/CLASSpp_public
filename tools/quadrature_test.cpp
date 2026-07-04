#include "quadrature.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr std::array<double, 16> kFermiDiracPowerMoments = {
    6.931471805599453094172321e-1,
    8.224670334241132182362076e-1,
    1.803085354739391428099607,
    5.682196976983475505459019,
    2.333087449072582334245572e1,
    1.182661309556922124918127e2,
    7.146675503444378035227065e2,
    5.021014329337345412105538e3,
    4.024316207687752693655701e4,
    3.625314565172607857391608e5,
    3.627047810325762864638089e6,
    3.990712756635035863264820e7,
    4.789434217892724115449366e8,
    6.226642012479158901428670e9,
    8.717563672678887863206590e10,
    1.307654444554356358920468e12,
};

void fermi_dirac_distribution(void* /*params*/, double q, double* f0) {
  *f0 = 1. / (std::exp(q) + 1.);
}

void quartic_test_function(void* /*params*/, double q, double* value) {
  *value = q * q * q * q;
}

void check_power_moments(const std::vector<double>& x, const std::vector<double>& w) {
  for (int power = 0; power < static_cast<int>(kFermiDiracPowerMoments.size()); ++power) {
    double integral = 0.;
    for (size_t i = 0; i < x.size(); ++i)
      integral += w[i] * std::pow(x[i], power);
    const double reference = kFermiDiracPowerMoments[power];
    assert(std::abs(integral - reference) / reference < 2.e-10);
  }
}

void test_fundamental_rule() {
  constexpr int N = 8;
  std::vector<double> x(N), w(N);
  assert(compute_FermiDirac(x.data(), w.data(), N));
  for (int i = 0; i < N; ++i) {
    assert(x[i] > 0.);
    assert(w[i] > 0.);
    if (i > 0)
      assert(x[i - 1] < x[i]);
  }
  check_power_moments(x, w);

  std::vector<double> x_max(17), w_max(17);
  assert(compute_FermiDirac(x_max.data(), w_max.data(), 17));
  double normalization = 0.;
  for (double weight : w_max)
    normalization += weight;
  assert(std::abs(normalization - kFermiDiracPowerMoments[0]) < 2.e-10);

  assert(!compute_FermiDirac(x.data(), w.data(), 18));
}

void test_manual_selection() {
  constexpr int N = 8;
  std::vector<double> x(N), w(N), dq(N);
  GBQuadParams gb;
  assert(get_qsampling_manual(x.data(),
                              w.data(),
                              dq.data(),
                              N,
                              0.,
                              qm_FermiDirac,
                              nullptr,
                              0,
                              fermi_dirac_distribution,
                              nullptr,
                              gb));
  check_power_moments(x, w);
}

void test_auto_selection() {
  std::array<double, _QUADRATURE_MAX_> x{};
  std::array<double, _QUADRATURE_MAX_> w{};
  int N = 0;
  get_qsampling(x.data(),
                w.data(),
                &N,
                _QUADRATURE_MAX_,
                1.e-10,
                nullptr,
                0,
                quartic_test_function,
                fermi_dirac_distribution,
                nullptr);

  // A three-point FD rule integrates the quartic exactly. The Laguerre and
  // adaptive candidates require more points at this tolerance, so this also
  // pins the automatic selector's new candidate.
  assert(N == 3);
  double integral = 0.;
  for (int i = 0; i < N; ++i)
    integral += w[i] * std::pow(x[i], 4);
  assert(std::abs(integral - kFermiDiracPowerMoments[4]) / kFermiDiracPowerMoments[4] < 1.e-10);
}

}  // namespace

int main() {
  test_fundamental_rule();
  test_manual_selection();
  test_auto_selection();
  std::printf("Fermi-Dirac quadrature tests passed\n");
  return 0;
}
