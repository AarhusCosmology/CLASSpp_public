// AcousticSpline (include/acoustic_spline.h): the l interpolation of the CMB spectra. It must pass
// through every sample, beat the cubic spline it replaces by a wide margin on an acoustic-like
// spectrum sampled as coarsely as CLASS samples l, and refuse what it cannot do.
// docs/superpowers/specs/2026-10-27-acoustic-l-interpolation-design.md
#include "acoustic_spline.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "arrays.h"
#include "common.h"

namespace {

/* the l list of the transfer module (flat space) */
std::vector<double> Multipoles(int l_linstep, int l_max) {
  std::vector<double> l{2.};
  int increment = std::max(static_cast<int>(l.back() * 0.12), 1);
  while (l.back() + increment < l_max && increment < l_linstep) {
    l.push_back(l.back() + increment);
    increment = std::max(static_cast<int>(l.back() * 0.12), 1);
  }
  while (l.back() + l_linstep <= l_max) {
    l.push_back(l.back() + l_linstep);
  }
  if (l.back() != l_max) {
    l.push_back(l_max);
  }
  return l;
}

/* an acoustic-like C_l: a damped, tilted envelope times peaks whose phase drifts with l, plus
   the odd/even modulation at half the frequency */
constexpr double kOmega = 0.0208;
double Spectrum(double l) {
  const double envelope = 5000. * std::pow(l / 200., 0.3) * std::exp(-std::pow(l / 1400., 2));
  const double phase    = kOmega * l + 0.4 + 2e-4 * l;
  const double D        = envelope * (1. + 0.45 * std::cos(phase) + 0.15 * std::cos(0.5 * phase));
  return D / (l * (l + 1.));
}

}  // namespace

int main() {
  const std::vector<double> l = Multipoles(50, 3500);
  const int n                 = static_cast<int>(l.size());
  std::vector<double> cl(n);
  for (int i = 0; i < n; i++) {
    cl[i] = Spectrum(l[i]);
  }

  assert(AcousticSpline::Applies(l.data(), n, cl.data(), 1, kOmega));
  const AcousticSpline spline(l.data(), n, cl.data(), 1, kOmega);

  /* through every sample */
  for (int i = 0; i < n; i++) {
    assert(std::fabs(spline(l[i]) - cl[i]) <= 1e-12 * std::fabs(cl[i]));
  }

  /* against CLASS's cubic spline (as the spectra module builds it) on the same samples */
  std::vector<double> ddcl(n);
  array_spline_table_lines(l.data(), n, cl.data(), 1, ddcl.data(), _SPLINE_EST_DERIV_);
  double worst_acoustic = 0.;
  double worst_cubic    = 0.;
  int last              = 0;
  for (int ll = 2; ll <= 3000; ll++) {
    double cubic = 0.;
    array_interpolate_spline(l.data(), n, cl.data(), ddcl.data(), 1, ll, &last, &cubic, 1);
    const double truth = Spectrum(ll);
    worst_acoustic     = std::max(worst_acoustic, std::fabs(spline(ll) / truth - 1.));
    worst_cubic        = std::max(worst_cubic, std::fabs(cubic / truth - 1.));
  }
  std::printf("largest relative error, l = 2..3000, sampled every 50: acoustic %.2e, cubic %.2e\n",
              worst_acoustic,
              worst_cubic);
  assert(worst_acoustic < 0.2 * worst_cubic);

  /* what it refuses: a zero spectrum, samples too sparse for omega, too few samples */
  const std::vector<double> zero(n, 0.);
  assert(!AcousticSpline::Applies(l.data(), n, zero.data(), 1, kOmega));
  assert(!AcousticSpline::Applies(l.data(), n, cl.data(), 1, 0.05));
  assert(!AcousticSpline::Applies(l.data(), 6, cl.data(), 1, kOmega));

  std::printf("acoustic spline tests passed\n");
  return 0;
}
