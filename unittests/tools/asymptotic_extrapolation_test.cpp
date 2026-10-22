// AsymptoticExtrapolation continues a source past the last k of its table as
// s(k) k^p = A ln k + B + C/k, fitted to the points in [k_max/3, k_max] through the last one.
#include "asymptotic_extrapolation.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

namespace {

bool Close(double a, double b, double tol) {
  return std::abs(a - b) <= tol * std::abs(b);
}

// A table of s(k) = f(k) / k^power on log-spaced k up to k_max.
std::vector<double> Table(const std::vector<double>& k, double (*f)(double), int power) {
  std::vector<double> s;
  for (double x : k) {
    s.push_back(f(x) / std::pow(x, power));
  }
  return s;
}

std::vector<double> LogGrid(double k_min, double k_max, int per_decade) {
  std::vector<double> k;
  for (double x = k_min; x <= k_max * (1. + 1e-12); x *= std::pow(10., 1. / per_decade)) {
    k.push_back(x);
  }
  return k;
}

double Asymptotic(double k) {
  return 2.0 * std::log(k) + 0.7 - 0.3 / k;
}

void TestExactForTheAsymptoticForm() {
  for (int power : {0, 2}) {
    const auto k = LogGrid(0.01, 2., 10);
    const auto s = Table(k, Asymptotic, power);
    const AsymptoticExtrapolation ext(k.data(), s.data(), static_cast<int>(k.size()), 1, power);
    for (double x : {2., 5., 20., 70.}) {
      assert(Close(ext(x), Asymptotic(x) / std::pow(x, power), 1e-12));
    }
  }
  std::printf("exact for A ln k + B + C/k: ok\n");
}

double Wiggly(double k) {
  return std::log(k) + 0.05 * std::sin(10. * k);
}

void TestJoinsTheLastPoint() {
  const auto k = LogGrid(0.01, 2., 10);
  const auto s = Table(k, Wiggly, 2);
  const AsymptoticExtrapolation ext(k.data(), s.data(), static_cast<int>(k.size()), 1, 2);
  assert(Close(ext(k.back()), s.back(), 1e-14));
  std::printf("passes through the last point: ok\n");
}

void TestStrideReadsOneRowOfATable() {
  // The source tables are [tau][k]; a fit reads one row with stride 1, a column with
  // stride tau_size.
  const auto k = LogGrid(0.01, 2., 10);
  const int n  = static_cast<int>(k.size());
  std::vector<double> interleaved(2 * k.size());
  for (std::size_t i = 0; i < k.size(); ++i) {
    interleaved[2 * i]     = -1.;  // the other column
    interleaved[2 * i + 1] = Asymptotic(k[i]);
  }
  const AsymptoticExtrapolation ext(k.data(), interleaved.data() + 1, n, 2, 0);
  assert(Close(ext(30.), Asymptotic(30.), 1e-12));
  std::printf("stride selects one column: ok\n");
}

void TestFewPoints() {
  // Two points: the 1/k term cannot be fitted, the slope in ln k still is.
  const std::array<double, 2> k = {1., 2.};
  const std::array<double, 2> s = {std::log(1.) + 3., std::log(2.) + 3.};
  const AsymptoticExtrapolation two(k.data(), s.data(), 2, 1, 0);
  assert(Close(two(8.), std::log(8.) + 3., 1e-12));
  // One point: constant.
  const AsymptoticExtrapolation one(k.data() + 1, s.data() + 1, 1, 1, 0);
  assert(Close(one(8.), s[1], 1e-14));
  std::printf("two points fit the slope, one point is constant: ok\n");
}

}  // namespace

int main() {
  TestExactForTheAsymptoticForm();
  TestJoinsTheLastPoint();
  TestStrideReadsOneRowOfATable();
  TestFewPoints();
  std::printf("asymptotic extrapolation test passed\n");
  return 0;
}
