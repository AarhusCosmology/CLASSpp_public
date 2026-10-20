#include "bisection.h"

#include <cassert>
#include <cmath>
#include <cstdio>

int main() {
  // Value-based: find root of f(x)=x-2 in [0,10] to tol 1e-9.
  // predicate(mid) is true when mid is on the upper (hi) side of the root.
  double root = bisect_value(0.0, 10.0, 1e-9, [](double x) { return x > 2.0; });
  assert(std::fabs(root - 2.0) < 1e-6);

  // Integer-bracket: smallest index where table[i] >= 5 in 0..9 (table[i]=i),
  // matching the (hi - lo) > 1 table-bracketing loops.
  int idx = bisect_index(0, 9, [](int i) { return i >= 5; });
  assert(idx == 5);

  // Regression (#395): tol below one ulp of the bracket must terminate rather
  // than spin. tol_tau_approx is an absolute time in Mpc, so tol=1e-12 against
  // tau ~ 1e4 (one ulp ~ 1.8e-12) is reachable from the input file, and the
  // approximation-switch bisection used to hang forever on it.
  double tau = bisect_value(1.0e4, 2.0e4, 1e-12, [](double t) { return t > 1.5e4; });
  assert(std::fabs(tau - 1.5e4) < 1e-8);

  // Degenerate tolerances must terminate too, at the best available bracket.
  double zero_tol = bisect_value(1.0e4, 2.0e4, 0.0, [](double t) { return t > 1.5e4; });
  assert(std::fabs(zero_tol - 1.5e4) < 1e-8);

  std::printf("bisection tests passed\n");
  return 0;
}
