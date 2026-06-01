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

  std::printf("bisection tests passed\n");
  return 0;
}
