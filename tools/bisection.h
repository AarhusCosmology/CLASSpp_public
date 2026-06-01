#pragma once

/**
 * Minimal bisection helpers, header-only, replacing hand-rolled loops across the
 * modules.
 *
 * bisect_value: continuous bracketing. `lo`,`hi` bracket the root; `pred(mid)`
 *   returns true when `mid` is on the `hi` side of the root. Iterates until
 *   (hi - lo) <= tol; returns the final midpoint.
 *
 * bisect_index: integer bracketing. `lo`,`hi` are indices with pred(lo)==false,
 *   pred(hi)==true; returns the smallest index in (lo,hi] where pred flips true,
 *   matching the `(hi - lo) > 1` table-bracketing loops.
 */
template <typename Pred>
double bisect_value(double lo, double hi, double tol, Pred pred) {
  while ((hi - lo) > tol) {
    double mid = 0.5 * (lo + hi);
    if (pred(mid)) {
      hi = mid;
    }
    else {
      lo = mid;
    }
  }
  return 0.5 * (lo + hi);
}

template <typename Pred>
int bisect_index(int lo, int hi, Pred pred) {
  while ((hi - lo) > 1) {
    int mid = (lo + hi) / 2;
    if (pred(mid)) {
      hi = mid;
    }
    else {
      lo = mid;
    }
  }
  return hi;
}
