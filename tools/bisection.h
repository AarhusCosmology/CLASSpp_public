#pragma once

/**
 * Minimal bisection helpers, header-only, replacing hand-rolled loops across the
 * modules.
 *
 * bisect_value: continuous bracketing. `lo`,`hi` bracket the root; `pred(mid)`
 *   returns true when `mid` is on the `hi` side of the root. Iterates until
 *   (hi - lo) <= tol; returns the final midpoint.
 *
 *   `tol` is an ABSOLUTE width, so a caller can legitimately ask for a
 *   tolerance below the spacing of doubles near the bracket (e.g. tol = 1e-12
 *   on tau ~ 1e4 Mpc, where one ulp is ~1.8e-12). The loop therefore also
 *   stops once the bracket can no longer be subdivided, and returns the best
 *   bracket instead of spinning forever.
 *
 *   The midpoint is 0.5*(lo + hi) rather than lo + 0.5*(hi - lo): the latter is
 *   the overflow-safe spelling, but it moves the midpoint by an ulp, and the
 *   evolver amplifies that into 6e-5 on P(k). Every caller brackets a physical
 *   quantity (conformal times of order 1e-5..1e5 Mpc, a momentum), so lo + hi
 *   cannot overflow; a caller passing values near DBL_MAX would have to switch
 *   to the other spelling and accept the drift.
 *
 * bisect_index: integer bracketing. `lo`,`hi` are indices with pred(lo)==false,
 *   pred(hi)==true; returns the smallest index in (lo,hi] where pred flips true,
 *   matching the `(hi - lo) > 1` table-bracketing loops.
 */
template <typename Pred>
double bisect_value(double lo, double hi, double tol, Pred pred) {
  while ((hi - lo) > tol) {
    double mid = 0.5 * (lo + hi);
    /* tol below one ulp of the bracket: mid can no longer separate from an
       endpoint, so every further iteration would repeat this one. */
    if (mid <= lo || mid >= hi)
      break;
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
