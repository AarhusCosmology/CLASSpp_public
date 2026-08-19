// Compiled WITHOUT -ffast-math (see CMakeLists.txt / Makefile). Do not add
// -ffast-math to this file and do not move this definition into the header: the
// whole reason it exists as a separate translation unit is that under
// -ffinite-math-only the compiler folds every non-finite test, including the
// bit-pattern one below. See include/nonfinite.h for the measurement.

#include "nonfinite.h"

#include <cstdint>
#include <cstring>

bool IsNonFinite(double x) {
  // memcpy is the standard-blessed type pun and compiles to a register move.
  // Exponent field all ones => Inf (zero mantissa) or NaN (non-zero mantissa);
  // both are wanted. Zero and denormals have a zero exponent => finite.
  std::uint64_t bits = 0;
  std::memcpy(&bits, &x, sizeof bits);
  return (bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL;
}

bool StepMakesNoProgress(double t, double h) {
  // Deliberately written as the rounding test itself rather than as
  // |h| < eps*|t|: this asks the machine the exact question the integrator cares
  // about, and stays correct at t = 0, across the denormal range and near the
  // exponent boundaries where a hand-rolled epsilon comparison drifts.
  return t + h == t;
}
