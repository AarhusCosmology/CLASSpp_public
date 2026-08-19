// rkdp45 must ABORT on a non-finite state, not grind against the step-size floor.
//
// Motivating failure (measured 2026-08-01, decaying-NCDM sector at Gamma = 1e8,
// dr_N_q = 30, k = 0.5): the perturbation state reached NaN at tau = 523.76 and the
// evolver then printed this line, identically, forever --
//
//   [rkdp45] x=5.23763103e+02 h=1.1654e-303 errmax=inf/1.0e-05 REJ worst=y[10]
//            (y=nan dy=-inf)
//
// 99.4% of steps rejected, h pinned. The only abort path was
// class_test(fabs(h) < hmin) with hmin = 100*DBL_MIN*|x| = 1.1654e-303 at that x --
// h had converged to exactly that value without going below it, so the guard never
// fired. Two multi-hour runs were burned on what should have been an immediate error.
//
// THE CHECK MUST BE BIT-LEVEL. classpp is compiled with -ffast-math
// (CMakeLists.txt), which implies -ffinite-math-only: the compiler may then assume
// NaN/Inf never occur and fold std::isnan/std::isfinite to a constant. The values
// still arise at runtime -- fast-math licenses the assumption, it does not prevent
// the arithmetic. So this test also pins that a NaN produced at RUN time (never a
// compile-time constant here, or the optimiser could fold the producer too) is
// actually caught.

#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "evolver_rkdp45.h"
#include "nonfinite.h"

namespace {

// A NaN the optimiser cannot constant-fold: both operands come from volatile
// storage, so 0.0/0.0 has to be evaluated at run time.
double RuntimeNaN() {
  volatile double zero_a = 0.0;
  volatile double zero_b = 0.0;
  return zero_a / zero_b;
}

double RuntimeInf() {
  volatile double one  = 1.0;
  volatile double zero = 0.0;
  return one / zero;
}

// y' = 1 until x > 0.5, then non-finite -- the shape of the real failure, where a
// healthy integration runs into a poisoned right-hand side partway through.
void derivs_goes_nan(double x, double* /*y*/, double* dy, void* /*p*/) {
  dy[0] = (x > 0.5) ? RuntimeNaN() : 1.0;
}

void derivs_finite(double /*x*/, double* /*y*/, double* dy, void* /*p*/) {
  dy[0] = 1.0;
}

// The evolver null-guards print_variables but NOT output (real callers always supply
// one), so a no-op is required rather than nullptr.
void output_noop(double /*x*/, double /*y*/[], double /*dy*/[], int /*index_x*/, void* /*p*/) {}

void test_bit_check_survives_fast_math() {
  // The load-bearing property: these must be true even where std::isnan is not
  // trustworthy. Integer inspection of the IEEE-754 exponent field cannot be
  // reasoned away by a floating-point assumption.
  assert(IsNonFinite(RuntimeNaN()));
  assert(IsNonFinite(RuntimeInf()));
  assert(IsNonFinite(-RuntimeInf()));
  assert(!IsNonFinite(0.0));
  assert(!IsNonFinite(-0.0));
  assert(!IsNonFinite(1.0));
  assert(!IsNonFinite(-1e300));
  // Denormals are finite: the stalled run sat at h ~ 1e-303, so a check that
  // rejected denormals would fire on healthy-but-tiny states.
  assert(!IsNonFinite(1.1654e-303));
  assert(!IsNonFinite(5e-324));
  std::printf("nonfinite: bit check ok (NaN, +/-Inf caught; denormals finite)\n");
}

void test_evolver_aborts_on_nonfinite() {
  double y               = 0.;
  int used               = 1;
  std::vector<double> xs = {0., 1.};
  bool threw             = false;
  std::string msg;
  try {
    evolver_rkdp45(derivs_goes_nan,
                   0.,
                   1.,
                   &y,
                   &used,
                   1,
                   nullptr,
                   1e-5,
                   1e-10,
                   nullptr,
                   0.,
                   xs.data(),
                   static_cast<int>(xs.size()),
                   output_noop,
                   nullptr,
                   nullptr);
  }
  catch (const std::runtime_error& e) {
    threw = true;
    msg   = e.what();
  }
  assert(threw);  // must ERROR, not hang and not silently return garbage
  // The message has to name the variable, or the next person is back to bisecting
  // a 1000-entry state vector by hand.
  assert(msg.find("non-finite") != std::string::npos);
  assert(msg.find("y[0]") != std::string::npos);
  std::printf("nonfinite: evolver aborted as required: %s\n", msg.c_str());
}

void test_finite_run_still_completes() {
  // Guard against a check that fires on healthy runs.
  double y               = 0.;
  int used               = 1;
  std::vector<double> xs = {0., 1.};
  evolver_rkdp45(derivs_finite,
                 0.,
                 1.,
                 &y,
                 &used,
                 1,
                 nullptr,
                 1e-5,
                 1e-10,
                 nullptr,
                 0.,
                 xs.data(),
                 static_cast<int>(xs.size()),
                 output_noop,
                 nullptr,
                 nullptr);
  assert(y > 0.999 && y < 1.001);  // y' = 1 over [0,1]
  std::printf("nonfinite: finite run unaffected (y=%.6f)\n", y);
}

// StepMakesNoProgress lives in the same -fno-fast-math TU and for the same reason:
// `t + h == t` is algebraically `h == 0`, which -ffast-math may substitute, at which
// point the test only catches h == 0 exactly -- a value the step controller never
// produces, because it floors h at minimum_variation.
//
// Motivating failure (measured 2026-08-02, same sector at Gamma = 1e8 under ndf15):
// 93% of 4.87M attempted steps sat at one tau = 6422.974 with h = 7.08e-16 against
// eps*tau = 1.4e-12. The steps were ACCEPTED (local error ~1e-12 vs rtol 1e-5), so
// ndf15's own "Step size too small" test never fired -- it lives inside the
// err > rtol branch. That floor is also ABSOLUTE (minimum_variation = DBL_EPSILON
// = 2.2e-16) where the question is RELATIVE to |t|: h was above the floor and still
// 2000x too small to move t.
void test_no_progress_check_survives_fast_math() {
  // Values from the real stall. h is far above DBL_EPSILON, so an absolute floor
  // misses it; t + h == t all the same.
  volatile double t_stall = 6422.974;
  volatile double h_stall = 7.08e-16;
  assert(StepMakesNoProgress(t_stall, h_stall));

  // A step that does advance t must NOT be flagged, including one that is tiny in
  // absolute terms but large relative to t.
  volatile double t_small = 1e-8;
  volatile double h_small = 1e-16;
  assert(!StepMakesNoProgress(t_small, h_small));
  assert(!StepMakesNoProgress(1.0, 1e-15));
  assert(!StepMakesNoProgress(6422.974, 1e-9));

  // h == 0 is the degenerate case the folded form WOULD still catch; it must stay
  // caught, so the guard is strictly stronger than the optimiser's rewrite.
  assert(StepMakesNoProgress(1.0, 0.0));

  // t == 0 is where a hand-rolled |h| < eps*|t| comparison breaks (it would flag
  // every step); the rounding test gets it right.
  assert(!StepMakesNoProgress(0.0, 5e-324));
}

}  // namespace

int main() {
  test_bit_check_survives_fast_math();
  test_no_progress_check_survives_fast_math();
  test_evolver_aborts_on_nonfinite();
  test_finite_run_still_completes();
  std::printf("evolver nonfinite test passed\n");
  return 0;
}
