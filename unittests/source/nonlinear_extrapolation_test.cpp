/**
 * Unit tests for the default continuation of a source beyond the computed k
 * (extrap_max_scaled), which halofit's sigma integrals and HMcode integrate over.
 *
 * The continuation is source ~ ln(a k), with a fitted through the last two points
 * of the table, (k_previous, source_previous) and (k_max, source_max). That is a
 * straight line in ln k through the two points, and these tests pin exactly that:
 * the function is finite, passes through both points and is affine in ln k.
 *
 * The cases that matter are those where the last two sources are nearly equal.
 * Then ln(a k_max) = source_max ln(k_max/k_previous) / (source_max - source_previous)
 * is large, and a = exp(ln(a k_max)) / k_max over- or underflows. That happens
 * when the table ends near an extremum of an oscillating source: dark acoustic
 * oscillations that the k grid does not resolve. The values below are those of
 * delta_m at z = 0 at three such points of an IDM-DRMD model (f_idm_drmd = 0.7,
 * delta_Neff_drmd near 0.0025), where HMcode stopped with x_min = nan and
 * halofit dropped the non-linear correction.
 */

#include <cassert>
#include <cmath>
#include <cstdio>

#include "nonlinear_module.h"

namespace {

/* The last two points of the k table at k_per_decade_for_pk = 10. */
constexpr double kMax      = 5.79305;
constexpr double kPrevious = 4.60165;  // kMax / 10^(1/10)

/* HMcode extrapolates to hmcode_max_k_extra = 1e6 /Mpc. */
constexpr double kFar = 1.e6;

int failures = 0;

void check(bool ok, const char* what, double got, double want) {
  if (!ok) {
    std::printf("  FAIL %-58s got %.16e want %.16e\n", what, got, want);
    failures++;
  }
}

double extrapolate(double k, double source_max, double source_previous) {
  return NonlinearModule::ExtrapolateMaxScaled(k, kMax, kPrevious, source_max, source_previous);
}

/* The continuation is finite at every k up to HMcode's reach, passes through the
   two table points and is affine in ln k. */
void check_is_the_line_through_the_two_points(double source_max, double source_previous) {
  const double scale = std::fabs(source_max);

  for (double k = kMax; k <= kFar; k *= 1.7) {
    const double s = extrapolate(k, source_max, source_previous);
    check(std::isfinite(s), "finite up to hmcode_max_k_extra", s, source_max);
  }

  const double at_max      = extrapolate(kMax, source_max, source_previous);
  const double at_previous = extrapolate(kPrevious, source_max, source_previous);
  check(std::fabs(at_max - source_max) <= 1.e-12 * scale,
        "passes through (k_max, source_max)",
        at_max,
        source_max);
  check(std::fabs(at_previous - source_previous) <= 1.e-12 * scale,
        "passes through (k_previous, source_previous)",
        at_previous,
        source_previous);

  /* Affine in ln k: the second difference on an even grid in ln k vanishes. */
  for (double u : {0.3, 2., 6.}) {
    const double s0 = extrapolate(kMax, source_max, source_previous);
    const double s1 = extrapolate(kMax * std::exp(u), source_max, source_previous);
    const double s2 = extrapolate(kMax * std::exp(2. * u), source_max, source_previous);
    check(std::fabs(s2 - 2. * s1 + s0) <= 1.e-12 * scale, "affine in ln k", s2 - 2. * s1 + s0, 0.);
  }
}

/* A CDM-like end of the table, |source| growing like ln k: ln(a k_max) = 5.3, so
   a is a moderate number and ln(a k) / ln(a k_max) can be evaluated as written. */
void test_cdm_like_source_is_source_max_ln_ak_over_ln_ak_max() {
  const double source_max      = -2.14962567e4;
  const double source_previous = -2.05652582e4;
  check_is_the_line_through_the_two_points(source_max, source_previous);

  const long double ln_ak_max = static_cast<long double>(source_max) * std::log(kMax / kPrevious) /
                                (static_cast<long double>(source_max) - source_previous);
  for (double k : {6., 17., 1.e3, kFar}) {
    const long double ln_ak = ln_ak_max + std::log(static_cast<long double>(k) / kMax);
    const double want       = static_cast<double>(source_max * ln_ak / ln_ak_max);
    const double got        = extrapolate(k, source_max, source_previous);
    check(std::fabs(got - want) <= 1.e-12 * std::fabs(want),
          "source_max ln(a k) / ln(a k_max)",
          got,
          want);
  }
}

/* |source| still grows, by 1.4 parts in 10^4: ln(a k_max) = +1672, so a overflows. */
void test_nearly_flat_growing_end_is_finite() {
  check_is_the_line_through_the_two_points(-2.0e4, -2.0e4 * (1. - 1.377e-4));
}

/* |source| falls, by 1.3 parts in 10^4: ln(a k_max) = -1826, so a underflows to 0. */
void test_nearly_flat_falling_end_is_finite() {
  check_is_the_line_through_the_two_points(-2.0e4, -2.0e4 * (1. + 1.261e-4));
}

/* ln(a k_max) = +40000. */
void test_almost_exactly_flat_end_is_finite() {
  check_is_the_line_through_the_two_points(-2.0e4, -2.0e4 * (1. - 5.765e-6));
}

/* Equal end values: the line is the constant source_max (the closed form was 0/0). */
void test_equal_end_values_continue_as_a_constant() {
  for (double k : {kMax, 17., kFar}) {
    const double s = extrapolate(k, -2.0e4, -2.0e4);
    check(s == -2.0e4, "equal end values give a constant", s, -2.0e4);
  }
}

/* The continuation is a continuous function of the table: as source_previous
   passes through source_max (ln(a k_max) through +-infinity), it does not jump. */
void test_continuous_through_equal_end_values() {
  const double source_max = -2.0e4;
  for (double eps : {1.e-6, 1.e-9, 1.e-12}) {
    const double above = extrapolate(kFar, source_max, source_max * (1. + eps));
    const double below = extrapolate(kFar, source_max, source_max * (1. - eps));
    /* The slope changes by 2 eps source_max / ln(k_max / k_previous), over
       ln(kFar / k_max) = 12: a change of order 100 eps in units of source_max. */
    check(std::fabs(above - below) <= 200. * eps * std::fabs(source_max),
          "continuous as source_previous crosses source_max",
          above,
          below);
  }
}

}  // namespace

int main() {
  test_cdm_like_source_is_source_max_ln_ak_over_ln_ak_max();
  test_nearly_flat_growing_end_is_finite();
  test_nearly_flat_falling_end_is_finite();
  test_almost_exactly_flat_end_is_finite();
  test_equal_end_values_continue_as_a_constant();
  test_continuous_through_equal_end_values();

  if (failures > 0) {
    std::printf("nonlinear_extrapolation_test: %d failure(s)\n", failures);
  }
  assert(failures == 0);
  std::printf("nonlinear_extrapolation_test: all passed\n");
  return 0;
}
