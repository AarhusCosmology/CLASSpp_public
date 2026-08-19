#ifndef __NONFINITE__
#define __NONFINITE__

/**
 * Non-finite (NaN or +/-Inf) test that survives -ffast-math.
 *
 * DECLARED HERE, DEFINED OUT OF LINE in tools/nonfinite.cpp, which CMake compiles
 * without fast-math. That separation IS the mechanism and must not be "simplified"
 * into an inline or header-only helper.
 *
 * classpp is built with -ffast-math, which implies -ffinite-math-only: the compiler
 * may ASSUME no NaN or Inf occurs. That is a licence to optimise, not a guarantee
 * about the arithmetic -- a decaying NCDM run at high Gamma genuinely reaches
 * errmax = inf with a NaN state vector. std::isnan folds to false under that
 * assumption, and -- measured, not assumed -- so does the IEEE-754 bit test through
 * memcpy, because the optimiser propagates "this value is finite" through the bitcast
 * and folds the integer comparison too. No check written in a fast-math TU can work,
 * however it is spelled. There is no LTO in this build, so an out-of-line call in a
 * TU compiled without the flag is opaque at the call site.
 *
 * Denormals and zero have a zero exponent field and read as FINITE, which is what
 * callers want: a stalled run sits at h ~ 1e-303 and a check that rejected denormals
 * would fire on healthy-but-tiny states.
 */
bool IsNonFinite(double x);

/**
 * True when advancing t by h cannot change t, i.e. `t + h == t` in IEEE-754.
 *
 * SAME REASON FOR LIVING HERE as IsNonFinite: it must be compiled WITHOUT
 * -ffast-math. `t + h == t` is algebraically `h == 0`, and -ffast-math licenses
 * exactly that substitution -- at which point the test stops catching the case it
 * exists for (h small but nonzero) and only catches h == 0, which never happens
 * because the step controller floors h. Do not move this back inline.
 *
 * WHY A RELATIVE TEST IS THE RIGHT ONE. ndf15's own floor is the ABSOLUTE
 * minimum_variation (= DBL_EPSILON) and rkdp45's is 100*DBL_MIN*|t|. Neither
 * expresses "h is too small to move t": at t = 6.4e3 that needs h > eps*|t|, some
 * four orders above the absolute floor, so a run can sit forever at one t with h
 * comfortably above every declared floor.
 */
bool StepMakesNoProgress(double t, double h);

#endif
