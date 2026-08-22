#ifndef __EVOERK__
#define __EVOERK__

/**
 * Shared core for the explicit embedded Runge-Kutta evolvers.
 *
 * evolver_rkdp45 (Dormand-Prince 5(4)) and evolver_tsit5 (Tsitouras 5(4)) differ
 * ONLY in their Butcher tableau and dense-output coefficients. Both run through
 * the one driver in evolver_erk_impl.h, so any comparison between them measures
 * the tableau and nothing else -- step controller, error norm, non-finite
 * handling, output emission and FSAL bookkeeping are literally the same code.
 *
 * The tableaux below are POLICY TYPES carrying constexpr coefficient arrays, and
 * the driver is a template over them, so every Butcher coefficient is a
 * compile-time constant inside the stage loops. That is not decoration: with the
 * coefficients behind a runtime pointer the compiler contracts the
 * multiply-adds differently, and the resulting ULP drift in the dense output is
 * amplified by the ODE into a ~2e-5 relative shift in C_l. Keeping them
 * constexpr reproduces the previous hand-written rkdp45 bit for bit.
 *
 * Layout conventions:
 *  - `A` is row-major stages*stages, strict lower triangle used; rest reads 0.
 *  - `b` are the propagating weights (the higher-order solution).
 *  - `e` are the ERROR weights, already differenced: e = b - bhat.
 *  - `R` is row-major stages*dense_deg holding the continuous extension
 *        b_i(theta) = sum_{m=1..dense_deg} R[i*dense_deg + (m-1)] * theta^m,
 *    so that y(t + theta*h) = y_n + h * sum_i b_i(theta) * k_i. b_i(0) = 0 is
 *    built into the parametrisation (no theta^0 column).
 *  - `fsal` means the last stage equals the propagating weights, i.e. row
 *    stages-1 of A equals b, so k_{s-1} of this step is k_0 of the next and a
 *    step costs stages-1 derivative evaluations rather than stages.
 */

/** Dormand-Prince 5(4) -- the historical CLASS `rkdp45`. */
struct ErkDormandPrince45 {
  static constexpr const char* kName = "rkdp45";
  static constexpr int kStages       = 7;
  static constexpr int kOrder        = 5;
  static constexpr int kDenseDeg     = 4;
  static constexpr bool kFsal        = true;

  // clang-format off
  static constexpr double c[7] = {0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0};

  static constexpr double A[49] = {
      /* row 0 */ 0., 0., 0., 0., 0., 0., 0.,
      /* row 1 */ 0.2, 0., 0., 0., 0., 0., 0.,
      /* row 2 */ 3.0 / 40.0, 9.0 / 40.0, 0., 0., 0., 0., 0.,
      /* row 3 */ 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0., 0., 0., 0.,
      /* row 4 */ 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0., 0.,
      0.,
      /* row 5 */ 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0,
      -5103.0 / 18656.0, 0., 0.,
      /* row 6 */ 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0,
      0.};

  static constexpr double b[7] =
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};

  static constexpr double e[7] = {71.0 / 57600.0,
                                  0.0,
                                  -71.0 / 16695.0,
                                  71.0 / 1920.0,
                                  -17253.0 / 339200.0,
                                  22.0 / 525.0,
                                  -1.0 / 40.0};

  /* Shampine's free 4th-order continuous extension (the one MATLAB's ode45
     uses), rewritten from the `i01..i03` + `ixx` spelling the hand-written
     rkdp45 carried. Same numbers. */
  static constexpr double R[28] = {
      /* i=0 */ 1.0, -183.0 / 64.0, 37.0 / 12.0, -145.0 / 128.0,
      /* i=1 */ 0.0, 0.0, 0.0, 0.0,
      /* i=2 */ 0.0, 1500.0 / 371.0, -1000.0 / 159.0, 1000.0 / 371.0,
      /* i=3 */ 0.0, -125.0 / 32.0, 125.0 / 12.0, -375.0 / 64.0,
      /* i=4 */ 0.0, 9477.0 / 3392.0, -729.0 / 106.0, 25515.0 / 6784.0,
      /* i=5 */ 0.0, -11.0 / 7.0, 11.0 / 3.0, -55.0 / 28.0,
      /* i=6 */ 0.0, 1.5, -4.0, 2.5};
  // clang-format on
};

/**
 * Tsitouras 5(4) (Ch. Tsitouras, Comput. Math. Appl. 62 (2011) 770-775).
 *
 * Same shape as Dormand-Prince: 7 stages, FSAL, free 4th-order dense output.
 * The pair minimises the principal 5th-order truncation-error norm subject to
 * only the FIRST-column simplifying assumption, where Dormand-Prince imposes the
 * full row-simplifying assumption; the extra freedom buys a smaller error
 * constant at the same stage count. It is the default non-stiff integrator in
 * OrdinaryDiffEq.jl, which is where its reputation comes from.
 */
struct ErkTsitouras54 {
  static constexpr const char* kName = "tsit5";
  static constexpr int kStages       = 7;
  static constexpr int kOrder        = 5;
  static constexpr int kDenseDeg     = 4;
  static constexpr bool kFsal        = true;

  // clang-format off
  static constexpr double c[7] = {0.0, 0.161, 0.327, 0.9, 0.9800255409045097, 1.0, 1.0};

  static constexpr double A[49] = {
      /* row 0 */ 0., 0., 0., 0., 0., 0., 0.,
      /* row 1 */ 0.161, 0., 0., 0., 0., 0., 0.,
      /* row 2 */ -0.008480655492356989, 0.335480655492357, 0., 0., 0., 0., 0.,
      /* row 3 */ 2.8971530571054935, -6.359448489975075, 4.3622954328695815, 0., 0., 0., 0.,
      /* row 4 */ 5.325864828439257, -11.748883564062828, 7.4955393428898365,
      -0.09249506636175525, 0., 0., 0.,
      /* row 5 */ 5.86145544294642, -12.92096931784711, 8.159367898576159, -0.071584973281401,
      -0.028269050394068383, 0., 0.,
      /* row 6 */ 0.09646076681806523, 0.01, 0.4798896504144996, 1.379008574103742,
      -3.290069515436081, 2.324710524099774, 0.};

  static constexpr double b[7] = {0.09646076681806523,
                                  0.01,
                                  0.4798896504144996,
                                  1.379008574103742,
                                  -3.290069515436081,
                                  2.324710524099774,
                                  0.0};

  static constexpr double e[7] = {-0.001780011052226,
                                  -0.000816434459657,
                                  0.007880878010262,
                                  -0.144711007173263,
                                  0.582357165452555,
                                  -0.458082105929187,
                                  0.015151515151515};

  static constexpr double R[28] = {
      /* i=0 */ 1.0, -2.763706197274826, 2.9132554618219126, -1.0530884977290216,
      /* i=1 */ 0.0, 0.13169999999999998, -0.2234, 0.1017,
      /* i=2 */ 0.0, 3.9302962368947516, -5.941033872131505, 2.490627285651253,
      /* i=3 */ 0.0, -12.411077166933676, 30.33818863028232, -16.548102889244902,
      /* i=4 */ 0.0, 37.50931341651104, -88.1789048947664, 47.37952196281928,
      /* i=5 */ 0.0, -27.896526289197286, 65.09189467479366, -34.8706578614966,
      /* i=6 */ 0.0, 1.5, -4.0, 2.5};
  // clang-format on
};

/**
 * Step-size controller shape. The historical CLASS behaviour is
 * {legacy, max_norm, safety 0.8, no growth cap, growth allowed right after a
 * rejection}; it is kept because every committed reference output used it.
 */
enum class ErkControllerKind {
  legacy, /* uncapped I controller that also grows h right after a rejection */
  capped, /* I controller with a growth cap and no growth right after a rejection */
  pi      /* Hairer/Gustafsson PI controller, capped */
};

enum class ErkErrorNorm {
  max_norm, /* infinity norm of the scaled error (MATLAB ode45, historical CLASS) */
  rms_norm  /* root-mean-square, as in Hairer's DOPRI5 */
};

/**
 * How output points are served.
 *
 *  dense  -- the step sequence ignores the output grid and every point inside a
 *            step is produced by the continuous extension. That extension is one
 *            order BELOW the step (4 against 5), so the sources CLASS actually
 *            stores are a full order less accurate than the solution the error
 *            controller is protecting.
 *  exact  -- h is truncated so that every output point is a step end. The
 *            interpolant is then never used and every stored source carries the
 *            5th-order solution, at the price of a step sequence dictated partly
 *            by the output grid rather than entirely by the error.
 *
 * Which one wins is not decidable from the orders alone: it depends on whether
 * the delivered error is dominated by the interpolation at the sampling point or
 * by the error accumulated up to it, and that differs between observables.
 */
enum class ErkOutputStepping { dense, exact };

/**
 * Which polynomial serves an output point that falls strictly inside a step.
 *
 *  tableau  -- the pair's own continuous extension over the seven stages. Free,
 *              and 4th order: ONE ORDER BELOW the step it interpolates.
 *  hermite3 -- a quintic Hermite through the last THREE step ends, using the
 *              value and the derivative at each. All six data are 5th-order
 *              accurate (the derivative at a step end is f(t, y) of an accepted
 *              solution, and FSAL means it is already in hand), so the
 *              interpolation error is O(h^6) instead of O(h^5) -- at the cost of
 *              two stored vectors and no extra derivative evaluations at all.
 *              The first step of each interval has no history and falls back to
 *              the tableau extension.
 *
 * This matters more than it sounds: measured on LCDM + Mnu, removing the
 * interpolation entirely (ErkOutputStepping::exact) improves linear P(k) by a
 * factor of 41 at the default tolerance. Almost all of the P(k) error CLASS
 * carries is interpolation error at the sampling points, not error in the
 * solution being sampled.
 */
enum class ErkInterpolant { tableau, hermite3 };

struct ErkControllerConfig {
  ErkControllerKind kind            = ErkControllerKind::legacy;
  ErkErrorNorm norm                 = ErkErrorNorm::max_norm;
  ErkOutputStepping output_stepping = ErkOutputStepping::dense;
  ErkInterpolant interpolant        = ErkInterpolant::tableau;
  double safety                     = 0.8;
  double fac_min                    = 0.1;  /* most h may shrink on an ACCEPTED step */
  double fac_max                    = 5.0;  /* most h may grow (ignored by `legacy`) */
  double pi_beta                    = 0.04; /* PI memory exponent; Hairer's DOPRI5 default */
};

/** Process-wide controller settings. Set once from input, read by every thread. */
void evolver_erk_configure(const ErkControllerConfig& config);
const ErkControllerConfig& evolver_erk_config();

/** Aggregate step statistics, accumulated across threads when enabled. */
struct ErkStats {
  long long steps_accepted = 0;
  long long steps_rejected = 0;
  long long derivs_calls   = 0;
  long long dense_points   = 0; /* output points served by the continuous extension */
  long long exact_points   = 0; /* output points that landed exactly on a step end */
};

void evolver_erk_stats_enable(bool on);
bool evolver_erk_stats_enabled();
void evolver_erk_stats_reset();
ErkStats evolver_erk_stats_get();

/**
 * Where the steps and the rejections actually happen.
 *
 * A rejection rate is a single number and does not say whether the controller
 * is overshooting everywhere or only at a few features. Two histograms answer
 * that: one over the accepted/rejected step's error ratio (errnorm / rtol,
 * which is 1 at the acceptance boundary) and one over the independent variable,
 * which for the perturbations is conformal time.
 */
struct ErkHistograms {
  static constexpr int kErrBins  = 32; /* log10(errnorm/rtol) in [-6, +2), 0.25/bin */
  static constexpr double kErrLo = -6.0, kErrStep = 0.25;
  static constexpr int kXBins  = 32; /* log10(x) in [-2, +6), 0.25/bin */
  static constexpr double kXLo = -2.0, kXStep = 0.25;
  long long err_accepted[kErrBins] = {0};
  long long err_rejected[kErrBins] = {0};
  long long x_accepted[kXBins]     = {0};
  long long x_rejected[kXBins]     = {0};
};

/* Histogram accumulation is a SEPARATE switch from the counters: it costs four
   more atomic increments per step, which is enough to move a wall-time
   measurement, and the counters are wanted in every benchmark while the
   histograms are wanted only when diagnosing the controller. */
void evolver_erk_histograms_enable(bool on);
ErkHistograms evolver_erk_histograms_get();

/* Counter hooks used by the templated driver; out of line so the header stays
   free of the atomics. */
namespace erk_detail {
void CountDerivs();
void CountAccepted(double x, double err_ratio);
void CountRejected(double x, double err_ratio);
void CountDense();
void CountExact();
}  // namespace erk_detail

#endif
