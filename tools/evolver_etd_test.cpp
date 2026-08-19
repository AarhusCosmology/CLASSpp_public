#include "evolver_etd.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

struct Problem {
  double lambda = 0.;
  double gain   = 0.;
  int derivs    = 0;
  int diagonals = 0;
  std::vector<double> sampled_y;
  std::vector<double> sampled_dy;
};

void linear_derivs(double /*x*/, double* y, double* dy, void* p) {
  auto* problem = static_cast<Problem*>(p);
  problem->derivs++;
  dy[0] = problem->lambda * y[0] + problem->gain;
}

void linear_diag(double /*x*/, double* /*y*/, double* diag, void* p) {
  auto* problem = static_cast<Problem*>(p);
  problem->diagonals++;
  diag[0] = problem->lambda;
}

void quadratic_derivs(double /*x*/, double* y, double* dy, void* p) {
  auto* problem = static_cast<Problem*>(p);
  problem->derivs++;
  dy[0] = y[0] * y[0];
}

void output_store(double /*x*/, double y[], double dy[], int index_x, void* p) {
  auto* problem = static_cast<Problem*>(p);
  if (index_x >= static_cast<int>(problem->sampled_y.size())) {
    problem->sampled_y.resize(index_x + 1);
    problem->sampled_dy.resize(index_x + 1);
  }
  problem->sampled_y[index_x]  = y[0];
  problem->sampled_dy[index_x] = dy[0];
}

double linear_exact(double y0, double lambda, double gain, double x) {
  const double ystar = -gain / lambda;
  return ystar + (y0 - ystar) * std::exp(lambda * x);
}

void test_stiff_constant_remainder_exact() {
  Problem problem;
  problem.lambda = -80.;
  problem.gain   = 3.;
  problem.sampled_y.resize(5);
  problem.sampled_dy.resize(5);

  double y               = 0.2;
  int used               = 1;
  std::vector<double> xs = {0., 0.25, 0.5, 0.75, 1.};

  evolver_etd(linear_derivs,
              0.,
              1.,
              &y,
              &used,
              1,
              &problem,
              1e-8,
              1e-10,
              nullptr,
              0.,
              xs.data(),
              static_cast<int>(xs.size()),
              output_store,
              nullptr,
              linear_diag);

  for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
    const double expected = linear_exact(0.2, problem.lambda, problem.gain, xs[i]);
    assert(std::fabs(problem.sampled_y[i] - expected) < 2e-12);
    /* dy comes out of the extension as L*y + N, so it is exact here too. */
    const double dexpected = problem.lambda * expected + problem.gain;
    assert(std::fabs(problem.sampled_dy[i] - dexpected) < 2e-10);
  }
  assert(std::fabs(y - linear_exact(0.2, problem.lambda, problem.gain, 1.)) < 2e-12);
  assert(problem.derivs < 40);
  std::printf("etd: stiff constant-remainder solve exact with dense samples (%d RHS, %d diag)\n",
              problem.derivs,
              problem.diagonals);
}

void test_zero_diagonal_rk4_fallback() {
  Problem problem;
  problem.sampled_y.resize(3);
  problem.sampled_dy.resize(3);

  double y               = 1.;
  int used               = 1;
  std::vector<double> xs = {0., 0.25, 0.5};

  evolver_etd(quadratic_derivs,
              0.,
              0.5,
              &y,
              &used,
              1,
              &problem,
              1e-8,
              1e-10,
              nullptr,
              0.,
              xs.data(),
              static_cast<int>(xs.size()),
              output_store,
              nullptr,
              nullptr);

  assert(std::fabs(y - 2.) < 2e-5);
  assert(std::fabs(problem.sampled_y[1] - 4. / 3.) < 2e-5);
  assert(std::fabs(problem.sampled_y[2] - 2.) < 2e-5);
  std::printf("etd: zero-diagonal RK4 fallback accurate (%d RHS)\n", problem.derivs);
}

/* ---------------------------------------------------------------- layer problem --
 *
 *   y' = -lam(t) (y - g(t)) + g'(t),   lam(t) = lam0 (1 + t),   g(t) = sin(w t)
 *   y(t) = g(t) + d0 exp(-lam0 (t + t^2/2))
 *
 * Stiff, with a diagonal that VARIES IN TIME (the constant-lambda tests above cannot
 * see a re-linearisation bug), and with a closed-form solution. The step is set by
 * resolving g, not by lam0, so the accepted step is far wider than the relaxation
 * layer of width ~1/lam0 -- which is exactly the regime the exponential scheme exists
 * for, and exactly the regime where a polynomial dense output fails.
 */
struct Layer {
  static constexpr double lam0 = 500.0;
  static constexpr double w    = 5.0;
  static constexpr double d0   = 1.0;
  int derivs                   = 0;
  std::vector<double> sampled_y;

  static double g(double t) {
    return std::sin(w * t);
  }
  static double gp(double t) {
    return w * std::cos(w * t);
  }
  static double lam(double t) {
    return lam0 * (1.0 + t);
  }
  static double exact(double t) {
    return g(t) + d0 * std::exp(-lam0 * (t + 0.5 * t * t));
  }
};

void layer_derivs(double t, double* y, double* dy, void* p) {
  static_cast<Layer*>(p)->derivs++;
  dy[0] = -Layer::lam(t) * (y[0] - Layer::g(t)) + Layer::gp(t);
}

void layer_diag(double t, double* /*y*/, double* diag, void* /*p*/) {
  diag[0] = -Layer::lam(t);
}

void layer_store(double /*x*/, double y[], double /*dy*/[], int index_x, void* p) {
  auto* l = static_cast<Layer*>(p);
  if (index_x >= static_cast<int>(l->sampled_y.size()))
    l->sampled_y.resize(index_x + 1);
  l->sampled_y[index_x] = y[0];
}

/* Returns the worst sampled error; fills `rhs` with the RHS count. */
double solve_layer(double rtol, int nsample, int* rhs) {
  Layer l;
  std::vector<double> xs(nsample);
  for (int i = 0; i < nsample; ++i)
    xs[i] = 1.0 * i / (nsample - 1);
  l.sampled_y.resize(nsample);

  double y = Layer::exact(0.);
  int used = 1;
  evolver_etd(layer_derivs,
              0.,
              1.,
              &y,
              &used,
              1,
              &l,
              rtol,
              0.,
              nullptr,
              0.,
              xs.data(),
              nsample,
              layer_store,
              nullptr,
              layer_diag);

  double worst = 0.;
  for (int i = 0; i < nsample; ++i)
    worst = std::max(worst, std::fabs(l.sampled_y[i] - Layer::exact(xs[i])));
  *rhs = l.derivs;
  return worst;
}

void test_dense_output_resolves_the_layer() {
  int rhs_coarse = 0, rhs_fine = 0;
  const double coarse = solve_layer(1e-6, 11, &rhs_coarse);
  const double fine   = solve_layer(1e-6, 1001, &rhs_fine);

  /* The bound that matters. A cubic Hermite through the step's endpoints scores ~6e-4
     on the fine grid here -- more than an order of magnitude over this threshold --
     because it interpolates straight across a layer the step jumped. The exponential
     extension carries the layer in e^{theta z} and stays with the trajectory. */
  assert(fine < 5e-5);
  assert(coarse < 5e-5);

  /* Dense output costs no RHS evaluations, so refining the output grid by 100x must
     not change the count at all. This is the entire cost argument for the extension:
     the sampled table is free, whatever its density. */
  assert(rhs_coarse == rhs_fine);

  std::printf(
      "etd: time-varying-diagonal layer resolved (worst %.2e at 11 samples, "
      "%.2e at 1001, %d RHS either way)\n",
      coarse,
      fine,
      rhs_fine);
}

void test_error_estimate_tracks_tolerance() {
  /* The embedded third-order companion has to actually control the step: tightening
     the tolerance must buy accuracy at close to the expected rate. Step doubling used
     to supply this estimate at ~11 RHS per attempted step; the companion costs none. */
  int rhs            = 0;
  const double loose = solve_layer(1e-5, 201, &rhs);
  int rhs_mid        = 0;
  const double mid   = solve_layer(1e-7, 201, &rhs_mid);
  int rhs_tight      = 0;
  const double tight = solve_layer(1e-9, 201, &rhs_tight);

  assert(loose > mid && mid > tight);
  assert(loose < 1e-3);
  assert(mid < 1e-5);
  assert(tight < 1e-7);
  /* Tightening by 1e-4 must not cost anything like 1e-4 in work; an order-4 accepted
     step with an O(h^4) estimate predicts ~10x, and the guard is deliberately loose. */
  assert(rhs_tight < 40 * rhs);
  std::printf("etd: embedded estimate controls the step (%.2e/%d RHS, %.2e/%d, %.2e/%d)\n",
              loose,
              rhs,
              mid,
              rhs_mid,
              tight,
              rhs_tight);
}

/* ------------------------------------------------------- nonlinear coupled system --
 *
 * Two components: one stiff and relaxing onto a state-dependent target, one slow and
 * driven by the first. The remainder N is genuinely nonlinear in y, so the constant-N
 * exactness that carries the tests above is unavailable and the order of the scheme
 * is what has to do the work. No closed form -- a tight run is the reference.
 */
struct Coupled {
  int derivs = 0;
};

void coupled_derivs(double t, double* y, double* dy, void* p) {
  static_cast<Coupled*>(p)->derivs++;
  dy[0] = -200. * (y[0] - y[1] * y[1]) + std::cos(t);
  dy[1] = -y[0] * y[1] + 0.5;
}

void coupled_diag(double /*t*/, double* y, double* diag, void* /*p*/) {
  diag[0] = -200.;
  diag[1] = -y[0];
}

void coupled_none(double /*x*/, double /*y*/[], double /*dy*/[], int /*i*/, void* /*p*/) {}

void test_nonlinear_coupled_converges() {
  auto solve = [](double rtol, double* out) {
    Coupled c;
    double y[2]            = {0.3, 0.8};
    int used[2]            = {1, 1};
    std::vector<double> xs = {0., 1., 2.};
    evolver_etd(coupled_derivs,
                0.,
                2.,
                y,
                used,
                2,
                &c,
                rtol,
                0.,
                nullptr,
                0.,
                xs.data(),
                3,
                coupled_none,
                nullptr,
                coupled_diag);
    out[0] = y[0];
    out[1] = y[1];
    return c.derivs;
  };

  double ref[2], loose[2], tighter[2];
  solve(1e-12, ref);
  const int rhs_loose = solve(1e-6, loose);
  solve(1e-8, tighter);

  const double e_loose   = std::max(std::fabs(loose[0] - ref[0]), std::fabs(loose[1] - ref[1]));
  const double e_tighter = std::max(std::fabs(tighter[0] - ref[0]), std::fabs(tighter[1] - ref[1]));

  assert(e_loose < 1e-5);
  assert(e_tighter < e_loose);
  assert(e_tighter < 1e-7);
  std::printf("etd: nonlinear coupled system converges (%.2e at 1e-6 / %d RHS, %.2e at 1e-8)\n",
              e_loose,
              rhs_loose,
              e_tighter);
}

}  // namespace

int main() {
  test_stiff_constant_remainder_exact();
  test_zero_diagonal_rk4_fallback();
  test_dense_output_resolves_the_layer();
  test_error_estimate_tracks_tolerance();
  test_nonlinear_coupled_converges();
  std::printf("evolver etd test passed\n");
  return 0;
}
