// The explicit-RK pairs, pinned at the level where they can silently go wrong.
//
// A mistyped Butcher or interpolation coefficient does not crash and does not
// obviously change an answer: it QUIETLY drops the order, and the only symptom
// is that the integrator becomes expensive at tight tolerances -- exactly the
// measurement one would then be making. So the tests here are convergence-order
// tests on the coefficients themselves, not smoke tests on the output.
//
// Reference behaviour for both pairs:
//   propagated solution   local error  O(h^6)   (order 5)
//   embedded difference   O(h^5)                (order 4, so the estimate is
//                                                asymptotically the true error
//                                                of the 4th-order solution)
//   dense output          O(h^5)                (order 4 continuous extension)

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "evolver_erk_impl.h"

namespace {

/* --- order conditions -------------------------------------------------

   A mistyped coefficient is caught here rather than by a convergence
   experiment, because the residual of an order condition is exact to roundoff
   while an empirical order fit depends on landing inside the asymptotic window
   -- for these pairs that window is only about one octave wide before the local
   error reaches 1e-16, which makes the fit a worse test than the algebra.

   The elementary weights below are the standard rooted-tree conditions through
   order 5, in the vector form b^T PHI(t) = 1/gamma(t).                     */

using Vec = std::vector<double>;

Vec apply_A(const double* A, int s, const Vec& v) {
  Vec out(s, 0.0);
  for (int i = 0; i < s; i++)
    for (int j = 0; j < s; j++)
      out[i] += A[i * s + j] * v[j];
  return out;
}
Vec mul(const Vec& a, const Vec& b) {
  Vec out(a.size());
  for (size_t i = 0; i < a.size(); i++)
    out[i] = a[i] * b[i];
  return out;
}
double dot(const Vec& w, const Vec& v) {
  double acc = 0.0;
  for (size_t i = 0; i < w.size(); i++)
    acc += w[i] * v[i];
  return acc;
}

struct Condition {
  const char* name;
  int order;
  Vec phi;      /* elementary weight vector */
  double gamma; /* the condition is b.phi == 1/gamma */
};

/* The 17 trees of order <= 5, built from c and A. */
template <class Tab>
std::vector<Condition> conditions() {
  constexpr int s = Tab::kStages;
  const double* A = Tab::A;
  Vec one(s, 1.0), c(Tab::c, Tab::c + s);
  Vec c2 = mul(c, c), c3 = mul(c2, c), c4 = mul(c3, c);
  Vec Ac = apply_A(A, s, c), Ac2 = apply_A(A, s, c2), Ac3 = apply_A(A, s, c3);
  Vec AAc = apply_A(A, s, Ac), AcAc = apply_A(A, s, mul(c, Ac)), AAc2 = apply_A(A, s, Ac2);
  Vec AAAc = apply_A(A, s, AAc);
  return {
      {"1", 1, one, 1.0},
      {"c", 2, c, 2.0},
      {"c^2", 3, c2, 3.0},
      {"Ac", 3, Ac, 6.0},
      {"c^3", 4, c3, 4.0},
      {"c(Ac)", 4, mul(c, Ac), 8.0},
      {"Ac^2", 4, Ac2, 12.0},
      {"AAc", 4, AAc, 24.0},
      {"c^4", 5, c4, 5.0},
      {"c^2(Ac)", 5, mul(c2, Ac), 10.0},
      {"(Ac)^2", 5, mul(Ac, Ac), 20.0},
      {"c(Ac^2)", 5, mul(c, Ac2), 15.0},
      {"Ac^3", 5, Ac3, 20.0},
      {"c(AAc)", 5, mul(c, AAc), 30.0},
      {"A(c Ac)", 5, AcAc, 40.0},
      {"AAc^2", 5, AAc2, 60.0},
      {"AAAc", 5, AAAc, 120.0},
  };
}

template <class Tab>
void test_tableau() {
  constexpr int s  = Tab::kStages;
  const auto conds = conditions<Tab>();

  /* Consistency of the stages: row sums of A are the abscissae. */
  double worst_row = 0.0;
  for (int i = 0; i < s; i++) {
    double row = 0.0;
    for (int j = 0; j < s; j++)
      row += Tab::A[i * s + j];
    worst_row = std::max(worst_row, fabs(row - Tab::c[i]));
  }
  assert(worst_row < 1e-14);

  /* FSAL: the last stage must BE the propagating solution, otherwise reusing
     k_{s-1} as the next k_0 silently integrates a different method. */
  if (Tab::kFsal) {
    double worst = 0.0;
    for (int j = 0; j < s; j++)
      worst = std::max(worst, fabs(Tab::A[(s - 1) * s + j] - Tab::b[j]));
    assert(worst < 1e-14);
  }

  Vec b(Tab::b, Tab::b + s), e(Tab::e, Tab::e + s), bhat(s);
  for (int j = 0; j < s; j++)
    bhat[j] = b[j] - e[j];

  double worst_b = 0.0, worst_bhat4 = 0.0, best_bhat5 = 1e300;
  for (const Condition& cond : conds) {
    const double rb = dot(b, cond.phi) - 1.0 / cond.gamma;
    worst_b         = std::max(worst_b, fabs(rb));
    const double rh = dot(bhat, cond.phi) - 1.0 / cond.gamma;
    if (cond.order <= 4)
      worst_bhat4 = std::max(worst_bhat4, fabs(rh));
    else
      best_bhat5 = std::min(best_bhat5, fabs(rh));
  }
  printf(
      "  %-8s order conditions: b through 5 -> %.1e, bhat through 4 -> %.1e, "
      "bhat at 5 -> %.1e (must NOT vanish)\n",
      Tab::kName,
      worst_b,
      worst_bhat4,
      best_bhat5);
  assert(worst_b < 1e-12);     /* the propagated solution is order 5 */
  assert(worst_bhat4 < 1e-12); /* the embedded solution is order 4 */
  /* ... and NOT order 5, or the "error estimate" would be a difference of two
     order-5 solutions and would not measure the local error at all. */
  assert(best_bhat5 > 1e-6);
}

/* --- continuous extension ------------------------------------------- */

/* b_i(theta) for the dense output. */
template <class Tab>
Vec dense_weights(double theta) {
  constexpr int s = Tab::kStages;
  Vec b(s, 0.0);
  for (int i = 0; i < s; i++) {
    double thp = theta;
    for (int m = 0; m < Tab::kDenseDeg; m++) {
      b[i] += Tab::R[i * Tab::kDenseDeg + m] * thp;
      thp  *= theta;
    }
  }
  return b;
}

template <class Tab>
void test_dense() {
  constexpr int s  = Tab::kStages;
  const auto conds = conditions<Tab>();

  /* At theta = 0 the extension must return y_n exactly, and at theta = 1 it
     must return the step the controller actually accepted -- otherwise the
     source functions are sampled off a different trajectory than the one being
     error-controlled. */
  const Vec b0 = dense_weights<Tab>(0.0), b1 = dense_weights<Tab>(1.0);
  double w0 = 0.0, w1 = 0.0;
  for (int i = 0; i < s; i++) {
    w0 = std::max(w0, fabs(b0[i]));
    w1 = std::max(w1, fabs(b1[i] - Tab::b[i]));
  }
  assert(w0 < 1e-14);
  assert(w1 < 1e-13);

  /* Continuous order conditions: sum_i b_i(theta) PHI_i(t) = theta^|t| / gamma(t)
     for every tree of order <= 4. Checked at several theta because they are
     polynomial identities, so a coefficient error that happens to cancel at one
     theta cannot cancel at all of them. */
  double worst4 = 0.0, worst5 = 0.0;
  for (double theta : {0.1, 0.25, 0.5, 0.75, 0.9}) {
    const Vec bt = dense_weights<Tab>(theta);
    double thp   = 1.0;
    for (const Condition& cond : conds) {
      thp             = pow(theta, cond.order);
      const double rr = dot(bt, cond.phi) - thp / cond.gamma;
      if (cond.order <= 4)
        worst4 = std::max(worst4, fabs(rr));
      else
        worst5 = std::max(worst5, fabs(rr));
    }
  }
  printf(
      "  %-8s dense output: continuous conditions through 4 -> %.1e, "
      "at 5 -> %.1e (4th-order extension, so nonzero as expected)\n",
      Tab::kName,
      worst4,
      worst5);
  assert(worst4 < 1e-12);
  /* The extension is one order BELOW the step it interpolates. That is the
     accuracy CLASS actually delivers at a source-sampling point, and it is a
     property worth having pinned rather than discovered later. */
  assert(worst5 > 1e-5);
}

/* A genuinely NONLINEAR system with a closed-form solution, used for the
   end-to-end checks below:
       y0' = -2 t y0^2 ,  y0(0) = 1  =>  y0 = 1/(1+t^2)
       y1' = y0 y1     ,  y1(0) = 1  =>  y1 = exp(atan t)                 */
void nl_derivs(double t, double* y, double* dy, void* p) {
  if (p != nullptr)
    (*static_cast<int*>(p))++;
  dy[0] = -2.0 * t * y[0] * y[0];
  dy[1] = y[0] * y[1];
}

void nl_exact(double t, double* y) {
  y[0] = 1.0 / (1.0 + t * t);
  y[1] = exp(atan(t));
}

/* --- end-to-end through the actual driver ---------------------------- */

struct Sampled {
  std::vector<double> t;
  std::vector<double> y0;
  int derivs = 0;
};

void store(double x, double y[], double /*dy*/[], int index_x, void* p) {
  auto* s        = static_cast<Sampled*>(p);
  s->t[index_x]  = x;
  s->y0[index_x] = y[0];
}

void counting_derivs(double t, double* y, double* dy, void* p) {
  nl_derivs(t, y, dy, &static_cast<Sampled*>(p)->derivs);
}

template <class Tab>
double run_to_tolerance(double rtol, int* rhs_calls) {
  const int n = 41;
  Sampled s;
  s.t.resize(n);
  s.y0.resize(n);
  std::vector<double> xs(n);
  for (int i = 0; i < n; i++)
    xs[i] = 0.05 * i; /* deliberately off the step ends: exercises dense output */

  double y[2] = {1.0, 1.0};
  int used[2] = {1, 1};
  evolver_erk_run<
      Tab>(counting_derivs, 0.0, 2.0, y, used, 2, &s, rtol, xs.data(), n, store, nullptr);

  double worst = 0.0;
  for (int i = 1; i < n; i++) {
    double ex[2];
    nl_exact(s.t[i], ex);
    worst = std::max(worst, fabs(s.y0[i] - ex[0]) / fabs(ex[0]));
  }
  *rhs_calls = s.derivs;
  return worst;
}

template <class Tab>
void test_end_to_end() {
  for (double rtol : {1e-4, 1e-6, 1e-8}) {
    int calls    = 0;
    double worst = run_to_tolerance<Tab>(rtol, &calls);
    printf("  %-8s rtol=%.0e -> worst dense-output relative error %.2e (%d RHS calls)\n",
           Tab::kName,
           rtol,
           worst,
           calls);
    /* The controller bounds the LOCAL error per step; the accumulated error at
       an output point is allowed to exceed rtol, but not by orders of magnitude
       on a problem this benign. */
    assert(worst < 50.0 * rtol);
    assert(calls > 0);
  }
}

/* The two pairs must agree with each other to within the tolerance they were
   asked for. Disagreeing by much more would mean one of them is not solving the
   stated problem, which no single-method test can catch. */
void test_pairs_agree() {
  const int n = 41;
  std::vector<double> xs(n);
  for (int i = 0; i < n; i++)
    xs[i] = 0.05 * i;
  Sampled a, b;
  a.t.resize(n);
  a.y0.resize(n);
  b.t.resize(n);
  b.y0.resize(n);
  double ya[2] = {1.0, 1.0}, yb[2] = {1.0, 1.0};
  int used[2] = {1, 1};
  evolver_erk_run<ErkDormandPrince45>(counting_derivs,
                                      0.0,
                                      2.0,
                                      ya,
                                      used,
                                      2,
                                      &a,
                                      1e-8,
                                      xs.data(),
                                      n,
                                      store,
                                      nullptr);
  evolver_erk_run<ErkTsitouras54>(counting_derivs,
                                  0.0,
                                  2.0,
                                  yb,
                                  used,
                                  2,
                                  &b,
                                  1e-8,
                                  xs.data(),
                                  n,
                                  store,
                                  nullptr);
  double worst = 0.0;
  for (int i = 1; i < n; i++)
    worst = std::max(worst, fabs(a.y0[i] - b.y0[i]) / fabs(a.y0[i]));
  printf("  rkdp45 vs tsit5 at rtol=1e-8: max relative difference %.2e\n", worst);
  assert(worst < 1e-5);
}

/* Exact output stepping must actually deliver what it claims: no output point
   may be served by the interpolant, and the answer must still be right. */
void test_exact_output_stepping() {
  const ErkControllerConfig saved = evolver_erk_config();
  ErkControllerConfig cfg;
  cfg.output_stepping = ErkOutputStepping::exact;
  evolver_erk_configure(cfg);
  evolver_erk_stats_reset();
  evolver_erk_stats_enable(true);
  int calls    = 0;
  double worst = run_to_tolerance<ErkTsitouras54>(1e-8, &calls);
  evolver_erk_stats_enable(false);
  const ErkStats st = evolver_erk_stats_get();
  printf(
      "  exact output stepping: error %.2e, %lld accepted, %lld dense points "
      "(only the interval start), %lld exact points\n",
      worst,
      st.steps_accepted,
      st.dense_points,
      st.exact_points);
  /* One sampling point coincides with the start of the interval. It is served
     through the extension at theta = 0, where b_i(0) = 0 and the result is y
     exactly, so it is not an interpolation in any meaningful sense; everything
     strictly inside a step must be an exact step end. */
  assert(st.dense_points <= 1);
  assert(st.exact_points == 40);
  assert(worst < 50.0 * 1e-8);
  evolver_erk_configure(saved);
}

/* The quintic Hermite interpolant must beat the tableau's own continuous
   extension on the same step sequence -- that is its entire reason to exist.
   Both are compared at the same tolerance so the step sequences are identical
   and the only difference measured is the polynomial. */
void test_hermite_interpolant() {
  const ErkControllerConfig saved = evolver_erk_config();
  double err_tab = 0.0, err_her = 0.0;
  int calls_tab = 0, calls_her = 0;
  {
    ErkControllerConfig cfg;
    cfg.interpolant = ErkInterpolant::tableau;
    evolver_erk_configure(cfg);
    err_tab = run_to_tolerance<ErkTsitouras54>(1e-6, &calls_tab);
  }
  {
    ErkControllerConfig cfg;
    cfg.interpolant = ErkInterpolant::hermite3;
    evolver_erk_configure(cfg);
    err_her = run_to_tolerance<ErkTsitouras54>(1e-6, &calls_her);
  }
  printf(
      "  at rtol=1e-6, same %d RHS calls: tableau extension %.2e, quintic Hermite %.2e "
      "(%.1fx better)\n",
      calls_tab,
      err_tab,
      err_her,
      err_tab / err_her);
  assert(calls_tab == calls_her); /* the interpolant must not change the stepping */
  assert(err_her < err_tab);
  evolver_erk_configure(saved);
}

/* Every controller must reach the requested accuracy. What differs is the COST
   and the fraction of attempts thrown away, which is the whole reason the
   alternatives exist -- so the counters are reported, not asserted on. */
void test_controllers() {
  const ErkControllerConfig saved = evolver_erk_config();
  struct Case {
    const char* name;
    ErkControllerKind kind;
    ErkErrorNorm norm;
  };
  const Case cases[] = {{"legacy", ErkControllerKind::legacy, ErkErrorNorm::max_norm},
                        {"capped", ErkControllerKind::capped, ErkErrorNorm::max_norm},
                        {"pi", ErkControllerKind::pi, ErkErrorNorm::max_norm},
                        {"pi+rms", ErkControllerKind::pi, ErkErrorNorm::rms_norm}};
  for (const Case& c : cases) {
    ErkControllerConfig cfg;
    cfg.kind = c.kind;
    cfg.norm = c.norm;
    evolver_erk_configure(cfg);
    evolver_erk_stats_reset();
    evolver_erk_stats_enable(true);
    int calls    = 0;
    double worst = run_to_tolerance<ErkTsitouras54>(1e-8, &calls);
    evolver_erk_stats_enable(false);
    const ErkStats st  = evolver_erk_stats_get();
    const long long at = st.steps_accepted + st.steps_rejected;
    printf(
        "  controller %-7s error %.2e, %lld accepted / %lld attempted (%.1f%% rejected), "
        "%lld dense points\n",
        c.name,
        worst,
        st.steps_accepted,
        at,
        at > 0 ? 100.0 * st.steps_rejected / at : 0.0,
        st.dense_points);
    assert(worst < 50.0 * 1e-8);
    assert(st.steps_accepted > 0);
    assert(st.derivs_calls == calls);
  }
  evolver_erk_configure(saved);
}

}  // namespace

int main() {
  printf("Butcher tableaux:\n");
  test_tableau<ErkDormandPrince45>();
  test_tableau<ErkTsitouras54>();
  printf("Continuous extensions:\n");
  test_dense<ErkDormandPrince45>();
  test_dense<ErkTsitouras54>();
  printf("End to end through the driver:\n");
  test_end_to_end<ErkDormandPrince45>();
  test_end_to_end<ErkTsitouras54>();
  printf("Cross-checks:\n");
  test_pairs_agree();
  printf("Output stepping:\n");
  test_exact_output_stepping();
  printf("Interpolant:\n");
  test_hermite_interpolant();
  printf("Step controllers:\n");
  test_controllers();
  printf("evolver_erk_test: all assertions passed\n");
  return 0;
}
