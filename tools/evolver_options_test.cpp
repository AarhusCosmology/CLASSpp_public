/** The options struct, and the rule that gives it its point: an option an
 *  evolver does not consume is an ERROR, not an ignore.
 *
 *  Before this, passing `timestep_over_timescale` to ndf15 was accepted and
 *  quietly did nothing, and the only record of which evolver honoured what was
 *  prose in five headers -- prose that had already drifted (evolver_ndf15.h said
 *  minimum_variation was unused while the code assigned it to hmin).
 */

#include "evolver_options.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include "evolver_etd.h"
#include "evolver_ndf15.h"
#include "evolver_rkck.h"
#include "evolver_rkdp45.h"

namespace {

/* y' = -y, y(0) = 1. Stiff enough for nothing, which is the point: these tests
   are about the interface, not the integration. */
void decay_derivs(double, double* y, double* dy, void*) {
  dy[0] = -y[0];
}
void decay_diagonal(double, double*, double* diag, void*) {
  diag[0] = -1.0;
}
void decay_timescale(double, void*, double* timescale) {
  *timescale = 1.0;
}
void no_output(double, double[], double[], int, void*) {}

/* One equation in every test here, so one flag. */
const int kUsedInOutput[1] = {1};

EvolverOptions Baseline(const std::vector<double>& sampling) {
  EvolverOptions options;
  options.rtol            = 1e-8;
  options.x_sampling      = sampling.data();
  options.x_sampling_size = static_cast<int>(sampling.size());
  options.used_in_output  = kUsedInOutput;
  options.output          = no_output;
  return options;
}

bool Rejects(const EvolverOptions& options,
             const char* name,
             std::initializer_list<EvolverFeature> honoured,
             double x_ini = 0.0) {
  try {
    EvolverOptionsCheck(options, name, x_ini, honoured);
  }
  catch (const std::exception&) {
    return true;
  }
  return false;
}

/* Every evolver requires these, so they are checked centrally rather than five
   times. x_sampling especially: ndf15 dereferences it before its first step, so
   a null used to be a segfault with no diagnostic. */
void test_universal_requirements() {
  const std::vector<double> sampling = {0.0, 1.0};

  EvolverOptions null_sampling = Baseline(sampling);
  null_sampling.x_sampling     = nullptr;
  assert(Rejects(null_sampling, "any", {}));

  EvolverOptions empty_sampling  = Baseline(sampling);
  empty_sampling.x_sampling_size = 0;
  assert(Rejects(empty_sampling, "any", {}));

  /* Also required: the explicit evolvers index it without checking. */
  EvolverOptions no_flags = Baseline(sampling);
  no_flags.used_in_output = nullptr;
  assert(Rejects(no_flags, "any", {}));

  /* And output, which every evolver calls without a null check. */
  EvolverOptions no_out = Baseline(sampling);
  no_out.output         = nullptr;
  assert(Rejects(no_out, "any", {}));

  /* One sampling point is not enough -- the contract asks for two endpoints. */
  EvolverOptions one_point  = Baseline(sampling);
  one_point.x_sampling_size = 1;
  assert(Rejects(one_point, "any", {}));

  /* A grid lying entirely below x_ini would send every evolver's opening
     `while (x_sampling[i] < x_ini) ++i` off the end of the array. Two points do
     not save it -- both can be below -- which is why reachability is checked
     rather than merely size. */
  assert(Rejects(Baseline(sampling), "any", {}, 99.0));

  for (double bad : {0.0, -1e-8}) {
    EvolverOptions bad_rtol = Baseline(sampling);
    bad_rtol.rtol           = bad;
    assert(Rejects(bad_rtol, "any", {}));

    EvolverOptions bad_abstol = Baseline(sampling);
    bad_abstol.abstol         = bad;
    assert(Rejects(bad_abstol, "any", {}));
  }

  assert(!Rejects(Baseline(sampling), "any", {}));
}

/* The table that used to live in five prose comments, now executable. */
void test_features_are_rejected_where_unsupported() {
  const std::vector<double> sampling = {0.0, 1.0};

  struct Case {
    const char* name;
    EvolverFeature feature;
    void (*apply)(EvolverOptions&);
  };
  const Case cases[] = {
      {"Timescale",
       EvolverFeature::Timescale,
       [](EvolverOptions& o) {
         o.evaluate_timescale      = decay_timescale;
         o.timestep_over_timescale = 0.1;
       }},
      {"Diagonal",
       EvolverFeature::Diagonal,
       [](EvolverOptions& o) { o.derivs_diagonal = decay_diagonal; }},
      {"MaxOrder", EvolverFeature::MaxOrder, [](EvolverOptions& o) { o.max_order = 2; }},
      {"ErkController",
       EvolverFeature::ErkController,
       [](EvolverOptions& o) { o.erk = ErkControllerConfig{}; }},
  };

  for (const Case& c : cases) {
    EvolverOptions options = Baseline(sampling);
    c.apply(options);
    /* Rejected by an evolver that honours nothing... */
    assert(Rejects(options, "honours-nothing", {}));
    /* ...and accepted by one that honours exactly this. */
    assert(!Rejects(options, "honours-it", {c.feature}));
  }

  /* An evolver that honours the timescale pair REQUIRES it: rk takes its whole
     step size from there and calls the function without a null check. So an
     unset pair fails too, not only a half-set one. */
  assert(Rejects(Baseline(sampling), "rk", {EvolverFeature::Timescale}));

  EvolverOptions half     = Baseline(sampling);
  half.evaluate_timescale = decay_timescale;
  assert(Rejects(half, "rk", {EvolverFeature::Timescale}));

  half.timestep_over_timescale = 0.1;
  assert(!Rejects(half, "rk", {EvolverFeature::Timescale}));
}

/* And the real evolvers must declare what they actually consume. */
void test_real_evolvers_reject_what_they_ignore() {
  const std::vector<double> sampling = {0.0, 1.0};
  double y                           = 1.0;

  const auto throws = [&](const EvolverOptions& options, auto evolver) {
    double state = 1.0;
    try {
      evolver(decay_derivs, 0.0, 1.0, &state, 1, nullptr, options);
    }
    catch (const std::exception&) {
      return true;
    }
    return false;
  };

  EvolverOptions diagonal  = Baseline(sampling);
  diagonal.derivs_diagonal = decay_diagonal;
  assert(throws(diagonal, evolver_ndf15));  /* ndf15 does not exponentiate it */
  assert(throws(diagonal, evolver_rkdp45)); /* nor does rkdp45 */
  assert(!throws(diagonal, evolver_etd));   /* etd is the one that does */

  EvolverOptions order = Baseline(sampling);
  order.max_order      = 2;
  assert(!throws(order, evolver_ndf15)); /* BDF order is an ndf15 notion */
  assert(throws(order, evolver_rkdp45));

  EvolverOptions controller = Baseline(sampling);
  controller.erk            = ErkControllerConfig{};
  assert(throws(controller, evolver_ndf15));
  assert(!throws(controller, evolver_rkdp45));

  (void) y;
}

/* Two capabilities the old interface could not express at all. */
void test_new_capabilities_reach_the_evolver() {
  const std::vector<double> sampling = {0.0, 0.5, 1.0};

  /* Statistics: ndf15 has collected these six counters into a local array and
     discarded them since it was written. */
  EvolverStats stats;
  EvolverOptions options = Baseline(sampling);
  options.stats          = &stats;
  double y               = 1.0;
  evolver_ndf15(decay_derivs, 0.0, 1.0, &y, 1, nullptr, options);

  assert(std::fabs(y - std::exp(-1.0)) < 1e-6);
  assert(stats.steps_accepted > 0);
  assert(stats.derivs_evaluations > 0);
  assert(stats.jacobians > 0);

  /* Maximum order: a lower cap is a real constraint, so it must cost steps. */
  EvolverStats capped_stats;
  EvolverOptions capped = Baseline(sampling);
  capped.stats          = &capped_stats;
  capped.max_order      = 1; /* backward Euler */
  capped.rtol           = 1e-10;
  double y_capped       = 1.0;
  evolver_ndf15(decay_derivs, 0.0, 1.0, &y_capped, 1, nullptr, capped);

  assert(std::fabs(y_capped - std::exp(-1.0)) < 1e-5);
  assert(capped_stats.steps_accepted > stats.steps_accepted);
}

/** Both evolver families report through the SAME struct.
 *
 *  They used to be asymmetric: ndf15 collected six counters into a local array
 *  and discarded them, while the explicit evolvers accumulated into process-wide
 *  atomics behind an enable/reset/get API. Each fills the fields it actually
 *  produces and leaves the rest at zero, which is the only asymmetry left and is
 *  a property of the methods rather than of the interface. */
void test_both_evolver_families_report_the_same_way() {
  const std::vector<double> sampling = {0.0, 0.25, 0.5, 0.75, 1.0};

  EvolverStats implicit_stats;
  EvolverOptions implicit_options = Baseline(sampling);
  implicit_options.stats          = &implicit_stats;
  double y_implicit               = 1.0;
  evolver_ndf15(decay_derivs, 0.0, 1.0, &y_implicit, 1, nullptr, implicit_options);

  EvolverStats explicit_stats;
  EvolverOptions explicit_options = Baseline(sampling);
  explicit_options.stats          = &explicit_stats;
  double y_explicit               = 1.0;
  evolver_rkdp45(decay_derivs, 0.0, 1.0, &y_explicit, 1, nullptr, explicit_options);

  /* Both solved the same problem... */
  assert(std::fabs(y_implicit - std::exp(-1.0)) < 1e-6);
  assert(std::fabs(y_explicit - std::exp(-1.0)) < 1e-6);

  /* ...and both reported the counters every method has. */
  for (const EvolverStats* st : {&implicit_stats, &explicit_stats}) {
    assert(st->steps_accepted > 0);
    assert(st->derivs_evaluations > 0);
  }

  /* Only the implicit method factorises a Jacobian. */
  assert(implicit_stats.jacobians > 0);
  assert(implicit_stats.lu_decompositions > 0);
  assert(explicit_stats.jacobians == 0);

  /* Only the explicit one has a continuous extension to serve output from. The
     sampling grid above deliberately falls off the step ends. */
  assert(explicit_stats.dense_points + explicit_stats.exact_points > 0);
  assert(implicit_stats.dense_points == 0);
}

/** Histograms are explicit-RK only, and cost a log10 per step -- hence a
 *  separate opt-in from the counters, which are free. */
void test_histograms() {
  const std::vector<double> sampling = {0.0, 0.5, 1.0};

  ErkHistograms histograms;
  EvolverOptions options = Baseline(sampling);
  options.histograms     = &histograms;
  double y               = 1.0;
  evolver_rkdp45(decay_derivs, 0.0, 1.0, &y, 1, nullptr, options);

  long long binned = 0;
  for (int i = 0; i < ErkHistograms::kErrBins; ++i) {
    binned += histograms.err_accepted[i] + histograms.err_rejected[i];
  }
  assert(binned > 0);

  /* ndf15 has no such notion, so asking it is an error rather than an empty
     result that looks like "no steps were taken". */
  EvolverOptions wrong = Baseline(sampling);
  wrong.histograms     = &histograms;
  double y_wrong       = 1.0;
  bool threw           = false;
  try {
    evolver_ndf15(decay_derivs, 0.0, 1.0, &y_wrong, 1, nullptr, wrong);
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);
}

/** abstol must actually reach EVERY evolver, not just ndf15.
 *
 *  Each of them hardwired 1e-15 independently, and the header claimed the field
 *  was honoured everywhere while two of them ignored it -- the exact silent
 *  ignore this struct exists to prevent, in the struct's own documentation.
 *  Raising the error-weight floor stops small components being resolved
 *  relatively, so it must cost fewer steps. */
void test_abstol_reaches_every_evolver() {
  /* Integrated far enough that y falls well below the raised floor. */
  const std::vector<double> sampling = {0.0, 30.0};

  const auto steps_at = [&](double abstol, auto evolver) {
    EvolverStats stats;
    EvolverOptions options = Baseline(sampling);
    options.rtol           = 1e-10;
    options.abstol         = abstol;
    options.stats          = &stats;
    double y               = 1.0;
    evolver(decay_derivs, 0.0, 30.0, &y, 1, nullptr, options);
    return stats.steps_accepted;
  };

  for (auto evolver : {evolver_ndf15, evolver_rkdp45}) {
    const long long tight = steps_at(1e-15, evolver);
    const long long loose = steps_at(1e-3, evolver);
    assert(tight > 0 && loose > 0);
    assert(loose < tight);
  }
}

/** Values that pass a bare "is it set" test but are still nonsense. */
void test_out_of_range_values() {
  const std::vector<double> sampling = {0.0, 1.0};

  /* max_order indexes five-element coefficient arrays, so 6 is an out-of-bounds
     read and a negative value would silently mean "default". */
  for (int bad : {-1, 6, 99}) {
    EvolverOptions options = Baseline(sampling);
    options.max_order      = bad;
    assert(Rejects(options, "ndf15", {EvolverFeature::MaxOrder}));
  }
  for (int good : {0, 1, 5}) {
    EvolverOptions options = Baseline(sampling);
    options.max_order      = good;
    assert(!Rejects(options, "ndf15", {EvolverFeature::MaxOrder}));
  }

  /* A negative step factor sends rk backwards while its loop still asks for
     x1 < x_end, i.e. away from the endpoint; NaN passes any `!= 0` test. */
  for (double bad : {-0.1, 0.0, std::nan("")}) {
    EvolverOptions options          = Baseline(sampling);
    options.evaluate_timescale      = decay_timescale;
    options.timestep_over_timescale = bad;
    assert(Rejects(options, "rk", {EvolverFeature::Timescale}));
  }
}

/** Non-finite inputs must be rejected in a build with -ffast-math.
 *
 *  This is not a hypothetical. EvolverOptionsCheck lives in classpp, which is
 *  compiled with -ffast-math (-ffinite-math-only), where the compiler may assume
 *  no NaN occurs and fold any test for one. Written as `!(x > 0.)` the NaN case
 *  passed under Apple clang and FAILED under the GCC on Linux CI -- so a NaN
 *  reached the evolver on one platform and not the other, and the local test
 *  suite was green throughout. The checks now go through IsNonFinite, declared in
 *  nonfinite.h and defined in the one translation unit CMake compiles without the
 *  flag.
 *
 *  This test binary is NOT built with -ffast-math, so the NaN it constructs is a
 *  real one crossing into a fast-math TU -- exactly the CI arrangement. */
void test_non_finite_inputs_are_rejected() {
  const std::vector<double> sampling = {0.0, 1.0};
  const double nan_value             = std::nan("");
  const double inf_value             = std::numeric_limits<double>::infinity();

  for (double bad : {nan_value, inf_value}) {
    EvolverOptions bad_rtol = Baseline(sampling);
    bad_rtol.rtol           = bad;
    assert(Rejects(bad_rtol, "any", {}));

    EvolverOptions bad_abstol = Baseline(sampling);
    bad_abstol.abstol         = bad;
    assert(Rejects(bad_abstol, "any", {EvolverFeature::AbsTol}));

    EvolverOptions bad_step          = Baseline(sampling);
    bad_step.evaluate_timescale      = decay_timescale;
    bad_step.timestep_over_timescale = bad;
    assert(Rejects(bad_step, "rk", {EvolverFeature::Timescale}));
  }

  /* And a non-finite integration start, which the reachability check compares
     against. */
  assert(Rejects(Baseline(sampling), "any", {}, nan_value));
}

/** abstol is feature-gated, not universal.
 *
 *  The legacy rk builds its error scale inside generic_integrator from
 *  _TINY_ = 1e-30, a different quantity that cannot be swapped for abstol
 *  without changing results. So rather than let abstol be silently ignored
 *  there -- the failure this whole struct exists to prevent -- setting a
 *  non-default value with that evolver is an error. */
void test_abstol_is_feature_gated() {
  const std::vector<double> sampling = {0.0, 1.0};

  EvolverOptions options          = Baseline(sampling);
  options.evaluate_timescale      = decay_timescale;
  options.timestep_over_timescale = 0.1;
  /* The default passes everywhere, so existing callers are unaffected. */
  assert(!Rejects(options, "rk", {EvolverFeature::Timescale}));

  options.abstol = 1e-8;
  assert(Rejects(options, "rk", {EvolverFeature::Timescale}));
  assert(!Rejects(options, "ndf15", {EvolverFeature::Timescale, EvolverFeature::AbsTol}));
}

}  // namespace

int main() {
  test_universal_requirements();
  test_features_are_rejected_where_unsupported();
  test_real_evolvers_reject_what_they_ignore();
  test_new_capabilities_reach_the_evolver();
  test_both_evolver_families_report_the_same_way();
  test_histograms();
  test_abstol_reaches_every_evolver();
  test_out_of_range_values();
  test_non_finite_inputs_are_rejected();
  test_abstol_is_feature_gated();
  std::printf("evolver options tests passed\n");
  return 0;
}
