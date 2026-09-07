#include "evolver_options.h"

#include <algorithm>

#include "errors.h"
#include "nonfinite.h"

namespace {

bool Honours(std::initializer_list<EvolverFeature> honoured, EvolverFeature feature) {
  return std::find(honoured.begin(), honoured.end(), feature) != honoured.end();
}

}  // namespace

void EvolverOptionsCheck(const EvolverOptions& options,
                         const char* evolver_name,
                         double x_ini,
                         std::initializer_list<EvolverFeature> honoured) {
  /* Required of every evolver. x_sampling is checked here rather than left to
     each evolver because ndf15 dereferences it before its first step, so a null
     is a segfault with no diagnostic. */
  class_test_severe(options.x_sampling == nullptr || options.x_sampling_size < 2,
                    "%s: x_sampling must be non-null with at least two points. Pass the "
                    "two endpoints if no dense output is wanted.",
                    evolver_name);

  /* Every evolver opens with `while (x_sampling[i] < x_ini) ++i` and no bound
     check, so a grid lying entirely below the start walks off the end of the
     array. Requiring two points does NOT prevent that -- both could be below --
     so reachability is checked directly. */
  class_test_severe(IsNonFinite(x_ini) || options.x_sampling[options.x_sampling_size - 1] < x_ini,
                    "%s: the whole x_sampling grid lies below x_ini = %g (last point %g). "
                    "Every evolver scans forward from the start of that array without a "
                    "bound check, so this would run off the end.",
                    evolver_name,
                    x_ini,
                    options.x_sampling[options.x_sampling_size - 1]);

  /* Called at least once, because x_sampling is non-empty, and invoked without a
     null check by every evolver. */
  class_test_severe(options.output == nullptr,
                    "%s: output must be non-null; it is called at every x_sampling point "
                    "without a null check. Pass a no-op if the samples are not wanted.",
                    evolver_name);
  class_test_severe(options.used_in_output == nullptr,
                    "%s: used_in_output must be non-null, one entry per equation. The "
                    "explicit evolvers index it without checking. Pass all-ones to mean "
                    "every component.",
                    evolver_name);
  /* IsNonFinite, not std::isnan or !(x > 0), and not because of taste. This
     translation unit is part of classpp, which is built with -ffast-math, i.e.
     -ffinite-math-only: the compiler may ASSUME no NaN occurs and fold any test
     for one. `!(x > 0.)` rejects NaN under Apple clang and does NOT under the
     GCC used on Linux CI -- which is how a NaN timestep reached the evolver on
     one platform and not the other. IsNonFinite is declared in nonfinite.h and
     defined in the single TU CMake compiles without the flag; see the comment
     there, which records that no in-TU spelling survives, the bit-pattern form
     included. */
  class_test_severe(IsNonFinite(options.rtol) || !(options.rtol > 0.),
                    "%s: rtol must be positive and finite, got %g",
                    evolver_name,
                    options.rtol);
  class_test_severe(IsNonFinite(options.abstol) || !(options.abstol > 0.),
                    "%s: abstol must be positive and finite, got %g",
                    evolver_name,
                    options.abstol);

  /* And the point of the struct: an option this evolver does not consume is an
     error. Accepting it and doing nothing is how a knob silently stops working
     (#395), which is far harder to notice than a message at startup. */
  const auto reject = [&](EvolverFeature feature, bool is_set, const char* what) {
    class_test_severe(is_set && !Honours(honoured, feature),
                      "%s does not honour %s; setting it would silently do nothing. "
                      "Either drop it or choose an evolver that consumes it.",
                      evolver_name,
                      what);
  };

  reject(EvolverFeature::Timescale,
         options.evaluate_timescale != nullptr || options.timestep_over_timescale != 0.,
         "evaluate_timescale / timestep_over_timescale");
  reject(EvolverFeature::Diagonal, options.derivs_diagonal != nullptr, "derivs_diagonal");
  reject(EvolverFeature::MaxOrder, options.max_order != 0, "max_order");
  reject(EvolverFeature::AbsTol, options.abstol != kEvolverDefaultAbsTol, "a non-default abstol");
  reject(EvolverFeature::ErkController, options.erk.has_value(), "erk");
  reject(EvolverFeature::Stats, options.stats != nullptr, "stats");
  reject(EvolverFeature::Histograms, options.histograms != nullptr, "histograms");

  /* An evolver that honours the timescale pair does not merely accept it, it
     REQUIRES it: the legacy rk takes its whole step size from there and calls
     the function without a null check. So both must be SET, not merely set
     together. */
  if (Honours(honoured, EvolverFeature::Timescale)) {
    class_test_severe(options.evaluate_timescale == nullptr,
                      "%s takes its step size from evaluate_timescale, so it is required.",
                      evolver_name);
    /* A negative factor would send the step backwards while the loop condition
       still asks for x1 < x_end, i.e. away from the endpoint, forever. NaN needs
       IsNonFinite rather than the comparison -- see the note on rtol above. */
    class_test_severe(IsNonFinite(options.timestep_over_timescale) ||
                          !(options.timestep_over_timescale > 0.),
                      "%s needs a positive, finite timestep_over_timescale, got %g",
                      evolver_name,
                      options.timestep_over_timescale);
  }

  /* Order is used to index five-element coefficient arrays. */
  if (Honours(honoured, EvolverFeature::MaxOrder)) {
    class_test_severe(options.max_order < 0 || options.max_order > 5,
                      "%s: max_order must be 0 (the evolver's own default) or 1..5, got %d",
                      evolver_name,
                      options.max_order);
  }
}
