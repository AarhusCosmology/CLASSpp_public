#ifndef CLASS_EVOLVER_OPTIONS_H
#define CLASS_EVOLVER_OPTIONS_H

/** @file evolver_options.h The one place the shared evolver interface is defined.
 *
 *  Design: docs/superpowers/specs/2026-09-06-evolver-options-design.md
 *
 *  Every evolver takes (derivs, range, state, workspace, options). Everything
 *  else lives here, so adding a knob touches this file and the evolvers that
 *  honour it, rather than every call site.
 *
 *  The rule that makes this more than tidying: **an option an evolver does not
 *  consume is an ERROR, not an ignore.** Each evolver-specific field therefore
 *  has an "unset" state distinguishable from any value a caller would pass on
 *  purpose, and each evolver declares what it honours through
 *  EvolverOptionsCheck. Silently ignoring a knob is the failure #395 punished.
 */

#include <initializer_list>
#include <optional>

#include "evolver_erk.h"  // ErkControllerConfig

/* Callback shapes. Named once so a call site cannot transpose two of them and
   still compile -- three of the five take the same (double, double*, double*,
   void*) prefix. */
using EvolverDerivs = void (*)(double x, double* y, double* dy, void* parameters_and_workspace);
using EvolverOutput =
    void (*)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace);
using EvolverPrint = void (*)(double x, double y[], double dy[], void* parameters_and_workspace);
using EvolverTimescale = void (*)(double x, void* parameters_and_workspace, double* timescale);
using EvolverDiagonal = void (*)(double x, double* y, double* diag, void* parameters_and_workspace);

/** Step statistics, filled when a caller asks for them.
 *
 *  ndf15 has collected exactly these six counters into a local array and thrown
 *  them away since it was written; this is how they get out. Not thread-shared:
 *  one EvolverStats belongs to one evolver call. */
/** The value every evolver used to hardwire. Also the sentinel for "unset":
 *  abstol is feature-gated like the rest, so it needs a distinguishable default. */
inline constexpr double kEvolverDefaultAbsTol = 1.0e-15;

struct EvolverStats {
  /* Filled by every evolver. */
  long long steps_accepted     = 0;
  long long steps_rejected     = 0;
  long long derivs_evaluations = 0;

  /* Implicit methods only; left at zero by the explicit ones. */
  long long jacobians         = 0;
  long long lu_decompositions = 0;
  long long linear_solves     = 0;

  /* Evolvers with a continuous extension only; left at zero by the others.
     `dense` counts output points served by interpolation, `exact` those that
     landed on a step end. */
  long long dense_points = 0;
  long long exact_points = 0;

  /** Accumulate another run's counters, for a caller summing many calls. */
  void Add(const EvolverStats& other) {
    steps_accepted     += other.steps_accepted;
    steps_rejected     += other.steps_rejected;
    derivs_evaluations += other.derivs_evaluations;
    jacobians          += other.jacobians;
    lu_decompositions  += other.lu_decompositions;
    linear_solves      += other.linear_solves;
    dense_points       += other.dense_points;
    exact_points       += other.exact_points;
  }
};

/** Optional features. An evolver names the ones it honours; anything else that
 *  is set makes the call fail rather than quietly doing nothing. */
enum class EvolverFeature {
  Timescale,     /**< evaluate_timescale + timestep_over_timescale */
  Diagonal,      /**< derivs_diagonal */
  MaxOrder,      /**< max_order */
  ErkController, /**< erk */
  AbsTol,        /**< abstol (non-default) */
  Stats,         /**< stats */
  Histograms,    /**< histograms */
};

struct EvolverOptions {
  /* ---- honoured by every evolver ---------------------------------------- */

  /** Relative tolerance on the local error. */
  double rtol = 1.0e-6;

  /** Floor on the error weight: the magnitude below which a variable is judged
   *  in absolute rather than relative terms.
   *
   *  ndf15, etd and the explicit Runge-Kutta evolvers each hardwired 1e-15
   *  independently; the default preserves that everywhere, and they now read it
   *  from here. The legacy `rk` does NOT honour it -- its error scale is built
   *  inside generic_integrator from _TINY_ = 1e-30, a different quantity -- so
   *  setting a non-default value with that evolver is an error rather than a
   *  silent no-op. That is the same rule as every other feature here. */
  double abstol = kEvolverDefaultAbsTol;

  /* ---- dense output, honoured by every evolver --------------------------- */

  /** Ascending grid at which `output` is called. REQUIRED, non-null, and at
   *  least two points: ndf15 dereferences it before its first step, so a null is
   *  a crash rather than "no dense output". Pass the two endpoints if no dense
   *  output is wanted.
   *
   *  At least one point must lie at or beyond x_ini. Every evolver opens with
   *  `while (x_sampling[i] < x_ini) ++i` and no bound check, so a grid entirely
   *  below the start walks off the end of the array. */
  const double* x_sampling = nullptr;
  int x_sampling_size      = 0;

  /** Per-variable flag, one per equation: is this component needed by `output`?
   *  REQUIRED and non-null, like x_sampling: the explicit evolvers index it
   *  without checking (evolver_erk_impl.h), so a null is a crash rather than a
   *  useful "all of them" default. Pass an all-ones array to mean all. */
  const int* used_in_output = nullptr;

  /** Called at each x_sampling point. REQUIRED: every evolver invokes it without
   *  a null check, and x_sampling is required, so there is always at least one
   *  call. Pass a no-op if the samples are not wanted. */
  EvolverOutput output         = nullptr;
  EvolverPrint print_variables = nullptr;

  /* ---- evolver-specific; unset at the default, validated ----------------- */

  /** Fixed-step control for the legacy `rk` evolver only. */
  EvolverTimescale evaluate_timescale = nullptr;
  double timestep_over_timescale      = 0.0;

  /** Analytic Jacobian diagonal. Consumed by etd, which exponentiates it. */
  EvolverDiagonal derivs_diagonal = nullptr;

  /** Maximum BDF order for ndf15. 0 means the evolver's own default (5).
   *  Spelled 0-means-default rather than defaulting to 5 so that "the caller
   *  asked for 5" is distinguishable from "the caller said nothing". */
  int max_order = 0;

  /** Step-controller settings for the explicit Runge-Kutta evolvers. Unset
   *  means their defaults. This used to be process-wide state that three
   *  modules each had to remember to establish before evolving. */
  std::optional<ErkControllerConfig> erk;

  /** Where to write step statistics, or null for none.
   *
   *  Null is how collection is turned OFF -- there is no separate enable flag,
   *  and no process-wide one. Because the counters are per call they need no
   *  synchronisation, so unlike the relaxed atomics they replace they are cheap
   *  enough to leave on while timing. */
  EvolverStats* stats = nullptr;

  /** Where to write step-acceptance histograms, or null for none.
   *
   *  Separate from `stats` because, unlike the counters, these are NOT free:
   *  each step costs a log10 to find its bin. Explicit Runge-Kutta only. */
  ErkHistograms* histograms = nullptr;
};

/** Reject options the named evolver does not honour, and check the ones every
 *  evolver requires. Throws (severe: these are run-configuration choices a
 *  sampler never varies). */
void EvolverOptionsCheck(const EvolverOptions& options,
                         const char* evolver_name,
                         double x_ini,
                         std::initializer_list<EvolverFeature> honoured);

#endif  // CLASS_EVOLVER_OPTIONS_H
