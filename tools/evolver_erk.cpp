#include "evolver_erk.h"

#include <cmath>

#include "evolver_erk_impl.h"

/**
 * The explicit embedded Runge-Kutta evolvers, each binding one tableau to the
 * templated driver in evolver_erk_impl.h, and the driver's out-of-line step
 * recording. Everything numerical lives in the header so the tableaux stay
 * constexpr at the point of use; tsit5 and rkdp45 therefore differ in their
 * coefficients and nothing else.
 */

namespace {

int bin_of(double value, double lo, double step, int nbins) {
  if (!(value > 0.0)) {
    return 0;
  }
  const int b = static_cast<int>((log10(value) - lo) / step);
  return b < 0 ? 0 : (b >= nbins ? nbins - 1 : b);
}

}  // namespace

namespace erk_detail {

void Record(
    long long* counter, ErkHistograms* histograms, bool accepted, double x, double err_ratio) {
  ++*counter;
  /* The counter above is a plain increment. The binning below is not -- it costs
     a log10 per step -- so it is gated on the caller having asked for it. */
  if (histograms == nullptr) {
    return;
  }
  long long* err = accepted ? histograms->err_accepted : histograms->err_rejected;
  long long* pos = accepted ? histograms->x_accepted : histograms->x_rejected;
  err[bin_of(err_ratio, ErkHistograms::kErrLo, ErkHistograms::kErrStep, ErkHistograms::kErrBins)] +=
      1;
  pos[bin_of(fabs(x), ErkHistograms::kXLo, ErkHistograms::kXStep, ErkHistograms::kXBins)] += 1;
}

}  // namespace erk_detail

static constexpr std::initializer_list<EvolverFeature> kErkFeatures =
    {EvolverFeature::ErkController,
     EvolverFeature::Stats,
     EvolverFeature::Histograms,
     EvolverFeature::AbsTol};

void evolver_rkdp45(EvolverDerivs derivs,
                    double x_ini,
                    double x_end,
                    double* y,
                    int y_size,
                    const EvolverOptions& options) {
  EvolverOptionsCheck(options, "evolver_rkdp45", x_ini, kErkFeatures);
  evolver_erk_run<ErkDormandPrince45>(derivs, x_ini, x_end, y, y_size, options);
}

void evolver_tsit5(EvolverDerivs derivs,
                   double x_ini,
                   double x_end,
                   double* y,
                   int y_size,
                   const EvolverOptions& options) {
  EvolverOptionsCheck(options, "evolver_tsit5", x_ini, kErkFeatures);
  evolver_erk_run<ErkTsitouras54>(derivs, x_ini, x_end, y, y_size, options);
}
