#ifndef __EVORKCK__
#define __EVORKCK__

#include "dei_rkck.h"

/**************************************************************/

/** Legacy fixed-step Runge-Kutta-Cash-Karp evolver.
 *
 *  The ONLY evolver that honours EvolverFeature::Timescale: it takes its step
 *  from a caller-supplied timescale rather than from an error estimate. */
#include "evolver_options.h"

void evolver_rk(EvolverDerivs derivs,
                double x_ini,
                double x_end,
                double* y,
                int y_size,
                const EvolverOptions& options);

/**************************************************************/

#endif
