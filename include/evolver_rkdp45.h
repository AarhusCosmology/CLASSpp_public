#ifndef __EVORKDP45__
#define __EVORKDP45__

#include "common.h"

/**************************************************************/

/** Dormand-Prince 5(4), adaptive. Honours EvolverFeature::ErkController. */
#include "evolver_options.h"

void evolver_rkdp45(EvolverDerivs derivs,
                    double x_ini,
                    double x_end,
                    double* y,
                    int y_size,
                    const EvolverOptions& options);

/**************************************************************/

#endif
