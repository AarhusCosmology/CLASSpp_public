#ifndef __EVOTSIT5__
#define __EVOTSIT5__

#include "common.h"

/**************************************************************/

/** Tsitouras 5(4), adaptive. Honours EvolverFeature::ErkController. */
#include "evolver_options.h"

void evolver_tsit5(EvolverDerivs derivs,
                   double x_ini,
                   double x_end,
                   double* y,
                   int y_size,
                   void* parameters_and_workspace_for_derivs,
                   const EvolverOptions& options);

/**************************************************************/

#endif
