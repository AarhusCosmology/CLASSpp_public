/** @file spectra.h Input data and shared types for power spectra. */

#ifndef __SPECTRA__
#define __SPECTRA__

#include "transfer.h"

/**
 * Input parameters for harmonic and matter power spectra.
 *
 * Computed spectra belong to SpectraModule; this structure holds the shared
 * configuration read by InputModule.
 */

struct spectra {
  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these parameters
      and the content of the 'background', 'perturbs', 'transfers' and
      'primordial' structures) */

  //@{

  double z_max_pk = 0.; /**< maximum z for P(k,z) output; set from perturbs by the reader */

  int non_diag = 0; /**< sets the number of cross-correlation spectra
                   that you want to calculate: 0 means only
                   auto-correlation, 1 means only adjacent bins,
                   and number of bins minus one means all
                   correlations */

  /** @name - technical parameters */

  //@{

  short spectra_verbose =
      0; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}
};

/*************************************************************************************************************/
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

/* end deprecated functions */

#ifdef __cplusplus
}
#endif

#endif
