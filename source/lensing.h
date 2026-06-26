/** @file lensing.h Input data and shared types for CMB lensing. */

#ifndef __LENSING__
#define __LENSING__

#include "spectra.h"

/**
 * Input parameters for CMB lensing.
 *
 * Computed lensed spectra belong to LensingModule; this structure holds the
 * shared configuration read by InputModule.
 */

struct lensing {
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  bool has_lensed_cls = false; /**< do we need to compute lensed \f$ C_l\f$'s at all ? */

  //@}

  /** @name - technical parameters */

  //@{

  short lensing_verbose =
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

#ifdef __cplusplus
}
#endif

#endif
