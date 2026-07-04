/** @file output.h Output configuration shared with OutputModule. */

#ifndef __OUTPUT__
#define __OUTPUT__

#include <string>

#include "common.h"
#include "lensing.h"

/**
 * Maximum number of values of redshift at which the spectra will be
 * written in output files
 */

#define _Z_PK_NUM_MAX_ 100

/**
 * Structure containing various informations on the output format,
 * all of them initialized by user in input module.
 *
 */

struct output {
  //@{

  std::string root = "output/"; /**< root for all file names */

  //@}

  /** @name - number and value(s) of redshift at which P(k,z) and T_i(k,z) should be written */

  //@{

  int z_pk_num = 1; /**< number of redshift at which P(k,z) and T_i(k,z) should be written */
  double z_pk[_Z_PK_NUM_MAX_] = {
      0.}; /**< value(s) of redshift at which P(k,z) and T_i(k,z) should be written */

  //@}

  /** @name - extra information on output */

  //@{

  bool write_header = true; /**< flag stating whether we should write a header in output files */

  /** which format for output files (definitions, order of columns, etc.) */
  file_format output_format = file_format::class_format;

  bool write_background     = false; /**< flag for outputing background evolution in file */
  bool write_thermodynamics = false; /**< flag for outputing thermodynamical evolution in file */
  bool write_perturbations =
      false; /**< flag for outputing perturbations of selected wavenumber(s) in file(s) */
  bool write_primordial =
      false; /**< flag for outputing scalar/tensor primordial spectra in files */

  //@}

  /** @name - technical parameters */

  //@{

  short output_verbose =
      0; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}
};

#endif
