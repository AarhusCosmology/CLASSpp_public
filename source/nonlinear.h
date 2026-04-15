/** @file nonlinear.h Documented includes for trg module */

#include <vector>

#include "primordial.h"
#include "trigonometric_integrals.h"

#ifndef __NONLINEAR__
#define __NONLINEAR__

#define _M_EV_TOO_BIG_FOR_HALOFIT_ \
  10. /**< above which value of non-CDM mass (in eV) do we stop trusting halofit? */

#define _M_SUN_ 1.98847e30 /**< Solar mass in Kg */

#define _MAX_NUM_EXTRAPOLATION_ 100000

enum non_linear_method { nl_none, nl_halofit, nl_HMcode };
enum pk_outputs { pk_linear, pk_nonlinear };

enum source_extrapolation {
  extrap_zero,
  extrap_only_max,
  extrap_only_max_units,
  extrap_max_scaled,
  extrap_hmcode,
  extrap_user_defined
};

enum halofit_integral_type { halofit_integral_one, halofit_integral_two, halofit_integral_three };

enum hmcode_baryonic_feedback_model {
  nl_emu_dmonly,
  nl_owls_dmonly,
  nl_owls_ref,
  nl_owls_agn,
  nl_owls_dblim,
  nl_user_defined
};
enum out_sigmas { out_sigma, out_sigma_prime, out_sigma_disp };

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct nonlinear {
  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

  enum non_linear_method method =
      nl_none; /**< method for computing non-linear corrections (none, Halogit, etc.) */

  enum source_extrapolation extrapolation_method =
      extrap_max_scaled; /**< method for analytical extrapolation of sources beyond pre-computed range */

  enum hmcode_baryonic_feedback_model feedback =
      nl_emu_dmonly; /** to choose between different baryonic feedback models
                                                in hmcode (dmonly, gas cooling, Agn or supernova feedback) */
  double c_min; /** for HMcode: minimum concentration in Bullock 2001 mass-concentration relation */
  double eta_0; /** for HMcode: halo bloating parameter */
  double z_infinity =
      10.; /** for HMcode: z value at which Dark Energy correction is evaluated; should be at early times (default: 10.) */

  //@}

  short has_pk_eq = _FALSE_; /**< flag: will we use the pk_eq method? */

  /** @name - technical parameters */

  //@{

  short nonlinear_verbose = 0; /**< amount of information written in standard output */

  //@}
};

/**
 * Structure containing variables used only internally in nonlinear module by various functions.
 *
 */

struct nonlinear_workspace {
  /** @name - quantitites used by HMcode */

  //@{

  double* rtab;   /** List of R values */
  double* stab;   /** List of Sigma Values */
  double* ddstab; /** Splined sigma */

  double* growtable;
  double* ztable;
  double* tautable;

  double** sigma_8;
  double** sigma_disp;
  double** sigma_disp_100;
  double** sigma_prime;

  std::vector<double> rtab_storage;
  std::vector<double> stab_storage;
  std::vector<double> ddstab_storage;

  std::vector<double> growtable_storage;
  std::vector<double> ztable_storage;
  std::vector<double> tautable_storage;

  std::vector<std::vector<double>> sigma_8_storage;
  std::vector<std::vector<double>> sigma_disp_storage;
  std::vector<std::vector<double>> sigma_disp_100_storage;
  std::vector<std::vector<double>> sigma_prime_storage;

  std::vector<double*> sigma_8_rows;
  std::vector<double*> sigma_disp_rows;
  std::vector<double*> sigma_disp_100_rows;
  std::vector<double*> sigma_prime_rows;

  double dark_energy_correction; /** this is the ratio [g_wcdm(z_infinity)/g_lcdm(z_infinity)]^1.5
                                  * (power comes from Dolag et al. (2004) correction)
                                  * it is 1, if has_fld == _FALSE_ */

  //@}
};

#endif
