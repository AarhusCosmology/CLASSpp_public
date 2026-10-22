/** @file nonlinear.h Input data and shared types for nonlinear corrections. */

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
  extrap_user_defined,
  extrap_asymptotic /**< A ln k + B + C/k fitted over [k_max/3, k_max] (AsymptoticExtrapolation) */
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

/** HMcode versions. v2020_unfitted is the halo model without the fitted
    parameters, used internally by the baryonic feedback response of 2020
    (Mead et al. 2021, arXiv:2009.01858). */
enum class HMcodeVersion { v2016, v2020, v2020_unfitted, v2020_baryonic_feedback };

/**
 * Input parameters for nonlinear corrections.
 *
 * Computed nonlinear spectra belong to NonlinearModule; this structure holds
 * the shared configuration read by InputModule.
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
      extrap_asymptotic; /**< method for analytical extrapolation of sources beyond pre-computed range; extrap_max_scaled until 2026-10-21 */

  enum hmcode_baryonic_feedback_model feedback =
      nl_emu_dmonly; /** to choose between different baryonic feedback models
                                                in hmcode (dmonly, gas cooling, Agn or supernova feedback) */
  double c_min =
      0.0; /** for HMcode: minimum concentration in Bullock 2001 mass-concentration relation */
  double eta_0 = 0.0; /** for HMcode: halo bloating parameter */
  double z_infinity =
      10.; /** for HMcode: z value at which Dark Energy correction is evaluated; should be at early times (default: 10.) */
  HMcodeVersion hmcode_version = HMcodeVersion::
      v2020; /** HMcode version: 2016 (Mead et al. 2016) or 2020 (Mead et al. 2021), with or without baryonic feedback */
  double log10T_heat_hmcode =
      7.8; /** for HMcode 2020 with baryonic feedback: log10 of the AGN heating temperature in K */
  int nk_wiggle =
      512; /** for HMcode 2020: number of wavenumbers on which the linear spectrum is split into wiggle and no-wiggle parts */

  //@}

  bool has_pk_eq = false; /**< flag: will we use the pk_eq method? */

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

  std::vector<double> rtab;   /** List of R values */
  std::vector<double> stab;   /** List of Sigma Values */
  std::vector<double> ddstab; /** Splined sigma */

  std::vector<double> growtable;
  std::vector<double> ztable;
  std::vector<double> tautable;

  std::vector<double*> sigma_8;        /* row pointers into sigma_8_storage */
  std::vector<double*> sigma_disp;     /* row pointers into sigma_disp_storage */
  std::vector<double*> sigma_disp_100; /* row pointers into sigma_disp_100_storage */
  std::vector<double*> sigma_prime;    /* row pointers into sigma_prime_storage */

  std::vector<std::vector<double>> sigma_8_storage;
  std::vector<std::vector<double>> sigma_disp_storage;
  std::vector<std::vector<double>> sigma_disp_100_storage;
  std::vector<std::vector<double>> sigma_prime_storage;

  double dark_energy_correction =
      0.0; /** this is the ratio [g(z_infinity)/g_lcdm(z_infinity)]^1.5 for
                                  * HMcode 2016 (power from Dolag et al. (2004) as in Mead et al. 2016)
                                  * and g(z_infinity)/g_lcdm(z_infinity) for 2020, with g_lcdm the
                                  * growth of a flat LCDM with the same Omega_m, in every cosmology.
                                  * For 2020, g is the growth of the cosmology without radiation below */

  /** HMcode 2020: linear growth in the cosmology without radiation of Mead et al. (2021),
      on a_norad (logarithmic in a). Columns of growth_norad: g/a, the integral of g/a over a,
      Omega_m(a); their splines in a are in dd_growth_norad. g_norad is g normalised to 1 today. */
  std::vector<double> a_norad;
  std::vector<double> growth_norad;
  std::vector<double> dd_growth_norad;
  std::vector<double> g_norad;
  std::vector<double> dd_g_norad; /**< spline of g_norad(a) */
  std::vector<double> dd_a_norad; /**< spline of a(g_norad), to find formation times */

  /** HMcode 2020: no-wiggle reference spectrum of Eisenstein & Hu on the extended k grid
      (ln P and its spline in ln k), and the Gaussian weights of the smoothing in ln k */
  std::vector<double> ln_pk_eh;
  std::vector<double> dd_ln_pk_eh;
  std::vector<double>
      smoothing_weights; /**< nk_wiggle x nk_wiggle, rows normalised; empty rows are not smoothed */

  /** HMcode 2020: wiggle part of the linear spectrum at the current time, as a dimensionless
      spectrum k^3/(2 pi^2) P_wiggle, on lnk_wiggle */
  std::vector<double> lnk_wiggle;
  std::vector<double> pk_wiggle;
  std::vector<double> ddpk_wiggle;

  //@}
};

#endif
