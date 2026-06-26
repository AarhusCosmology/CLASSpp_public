/** @file transfer.h Input data and shared types for transfer functions. */

#ifndef __TRANSFER__
#define __TRANSFER__

#include <string>
#include <vector>

#include "hyperspherical.h"
#include "nonlinear.h"

/* macro: test if index_tt is in the range between index and index+num, while the flag is true */
#define _index_tt_in_range_(index, num, flag) \
  (flag == _TRUE_) && (index_tt >= index) && (index_tt < index + num)
/* macro: test if index_tt corresponds to an integrated nCl/sCl contribution */
#define _integrated_ncl_                                                                         \
  (_index_tt_in_range_(index_tt_lensing_, ppt->selection_num, ppt->has_cl_lensing_potential)) || \
      (_index_tt_in_range_(index_tt_nc_lens_, ppt->selection_num, ppt->has_nc_lens)) ||          \
      (_index_tt_in_range_(index_tt_nc_g4_, ppt->selection_num, ppt->has_nc_gr)) ||              \
      (_index_tt_in_range_(index_tt_nc_g5_, ppt->selection_num, ppt->has_nc_gr))
/* macro: test if index_tt corresponds to an non-integrated nCl/sCl contribution */
#define _nonintegrated_ncl_                                                            \
  (_index_tt_in_range_(index_tt_density_, ppt->selection_num, ppt->has_nc_density)) || \
      (_index_tt_in_range_(index_tt_rsd_, ppt->selection_num, ppt->has_nc_rsd)) ||     \
      (_index_tt_in_range_(index_tt_d0_, ppt->selection_num, ppt->has_nc_rsd)) ||      \
      (_index_tt_in_range_(index_tt_d1_, ppt->selection_num, ppt->has_nc_rsd)) ||      \
      (_index_tt_in_range_(index_tt_nc_g1_, ppt->selection_num, ppt->has_nc_gr)) ||    \
      (_index_tt_in_range_(index_tt_nc_g2_, ppt->selection_num, ppt->has_nc_gr)) ||    \
      (_index_tt_in_range_(index_tt_nc_g3_, ppt->selection_num, ppt->has_nc_gr))
/* macro: bin number associated to particular redshift bin and selection function for non-integrated contributions*/
#define _get_bin_nonintegrated_ncl_(index_tt)                                          \
  if (_index_tt_in_range_(index_tt_density_, ppt->selection_num, ppt->has_nc_density)) \
    bin = index_tt - index_tt_density_;                                                \
  if (_index_tt_in_range_(index_tt_rsd_, ppt->selection_num, ppt->has_nc_rsd))         \
    bin = index_tt - index_tt_rsd_;                                                    \
  if (_index_tt_in_range_(index_tt_d0_, ppt->selection_num, ppt->has_nc_rsd))          \
    bin = index_tt - index_tt_d0_;                                                     \
  if (_index_tt_in_range_(index_tt_d1_, ppt->selection_num, ppt->has_nc_rsd))          \
    bin = index_tt - index_tt_d1_;                                                     \
  if (_index_tt_in_range_(index_tt_nc_g1_, ppt->selection_num, ppt->has_nc_gr))        \
    bin = index_tt - index_tt_nc_g1_;                                                  \
  if (_index_tt_in_range_(index_tt_nc_g2_, ppt->selection_num, ppt->has_nc_gr))        \
    bin = index_tt - index_tt_nc_g2_;                                                  \
  if (_index_tt_in_range_(index_tt_nc_g3_, ppt->selection_num, ppt->has_nc_gr))        \
    bin = index_tt - index_tt_nc_g3_;
/* macro: bin number associated to particular redshift bin and selection function for integrated contributions*/
#define _get_bin_integrated_ncl_(index_tt)                                                       \
  if (_index_tt_in_range_(index_tt_lensing_, ppt->selection_num, ppt->has_cl_lensing_potential)) \
    bin = index_tt - index_tt_lensing_;                                                          \
  if (_index_tt_in_range_(index_tt_nc_lens_, ppt->selection_num, ppt->has_nc_lens))              \
    bin = index_tt - index_tt_nc_lens_;                                                          \
  if (_index_tt_in_range_(index_tt_nc_g4_, ppt->selection_num, ppt->has_nc_gr))                  \
    bin = index_tt - index_tt_nc_g4_;                                                            \
  if (_index_tt_in_range_(index_tt_nc_g5_, ppt->selection_num, ppt->has_nc_gr))                  \
    bin = index_tt - index_tt_nc_g5_;
/**
 * Input parameters for transfer functions.
 *
 * Computed transfer-function tables belong to TransferModule; this structure
 * holds the shared configuration read by InputModule.
 */

struct transfers {
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of previous structures) */

  //@{

  double lcmb_rescale = 1.;  /**< normally set to one, can be used
                          exceptionally to rescale by hand the CMB
                          lensing potential */
  double lcmb_tilt    = 0.;  /**< normally set to zero, can be used
                          exceptionally to tilt by hand the CMB
                          lensing potential */
  double lcmb_pivot   = 0.1; /**< if lcmb_tilt non-zero, corresponding pivot
                          scale */

  double selection_bias[_SELECTION_NUM_MAX_] = {
      1.}; /**< light-to-mass bias; element [0] default 1, rest 0 */
  double selection_magnification_bias[_SELECTION_NUM_MAX_] = {
      0.}; /**< magnification bias; default 0 */

  bool has_nz_file = false; /**< Has dN/dz (selection function) input file? */
  bool has_nz_analytic =
      false;                /**< Use analytic form for dN/dz (selection function) distribution? */
  std::string nz_file_name; /**< dN/dz (selection function) input file name */

  bool has_nz_evo_file = false; /**< Has dN/dz (evolution function) input file? */
  bool has_nz_evo_analytic =
      false; /**< Use analytic form for dN/dz (evolution function) distribution? */
  std::string nz_evo_file_name; /**< dN/dz (evolution function) input file name */

  //@}

  /** @name - technical parameters */

  //@{

  bool initialise_HIS_cache =
      false; /**< only true if we are using CLASS for setting up a cache of HIS structures */

  short transfer_verbose =
      0; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}
};

/**
 * Structure containing all the quantities that each thread needs to
 * know for computing transfer functions (but that can be forgotten
 * once the transfer functions are known, otherwise they would be
 * stored in the transfer module)
*/

struct transfer_workspace {
  /** @name - quantities related to Bessel functions */

  //@{

  HyperInterpStruct
      HIS; /**< structure containing all hyperspherical bessel functions (flat case) or all hyperspherical bessel functions for a given value of beta=q/sqrt(|K|) (non-flat case). HIS = Hyperspherical Interpolation Structure. */

  HyperInterpStruct*
      pBIS; /**< pointer to structure containing all the spherical bessel functions of the flat case (used even in the non-flat case, for approximation schemes). pBIS = pointer to Bessel Interpolation Structure. */

  int l_size; /**< number of l values */

  //@}

  /** @name - quantities related to the integrand of the transfer functions (most of them are arrays of time) */

  //@{

  int tau_size;     /**< number of discrete time values for a given type */
  int tau_size_max; /**< maximum number of discrete time values for all types */
  std::vector<double> interpolated_sources; /**< interpolated_sources[index_tau]:
                                    sources interpolated from the
                                    perturbation module at the right
                                    value of k */
  std::vector<double> sources;              /**< sources[index_tau]: sources
                                    used in transfer module, possibly
                                    differing from those in the
                                    perturbation module by some
                                    resampling or rescaling */
  std::vector<double> tau0_minus_tau; /**< tau0_minus_tau[index_tau]: values of (tau0 - tau) */
  std::vector<double>
      w_trapz; /**< w_trapz[index_tau]: values of weights in trapezoidal integration (related to time steps) */
  std::vector<double> chi;     /**< chi[index_tau]: value of argument of bessel
                                    function: k(tau0-tau) (flat case)
                                    or sqrt(|K|)(tau0-tau) (non-flat
                                    case) */
  std::vector<double> cscKgen; /**< cscKgen[index_tau]: useful trigonometric function */
  std::vector<double> cotKgen; /**< cotKgen[index_tau]: useful trigonometric function */

  /** Pre-allocated temporary buffers for transfer_radial_function and transfer_integrate.
   *  Size tau_size_max. Avoids repeated malloc/free in the hot per-(k,l) loop. */
  std::vector<double> Phi;        /**< Phi[index_tau]: Bessel function values */
  std::vector<double> dPhi;       /**< dPhi[index_tau]: first derivative of Bessel function */
  std::vector<double> d2Phi;      /**< d2Phi[index_tau]: second derivative of Bessel function */
  std::vector<double> chireverse; /**< chireverse[index_tau]: reversed chi grid */
  std::vector<double> rescale_function; /**< rescale_function[index_tau]: amplitude rescaling */
  std::vector<double>
      radial_function; /**< radial_function[index_tau]: output of transfer_radial_function */
  std::vector<double> chi_full_reverse; /**< workspace: reversed chi array */

  //@}

  /** @name - parameters defining the spatial curvature (copied from background structure) */

  //@{

  double K; /**< curvature parameter (see background module for details) */
  int sgnK; /**< 0 (flat), 1 (positive curvature, spherical, closed), -1 (negative curvature, hyperbolic, open) */

  //@}

  double
      tau0_minus_tau_cut; /**< critical value of (tau0-tau) in time cut approximation for the wavenumber at hand */
  short
      neglect_late_source; /**< flag stating whether we use the time cut approximation for the wavenumber at hand */
};

/**
 * enumeration of possible source types. This looks redundant with
 * respect to the definition of indices index_tt_... This definition is however
 * convenient and time-saving: it allows to use a "case" statement in
 * transfer_radial_function()
 */

typedef enum {
  SCALAR_TEMPERATURE_0,
  SCALAR_TEMPERATURE_1,
  SCALAR_TEMPERATURE_2,
  SCALAR_POLARISATION_E,
  VECTOR_TEMPERATURE_1,
  VECTOR_TEMPERATURE_2,
  VECTOR_POLARISATION_E,
  VECTOR_POLARISATION_B,
  TENSOR_TEMPERATURE_2,
  TENSOR_POLARISATION_E,
  TENSOR_POLARISATION_B,
  NC_RSD
} radial_function_type;

enum Hermite_Interpolation_Order { HERMITE3, HERMITE4, HERMITE6 };

#endif
