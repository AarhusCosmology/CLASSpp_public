/** @file transfer.h Input data and shared types for transfer functions. */

#ifndef __TRANSFER__
#define __TRANSFER__

#include <string>
#include <vector>

#include "hyperspherical.h"
#include "nonlinear.h"

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

  std::string
      selection_file_name; /**< for selection=tabulated: file with the redshift z in the first column and the selection function dN/dz of each bin in the next ones */
  std::vector<double> selection_file_z; /**< for selection=tabulated: redshifts of this file */
  std::vector<std::vector<double>>
      selection_file_dNdz; /**< for selection=tabulated: selection function of each bin, selection_file_dNdz[bin][index_z] */
  std::vector<double>
      selection_z_min; /**< for selection=tabulated: smallest redshift at which the selection function of each bin is non-zero (one tabulated point before the first non-zero value) */
  std::vector<double>
      selection_z_max; /**< for selection=tabulated: largest redshift at which the selection function of each bin is non-zero (one tabulated point after the last non-zero value) */

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
 * enumeration of possible radial function types (Bessel function or
 * combination of its derivatives) convolved with the transfer source,
 * see transfer_radial_function()
 */

using radial_function_type = enum {
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
};

/**
 * The physical contribution carried by a transfer type: a CMB term, or
 * one term of the number counts (see the CLASSgal paper arXiv:1307.1459)
 * or of the galaxy lensing potential in a given redshift bin.
 */
enum class TransferTerm {
  cmb_t0,
  cmb_t1,
  cmb_t2,
  cmb_e,
  cmb_b,
  cmb_lensing,
  /* local in a redshift bin */
  nc_density,
  nc_rsd,
  nc_d0,
  nc_d1,
  nc_g1,
  nc_g2,
  nc_g3,
  /* integrated along the line of sight from a redshift bin to us */
  nc_lens,
  nc_g4,
  nc_g5,
  galaxy_lensing
};

/** The window function W(tau) multiplying the perturbation source of a transfer type. */
enum class TransferWindow {
  none,        /**< W=1 (CMB temperature and polarisation) */
  cmb_lensing, /**< CMB lensing kernel, with support [tau_rec, tau0] */
  selection, /**< selection function of a bin times a function of the background, with the support of the bin */
  integrated /**< selection function of a bin convolved with a kernel, with support [tau_min of the bin, tau0] */
};

/** The ways of sampling a line-of-sight integral in time. */
enum class TransferSampling {
  from_perturbations, /**< the times of the perturbation module larger than tau_min */
  uniform, /**< tau_size evenly spaced times between tau_min and tau_max (just tau_min if tau_size=1) */
  logarithmic /**< the times of the perturbation module larger than tau_min as long as their step is smaller than (tau0-tau)*log_step; then this logarithmic step until tau_max; then tau0 */
};

/** The ways of computing a transfer function at a given (q, l). */
enum class TransferMethod {
  zero,         /**< it vanishes or can be neglected */
  limber,       /**< Limber approximation */
  line_of_sight /**< line-of-sight integral of the transfer source times a radial function */
};

/**
 * Everything that distinguishes a transfer type X, defined once for all
 * in TransferModule::transfer_type_properties(). The transfer source is
 * \f$ S^{Y(X)}(k(q),\tau) \, W^X(\tau) \, (k/k_{pivot})^{k_{power}} \f$ for every type.
 */
struct TransferType {
  TransferTerm term{};                   /**< physical contribution */
  int index_tp                     = -1; /**< index of the perturbation source Y(X) */
  int bin                          = -1; /**< redshift bin (-1 for CMB types) */
  int l_max                        = 0;  /**< largest multipole needed */
  radial_function_type radial_type = SCALAR_TEMPERATURE_0; /**< radial function */
  TransferWindow window_type       = TransferWindow::none; /**< kind of window W(tau) */
  double k_power                   = 0.;                   /**< power of the factor (k/k_pivot) */
  double k_pivot                   = 1.;                   /**< pivot of the factor (k/k_pivot) */
  double l_limber        = _HUGE_; /**< the Limber approximation is used for l >= l_limber */
  double delta_k_neglect = _HUGE_; /**< Delta_l(q)=0 for l < (q-delta_k_neglect) r_a(tau_rec) */
  double l_half_neglect =
      0.; /**< if > 0, the margin at l is delta_k_neglect sqrt(l/(l+l_half_neglect)) */
  bool late_source_cut =
      false; /**< sources at tau > tau_cut can be neglected for l > ppr->transfer_neglect_late_source */
  bool full_limber = false;   /**< also computed on the full-Limber grid q_limber */
  int index_ts     = -1;      /**< index of its time sampling */
  std::vector<double> window; /**< W(tau) on this time sampling */
};

/**
 * A time sampling of the line-of-sight integrals, shared by all transfer
 * types with an identical sampling, defined once for all in
 * TransferModule::transfer_time_samplings().
 */
struct TransferTimeSampling {
  std::vector<double> tau0_minus_tau; /**< sampled values of (tau0-tau), in decreasing order */
  std::vector<double> w_trapz;        /**< weights of the trapezoidal integration rule */
  std::vector<int>
      index_tau_pt; /**< for each sampled time, index of the perturbation time just before it, for linear interpolation of the perturbation sources */
  std::vector<double> weight_tau_pt; /**< corresponding weight of the next perturbation time */
  double q_max_sampled =
      _HUGE_; /**< largest q at which the sampling resolves the Bessel functions with ppr->selection_sampling_bessel_min points per radian (evenly spaced samplings only); above it, the Limber approximation is used */
  std::vector<int> types; /**< the transfer types using this sampling */
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

  const HyperInterpStruct*
      pBIS; /**< pointer to structure containing all the spherical bessel functions of the flat case (used even in the non-flat case, for approximation schemes). pBIS = pointer to Bessel Interpolation Structure. */

  //@}

  /** @name - transfer sources of the current wavenumber and time sampling */

  //@{

  std::vector<std::vector<double>>
      interpolated_sources; /**< interpolated_sources[index_ic*tp_size+index_tp][index_tau]: perturbation sources interpolated at the current k, on the time sampling of the perturbation module (only allocated for the sources used by some type) */
  std::vector<char>
      is_interpolated; /**< has the previous source been interpolated at the current k? */
  std::vector<double>
      sources; /**< sources[index_n*tau_size_max+index_tau]: transfer source of the n-th (initial condition, type) pair using the current time sampling */
  std::vector<TransferMethod> method; /**< method[index_n] for the current l */
  std::vector<int>
      index_tau_max; /**< index_tau_max[index_n]: last sampled time contributing to the line-of-sight integral */
  std::vector<char>
      bessel_truncation; /**< bessel_truncation[index_n]: is the integral truncated by the lower cut of the Bessel function? */

  std::vector<double> chi;     /**< chi[index_tau]: value of argument of bessel
                                    function: k(tau0-tau) (flat case)
                                    or sqrt(|K|)(tau0-tau) (non-flat
                                    case) */
  std::vector<double> cscKgen; /**< cscKgen[index_tau]: useful trigonometric function */
  std::vector<double> cotKgen; /**< cotKgen[index_tau]: useful trigonometric function */

  /** Bessel functions of the current (q, l) at the first bessel_size times of the current
   *  sampling, interpolated once for all types by transfer_bessel_interpolate(), in reverse
   *  order (order of increasing argument): the value at index_tau is at bessel_size-1-index_tau.
   *  Buffers of size tau_size_max, allocated once per wavenumber. */
  int bessel_size; /**< number of times at which they have been interpolated */
  double
      rescale_argument; /**< rescaling of their argument in the flat rescaling approximation (curved space) */
  std::vector<double> Phi;        /**< Bessel function */
  std::vector<double> dPhi;       /**< first derivative of Bessel function */
  std::vector<double> d2Phi;      /**< second derivative of Bessel function */
  std::vector<double> chireverse; /**< their argument (curved space) */
  std::vector<double>
      rescale_function; /**< amplitude rescaling of the flat rescaling approximation (curved space) */
  std::vector<double>
      radial_function; /**< radial_function[index_tau]: output of transfer_radial_function */
  std::vector<double> chi_full_reverse; /**< workspace: reversed chi array */
  std::vector<int>
      last_nonzero_source; /**< per (initial condition, type) of a time sampling: the last
                                           time index at or before i with a nonzero source, or -1 */
  HermiteLocations
      bessel_location; /**< flat case: where each chi_full_reverse falls in the Bessel table, found once per
                          wavenumber and time sampling and shared by every multipole */

  //@}

  /** @name - parameters defining the spatial curvature (copied from background structure) */

  //@{

  double K; /**< curvature parameter (see background module for details) */
  int sgnK; /**< 0 (flat), 1 (positive curvature, spherical, closed), -1 (negative curvature, hyperbolic, open) */

  //@}

  double
      tau0_minus_tau_cut; /**< critical value of (tau0-tau) in time cut approximation for the wavenumber at hand */
};

enum Hermite_Interpolation_Order { HERMITE4, HERMITE6 };

#endif
