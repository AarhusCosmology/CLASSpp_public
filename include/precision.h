/** @file precision.h Precision parameters and the enums they use. */
#ifndef CLASS_PRECISION_H
#define CLASS_PRECISION_H

#include <float.h>  // DBL_EPSILON
#include <string>

#include "constants.h"  // physics constants used in precision defaults
#include "errors.h"     // status codes

class FileContent;  // forward decl for precision::parse

/** parameters related to the precision of the code and to the method of calculation */

/**
 * list of evolver types for integrating perturbations over time
 */
enum class evolver_type {
  rk,     /* Runge-Kutta integrator */
  ndf15,  /* stiff integrator */
  rkdp45, /* Dormand-Prince 4(5) explicit adaptive integrator */
  etd     /* exponential Rosenbrock (ETDRK4, Cox-Matthews); integrates a species-reported
             Jacobian DIAGONAL exactly, so a diagonal stiffness stops limiting
             the step at explicit cost. Available on the background and the
             perturbations; species that report no diagonal get the explicit
             scheme, since at diag = 0 the phi functions reduce ETDRK4 to
             classical RK4. Thermodynamics passes no diagonal.

             Worth it where a diagonal stiffness dominates, which on the dncdm_inv
             perturbations it does across the measured range (Gamma = 1e5 .. 1e9,
             1.2x to 4.9x faster than rkdp45). Elsewhere the diagonal is zero and
             this is RK4 with a per-step overhead, so it is opt-in per module. */
};

/**
 * List of ways in which matter power spectrum P(k) can be defined.
 * The standard definition is the first one (delta_m_squared) but
 * alternative definitions can be useful in some projects.
 *
 */
enum class pk_def {
  delta_m_squared, /**< normal definition (delta_m includes all non-relativistic species at late times) */
  delta_tot_squared, /**< delta_tot includes all species contributions to (delta rho), and only non-relativistic contributions to rho */
  delta_bc_squared, /**< delta_bc includes contribution of baryons and cdm only to (delta rho) and to rho */
  delta_tot_from_poisson_squared /**< use delta_tot inferred from gravitational potential through Poisson equation */
};
/**
 * Different ways to present output files
 */

enum class file_format { class_format, camb_format };

/** Approximation scheme labels used as precision parameter defaults */

enum class tca_method {
  first_order_MB,
  first_order_CAMB,
  first_order_CLASS,
  second_order_CRS,
  second_order_CLASS,
  compromise_CLASS
};
enum class rsa_method { rsa_null, rsa_MD, rsa_MD_with_reio, rsa_none };
enum class rsa_idr_method { rsa_idr_none, rsa_idr_MD };
enum class ufa_method { ufa_mb, ufa_hu, ufa_CLASS, ufa_none };
enum class ncdmfa_method { ncdmfa_mb, ncdmfa_hu, ncdmfa_CLASS, ncdmfa_none };

/**
 * All precision parameters.
 *
 * Includes integrations
 * steps, flags telling how the computation is to be performed, etc.
 */
struct precision {
  std::string class_dir;

  /*
   * Background Quantities
   * */

  /**
   * Default initial value of scale factor used in the integration of background quantities.
   * For models like ncdm, the code may decide to start the integration earlier.
   */
  double a_ini_over_a_today_default = 1.e-14;
  /**
   * Default stepsize in conformal time for the background integration,
   * in units for the conformal Hubble time. dtau = back_integration_stepsize/aH
   */
  double back_integration_stepsize = 7.e-3;
  /**
   * Tolerance of the background integration, giving the allowed relative integration error.
   * Passed to the evolver in background_solve_evolver(). The default is 1e-6 (not the
   * historical 1e-2): the removed legacy RK solver had no dense output, so the requested
   * output points capped the step size and the step count acted as the effective precision
   * knob; the evolver instead relies on this tolerance directly, so it must be tight.
   */
  double tol_background_integration = 1.e-6;
  /**
   * Tolerance of the deviation of \f$ \Omega_r \f$ from 1 for which to start integration:
   * The starting point of integration will be chosen,
   * such that the Omega of radiation at that point is close to 1 within tolerance.
   * (Class starts background integration during complete radiation domination)
   */
  double tol_initial_Omega_r = 1.e-4;
  /**
   * Tolerance of relative deviation of the used non-cold dark matter mass compared to that which would give the correct density.
   * The dark matter mass is estimated from the dark matter density using a Newton-Method.
   * In the nonrelativistic limit, this could be estimated using M=density/number density
   */
  double tol_M_ncdm = 1.e-7;
  /**
   * Tolerance on the relative precision of the integration over
   * non-cold dark matter phase-space distributions in the synchronous gauge.
   * The convenience input key "tol_ncdm" sets both this and tol_ncdm_newtonian.
   */
  double tol_ncdm_synchronous = 3.4e-3;
  /**
   * Tolerance on the relative precision of the integration over
   * non-cold dark matter phase-space distributions in the newtonian gauge.
   * Newtonian gauge needs a finer ncdm momentum grid than synchronous for the
   * same P(k) accuracy, hence the tighter default.
   */
  double tol_ncdm_newtonian = 9.05e-5;
  /**
   * Tolerance on the relative precision of the integration over
   * non-cold dark matter phase-space distributions during the background evolution.
   */
  double tol_ncdm_bg = 1.e-5;
  /**
   * Tolerance on the initial deviation of non-cold dark matter from being fully relativistic.
   * Using w = pressure/density, this quantifies the maximum deviation from 1/3. (for relativistic species)
   */
  double tol_ncdm_initial_w = 1.e-3;
  /**
   * Tolerance on the deviation of the conformal time of equality from the true value in 1/Mpc.
   */
  double tol_tau_eq = 1.e-6;
  /**
   * Minimum amount of cdm to allow calculations in synchronous gauge comoving with cdm.
   */
  double Omega0_cdm_min_synchronous = 1.e-10;
  /*
   * Currently unused parameter.
   */
  //class_precision_parameter(safe_phi_scf,double,0.0)
  /**
   * Big Bang Nucleosynthesis file path. The file specifies the predictions for
   * \f$ Y_\mathrm{He} \f$ for given \f$ \omega_b \f$ and \f$ N_\mathrm{eff} \f$.
   */
  std::string sBBN_file = "/bbn/sBBN_2017.dat";

  /*
   *  Thermodynamical quantities
   * */

  /**
   * The initial z for the recfast calculation of the recombination history, e.g. 10^4
   */
  double recfast_z_initial = 1.0e4;
  /**
   * Number of recfast integration steps, e.g. if this is 1.10^4 and the previous one is 10^4, the step will be Delta z = 0.5
   */
  int recfast_Nz0 = 20000;
  /**
   * If there is interacting DM, we want the thermodynamics table to
   * start at a much larger z, in order to capture the possible
   * non-trivial behavior of the dark matter interaction rate at early
   * times:
   *
   * - The new initial redshift will be thermo_z_initial_idm_dr
   *
   * - the highest redhsift will be sampled with thermo_Nz1_idm_dr values, and the step will be
   * Delta z = (thermo_z_initial_idm_dr-recfast_z_initial)/thermo_Nz1_idm_dr
   * For instance, if the previous value is 10^9 and this value is 10^4, then Delta z simeq 10^5
   *
   * - But the first interval after recfast_z_initial will be better
   * sampled with thermo_Nz2_idm_dr values, in order to ensure a smoother
   * transition from a small step to a large step. The intermediate
   * stepsize will then be
   * Delta z = (thermo_z_initial_idm_dr-recfast_z_initial)/thermo_Nz1_idm_dr/thermo_Nz1_idm_dr.
   * For instance, if the three values are (10^9, 10^4, 10^2), then the intermediate timestep is Delta z simeq 10^3
  */
  double thermo_z_initial_idm_dr = 1.0e9;
  int thermo_Nz1_idm_dr          = 10000;
  int thermo_Nz2_idm_dr          = 100;
  /**
   * Tolerance for thermodynamical integration. RECFAST now uses the shared evolver
   * interface with dense output, so this is a true ODE tolerance rather than a
   * legacy per-redshift-bin RK correction tolerance.
   */
  double tol_thermo_integration = 1.0e-7;
  /*
   * Recfast 1.4 switch parameters
   */
  int recfast_Heswitch =
      6; /**< from recfast 1.4, specifies how accurate the Helium recombination should be handled */
  double recfast_fudge_He =
      0.86; /**< from recfast 1.4, fugde factor for Peeble's equation coefficient of Helium */

  /*
   * Recfast 1.5 parameters
   */
  int recfast_Hswitch =
      1; /**< from recfast 1.5, specifies how accurate the Hydrogen recombination should be handled */
  double recfast_fudge_H =
      1.14; /**< from recfast 1.4, fudge factor for Peeble's equation coeffient of Hydrogen */
  double recfast_delta_fudge_H =
      -0.015; /**< from recfast 1.5.2, increasing Hydrogen fudge factor if Hswitch is enabled */
  double recfast_AGauss1 = -0.14; /**< from recfast 1.5, Gaussian Peeble prefactor fit, amplitude */
  double recfast_AGauss2 =
      0.079; /**< from recfast 1.5.2, Gaussian Peeble prefactor fit, amplitude */
  double recfast_zGauss1 = 7.28; /**< from recfast 1.5, Gaussian Peeble prefactor fit, center */
  double recfast_zGauss2 = 6.73; /**< from recfast 1.5.2, Gaussian Peeble prefactor fit, center */
  double recfast_wGauss1 = 0.18; /**< from recfast 1.5, Gaussian Peeble prefactor fit, width */
  double recfast_wGauss2 = 0.33; /**< from recfast 1.5, Gaussian Peeble prefactor fit, width */

  double recfast_z_He_1 = 8000.0; /**< from recfast 1.4, Starting value of Helium recombination 1 */
  double recfast_delta_z_He_1 =
      50.0; /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
  double recfast_z_He_2 = 5000.0; /**< from recfast 1.4, Ending value of Helium recombination 1 */
  double recfast_delta_z_He_2 =
      100.0; /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
  double recfast_z_He_3 = 3500.0; /**< from recfast 1.4, Starting value of Helium recombination 2 */
  double recfast_delta_z_He_3 =
      50.0; /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
  double recfast_x_He0_trigger =
      0.995; /**< Switch for Helium full calculation during reco, raised from 0.99 to 0.995 for smoother Helium */
  double recfast_x_He0_trigger_delta =
      0.05; /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
  double recfast_x_H0_trigger =
      0.995; /**< Switch for Hydrogen full calculation during reco, raised from 0.99 to 0.995 for smoother Hydrogen */
  double recfast_x_H0_trigger2 =
      0.995; /**< Switch for Hydrogen full calculation during reco, for changing Hydrogen flag, raised from 0.98 to same as previous one for smoother Hydrogen */
  double recfast_x_H0_trigger_delta =
      0.05; /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */

  double recfast_H_frac =
      1.0e-3; /**< from recfast 1.4, specifies the time at which the temperature evolution is calculated by the more precise equation */

  double reionization_z_start_max = 50.0;   /**< Maximum starting value in z for reionization */
  double reionization_sampling    = 5.0e-2; /**< Sampling density in z during reionization */
  double reionization_optical_depth_tol =
      1.0e-4; /**< Relative tolerance on finding the user-given optical depth of reionization given a certain redshift of reionization */
  double reionization_start_factor =
      8.0; /**< Searching optical depth corresponding to the redshift is started from an initial offset beyond z_reionization_start, multiplied by reionization_width */

  int thermo_rate_smoothing_radius =
      50; /**< Smoothing in redshift of the variation rate of \f$ \exp(-\kappa) \f$, g, and \f$ \frac{dg}{d\tau} \f$ that is used as a timescale afterwards */

  std::string hyrec_Alpha_inf_file =
      "/hyrec/Alpha_inf.dat"; /**< File containing the alpha parameter of hyrec */
  std::string hyrec_R_inf_file =
      "/hyrec/R_inf.dat"; /**< File containing the R_inf parameter of hyrec */
  std::string hyrec_two_photon_tables_file =
      "/hyrec/two_photon_tables.dat"; /**< File containing the two-photon interaction parameter of hyrec */

  double k_min_tau0 =
      0.1; /**< number defining k_min for the computation of Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */

  double k_max_tau0_over_l_max =
      1.8; /**< number defining k_max for the computation of Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two. Since v3.2.2, the separate full-Limber grid keeps CMB lensing accurate at high l, allowing the standard transfer grid to stop at 1.8 instead of 2.4. */
  double k_step_sub =
      0.05; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_step_super =
      0.002; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */
  double k_step_transition =
      0.2; /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */
  double k_step_super_reduction =
      0.1; /**< the step k_step_super is reduced by this amount in the k-->0 limit (below scale of Hubble and/or curvature radius) */

  double k_per_decade_for_pk =
      10.0; /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade outside the BAO region*/

  double idmdr_boost_k_per_decade_for_pk =
      1.0; /**< boost factor for the case of DAO in idm-idr models */

  double k_per_decade_for_bao =
      70.0; /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)*/

  double k_bao_center =
      3.0; /**< in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak) */

  double k_bao_width =
      4.0; /**< in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak) */

  double start_small_k_at_tau_c_over_tau_h =
      0.0015; /**< largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_largek_at_tau_c_over_tau_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time. */

  double start_large_k_at_tau_h_over_tau_k =
      0.07; /**< largest wavelengths start being sampled when mode is sufficiently outside Hubble scale. This is quantified in terms of the ratio of hubble time scale to wavenumber time scale, \f$ \tau_h/\tau_k \f$ which is roughly equal to (k*tau). Start when this ratio equals start_large_k_at_tau_k_over_tau_h. Decrease this value to start integrating the wavenumbers earlier in time. */

  /**
   * when to switch off tight-coupling approximation: first condition:
   * \f$ \tau_c/\tau_H \f$ > tight_coupling_trigger_tau_c_over_tau_h.
   * Decrease this value to switch off earlier in time.  If this
   * number is larger than start_sources_at_tau_c_over_tau_h, the code
   * returns an error, because the source computation requires
   * tight-coupling to be switched off.
   */
  double tight_coupling_trigger_tau_c_over_tau_h = 0.015;

  /**
   * when to switch off tight-coupling approximation:
   * second condition: \f$ \tau_c/\tau_k \equiv k \tau_c \f$ <
   * tight_coupling_trigger_tau_c_over_tau_k.
   * Decrease this value to switch off earlier in time.
   */
  double tight_coupling_trigger_tau_c_over_tau_k = 0.01;

  double start_sources_at_tau_c_over_tau_h =
      0.008; /**< sources start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_sources_at_tau_c_over_tau_h equals this ratio. Decrease this value to start sampling the sources earlier in time. */

  int tight_coupling_approximation = static_cast<int>(
      tca_method::compromise_CLASS); /**< method for tight coupling approximation */

  double idm_dr_tight_coupling_trigger_tau_c_over_tau_k =
      0.01; /**< when to switch off the dark-tight-coupling approximation, first condition (see normal tca for full definition) */
  double idm_dr_tight_coupling_trigger_tau_c_over_tau_h =
      0.015; /**< when to switch off the dark-tight-coupling approximation, second condition (see normal tca for full definition) */

  double idm_drmd_tight_coupling_trigger_G_over_aH =
      100000; /**< when to switch off the dark-tight-coupling approximation in DRMD, should be larger than at least 100 (currently set to a very high number as the code runs perfectly fine without the approximation.) */
  int l_max_g =
      12; /**< number of momenta in Boltzmann hierarchy for photon temperature (scalar), at least 4 */
  int l_max_pol_g =
      10; /**< number of momenta in Boltzmann hierarchy for photon polarization (scalar), at least 4 */
  int l_max_dr =
      17; /**< number of momenta in Boltzmann hierarchy for decay radiation, at least 4 */
  int l_max_dr_col =
      17; /**< number of collision terms in Boltzmann hierarchy for decay radiation, at least 2 */
  int l_max_dncdm_col =
      17; /**< max multipole receiving the DNCDM inverse-decay collision term; must be <= l_max_ncdm */
  int l_max_ur =
      17; /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_idr =
      17; /**< number of momenta in Boltzmann hierarchy for interacting dark radiation */
  int l_max_ncdm =
      17; /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_g_ten =
      5; /**< number of momenta in Boltzmann hierarchy for photon temperature (tensor), at least 4 */
  int l_max_pol_g_ten =
      5; /**< number of momenta in Boltzmann hierarchy for photon polarization (tensor), at least 4 */

  double curvature_ini = 1.0; /**< initial condition for curvature for adiabatic */
  double entropy_ini   = 1.0; /**< initial condition for entropy perturbation for isocurvature */
  double gw_ini        = 1.0; /**< initial condition for tensor metric perturbation h */

  /**
   * default step \f$ d \tau \f$ in perturbation integration, in units of the timescale involved in the equations (usually, the min of \f$ 1/k \f$, \f$ 1/aH \f$, \f$ 1/\dot{\kappa} \f$)
   */
  double perturb_integration_stepsize = 0.5;

  /**
   * default step \f$ d \tau \f$ for sampling the source function, in units of the timescale involved in the sources: \f$ (\dot{\kappa}- \ddot{\kappa}/\dot{\kappa})^{-1} \f$
   */
  double perturb_sampling_stepsize = 0.1;
  /**
   * Age fraction above which source sampling is twice as fine. This improves
   * the low-l CMB lensing line-of-sight integral; 1.0 disables the boost.
   */
  double perturbations_sampling_boost_above_age_fraction = 0.9;

  /**
   * control parameter for the precision of the perturbation integration,
   * IMPORTANT FOR SETTING THE STEPSIZE OF NDF15
   */
  double tol_perturb_integration = 1.0e-5;

  /**
   * cutoff relevant for controlling stiffness in the PPF scheme. It is
   * neccessary for the Runge-Kutta evolver, but not for ndf15. However,
   * the approximation is excellent for a cutoff value of 1000, so we
   * leave it on for both evolvers. (CAMB uses a cutoff value of 30.)
   */
  double c_gamma_k_H_square_max = 1.0e3;

  /**
   * precision with which the code should determine (by bisection) the
   * times at which sources start being sampled, and at which
   * approximations must be switched on/off (units of Mpc)
   */
  double tol_tau_approx = 1.0e-10;

  /**
   * method for switching off photon perturbations
   */
  int radiation_streaming_approximation = static_cast<int>(rsa_method::rsa_MD_with_reio);

  /**
   * when to switch off photon perturbations, ie when to switch
   * on photon free-streaming approximation (keep density and thtau, set
   * shear and higher momenta to zero):
   * first condition: \f$ k \tau \f$ > radiation_streaming_trigger_tau_h_over_tau_k
   */
  double radiation_streaming_trigger_tau_over_tau_k = 45.0;

  /**
   * when to switch off photon perturbations, ie when to switch
   * on photon free-streaming approximation (keep density and theta, set
   * shear and higher momenta to zero):
   * second condition:
   */
  double radiation_streaming_trigger_tau_c_over_tau = 5.0;

  int idr_streaming_approximation = static_cast<int>(
      rsa_idr_method::rsa_idr_none); /**< method for dark radiation free-streaming approximation */
  double idr_streaming_trigger_tau_over_tau_k =
      50.0; /**< when to switch on dark radiation (idr) free-streaming approximation, first condition */
  double idr_streaming_trigger_tau_c_over_tau =
      10.0; /**< when to switch on dark radiation (idr) free-streaming approximation, second condition */

  int ur_fluid_approximation = static_cast<int>(
      ufa_method::ufa_CLASS); /**< method for ultra relativistic fluid approximation */

  /**
   * when to switch off ur (massless neutrinos / ultra-relativistic
   * relics) fluid approximation
   */
  double ur_fluid_trigger_tau_over_tau_k = 30.0;

  int ncdm_fluid_approximation = static_cast<int>(
      ncdmfa_method::ncdmfa_CLASS); /**< method for non-cold dark matter fluid approximation */

  /**
   * when to switch off ncdm (massive neutrinos / non-cold
   * relics) fluid approximation
   */
  double ncdm_fluid_trigger_tau_over_tau_k = 31.0;

  /**
   * whether CMB source functions can be approximated as zero when
   * visibility function g(tau) is tiny
   */
  double neglect_CMB_sources_below_visibility = 1.0e-3;

  /**
   * The type of evolver to use: options are ndf15, rk, rkdp45 or etd.
   *
   * PER-MODULE, because the three callers do not want the same integrator. The
   * decaying-NCDM inverse-decay sector is the case that forced this: its
   * BACKGROUND Jacobian is genuinely dense (every parent momentum bin couples to
   * every daughter bin through the collision integral, so numjac abandons
   * sparsity and ndf15 costs ~84x more than rkdp45 for answers agreeing to 4e-9
   * of the energy budget), while its PERTURBATION Jacobian stays sparse and
   * stiffness-limited, where an implicit method is the natural choice. With one
   * global setting that model cannot express its own best configuration.
   *
   * Input is an XOR: either the blanket `evolver` key, or any subset of
   * `evolver_background` / `evolver_thermodynamics` / `evolver_perturbations`
   * (unset members keep the default). Mixing the two is rejected rather than
   * silently resolved by a precedence rule -- see precision::parse.
   */
  evolver_type evolver_background     = evolver_type::ndf15;
  evolver_type evolver_thermodynamics = evolver_type::ndf15;
  evolver_type evolver_perturbations  = evolver_type::ndf15;

  /*
   * Primordial parameters
   * */

  double k_per_decade_primordial =
      10.0; /**< logarithmic sampling for primordial spectra (number of points per decade in k space) */

  double primordial_inflation_ratio_min =
      100.0; /**< for each k, start following wavenumber when aH = k/primordial_inflation_ratio_min */
  double primordial_inflation_ratio_max =
      1.0 /
      50.0; /**< for each k, stop following wavenumber, at the latest, when aH = k/primordial_inflation_ratio_max */
  int primordial_inflation_phi_ini_maxit =
      10000; /**< maximum number of iteration when searching a suitable initial field value phi_ini (value reached when no long-enough slow-roll period before the pivot scale) */
  double primordial_inflation_pt_stepsize =
      0.01; /**< controls the integration timestep for inflaton perturbations */
  double primordial_inflation_bg_stepsize =
      0.005; /**< controls the integration timestep for inflaton background */
  double primordial_inflation_tol_integration =
      1.0e-3; /**< controls the precision of the ODE integration during inflation */
  double primordial_inflation_attractor_precision_pivot =
      0.001; /**< targeted precision when searching attractor solution near phi_pivot */
  double primordial_inflation_attractor_precision_initial =
      0.1; /**< targeted precision when searching attractor solution near phi_ini */
  int primordial_inflation_attractor_maxit =
      10; /**< maximum number of iteration when searching attractor solution */
  double primordial_inflation_tol_curvature =
      1.0e-3; /**< for each k, stop following wavenumber, at the latest, when curvature perturbation R is stable up to to this tolerance */
  double primordial_inflation_aH_ini_target =
      0.9; /**< control the step size in the search for a suitable initial field value */
  double primordial_inflation_end_dphi =
      1.0e-10; /**< first bracketing width, when trying to bracket the value phi_end at which inflation ends naturally */
  double primordial_inflation_end_logstep =
      10.0; /**< logarithmic step for updating the bracketing width, when trying to bracket the value phi_end at which inflation ends naturally */
  double primordial_inflation_small_epsilon =
      0.1; /**< value of slow-roll parameter epsilon used to define a field value phi_end close to the end of inflation (doesn't need to be exactly at the end): epsilon(phi_end)=small_epsilon (should be smaller than one) */
  double primordial_inflation_small_epsilon_tol = 0.01; /**< tolerance in the search for phi_end */
  double primordial_inflation_extra_efolds =
      2.0; /**< a small number of efolds, irrelevant at the end, used in the search for the pivot scale (backward from the end of inflation) */

  /*
   * Transfer function parameters
   * */

  int l_linstep =
      40; /**< factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled */

  double l_logstep =
      1.12; /**< maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point) */

  double hyper_x_min =
      1.0e-5; /**< flat case: lower bound on the smallest value of x at which we sample \f$ \Phi_l^{\nu}(x)\f$ or \f$ j_l(x)\f$ */
  double hyper_sampling_flat =
      8.0; /**< flat case: number of sampled points x per approximate wavelength \f$ 2\pi \f$, should remain >7.5 */
  double hyper_sampling_curved_low_nu =
      7.0; /**< open/closed cases: number of sampled points x per approximate wavelength \f$ 2\pi/\nu\f$, when \f$ \nu \f$ smaller than hyper_nu_sampling_step */
  double hyper_sampling_curved_high_nu =
      3.0; /**< open/closed cases: number of sampled points x per approximate wavelength \f$ 2\pi/\nu\f$, when \f$ \nu \f$ greater than hyper_nu_sampling_step */
  double hyper_nu_sampling_step =
      1000.0; /**< open/closed cases: value of nu at which sampling changes  */
  double hyper_phi_min_abs =
      1.0e-10; /**< small value of Bessel function used in calculation of first point x (\f$ \Phi_l^{\nu}(x) \f$ equals hyper_phi_min_abs) */
  double hyper_x_tol = 1.0e-4; /**< tolerance parameter used to determine first value of x */
  double hyper_flat_approximation_nu =
      4000.0; /**< value of nu below which the flat approximation is used to compute Bessel function */

  double q_linstep = 0.45; /**< asymptotic linear sampling step in q
  // UNHANDLED: space, in units of \f$ 2\pi/r_a(\tau_rec) \f$
  // UNHANDLED: (comoving angular diameter distance to
  // UNHANDLED: recombination), very important for CMB */

  double q_logstep_spline = 170.0; /**< initial logarithmic sampling step in q
  // UNHANDLED: space, in units of \f$ 2\pi/r_a(\tau_{rec})\f$
  // UNHANDLED: (comoving angular diameter distance to
  // UNHANDLED: recombination), very important for CMB and LSS */

  double q_logstep_open = 6.0; /**< in open models, the value of
  // UNHANDLED: q_logstep_spline must be decreased
  // UNHANDLED: according to curvature. Increasing
  // UNHANDLED: this number will make the calculation
  // UNHANDLED: more accurate for large positive
  // UNHANDLED: Omega_k */

  double q_logstep_trapzd = 20.0; /**< initial logarithmic sampling step in q
  // UNHANDLED: space, in units of \f$ 2\pi/r_a(\tau_{rec}) \f$
  // UNHANDLED: (comoving angular diameter distance to
  // UNHANDLED: recombination), in the case of small
  // UNHANDLED: q's in the closed case, for which one
  // UNHANDLED: must used trapezoidal integration
  // UNHANDLED: instead of spline (the number of q's
  // UNHANDLED: for which this is the case decreases
  // UNHANDLED: with curvature and vanishes in the
  // UNHANDLED: flat limit) */

  double q_numstep_transition = 250.0; /**< number of steps for the transition
  // UNHANDLED: from q_logstep_trapzd steps to
  // UNHANDLED: q_logstep_spline steps (transition
  // UNHANDLED: must be smooth for spline) */

  double q_logstep_limber =
      1.025; /**< logarithmic q-step ratio for the separate full-Limber CMB lensing grid */
  double k_max_limber_over_l_max_scalars =
      0.001; /**< full-Limber perturbation source cutoff k_max/l_max_scalars in 1/Mpc */

  double transfer_neglect_delta_k_S_t0 =
      0.15; /**< for temperature source function T0 of scalar mode, range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero */
  double transfer_neglect_delta_k_S_t1 =
      0.04; /**< same for temperature source function T1 of scalar mode */
  double transfer_neglect_delta_k_S_t2 =
      0.15; /**< same for temperature source function T2 of scalar mode */
  double transfer_neglect_delta_k_S_e =
      0.11; /**< same for polarization source function E of scalar mode */
  double transfer_neglect_delta_k_V_t1 =
      1.0; /**< same for temperature source function T1 of vector mode */
  double transfer_neglect_delta_k_V_t2 =
      1.0; /**< same for temperature source function T2 of vector mode */
  double transfer_neglect_delta_k_V_e =
      1.0; /**< same for polarization source function E of vector mode */
  double transfer_neglect_delta_k_V_b =
      1.0; /**< same for polarization source function B of vector mode */
  double transfer_neglect_delta_k_T_t2 =
      0.2; /**< same for temperature source function T2 of tensor mode */
  double transfer_neglect_delta_k_T_e =
      0.25; /**< same for polarization source function E of tensor mode */
  double transfer_neglect_delta_k_T_b =
      0.1; /**< same for polarization source function B of tensor mode */

  double transfer_neglect_late_source =
      400.0; /**< value of l below which the CMB source functions can be neglected at late time, excepted when there is a Late ISW contribution */

  double l_switch_limber =
      10.; /**< when to use the Limber approximation for project gravitational potential cl's */
  // For density Cl, we recommend not to use the Limber approximation
  // at all, and hence to put here a very large number (e.g. 10000); but
  // if you have wide and smooth selection functions you may wish to
  // use it; then 100 might be OK
  double l_switch_limber_for_nc_local_over_z =
      100.0; /**< when to use the Limber approximation for local number count contributions to cl's (relative to central redshift of each bin) */
  // For terms integrated along the line-of-sight involving spherical
  // Bessel functions (but not their derivatives), Limber
  // approximation works well. High precision can be reached with 2000
  // only. But if you have wide and smooth selection functions you may
  // reduce to e.g. 30.
  double l_switch_limber_for_nc_los_over_z =
      30.0; /**< when to use the Limber approximation for number count contributions to cl's integrated along the line-of-sight (relative to central redshift of each bin) */

  double selection_cut_at_sigma =
      5.0; /**< in sigma units, where to cut gaussian selection functions */
  double selection_sampling =
      50.0; /**< controls sampling of integral over time when selection functions vary quicker than Bessel functions. Increase for better sampling. */
  double selection_sampling_bessel =
      20.0; /**< controls sampling of integral over time when selection functions vary slower than Bessel functions. Increase for better sampling. IMPORTANT for lensed contributions. */
  double selection_sampling_bessel_los =
      20.0; /**< controls sampling of integral over time when selection functions vary slower than Bessel functions. This parameter is specific to number counts contributions to Cl integrated along the line of sight. Increase for better sampling */
  double selection_tophat_edge =
      0.1; /**< controls how smooth are the edge of top-hat window function (<<1 for very sharp, 0.1 for sharp) */

  /*
   * Nonlinear module precision parameters
   * */

  double sigma_k_per_decade =
      80.; /**< logarithmic stepsize controlling the precision of integrals for sigma(R,k) and similar quantitites */

  double nonlinear_min_k_max = 5.0; /**< when
  // UNHANDLED: using an algorithm to compute nonlinear
  // UNHANDLED: corrections, like halofit or hmcode,
  // UNHANDLED: k_max must be at least equal to this
  // UNHANDLED: value. Calculations are done internally
  // UNHANDLED: until this k_max, but the P(k,z) output
  // UNHANDLED: is still controlled by P_k_max_1/Mpc or
  // UNHANDLED: P_k_max_h/Mpc even if they are
  // UNHANDLED: smaller */

  /** parameters relevant for HALOFIT computation */

  double halofit_min_k_nonlinear =
      1.0e-4; /**< value of k in 1/Mpc below which non-linear corrections will be neglected */

  double halofit_min_k_max = 5.0; /**< DEPRECATED: should use instead nonlinear_min_k_max */

  double halofit_k_per_decade = 80.0; /**< halofit needs to evalute integrals
  // UNHANDLED: (linear power spectrum times some
  // UNHANDLED: kernels). They are sampled using
  // UNHANDLED: this logarithmic step size. */

  double halofit_sigma_precision = 0.05; /**< a smaller value will lead to a
  // UNHANDLED: more precise halofit result at the *highest*
  // UNHANDLED: redshift at which halofit can make computations,
  // UNHANDLED: at the expense of requiring a larger k_max; but
  // UNHANDLED: this parameter is not relevant for the
  // UNHANDLED: precision on P_nl(k,z) at other redshifts, so
  // UNHANDLED: there is normally no need to change it */

  double halofit_tol_sigma = 1.0e-6; /**< tolerance required on sigma(R) when
  // UNHANDLED: matching the condition sigma(R_nl)=1,
  // UNHANDLED: whcih defines the wavenumber of
  // UNHANDLED: non-linearity, k_nl=1./R_nl */

  double pk_eq_z_max = 5.0;    /**< Maximum z for the pk_eq method */
  double pk_eq_tol   = 1.0e-7; /**< Tolerance on the pk_eq method for finding the pk */

  /** Parameters relevant for HMcode computation */

  double hmcode_max_k_extra = 1.e6; /**< parameter specifying the maximum k value for
  // UNHANDLED: the extrapolation of the linear power spectrum
  // UNHANDLED: (needed for the sigma computation) */

  double hmcode_min_k_max = 5.; /**< DEPRECATED: should use instead nonlinear_min_k_max */

  double hmcode_tol_sigma = 1.e-6; /**< tolerance required on sigma(R) when matching the
  // UNHANDLED: condition sigma(R_nl)=1, which defines the wavenumber
  // UNHANDLED: of non-linearity, k_nl=1./R_nl */

  /**
   * parameters controlling stepsize and min/max r & a values for
   * sigma(r) & grow table
   */
  int n_hmcode_tables      = 64;
  double rmin_for_sigtab   = 1.e-5;
  double rmax_for_sigtab   = 1.e3;
  double ainit_for_growtab = 1.e-3;
  double amax_for_growtab  = 1.;

  /**
   * parameters controlling stepsize and min/max halomass values for the
   * 1-halo-power integral
   */
  int nsteps_for_p1h_integral  = 256;
  double mmin_for_p1h_integral = 1.e3;
  double mmax_for_p1h_integral = 1.e18;

  /*
   * Lensing precision parameters
   * */

  int accurate_lensing =
      0; /**< switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles */
  int num_mu_minus_lmax =
      70;                /**< difference between num_mu and l_max, increase for more precision */
  int delta_l_max = 500; /**< difference between l_max in unlensed and lensed spectra */
  double tol_gauss_legendre =
      DBL_EPSILON; /**< tolerance with which quadrature points are found: must be very small for an accurate integration (if not entered manually, set automatically to match machine precision) */

  /** @name - general precision parameters */
  //@{
  double smallest_allowed_variation =
      DBL_EPSILON; /**< machine-dependent, assigned automatically by the code */
  //@}

  /** Prepend the runtime class_dir to the relative data-file path defaults
   *  (sBBN_file, hyrec_*). Call once after class_dir is set and before parse(),
   *  so a user-supplied absolute override in parse() still replaces verbatim. */
  void ResolveDataPaths();

  /** Parse precision parameters from a configuration file. */
  void parse(const FileContent& fc);
};

#endif  // CLASS_PRECISION_H
