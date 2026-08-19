#ifndef BACKGROUND_MODULE_H
#define BACKGROUND_MODULE_H

#include "base_module.h"
#include "input_module.h"

class BackgroundModule : public BaseModule {
 public:
  /** Look up a species by key (CDM, Lambda, Fluid, ScalarField, UR, IDM_DR_IDR,
   *  IDM_DRMD_IDR_DRMD, DCDM_DR, NCDM, DNCDM_DR, ...) and return its present-day
   *  Omega0; returns 0 if the species is absent.  For composite sub-species
   *  (e.g. IDR inside IDM_DR_IDR) use the dedicated sub_key form. */
  double GetOmega0Species(const std::string& key) const;
  BackgroundModule(InputModulePtr input_module);
  ~BackgroundModule();
  void background_output_titles(std::string& titles) const;
  void background_output_data(int number_of_titles, double* data) const;
  void background_at_tau(
      double tau, short return_format, short inter_mode, int* last_index, double* pvecback) const;
  void background_tau_of_z(double z, double* tau) const;
  void background_w_fld(double a,
                        double* w_fld,
                        double* dw_over_da_fld,
                        double* integral_fld) const;
  void background_idm_drmd(
      double a, double rho_idm_over_rho_idr, double* Rint, double* csp2, double* Gint) const;

  /** Sum of GetOmega0() across all NCDM-family species (NCDM, NCDMInteracting,
   *  and the DNCDMSpecies child of every DNCDM_DR_Species composite). */
  double GetOmega0NcdmTot() const;

  /** Count of DR-bearing composites: 1 if DCDM_DR is present, plus the
   *  number of DNCDM_DR_Species composites in all_species_. */
  int GetNDecayDr() const;

  /** Generic look-up: find species by all_species_ key, then call
   *  GetParam(param).  Returns 0.0 if either is absent.  Used by the
   *  Python wrapper to read species-specific quantities without
   *  downcasting. */
  double GetSpeciesParam(const std::string& key, const std::string& param) const;

  /** @name - all indices for the vector of background (=bg) quantities stored in table */

  //@{

  int index_bg_a_;       /**< scale factor */
  int index_bg_H_;       /**< Hubble parameter in \f$Mpc^{-1}\f$ */
  int index_bg_H_prime_; /**< its derivative w.r.t. conformal time */

  /* end of vector in short format, now quantities in normal format */

  int index_bg_rho_tot_;     /**< Total density */
  int index_bg_p_tot_;       /**< Total pressure */
  int index_bg_p_tot_prime_; /**< Conf. time derivative of total pressure */

  int index_bg_Omega_r_; /**< relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) */

  /* end of vector in normal format, now quantities in long format */

  int index_bg_rho_crit_; /**< critical density */
  int index_bg_Omega_m_; /**< non-relativistic density fraction (\f$ \Omega_b + \Omega_cdm + \Omega_{\nu nr} \f$) */
  int index_bg_conf_distance_; /**< conformal distance (from us) in Mpc */
  int index_bg_ang_distance_;  /**< angular diameter distance in Mpc */
  int index_bg_lum_distance_;  /**< luminosity distance in Mpc */
  int index_bg_time_;          /**< proper (cosmological) time in Mpc */
  int index_bg_rs_;            /**< comoving sound horizon in Mpc */

  int index_bg_D_; /**< scale independent growth factor D(a) for CDM perturbations */
  int index_bg_f_; /**< corresponding velocity growth factor [dlnD]/[dln a] */

  int bg_size_short_;  /**< size of background vector in the "short format" */
  int bg_size_normal_; /**< size of background vector in the "normal format" */
  int bg_size_;        /**< size of background vector in the "long format" */

  //@}
  int bt_size_; /**< number of lines (i.e. time-steps) in the array */
  std::vector<double>
      tau_table_; /**< vector tau_table[index_tau] with values of \f$ \tau \f$ (conformal time) */
  std::vector<double>
      background_table_; /**< table background_table[index_tau*bg_size_+pba->index_bg] with all other quantities (array of size bg_size*bt_size) **/

  double conformal_age_; /**< conformal age in Mpc */
  double
      Neff_; /**< so-called "effective neutrino number", computed at earliest time in interpolation table */
  double a_eq_;     /**< scale factor at radiation/matter equality */
  double H_eq_;     /**< Hubble rate at radiation/matter equality [Mpc^-1] */
  double Omega0_m_; /**< total non-relativistic matter today */
  double
      Omega0_de_; /**< total dark energy density today, currently defined as 1 - Omega0_m - Omega0_r - Omega0_k */

  double age_;         /**< age in Gyears */
  double Omega0_r_;    /**< total ultra-relativistic radiation today */
  double z_eq_;        /**< redshift at radiation/matter equality */
  double tau_eq_;      /**< conformal time at radiation/matter equality [Mpc] */
  double Omega0_dcdm_; /**< \f$ \Omega_{0 dcdm} \f$: decaying cold dark matter */
  double Omega0_dr_;   /**< \f$ \Omega_{0 dr} \f$: decay radiation */

  double Omega0_idr_drmd_;
  double Omega0_idm_drmd_;
  double Omega0_idm_;

  // DRMD diagnostics exposed to the Python wrapper. The data now lives on
  // IDM_DRMD_IDR_DRMD_Species (#319); these delegate to it (neutral values if
  // no DRMD species is present).
  double f_idr_drmd() const;
  double z_dec_drmd() const;

  /** True while background_functions is filling a STORED row of background_table_
   *  (return_format == long_info), false during the trial evaluations the evolver
   *  makes for derivatives, initial conditions and rejected steps.
   *
   *  A species self-check must consult this before reporting anything. Trial states
   *  are allowed to be unphysical -- the evolver's whole job is to propose a step,
   *  measure it and throw it away -- so a check that fires on them reports failures
   *  that never reached the solution. That is not hypothetical: the first version of
   *  DNCDMSpecies's floor guard did exactly that, and its false positives said the
   *  lasing regime was broken when the accepted solution never came within 75
   *  decades of the floor.
   *
   *  Safe as plain state because the background is solved on ONE thread; the
   *  perturbation modules read the finished table through background_at_tau and
   *  never re-enter background_functions. */
  bool StoringBackgroundTable() const {
    return storing_background_table_;
  }

 private:
  void background_functions(double* pvecback_B, short return_format, double* pvecback);
  bool storing_background_table_ = false;

  void background_init();
  void background_indices();
  void background_solve_evolver();
  void background_initial_conditions(double* pvecback, double* pvecback_integration);
  void background_find_equality();
  void background_derivs_member(double z, double* y, double* dy, void* parameters_and_workspace);
  void background_derivs_loga_member(double loga,
                                     double* y,
                                     double* dy,
                                     void* parameters_and_workspace);
  static void background_derivs_loga(double loga,
                                     double* y,
                                     double* dy,
                                     void* parameters_and_workspace);
  /** d(dy_i)/d(y_i) of the loga RHS, for the exponential evolver. Same shape and
   *  units as background_derivs_loga, including the 1/(a H) conversion. */
  void background_derivs_diagonal_loga_member(double loga,
                                              double* y,
                                              double* diag,
                                              void* parameters_and_workspace);
  static void background_derivs_diagonal_loga(double loga,
                                              double* y,
                                              double* diag,
                                              void* parameters_and_workspace);
  static void background_timescale(double loga, void* parameters_and_workspace, double* timescale);
  void background_add_line_to_bg_table_member(
      double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace);
  static void background_add_line_to_bg_table(
      double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace);
  void background_output_budget();
  static void background_print_variables(double loga,
                                         double* y,
                                         double* dy,
                                         void* parameters_and_workspace);

  /** @name - all indices for the vector of background quantities to be integrated (=bi)
   *
   * Most background quantities can be immediately inferred from the
   * scale factor. Only few of them require an integration with
   * respect to conformal time (in the minimal case, only one quantity needs to
   * be integrated with time: the scale factor, using the Friedmann
   * equation). These indices refer to the vector of
   * quantities to be integrated with time.
   * {B} quantities are needed by background_functions() while {C} quantities are not.
   */

  //@{

  int index_bi_a_; /**< {B} scale factor */

  int index_bi_time_;    /**< {C} proper (cosmological) time in Mpc */
  int index_bi_rs_;      /**< {C} sound horizon */
  int index_bi_tau_;     /**< {C} conformal time in Mpc */
  int index_bi_D_;       /**< {C} scale independent growth factor D(a) for CDM perturbations. */
  int index_bi_D_prime_; /**< {C} D satisfies \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */

  int bi_B_size_; /**< Number of {B} parameters */
  int bi_size_;   /**< Number of {B}+{C} parameters */

  //@}

  /** @name - background interpolation tables */

  //@{

  std::vector<double>
      z_table_; /**< vector z_table[index_tau] with values of \f$ z \f$ (redshift) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  std::vector<double>
      d2tau_dz2_table_; /**< vector d2tau_dz2_table[index_tau] with values of \f$ d^2 \tau / dz^2 \f$ (conformal time) */
  std::vector<double>
      d2background_dtau2_table_; /**< table d2background_dtau2_table[index_tau*bg_size_+pba->index_bg] with values of \f$ d^2 b_i / d\tau^2 \f$ (conformal time) */

  //@}
};

/**
 * temporary parameters and workspace passed to the background_derivs_loga / background_derivs evolver functions
 */

struct background_parameters_and_workspace {
  background_parameters_and_workspace(BackgroundModule* b_m) : background_module(b_m) {}
  BackgroundModule* const background_module;
  double* pvecback;
};

#endif  //BACKGROUND_MODULE_H
