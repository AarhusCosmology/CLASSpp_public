#ifndef PERTURBATIONS_MODULE_H
#define PERTURBATIONS_MODULE_H

#include <vector>

#include "base_module.h"
#include "input_module.h"

class FluidSpecies;

class PerturbationsModule : public BaseModule {
 public:
  // The single PPF fluid, or nullptr. PPF is a modified-gravity closure defined
  // relative to the whole universe (intrinsically singular), so the module owns
  // it; non-PPF fluids are ordinary species. Resolved once in the ctor.
  FluidSpecies* ppf_fluid() const {
    return ppf_fluid_;
  }

  PerturbationsModule(InputModulePtr input_module,
                      BackgroundModulePtr background_module,
                      ThermodynamicsModulePtr thermodynamics_module);
  ~PerturbationsModule();
  void perturb_output_data(enum file_format output_format,
                           double z,
                           int number_of_titles,
                           double* data) const;
  void perturb_output_titles(enum file_format output_format, std::string& titles) const;
  void perturb_output_firstline_and_ic_suffix(int index_ic,
                                              std::string& first_line,
                                              std::string& ic_suffix) const;

  /** Accessors that expose protected BaseModule pointers to species code. */
  const background* GetBackground() const noexcept {
    return pba;
  }
  const perturbs* GetPerturbs() const noexcept {
    return ppt;
  }
  const precision* GetPrecision() const noexcept {
    return ppr;
  }
  BackgroundModulePtr GetBackgroundModule() const {
    return background_module_;
  }
  ThermodynamicsModulePtr GetThermodynamicsModule() const {
    return thermodynamics_module_;
  }

  /** Write a value into the source table at (mode, ic, type, tau, k). */
  void SetSourceValue(
      int index_md, int index_ic, int index_tp, int index_tau, int index_k, double value) {
    sources_[index_md][index_ic * tp_size_[index_md] + index_tp]
            [index_tau * k_size_[index_md] + index_k] = value;
  }

  /** @name - indices running on modes (scalar, vector, tensor) */
  //@{
  int index_md_scalars_; /**< index value for scalars */
  int index_md_tensors_; /**< index value for tensors */
  int index_md_vectors_; /**< index value for vectors */
  int md_size_;          /**< number of modes included in computation */
  //@}
  /** @name - indices running on initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) */
  //@{
  int index_ic_ad_  = -1; /**< index value for adiabatic */
  int index_ic_cdi_ = -1; /**< index value for CDM isocurvature */
  int index_ic_bi_  = -1; /**< index value for baryon isocurvature */
  int index_ic_nid_ = -1; /**< index value for neutrino density isocurvature */
  int index_ic_niv_ = -1; /**< index value for neutrino velocity isocurvature */
  int index_ic_ten_ = -1; /**< index value for unique possibility for tensors */
  std::vector<int>
      ic_size_; /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  //@}
  /** @name - flags and indices running on types (temperature, polarization, lensing, ...) */

  //@{

  /* remember that the temperature source function includes three
     terms that we call 0,1,2 (since the strategy in class v > 1.7 is
     to avoid the integration by part that would reduce the source to
     a single term) */
  int index_tp_t0_;        /**< index value for temperature (j=0 term) */
  int index_tp_t1_;        /**< index value for temperature (j=1 term) */
  int index_tp_t2_;        /**< index value for temperature (j=2 term) */
  int index_tp_p_;         /**< index value for polarization */
  int index_tp_delta_tot_; /**< index value for total density fluctuation */
  int index_tp_perturbed_recombination_delta_temp_; /**< Gas temperature perturbation */
  int index_tp_perturbed_recombination_delta_chi_;  /**< Inionization fraction perturbation */

  int index_tp_theta_m_;   /**< index value for matter velocity fluctuation */
  int index_tp_theta_cb_;  /**< index value for theta cb */
  int index_tp_theta_tot_; /**< index value for total velocity fluctuation */

  int index_tp_phi_;          /**< index value for metric fluctuation phi */
  int index_tp_phi_prime_;    /**< index value for metric fluctuation phi' */
  int index_tp_phi_plus_psi_; /**< index value for metric fluctuation phi+psi */
  int index_tp_psi_;          /**< index value for metric fluctuation psi */
  int index_tp_h_;            /**< index value for metric fluctuation h */
  int index_tp_h_prime_;      /**< index value for metric fluctuation h' */
  int index_tp_eta_;          /**< index value for metric fluctuation eta */
  int index_tp_eta_prime_;    /**< index value for metric fluctuation eta' */
  int index_tp_H_T_Nb_prime_; /**< index value for metric fluctuation H_T_Nb' */
  int index_tp_k2gamma_Nb_; /**< index value for metric fluctuation gamma times k^2 in Nbody gauge */

  int index_tp_delta_m_;  /**< index value for matter density fluctuation */
  int index_tp_delta_cb_; /**< index value for delta cb */
  std::vector<int>
      tp_size_; /**< number of types tp_size[index_md] included in computation for each mode */

  bool has_source_t_;            /**< do we need source for CMB temperature? */
  bool has_source_p_;            /**< do we need source for CMB polarization? */
  bool has_source_delta_m_;      /**< do we need source for delta of total matter? */
  bool has_source_delta_cb_;     /**< do we ALSO need source for delta of ONLY cdm and baryon? */
  bool has_source_delta_tot_;    /**< do we need source for delta total? */
  bool has_source_theta_m_;      /**< do we need source for theta of total matter? */
  bool has_source_theta_cb_;     /**< do we ALSO need source for theta of ONLY cdm and baryon? */
  bool has_source_theta_tot_;    /**< do we need source for theta total? */
  bool has_source_phi_;          /**< do we need source for metric fluctuation phi? */
  bool has_source_phi_prime_;    /**< do we need source for metric fluctuation phi'? */
  bool has_source_phi_plus_psi_; /**< do we need source for metric fluctuation (phi+psi)? */
  bool has_source_psi_;          /**< do we need source for metric fluctuation psi? */
  bool has_source_h_;            /**< do we need source for metric fluctuation h? */
  bool has_source_h_prime_;      /**< do we need source for metric fluctuation h'? */
  bool has_source_eta_;          /**< do we need source for metric fluctuation eta? */
  bool has_source_eta_prime_;    /**< do we need source for metric fluctuation eta'? */
  bool has_source_H_T_Nb_prime_; /**< do we need source for metric fluctuation H_T_Nb'? */
  bool
      has_source_k2gamma_Nb_; /**< do we need source for metric fluctuation gamma in Nbody gauge? */

  /** @name - arrays storing the evolution of all sources for given k values, passed as k_output_values */
  //@{
  std::vector<int>
      index_k_output_values_; /**< List of indices corresponding to k-values close to k_output_values for each mode.
                                     index_k_output_values[index_md*k_output_values_num+ik]*/

  std::string
      scalar_titles_; /**< _DELIMITER_ separated string of titles for scalar perturbation output files. */
  std::string
      vector_titles_; /**< _DELIMITER_ separated string of titles for vector perturbation output files. */
  std::string
      tensor_titles_; /**< _DELIMITER_ separated string of titles for tensor perturbation output files. */

  std::vector<double> scalar_perturbations_data_
      [_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for scalars */
  std::vector<double> vector_perturbations_data_
      [_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for vectors */
  std::vector<double> tensor_perturbations_data_
      [_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for tensors */

  int size_scalar_perturbation_data_
      [_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of scalar double pointers  */
  int size_vector_perturbation_data_
      [_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of vector double pointers  */
  int size_tensor_perturbation_data_
      [_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of tensor double pointers  */
  //@}

  /** @name - source functions interpolation table */

  //@{

  std::vector<std::vector<std::vector<double>>> sources_; /**< Source interpolation table
                         sources[index_md]
                         [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                         [index_tau * perturbations_module_->k_size_ + index_k] */

  //@}
  /** @name - arrays related to the interpolation table for sources at late times, corresponding to z < z_max_pk (used for Fourier transfer function and spectra output) */
  //@{

  std::vector<double>
      ln_tau_; /**< log of the arrau tau_sampling, covering only the final time range required for the output of
                            Fourier transfer functions (used for interpolations) */
  int ln_tau_size_; /**< number of values in this array */

  //@}

  std::vector<double> tau_sampling_; /**< array of tau values */
  int tau_size_;                     /**< number of values in this array */

  std::vector<int> k_size_cl_;         /**< k_size_cl[index_md] number of k values used
                       for non-CMB \f$ C_l \f$ calculations, requiring a coarse
                       sampling in k-space. */
  std::vector<int> k_size_;            /**< k_size[index_md] = total number of k
                       values, including those needed for P(k) but not
                       for \f$ C_l \f$'s */
  std::vector<std::vector<double>> k_; /**< k[index_md][index_k] = list of values */
  double k_min_;                       /**< minimum value (over all modes) */
  double k_max_;                       /**< maximum value (over all modes) */

 private:
  void perturb_sources_at_tau(
      int index_md, int index_ic, int index_tp, double tau, double* pvecsources) const;
  void perturb_init();
  void perturb_indices_of_perturbs();
  void perturb_timesampling_for_sources();
  void perturb_get_k_list();
  void perturb_workspace_init(int index_md, perturb_workspace* ppw);
  void perturb_solve(int index_md, int index_ic, int index_k, perturb_workspace* ppw);
  void perturb_prepare_k_output();
  void perturb_find_approximation_number(int index_md,
                                         double k,
                                         perturb_workspace* ppw,
                                         double tau_ini,
                                         double tau_end,
                                         int* interval_number,
                                         int* interval_number_of);
  void perturb_find_approximation_switches(int index_md,
                                           double k,
                                           perturb_workspace* ppw,
                                           double tau_ini,
                                           double tau_end,
                                           double precision,
                                           int interval_number,
                                           int* interval_number_of,
                                           double* interval_limit,
                                           int** interval_approx);
  void perturb_vector_init(
      int index_md, int index_ic, double k, double tau, perturb_workspace* ppw, int* pa_old);
  void perturb_initial_conditions(
      int index_md, int index_ic, double k, double tau, perturb_workspace* ppw);
  void perturb_approximations(int index_md, double k, double tau, perturb_workspace* ppw);
  void perturb_einstein(int index_md, double k, double tau, double* y, perturb_workspace* ppw);
  void perturb_total_stress_energy(int index_md, double k, double* y, perturb_workspace* ppw);
  int perturb_timescale_member(double tau, void* parameters_and_workspace, double* timescale);
  static int perturb_timescale(double tau, void* parameters_and_workspace, double* timescale);
  int perturb_sources_member(double tau,
                             double* pvecperturbations,
                             double* pvecderivs,
                             int index_tau,
                             void* parameters_and_workspace);
  static int perturb_sources(double tau,
                             double* pvecperturbations,
                             double* pvecderivs,
                             int index_tau,
                             void* parameters_and_workspace);
  int perturb_print_variables_member(double tau,
                                     double* y,
                                     double* dy,
                                     void* parameters_and_workspace);
  static int perturb_print_variables(double tau,
                                     double* y,
                                     double* dy,
                                     void* parameters_and_workspace);
  int perturb_derivs_member(double tau, double* y, double* dy, void* parameters_and_workspace);
  static int perturb_derivs(double tau, double* y, double* dy, void* parameters_and_workspace);
  void perturb_tca_slip_and_shear(double* y, void* parameters_and_workspace);
  void perturb_rsa_delta_and_theta(
      double k, double* y, double a_prime_over_a, double* pvecthermo, perturb_workspace* ppw);
  void perturb_rsa_idr_delta_and_theta(
      double k, double* y, double a_prime_over_a, double* pvecthermo, perturb_workspace* ppw);

  BackgroundModulePtr background_module_;
  ThermodynamicsModulePtr thermodynamics_module_;
  FluidSpecies* ppf_fluid_ = nullptr;

  short
      evolve_tensor_ur_; /**< will we evolve ur tensor perturbations (either because we have ur species, or we have ncdm species with massless approximation) ? */
  short
      evolve_tensor_ncdm_; /**< will we evolve ncdm tensor perturbations (if we have ncdm species and we use the exact method) ? */

  /** @name - useful flags  */

  //@{

  bool has_cmb_; /**< do we need CMB-related sources (temperature, polarization) ? */
  bool has_lss_; /**< do we need LSS-related sources (lensing potential, ...) ? */

  //@}

  std::vector<std::vector<double*>> late_sources_; /**< Pointer into the source interpolation table
                                late_sources[index_md]
                                            [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                Note: these are non-owning pointers into sources_ memory (alias) */

  std::vector<std::vector<std::vector<double>>>
      ddlate_sources_; /**< Splined source interpolation table with second derivatives with respect to time
                              ddlate_sources[index_md]
                                            [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                            [index_tau * perturbations_module_->k_size_ + index_k] */

  //@}
  int number_of_scalar_titles_; /**< number of titles/columns in scalar perturbation output files */
  int number_of_vector_titles_; /**< number of titles/columns in vector perturbation output files*/
  int number_of_tensor_titles_; /**< number of titles/columns in tensor perturbation output files*/

  /** @name - list of k values for each mode */

  //@{

  std::vector<int> k_size_cmb_; /**< k_size_cmb[index_md] number of k values used
                        for CMB calculations, requiring a fine
                        sampling in k-space */

  //@}
};

/**
 * Structure pointing towards all what the function that perturb_derivs
 * needs to know: fixed input parameters and indices contained in the
 * various structures, workspace, etc.
 */

struct perturb_parameters_and_workspace {
  perturb_parameters_and_workspace(PerturbationsModule* p_m) : perturbations_module(p_m) {}
  PerturbationsModule* const perturbations_module;
  int index_md;           /**< index of mode (scalar/.../vector/tensor) */
  int index_ic;           /**< index of initial condition (adiabatic/isocurvature(s)/...) */
  int index_k;            /**< index of wavenumber */
  double k;               /**< current value of wavenumber in 1/Mpc */
  perturb_workspace* ppw; /**< workspace defined above */
};

#endif  //PERTURBATIONS_MODULE_H
