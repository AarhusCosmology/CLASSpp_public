#ifndef THERMODYNAMICS_MODULE_H
#define THERMODYNAMICS_MODULE_H

#include <vector>

#include "base_module.h"
#include "input_module.h"
#include "recombination_model.h"

class ThermodynamicsModule : public BaseModule {
 public:
  ThermodynamicsModule(InputModulePtr input_module, BackgroundModulePtr background_module);
  ~ThermodynamicsModule();
  void thermodynamics_output_titles(std::string& titles) const;
  void thermodynamics_output_data(int number_of_titles, double* data) const;
  void thermodynamics_at_z(
      double z, short inter_mode, int* last_index, double* pvecback, double* pvecthermo) const;

  const thermo* GetThermodynamics() const {
    return pth;
  }

  double
      tau_ini_; /**< initial conformal time at which thermodynamical variables have been be integrated */
  double YHe_;
  /**
   *@name - some flags needed for thermodynamics functions
   */
  //@{
  short
      inter_normal_; /**< flag for calling thermodynamics_at_z and find position in interpolation table normally */
  short
      inter_closeby_; /**< flag for calling thermodynamics_at_z and find position in interpolation table starting from previous position in previous call */
  //@}

  int tt_size_;  /**< number of lines (redshift steps) in the tables */
  double z_rec_; /**< z at which the visibility reaches its maximum (= recombination redshift) */
  int th_size_;  /**< size of thermodynamics vector */
  double
      tau_rec_; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double
      angular_rescaling_; /**< [ratio ra_rec / (tau0-tau_rec)]: gives CMB rescaling in angular space relative to flat model (=1 for curvature K=0) */
  double
      tau_free_streaming_; /**< minimum value of tau at which free-streaming approximation can be switched on */
  double
      tau_idr_free_streaming_; /**< trigger for dark radiation free streaming approximation (idm-idr) */
  double
      tau_cut_; /**< at at which the visibility goes below a fixed fraction of the maximum visibility, used for an approximation in perturbation module */

  double tau_reionization_;
  double z_reionization_;
  double n_e_; /**< total number density of electrons today (free or not) */

  /** @name - characteristic quantities like redshift, conformal time and sound horizon at recombination */

  //@{

  double ds_rec_; /**< physical sound horizon at recombination */
  double da_rec_; /**< physical angular diameter distance to recombination */
  double rd_rec_; /**< comoving photon damping scale at recombination */
  double rs_rec_; /**< comoving sound horizon at recombination */
  double ra_rec_; /**< conformal angular diameter distance to recombination */

  double rs_star_;  /**< comoving sound horizon at z_star */
  double ra_star_;  /**< conformal angular diameter distance to z_star */
  double z_star_;   /**< redshift at which photon optical depth crosses one */
  double tau_star_; /**< conformal time at which photon optical depth crosses one */
  double ds_star_;  /**< physical sound horizon at z_star */
  double da_star_;  /**< physical angular diameter distance to z_star */
  double rd_star_;  /**< comoving photon damping scale at z_star */

  double z_d_;   /**< baryon drag redshift */
  double tau_d_; /**< baryon drag time */
  double rs_d_;  /**< comoving sound horizon at baryon drag */
  double ds_d_;  /**< physical sound horizon at baryon drag */

  //@}

  /** @name - all indices for the vector of thermodynamical (=th) quantities stored in table */

  //@{

  int index_th_xe_;   /**< ionization fraction \f$ x_e \f$ */
  int index_th_rate_; /**< maximum variation rate of \f$ exp^{-\kappa}\f$, g and \f$ (d g / d \tau) \f$, used for computing integration step in perturbation module */
  int index_th_tau_d_;       /**< Baryon drag optical depth */
  int index_th_dkappa_;      /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_ddkappa_;     /**< scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_dddkappa_;    /**< scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_exp_m_kappa_; /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g_;   /**< visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_dg_;  /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_ddg_; /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */
  int index_th_dmu_idm_dr_; /**< scattering rate of idr with idm_dr (i.e. idr opacity to idm_dr scattering) (units 1/Mpc) */
  int index_th_ddmu_idm_dr_;  /**< derivative of this scattering rate */
  int index_th_dddmu_idm_dr_; /**< second derivative of this scattering rate */
  int index_th_dmu_idr_;      /**< idr self-interaction rate */
  int index_th_tau_idm_dr_;   /**< optical depth of idm_dr (due to interactions with idr) */
  int index_th_tau_idr_;      /**< optical depth of idr (due to self-interactions) */
  int index_th_g_idm_dr_;     /**< visibility function of idm_idr */
  int index_th_cidm_dr2_;     /**< interacting dark matter squared sound speed \f$ c_{dm}^2 \f$ */
  int index_th_Tidm_dr_;      /**< temperature of DM interacting with DR \f$ T_{idm_dr} \f$ */
  int index_th_Tb_;           /**< baryon temperature \f$ T_b \f$ */
  int index_th_wb_;           /**< baryon equation of state parameter \f$ w_b = k_B T_b / \mu \f$ */
  int index_th_cb2_;          /**< squared baryon adiabatic sound speed \f$ c_b^2 \f$ */
  int index_th_dcb2_; /**< derivative wrt conformal time of squared baryon sound speed \f$ d [c_b^2] / d \tau \f$ (only computed if some non-minimal tight-coupling schemes is requested) */
  int index_th_ddcb2_; /**< second derivative wrt conformal time of squared baryon sound speed  \f$ d^2 [c_b^2] / d \tau^2 \f$ (only computed if some non0-minimal tight-coupling schemes is requested) */
  int index_th_r_d_;   /**< simple analytic approximation to the photon comoving damping scale */

  //@}

 private:
  void thermodynamics_init();
  void thermodynamics_indices(recombination* preco, reionization* preio);
  void thermodynamics_helium_from_bbn();
  void thermodynamics_onthespot_energy_injection(recombination* preco,
                                                 double z,
                                                 double* energy_rate);
  void thermodynamics_energy_injection(recombination* preco, double z, double* energy_rate);
  void thermodynamics_reionization_function(double z, reionization* preio, double* xe);
  void thermodynamics_reionization(recombination* preco, reionization* preio, double* pvecback);
  void thermodynamics_reionization_sample(recombination* preco,
                                          reionization* preio,
                                          double* pvecback);
  void thermodynamics_get_xe_before_reionization(recombination* preco, double z, double* xe);
  void thermodynamics_recombination(recombination* preco, double* pvecback);
  void thermodynamics_recombination_integrate(recombination* prec, double* pvecback);
  double thermodynamics_recfast_hydrogen_saha_xH(const recombination* preco, double z) const;
  double thermodynamics_recfast_helium_first_saha_xe(const recombination* preco, double z) const;
  double thermodynamics_recfast_helium_second_saha_xe(const recombination* preco, double z) const;
  void thermodynamics_recfast_store_row(
      recombination* preco, int sample_index, double z, double xe, double Tb, double dTbdz) const;
  double thermodynamics_recfast_xe_after_helium_ode(const recombination* preco,
                                                    double z,
                                                    const double* y) const;
  double thermodynamics_recfast_xe_after_full_ode(const recombination* preco,
                                                  double z,
                                                  const double* y) const;
  void thermodynamics_recombination_derivs_member(double z,
                                                  double* y,
                                                  double* dy,
                                                  void* fixed_parameters);
  static void thermodynamics_recombination_derivs(double z,
                                                  double* y,
                                                  double* dy,
                                                  void* fixed_parameters);
  static void thermodynamics_recfast_derivs(double minus_z,
                                            double* y,
                                            double* dy,
                                            void* fixed_parameters);
  static void thermodynamics_recfast_output(
      double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace);
  static void thermodynamics_recfast_output_none(
      double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace);
  static void thermodynamics_recfast_timescale(double minus_z,
                                               void* parameters_and_workspace,
                                               double* timescale);
  void thermodynamics_merge_reco_and_reio(recombination* preco, reionization* preio);
  void thermodynamics_tanh(
      double x, double center, double before, double after, double width, double* result);

  BackgroundModulePtr background_module_;

  /** @name - thermodynamics interpolation tables */

  //@{
  std::vector<double>
      z_table_; /**< vector z_table[index_z] with values of redshift (vector of size tt_size) */
  std::vector<double>
      thermodynamics_table_; /**< table thermodynamics_table[index_z*tt_size_+pba->index_th] with all other quantities (array of size th_size*tt_size) */
  std::vector<double>
      d2thermodynamics_dz2_table_; /**< table d2thermodynamics_dz2_table[index_z*tt_size_+pba->index_th] with values of \f$ d^2 t_i / dz^2 \f$ (array of size th_size*tt_size) */
  //@}
};

enum class RecombinationPhase { analytic, helium, full };

/**
 * temporary  parameters and workspace passed to the thermodynamics_derivs function
 */

struct thermodynamics_parameters_and_workspace {
  thermodynamics_parameters_and_workspace(ThermodynamicsModule* p_m) : thermodynamics_module(p_m) {}
  ThermodynamicsModule* const thermodynamics_module;
  /* structures containing fixed input parameters (indices, ...) */
  recombination* preco;

  /* which atomic physics supplies dx_H/dz and dx_He/dz */
  const RecombinationModel* recombination_model = nullptr;

  /* workspace */
  double* pvecback;
  int recfast_output_index_offset        = 0;
  RecombinationPhase recombination_phase = RecombinationPhase::full;
};

#endif  //THERMODYNAMICS_MODULE_H
