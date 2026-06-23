/** @file perturbations.h Documented includes for perturbation module */

#ifndef __PERTURBATIONS__
#define __PERTURBATIONS__

#include <memory>
#include <vector>

#include "../species/base_species.h"
#include "evolver_ndf15.h"
#include "evolver_rkck.h"
#include "evolver_rkdp45.h"
#include "thermodynamics.h"

#define _scalars_ (ppt->has_scalars && (index_md == index_md_scalars_))
#define _vectors_ (ppt->has_vectors && (index_md == index_md_vectors_))
#define _tensors_ (ppt->has_tensors && (index_md == index_md_tensors_))

//TODO: Remove those when possible!
#define _scalarsEXT_ (ppt->has_scalars && (index_md == perturbations_module_->index_md_scalars_))
#define _vectorsEXT_ (ppt->has_vectors && (index_md == perturbations_module_->index_md_vectors_))
#define _tensorsEXT_ (ppt->has_tensors && (index_md == perturbations_module_->index_md_tensors_))

#define _set_source_(index) \
  sources_[index_md][index_ic * tp_size_[index_md] + index][index_tau * k_size_[index_md] + index_k]

/**
 * flags for various approximation schemes
 * (tca = tight-coupling approximation,
 *  rsa = radiation streaming approximation,
 *  ufa = massless neutrinos / ultra-relativistic relics fluid approximation)
 *
 * CAUTION: must be listed below in chronological order, and cannot be
 * reversible. When integrating equations for a given mode, it is only
 * possible to switch from left to right in the lists below.
 */

//@{

enum tca_flags { tca_on, tca_off };
enum rsa_flags { rsa_off, rsa_on };
enum tca_idm_dr_flags { tca_idm_dr_on, tca_idm_dr_off };
enum tca_idm_drmd_flags { tca_idm_drmd_on, tca_idm_drmd_off };
enum rsa_idr_flags { rsa_idr_off, rsa_idr_on };
enum ufa_flags { ufa_off, ufa_on };
enum ncdmfa_flags { ncdmfa_off, ncdmfa_on };

//@}

/**
 * labels for the way in which each approximation scheme is implemented
 */

//@{

enum idr_method { idr_free_streaming, idr_fluid }; /* for the idm-idr case */
enum tensor_methods { tm_photons_only, tm_massless_approximation, tm_exact };

//@}

/**
 * List of coded gauges. More gauges can in principle be defined.
 */

//@{

enum class possible_gauges {
  newtonian,  /**< newtonian (or longitudinal) gauge */
  synchronous /**< synchronous gauge with \f$ \theta_{cdm} = 0 \f$ by convention */
};

//@}

//@{

/**
 * maximum number and types of selection function (for bins of matter density or cosmic shear)
 */
#define _SELECTION_NUM_MAX_ 100
enum selection_type { gaussian, tophat, dirac };

//@}

//@{

/**
 * maximum number of k-values for perturbation output
 */
#define _MAX_NUMBER_OF_K_FILES_ 30

//@}

/**
 * Structure containing everything about perturbations that other
 * modules need to know, in particular tabled values of the source
 * functions \f$ S(k, \tau) \f$ for all requested modes
 * (scalar/vector/tensor), initial conditions, types (temperature,
 * E-polarization, B-polarization, lensing potential, etc), multipole
 * l and wavenumber k.
 *
 */

struct perturbs {
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  bool has_perturbations; /**< do we need to compute perturbations at all ? */

  bool
      has_cls; /**< do we need any harmonic space spectrum \f$ C_l \f$ (and hence Bessel functions, transfer functions, ...)? */

  bool has_scalars = true;  /**< do we need scalars? */
  bool has_vectors = false; /**< do we need vectors? */
  bool has_tensors = false; /**< do we need tensors? */

  bool has_ad  = true;  /**< do we need adiabatic mode? */
  bool has_bi  = false; /**< do we need isocurvature bi mode? */
  bool has_cdi = false; /**< do we need isocurvature cdi mode? */
  bool has_nid = false; /**< do we need isocurvature nid mode? */
  bool has_niv = false; /**< do we need isocurvature niv mode? */

  /* perturbed recombination */
  /** Do we want to consider perturbed temperature and ionization fraction? */
  bool has_perturbed_recombination = false;
  /** Neutrino contribution to tensors */
  enum tensor_methods tensor_method =
      tm_massless_approximation; /**< way to treat neutrinos in tensor perturbations(neglect, approximate as massless, take exact equations) */

  bool has_cl_cmb_temperature  = false; /**< do we need \f$ C_l \f$'s for CMB temperature? */
  bool has_cl_cmb_polarization = false; /**< do we need \f$ C_l \f$'s for CMB polarization? */
  bool has_cl_cmb_lensing_potential =
      false; /**< do we need \f$ C_l \f$'s for CMB lensing potential? */
  bool has_cl_lensing_potential =
      false;                        /**< do we need \f$ C_l \f$'s for galaxy lensing potential? */
  bool has_cl_number_count = false; /**< do we need \f$ C_l \f$'s for density number count? */
  bool has_pk_matter       = false; /**< do we need matter Fourier spectrum? */
  bool has_density_transfers =
      false; /**< do we need to output individual matter density transfer functions? */
  bool has_velocity_transfers =
      false; /**< do we need to output individual matter velocity transfer functions? */
  bool has_metricpotential_transfers =
      false; /**< do we need to output individual transfer functions for scalar metric perturbations? */
  bool has_Nbody_gauge_transfers =
      false; /**< should we convert density and velocity transfer functions to Nbody gauge? */

  bool has_nl_corrections_based_on_delta_m =
      false; /**< do we want to compute non-linear corrections with an algorithm relying on delta_m (like halofit)? */

  bool has_nc_density = false; /**< in dCl, do we want density terms ? */
  bool has_nc_rsd     = false; /**< in dCl, do we want redshift space distortion terms ? */
  bool has_nc_lens    = false; /**< in dCl, do we want lensing terms ? */
  bool has_nc_gr      = false; /**< in dCl, do we want gravity terms ? */

  int l_scalar_max = 2500; /**< maximum l value for CMB scalars \f$ C_l \f$'s */
  int l_vector_max = 500;  /**< maximum l value for CMB vectors \f$ C_l \f$'s */
  int l_tensor_max = 500;  /**< maximum l value for CMB tensors \f$ C_l \f$'s */
  int l_lss_max =
      300; /**< maximum l value for LSS \f$ C_l \f$'s (density and lensing potential in  bins) */
  double k_max_for_pk =
      1.; /**< maximum value of k in 1/Mpc in P(k) (if \f$ C_l \f$'s also requested, overseeded by value kmax inferred from l_scalar_max if it is bigger) */

  int selection_num                          = 1;        /**< number of selection functions
                                                   (i.e. bins) for matter density \f$ C_l \f$'s */
  enum selection_type selection              = gaussian; /**< type of selection functions */
  double selection_mean[_SELECTION_NUM_MAX_] = {
      1.}; /**< centers of selection functions; element [0] default 1 */
  double selection_width[_SELECTION_NUM_MAX_] = {
      0.1}; /**< widths of selection functions; element [0] default 0.1 */

  int switch_sw =
      1; /**< in temperature calculation, do we want to include the intrinsic temperature + Sachs Wolfe term? */
  int switch_eisw =
      1; /**< in temperature calculation, do we want to include the early integrated Sachs Wolfe term? */
  int switch_lisw =
      1; /**< in temperature calculation, do we want to include the late integrated Sachs Wolfe term? */
  int switch_dop = 1; /**< in temperature calculation, do we want to include the Doppler term? */
  int switch_pol =
      1; /**< in temperature calculation, do we want to include the polarization-related term? */
  double eisw_lisw_split_z =
      120; /**< at which redshift do we define the cut between eisw and lisw ?*/

  int store_perturbations = _FALSE_; /**< Do we want to store perturbations? */
  int k_output_values_num = 0;       /**< Number of perturbation outputs (default=0) */
  double k_output_values
      [_MAX_NUMBER_OF_K_FILES_]; /**< List of k values where perturbation output is requested. */

  double three_ceff2_ur =
      1.; /**< 3 x effective squared sound speed for the ultrarelativistic perturbations */
  double three_cvis2_ur =
      1.; /**< 3 x effective viscosity parameter for the ultrarelativistic perturbations */
  double G_eff_ur = 0.;

  double z_max_pk =
      0.; /**< when we compute only the matter spectrum / transfer functions, but not the CMB, we are sometimes interested to sample source functions at very high redshift, way before recombination. This z_max_pk will then fix the initial sampling time of the sources. */

  //@}

  /** @name - gauge in which to perform the calculation */

  //@{

  possible_gauges gauge =
      possible_gauges::synchronous; /**< gauge in which to perform this calculation */

  //@}

  /** @name - list of conformal time values in the source table
      (common to all modes and types) */

  //@{

  double
      selection_min_of_tau_min; /**< used in presence of selection functions (for matter density, cosmic shear...) */
  double
      selection_max_of_tau_max; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double
      selection_delta_tau; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double*
      selection_tau_min; /**< value of conformal time below which W(tau) is considered to vanish for each bin */
  double*
      selection_tau_max; /**< value of conformal time above which W(tau) is considered to vanish for each bin */
  double* selection_tau; /**< value of conformal time at the center of each bin */
  double*
      selection_function; /**< selection function W(tau), normalized to \f$ \int W(tau) dtau=1 \f$, stored in selection_function[bin*tau_size_+index_tau] */

  //@}

  /** @name - technical parameters */

  //@{

  short perturbations_verbose =
      0; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}
};

/**
 * Structure containing the indices and the values of the perturbation
 * variables which are integrated over time (as well as their
 * time-derivatives). For a given wavenumber, the size of these
 * vectors changes when the approximation scheme changes.
 */

struct perturb_vector {
  // Per-species perturbation layouts, parallel to all_species_ (lex-key order).
  // Each thread allocates its own pv, so this storage is per-thread.
  std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> species_layouts;

  // Non-owning view into species_layouts: the species that registered ≥1
  // perturbation variable in THIS pv's mode (no-ops like Lambda excluded), in
  // lex-key order. Consumed by the scalar/vector/tensor RHS dispatch loops.
  struct ActiveSpecies {
    const BaseSpecies* species;
    const BaseSpecies::PerturbLayout* layout;
    bool clusters_as_matter;  // ClustersAsMatter()    (delta_m/theta_m membership)
    bool is_cold;             // IsColdMatterSpecies() (delta_cb bucket)
  };
  std::vector<ActiveSpecies> active_species;

  // Non-owning view into species_layouts: the species that emit a per-species
  // transfer source (delta_g, delta_cdm, delta_ncdm, ...), in lex-key order.
  // Source-column membership is k-independent (fixed in
  // perturb_indices_of_perturbs); this list materialises it per pv so
  // perturb_sources_member dispatches FillSources by pointer. Empty in runs
  // without dTk/mTk/k_output, so that per-sample loop is skipped entirely.
  struct SourceSpecies {
    const BaseSpecies* species;
    const BaseSpecies::PerturbLayout* layout;
  };
  std::vector<SourceSpecies> source_species;

  // Always-present species, resolved once per pv (non-owning, into
  // species_layouts). Base-typed because perturbations.h cannot see the nested
  // PhotonsSpecies/BaryonsSpecies::PerturbLayout types (photons.h/baryons.h
  // include perturbations.h). Read sites static_cast to the concrete type.
  BaseSpecies::PerturbLayout* photon_layout = nullptr;
  BaseSpecies::PerturbLayout* baryon_layout = nullptr;
  BaseSpecies::PerturbLayout* cdm_layout    = nullptr;

  /* Photons bare fields removed — use pv->species_layouts[g_i] (PhotonsSpecies::PerturbLayout). */
  /* Baryons bare fields removed — use pv->species_layouts[b_i] (BaryonsSpecies::PerturbLayout). */
  /* CDM bare fields removed — use pv->species_layouts[cdm_i] (CDMSpecies::PerturbLayout). */
  /* IDM_DR bare fields removed — use pv->species_layouts[idm_dr_idr_i]
     (IDM_DR_IDR_Species::PerturbLayout::idm_dr). */
  /* IDM_DRMD bare fields removed — use pv->species_layouts[idm_drmd_idr_drmd_i]
     (IDM_DRMD_IDR_DRMD_Species::PerturbLayout::idm_drmd). */
  /* DCDM bare fields removed — use pv->species_layouts[dcdm_dr_i] (DCDM_DR_Species::PerturbLayout::dcdm). */
  /* UR bare fields removed — use pv->species_layouts[ur_i] (UltraRelativisticSpecies::PerturbLayout). */
  /* IDR bare fields removed — use pv->species_layouts[idm_dr_idr_i]
     (IDM_DR_IDR_Species::PerturbLayout::idr). */
  /* IDR_DRMD bare fields removed — use pv->species_layouts[idm_drmd_idr_drmd_i]
     (IDM_DRMD_IDR_DRMD_Species::PerturbLayout::idr_drmd). */

  /* perturbed recombination */
  int index_pt_perturbed_recombination_delta_temp; /**< Gas temperature perturbation */
  int index_pt_perturbed_recombination_delta_chi;  /**< Inionization fraction perturbation */

  int index_pt_eta;      /**< synchronous gauge metric perturbation eta*/
  int index_pt_phi;      /**< newtonian gauge metric perturbation phi */
  int index_pt_hv_prime; /**< vector metric perturbation h_v' in synchronous gauge */
  int index_pt_V;        /**< vector metric perturbation V in Newtonian gauge */

  int index_pt_gw;    /**< tensor metric perturbation h (gravitational waves) */
  int index_pt_gwdot; /**< its time-derivative */

  /** Tensor "relativistic neutrino" Boltzmann hierarchy. Unlike the scalar/vector
      hierarchies this is NOT owned by a species: it is a single hierarchy sourced
      by the gravitational waves (index_pt_gwdot) and exists whenever
      evolve_tensor_ur_ is set — i.e. when any massless relativistic species
      contributes, be it the ultra-relativistic neutrinos or massive ncdm under the
      massless approximation. The perturb_vector owns it and the perturbations
      module registers/evolves it directly (formerly the bare index_pt_*_ur fields). */
  struct TensorRelativisticLayout {
    int idx_delta = -1; /**< was index_pt_delta_ur */
    int idx_theta = -1; /**< was index_pt_theta_ur */
    int idx_shear = -1; /**< was index_pt_shear_ur */
    int idx_l3    = -1; /**< was index_pt_l3_ur (only registered when l_max >= 3) */
    int l_max     = -1; /**< was l_max_ur */
  } tensor_ur_layout;

  int pt_size; /**< size of perturbation vector */

  std::vector<double> y;  /**< vector of perturbations to be integrated */
  std::vector<double> dy; /**< time-derivative of the same vector */

  std::vector<int> used_in_sources; /**< boolean array specifying which
                           perturbations enter in the calculation of
                           source functions */
};

/**
 * Workspace containing, among other things, the value at a given time
 * of all background/perturbed quantities, as well as their indices.
 * There will be one such structure created for each mode
 * (scalar/.../tensor) and each thread (in case of parallel computing)
 */

struct perturb_workspace {
  /** @name - all possible useful indices for those metric
      perturbations which are not integrated over time, but just
      inferred from Einstein equations. "_mt_" stands for "metric".*/

  //@{

  int index_mt_psi;           /**< psi in longitudinal gauge */
  int index_mt_phi_prime;     /**< (d phi/d conf.time) in longitudinal gauge */
  int index_mt_h_prime;       /**< h' (wrt conf. time) in synchronous gauge */
  int index_mt_h_prime_prime; /**< h'' (wrt conf. time) in synchronous gauge */
  int index_mt_eta_prime;     /**< eta' (wrt conf. time) in synchronous gauge */
  int index_mt_alpha;         /**< \f$ \alpha = (h' + 6 \eta') / (2 k^2) \f$ in synchronous gauge */
  int index_mt_alpha_prime;   /**< \f$ \alpha'\f$ wrt conf. time) in synchronous gauge */
  int index_mt_gw_prime_prime; /**< second derivative wrt conformal time of gravitational wave field, often called h */
  int index_mt_V_prime; /**< derivative of Newtonian gauge vector metric perturbation V */
  int index_mt_hv_prime_prime; /**< Second derivative of Synchronous gauge vector metric perturbation \f$ h_v\f$ */
  int mt_size;                 /**< size of metric perturbation vector */

  //@}

  /** @name - value at a given time of all background/perturbed
      quantities
  */

  //@{

  std::vector<double> pvecback;   /**< background quantities */
  std::vector<double> pvecthermo; /**< thermodynamics quantities */
  std::vector<double> pvecmetric; /**< metric quantities */
  /* Single owner of the perturb_vector; species methods receive a
     non-owning raw pointer via pv.get() (used transiently while the
     workspace is alive, never stored or shared). */
  std::unique_ptr<perturb_vector> pv; /**< owning pointer to vector of
                                 integrated perturbations and their
                                 time-derivatives */

  double delta_rho;        /**< total density perturbation (gives delta Too) */
  double rho_plus_p_theta; /**< total (rho+p)*theta perturbation (gives delta Toi) */
  double rho_plus_p_shear; /**< total (rho+p)*shear (gives delta Tij) */
  double delta_p;          /**< total pressure perturbation (gives Tii) */

  double rho_plus_p_tot; /**< total (rho+p) (used to infer theta_tot from rho_plus_p_theta) */

  double
      gw_source; /**< stress-energy source term in Einstein's tensor equations (gives Tij[tensor]) */
  double vector_source_pi; /**< first stress-energy source term in Einstein's vector equations */
  double vector_source_v;  /**< second stress-energy source term in Einstein's vector equations */

  /* Module-owned approximation closures. Computed by the perturbations module
     (perturb_tca_slip_and_shear and perturb_rsa_* / perturb_rsa_idr_*, the
     latter run inside perturb_einstein after the metric solve because they
     depend on h'/phi') and then read — read-only — by the relevant species
     (photons, ultra-relativistic, IDM_DR). These are not per-species private
     state despite the species-suffixed names: they live here because the module
     owns the approximation machinery and several consumers cross species. */
  double tca_shear_g;  /**< photon shear in tight-coupling approximation */
  double tca_slip;     /**< photon-baryon slip in tight-coupling approximation */
  double rsa_delta_g;  /**< photon density in radiation streaming approximation */
  double rsa_theta_g;  /**< photon velocity in radiation streaming approximation */
  double rsa_delta_ur; /**< photon density in radiation streaming approximation */
  double rsa_theta_ur; /**< photon velocity in radiation streaming approximation */
  double
      rsa_delta_idr; /**< interacting dark radiation density in dark radiation streaming approximation */
  double
      rsa_theta_idr; /**< interacting dark radiation velocity in dark radiation streaming approximation */

  double delta_m; /**< relative density perturbation of all non-relativistic species */
  double theta_m; /**< velocity divergence theta of all non-relativistic species */

  double delta_cb; /**< relative density perturbation of only cdm and baryon */
  double theta_cb; /**< velocity divergence theta of only cdm and baryon */

  /* Module-owned PPF dark-energy closure. Written by FluidSpecies::ComputePpf
     during stress-energy assembly (it depends on the totals of every other
     species), then read back by the module totals and by the fluid's own
     derivs/sources. PPF fluid is defined relative to the global stress-energy,
     so this is genuine module-level coupling state, not fluid-private scratch.
     (The S quantity of eq. 15 in 0808.3125 is a ComputePpf-internal temporary,
     no longer stored here.) */
  double delta_rho_fld; /**< density perturbation of fluid, not so trivial in PPF scheme */
  double delta_p_fld;   /**< pressure perturbation of fluid, very non-trivial in PPF scheme */
  double rho_plus_p_theta_fld; /**< velocity divergence of fluid, not so trivial in PPF scheme */
  double Gamma_prime_fld;      /**< Gamma_prime in PPF scheme (equivalent to eq. 14 in 0808.3125) */

  FILE* perturb_output_file; /**< filepointer to output file*/
  int index_ikout;           /**< index for output k value (when k_output_values is set) */

  //@}

  /** @name - pre-computed scalar perturbation context, populated by
   *  perturb_derivs_member() before calling species::PerturbDerivs(). */
  //@{
  PerturbScalarContext scalar_ctx;
  //@}

  /** @name - indices useful for searching background/thermo quantities in tables */

  //@{

  short
      inter_mode; /**< flag defining the method used for interpolation background/thermo quantities tables */

  int last_index_back; /**< the background interpolation function background_at_tau() keeps memory of the last point called through this index */
  int last_index_thermo; /**< the thermodynamics interpolation function thermodynamics_at_z() keeps memory of the last point called through this index */

  //@}

  /** @name - approximations used at a given time */

  //@{

  int index_ap_tca;          /**< index for tight-coupling approximation */
  int index_ap_rsa;          /**< index for radiation streaming approximation */
  int index_ap_tca_idm_dr;   /**< index for dark tight-coupling approximation (idm-idr) */
  int index_ap_tca_idm_drmd; /**< index for dark tight-coupling approximation (DRMD) */
  int index_ap_rsa_idr;      /**< index for dark radiation streaming approximation */
  int index_ap_ufa;          /**< index for ur fluid approximation */
  int index_ap_ncdmfa;       /**< index for ncdm fluid approximation */
  int ap_size;               /**< number of relevant approximations for a given mode */

  std::vector<int>
      approx; /**< array of approximation flags holding at a given time: approx[index_ap] */

  //@}

  /** @name - approximations used at a given time */

  //@{

  int max_l_max; /**< maximum l_max for any multipole */
  std::vector<double>
      s_l; /**< array of freestreaming coefficients \f$ s_l = \sqrt{1-K*(l^2-1)/k^2} \f$*/
  double
      cotKgen; /**< generalised cot(sqrt(|K|)*tau)/k, for closing free-streaming hierarchies; valid for all modes */

  //@}
};

#endif
