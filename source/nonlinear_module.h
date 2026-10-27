#ifndef NONLINEAR_MODULE_H
#define NONLINEAR_MODULE_H

#include <vector>

#include "base_module.h"
#include "input_module.h"

class NonlinearModule : public BaseModule {
 public:
  NonlinearModule(InputModulePtr input_module,
                  BackgroundModulePtr background_module,
                  PerturbationsModulePtr perturbations_module,
                  PrimordialModulePtr primordial_module);
  ~NonlinearModule();

  /* external functions (meant to be called from other modules) */
  void nonlinear_pk_at_z(enum linear_or_logarithmic mode,
                         enum pk_outputs pk_output,
                         double z,
                         int index_pk,
                         double* out_pk,
                         double* out_pk_ic) const;
  void nonlinear_pks_at_z(enum linear_or_logarithmic mode,
                          enum pk_outputs pk_output,
                          double z,
                          double* out_pk,
                          double* out_pk_ic,
                          double* out_pk_cb,
                          double* out_pk_cb_ic) const;
  void nonlinear_pk_at_k_and_z(enum pk_outputs pk_output,
                               double k,
                               double z,
                               int index_pk,
                               double* out_pk,
                               double* out_pk_ic) const;
  void nonlinear_pks_at_k_and_z(enum pk_outputs pk_output,
                                double k,
                                double z,
                                double* out_pk,
                                double* out_pk_ic,
                                double* out_pk_cb,
                                double* out_pk_cb_ic) const;
  void nonlinear_pks_at_kvec_and_zvec(enum pk_outputs pk_output,
                                      double* kvec,
                                      int kvec_size,
                                      double* zvec,
                                      int zvec_size,
                                      double* out_pk,
                                      double* out_pk_cb) const;
  /* Hand back the stored P(k,z) table on CLASS's own k sampling, so that a
     caller (the Python wrapper) never has to know this module's array layout.
     These are the tabulated nodes, not an interpolation.
     Design: docs/superpowers/specs/2026-09-16-classy-lss-methods-design.md */
  void GetPkGrid(enum pk_outputs pk_output,
                 int index_pk,
                 std::vector<double>& k,
                 std::vector<double>& pk) const;
  double GetSigma8(int index_pk) const;
  void nonlinear_sigmas_at_z(
      double R, double z, int index_pk, enum out_sigmas sigma_output, double* result) const;
  void nonlinear_pk_tilt_at_k_and_z(
      enum pk_outputs pk_output, double k, double z, int index_pk, double* pk_tilt) const;
  void nonlinear_k_nl_at_z(double z, double* k_nl, double* k_nl_cb) const;

  /* The default continuation of a source beyond the computed k (extrap_max_scaled):
     source ~ ln(a k), with a fitted through the last two points of the table. */
  static double ExtrapolateMaxScaled(
      double k, double k_max, double k_previous, double source_max, double source_previous);

  // Deprecated:
  void nonlinear_sigma_at_z(
      double R, double z, int index_pk, double k_per_decade, double* result) const;

  int k_size_ = 0;           /**< k_size = total number of k values */
  std::vector<double> ln_k_; /**< ln_k[index_k] = list of log(k) values */

  std::vector<std::vector<double>>
      nl_corr_density_; /**< nl_corr_density[index_pk][index_tau * perturbations_module_->k_size_ + index_k] */

  /** sqrt(P_NL/P_L) at a k beyond the last perturbation mode, up to
   *  perturbs::k_max_lensing_extrapolated (held constant beyond), for the CMB lensing source
   *  continued there; 1 where the correction is not computed. */
  [[nodiscard]] double NlCorrDensityBeyondKMax(int index_pk, int index_tau, double k) const;

  std::vector<bool>
      is_non_zero_; /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */
  int ic_size_ =
      0; /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int ic_ic_size_ =
      0; /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  bool has_pk_m_  = false; /**< do we want spectra for total matter? */
  bool has_pk_cb_ = false; /**< do we want spectra for cdm+baryons? */

  int index_pk_m_  = -1; /**< index of pk for matter (meaningful only when has_pk_m_) */
  int index_pk_cb_ = -1; /**< index of pk for cdm+baryons (meaningful only when has_pk_cb_) */

  /* and two redundent but useful indices: */

  int index_pk_total_   = -1; /**< always equal to index_pk_m
                                  (always defined, useful e.g. for weak lensing spectrum) */
  int index_pk_cluster_ = -1; /**< equal to index_pk_cb if it exists, otherwise to index_pk_m
                                  (always defined, useful e.g. for galaxy clustering spectrum).
                                  delta_cb is the cold part of delta_m, so without warm matter
                                  P_cb is P_m and is stored only once */

  int pk_size_ = 0;            /**< k_size = total number of pk */
  std::vector<double> sigma8_; /**< sigma8[index_pk] */

 private:
  /* internal functions */
  void nonlinear_init();
  void nonlinear_indices();
  void nonlinear_get_k_list();
  void nonlinear_get_tau_list();
  void nonlinear_get_source(int index_k,
                            int index_ic,
                            int index_tp,
                            int index_tau,
                            const std::vector<std::vector<double>>& sources,
                            double* source);
  void nonlinear_pk_linear(int index_pk, int index_tau, int k_size, double* lnpk, double* lnpk_ic);
  void nonlinear_sigmas(double R,
                        double* lnpk_l,
                        double* ddlnpk_l,
                        int k_size,
                        double k_per_decade,
                        enum out_sigmas sigma_output,
                        double* result) const;
  void nonlinear_halofit(int index_pk,
                         double tau,
                         double* pk_nl,
                         double* lnpk_l,
                         double* ddlnpk_l,
                         double* k_nl,
                         bool* nl_corr_not_computable_at_this_k);
  void nonlinear_halofit_integrate(const double* integrand_array,
                                   int integrand_size,
                                   int ia_size,
                                   int index_ia_k,
                                   int index_ia_pk,
                                   const double* integration_weights,
                                   double R,
                                   enum halofit_integral_type type,
                                   double* sum);
  void nonlinear_hmcode(int index_pk,
                        double tau,
                        double* pk_nl,
                        std::vector<std::vector<double>>& lnpk_l,
                        std::vector<std::vector<double>>& ddlnpk_l,
                        double* k_nl,
                        bool* nl_corr_not_computable_at_this_k,
                        const nonlinear_workspace* pnw,
                        nonlinear_time_workspace* pntw);
  void nonlinear_hmcode_compute(HMcodeVersion version,
                                int index_pk,
                                double tau,
                                double* pk_nl,
                                std::vector<std::vector<double>>& lnpk_l,
                                std::vector<std::vector<double>>& ddlnpk_l,
                                double* k_nl,
                                bool* nl_corr_not_computable_at_this_k,
                                const nonlinear_workspace* pnw,
                                const nonlinear_time_workspace* pntw);
  void nonlinear_hmcode_noradiation_growth(nonlinear_workspace* pnw);
  double nonlinear_hmcode_growth_norad(const nonlinear_workspace* pnw, double a) const;
  void nonlinear_hmcode_eisenstein_hu(nonlinear_workspace* pnw);
  void nonlinear_hmcode_wiggle_init(nonlinear_workspace* pnw);
  void nonlinear_hmcode_nowiggle(const double* lnpk_l,
                                 const double* ddlnpk_l,
                                 const nonlinear_workspace* pnw,
                                 nonlinear_time_workspace* pntw);
  void nonlinear_hmcode_workspace_init(nonlinear_workspace* pnw);
  void nonlinear_hmcode_dark_energy_correction(nonlinear_workspace* pnw);
  void nonlinear_hmcode_baryonic_feedback();
  void nonlinear_hmcode_fill_sigtab(double* lnpk_l,
                                    double* ddlnpk_l,
                                    nonlinear_time_workspace* pntw);
  void nonlinear_hmcode_fill_growtab(nonlinear_workspace* pnw);
  void nonlinear_hmcode_growint(
      double a, double w, double wa, double Omega_k, double Omega_de, double* growth);
  void nonlinear_hmcode_window_nfw(double k, double rv, double c, double* window_nfw);
  void nonlinear_hmcode_halomassfunction(double nu, double* hmf);
  void prepare_pk_eq();

  BackgroundModulePtr background_module_;
  PerturbationsModulePtr perturbations_module_;
  PrimordialModulePtr primordial_module_;

  /** @name - arrays for the Fourier power spectra P(k,tau) */

  //@{

  bool has_pk_matter_ = false; /**< do we need matter Fourier spectrum? */

  std::vector<double> k_; /**< k[index_k] = list of k values */

  std::vector<double> ln_tau_; /**< log(tau) array, only needed if user wants
                        some output at z>0, instead of only z=0.  This
                        array only covers late times, used for the
                        output of P(k) or T(k), and matching the
                        condition z(tau) < z_max_pk */

  int ln_tau_size_ = 0; /**< number of values in this array */

  std::vector<std::vector<double>> ln_pk_ic_l_; /**< Matter power spectrum (linear).
                             Depends on indices index_pk, index_ic1_ic2, index_k, index_tau as:
                             ln_pk_ic_l[index_pk][(index_tau * k_size_ + index_k)* ic_ic_size_ + index_ic1_ic2]
                             where index-pk labels P(k) types (_m = total matter, _cb = baryons+CDM),
                             while index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since
                             the primordial spectrum is symmetric in (index_ic1, index_ic2)).
                             - for diagonal elements (index_ic1 = index_ic2) this arrays contains
                             ln[P(k)] where P(k) is positive by construction.
                             - for non-diagonal elements this arrays contains the k-dependent
                             cosine of the correlation angle, namely
                             P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                             This choice is convenient since the sign of the non-diagonal cross-correlation
                             could be negative. For fully correlated or anti-correlated initial conditions,
                             this non-diagonal element is independent on k, and equal to +1 or -1.
                          */

  std::vector<std::vector<double>>
      ddln_pk_ic_l_; /**< second derivative of above array with respect to log(tau), for spline interpolation. So:
                             - for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is
                             good since this function is usually smooth.
                             - for non-diagonal coefficients, we spline
                             P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                             vs. ln(k), which is fine since this quantity is often assumed to be
                             constant (e.g for fully correlated/anticorrelated initial conditions)
                             or nearly constant, and with arbitrary sign.
                          */

  std::vector<std::vector<double>>
      ln_pk_l_; /**< Total matter power spectrum summed over initial conditions (linear).
                          Only depends on indices index_pk,index_k, index_tau as:
                          ln_pk[index_pk][index_tau * k_size_ + index_k]
                       */

  std::vector<std::vector<double>>
      ddln_pk_l_; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  std::vector<std::vector<double>>
      ln_pk_nl_; /**< Total matter power spectrum summed over initial conditions (nonlinear).
                          Only depends on indices index_pk,index_k, index_tau as:
                          ln_pk[index_pk][index_tau * k_size_ + index_k]
                       */

  std::vector<std::vector<double>>
      ddln_pk_nl_; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  int k_size_extra_ = 0; /** total number of k values of extrapolated k array (high k)*/
  int k_size_nl_ =
      0; /**< P_NL and the correction are computed for index_k < k_size_nl_: k_size_, and beyond for the CMB lensing source continued past the last mode */
  double k_max_trusted_ =
      0.; /**< largest k of the linear spectrum that halofit and HMcode trust: the last mode, or k_max_extrapolated if larger */
  std::vector<std::vector<double>>
      nl_corr_density_extra_; /**< nl_corr_density_extra_[index_pk][index_tau * (k_size_nl_ - k_size_) + index_k - k_size_] */

  int tau_size_ = 0;        /**< tau_size = number of values */
  std::vector<double> tau_; /**< tau[index_tau] = list of time values, covering
                      all the values of the perturbation module */

  std::vector<std::vector<double>>
      k_nl_; /**< wavenumber at which non-linear corrections become important,
                                    defined differently by different non_linear_method's */
  int index_tau_min_nl_ =
      -1; /**< index of smallest value of tau at which nonlinear corrections have been computed
                                    (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - information on number of modes and pairs of initial conditions */

  //@{

  int index_md_scalars_ = -1; /**< set equal to psp->index_md_scalars
                           (useful since this module only deals with
                           scalars) */

  //@}

  /** @name - information on the type of power spectra (_cb, _m...) */

  //@{

  double c_min_ =
      0.0; /** for HMcode: minimum concentration in Bullock 2001 mass-concentration relation */
  double eta_0_ = 0.0; /** for HMcode: halo bloating parameter */

  //@}

  /** @name - parameters for the pk_eq method */

  //@{

  int index_pk_eq_w_       = -1;  /**< index of w in table pk_eq_w_and_Omega */
  int index_pk_eq_Omega_m_ = -1;  /**< index of Omega_m in table pk_eq_w_and_Omega */
  int pk_eq_size_          = 0;   /**< number of indices in table pk_eq_w_and_Omega */
  int pk_eq_tau_size_      = 0;   /**< number of times (and raws in table pk_eq_w_and_Omega) */
  std::vector<double> pk_eq_tau_; /**< table of time values */
  std::vector<double> pk_eq_w_and_Omega_;     /**< table of background quantites */
  std::vector<double> pk_eq_ddw_and_ddOmega_; /**< table of second derivatives */

  //@}
};

#endif  //NONLINEAR_MODULE_H
