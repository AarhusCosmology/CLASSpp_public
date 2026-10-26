#ifndef TRANSFER_MODULE_H
#define TRANSFER_MODULE_H

#include <vector>

#include "asymptotic_extrapolation.h"
#include "base_module.h"
#include "bessel_table_cache.h"
#include "input_module.h"

class TransferModule : public BaseModule {
 public:
  TransferModule(InputModulePtr input_module,
                 BackgroundModulePtr background_module,
                 ThermodynamicsModulePtr thermodynamics_module,
                 PerturbationsModulePtr perturbations_module,
                 NonlinearModulePtr nonlinear_module,
                 std::shared_ptr<BesselTableCache> bessel_cache);
  ~TransferModule();

  /** @name - number of modes and transfer function types */
  //@{

  /* class_define_index assigns only the requested types; the others stay -1, and are still
     read (transfer_type_properties passes every one of them by value) */
  int index_tt_t0_      = -1; /**< index for transfer type = temperature (j=0 term) */
  int index_tt_t1_      = -1; /**< index for transfer type = temperature (j=1 term) */
  int index_tt_t2_      = -1; /**< index for transfer type = temperature (j=2 term) */
  int index_tt_e_       = -1; /**< index for transfer type = E-polarization */
  int index_tt_b_       = -1; /**< index for transfer type = B-polarization */
  int index_tt_lcmb_    = -1; /**< index for transfer type = CMB lensing */
  int index_tt_density_ = -1; /**< index for first bin of transfer type = matter density */
  int index_tt_lensing_ = -1; /**< index for first bin of transfer type = galaxy lensing */

  int index_tt_rsd_ =
      -1; /**< index for first bin of transfer type = redshift space distortion of number count */
  int index_tt_d0_ =
      -1; /**< index for first bin of transfer type = doppler effect for of number count (j=0 term) */
  int index_tt_d1_ =
      -1; /**< index for first bin of transfer type = doppler effect for of number count (j=1 term) */
  int index_tt_nc_lens_ =
      -1; /**< index for first bin of transfer type = lensing for of number count */
  int index_tt_nc_g1_ =
      -1; /**< index for first bin of transfer type = gravity term G1 for of number count */
  int index_tt_nc_g2_ =
      -1; /**< index for first bin of transfer type = gravity term G2 for of number count */
  int index_tt_nc_g3_ =
      -1; /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g4_ =
      -1; /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g5_ =
      -1; /**< index for first bin of transfer type = gravity term G3 for of number count */

  std::vector<int>
      tt_size_; /**< number of requested transfer types tt_size[index_md] for each mode */
  //@}

  /** @name - number and list of multipoles */
  //@{
  int l_size_max_ = 0; /**< greatest of all l_size[index_md] */
  std::vector<std::vector<int>>
      l_size_tt_; /**< number of multipole values for which we effectively compute the transfer function,l_size_tt[index_md][index_tt] */
  std::vector<int>
      l_size_;         /**< number of multipole values for each requested mode, l_size[index_md] */
  std::vector<int> l_; /**< list of multipole values l[index_l] */
  //@}

  /** @name - number and list of wavenumbers */
  //@{
  int q_size_ = 0;        /**< number of wavenumber values */
  std::vector<double> q_; /**< list of wavenumber values, q[index_q] */
  std::vector<std::vector<double>>
      k_; /**< list of wavenumber values for each requested mode, k[index_md][index_q]. In flat universes k=q. In non-flat universes q and k differ through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector, tensor. q should be used throughout the transfer module, excepted when interpolating or manipulating the source functions S(k,tau): for a given value of q this should be done in k(q). */
  int index_q_flat_approximation_ =
      -1; /**< index of the first q value using the flat rescaling approximation */
  //@}

  /** @name - full Limber scheme for CMB lensing */
  //@{
  bool do_lcmb_full_limber_ = false;
  /** The linear CMB lensing source past the last perturbation mode, per initial condition
   *  and time: lensing_source_extrapolation_[index_ic * tau_size + index_tau]. Empty unless
   *  the full-Limber grid reaches beyond that mode. */
  std::vector<AsymptoticExtrapolation> lensing_source_extrapolation_;
  int q_size_limber_ = 0;
  std::vector<double> q_limber_;
  std::vector<std::vector<double>> k_limber_;
  std::vector<std::vector<double>> transfer_limber_;
  //@}

  /** @name - transfer functions */
  //@{
  std::vector<std::vector<double>>
      transfer_; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * transfer_module_->tt_size_[index_md] + index_tt) * transfer_module_->l_size_[index_md] + index_l) * transfer_module_->q_size_ + index_q] */
                 //@}

  /** @name - the Limber delta */
  //@{
  /**
   * Where the Limber delta sits, and the curvature factor that belongs there.
   *
   * The prescription is j_l(x) -> Int_0^inf j_l dx * delta(x - L) at the
   * Langer-corrected multipole L = l + 1/2. Substituting U = Phi sin_K(chi)
   * turns the hyperspherical radial equation into the flat one with
   * chi -> sin_K(chi), so the delta sits at sin_K(chi_L) = L/q and its weight
   * picks up the WKB factor cos_K(chi_L)^(-1/2) = (1 - K L^2/q^2)^(-1/4),
   * evaluated at that same L. Issue #423; see
   * docs/superpowers/specs/2026-09-08-limber-effective-multipole-design.md.
   */
  struct LimberPoint {
    double tau0_minus_tau; /**< chi_L: the distance at which the source is read */
    double amplitude;      /**< cos_K(chi_L)^(-1/2); exactly 1 in flat space */
  };
  static LimberPoint LimberDeltaPoint(double l, double q, int sgnK, double K);

  /**
   * True when the source support reaches past the equator of a closed universe.
   *
   * `LimberDeltaPoint` places the delta at the near turning point, the only one
   * `asin` can return. A closed universe has a second one at `pi R - chi_L`, so
   * the approximation is complete only while the sources stay inside the first
   * quarter circumference. Flat and open geometries have no equator.
   */
  static bool LimberSourceCrossesEquator(double tau0_minus_tau_max, int sgnK, double K);

  /**
   * Lower anchor of the separate full-Limber CMB lensing quadrature grid.
   *
   * That grid is a geometric progression, so its anchor fixes every node. A
   * closed universe has modes down to `nu = 3` (`q = 3 sqrt(K)`), and anchoring
   * there made the whole grid slide with `K`, so it never converged to the flat
   * grid as `K -> 0` and `C_l^phiphi` inherited the quadrature-phase difference.
   * Below the flat truncation scale those nodes carry no integrand anyway -- the
   * Limber transfer vanishes for `q < l_switch_limber/tau0` -- so the closed
   * anchor is clamped by the flat one. Flat and open are unchanged: open's
   * `sqrt(k_min^2 + K)` already collapses to the flat anchor.
   */
  static double LimberGridQMin(double k_min, double k_min_flat, int sgnK, double K);
  //@}

 private:
  void transfer_functions_at_q(
      int index_md, int index_ic, int index_tt, int index_l, double q, double* transfer_function);
  void transfer_init();
  void transfer_indices_of_transfers(double q_period, double K, int sgnK);
  void transfer_type_properties();
  void transfer_perturbation_copy_sources_and_nl_corrections(
      std::vector<std::vector<double*>>& sources,
      std::vector<std::vector<std::vector<double>>>& owned_sources);
  void transfer_lensing_source_extrapolation();
  void transfer_extrapolate_lensing_source(double k,
                                           int index_ic,
                                           double* interpolated_sources) const;
  void transfer_perturbation_source_spline(
      const std::vector<std::vector<double*>>& sources,
      std::vector<std::vector<double*>>& sources_spline_ptrs,
      std::vector<std::vector<std::vector<double>>>& sources_spline_storage);
  void transfer_get_l_list();
  void transfer_get_q_list(double q_period, double K, int sgnK);
  void transfer_get_q_limber_list(double K, int sgnK);
  void transfer_get_k_list(double K);
  void transfer_selection_times();
  void transfer_time_samplings(double tau_rec);
  std::vector<double> transfer_logarithmic_sampling(double tau_min,
                                                    double tau_max,
                                                    double log_step);
  void transfer_precompute_windows(double tau_rec);
  void transfer_compute_for_each_q(int index_q,
                                   double tau_rec,
                                   double** const* pert_sources,
                                   double** const* pert_sources_spline,
                                   struct transfer_workspace* ptw,
                                   bool use_full_limber);
  void transfer_radial_coordinates(struct transfer_workspace* ptw,
                                   const double* tau0_minus_tau,
                                   int tau_size,
                                   double k);
  void transfer_interpolate_sources(double k,
                                    int index_md,
                                    const double* pert_source,
                                    const double* pert_source_spline,
                                    double* interpolated_sources);
  void transfer_sources(
      const double* interpolated_sources, double k, int index_md, int index_tt, double* sources);
  void transfer_selection_function(int bin, double z, double* selection);
  void transfer_dNdz_analytic(double z, double* dNdz, double* dln_dNdz_dz);
  void transfer_selection_compute(double* selection,
                                  const double* tau0_minus_tau,
                                  const double* w_trapz,
                                  int tau_size,
                                  double* pvecback,
                                  double tau0,
                                  int bin);
  TransferMethod transfer_select_method(struct transfer_workspace* ptw,
                                        int index_q,
                                        int index_md,
                                        int index_tt,
                                        int index_l,
                                        double q,
                                        double ra_rec,
                                        double q_max_bessel,
                                        const TransferTimeSampling& sampling,
                                        const double* sources,
                                        bool use_full_limber,
                                        double tau0_minus_tau_min_bessel,
                                        int* index_tau_max,
                                        char* bessel_truncation);
  void transfer_integrate(struct transfer_workspace* ptw,
                          int index_l,
                          double k,
                          radial_function_type radial_type,
                          const TransferTimeSampling& sampling,
                          const double* sources,
                          int index_tau_max,
                          bool bessel_truncation,
                          double tau0_minus_tau_min_bessel,
                          double* trsf);
  void transfer_limber(struct transfer_workspace* ptw,
                       double l,
                       double q,
                       radial_function_type radial_type,
                       const std::vector<double>& tau0_minus_tau,
                       const double* sources,
                       double* trsf);
  void transfer_limber_interpolate(const double* tau0_minus_tau,
                                   const double* sources,
                                   int tau_size,
                                   double tau0_minus_tau_limber,
                                   double* S);
  /** which of the Bessel function and its first two derivatives a radial function needs */
  struct BesselNeeds {
    bool Phi;
    bool dPhi;
    bool d2Phi;
  };
  static BesselNeeds transfer_radial_function_needs(radial_function_type radial_type);
  double transfer_bessel_cut(struct transfer_workspace* ptw, int index_q, int index_l, double k);
  void transfer_bessel_interpolate(struct transfer_workspace* ptw,
                                   int index_q,
                                   int index_l,
                                   int tau_size,
                                   int x_size,
                                   BesselNeeds needs);
  void transfer_radial_function(struct transfer_workspace* ptw,
                                double k,
                                int index_l,
                                int x_size,
                                double* radial_function,
                                radial_function_type radial_type);
  void transfer_global_selection_read();
  void transfer_workspace_init(struct transfer_workspace* ptw,
                               double K,
                               int sgnK,
                               double tau0_minus_tau_cut,
                               const HyperInterpStruct* pBIS);
  void transfer_update_HIS(struct transfer_workspace* ptw, int index_q, double tau0);
  void transfer_get_lmax(void (*get_xmin_generic)(int sgnK,
                                                  int l,
                                                  double nu,
                                                  double xtol,
                                                  double phiminabs,
                                                  double* x_nonzero,
                                                  int* fevals),
                         int sgnK,
                         double nu,
                         int* lvec,
                         int lsize,
                         double phiminabs,
                         double xmax,
                         double xtol,
                         int* index_l_left,
                         int* index_l_right);
  void transfer_f_evo(const double* pvecback, int last_index, double cotKgen, double* f_evo);

  BackgroundModulePtr background_module_;
  ThermodynamicsModulePtr thermodynamics_module_;
  PerturbationsModulePtr perturbations_module_;
  NonlinearModulePtr nonlinear_module_;
  /** where the flat Bessel table comes from, possibly kept from an earlier cosmology */
  std::shared_ptr<BesselTableCache> bessel_cache_;
  /** largest x = q (tau0 - tau) at which this cosmology needs the flat Bessel functions. The
      table may reach further; decisions use this, so they do not depend on what was cached. */
  double bessel_x_max_ = 0.;

  bool has_cls_ = false; /**< copy of same flag in perturbation structure */
  int md_size_  = 0;     /**< number of modes included in computation */

  /** @name - transfer types and time samplings, defined once before the loop over wavenumbers
   *  (see docs/superpowers/specs/2026-10-08-transfer-restructure-design.md) */
  //@{
  std::vector<std::vector<TransferType>>
      types_; /**< properties of each type, types_[index_md][index_tt] */
  std::vector<std::vector<TransferTimeSampling>>
      samplings_; /**< distinct time samplings of each mode, samplings_[index_md][index_ts] */
  int tau_size_max_  = 0; /**< largest number of sampled times over all samplings */
  int tt_per_ts_max_ = 0; /**< largest number of types sharing a time sampling */
  std::vector<double>
      selection_tau_min_; /**< smallest time at which the selection function of each bin is non-zero */
  std::vector<double> selection_tau_mean_; /**< time at the mean redshift of each bin */
  std::vector<double>
      selection_tau_max_; /**< largest time at which the selection function of each bin is non-zero */
  //@}

  int nz_size_ = 0;          /**< number of redshift values in input tabulated selection function */
  std::vector<double> nz_z_; /**< redshift values in input tabulated selection function */
  std::vector<double> nz_nz_;   /**< input tabulated values of selection function */
  std::vector<double> nz_ddnz_; /**< second derivatives in splined selection function*/

  int nz_evo_size_ = 0; /**< number of redshift values in input tabulated evolution function */
  std::vector<double> nz_evo_z_;       /**< redshift values in input tabulated evolution function */
  std::vector<double> nz_evo_nz_;      /**< input tabulated values of evolution function */
  std::vector<double> nz_evo_dlog_nz_; /**< log of tabulated values of evolution function */
  std::vector<double>
      nz_evo_dd_dlog_nz_; /**< second derivatives in splined log of evolution function */

  std::vector<std::vector<double>>
      selection_file_ddNdz_; /**< for selection=tabulated: second derivatives of the selection function of each bin with respect to z, selection_file_ddNdz_[bin][index_z] */
};

#endif  //TRANSFER_MODULE_H
