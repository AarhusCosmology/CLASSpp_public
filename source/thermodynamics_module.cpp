/** @file thermodynamics_module.cpp
 * Computes recombination and reionization histories and exposes interpolated
 * thermodynamic quantities for the lifetime of a ThermodynamicsModule instance.
 */

#include "thermodynamics_module.h"

#include <algorithm>
#include <fstream>

#include "../species/fluid.h"
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
#include "background_module.h"
#include "evolver_erk.h"
#include "evolver_ndf15.h"
#include "evolver_rkdp45.h"
#include "evolver_tsit5.h"
#include "hyrec_model.h"
#include "recfast_model.h"

ThermodynamicsModule::ThermodynamicsModule(InputModulePtr input_module,
                                           BackgroundModulePtr background_module)
    : BaseModule(std::move(input_module)), background_module_(std::move(background_module)) {
  thermodynamics_init();
}

ThermodynamicsModule::~ThermodynamicsModule() {}

// Wrapper functions to pass non-static member functions
void ThermodynamicsModule::thermodynamics_recombination_derivs(double z,
                                                               double* y,
                                                               double* dy,
                                                               void* fixed_parameters) {
  auto tppaw = static_cast<thermodynamics_parameters_and_workspace*>(fixed_parameters);
  tppaw->thermodynamics_module->thermodynamics_recombination_derivs_member(z,
                                                                           y,
                                                                           dy,
                                                                           fixed_parameters);
}

void ThermodynamicsModule::thermodynamics_recfast_timescale(double /*minus_z*/,
                                                            void* /*parameters_and_workspace*/,
                                                            double* timescale) {
  *timescale = 1.;
}

// Single evolver derivative shim: integrates in minus_z = -z, dispatches on phase.
void ThermodynamicsModule::thermodynamics_recfast_derivs(double minus_z,
                                                         double* y,
                                                         double* dy,
                                                         void* parameters_and_workspace) {
  auto ws        = static_cast<thermodynamics_parameters_and_workspace*>(parameters_and_workspace);
  auto* module   = ws->thermodynamics_module;
  const double z = -minus_z;

  if (ws->recombination_phase == RecombinationPhase::helium) {
    // State is {x_He, Tmat}; hydrogen follows its Saha branch.
    double y_full[3] = {module->thermodynamics_recfast_hydrogen_saha_xH(ws->preco, z), y[0], y[1]};
    double dy_full_dz[3];
    module->thermodynamics_recombination_derivs_member(z,
                                                       y_full,
                                                       dy_full_dz,
                                                       parameters_and_workspace);
    dy[0] = -dy_full_dz[1];
    dy[1] = -dy_full_dz[2];
  }
  else {
    // State is {x_H, x_He, Tmat}.
    double dy_dz[_RECFAST_INTEG_SIZE_];
    module->thermodynamics_recombination_derivs_member(z, y, dy_dz, parameters_and_workspace);
    for (int index_y = 0; index_y < _RECFAST_INTEG_SIZE_; index_y++)
      dy[index_y] = -dy_dz[index_y];
  }
}

// Single dense-output shim: stores one recombination-table row, dispatching on phase.
void ThermodynamicsModule::thermodynamics_recfast_output(
    double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace) {
  auto ws      = static_cast<thermodynamics_parameters_and_workspace*>(parameters_and_workspace);
  auto* module = ws->thermodynamics_module;
  recombination* preco   = ws->preco;
  const double z         = -minus_z;
  const int sample_index = ws->recfast_output_index_offset + index_x;

  if (ws->recombination_phase == RecombinationPhase::helium) {
    const double y_full[3] = {module->thermodynamics_recfast_hydrogen_saha_xH(preco, z),
                              y[0],
                              y[1]};
    const double xe        = module->thermodynamics_recfast_xe_after_helium_ode(preco, z, y_full);
    module->thermodynamics_recfast_store_row(preco, sample_index, z, xe, y[1], -dy[1]);
  }
  else {
    const double xe = module->thermodynamics_recfast_xe_after_full_ode(preco, z, y);
    module->thermodynamics_recfast_store_row(preco, sample_index, z, xe, y[2], -dy[2]);
  }
}

// No-op dense-output shim: used when the helium ODE runs with no requested samples.
void ThermodynamicsModule::thermodynamics_recfast_output_none(double /*minus_z*/,
                                                              double /*y*/[],
                                                              double /*dy*/[],
                                                              int /*index_x*/,
                                                              void* /*parameters_and_workspace*/) {}

/**
 * Thermodynamics quantities at given redshift z.
 *
 * Evaluates all thermodynamics quantities at a given value of
 * the redshift by reading the pre-computed table and interpolating.
 *
 * @param z          Input: redshift
 * @param inter_mode Input: interpolation mode (normal or growing_closeby)
 * @param last_index Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback   Input: vector of background quantities (used only in case z>z_initial for getting ddkappa and dddkappa; in that case, should be already allocated and filled, with format short_info or larger; in other cases, will be ignored)
 * @param pvecthermo Output: vector of thermodynamics quantities (assumed to be already allocated)
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_at_z(
    double z, short inter_mode, int* last_index, double* pvecback, double* pvecthermo) const {
  /** Summary: */

  /** - define local variables */

  /* - the fact that z is in the pre-computed range 0 <= z <= z_initial
     will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z
     > z_intial: then, all relevant quantities can be extrapolated
     using simple analytic approximations */

  if (z >= z_table_[tt_size_ - 1]) {
    /* Seed the caller's interpolation cursor even though nothing is
       interpolated here. inter_normal is output-only for last_index, so callers
       legitimately use this call to initialise a cursor they then sweep with
       inter_closeby; leaving it unwritten hands the next hunt whatever was on
       the stack (#380). The last row is the one this branch extrapolates from,
       and the sweep proceeds downward in z from there --
       array_hunt_growing_closeby documents n_lines-1 as a valid hint. */
    *last_index = tt_size_ - 1;

    /* ionization fraction assumed to remain constant at large z */
    double x0                = thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_xe_];
    pvecthermo[index_th_xe_] = x0;

    /* Calculate dkappa/dtau (dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo[index_th_dkappa_] = (1. + z) * (1. + z) * n_e_ * x0 * _sigma_ * _Mpc_over_m_;

    /* tau_d scales like (1+z)**2 */
    pvecthermo[index_th_tau_d_] =
        thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_d_] *
        pow((1 + z) / (1. + z_table_[tt_size_ - 1]), 2);

    if (pth->compute_damping_scale) {
      /* r_d scales like (1+z)**-3/2 */
      pvecthermo[index_th_r_d_] = thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_r_d_] *
                                  pow((1 + z) / (1. + z_table_[tt_size_ - 1]), -1.5);
    }

    /* Calculate d2kappa/dtau2 = dz/dtau d/dz[dkappa/dtau] given that [dkappa/dtau] proportional to (1+z)^2 and dz/dtau = -H */
    pvecthermo[index_th_ddkappa_] = -pvecback[background_module_->index_bg_H_] * 2. / (1. + z) *
                                    pvecthermo[index_th_dkappa_];

    /* Calculate d3kappa/dtau3 given that [dkappa/dtau] proportional to (1+z)^2 */
    pvecthermo[index_th_dddkappa_] = (pvecback[background_module_->index_bg_H_] *
                                          pvecback[background_module_->index_bg_H_] / (1. + z) -
                                      pvecback[background_module_->index_bg_H_prime_]) *
                                     2. / (1. + z) * pvecthermo[index_th_dkappa_];

    /* \f$ exp^{-\kappa}, g, g', g'' \f$ can be set to zero: they are
       used only for computing the source functions in the
       perturbation module; but source functions only need to be
       sampled below z_initial (the condition that
       z_start_sources<z_initial is checked in the perturbation
       module) */
    pvecthermo[index_th_exp_m_kappa_] = 0.;
    pvecthermo[index_th_g_]           = 0.;
    pvecthermo[index_th_dg_]          = 0.;
    pvecthermo[index_th_ddg_]         = 0.;

    /* Calculate Tb */
    pvecthermo[index_th_Tb_] = pba->T_cmb * (1. + z);

    /* Calculate baryon equation of state parameter wb = (k_B/mu) Tb */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[index_th_wb_] = _k_B_ / (_c_ * _c_ * _m_H_) *
                               (1. + (1. / _not4_ - 1.) * YHe_ + x0 * (1. - YHe_)) * pba->T_cmb *
                               (1. + z);

    /* Calculate baryon adiabatic sound speed cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[index_th_cb2_] = pvecthermo[index_th_wb_] * 4. / 3.;

    /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
    if (pth->compute_cb2_derivatives) {
      /* since cb2 proportional to (1+z) or 1/a, its derivative wrt conformal time is given by dcb2 = - a H cb2 */
      pvecthermo[index_th_dcb2_] = -pvecback[background_module_->index_bg_H_] *
                                   pvecback[background_module_->index_bg_a_] *
                                   pvecthermo[index_th_cb2_];

      /* then its second derivative is given by ddcb2 = - a H' cb2 */
      pvecthermo[index_th_ddcb2_] = -pvecback[background_module_->index_bg_H_prime_] *
                                    pvecback[background_module_->index_bg_a_] *
                                    pvecthermo[index_th_cb2_];
    }

    /* in this regime, variation rate = dkappa/dtau */
    pvecthermo[index_th_rate_] = pvecthermo[index_th_dkappa_];

    /* quantities related to DM interacting with DR */
    if (all_species_.count("IDM_DR_IDR") > 0) {
      const auto& idm_dr_comp = static_cast<const IDM_DR_IDR_Species&>(
          *all_species_.at("IDM_DR_IDR"));
      const double Omega0_idm_dr_thermo = idm_dr_comp.idm_dr().GetOmega0();
      const double Omega0_idr_thermo    = idm_dr_comp.idr().GetOmega0();
      const double T_idr                = idm_dr_comp.idr().T_idr();
      /* calculate dmu_idm_dr and approximate its derivatives as zero */
      pvecthermo[index_th_dmu_idm_dr_]  = idm_dr_comp.idm_dr().a_idm_dr() *
                                          pow((1. + z) / 1.e7,
                                              idm_dr_comp.idm_dr().nindex_idm_dr()) *
                                          Omega0_idm_dr_thermo * pow(pba->h, 2);
      pvecthermo[index_th_ddmu_idm_dr_] = -pvecback[background_module_->index_bg_H_] *
                                          idm_dr_comp.idm_dr().nindex_idm_dr() / (1 + z) *
                                          pvecthermo[index_th_dmu_idm_dr_];
      pvecthermo[index_th_dddmu_idm_dr_] =
          (pvecback[background_module_->index_bg_H_] * pvecback[background_module_->index_bg_H_] /
               (1. + z) -
           pvecback[background_module_->index_bg_H_prime_]) *
          idm_dr_comp.idm_dr().nindex_idm_dr() / (1. + z) * pvecthermo[index_th_dmu_idm_dr_];

      /* calculate dmu_idr (self interaction) */
      pvecthermo[index_th_dmu_idr_] = idm_dr_comp.idr().b_idr() *
                                      pow((1. + z) / 1.e7, idm_dr_comp.idm_dr().nindex_idm_dr()) *
                                      Omega0_idr_thermo * pow(pba->h, 2);

      /* extrapolate optical depth of idm_dr and idr */
      pvecthermo[index_th_tau_idm_dr_] =
          thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idm_dr_] +
          (thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idm_dr_] -
           thermodynamics_table_[(tt_size_ - 2) * th_size_ + index_th_tau_idm_dr_]) *
              (z - z_table_[tt_size_ - 1]) / (z_table_[tt_size_ - 1] - z_table_[tt_size_ - 2]);

      pvecthermo[index_th_tau_idr_] =
          thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idr_] +
          (thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idr_] -
           thermodynamics_table_[(tt_size_ - 2) * th_size_ + index_th_tau_idr_]) *
              (z - z_table_[tt_size_ - 1]) / (z_table_[tt_size_ - 1] - z_table_[tt_size_ - 2]);

      /* extrapolate idm_dr visibility function */
      pvecthermo[index_th_g_idm_dr_] =
          thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_g_idm_dr_];

      /* calculate interacting dark matter sound speed */
      pvecthermo[index_th_cidm_dr2_] = 4 * _k_B_ * T_idr * (1. + z) / _eV_ / 3. /
                                       idm_dr_comp.idm_dr().m_idm();

      /* calculate interacting dark matter temperature (equal to idr temperature at this redhsift) */
      pvecthermo[index_th_Tidm_dr_] = T_idr * (1. + z);
    }
  }

  /** - interpolate in table with array_interpolate_spline() (normal
      mode) or array_interpolate_spline_growing_closeby() (closeby
      mode) */

  else {
    /* some very specific cases require linear interpolation because of a break in the derivative of the functions */
    if (((pth->reio_parametrization == reio_half_tanh) && (z < 2 * z_reionization_)) ||
        ((pth->reio_parametrization == reio_inter) && (z < 50.))) {
      array_interpolate_linear(z_table_.data(),
                               tt_size_,
                               thermodynamics_table_.data(),
                               th_size_,
                               z,
                               last_index,
                               pvecthermo,
                               th_size_);
    }

    /* in the "normal" case, use spline interpolation */
    else {
      if (inter_mode == inter_normal_) {
        array_interpolate_spline(z_table_.data(),
                                 tt_size_,
                                 thermodynamics_table_.data(),
                                 d2thermodynamics_dz2_table_.data(),
                                 th_size_,
                                 z,
                                 last_index,
                                 pvecthermo,
                                 th_size_);
      }

      if (inter_mode == inter_closeby_) {
        array_interpolate_spline_growing_closeby(z_table_.data(),
                                                 tt_size_,
                                                 thermodynamics_table_.data(),

                                                 d2thermodynamics_dz2_table_.data(),
                                                 th_size_,
                                                 z,
                                                 last_index,
                                                 pvecthermo,
                                                 th_size_);
      }
    }
  }
}

/**
 * Initialize the thermo structure, and in particular the
 * thermodynamics interpolation table.
 *
 * @return the error status
 */
void ThermodynamicsModule::thermodynamics_init() {
  /** Summary: */

  /** - define local variables */

  /* index running over time*/
  int index_tau;
  /* vector of background values for calling background_at_tau() */
  std::vector<double> pvecback;
  /* index for calling background_at_tau() */
  int last_index_back;
  /* temporary table of values of tau associated with z values in z_table_ */
  std::vector<double> tau_table;
  /* same ordered in growing time rather than growing redshift */
  std::vector<double> tau_table_growing;
  /* R = (3./4.)*(rho_b/rho_g) */
  double R;
  /* structures for storing temporarily information on recombination
     and reionization */
  struct recombination reco;
  struct reionization reio;
  struct recombination* preco;
  struct reionization* preio;

  double tau;

  tau_reionization_ = pth->tau_reio;
  z_reionization_   = pth->z_reio;
  YHe_              = pth->YHe;

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics");

  /** - compute and check primordial Helium fraction  */

  /* Y_He */
  if (YHe_ == _BBN_) {
    thermodynamics_helium_from_bbn();
    if (pth->thermodynamics_verbose > 0)
      printf(" with Y_He=%.4f\n", YHe_);
  }
  else {
    if (pth->thermodynamics_verbose > 0)
      printf("\n");
  }

  class_test((YHe_ < _YHE_SMALL_) || (YHe_ > _YHE_BIG_),
             "Y_He=%g out of bounds (%g<Y_He<%g)",
             YHe_,
             _YHE_SMALL_,
             _YHE_BIG_);

  /** - check energy injection parameters */

  class_test((pth->annihilation < 0), "annihilation parameter cannot be negative");

  class_test((pth->annihilation > 1.e-4),
             "annihilation parameter suspiciously large (%e, while typical bounds are in the range "
             "of 1e-7 to 1e-6)",
             pth->annihilation);

  class_test((pth->annihilation_variation > 0),
             "annihilation variation parameter must be negative (decreasing annihilation rate)");

  class_test((pth->annihilation_z < 0), "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation_zmin < 0),
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation_zmax < 0),
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation > 0) &&
                 (all_species_.count("CDM") == 0 && all_species_.count("IDM_DR_IDR") == 0),
             "CDM annihilation effects require the presence of CDM or IDM!");

  class_test((pth->annihilation_f_halo > 0) && (pth->recombination == recfast),
             "Switching on DM annihilation in halos requires using HyRec instead of RECFAST. "
             "Otherwise some values go beyond their range of validity in the RECFAST fits, and the "
             "thermodynamics module fails. Two solutions: add 'recombination = HyRec' to your "
             "input, or set 'annihilation_f_halo = 0.' (default).");

  class_test((pth->annihilation_f_halo < 0),
             "Parameter for DM annihilation in halos cannot be negative");

  class_test((pth->annihilation_z_halo < 0),
             "Parameter for DM annihilation in halos cannot be negative");

  if (pth->thermodynamics_verbose > 0)
    if ((pth->annihilation > 0) && (pth->reio_parametrization == reio_none) &&
        (ppr->recfast_Heswitch >= 3) && (pth->recombination == recfast))
      printf(
          "Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch "
          ">= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, "
          "producing nan's. This is however masked by reionization if you are not in reio_none "
          "mode.");

  class_test((pth->decay < 0), "decay parameter cannot be negative");

  class_test((pth->decay > 0) &&
                 (all_species_.count("CDM") == 0 && all_species_.count("IDM_DR_IDR") == 0),
             "CDM decay effects require the presence of CDM or IDM!");

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0., "stop to avoid division by zero");
  class_test(YHe_ == 1., "stop to avoid division by zero");

  /** - initialize pointers */

  preco = &reco;
  preio = &reio;

  /** - assign values to all indices in the structures with thermodynamics_indices()*/

  thermodynamics_indices(preco, preio);

  /** - allocate background vector */

  pvecback.resize(background_module_->bg_size_);

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_recombination() */

  thermodynamics_recombination(preco, pvecback.data());

  /** - if there is reionization, solve reionization and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_reionization()*/

  if (pth->reio_parametrization != reio_none) {
    thermodynamics_reionization(preco, preio, pvecback.data());
  }
  else {
    preio->rt_size                    = 0;
    preio->index_reco_when_reio_start = -1;
  }

  /** - merge tables in recombination and reionization structures into
      a single table in thermo structure */

  thermodynamics_merge_reco_and_reio(preco, preio);

  /** - compute table of corresponding conformal times */

  tau_table.resize(tt_size_);

  for (index_tau = 0; index_tau < tt_size_; index_tau++) {
    background_module_->background_tau_of_z(z_table_[index_tau], tau_table.data() + index_tau);
  }

  /** - store initial value of conformal time in the structure */

  tau_ini_ = tau_table[tt_size_ - 1];

  /** - fill missing columns (quantities not computed previously but related) */

  /** - --> minus the baryon drag interaction rate, -dkappa_d/dtau = -[1/R * kappa'], with R = 3 rho_b / 4 rho_gamma, stored temporarily in column ddkappa */

  // bt_size_ (ROWS = time steps), not bg_size_ (COLUMNS = entries per row). This is a
  // starting guess for the closeby hunt below, which walks a ROW index of the background
  // table, so it must lie in [0, bt_size_-1]; tau_table here is DECREASING (tau_ini_ is
  // its last element), so the sweep starts at the latest time and the top row is the
  // right guess.
  //
  // It read bg_size_-1 for years and only ever mattered once bg_size_ grew past bt_size_:
  // background columns scale with the momentum grids, so a decaying-NCDM sector with
  // dr_N_q=1600 pushed bg_size_ to 6481 against bt_size_=4605 and the run died in
  // array_hunt_growing_closeby with "*last_index=6480 out of range [0:4604]". Below that
  // threshold the wrong value is merely a poor guess that the hunt walks away from, which
  // is why no result ever changed -- and why this must NOT be "fixed" by clamping inside
  // the hunt, which would hide the next caller that confuses the two sizes.
  last_index_back = background_module_->bt_size_ - 1;

  for (index_tau = 0; index_tau < tt_size_; index_tau++) {
    background_module_->background_at_tau(tau_table[index_tau],
                                          pba->normal_info,
                                          pba->inter_closeby,
                                          &last_index_back,
                                          pvecback.data());

    R = 3. / 4. * all_species_.baryons().Rho(pvecback.data()) /
        all_species_.photons().Rho(pvecback.data());

    thermodynamics_table_[index_tau * th_size_ + index_th_ddkappa_] =
        -1. / R * thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_];

    if (all_species_.count("IDM_DR_IDR") > 0) {
      auto& idm_idr = static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
      const double Omega0_idm_dr_table = idm_idr.idm_dr().GetOmega0();
      const double Omega0_idr_table    = idm_idr.idr().GetOmega0();
      /* - --> idr interaction rate with idm_dr (i.e. idr opacity to idm_dr scattering) */
      thermodynamics_table_[index_tau * th_size_ + index_th_dmu_idm_dr_] =
          idm_idr.idm_dr().a_idm_dr() *
          pow((1. + z_table_[index_tau]) / 1.e7, idm_idr.idm_dr().nindex_idm_dr()) *
          Omega0_idm_dr_table * pow(pba->h, 2);

      /* - --> idm_dr interaction rate with idr (i.e. idm_dr opacity
               to idr scattering), [Sinv*dmu_idm_dr] with Sinv = (4
               rho_idr) / (3 rho_idm_dr), stored temporarily in
               ddmu_idm_dr */
      {
        const double rho_idr    = idm_idr.idr().Rho(pvecback.data());
        const double rho_idm_dr = idm_idr.idm_dr().Rho(pvecback.data());
        if (rho_idr > 0. && rho_idm_dr > 0.) {
          thermodynamics_table_[index_tau * th_size_ + index_th_ddmu_idm_dr_] =
              4. / 3. * rho_idr / rho_idm_dr *
              thermodynamics_table_[index_tau * th_size_ + index_th_dmu_idm_dr_];
        }
        else {
          thermodynamics_table_[index_tau * th_size_ + index_th_ddmu_idm_dr_] = 0.;
        }
      }

      /* - --> idr self-interaction rate */
      thermodynamics_table_[index_tau * th_size_ + index_th_dmu_idr_] =
          idm_idr.idr().b_idr() *
          pow((1. + z_table_[index_tau]) / 1.e7, idm_idr.idm_dr().nindex_idm_dr()) *
          Omega0_idr_table * pow(pba->h, 2);
    }
  }

  /** - --> second derivative of this rate, -[1/R * kappa']'', stored temporarily in column dddkappa */
  array_spline_table_line_to_line(tau_table.data(),
                                  tt_size_,
                                  thermodynamics_table_.data(),
                                  th_size_,
                                  index_th_ddkappa_,
                                  index_th_dddkappa_,
                                  _SPLINE_EST_DERIV_);

  /** - --> compute tau_d = [int_{tau_today}^{tau} dtau -dkappa_d/dtau] */
  array_integrate_spline_table_line_to_line(tau_table.data(),
                                            tt_size_,
                                            thermodynamics_table_.data(),
                                            th_size_,
                                            index_th_ddkappa_,
                                            index_th_dddkappa_,
                                            index_th_tau_d_);

  /* the temporary quantities stored in columns ddkappa and dddkappa
     will not be used anymore, so they can be overwritten by other
     intermediate steps of other computations */

  if (all_species_.count("IDM_DR_IDR") > 0) {
    /** --> second derivative of idm_dr interaction rate (with idr), [Sinv*dmu_idm_dr]'', stored temporarily in column dddmu */
    array_spline_table_line_to_line(tau_table.data(),
                                    tt_size_,
                                    thermodynamics_table_.data(),
                                    th_size_,
                                    index_th_ddmu_idm_dr_,
                                    index_th_dddmu_idm_dr_,
                                    _SPLINE_EST_DERIV_);

    /** - --> compute optical depth of idm, tau_idm_dr = [int_{tau_today}^{tau} dtau [Sinv*dmu_idm_dr] ].
              This step gives -tau_idm_dr. The resulty is mutiplied by -1 later on. */
    array_integrate_spline_table_line_to_line(tau_table.data(),
                                              tt_size_,
                                              thermodynamics_table_.data(),
                                              th_size_,
                                              index_th_ddmu_idm_dr_,
                                              index_th_dddmu_idm_dr_,
                                              index_th_tau_idm_dr_);

    /** - --> second derivative of idr interaction rate (with idm_dr), [dmu_idm_idr]'', stored temporarily in column dddmu */
    array_spline_table_line_to_line(tau_table.data(),
                                    tt_size_,
                                    thermodynamics_table_.data(),
                                    th_size_,
                                    index_th_dmu_idm_dr_,
                                    index_th_dddmu_idm_dr_,
                                    _SPLINE_EST_DERIV_);

    /** - --> compute optical depth of idr, tau_idr = [int_{tau_today}^{tau} dtau [dmu_idm_idr] ].
              This step gives -tau_idr. The resulty is mutiplied by -1 later on. */
    array_integrate_spline_table_line_to_line(tau_table.data(),
                                              tt_size_,
                                              thermodynamics_table_.data(),
                                              th_size_,
                                              index_th_dmu_idm_dr_,
                                              index_th_dddmu_idm_dr_,
                                              index_th_tau_idr_);
  }

  /** - --> compute damping scale:

      r_d = 2pi/k_d = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') 1/6 (R^2+16/15(1+R))/(1+R)^2]^1/2
      = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') 1/6 (R^2/(1+R)+16/15)/(1+R)]^1/2

      which is like in CosmoTherm (CT), but slightly
      different from Wayne Hu (WH)'s thesis eq. (5.59):
      the factor 16/15 in CT is 4/5 in WH */

  if (pth->compute_damping_scale) {
    tau_table_growing.resize(tt_size_);

    /* compute integrand and store temporarily in column "ddkappa" */
    for (index_tau = 0; index_tau < tt_size_; index_tau++) {
      tau_table_growing[index_tau] = tau_table[tt_size_ - 1 - index_tau];

      background_module_->background_at_tau(tau_table_growing[index_tau],
                                            pba->normal_info,
                                            pba->inter_closeby,
                                            &last_index_back,
                                            pvecback.data());

      R = 3. / 4. * all_species_.baryons().Rho(pvecback.data()) /
          all_species_.photons().Rho(pvecback.data());

      thermodynamics_table_[index_tau * th_size_ + index_th_ddkappa_] =
          1. / 6. /
          thermodynamics_table_[(tt_size_ - 1 - index_tau) * th_size_ + index_th_dkappa_] *
          (R * R / (1 + R) + 16. / 15.) / (1. + R);
    }

    /* compute second derivative of integrand, and store temporarily in column "dddkappa" */
    array_spline_table_line_to_line(tau_table_growing.data(),
                                    tt_size_,
                                    thermodynamics_table_.data(),
                                    th_size_,
                                    index_th_ddkappa_,
                                    index_th_dddkappa_,
                                    _SPLINE_EST_DERIV_);

    /* compute integral and store temporarily in column "g" */
    array_integrate_spline_table_line_to_line(tau_table_growing.data(),
                                              tt_size_,
                                              thermodynamics_table_.data(),
                                              th_size_,
                                              index_th_ddkappa_,
                                              index_th_dddkappa_,
                                              index_th_g_);

    /* we could now write the result as r_d = 2pi * sqrt(integral),
       but we will first better acount for the contribution frokm the tau_ini boundary.
       Close to this boundary, R=0 and the integrand is just 16/(15*6)/kappa'
       Using kappa' propto 1/a^2 and tau propro a during RD, we get the analytic result:
       int_0^{tau_ini} dtau / kappa' = tau_ini / 3 / kappa'_ini
       Thus r_d = 2pi * sqrt( 16/(15*6*3) * (tau_ini/ kappa'_ini) * integral) */

    double tau_ini    = tau_table[tt_size_ - 1];
    double dkappa_ini = thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_dkappa_];

    for (index_tau = 0; index_tau < tt_size_; index_tau++) {
      thermodynamics_table_[index_tau * th_size_ + index_th_r_d_] =
          2. * _PI_ *
          sqrt(16. / (15. * 6. * 3.) * tau_ini / dkappa_ini +
               thermodynamics_table_[(tt_size_ - 1 - index_tau) * th_size_ + index_th_g_]);
    }

  }  // end of damping scale calculation

  /** - --> second derivative with respect to tau of dkappa (in view of spline interpolation) */
  array_spline_table_line_to_line(tau_table.data(),
                                  tt_size_,
                                  thermodynamics_table_.data(),
                                  th_size_,
                                  index_th_dkappa_,
                                  index_th_dddkappa_,
                                  _SPLINE_EST_DERIV_);

  /** - --> first derivative with respect to tau of dkappa (using spline interpolation) */
  array_derive_spline_table_line_to_line(tau_table.data(),
                                         tt_size_,
                                         thermodynamics_table_.data(),
                                         th_size_,
                                         index_th_dkappa_,
                                         index_th_dddkappa_,
                                         index_th_ddkappa_);

  /** - --> compute -kappa = [int_{tau_today}^{tau} dtau dkappa/dtau], store temporarily in column "g" */
  array_integrate_spline_table_line_to_line(tau_table.data(),
                                            tt_size_,
                                            thermodynamics_table_.data(),
                                            th_size_,
                                            index_th_dkappa_,
                                            index_th_dddkappa_,
                                            index_th_g_);

  /** - --> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives) {
    /** - ---> second derivative with respect to tau of cb2 */
    array_spline_table_line_to_line(tau_table.data(),
                                    tt_size_,
                                    thermodynamics_table_.data(),
                                    th_size_,
                                    index_th_cb2_,
                                    index_th_ddcb2_,
                                    _SPLINE_EST_DERIV_);

    /** - ---> first derivative with respect to tau of cb2 (using spline interpolation) */
    array_derive_spline_table_line_to_line(tau_table.data(),
                                           tt_size_,
                                           thermodynamics_table_.data(),
                                           th_size_,
                                           index_th_cb2_,
                                           index_th_ddcb2_,
                                           index_th_dcb2_);
  }

  /** - --> compute visibility: \f$ g= (d \kappa/d \tau) e^{- \kappa} \f$ */

  /* loop on z (decreasing z, increasing time) */
  for (index_tau = tt_size_ - 1; index_tau >= 0; index_tau--) {
    /** - ---> compute g */
    double g = thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_] *
               exp(thermodynamics_table_[index_tau * th_size_ + index_th_g_]);

    /** - ---> compute exp(-kappa) */
    thermodynamics_table_[index_tau * th_size_ + index_th_exp_m_kappa_] = exp(
        thermodynamics_table_[index_tau * th_size_ + index_th_g_]);

    /** - ---> compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    thermodynamics_table_[index_tau * th_size_ + index_th_dg_] =
        (thermodynamics_table_[index_tau * th_size_ + index_th_ddkappa_] +
         thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_] *
             thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_]) *
        exp(thermodynamics_table_[index_tau * th_size_ + index_th_g_]);

    /** - ---> compute g''  */
    thermodynamics_table_[index_tau * th_size_ + index_th_ddg_] =
        (thermodynamics_table_[index_tau * th_size_ + index_th_dddkappa_] +
         thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_] *
             thermodynamics_table_[index_tau * th_size_ + index_th_ddkappa_] * 3. +
         thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_] *
             thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_] *
             thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_]) *
        exp(thermodynamics_table_[index_tau * th_size_ + index_th_g_]);

    /** - ---> store g */
    thermodynamics_table_[index_tau * th_size_ + index_th_g_] = g;

    /** - ---> compute variation rate */
    class_test(thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_] == 0.,
               "variation rate diverges");

    thermodynamics_table_[index_tau * th_size_ + index_th_rate_] = sqrt(
        pow(thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_], 2) +
        pow(thermodynamics_table_[index_tau * th_size_ + index_th_ddkappa_] /
                thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_],
            2) +
        fabs(thermodynamics_table_[index_tau * th_size_ + index_th_dddkappa_] /
             thermodynamics_table_[index_tau * th_size_ + index_th_dkappa_]));

    /* - ---> restore correct sign for idm_dr and idr optical depth, and calculate idm_dr visibility function */
    if (all_species_.count("IDM_DR_IDR") > 0) {
      /* restore the correct sign for tau_idm_dr */
      thermodynamics_table_[index_tau * th_size_ + index_th_tau_idm_dr_] *= -1.;

      /* restore the correct sign for tau_idr */
      thermodynamics_table_[index_tau * th_size_ + index_th_tau_idr_] *= -1.;

      /* visibility function for idm_dr : g_idm_dr = [Sinv*dmu_idm_dr] * exp(-tau_idm_dr) */
      thermodynamics_table_[index_tau * th_size_ + index_th_g_idm_dr_] =
          thermodynamics_table_[index_tau * th_size_ + index_th_ddmu_idm_dr_] *
          exp(-thermodynamics_table_[index_tau * th_size_ + index_th_tau_idm_dr_]);
    }
  }

  /** - smooth the rate (details of smoothing unimportant: only the
      order of magnitude of the rate matters) */
  array_smooth(thermodynamics_table_.data(),
               th_size_,
               tt_size_,
               index_th_rate_,
               ppr->thermo_rate_smoothing_radius);

  /* - ---> fill columns for ddmu_idm_dr and dddmu_idm_dr with true values, and compute idm_dr temperature and sound speed */
  if (all_species_.count("IDM_DR_IDR") > 0) {
    const auto& comp_ref = static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
    const double Omega0_idr_local = comp_ref.idr().GetOmega0();
    const double T_idr_param      = comp_ref.idr().T_idr();
    double Gamma_heat_idm_dr, dTdz_idm_dr, T_idm_dr, T_idr, dz, T_adia, z_adia;
    double z;

    /** - --> second derivative with respect to tau of dmu_idm_dr (in view of spline interpolation) */
    array_spline_table_line_to_line(tau_table.data(),
                                    tt_size_,
                                    thermodynamics_table_.data(),
                                    th_size_,
                                    index_th_dmu_idm_dr_,
                                    index_th_dddmu_idm_dr_,
                                    _SPLINE_EST_DERIV_);

    /** - --> first derivative with respect to tau of dmu_idm_dr (using spline interpolation) */
    array_derive_spline_table_line_to_line(tau_table.data(),
                                           tt_size_,
                                           thermodynamics_table_.data(),
                                           th_size_,
                                           index_th_dmu_idm_dr_,
                                           index_th_dddmu_idm_dr_,
                                           index_th_ddmu_idm_dr_);

    /** - --> now compute idm_dr temperature and sound speed in various regimes */

    /* (A) - initial value of T_idm_dr at the maximum z (minimum tau) */

    z = z_table_[tt_size_ - 1];

    background_module_->background_tau_of_z(z, &(tau));

    background_module_->background_at_tau(tau,
                                          pba->short_info,
                                          pba->inter_normal,
                                          &last_index_back,
                                          pvecback.data());

    Gamma_heat_idm_dr = 2. * Omega0_idr_local * pow(pba->h, 2) * comp_ref.idm_dr().a_idm_dr() *
                        pow((1. + z), (comp_ref.idm_dr().nindex_idm_dr() + 1.)) /
                        pow(1.e7, comp_ref.idm_dr().nindex_idm_dr());

    /* (A1) --> if Gamma is not much smaller than H, set T_idm_dr to T_idm_dr = T_idr = xi*T_gamma (tight coupling solution) */
    if (Gamma_heat_idm_dr > 1.e-3 * pvecback[background_module_->index_bg_a_] *
                                pvecback[background_module_->index_bg_H_]) {
      T_idr       = T_idr_param * (1. + z);
      T_idm_dr    = T_idr;
      dTdz_idm_dr = T_idr_param;
    }

    /* (A2) --> otherwise, if Gamma << H, set initial T_idm_dr to the
       approximate analytic solution (Gamma/aH)/(1+(Gamma/aH)*T_idr)
       (eq. (A62) in ETHOS I ) */
    else {
      T_idr       = T_idr_param * (1. + z);
      T_idm_dr    = Gamma_heat_idm_dr /
                    (pvecback[background_module_->index_bg_a_] *
                     pvecback[background_module_->index_bg_H_]) /
                    (1. + Gamma_heat_idm_dr / (pvecback[background_module_->index_bg_a_] *
                                               pvecback[background_module_->index_bg_H_])) *
                    T_idr;
      dTdz_idm_dr = 2. * T_idm_dr - Gamma_heat_idm_dr / pvecback[background_module_->index_bg_H_] *
                                        (T_idr - T_idm_dr);
    }

    thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_Tidm_dr_] = T_idm_dr;
    thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_cidm_dr2_] =
        _k_B_ * T_idm_dr / _eV_ / comp_ref.idm_dr().m_idm() * (1. + dTdz_idm_dr / 3. / T_idm_dr);

    /* T_adia and z_adia will be used later. They are defined as "the
       last T_idm_dr(z) at which the temperature was evaluated
       explicitely, rather than scaled like a^{-2} (decoupled DM
       regime)". Here we just initialize them. They will be updated
       each time that we recompte T_idm_dr explicitely. */
    T_adia = T_idm_dr;
    z_adia = z;

    /* (B) - iterate over growing tau / decreasing z to find other
       values. At each new z we need to compute the following
       quantities: T_idr, T_idm_dr, Gamma_heat_idm_dr, a, H, dT_idm_dr,/dz,
       c_s_idm_dr^2. They all needed to be known from step to step, even
       if the final goal is only to store T_idm_dr, c_s_idm^2 */
    for (index_tau = tt_size_ - 2; index_tau >= 0; index_tau--) {
      /* (B1) --> tight-coupling solution: Gamma >> H implies T_idm_dr=T_idr=xi*T_gamma */
      if (Gamma_heat_idm_dr > 1.e3 * pvecback[background_module_->index_bg_a_] *
                                  pvecback[background_module_->index_bg_H_]) {
        z                 = z_table_[index_tau];
        T_idr             = T_idr_param * (1. + z);
        T_idm_dr          = T_idr;
        Gamma_heat_idm_dr = 2. * Omega0_idr_local * pow(pba->h, 2) * comp_ref.idm_dr().a_idm_dr() *
                            pow((1. + z), (comp_ref.idm_dr().nindex_idm_dr() + 1.)) /
                            pow(1.e7, comp_ref.idm_dr().nindex_idm_dr());
        background_module_->background_tau_of_z(z, &(tau));
        background_module_->background_at_tau(tau,
                                              pba->short_info,
                                              pba->inter_normal,
                                              &last_index_back,
                                              pvecback.data());
        dTdz_idm_dr = T_idr_param;
      }

      /* (B2) --> intermediate solution: integrate differential equation equation dT_idm_dr/dz = 2 a T_DM - Gamma/H (T_idr - T_idm_dr) */
      else if (Gamma_heat_idm_dr > 1.e-3 * pvecback[background_module_->index_bg_a_] *
                                       pvecback[background_module_->index_bg_H_]) {
        dz = z_table_[index_tau + 1] - z_table_[index_tau];

        /* (B2a) ----> if dz << H/Gamma the equation is not too stiff and the traditional forward Euler method converges */
        if (dz < pvecback[background_module_->index_bg_H_] / Gamma_heat_idm_dr / 10.) {
          z                  = z_table_[index_tau];
          T_idr              = T_idr_param * (1. + z);
          T_idm_dr          -= dTdz_idm_dr * dz;
          Gamma_heat_idm_dr  = 2. * Omega0_idr_local * pow(pba->h, 2) *
                               comp_ref.idm_dr().a_idm_dr() *
                               pow((1. + z), (comp_ref.idm_dr().nindex_idm_dr() + 1.)) /
                               pow(1.e7, comp_ref.idm_dr().nindex_idm_dr());
          background_module_->background_tau_of_z(z, &(tau));
          background_module_->background_at_tau(tau,
                                                pba->short_info,
                                                pba->inter_normal,
                                                &last_index_back,
                                                pvecback.data());
          dTdz_idm_dr = 2. * pvecback[background_module_->index_bg_a_] * T_idm_dr -
                        Gamma_heat_idm_dr / (pvecback[background_module_->index_bg_H_]) *
                            (T_idr - T_idm_dr);
        }

        /* (B2b) ----> otherwise, the equation is too stiff and the
           traditional forward Euler method diverges with this
           stepsize. But we can just decreasee dz to bring it back
           well within the convergence radius H/Gamma of the
           equation. */
        else {
          int N_sub_steps =
              (int) (dz / (pvecback[background_module_->index_bg_H_] / Gamma_heat_idm_dr / 10.)) +
              1;
          double dz_sub_step = dz / N_sub_steps;

          /* loop over sub-steps */
          for (int n = 0; n < N_sub_steps; n++) {
            /* evolve quantities over  sub-step wioth forward Euler method */

            z -= dz_sub_step;
            /* final redshift last sub-step overwritten to avoid small rounding error */
            if (n == (N_sub_steps - 1))
              z = z_table_[index_tau];

            T_idr              = T_idr_param * (1. + z);
            T_idm_dr          -= dTdz_idm_dr * dz_sub_step;
            Gamma_heat_idm_dr  = 2. * Omega0_idr_local * pow(pba->h, 2) *
                                 comp_ref.idm_dr().a_idm_dr() *
                                 pow((1. + z), (comp_ref.idm_dr().nindex_idm_dr() + 1.)) /
                                 pow(1.e7, comp_ref.idm_dr().nindex_idm_dr());
            background_module_->background_tau_of_z(z, &(tau));
            background_module_->background_at_tau(tau,
                                                  pba->short_info,
                                                  pba->inter_normal,
                                                  &last_index_back,
                                                  pvecback.data());
            dTdz_idm_dr = 2. * pvecback[background_module_->index_bg_a_] * T_idm_dr -
                          Gamma_heat_idm_dr / (pvecback[background_module_->index_bg_H_]) *
                              (T_idr - T_idm_dr);
          }
        }

        /* update T_adia, z_adia */
        T_adia = T_idm_dr;
        z_adia = z;
      }

      /* (B3) --> decoupled solution: T_idm_dr scales like a^-2 */
      else {
        z                 = z_table_[index_tau];
        T_idr             = T_idr_param * (1. + z);
        T_idm_dr          = T_adia * pow((1. + z) / (1. + z_adia), 2);
        Gamma_heat_idm_dr = 2. * Omega0_idr_local * pow(pba->h, 2) * comp_ref.idm_dr().a_idm_dr() *
                            pow((1. + z), (comp_ref.idm_dr().nindex_idm_dr() + 1.)) /
                            pow(1.e7, comp_ref.idm_dr().nindex_idm_dr());
        background_module_->background_tau_of_z(z, &(tau));
        background_module_->background_at_tau(tau,
                                              pba->short_info,
                                              pba->inter_normal,
                                              &last_index_back,
                                              pvecback.data());
        dTdz_idm_dr = 2. / (1 + z) * T_idm_dr;
      }

      thermodynamics_table_[index_tau * th_size_ + index_th_Tidm_dr_] = T_idm_dr;
      thermodynamics_table_[index_tau * th_size_ + index_th_cidm_dr2_] =
          _k_B_ * T_idm_dr / _eV_ / comp_ref.idm_dr().m_idm() * (1. + dTdz_idm_dr / 3. / T_idm_dr);
    }
  }

  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */

  array_spline_table_lines(z_table_.data(),
                           tt_size_,
                           thermodynamics_table_.data(),
                           th_size_,
                           d2thermodynamics_dz2_table_.data(),
                           _SPLINE_EST_DERIV_);

  /** - find maximum of g */

  index_tau = tt_size_ - 1;
  while (z_table_[index_tau] > _Z_REC_MAX_) {
    index_tau--;
  }

  class_test(thermodynamics_table_[(index_tau + 1) * th_size_ + index_th_g_] >
                 thermodynamics_table_[index_tau * th_size_ + index_th_g_],
             "found a recombination redshift greater or equal to the maximum value imposed in "
             "thermodynamics.h, z_rec_max=%g",
             _Z_REC_MAX_);

  while (thermodynamics_table_[(index_tau + 1) * th_size_ + index_th_g_] <
         thermodynamics_table_[index_tau * th_size_ + index_th_g_]) {
    index_tau--;
  }

  double g_max      = thermodynamics_table_[index_tau * th_size_ + index_th_g_];
  int index_tau_max = index_tau;

  /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
  z_rec_ = z_table_[index_tau + 1] +
           0.5 * (z_table_[index_tau + 1] - z_table_[index_tau]) *
               (thermodynamics_table_[(index_tau) *th_size_ + index_th_g_] -
                1. * thermodynamics_table_[(index_tau + 2) * th_size_ + index_th_g_]) /
               (thermodynamics_table_[(index_tau) *th_size_ + index_th_g_] -
                2. * thermodynamics_table_[(index_tau + 1) * th_size_ + index_th_g_] +
                thermodynamics_table_[(index_tau + 2) * th_size_ + index_th_g_]);

  class_test(z_rec_ + ppr->smallest_allowed_variation >= _Z_REC_MAX_,
             "found a recombination redshift greater or equal to the maximum value imposed in "
             "thermodynamics.h, z_rec_max=%g",
             _Z_REC_MAX_);

  class_test(z_rec_ - ppr->smallest_allowed_variation <= _Z_REC_MIN_,
             "found a recombination redshift smaller or equal to the maximum value imposed in "
             "thermodynamics.h, z_rec_min=%g",
             _Z_REC_MIN_);

  /** - find conformal recombination time using background_tau_of_z() **/

  background_module_->background_tau_of_z(z_rec_, &(tau_rec_));

  background_module_->background_at_tau(tau_rec_,
                                        pba->long_info,
                                        pba->inter_normal,
                                        &last_index_back,
                                        pvecback.data());

  rs_rec_            = pvecback[background_module_->index_bg_rs_];
  ds_rec_            = rs_rec_ / (1. + z_rec_);
  da_rec_            = pvecback[background_module_->index_bg_ang_distance_];
  ra_rec_            = da_rec_ * (1. + z_rec_);
  angular_rescaling_ = ra_rec_ / (background_module_->conformal_age_ - tau_rec_);

  /** - find damping scale at recombination (using linear interpolation) */

  if (pth->compute_damping_scale) {
    rd_rec_ = (z_table_[index_tau + 1] - z_rec_) / (z_table_[index_tau + 1] - z_table_[index_tau]) *
                  thermodynamics_table_[(index_tau) *th_size_ + index_th_r_d_] +
              (z_rec_ - z_table_[index_tau]) / (z_table_[index_tau + 1] - z_table_[index_tau]) *
                  thermodynamics_table_[(index_tau + 1) * th_size_ + index_th_r_d_];
  }

  /** - find time (always after recombination) at which tau_c/tau
      falls below some threshold, defining tau_free_streaming */

  background_module_->background_tau_of_z(z_table_[index_tau], &tau);

  while ((1. / thermodynamics_table_[(index_tau) *th_size_ + index_th_dkappa_] / tau <
          ppr->radiation_streaming_trigger_tau_c_over_tau) &&
         (index_tau > 0)) {
    index_tau--;

    background_module_->background_tau_of_z(z_table_[index_tau], &tau);
  }

  tau_free_streaming_ = tau;

  /** - Find interacting dark radiation free-streaming time */
  int index_tau_fs = index_tau;

  double tau_idm_dr_fs = 0.;

  if (all_species_.count("IDM_DR_IDR") > 0) {
    const auto& comp = static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
    if (comp.has_idm_dr()) {
      if (comp.idm_dr().nindex_idm_dr() >= 2) {
        index_tau = index_tau_fs - 1;
        /* comment: using index_tau_max (index_tau_fs) instead of tt_size_ - 1 ensures that the switch is always after recombination (free streaming) */
      }
      else {
        index_tau = 0;
      }

      background_module_->background_tau_of_z(z_table_[index_tau], &tau);

      while ((1. / thermodynamics_table_[(index_tau) *th_size_ + index_th_dmu_idm_dr_] / tau <
              ppr->idr_streaming_trigger_tau_c_over_tau) &&
             ((comp.idm_dr().nindex_idm_dr() >= 2 && index_tau > 0) ||
              (comp.idm_dr().nindex_idm_dr() < 2 && index_tau < tt_size_ - 1))) {
        if (comp.idm_dr().nindex_idm_dr() >= 2) {
          index_tau--;
        }
        else {
          index_tau++;
        }

        background_module_->background_tau_of_z(z_table_[index_tau], &tau);
      }

      tau_idm_dr_fs           = tau;
      tau_idr_free_streaming_ = tau;
    }

    /* case of idr alone without idm_dr */
    else {
      index_tau = index_tau_fs - 1;
      background_module_->background_tau_of_z(z_table_[index_tau], &tau);
      tau_idm_dr_fs           = tau;
      tau_idr_free_streaming_ = tau;
    }
  }

  /** - find z_star (when optical depth kappa crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau = 0;
  while ((thermodynamics_table_[(index_tau) *th_size_ + index_th_exp_m_kappa_] > 1. / _E_) &&
         (index_tau < tt_size_))
    index_tau++;

  z_star_ = z_table_[index_tau - 1] +
            (1. / _E_ - thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_exp_m_kappa_]) /
                (thermodynamics_table_[(index_tau) *th_size_ + index_th_exp_m_kappa_] -
                 thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_exp_m_kappa_]) *
                (z_table_[index_tau] - z_table_[index_tau - 1]);

  background_module_->background_tau_of_z(z_star_, &(tau_star_));

  background_module_->background_at_tau(tau_star_,
                                        pba->long_info,
                                        pba->inter_normal,
                                        &last_index_back,
                                        pvecback.data());

  rs_star_ = pvecback[background_module_->index_bg_rs_];
  ds_star_ = rs_star_ / (1. + z_star_);
  da_star_ = pvecback[background_module_->index_bg_ang_distance_];
  ra_star_ = da_star_ * (1. + z_star_);

  if (pth->compute_damping_scale) {
    rd_star_ = (z_table_[index_tau + 1] - z_star_) /
                   (z_table_[index_tau + 1] - z_table_[index_tau]) *
                   thermodynamics_table_[(index_tau) *th_size_ + index_th_r_d_] +
               (z_star_ - z_table_[index_tau]) / (z_table_[index_tau + 1] - z_table_[index_tau]) *
                   thermodynamics_table_[(index_tau + 1) * th_size_ + index_th_r_d_];
  }

  /** - find baryon drag time (when tau_d crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau = 0;
  while ((thermodynamics_table_[(index_tau) *th_size_ + index_th_tau_d_] < 1.) &&
         (index_tau < tt_size_))
    index_tau++;

  z_d_ = z_table_[index_tau - 1] +
         (1. - thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_tau_d_]) /
             (thermodynamics_table_[(index_tau) *th_size_ + index_th_tau_d_] -
              thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_tau_d_]) *
             (z_table_[index_tau] - z_table_[index_tau - 1]);

  background_module_->background_tau_of_z(z_d_, &(tau_d_));

  background_module_->background_at_tau(tau_d_,
                                        pba->long_info,
                                        pba->inter_normal,
                                        &last_index_back,
                                        pvecback.data());

  rs_d_ = pvecback[background_module_->index_bg_rs_];
  ds_d_ = rs_d_ / (1. + z_d_);

  /** - find idm_dr and idr drag times */
  double tau_idm_dr = 0.0, tau_idr = 0.0;
  if (all_species_.count("IDM_DR_IDR") > 0) {
    if ((thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idm_dr_] > 1.) &&
        (thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idr_] > 1.)) {
      index_tau = 0;

      while ((thermodynamics_table_[(index_tau) *th_size_ + index_th_tau_idm_dr_] < 1.) &&
             (index_tau < tt_size_ - 1))
        index_tau++;

      double z_idm_dr =
          z_table_[index_tau - 1] +
          (1. - thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_tau_idm_dr_]) /
              (thermodynamics_table_[(index_tau) *th_size_ + index_th_tau_idm_dr_] -
               thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_tau_idm_dr_]) *
              (z_table_[index_tau] - z_table_[index_tau - 1]);

      background_module_->background_tau_of_z(z_idm_dr, &(tau_idm_dr));

      index_tau = 0;

      while ((thermodynamics_table_[(index_tau) *th_size_ + index_th_tau_idr_] < 1.) &&
             (index_tau < tt_size_ - 1))
        index_tau++;

      double z_idr = z_table_[index_tau - 1] +
                     (1. - thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_tau_idr_]) /
                         (thermodynamics_table_[(index_tau) *th_size_ + index_th_tau_idr_] -
                          thermodynamics_table_[(index_tau - 1) * th_size_ + index_th_tau_idr_]) *
                         (z_table_[index_tau] - z_table_[index_tau - 1]);

      background_module_->background_tau_of_z(z_idr, &(tau_idr));
    }
  }

  /** - find time above which visibility falls below a given fraction of its maximum */

  index_tau = index_tau_max;
  while ((thermodynamics_table_[(index_tau) *th_size_ + index_th_g_] >
          g_max * ppr->neglect_CMB_sources_below_visibility) &&
         (index_tau > 0))
    index_tau--;

  background_module_->background_tau_of_z(z_table_[index_tau], &(tau_cut_));

  /** - if verbose flag set to next-to-minimum value, print the main results */

  if (pth->thermodynamics_verbose > 0) {
    if (all_species_.count("IDM_DR_IDR") > 0) {
      if ((thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idm_dr_] > 1.) &&
          (thermodynamics_table_[(tt_size_ - 1) * th_size_ + index_th_tau_idr_] > 1.)) {
        printf(" -> idr decouples at tau_idr = %e Mpc\n", tau_idr);
        printf(" -> idm_dr decouples at tau_idm_dr = %e Mpc\n", tau_idm_dr);
      }
      else {
        printf(
            " -> computation of decoupling time of idm_dr and idr skipped, because z would not be "
            "in z_table\n");
      }
    }
    printf(" -> recombination at z = %f (max of visibility function)\n", z_rec_);
    printf("    corresponding to conformal time = %f Mpc\n", tau_rec_);
    printf("    with comoving sound horizon = %f Mpc\n", rs_rec_);
    printf("    angular diameter distance = %f Mpc\n", da_rec_);
    printf("    and sound horizon angle 100*theta_s = %f\n", 100. * rs_rec_ / ra_rec_);
    if (pth->compute_damping_scale) {
      printf("    and with comoving photon damping scale = %f Mpc\n", rd_rec_);
      printf("    or comoving damping wavenumber k_d = %f 1/Mpc\n", 2. * _PI_ / rd_rec_);
    }
    printf("    Thomson optical depth crosses one at z_* = %f\n", z_star_);
    printf("    giving an angle 100*theta_* = %f\n", 100. * rs_star_ / ra_star_);
    printf(" -> baryon drag stops at z = %f\n", z_d_);
    printf("    corresponding to conformal time = %f Mpc\n", tau_d_);
    printf("    with comoving sound horizon rs = %f Mpc\n", rs_d_);
    if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
      if (pth->reio_z_or_tau == reio_tau)
        printf(" -> reionization  at z = %f\n", z_reionization_);
      if (pth->reio_z_or_tau == reio_z)
        printf(" -> reionization with optical depth = %f\n", tau_reionization_);
      double tau_reio;
      background_module_->background_tau_of_z(z_reionization_, &tau_reio);
      printf("    corresponding to conformal time = %f Mpc\n", tau_reio);
    }
    if (pth->reio_parametrization == reio_bins_tanh) {
      printf(" -> binned reionization gives optical depth = %f\n", tau_reionization_);
    }
    if (pth->reio_parametrization == reio_many_tanh) {
      printf(" -> many-step reionization gives optical depth = %f\n", tau_reionization_);
    }
    if (pth->reio_parametrization == reio_inter) {
      printf(" -> interpolated reionization history gives optical depth = %f\n", tau_reionization_);
    }
    if (pth->thermodynamics_verbose > 1) {
      printf(" -> free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
             tau_free_streaming_);
    }
    if ((all_species_.count("IDM_DR_IDR") > 0) && (pth->thermodynamics_verbose > 1)) {
      printf(" -> dark free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
             tau_idm_dr_fs);
    }
  }
}

/**
 * Assign value to each relevant index in vectors of thermodynamical quantities,
 * as well as in vector containing reionization parameters.
 *
 *
 * @param preco Input/Output: pointer to recombination structure
 * @param preio Input/Output: pointer to reionization structure
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_indices(recombination* preco, reionization* preio) {
  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - initialization of all indices and flags in thermo structure */
  index = 0;

  index_th_xe_ = index;
  index++;
  index_th_dkappa_ = index;
  index++;
  index_th_tau_d_ = index;
  index++;
  index_th_ddkappa_ = index;
  index++;
  index_th_dddkappa_ = index;
  index++;
  index_th_exp_m_kappa_ = index;
  index++;
  index_th_g_ = index;
  index++;
  index_th_dg_ = index;
  index++;
  index_th_ddg_ = index;
  index++;
  index_th_Tb_ = index;
  index++;
  index_th_wb_ = index;
  index++;
  index_th_cb2_ = index;
  index++;

  if (all_species_.count("IDM_DR_IDR") > 0) {
    index_th_dmu_idm_dr_ = index;
    index++;
    index_th_ddmu_idm_dr_ = index;
    index++;
    index_th_dddmu_idm_dr_ = index;
    index++;
    index_th_tau_idm_dr_ = index;
    index++;
    index_th_tau_idr_ = index;
    index++;
    index_th_g_idm_dr_ = index;
    index++;
    index_th_cidm_dr2_ = index;
    index++;
    index_th_Tidm_dr_ = index;
    index++;
    index_th_dmu_idr_ = index;
    index++;
  }

  /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives) {
    index_th_dcb2_ = index;
    index++;
    index_th_ddcb2_ = index;
    index++;
  }

  index_th_rate_ = index;
  index++;

  if (pth->compute_damping_scale) {
    index_th_r_d_ = index;
    index++;
  }

  /* end of indices */
  th_size_ = index;

  /** - initialization of all indices and flags in recombination structure */
  index = 0;

  preco->index_re_z = index;
  index++;
  preco->index_re_xe = index;
  index++;
  preco->index_re_dkappadtau = index;
  index++;
  preco->index_re_Tb = index;
  index++;
  preco->index_re_wb = index;
  index++;
  preco->index_re_cb2 = index;
  index++;

  /* end of indices */
  preco->re_size = index;

  /** - initialization of all indices and flags in reionization structure */
  index = 0;

  preio->index_re_z = index;
  index++;
  preio->index_re_xe = index;
  index++;
  preio->index_re_Tb = index;
  index++;
  preio->index_re_wb = index;
  index++;
  preio->index_re_cb2 = index;
  index++;
  preio->index_re_dkappadtau = index;
  index++;
  preio->index_re_dkappadz = index;
  index++;
  preio->index_re_d3kappadz3 = index;
  index++;

  /* end of indices */
  preio->re_size = index;

  /** - same with parameters of the function \f$ X_e(z)\f$ */

  index = 0;

  preio->index_reio_start = index;
  index++;

  /* case where x_e(z) taken like in CAMB (other cases can be added) */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    preio->index_reio_redshift = index;
    index++;
    preio->index_reio_exponent = index;
    index++;
    preio->index_reio_width = index;
    index++;
    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;
  }

  /* case where x_e(z) is binned */
  if (pth->reio_parametrization == reio_bins_tanh) {
    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the bin center array */

    preio->reio_num_z = pth->binned_reio_num + 2; /* add two values: beginning and end of reio */

    preio->index_reio_first_z         = index;
    index                            += preio->reio_num_z;
    preio->index_reio_first_xe        = index;
    index                            += preio->reio_num_z;
    preio->index_reio_step_sharpness  = index;
    index++;
  }

  /* case where x_e(z) has many tanh jumps */
  if (pth->reio_parametrization == reio_many_tanh) {
    /* the code will not only copy here the "jump centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the jump center array */

    preio->reio_num_z = pth->many_tanh_num + 2; /* add two values: beginning and end of reio */

    preio->index_reio_first_z         = index;
    index                            += preio->reio_num_z;
    preio->index_reio_first_xe        = index;
    index                            += preio->reio_num_z;
    preio->index_reio_step_sharpness  = index;
    index++;
  }

  /* case where x_e(z) must be interpolated */
  if (pth->reio_parametrization == reio_inter) {
    preio->reio_num_z = pth->reio_inter_num;

    preio->index_reio_first_z   = index;
    index                      += preio->reio_num_z;
    preio->index_reio_first_xe  = index;
    index                      += preio->reio_num_z;
  }

  preio->reio_num_params = index;

  /* flags for calling the interpolation routine */

  inter_normal_  = 0;
  inter_closeby_ = 1;
}

/**
 * Infer the primordial helium fraction from standard BBN, as a
 * function of the baryon density and expansion rate during BBN.
 *
 * This module is simpler then the one used in arXiv:0712.2826 because
 * it neglects the impact of a possible significant chemical
 * potentials for electron neutrinos. The full code with xi_nu_e could
 * be introduced here later.
 *
 * @return the error status
 */
void ThermodynamicsModule::thermodynamics_helium_from_bbn() {
  std::string line;
  const char* left;

  int num_omegab = 0;
  int num_deltaN = 0;

  std::vector<double> omegab;
  std::vector<double> deltaN;
  std::vector<double> YHe;
  std::vector<double> ddYHe;
  std::vector<double> YHe_at_deltaN;
  std::vector<double> ddYHe_at_deltaN;

  int array_line = 0;
  double DeltaNeff;
  double omega_b;
  int last_index;
  std::vector<double> pvecback;

  /**Summary: */
  /** - Infer effective number of neutrinos at the time of BBN */
  pvecback.resize(background_module_->bg_size_);

  /** - 8.6173e-11 converts from Kelvin to MeV. We randomly choose 0.1 MeV to be the temperature of BBN */
  double z_bbn = 0.1 / (8.6173e-11 * pba->T_cmb) - 1.0;

  double tau_bbn;
  background_module_->background_tau_of_z(z_bbn, &tau_bbn);

  background_module_->background_at_tau(tau_bbn,
                                        pba->long_info,
                                        pba->inter_normal,
                                        &last_index,
                                        pvecback.data());

  const double rho_g = all_species_.photons().Rho(pvecback.data());
  double Neff_bbn    = (pvecback[background_module_->index_bg_Omega_r_] *
                            pvecback[background_module_->index_bg_rho_crit_] -
                        rho_g) /
                       (7. / 8. * pow(4. / 11., 4. / 3.) * rho_g);

  /**DRMD**/
  if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
    auto& drmd  = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
    Neff_bbn   -= drmd.idr_drmd().Rho(pvecback.data()) / (7. / 8. * pow(4. / 11., 4. / 3.) * rho_g);
  }

  //  printf("Neff early = %g, Neff at bbn: %g\n", background_module_->Neff_, Neff_bbn);

  /** - compute Delta N_eff as defined in bbn file, i.e. \f$ \Delta N_{eff}=0\f$ means \f$ N_{eff}=3.046\f$ */
  DeltaNeff = Neff_bbn - 3.046;

  /* the following file is assumed to contain (apart from comments and blank lines):
     - the two numbers (num_omegab, num_deltaN) = number of values of BBN free parameters
     - three columns (omegab, deltaN, YHe) where omegab = Omega0_b h^2 and deltaN = Neff-3.046 by definition
     - omegab and deltaN are assumed to be arranged as:
     omegab1 deltaN1 YHe
     omegab2 deltaN1 YHe
     .....
     omegab1 delatN2 YHe
     omegab2 deltaN2 YHe
     .....
  */

  std::ifstream bbn_file(ppr->sBBN_file);
  class_test(!bbn_file.is_open(), "could not open BBN file with name %s", ppr->sBBN_file.c_str());

  /* go through each line */
  while (std::getline(bbn_file, line)) {
    /* eliminate blank spaces at beginning of line */
    left = line.c_str();
    while (left[0] == ' ') {
      left++;
    }

    /* check that the line is neither blank neither a comment. In
       ASCII, left[0]>39 means that first non-blank character might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */
    if (left[0] > 39) {
      /* if the line contains data, we must interpret it. If
         (num_omegab, num_deltaN)=(0,0), the current line must contain
         their values. Otherwise, it must contain (omegab, delatN,
         YHe). */
      if ((num_omegab == 0) && (num_deltaN == 0)) {
        /* read (num_omegab, num_deltaN), infer size of arrays and allocate them */
        class_test(sscanf(line.c_str(), "%d %d", &num_omegab, &num_deltaN) != 2,
                   "could not read value of parameters (num_omegab,num_deltaN) in file %s\n",
                   ppr->sBBN_file.c_str());

        class_test(num_omegab <= 0 || num_deltaN <= 0,
                   "read num_omegab=%d, num_deltaN=%d in file %s, expected positive values\n",
                   num_omegab,
                   num_deltaN,
                   ppr->sBBN_file.c_str());

        omegab.resize(num_omegab);
        deltaN.resize(num_deltaN);
        YHe.resize(num_omegab * num_deltaN);
        ddYHe.resize(num_omegab * num_deltaN);
        YHe_at_deltaN.resize(num_omegab);
        ddYHe_at_deltaN.resize(num_omegab);
        array_line = 0;
      }
      else {
        /* read (omegab, deltaN, YHe) */
        class_test(sscanf(line.c_str(),
                          "%lg %lg %lg",
                          &(omegab[array_line % num_omegab]),
                          &(deltaN[array_line / num_omegab]),
                          &(YHe[array_line])) != 3,
                   "could not read value of parameters (omegab,deltaN,YHe) in file %s\n",
                   ppr->sBBN_file.c_str());
        array_line++;
      }
    }
  }

  /** - spline in one dimension (along deltaN) */
  array_spline_table_lines(deltaN.data(),
                           num_deltaN,
                           YHe.data(),
                           num_omegab,
                           ddYHe.data(),
                           _SPLINE_NATURAL_);

  omega_b = pba->Omega0_b * pba->h * pba->h;

  class_test(omega_b < omegab[0],

             "You have asked for an unrealistic small value omega_b = %e. The corresponding value "
             "of the primordial helium fraction cannot be found in the interpolation table. If you "
             "really want this value, you should fix YHe to a given value rather than to BBN",
             omega_b);

  class_test(omega_b > omegab[num_omegab - 1],

             "You have asked for an unrealistic high value omega_b = %e. The corresponding value "
             "of the primordial helium fraction cannot be found in the interpolation table. If you "
             "really want this value, you should fix YHe to a given value rather than to BBN",
             omega_b);

  class_test(DeltaNeff < deltaN[0],

             "You have asked for an unrealistic small value of Delta N_eff = %e. The corresponding "
             "value of the primordial helium fraction cannot be found in the interpolation table. "
             "If you really want this value, you should fix YHe to a given value rather than to "
             "BBN",
             DeltaNeff);

  class_test(DeltaNeff > deltaN[num_deltaN - 1],

             "You have asked for an unrealistic high value of Delta N_eff = %e. The corresponding "
             "value of the primordial helium fraction cannot be found in the interpolation table. "
             "If you really want this value, you should fix YHe to a given value rather than to "
             "BBN",
             DeltaNeff);

  /** - interpolate in one dimension (along deltaN) */
  array_interpolate_spline(deltaN.data(),
                           num_deltaN,
                           YHe.data(),
                           ddYHe.data(),
                           num_omegab,
                           DeltaNeff,
                           &last_index,
                           YHe_at_deltaN.data(),
                           num_omegab);

  /** - spline in remaining dimension (along omegab) */
  array_spline_table_lines(omegab.data(),
                           num_omegab,
                           YHe_at_deltaN.data(),
                           1,
                           ddYHe_at_deltaN.data(),
                           _SPLINE_NATURAL_);

  /** - interpolate in remaining dimension (along omegab) */
  array_interpolate_spline(omegab.data(),
                           num_omegab,
                           YHe_at_deltaN.data(),
                           ddYHe_at_deltaN.data(),
                           1,
                           omega_b,
                           &last_index,
                           &(YHe_),
                           1);

  /** - deallocate arrays */
}

/**
 * In case of non-minimal cosmology, this function determines the
 * energy rate injected in the IGM at a given redshift z (= on-the-spot
 * annihilation). This energy injection may come e.g. from dark matter
 * annihilation or decay.
 *
 * @param preco Input: pointer to recombination structure
 * @param z Input: redshift
 * @param energy_rate Output: energy density injection rate
 * @param error_message Output: error message
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_onthespot_energy_injection(recombination* preco,
                                                                     double z,
                                                                     double* energy_rate) {
  /*redshift-dependent annihilation parameter*/

  double annihilation_at_z;
  if (z > preco->annihilation_zmax) {
    annihilation_at_z = preco->annihilation *
                        exp(-preco->annihilation_variation *
                            pow(log((preco->annihilation_z + 1.) / (preco->annihilation_zmax + 1.)),
                                2));
  }
  else if (z > preco->annihilation_zmin) {
    annihilation_at_z =
        preco->annihilation *
        exp(preco->annihilation_variation *
            (-pow(log((preco->annihilation_z + 1.) / (preco->annihilation_zmax + 1.)), 2) +
             pow(log((z + 1.) / (preco->annihilation_zmax + 1.)), 2)));
  }
  else {
    annihilation_at_z =
        preco->annihilation *
        exp(preco->annihilation_variation *
            (-pow(log((preco->annihilation_z + 1.) / (preco->annihilation_zmax + 1.)), 2) +
             pow(log((preco->annihilation_zmin + 1.) / (preco->annihilation_zmax + 1.)), 2)));
  }

  const double Omega0_idm_dr_ann = all_species_.count("IDM_DR_IDR")
                                       ? static_cast<const IDM_DR_IDR_Species&>(
                                             *all_species_.at("IDM_DR_IDR"))
                                             .idm_dr()
                                             .GetOmega0()
                                       : 0.;
  const double Omega0_cdm_ann    = all_species_.count("CDM") ? all_species_.at("CDM")->GetOmega0()
                                                             : 0.;
  double rho_cdm_today           = pow(pba->H0 * _c_ / _Mpc_over_m_, 2) * 3 / 8. / _PI_ / _G_ *
                                   (Omega0_idm_dr_ann + Omega0_cdm_ann) * _c_ *
                                   _c_; /* energy density in J/m^3 */

  double u_min = (1 + z) / (1 + preco->annihilation_z_halo);

  double erfc = pow(1. + 0.278393 * u_min + 0.230389 * u_min * u_min +
                        0.000972 * u_min * u_min * u_min + 0.078108 * u_min * u_min * u_min * u_min,
                    -4);

  *energy_rate = pow(rho_cdm_today, 2) / _c_ / _c_ * pow((1 + z), 3) *
                     (pow((1. + z), 3) * annihilation_at_z + preco->annihilation_f_halo * erfc) +
                 rho_cdm_today * pow((1 + z), 3) * preco->decay;
  /* energy density rate in J/m^3/s (remember that annihilation_at_z is in m^3/s/Kg and decay in s^-1) */
}

/**
 * In case of non-minimal cosmology, this function determines the
 * effective energy rate absorbed by the IGM at a given redshift
 * (beyond the on-the-spot annihilation). This energy injection may
 * come e.g. from dark matter annihilation or decay.
 *
 * @param preco Input: pointer to recombination structure
 * @param z Input: redshift
 * @param energy_rate Output: energy density injection rate
 * @param error_message Output: error message
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_energy_injection(recombination* preco,
                                                           double z,
                                                           double* energy_rate) {
  if (preco->annihilation > 0) {
    double result;

    if (!preco->has_on_the_spot) {
      /* number of hydrogen nuclei today in m**-3 */
      double nH0 = 3. * preco->H0 * preco->H0 * pba->Omega0_b / (8. * _PI_ * _G_ * _m_H_) *
                   (1. - preco->YHe);

      const double Omega0_cdm_inj = all_species_.count("CDM") ? all_species_.at("CDM")->GetOmega0()
                                                              : 0.;
      const double Omega0_idm_dr_inj = all_species_.count("IDM_DR_IDR")
                                           ? static_cast<const IDM_DR_IDR_Species&>(
                                                 *all_species_.at("IDM_DR_IDR"))
                                                 .idm_dr()
                                                 .GetOmega0()
                                           : 0.;
      /* factor = c sigma_T n_H(0) / (H(0) \sqrt(Omega_m)) (dimensionless) */
      double factor = _sigma_ * nH0 / pba->H0 * _Mpc_over_m_ /
                      sqrt(pba->Omega0_b + Omega0_cdm_inj + Omega0_idm_dr_inj);

      /* integral over z'(=zp) with step dz */
      double dz = 1.;

      /* first point in trapezoidal integral */
      double zp = z;
      double onthespot;
      thermodynamics_onthespot_energy_injection(preco, zp, &onthespot);
      double first_integrand =
          factor * pow(1 + z, 8) / pow(1 + zp, 7.5) *
          exp(2. / 3. * factor * (pow(1 + z, 1.5) - pow(1 + zp, 1.5))) *
          onthespot;  // beware: versions before 2.4.3, there were wrong exponents: 6 and 5.5 instead of 8 and 7.5
      result = 0.5 * dz * first_integrand;

      /* other points in trapezoidal integral */
      double integrand;
      do {
        zp += dz;
        thermodynamics_onthespot_energy_injection(preco, zp, &onthespot);
        integrand =
            factor * pow(1 + z, 8) / pow(1 + zp, 7.5) *
            exp(2. / 3. * factor * (pow(1 + z, 1.5) - pow(1 + zp, 1.5))) *
            onthespot;  // beware: versions before 2.4.3, there were wrong exponents: 6 and 5.5 instead of 8 and 7.5
        result += dz * integrand;

      } while (integrand / first_integrand > 0.02);

      /* uncomment these lines if you also want to compute the on-the-spot for comparison */
      thermodynamics_onthespot_energy_injection(preco, z, &onthespot);
    }
    else {
      thermodynamics_onthespot_energy_injection(preco, z, &result);
    }

    /* these test lines print the energy rate rescaled by (1+z)^6 in J/m^3/s, with or without the on-the-spot approximation */
    /*
      fprintf(stdout,"%e  %e  %e \n",
      1.+z,
      result/pow(1.+z,6),
      onthespot/pow(1.+z,6));
    */

    /* effective energy density rate in J/m^3/s  */
    *energy_rate = result;
  }
  else {
    *energy_rate = 0.;
  }
}

/**
 * This subroutine contains the reionization function \f$ X_e(z) \f$
 * (one for each scheme; so far, only the function corresponding to
 * the reio_camb scheme is coded)
 *
 * @param z     Input: redshift
 * @param preio Input: pointer to reionization structure, containing the parameters of the function \f$ X_e(z) \f$
 * @param xe    Output: \f$ X_e(z) \f$
 */

void ThermodynamicsModule::thermodynamics_reionization_function(double z,
                                                                reionization* preio,
                                                                double* xe) {
  /** Summary: */

  /** - define local variables */
  double argument;

  /** - implementation of ionization function similar to the one in CAMB */

  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];
    }

    else {
      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */

      argument =
          (pow((1. + preio->reionization_parameters[preio->index_reio_redshift]),
               preio->reionization_parameters[preio->index_reio_exponent]) -
           pow((1. + z), preio->reionization_parameters[preio->index_reio_exponent])) /
          (preio->reionization_parameters[preio->index_reio_exponent]
           /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
           * pow((1. + preio->reionization_parameters[preio->index_reio_redshift]),
                 (preio->reionization_parameters[preio->index_reio_exponent] - 1.))) /
          preio->reionization_parameters[preio->index_reio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */

      if (pth->reio_parametrization == reio_camb) {
        *xe = (preio->reionization_parameters[preio->index_reio_xe_after] -
               preio->reionization_parameters[preio->index_reio_xe_before]) *
                  (tanh(argument) + 1.) / 2. +
              preio->reionization_parameters[preio->index_reio_xe_before];
      }
      else {
        *xe = (preio->reionization_parameters[preio->index_reio_xe_after] -
               preio->reionization_parameters[preio->index_reio_xe_before]) *
                  tanh(argument) +
              preio->reionization_parameters[preio->index_reio_xe_before];
      }

      /** - --> case z < z_reio_start: helium contribution (tanh of simpler argument) */

      if (pth->reio_parametrization == reio_camb) {
        argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z) /
                   preio->reionization_parameters[preio->index_helium_fullreio_width];
        /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
        *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction] *
               (tanh(argument) + 1.) / 2.;
      }
    }

    return;
  }

  /** - implementation of binned ionization function similar to astro-ph/0606552 */

  if (pth->reio_parametrization == reio_bins_tanh) {
    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1];
    }

    else if (z < preio->reionization_parameters[preio->index_reio_first_z]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    else {
      int i = 0;
      while (preio->reionization_parameters[preio->index_reio_first_z + i + 1] < z)
        i++;

      /* This is the expression of the tanh-like jumps of the
         reio_bins_tanh scheme until the 10.06.2015. It appeared to be
         not robust enough. It could lead to a kink in xe(z) near the
         maximum value of z at which reionisation is sampled. It has
         been replaced by the simpler and more robust expression
         below.

         *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
         +0.5*(tanh((2.*(z-preio->reionization_parameters[preio->index_reio_first_z+i])
         /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
         -preio->reionization_parameters[preio->index_reio_first_z+i])-1.)
         /preio->reionization_parameters[preio->index_reio_step_sharpness])
         /tanh(1./preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
         *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
         -preio->reionization_parameters[preio->index_reio_first_xe+i]);
      */

      /* compute the central redshift value of the tanh jump */

      double z_jump;
      if (i == preio->reio_num_z - 2) {
        z_jump = preio->reionization_parameters[preio->index_reio_first_z + i] +
                 0.5 * (preio->reionization_parameters[preio->index_reio_first_z + i] -
                        preio->reionization_parameters[preio->index_reio_first_z + i - 1]);
      }
      else {
        z_jump = 0.5 * (preio->reionization_parameters[preio->index_reio_first_z + i + 1] +
                        preio->reionization_parameters[preio->index_reio_first_z + i]);
      }

      /* implementation of the tanh jump */

      *xe = preio->reionization_parameters[preio->index_reio_first_xe + i] +
            0.5 *
                (tanh((z - z_jump) /
                      preio->reionization_parameters[preio->index_reio_step_sharpness]) +
                 1.) *
                (preio->reionization_parameters[preio->index_reio_first_xe + i + 1] -
                 preio->reionization_parameters[preio->index_reio_first_xe + i]);
    }

    return;
  }

  /** - implementation of many tanh jumps */

  if (pth->reio_parametrization == reio_many_tanh) {
    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1];
    }

    else if (z > preio->reionization_parameters[preio->index_reio_first_z]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1];

      for (int jump = 1; jump < preio->reio_num_z - 1; jump++) {
        double center = preio->reionization_parameters[preio->index_reio_first_z +
                                                       preio->reio_num_z - 1 - jump];
        // before and after are meant with respect to growing z, not growing time
        double before =
            preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1 -
                                           jump] -
            preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - jump];
        double after = 0.;
        double width = preio->reionization_parameters[preio->index_reio_step_sharpness];

        double one_jump;
        thermodynamics_tanh(z, center, before, after, width, &one_jump);

        *xe += one_jump;
      }
    }

    else {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    return;
  }

  /** - implementation of reio_inter */

  if (pth->reio_parametrization == reio_inter) {
    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1];
      class_stop("Check: is it normal that we are here?");
    }

    else {
      int i = 0;
      while (preio->reionization_parameters[preio->index_reio_first_z + i + 1] < z)
        i++;

      double z_min = preio->reionization_parameters[preio->index_reio_first_z + i];
      double z_max = preio->reionization_parameters[preio->index_reio_first_z + i + 1];

      class_test(z < z_min, "");

      class_test(z > z_max, "");

      double x = (z - preio->reionization_parameters[preio->index_reio_first_z + i]) /
                 (preio->reionization_parameters[preio->index_reio_first_z + i + 1] -
                  preio->reionization_parameters[preio->index_reio_first_z + i]);

      *xe = preio->reionization_parameters[preio->index_reio_first_xe + i] +
            x * (preio->reionization_parameters[preio->index_reio_first_xe + i + 1] -
                 preio->reionization_parameters[preio->index_reio_first_xe + i]);

      class_test(*xe < 0.,
                 "%e %e %e\n",
                 x,
                 preio->reionization_parameters[preio->index_reio_first_xe + i],
                 preio->reionization_parameters[preio->index_reio_first_xe + i + 1]);
    }

    return;
  }

  class_test(0 == 0, "value of reio_parametrization=%d unclear", pth->reio_parametrization);
}

/**
 * This subroutine reads \f$ X_e(z) \f$ in the recombination table at
 * the time at which reionization starts. Hence it provides correct
 * initial conditions for the reionization function.
 *
 * @param preco Input: pointer to recombination structure
 * @param z     Input: redshift z_reio_start
 * @param xe    Output: \f$ X_e(z) \f$ at z
 */

void ThermodynamicsModule::thermodynamics_get_xe_before_reionization(recombination* preco,
                                                                     double z,
                                                                     double* xe) {
  int last_index = 0;

  array_interpolate_one_growing_closeby(preco->recombination_table.data(),
                                        preco->re_size,
                                        preco->rt_size,
                                        preco->index_re_z,
                                        z,
                                        &last_index,
                                        preco->index_re_xe,
                                        xe);
}

/**
 * This routine computes the reionization history. In the reio_camb
 * scheme, this is straightforward if the input parameter is the
 * reionization redshift. If the input is the optical depth, need to
 * find z_reio by dichotomy (trying several z_reio until the correct
 * tau_reio is approached).
 *
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantities (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_reionization(recombination* preco,
                                                       reionization* preio,
                                                       double* pvecback) {
  /** Summary: */

  /** - define local variables */

  double z_sup, z_inf;
  double tau_sup, tau_inf;

  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */

  preio->reionization_parameters.resize(preio->reio_num_params);

  /** - (a) if reionization implemented like in CAMB */

  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    /** - --> set values of these parameters, excepted those depending on the reionization redshift */

    if (pth->reio_parametrization == reio_camb) {
      preio->reionization_parameters[preio->index_reio_xe_after] =
          1. +
          YHe_ /
              (_not4_ *
               (1. -
                YHe_)); /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    }
    if (pth->reio_parametrization == reio_half_tanh) {
      preio->reionization_parameters[preio->index_reio_xe_after] =
          1.; /* xe_after_reio: neglect He ionization */
      //+ 2*YHe_/(_not4_*(1.-YHe_));    /* xe_after_reio: H + fully ionized He */
    }
    preio->reionization_parameters[preio->index_reio_exponent] =
        pth->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_reio_width] =
        pth->reionization_width; /* reio_width */
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] =
        YHe_ /
        (_not4_ *
         (1. -
          YHe_)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] =
        pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] =
        pth->helium_fullreio_width; /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_reio_exponent] == 0,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_reio_width] == 0,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width] == 0,
               "stop to avoid division by zero");

    /** - --> if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {
      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_reionization_;

      /* infer starting redshift for hydrogen */

      if (pth->reio_parametrization == reio_camb) {
        preio->reionization_parameters[preio->index_reio_start] =
            preio->reionization_parameters[preio->index_reio_redshift] +
            ppr->reionization_start_factor * pth->reionization_width;

        /* if starting redshift for helium is larger, take that one
           (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_reio_start] <
            pth->helium_fullreio_redshift +
                ppr->reionization_start_factor * pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_reio_start] = pth->helium_fullreio_redshift +
                                                                    ppr->reionization_start_factor *
                                                                        pth->helium_fullreio_width;
      }
      else {
        preio->reionization_parameters[preio->index_reio_start] = z_reionization_;
      }

      class_test(preio->reionization_parameters[preio->index_reio_start] >
                     ppr->reionization_z_start_max,
                 "starting redshift for reionization > reionization_z_start_max = %e\n",
                 ppr->reionization_z_start_max);

      /* infer xe_before_reio */
      thermodynamics_get_xe_before_reionization(
          preco,
          preio->reionization_parameters[preio->index_reio_start],
          &(preio->reionization_parameters[preio->index_reio_xe_before]));

      /* fill reionization table */
      thermodynamics_reionization_sample(preco, preio, pvecback);

      tau_reionization_ = preio->reionization_optical_depth;
    }

    /** - --> if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {
      /* upper value */

      z_sup = ppr->reionization_z_start_max -
              ppr->reionization_start_factor * pth->reionization_width;
      class_test(z_sup < 0.,
                 "parameters are such that reionization cannot take place before today while "
                 "starting after z_start_max; need to increase z_start_max");

      /* maximum possible reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup;
      /* maximum possible starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */
      thermodynamics_get_xe_before_reionization(
          preco,
          preio->reionization_parameters[preio->index_reio_start],
          &(preio->reionization_parameters[preio->index_reio_xe_before]));

      /* fill reionization table */
      thermodynamics_reionization_sample(preco, preio, pvecback);

      tau_sup = preio->reionization_optical_depth;

      if (tau_sup < tau_reionization_) {
        class_stop("parameters are such that reionization cannot start after z_start_max");
      }

      /* lower value */

      z_inf   = 0.;
      tau_inf = 0.;

      /* try intermediate values */

      int counter = 0;
      /* bisection kept inline: bisects z but stops on (tau_sup - tau_inf) and the reionization table from the final iteration is read afterward */
      while ((tau_sup - tau_inf) > tau_reionization_ * ppr->reionization_optical_depth_tol) {
        double z_mid = 0.5 * (z_sup + z_inf);

        /* reionization redshift */
        preio->reionization_parameters[preio->index_reio_redshift] = z_mid;
        /* infer starting redshift for hygrogen */
        preio->reionization_parameters[preio->index_reio_start] =
            preio->reionization_parameters[preio->index_reio_redshift] +
            ppr->reionization_start_factor * pth->reionization_width;
        /* if starting redshift for helium is larger, take that one
           (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_reio_start] <
            pth->helium_fullreio_redshift +
                ppr->reionization_start_factor * pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_reio_start] = pth->helium_fullreio_redshift +
                                                                    ppr->reionization_start_factor *
                                                                        pth->helium_fullreio_width;

        class_test(preio->reionization_parameters[preio->index_reio_start] >
                       ppr->reionization_z_start_max,
                   "starting redshift for reionization > reionization_z_start_max = %e",
                   ppr->reionization_z_start_max);

        /* infer xe_before_reio */
        thermodynamics_get_xe_before_reionization(
            preco,
            preio->reionization_parameters[preio->index_reio_start],
            &(preio->reionization_parameters[preio->index_reio_xe_before]));

        /* clean and fill reionization table */

        thermodynamics_reionization_sample(preco, preio, pvecback);

        double tau_mid = preio->reionization_optical_depth;

        /* trial */

        if (tau_mid > tau_reionization_) {
          z_sup   = z_mid;
          tau_sup = tau_mid;
        }
        else {
          z_inf   = z_mid;
          tau_inf = tau_mid;
        }

        counter++;
        class_test(counter > _MAX_IT_,
                   "while searching for reionization_optical_depth, maximum number of iterations "
                   "exceeded");
      }

      /* store z_reionization in thermodynamics structure */
      z_reionization_ = preio->reionization_parameters[preio->index_reio_redshift];
    }

    return;
  }

  /** - (b) if reionization implemented with reio_bins_tanh scheme */

  if (pth->reio_parametrization == reio_bins_tanh) {
    /* this algorithm requires at least two bin centers (i.e. at least
       4 values in the (z,xe) array, counting the edges). */
    class_test(pth->binned_reio_num < 2,
               "current implementation of binned reio requires at least two bin centers");

    /* check that this input can be interpreted by the code */
    for (int bin = 1; bin < pth->binned_reio_num; bin++) {
      class_test(pth->binned_reio_z[bin - 1] >= pth->binned_reio_z[bin],
                 "value of reionization bin centers z_i expected to be passed in growing order: "
                 "%e, %e",
                 pth->binned_reio_z[bin - 1],
                 pth->binned_reio_z[bin]);
    }

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */

    for (int bin = 1; bin < preio->reio_num_z - 1; bin++) {
      preio->reionization_parameters[preio->index_reio_first_z + bin] = pth->binned_reio_z[bin - 1];
      preio->reionization_parameters[preio->index_reio_first_xe + bin] =
          pth->binned_reio_xe[bin - 1];
    }

    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + 2*(the distance between z_(i_max) and z_(i_max-1)). E.g. if
       the bins are in 10,12,14, the largest z will be 18. */
    preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1] =
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 2] +
        2. * (preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 2] -
              preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 3]);

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] =
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] >
                   ppr->reionization_z_start_max,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must "
               "change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose
       to define it as z_0 - (the distance between z_1 and z_0). E.g. if
       the bins are in 10,12,14, the stop redshift will be 8. */

    preio->reionization_parameters[preio->index_reio_first_z] =
        2. * preio->reionization_parameters[preio->index_reio_first_z + 1] -
        preio->reionization_parameters[preio->index_reio_first_z + 2];

    /* check it's not too small */
    /* 6.06.2015: changed this test to simply imposing that the first z is at least zero */
    /*
      class_test(preio->reionization_parameters[preio->index_reio_first_z] < 0,
      "final redshift for reionization = %e, you must change the binning or redefine the way in which the code extrapolates below the first value of z_i",preio->reionization_parameters[preio->index_reio_first_z]);
    */
    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe before reio */
    thermodynamics_get_xe_before_reionization(
        preco,
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1],
        &(preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1]));

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_reio_first_xe] =
        1. +
        YHe_ /
            (_not4_ *
             (1. -
              YHe_)); /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_reio_step_sharpness] =
        pth->binned_reio_step_sharpness;

    /* fill reionization table */
    thermodynamics_reionization_sample(preco, preio, pvecback);

    tau_reionization_ = preio->reionization_optical_depth;

    return;
  }

  /** - (c) if reionization implemented with reio_many_tanh scheme */

  if (pth->reio_parametrization == reio_many_tanh) {
    /* this algorithm requires at least one jump centers */
    class_test(pth->many_tanh_num < 1,
               "current implementation of reio_many_tanh requires at least one jump center");

    /* check that z input can be interpreted by the code */
    for (int bin = 1; bin < pth->many_tanh_num; bin++) {
      class_test(pth->many_tanh_z[bin - 1] >= pth->many_tanh_z[bin],
                 "value of reionization bin centers z_i expected to be passed in growing order: "
                 "%e, %e",
                 pth->many_tanh_z[bin - 1],
                 pth->many_tanh_z[bin]);
    }

    /* the code will not only copy here the "jump centers" passed in
       input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */

    for (int bin = 1; bin < preio->reio_num_z - 1; bin++) {
      preio->reionization_parameters[preio->index_reio_first_z + bin] = pth->many_tanh_z[bin - 1];

      /* check that xe input can be interpreted by the code */
      double xe_input = pth->many_tanh_xe[bin - 1];
      double xe_actual;
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input < -0.9) && (xe_input > -1.1)) {
        xe_actual = 1. + YHe_ / (_not4_ * (1. - YHe_));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input < -1.9) && (xe_input > -2.1)) {
        xe_actual = 1. + 2. * YHe_ / (_not4_ * (1. - YHe_));
      }
      //other negative number is nonsense
      else {
        class_stop(
            "Your entry for many_tanh_xe[%d] is %e, this makes no sense (either positive or "
            "0,-1,-2)",
            bin - 1,
            pth->many_tanh_xe[bin - 1]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe + bin] = xe_actual;
    }

    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + ppr->reionization_start_factor*step_sharpness. */
    preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1] =
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 2] +
        ppr->reionization_start_factor * pth->many_tanh_width;

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] =
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] >
                   ppr->reionization_z_start_max,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must "
               "change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose
       to define it as z_0 - ppr->reionization_start_factor*step_sharpness, but at least zero. */

    preio->reionization_parameters[preio->index_reio_first_z] =
        preio->reionization_parameters[preio->index_reio_first_z + 1] -
        ppr->reionization_start_factor * pth->many_tanh_width;

    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe before reio */
    thermodynamics_get_xe_before_reionization(
        preco,
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1],
        &(preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1]));

    /* infer xe after reio */

    preio->reionization_parameters[preio->index_reio_first_xe] =
        preio->reionization_parameters[preio->index_reio_first_xe + 1];

    /* if we want to model only hydrogen reionization and neglect both helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1.;

    /* if we want to model only hydrogen + first helium reionization and neglect second helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + YHe_/(_not4_*(1.-YHe_));

    /* if we want to model hydrogen + two helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + 2.*YHe_/(_not4_*(1.-YHe_));

    /* pass step sharpness parameter */
    class_test(pth->many_tanh_width <= 0,
               "many_tanh_width must be strictly positive, you passed %e",
               pth->many_tanh_width);

    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->many_tanh_width;

    /* fill reionization table */
    thermodynamics_reionization_sample(preco, preio, pvecback);

    tau_reionization_ = preio->reionization_optical_depth;

    return;
  }

  /** - (d) if reionization implemented with reio_inter scheme */

  if (pth->reio_parametrization == reio_inter) {
    /* this parametrization requires at least one point (z,xe) */
    class_test(pth->reio_inter_num < 1,
               "current implementation of reio_inter requires at least one point (z,xe)");

    /* this parametrization requires that the first z value is zero */
    class_test(pth->reio_inter_z[0] != 0.,
               "For reio_inter scheme, the first value of reio_inter_z[...]  should always be "
               "zero, you passed %e",
               pth->reio_inter_z[0]);

    /* check that z input can be interpreted by the code */
    for (int point = 1; point < pth->reio_inter_num; point++) {
      class_test(pth->reio_inter_z[point - 1] >= pth->reio_inter_z[point],
                 "value of reionization bin centers z_i expected to be passed in growing order, "
                 "unlike: %e, %e",
                 pth->reio_inter_z[point - 1],
                 pth->reio_inter_z[point]);
    }

    /* this parametrization requires that the last x_i value is zero
       (the code will substitute it with the value that one would get in
       absence of reionization, as compute by the recombination code) */
    class_test(pth->reio_inter_xe[pth->reio_inter_num - 1] != 0.,
               "For reio_inter scheme, the last value of reio_inter_xe[...]  should always be "
               "zero, you passed %e",
               pth->reio_inter_xe[pth->reio_inter_num - 1]);

    /* copy here the (z,xe) values passed in input. */

    for (int point = 0; point < preio->reio_num_z; point++) {
      preio->reionization_parameters[preio->index_reio_first_z + point] = pth->reio_inter_z[point];

      /* check that xe input can be interpreted by the code */
      double xe_input = pth->reio_inter_xe[point];
      double xe_actual;
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input < -0.9) && (xe_input > -1.1)) {
        xe_actual = 1. + YHe_ / (_not4_ * (1. - YHe_));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input < -1.9) && (xe_input > -2.1)) {
        xe_actual = 1. + 2. * YHe_ / (_not4_ * (1. - YHe_));
      }
      //other negative number is nonsense
      else {
        class_stop(
            "Your entry for reio_inter_xe[%d] is %e, this makes no sense (either positive "
            "or 0,-1,-2)",
            point,
            pth->reio_inter_xe[point]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe + point] = xe_actual;
    }

    /* copy highest redshift in reio_start */
    preio->reionization_parameters[preio->index_reio_start] =
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] >
                   ppr->reionization_z_start_max,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must "
               "change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* infer xe before reio */
    thermodynamics_get_xe_before_reionization(
        preco,
        preio->reionization_parameters[preio->index_reio_first_z + preio->reio_num_z - 1],
        &(preio->reionization_parameters[preio->index_reio_first_xe + preio->reio_num_z - 1]));

    /* fill reionization table */
    thermodynamics_reionization_sample(preco, preio, pvecback);

    tau_reionization_ = preio->reionization_optical_depth;

    return;
  }

  class_test(0 == 0, "value of reio_z_or_tau=%d unclear", pth->reio_z_or_tau);
}

/**
 * For fixed input reionization parameters, this routine computes the
 * reionization history and fills the reionization table.
 *
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantities (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_reionization_sample(recombination* preco,
                                                              reionization* preio,
                                                              double* pvecback) {
  /** Summary: */

  double Yp = YHe_;

  /** - (a) allocate vector of values related to reionization */
  std::vector<double> reio_vector;
  reio_vector.resize(preio->re_size);

  /** - (b) create a growing buffer holding the reionization table rows */
  std::vector<double> reio_buffer;

  /** - (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - --> look where to start in current thermodynamics table */
  int i = 0;
  while (preco->recombination_table[i * preco->re_size + preco->index_re_z] <
         preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    class_test(i == ppr->recfast_Nz0,
               "reionization_z_start_max = %e > largest redshift in thermodynamics table",
               ppr->reionization_z_start_max);
  }

  /** - --> get redshift */
  double z = preco->recombination_table[i * preco->re_size + preco->index_re_z];
  reio_vector[preio->index_re_z]    = z;
  preio->index_reco_when_reio_start = i;

  /** - --> get \f$ X_e \f$ */
  double xe;
  thermodynamics_reionization_function(z, preio, &xe);

  reio_vector[preio->index_re_xe] = xe;

  /** -  --> get \f$ d \kappa / d z = (d \kappa / d \tau) * (d \tau / d z) = - (d \kappa / d \tau) / H \f$ */

  double tau;
  background_module_->background_tau_of_z(z, &tau);

  int last_index_back;
  background_module_->background_at_tau(tau,
                                        pba->short_info,
                                        pba->inter_normal,
                                        &last_index_back,
                                        pvecback);

  reio_vector[preio->index_re_dkappadtau] = (1. + z) * (1. + z) * n_e_ * xe * _sigma_ *
                                            _Mpc_over_m_;

  class_test(pvecback[background_module_->index_bg_H_] == 0., "stop to avoid division by zero");

  reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadtau] /
                                          pvecback[background_module_->index_bg_H_];

  double dkappadz   = reio_vector[preio->index_re_dkappadz];
  double dkappadtau = reio_vector[preio->index_re_dkappadtau];

  /** - --> get baryon temperature **/
  double Tb = preco->recombination_table[i * preco->re_size + preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - --> after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  //Tba2 = Tb/(1+z)/(1+z);

  /** - --> get baryon equation of state */
  reio_vector[preio->index_re_wb] = _k_B_ / (_c_ * _c_ * _m_H_) *
                                    (1. + (1. / _not4_ - 1.) * Yp + xe * (1. - Yp)) * Tb;

  /** - --> get baryon adiabatic sound speed */
  reio_vector[preio->index_re_cb2] = 5. / 3. * reio_vector[preio->index_re_wb];

  /** - --> store these values in growing buffer */
  reio_buffer.insert(reio_buffer.end(), reio_vector.begin(), reio_vector.end());

  int number_of_redshifts = 1;

  /** - (d) set the maximum step value (equal to the step in thermodynamics table) */
  double dz_max = preco->recombination_table[i * preco->re_size + preco->index_re_z] -
                  preco->recombination_table[(i - 1) * preco->re_size + preco->index_re_z];

  /** - (e) loop over redshift values in order to find values of z, x_e, kappa' (Tb and cb2 found later by integration). The sampling in z space is found here. */

  /* initial step */
  double dz = dz_max;

  while (z > 0.) {
    class_test(dz < ppr->smallest_allowed_variation,
               "stuck in the loop for reionization sampling, as if you were trying to impose a "
               "discontinuous evolution for xe(z)");

    /* - try next step */
    double z_next = z - dz;
    if (z_next < 0.)
      z_next = 0.;

    double xe_next;
    thermodynamics_reionization_function(z_next, preio, &xe_next);

    background_module_->background_tau_of_z(z_next, &tau);

    background_module_->background_at_tau(tau,
                                          pba->short_info,
                                          pba->inter_normal,
                                          &last_index_back,
                                          pvecback);

    class_test(pvecback[background_module_->index_bg_H_] == 0., "stop to avoid division by zero");

    double dkappadz_next = (1. + z_next) * (1. + z_next) * n_e_ * xe_next * _sigma_ * _Mpc_over_m_ /
                           pvecback[background_module_->index_bg_H_];

    double dkappadtau_next = (1. + z_next) * (1. + z_next) * n_e_ * xe_next * _sigma_ *
                             _Mpc_over_m_;

    class_test((dkappadz == 0.) || (dkappadtau == 0.), "stop to avoid division by zero");

    double relative_variation = fabs((dkappadz_next - dkappadz) / dkappadz) +
                                fabs((dkappadtau_next - dkappadtau) / dkappadtau);

    if (relative_variation < ppr->reionization_sampling) {
      /* accept the step: get \f$ z, X_e, d kappa / d z \f$ and store in growing table */

      z          = z_next;
      xe         = xe_next;
      dkappadz   = dkappadz_next;
      dkappadtau = dkappadtau_next;

      class_test((dkappadz == 0.) || (dkappadtau == 0.),
                 "dkappadz=%e, dkappadtau=%e, stop to avoid division by zero",
                 dkappadz,
                 dkappadtau);

      reio_vector[preio->index_re_z]          = z;
      reio_vector[preio->index_re_xe]         = xe;
      reio_vector[preio->index_re_dkappadz]   = dkappadz;
      reio_vector[preio->index_re_dkappadtau] = dkappadz *
                                                pvecback[background_module_->index_bg_H_];

      reio_buffer.insert(reio_buffer.end(), reio_vector.begin(), reio_vector.end());

      number_of_redshifts++;

      dz = std::min(0.9 * (ppr->reionization_sampling / relative_variation), 5.) * dz;
      dz = std::min(dz, dz_max);
    }
    else {
      /* do not accept the step and update dz */
      dz = 0.9 * (ppr->reionization_sampling / relative_variation) * dz;
    }
  }

  /** - (f) allocate reionization_table with correct size */
  preio->reionization_table.resize(preio->re_size * number_of_redshifts);

  preio->rt_size = number_of_redshifts;

  /** - (g) retrieve data stored in the growing buffer */
  double* pData = reio_buffer.data();

  /** - (h) copy buffer to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (int i = 0; i < preio->rt_size; i++) {
    void* memcopy_result = memcpy(preio->reionization_table.data() + i * preio->re_size,
                                  pData + (preio->rt_size - i - 1) * preio->re_size,
                                  preio->re_size * sizeof(double));
    class_test(memcopy_result != preio->reionization_table.data() + i * preio->re_size,
               "cannot copy data back to reionization_temporary_table");
  }

  /** - (i) release the now-finished buffer (matches the previous gt_free() lifetime;
   *        the data lives on in reionization_table) */
  reio_buffer = std::vector<double>();

  /** - (j) another loop on z, to integrate equation for Tb and to compute cb2 */
  for (int i = preio->rt_size - 1; i > 0; i--) {
    z = preio->reionization_table[i * preio->re_size + preio->index_re_z];

    background_module_->background_tau_of_z(z, &tau);

    background_module_->background_at_tau(tau,
                                          pba->normal_info,
                                          pba->inter_normal,
                                          &last_index_back,
                                          pvecback);

    dz = (preio->reionization_table[i * preio->re_size + preio->index_re_z] -
          preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_z]);

    double opacity = (1. + z) * (1. + z) * n_e_ *
                     preio->reionization_table[i * preio->re_size + preio->index_re_xe] * _sigma_ *
                     _Mpc_over_m_;

    double mu = _m_H_ /
                (1. + (1. / _not4_ - 1.) * YHe_ +
                 preio->reionization_table[i * preio->re_size + preio->index_re_xe] * (1. - YHe_));

    /** - --> derivative of baryon temperature */

    double dTdz = 2. / (1 + z) *
                      preio->reionization_table[i * preio->re_size + preio->index_re_Tb] -
                  2. * mu / _m_e_ * 4. * all_species_.photons().Rho(pvecback) / 3. /
                      all_species_.baryons().Rho(pvecback) * opacity *
                      (pba->T_cmb * (1. + z) -
                       preio->reionization_table[i * preio->re_size + preio->index_re_Tb]) /
                      pvecback[background_module_->index_bg_H_];

    if (preco->annihilation > 0) {
      double energy_rate;
      thermodynamics_energy_injection(preco, z, &energy_rate);

      xe = preio->reionization_table[i * preio->re_size + preio->index_re_xe];
      const EnergyDeposition dep = energy_deposition_fractions(pth->chi_type, xe);

      dTdz += -2. / (3. * _k_B_) * energy_rate * dep.heat / (preco->Nnow * pow(1. + z, 3)) /
              (1. + preco->fHe +
               preio->reionization_table[i * preio->re_size + preio->index_re_xe]) /
              (pvecback[background_module_->index_bg_H_] * _c_ / _Mpc_over_m_ *
               (1. + z)); /* energy injection */
    }

    /** - --> increment baryon temperature */

    preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_Tb] =
        preio->reionization_table[i * preio->re_size + preio->index_re_Tb] - dTdz * dz;

    /** - --> get baryon equation of state */

    preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_wb] =
        _k_B_ / (_c_ * _c_ * mu) *
        preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_Tb];

    /** - --> get baryon adiabatic sound speed */

    preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_cb2] =
        preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_wb] *
        (1. + (1 + z) / 3. * dTdz /
                  preio->reionization_table[(i - 1) * preio->re_size + preio->index_re_Tb]);
  }

  /** - --> spline \f$ d \tau / dz \f$ with respect to z in view of integrating for optical depth */
  array_spline(preio->reionization_table.data(),
               preio->re_size,
               preio->rt_size,
               preio->index_re_z,
               preio->index_re_dkappadz,
               preio->index_re_d3kappadz3,
               _SPLINE_EST_DERIV_);

  /** - --> integrate for optical depth */
  array_integrate_all_spline(preio->reionization_table.data(),
                             preio->re_size,
                             preio->rt_size,
                             preio->index_re_z,
                             preio->index_re_dkappadz,
                             preio->index_re_d3kappadz3,
                             &(preio->reionization_optical_depth));
}

/**
 * Integrate thermodynamics with your favorite recombination code.
 *
 */

void ThermodynamicsModule::thermodynamics_recombination(recombination* preco, double* pvecback) {
  /* Both recombination codes go through the same integration: they differ only
     in the atomic physics supplying dx_H/dz and dx_He/dz, which is selected
     inside. */
  thermodynamics_recombination_integrate(preco, pvecback);
}

double ThermodynamicsModule::thermodynamics_recfast_hydrogen_saha_xH(const recombination* preco,
                                                                     double z) const {
  const double rhs = exp(1.5 * log(preco->CR * preco->Tnow / (1. + z)) -
                         preco->CB1 / (preco->Tnow * (1. + z))) /
                     preco->Nnow;
  return 0.5 * (sqrt(rhs * rhs + 4. * rhs) - rhs);
}

double ThermodynamicsModule::thermodynamics_recfast_helium_first_saha_xe(const recombination* preco,
                                                                         double z) const {
  const double rhs = exp(1.5 * log(preco->CR * preco->Tnow / (1. + z)) -
                         preco->CB1_He2 / (preco->Tnow * (1. + z))) /
                     preco->Nnow;
  return 0.5 * (sqrt(pow(rhs - 1. - preco->fHe, 2) + 4. * (1. + 2. * preco->fHe) * rhs) -
                (rhs - 1. - preco->fHe));
}

double ThermodynamicsModule::thermodynamics_recfast_helium_second_saha_xe(
    const recombination* preco, double z) const {
  const double rhs = 4. *
                     exp(1.5 * log(preco->CR * preco->Tnow / (1. + z)) -
                         preco->CB1_He1 / (preco->Tnow * (1. + z))) /
                     preco->Nnow;
  return 0.5 * (sqrt(pow(rhs - 1., 2) + 4. * (1. + preco->fHe) * rhs) - (rhs - 1.));
}

void ThermodynamicsModule::thermodynamics_recfast_store_row(
    recombination* preco, int sample_index, double z, double xe, double Tb, double dTbdz) const {
  const int table_index = preco->rt_size - sample_index - 1;

  preco->recombination_table[table_index * preco->re_size + preco->index_re_z]  = z;
  preco->recombination_table[table_index * preco->re_size + preco->index_re_xe] = xe;
  preco->recombination_table[table_index * preco->re_size + preco->index_re_Tb] = Tb;

  preco->recombination_table[table_index * preco->re_size + preco->index_re_wb] =
      _k_B_ / (_c_ * _c_ * _m_H_) *
      (1. + (1. / _not4_ - 1.) * preco->YHe + xe * (1. - preco->YHe)) * Tb;

  preco->recombination_table[table_index * preco->re_size + preco->index_re_cb2] =
      preco->recombination_table[table_index * preco->re_size + preco->index_re_wb] *
      (1. + (1. + z) * dTbdz / Tb / 3.);

  preco->recombination_table[table_index * preco->re_size + preco->index_re_dkappadtau] =
      (1. + z) * (1. + z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;
}

double ThermodynamicsModule::thermodynamics_recfast_xe_after_helium_ode(const recombination* preco,
                                                                        double z,
                                                                        const double* y) const {
  if (ppr->recfast_x_He0_trigger - y[1] < ppr->recfast_x_He0_trigger_delta) {
    const double x0_previous = thermodynamics_recfast_helium_second_saha_xe(preco, z);
    const double x0_new      = y[0] + preco->fHe * y[1];
    const double s      = (ppr->recfast_x_He0_trigger - y[1]) / ppr->recfast_x_He0_trigger_delta;
    const double weight = f2(s);
    return weight * x0_new + (1. - weight) * x0_previous;
  }

  return y[0] + preco->fHe * y[1];
}

double ThermodynamicsModule::thermodynamics_recfast_xe_after_full_ode(const recombination* preco,
                                                                      double z,
                                                                      const double* y) const {
  if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta) {
    const double x_H0   = thermodynamics_recfast_hydrogen_saha_xH(preco, z);
    const double s      = (ppr->recfast_x_H0_trigger - y[0]) / ppr->recfast_x_H0_trigger_delta;
    const double weight = f2(s);
    return weight * y[0] + (1. - weight) * x_H0 + preco->fHe * y[1];
  }

  return y[0] + preco->fHe * y[1];
}

/**
 * Integrate thermodynamics with RECFAST.
 *
 * Integrate thermodynamics with RECFAST, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination, from thermodynamics_init().
 *
 *
 *******************************************************************************
 * RECFAST is an integrator for Cosmic Recombination of Hydrogen and Helium,
 * developed by Douglas Scott (dscott@astro.ubc.ca)
 * based on calculations in the paper Seager, Sasselov & Scott
 * (ApJ, 523, L1, 1999).
 * and "fudge" updates in Wong, Moss & Scott (2008).
 *
 * Permission to use, copy, modify and distribute without fee or royalty at
 * any tier, this software and its documentation, for any purpose and without
 * fee or royalty is hereby granted, provided that you agree to comply with
 * the following copyright notice and statements, including the disclaimer,
 * and that the same appear on ALL copies of the software and documentation,
 * including modifications that you make for internal use or for distribution:
 *
 * Copyright 1999-2010 by University of British Columbia.  All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO
 * REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
 * BY WAY OF EXAMPLE, BUT NOT LIMITATION,
 * U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF
 * MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
 * THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
 * ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
 *******************************************************************************
 *
 * Version 1.5: includes extra fitting function from
 *              Rubino-Martin et al. arXiv:0910.4383v1 [astro-ph.CO]
 *
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_recombination_integrate(recombination* preco,
                                                                  double* pvecback) {
  /** Summary: */

  /** - define local variables */

  /* vector of variables to be integrated: x_H, x_He, Tmat */
  double y[3], dy[3];

  /* contains all fixed parameters which should be passed to thermodynamics_recombination_derivs */
  thermodynamics_parameters_and_workspace tpaw{this};

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  preco->recombination_table.resize(preco->re_size * preco->rt_size);

  /** - read a few precision/cosmological parameters */

  /* Nz */
  int Nz = ppr->recfast_Nz0;

  /* preco->H0 is H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;

  /* Omega_b */
  double OmegaB = pba->Omega0_b;

  /* Yp */
  preco->YHe = YHe_;

  /* Tnow */
  preco->Tnow = pba->T_cmb;

  /* z_initial */
  double zinitial = ppr->recfast_z_initial;

  /* H_frac */
  preco->H_frac = ppr->recfast_H_frac;

  /* H fudging */
  class_test((ppr->recfast_Hswitch != 0) && (ppr->recfast_Hswitch != 1),
             "RECFAST error: unknown H fudging scheme");
  preco->fu = ppr->recfast_fudge_H;
  if (ppr->recfast_Hswitch)
    preco->fu += ppr->recfast_delta_fudge_H;

  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
             "RECFAST error: unknown He fudging scheme");

  /* related quantities */
  double z    = zinitial;
  double mu_H = 1. / (1. - preco->YHe);
  //mu_T = _not4_ /(_not4_ - (_not4_-1.)*preco->YHe); /* recfast 1.4*/
  preco->fHe  = preco->YHe / (_not4_ * (1. - preco->YHe)); /* recfast 1.4 */
  preco->Nnow = 3. * preco->H0 * preco->H0 * OmegaB / (8. * _PI_ * _G_ * mu_H * _m_H_);
  n_e_        = preco->Nnow;

  /* energy injection parameters */
  preco->annihilation           = pth->annihilation;
  preco->has_on_the_spot        = pth->has_on_the_spot;
  preco->annihilation_variation = pth->annihilation_variation;
  preco->annihilation_z         = pth->annihilation_z;
  preco->annihilation_zmax      = pth->annihilation_zmax;
  preco->annihilation_zmin      = pth->annihilation_zmin;
  preco->decay                  = pth->decay;
  preco->annihilation_f_halo    = pth->annihilation_f_halo;
  preco->annihilation_z_halo    = pth->annihilation_z_halo;

  /* quantities related to constants defined in thermodynamics.h */
  //n = preco->Nnow * pow((1.+z),3);
  double Lalpha    = 1. / _L_H_alpha_;
  double Lalpha_He = 1. / _L_He_2p_;
  double DeltaB    = _h_P_ * _c_ * (_L_H_ion_ - _L_H_alpha_);
  preco->CDB       = DeltaB / _k_B_;
  double DeltaB_He = _h_P_ * _c_ * (_L_He1_ion_ - _L_He_2s_);
  preco->CDB_He    = DeltaB_He / _k_B_;
  preco->CB1       = _h_P_ * _c_ * _L_H_ion_ / _k_B_;
  preco->CB1_He1   = _h_P_ * _c_ * _L_He1_ion_ / _k_B_;
  preco->CB1_He2   = _h_P_ * _c_ * _L_He2_ion_ / _k_B_;
  preco->CR        = 2. * _PI_ * (_m_e_ / _h_P_) * (_k_B_ / _h_P_);
  preco->CK        = pow(Lalpha, 3) / (8. * _PI_);
  preco->CK_He     = pow(Lalpha_He, 3) / (8. * _PI_);
  preco->CL        = _c_ * _h_P_ / (_k_B_ * Lalpha);
  preco->CL_He     = _c_ * _h_P_ / (_k_B_ / _L_He_2s_);
  preco->CT        = (8. / 3.) * (_sigma_ / (_m_e_ * _c_)) *
                     (8. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 3));

  preco->Bfact = _h_P_ * _c_ * (_L_He_2p_ - _L_He_2s_) / _k_B_;

  /** - define the fields of the 'thermodynamics parameter and workspace' structure */
  tpaw.preco    = preco;
  tpaw.pvecback = pvecback;

  /* The atomic physics behind dx_H/dz and dx_He/dz. Everything else about the
     integration -- the Saha phases, the baryon temperature, the evolver, the
     sampling -- is shared, so the two recombination codes differ in exactly this
     one place. In particular both are handed CLASS's own H(z) at every step. */
  std::unique_ptr<RecombinationModel> recombination_model;
  if (pth->recombination == hyrec) {
    /* SWIFT's correction function is calibrated against a fiducial cosmology and
       takes the NON-free-streaming matter density: total matter less the massive
       neutrinos, which stream out of the perturbations it was fitted to. Adding
       them here is exactly the double-count that made the old HyRec path wrong
       (#396), reached through a different term. */
    recombination_model = std::make_unique<HyrecModel>(ppr->hyrec_path,
                                                       pba->T_cmb,
                                                       pba->Omega0_b * pba->h * pba->h,
                                                       (background_module_->Omega0_m_ -
                                                        background_module_->GetOmega0NcdmTot()) *
                                                           pba->h * pba->h,
                                                       pba->Omega0_k * pba->h * pba->h,
                                                       background_module_->Neff_,
                                                       YHe_,
                                                       preco->fHe,
                                                       preco->Nnow,
                                                       ppr->recfast_z_initial);
  }
  else {
    recombination_model = std::make_unique<RecfastModel>(ppr, preco);
  }
  tpaw.recombination_model = recombination_model.get();

  if (pth->thermodynamics_verbose > 0)
    printf(" -> recombination computed with %s\n", recombination_model->Name());

  /** - impose initial conditions at early times */

  class_test(zinitial < ppr->recfast_z_He_3,
             "increase zinitial, otherwise should get initial conditions from recfast's get_init "
             "routine (less precise anyway)");

  y[0] = 1.;
  y[1] = 1.;
  double x0;
  y[2] = preco->Tnow * (1. + z);

  auto sample_z = [zinitial, Nz](int index_z) {
    return zinitial * (double) (Nz - index_z - 1) / (double) Nz;
  };

  auto first_sample_below = [&](int first_index, double z_boundary) {
    int index_z = first_index;
    while ((index_z < Nz) && (sample_z(index_z) >= z_boundary))
      index_z++;
    return index_z;
  };

  auto sample_count_down_to = [&](int first_index, double z_boundary) {
    int count = 0;
    while ((first_index + count < Nz) && (sample_z(first_index + count) >= z_boundary))
      count++;
    return count;
  };

  auto bisection_root = [&](auto&& function, double z_low, double z_high) {
    double f_low  = function(z_low);
    double f_high = function(z_high);

    class_test(f_low * f_high > 0.,
               "RECFAST boundary is not bracketed: f(%e)=%e, f(%e)=%e",
               z_low,
               f_low,
               z_high,
               f_high);

    for (int iter = 0; iter < 200; iter++) {
      const double z_mid = 0.5 * (z_low + z_high);
      const double f_mid = function(z_mid);

      if (fabs(z_high - z_low) <= ppr->smallest_allowed_variation * std::max(1., z_mid))
        return z_mid;

      if (f_mid * f_low <= 0.) {
        z_high = z_mid;
        f_high = f_mid;
      }
      else {
        z_low = z_mid;
        f_low = f_mid;
      }
    }

    return 0.5 * (z_low + z_high);
  };

  auto recfast_analytic_second_helium_xHe = [&](double z_analytic) {
    double x0_analytic = thermodynamics_recfast_helium_second_saha_xe(preco, z_analytic);

    if (z_analytic > ppr->recfast_z_He_3 - ppr->recfast_delta_z_He_3) {
      const double x0_previous = 1. + preco->fHe;
      const double s           = (ppr->recfast_z_He_3 - z_analytic) / ppr->recfast_delta_z_He_3;
      const double weight      = f1(s);
      x0_analytic              = weight * x0_analytic + (1. - weight) * x0_previous;
    }

    return (x0_analytic - 1.) / preco->fHe;
  };

  const double z_helium_ode_start = bisection_root(
      [&](double z_root) {
        return recfast_analytic_second_helium_xHe(z_root) - ppr->recfast_x_He0_trigger;
      },
      0.,
      ppr->recfast_z_He_3 + ppr->recfast_delta_z_He_3);

  const double z_hydrogen_ode_start = bisection_root(
      [&](double z_root) {
        return thermodynamics_recfast_hydrogen_saha_xH(preco, z_root) - ppr->recfast_x_H0_trigger;
      },
      0.,
      z_helium_ode_start);

  std::vector<double> minus_z_sampling(Nz);
  for (int i = 0; i < Nz; i++) {
    const double zend   = sample_z(i);
    minus_z_sampling[i] = -zend;
  }

  std::vector<int> used_in_output_full(_RECFAST_INTEG_SIZE_, 1);
  std::vector<int> used_in_output_helium(2, 1);

  /* RECFAST's recombination ODE is stiff. The implicit ndf15 (default) and the adaptive
     rkdp45 both handle it; the explicit Cash-Karp rk does not (it yields an unphysical
     history and trips the z_rec sanity check). When rk is requested, fall back to rkdp45
     for the recombination integration only, leaving the user's evolver choice intact for
     the background and perturbation modules. */
  /* Configure the shared explicit-RK controller from THIS module's ppr, here
     rather than at parse time: Cosmology is lazy, so a second object parsed in
     between would otherwise be the one whose settings this run reads. */
  evolver_erk_configure(ppr->erk_controller_config());

  auto generic_evolver = &evolver_ndf15;
  if (ppr->evolver_thermodynamics == evolver_type::rkdp45) {
    generic_evolver = &evolver_rkdp45;
  }
  else if (ppr->evolver_thermodynamics == evolver_type::tsit5) {
    generic_evolver = &evolver_tsit5;
  }
  else if (ppr->evolver_thermodynamics == evolver_type::rk) {
    generic_evolver = &evolver_rkdp45;
    printf(
        "\nWarning: evolver=rk cannot integrate the stiff RECFAST system; using rkdp45 for "
        "recombination instead.\n");
  }
  else if (ppr->evolver_thermodynamics == evolver_type::etd) {
    // RECFAST's stiffness is not diagonal in the sense ETD exploits, and no species
    // reports a diagonal here, so ETD would be explicit Heun on a stiff system.
    generic_evolver = &evolver_rkdp45;
    printf(
        "\nWarning: evolver=etd applies to the background only; using rkdp45 for "
        "recombination instead.\n");
  }

  /** - Fill all analytic RECFAST regimes directly on the requested sampling. */

  tpaw.recombination_phase = RecombinationPhase::analytic;

  int i = 0;
  for (; (i < Nz) && (sample_z(i) > z_helium_ode_start); i++) {
    z = sample_z(i);

    double x_H0        = 1.;
    double x_He0       = 1.;
    double x0_previous = 0.;
    double x0_new      = 0.;
    double s           = 0.;
    double weight      = 0.;

    /** - --> first approximation: H and Helium fully ionized */

    if (z > ppr->recfast_z_He_1 + ppr->recfast_delta_z_He_1) {
      x0 = 1. + 2. * preco->fHe;
    }

    /** - --> second approximation: first Helium recombination (analytic approximation) */

    else if (z > ppr->recfast_z_He_2 + ppr->recfast_delta_z_He_2) {
      if (z > ppr->recfast_z_He_1 - ppr->recfast_delta_z_He_1) {
        x0_previous = 1. + 2. * preco->fHe;
        x0_new      = thermodynamics_recfast_helium_first_saha_xe(preco, z);
        s           = (ppr->recfast_z_He_1 - z) / ppr->recfast_delta_z_He_1;
        weight      = f1(s);
        x0          = weight * x0_new + (1. - weight) * x0_previous;
      }
      else {
        x0 = thermodynamics_recfast_helium_first_saha_xe(preco, z);
      }
    }

    /** - --> third approximation: first Helium recombination completed */

    else if (z > ppr->recfast_z_He_3 + ppr->recfast_delta_z_He_3) {
      if (z > ppr->recfast_z_He_2 - ppr->recfast_delta_z_He_2) {
        x0_previous = thermodynamics_recfast_helium_first_saha_xe(preco, z);
        x0_new      = 1. + preco->fHe;
        s           = (ppr->recfast_z_He_2 - z) / ppr->recfast_delta_z_He_2;
        weight      = f1(s);
        x0          = weight * x0_new + (1. - weight) * x0_previous;
      }
      else {
        x0 = 1. + preco->fHe;
      }
    }

    /** - --> fourth approximation: second Helium recombination starts (analytic approximation) */

    else if (y[1] > ppr->recfast_x_He0_trigger) {
      x0_new = thermodynamics_recfast_helium_second_saha_xe(preco, z);

      if (z > ppr->recfast_z_He_3 - ppr->recfast_delta_z_He_3) {
        x0_previous = 1. + preco->fHe;
        s           = (ppr->recfast_z_He_3 - z) / ppr->recfast_delta_z_He_3;
        weight      = f1(s);
        x0          = weight * x0_new + (1. - weight) * x0_previous;
      }
      else {
        x0 = x0_new;
      }

      x_He0 = (x0 - 1.) / preco->fHe;
    }

    else {
      class_stop("RECFAST analytic sampling crossed the Helium ODE boundary unexpectedly");
    }

    y[0] = x_H0;
    y[1] = x_He0;
    y[2] = preco->Tnow * (1. + z);

    thermodynamics_recombination_derivs(z, y, dy, &tpaw);
    thermodynamics_recfast_store_row(preco, i, z, x0, y[2], dy[2]);
  }

  /** - Evolve Helium and baryon temperature while Hydrogen follows its Saha branch. */

  y[0] = thermodynamics_recfast_hydrogen_saha_xH(preco, z_helium_ode_start);
  y[1] = ppr->recfast_x_He0_trigger;
  y[2] = preco->Tnow * (1. + z_helium_ode_start);

  const int first_helium_sample = first_sample_below(i, z_helium_ode_start);
  const int helium_sample_count = sample_count_down_to(first_helium_sample, z_hydrogen_ode_start);

  /* The analytic loop filled rows [0, i) and the helium phase begins at
     first_helium_sample; they must meet exactly. They only differ if a requested sample
     lands exactly on z_helium_ode_start (the analytic loop's strict '>' skips it while
     first_sample_below's '>=' also skips it), which would leave that row unwritten. This
     is unreachable in practice but cheap to guard. */
  class_test(first_helium_sample != i,
             "RECFAST sampling gap at the analytic/helium boundary (a sample coincided with the "
             "helium transition redshift); nudge recfast_Nz0 or recfast_z_He_3");

  /* The helium ODE must run even when no requested sample falls inside its redshift
     window (only at pathologically small recfast_Nz0), because the full phase needs the
     evolved {x_He, Tmat} handed off at z_hydrogen_ode_start. In that case feed the
     evolver a single boundary sample and store nothing. */
  double y_helium[2]              = {y[1], y[2]};
  double helium_boundary_sample   = -z_hydrogen_ode_start;
  double* helium_sampling         = (helium_sample_count > 0)
                                        ? minus_z_sampling.data() + first_helium_sample
                                        : &helium_boundary_sample;
  const int helium_sampling_count = (helium_sample_count > 0) ? helium_sample_count : 1;
  auto helium_output              = (helium_sample_count > 0) ? &thermodynamics_recfast_output
                                                              : &thermodynamics_recfast_output_none;

  tpaw.recombination_phase         = RecombinationPhase::helium;
  tpaw.recfast_output_index_offset = first_helium_sample;

  generic_evolver(thermodynamics_recfast_derivs,
                  -z_helium_ode_start,
                  -z_hydrogen_ode_start,
                  y_helium,
                  used_in_output_helium.data(),
                  2,
                  &tpaw,
                  ppr->tol_thermo_integration,
                  ppr->smallest_allowed_variation,
                  thermodynamics_recfast_timescale,
                  z_helium_ode_start - z_hydrogen_ode_start,
                  helium_sampling,
                  helium_sampling_count,
                  helium_output,
                  nullptr,
                  // No species reports a diagonal on the thermodynamics path.
                  nullptr);

  y[0] = ppr->recfast_x_H0_trigger;
  y[1] = y_helium[0];
  y[2] = y_helium[1];

  /** - Evolve Hydrogen, Helium, and baryon temperature with dense-output sampling. */

  const int first_full_sample = first_sample_below(first_helium_sample + helium_sample_count,
                                                   z_hydrogen_ode_start);

  if (first_full_sample < Nz) {
    tpaw.recombination_phase         = RecombinationPhase::full;
    tpaw.recfast_output_index_offset = first_full_sample;

    generic_evolver(thermodynamics_recfast_derivs,
                    -z_hydrogen_ode_start,
                    minus_z_sampling[Nz - 1],
                    y,
                    used_in_output_full.data(),
                    _RECFAST_INTEG_SIZE_,
                    &tpaw,
                    ppr->tol_thermo_integration,
                    ppr->smallest_allowed_variation,
                    thermodynamics_recfast_timescale,
                    z_hydrogen_ode_start,
                    minus_z_sampling.data() + first_full_sample,
                    Nz - first_full_sample,
                    thermodynamics_recfast_output,
                    nullptr,
                    nullptr);
  }
}

/**
 * Subroutine evaluating the derivative with respect to redshift of
 * thermodynamical quantities (from RECFAST version 1.4).
 *
 * Computes derivatives of the three variables to integrate: \f$ d x_H
 * / dz, d x_{He} / dz, d T_{mat} / dz \f$.
 *
 * This is one of the few functions in the code which are passed to
 * the generic evolver routines. Since the evolver interface should work with
 * functions passed from various modules, the argument format is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic
 *   pointer. Here, this pointer contains the precision, background
 *   and recombination structures, plus a background vector, but the evolver
 *   does not know its fine structure.
 *
 * @param z                        Input: redshift
 * @param y                        Input: vector of variable to integrate
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices) and workspace (already allocated)
 */

void ThermodynamicsModule::thermodynamics_recombination_derivs_member(
    double z, double* y, double* dy, void* parameters_and_workspace) {
  /* define local variables */

  struct thermodynamics_parameters_and_workspace* ptpaw =
      (struct thermodynamics_parameters_and_workspace*) parameters_and_workspace;
  struct recombination* preco = ptpaw->preco;
  double* pvecback            = ptpaw->pvecback;

  /* Approximation switches are owned by the driver phase, not re-derived from the
     state here (mirrors the perturbations module): hydrogen stays frozen on its Saha
     branch until the full phase, and the RECFAST 1.4 helium corrections run once
     helium is recombining (helium and full phases). */
  const RecombinationPhase phase = ptpaw->recombination_phase;
  const bool hydrogen_frozen     = (phase != RecombinationPhase::full);
  const bool helium_corrections  = (phase != RecombinationPhase::analytic);

  double x_H  = y[0];
  double x_He = y[1];
  double x    = x_H + preco->fHe * x_He;
  double Tmat = y[2];

  double n    = preco->Nnow * (1. + z) * (1. + z) * (1. + z);
  double n_He = preco->fHe * n;
  double Trad = preco->Tnow * (1. + z);

  double tau;
  background_module_->background_tau_of_z(z, &tau);

  int last_index_back;
  background_module_->background_at_tau(tau,
                                        pba->short_info,
                                        pba->inter_normal,
                                        &last_index_back,
                                        pvecback);

  double energy_rate;
  thermodynamics_energy_injection(preco, z, &energy_rate);

  /* Fractions of that rate reaching each deposition channel. Skipped entirely
     when nothing is being injected -- the default -- so the usual path does not
     pay for a spline evaluation per derivative call. */
  const EnergyDeposition dep = (energy_rate > 0.) ? energy_deposition_fractions(pth->chi_type, x)
                                                  : EnergyDeposition{0., 0., 0., 0., 0.};

  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  double Hz = pvecback[background_module_->index_bg_H_] * _c_ / _Mpc_over_m_;

  const RecombinationState state =
      {z, x_H, x_He, x, n, Hz, Tmat, Trad, hydrogen_frozen, helium_corrections};

  const IonisationDerivatives dx = ptpaw->recombination_model->Derivatives(state, dep, energy_rate);
  dy[0]                          = dx.dx_H_dz;
  dy[1]                          = dx.dx_He_dz;

  double timeTh = (1. / (preco->CT * pow(Trad, 4))) * (1. + x + preco->fHe) / x;
  double timeH  = 2. / (3. * preco->H0 * pow(1. + z, 1.5));

  if (timeTh < preco->H_frac * timeH) {
    /*   dy[2]=Tmat/(1.+z); */
    /* v 1.5: like in camb, add here a smoothing term as suggested by Adam Moss */
    double dHdz    = -pvecback[background_module_->index_bg_H_prime_] /
                     pvecback[background_module_->index_bg_H_] * _c_ / _Mpc_over_m_;
    double epsilon = Hz * (1. + x + preco->fHe) / (preco->CT * pow(Trad, 3) * x);
    dy[2]          = preco->Tnow +
                     epsilon * ((1. + preco->fHe) / (1. + preco->fHe + x)) *
                         ((dy[0] + preco->fHe * dy[1]) / x) -
                     epsilon * dHdz / Hz + 3. * epsilon / (1. + z);
  }
  else {
    /* equations modified to take into account energy injection from dark matter */

    dy[2] = preco->CT * pow(Trad, 4) * x / (1. + x + preco->fHe) * (Tmat - Trad) / (Hz * (1. + z)) +
            2. * Tmat / (1. + z) -
            2. / (3. * _k_B_) * energy_rate * dep.heat / n / (1. + preco->fHe + x) /
                (Hz * (1. + z)); /* energy injection */
  }
}

/**
 * This routine merges the two tables 'recombination_table' and
 * 'reionization_table' inside the table 'thermodynamics_table', and
 * frees the temporary structures 'recombination' and 'reionization'.
 *
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input: pointer to reionization structure
 * @return the error status
 */

void ThermodynamicsModule::thermodynamics_merge_reco_and_reio(recombination* preco,
                                                              reionization* preio) {
  /** Summary: */

  /** - define local variables */

  int index_th, index_re;

  /** - first, a little check that the two tables match each other and can be merged */

  if (pth->reio_parametrization != reio_none) {
    class_test(preco->recombination_table[preio->index_reco_when_reio_start * preco->re_size +
                                          preco->index_re_z] !=
                   preio->reionization_table[(preio->rt_size - 1) * preio->re_size +
                                             preio->index_re_z],
               "mismatch which should never happen");
  }

  /** - find number of redshift in full table = number in reco + number in reio - overlap */

  tt_size_ = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;

  /** - add  more points to start earlier in presence of interacting DM */

  if (all_species_.count("IDM_DR_IDR") > 0)
    tt_size_ += ppr->thermo_Nz1_idm_dr + ppr->thermo_Nz2_idm_dr - 1;

  /** - allocate arrays in thermo structure */

  z_table_.resize(tt_size_);
  thermodynamics_table_.resize(th_size_ * tt_size_);
  d2thermodynamics_dz2_table_.resize(th_size_ * tt_size_);

  /** - fill these arrays */

  for (int i = 0; i < preio->rt_size; i++) {
    z_table_[i] = preio->reionization_table[i * preio->re_size + preio->index_re_z];
    thermodynamics_table_[i * th_size_ + index_th_xe_] =
        preio->reionization_table[i * preio->re_size + preio->index_re_xe];
    thermodynamics_table_[i * th_size_ + index_th_dkappa_] =
        preio->reionization_table[i * preio->re_size + preio->index_re_dkappadtau];
    thermodynamics_table_[i * th_size_ + index_th_Tb_] =
        preio->reionization_table[i * preio->re_size + preio->index_re_Tb];
    thermodynamics_table_[i * th_size_ + index_th_wb_] =
        preio->reionization_table[i * preio->re_size + preio->index_re_wb];
    thermodynamics_table_[i * th_size_ + index_th_cb2_] =
        preio->reionization_table[i * preio->re_size + preio->index_re_cb2];
  }
  for (int i = 0; i < ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1; i++) {
    index_th           = i + preio->rt_size;
    index_re           = i + preio->index_reco_when_reio_start + 1;
    z_table_[index_th] = preco->recombination_table[index_re * preco->re_size + preco->index_re_z];
    thermodynamics_table_[index_th * th_size_ + index_th_xe_] =
        preco->recombination_table[index_re * preco->re_size + preco->index_re_xe];
    thermodynamics_table_[index_th * th_size_ + index_th_dkappa_] =
        preco->recombination_table[index_re * preco->re_size + preco->index_re_dkappadtau];
    thermodynamics_table_[index_th * th_size_ + index_th_Tb_] =
        preco->recombination_table[index_re * preco->re_size + preco->index_re_Tb];
    thermodynamics_table_[index_th * th_size_ + index_th_wb_] =
        preco->recombination_table[index_re * preco->re_size + preco->index_re_wb];
    thermodynamics_table_[index_th * th_size_ + index_th_cb2_] =
        preco->recombination_table[index_re * preco->re_size + preco->index_re_cb2];
  }

  /** - add more points at larger redshift in presence of interacting
        DM. This is necessary because the value of integrated
        quantitites like tau_idm_dr or tau_idr will then be computed
        exactly up to high redshift. With extrapolations in
        thermodynamics_at_z() we could not obtain this. */

  if (all_species_.count("IDM_DR_IDR") > 0) {
    for (int i = 0; i < ppr->thermo_Nz2_idm_dr + ppr->thermo_Nz1_idm_dr - 1; i++) {
      /* with an intermediate step Delta z = (thermo_z_initial_idm_dr-recfast_z_initial)/thermo_Nz1_idm_dr/thermo_Nz1_idm_dr */
      if (i < ppr->thermo_Nz2_idm_dr - 1) {
        index_th = i + preio->rt_size + ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1;
        z_table_[index_th] = ppr->recfast_z_initial +
                             ((double) i + 1.) *
                                 (ppr->thermo_z_initial_idm_dr - ppr->recfast_z_initial) /
                                 (double) ppr->thermo_Nz1_idm_dr / (double) ppr->thermo_Nz2_idm_dr;
      }
      /* with a large step Delta z  = (thermo_z_initial_idm_dr-recfast_z_initial)/thermo_Nz1_idm_dr */
      else {
        index_th           = (i - ppr->thermo_Nz2_idm_dr + 1) + preio->rt_size + ppr->recfast_Nz0 -
                             preio->index_reco_when_reio_start - 1 + ppr->thermo_Nz2_idm_dr - 1;
        z_table_[index_th] = ppr->recfast_z_initial +
                             ((double) (i - ppr->thermo_Nz2_idm_dr + 1) + 1.) *
                                 (ppr->thermo_z_initial_idm_dr - ppr->recfast_z_initial) /
                                 (double) ppr->thermo_Nz1_idm_dr;
      }
      /* same extrapolation formulas as in thermodynamics_at_z() */
      double x0 = thermodynamics_table_[(preio->rt_size + ppr->recfast_Nz0 -
                                         preio->index_reco_when_reio_start - 2) *
                                            th_size_ +
                                        index_th_xe_];
      thermodynamics_table_[index_th * th_size_ + index_th_xe_]     = x0;
      thermodynamics_table_[index_th * th_size_ + index_th_dkappa_] = (1. + z_table_[index_th]) *
                                                                      (1. + z_table_[index_th]) *
                                                                      n_e_ * x0 * _sigma_ *
                                                                      _Mpc_over_m_;
      thermodynamics_table_[index_th * th_size_ + index_th_Tb_]     = pba->T_cmb *
                                                                      (1. + z_table_[index_th]);
      thermodynamics_table_[index_th * th_size_ + index_th_wb_] =
          _k_B_ / (_c_ * _c_ * _m_H_) * (1. + (1. / _not4_ - 1.) * YHe_ + x0 * (1. - YHe_)) *
          pba->T_cmb * (1. + z_table_[index_th]);
      thermodynamics_table_[index_th * th_size_ + index_th_cb2_] =
          thermodynamics_table_[index_th * th_size_ + index_th_wb_] * 4. / 3.;
    }
  }

  /** - release the temporary recombination/reionization tables */

  preco->recombination_table.clear();
  preco->recombination_table.shrink_to_fit();

  if (pth->reio_parametrization != reio_none) {
    preio->reionization_table.clear();
    preio->reionization_table.shrink_to_fit();
  }
}

/**
 * Subroutine for formatting thermodynamics output
 */

void ThermodynamicsModule::thermodynamics_output_titles(std::string& titles) const {
  class_store_columntitle(titles, "z", true);
  class_store_columntitle(titles, "conf. time [Mpc]", true);
  class_store_columntitle(titles, "x_e", true);
  class_store_columntitle(titles, "kappa' [Mpc^-1]", true);
  class_store_columntitle(titles, "exp(-kappa)", true);
  class_store_columntitle(titles, "g [Mpc^-1]", true);
  class_store_columntitle(titles, "Tb [K]", true);
  class_store_columntitle(titles, "w_b", true);
  class_store_columntitle(titles, "c_b^2", true);
  class_store_columntitle(titles, "tau_d", true);
  class_store_columntitle(titles, "r_d", pth->compute_damping_scale);

  if (all_species_.count("IDM_DR_IDR") > 0) {
    class_store_columntitle(titles, "dmu_idm_dr", true);
    class_store_columntitle(titles, "tau_idm_dr", true);
    class_store_columntitle(titles, "tau_idr", true);
    class_store_columntitle(titles, "g_idm_dr [Mpc^-1]", true);
    class_store_columntitle(titles, "c_idm_dr^2", true);
    class_store_columntitle(titles, "T_idm_dr", true);
    class_store_columntitle(titles, "dmu_idr", true);
  }
}

void ThermodynamicsModule::thermodynamics_output_data(int number_of_titles, double* data) const {
  int storeidx;
  double* dataptr;
  double z, tau;

  /* Store quantities: */
  for (int index_z = 0; index_z < tt_size_; index_z++) {
    dataptr                  = data + index_z * number_of_titles;
    const double* pvecthermo = thermodynamics_table_.data() + index_z * th_size_;
    z                        = z_table_[index_z];
    storeidx                 = 0;

    background_module_->background_tau_of_z(z, &tau);

    class_store_double(dataptr, z, true, storeidx);
    class_store_double(dataptr, tau, true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_xe_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_dkappa_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_exp_m_kappa_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_g_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_Tb_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_wb_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_cb2_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_tau_d_], true, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_r_d_], pth->compute_damping_scale, storeidx);

    if (all_species_.count("IDM_DR_IDR") > 0) {
      class_store_double(dataptr, pvecthermo[index_th_dmu_idm_dr_], true, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_tau_idm_dr_], true, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_tau_idr_], true, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_g_idm_dr_], true, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_cidm_dr2_], true, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_Tidm_dr_], true, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_dmu_idr_], true, storeidx);
    }
  }
}

void ThermodynamicsModule::thermodynamics_tanh(
    double x, double center, double before, double after, double width, double* result) {
  *result = before + (after - before) * (tanh((x - center) / width) + 1.) / 2.;
}
