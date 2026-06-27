/** @file background_module.cpp
 * Integrates the cosmological background and provides interpolated background
 * quantities for the lifetime of a BackgroundModule instance.
 */

#include "background_module.h"

#include <algorithm>

#include "../species/background_ic_context.h"
#include "../species/dcdm_dr_species.h"
#include "../species/dncdm_dr_species.h"
#include "../species/fluid.h"
#include "../species/idm_dr_idr_species.h"
#include "../species/idm_drmd_idr_drmd_species.h"
#include "../species/ncdm_species.h"
#include "../species/scalar_field.h"
#include "../species/type3_species.h"
#include "bisection.h"

/**
 * Return all NCDMBaseSpecies pointers from all_species_ in deterministic order.
 * For DNCDM_DR_Species composites, returns the wrapped DNCDMSpecies.
 * Since all_species_ is a std::map (sorted by key), "NCDM_0" < "NCDM_1" < ...
 * so the order is automatically deterministic.
 */
static std::vector<NCDMBaseSpecies*> GetNcdmSpecies(const SpeciesCollection& all_species) {
  std::vector<NCDMBaseSpecies*> result;
  for (const auto& sp : all_species) {
    if (auto* ncdm = dynamic_cast<NCDMBaseSpecies*>(sp.get()))
      result.push_back(ncdm);
    else if (auto* dncdm_dr = dynamic_cast<DNCDM_DR_Species*>(sp.get()))
      result.push_back(&dncdm_dr->dncdm());
  }
  return result;
}

BackgroundModule::BackgroundModule(InputModulePtr input_module) : BaseModule(input_module) {
  background_init();
}

BackgroundModule::~BackgroundModule() {}

double BackgroundModule::GetOmega0Species(const std::string& key) const {
  // Top-level species first.
  if (auto* ptr = all_species_.find(key))
    return (*ptr)->GetOmega0();
  // Sub-species of composites: IDR / IDM_DR live in IDM_DR_IDR, etc.
  if (auto* ptr = all_species_.find("IDM_DR_IDR")) {
    const auto& comp = static_cast<const IDM_DR_IDR_Species&>(**ptr);
    if (key == "IDR")
      return comp.idr().GetOmega0();
    if (key == "IDM_DR")
      return comp.idm_dr().GetOmega0();
  }
  if (auto* ptr = all_species_.find("IDM_DRMD_IDR_DRMD")) {
    const auto& comp = static_cast<const IDM_DRMD_IDR_DRMD_Species&>(**ptr);
    if (key == "IDR_DRMD")
      return comp.idr_drmd().GetOmega0();
    if (key == "IDM_DRMD")
      return comp.idm_drmd().GetOmega0();
  }
  if (auto* ptr = all_species_.find("DCDM_DR")) {
    const auto& comp = static_cast<const DCDM_DR_Species&>(**ptr);
    if (key == "DCDM")
      return comp.dcdm().GetOmega0();
  }
  return 0.;
}

double BackgroundModule::GetSpeciesParam(const std::string& key, const std::string& param) const {
  if (auto* ptr = all_species_.find(key)) {
    if (auto val = (*ptr)->GetParam(param))
      return *val;
  }
  return 0.;
}

// Wrapper functions to pass non-static member functions
int BackgroundModule::background_derivs_loga(double loga,
                                             double* y,
                                             double* dy,
                                             void* parameters_and_workspace) {
  auto pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  return pbpaw->background_module->background_derivs_loga_member(loga,
                                                                 y,
                                                                 dy,
                                                                 parameters_and_workspace);
}
int BackgroundModule::background_add_line_to_bg_table(
    double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace) {
  auto pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  return pbpaw->background_module->background_add_line_to_bg_table_member(loga,
                                                                          y,
                                                                          dy,
                                                                          index_loga,
                                                                          parameters_and_workspace);
}

/**
 * Evolution timescale for the background integration. This only exists to feed
 * the legacy explicit Runge-Kutta evolver (evolver=rk), which asks the caller
 * for a step size rather than judging the evolution speed itself; the stiff
 * ndf15 default and the modern rkdp45 evolver size their own steps and ignore
 * it. The background is integrated in loga (see background_derivs_loga_member),
 * where a unit step is the natural variation scale, so the timescale is simply
 * constant. evolver_rk then advances in steps of timestep_over_timescale, which
 * exceeds the requested output spacing, so it effectively integrates from one
 * output point to the next under the control of tol_background_integration.
 */
int BackgroundModule::background_timescale(double loga,
                                           void* parameters_and_workspace,
                                           double* timescale) {
  (void) loga;
  (void) parameters_and_workspace;
  *timescale = 1.;
  return _SUCCESS_;
}

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

void BackgroundModule::background_at_tau(
    double tau,
    short return_format,
    short intermode,
    int* last_index,
    double*
        pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
) const {
  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < tau_table_[0],
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter "
             "a_ini_over_a_today_default\n",
             tau,
             tau_table_[0]);

  class_test(tau > tau_table_[bt_size_ - 1],
             "out of range: tau=%e > tau_max=%e\n",
             tau,
             tau_table_[bt_size_ - 1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size = bg_size_normal_;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size = bg_size_short_;
    }
    else {
      pvecback_size = bg_size_;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    array_interpolate_spline(const_cast<double*>(tau_table_.data()),
                             bt_size_,
                             const_cast<double*>(background_table_.data()),
                             const_cast<double*>(d2background_dtau2_table_.data()),
                             bg_size_,
                             tau,
                             last_index,
                             pvecback,
                             pvecback_size);
  }
  if (intermode == pba->inter_closeby) {
    array_interpolate_spline_growing_closeby(const_cast<double*>(tau_table_.data()),
                                             bt_size_,
                                             const_cast<double*>(background_table_.data()),
                                             const_cast<double*>(d2background_dtau2_table_.data()),
                                             bg_size_,
                                             tau,
                                             last_index,
                                             pvecback,
                                             pvecback_size);
  }
}

/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

void BackgroundModule::background_tau_of_z(double z, double* tau) const {
  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < z_table_[bt_size_ - 1],
             "out of range: z=%e < z_min=%e\n",
             z,
             z_table_[bt_size_ - 1]);

  class_test(z > z_table_[0], "out of range: a=%e > a_max=%e\n", z, z_table_[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  array_interpolate_spline(const_cast<double*>(z_table_.data()),
                           bt_size_,
                           const_cast<double*>(tau_table_.data()),
                           const_cast<double*>(d2tau_dz2_table_.data()),
                           1,
                           z,
                           &last_index,
                           tau,
                           1);
}

/**
 * Background quantities at given \f$ a \f$.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the beginning of this file). In extended
 * cosmological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pvecback_B    Input: vector containing all {B} type quantities (scale factor, ...)
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assumed to be already allocated)
 * @return the error status
 */

void BackgroundModule::background_functions(
    double* pvecback_B, /* Vector containing all {B} quantities. */
    short return_format,
    double*
        pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
) {
  /** Summary: */

  /** - initialize local variables */
  double a       = pvecback_B[index_bi_a_];
  double rho_tot = 0.;
  double p_tot   = 0.;
  double rho_r   = 0.;
  double rho_m   = 0.;

  class_test(a <= 0., "a = %e instead of strictly positive", a);

  /** - pass value of \f$ a\f$ to output */
  pvecback[index_bg_a_] = a;

  /** - compute each component's density and pressure */

  /* Helper: accumulate rho/p/dp into totals and into rho_r or rho_m based on species type. */
  auto accumulate = [&](const BaseSpecies& sp) {
    const double rho  = sp.Rho(pvecback);
    const double p    = sp.P(pvecback);
    rho_tot          += rho;
    p_tot            += p;
    switch (sp.energy_type()) {
      case BaseSpecies::EnergyType::Radiation:
        rho_r += rho;
        break;
      case BaseSpecies::EnergyType::Matter:
        rho_m += rho;
        break;
      case BaseSpecies::EnergyType::Other:
        rho_r += 3. * p;
        rho_m += rho - 3. * p;
        break;
      default:
        break; /* DarkEnergy: no rho_r / rho_m contribution */
    }
  };

  /* Compute background for all species. Each species writes its own pvecback
     slots; the Fluid computes its own w(a) in ComputeBackground. */
  for (const auto& [name, sp] : all_species_) {
    sp->ComputeBackground(a, pvecback_B, pvecback);
    accumulate(*sp);
  }

  /** - compute expansion rate H from Friedmann equation: this is the
      only place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of \f$ [3c^2/8\pi G] \f$, ie
      \f$ \rho_{class} = [8 \pi G \rho_{physical} / 3 c^2]\f$ */
  pvecback[index_bg_H_] = sqrt(rho_tot - pba->K / a / a);
  const double H        = pvecback[index_bg_H_];

  /** - compute derivative of H with respect to conformal time */
  pvecback[index_bg_H_prime_] = -3. / 2. * (rho_tot + p_tot) * a + pba->K / a;

  /* Total energy density*/
  pvecback[index_bg_rho_tot_] = rho_tot;

  /* Total pressure */
  pvecback[index_bg_p_tot_] = p_tot;

  /* Derivative of total pressure w.r.t. conformal time: accumulated post-H so
     each species applies a' / a = a*H itself (the scalar field's p' is not of
     the a*H*dp/dlna form). FinalizeBackground writes any H-dependent owned slot
     (e.g. scf p_prime_scf). */
  double p_tot_prime = 0.;
  for (const auto& [name, sp] : all_species_) {
    p_tot_prime += sp->PPrime(a, H, pvecback_B, pvecback);
    sp->FinalizeBackground(a, H, pvecback_B, pvecback);
  }
  pvecback[index_bg_p_tot_prime_] = p_tot_prime;

  /** - compute critical density */
  double rho_crit = rho_tot - pba->K / a / a;
  class_test(rho_crit <= 0., "rho_crit = %e instead of strictly positive", rho_crit);

  /** - compute relativistic density to total density ratio */
  pvecback[index_bg_Omega_r_] = rho_r / rho_crit;

  /** - compute other quantities in the exhaustive, redundant format */
  if (return_format == pba->long_info) {
    /** - store critical density */
    pvecback[index_bg_rho_crit_] = rho_crit;

    /** - compute Omega_m */
    pvecback[index_bg_Omega_m_] = rho_m / rho_crit;

    /* one can put other variables here */
    /*  */
    /*  */
  }
}

/**
 * Single place where the fluid equation of state is
 * defined. Parameters of the function are passed through the
 * background structure. Generalisation to arbitrary functions should
 * be simple.
 *
 * @param a              Input: current value of scale factor
 * @param w_fld          Output: equation of state parameter w_fld(a)
 * @param dw_over_da_fld Output: function dw_fld/da
 * @param integral_fld   Output: function \f$ \int_{a}^{a_0} da 3(1+w_{fld})/a \f$
 * @return the error status
 */

int BackgroundModule::background_w_fld(double a,
                                       double* w_fld,
                                       double* dw_over_da_fld,
                                       double* integral_fld) const {
  if (all_species_.count("Fluid")) {
    return static_cast<FluidSpecies&>(*all_species_.at("Fluid"))
        .ComputeWFld(a, w_fld, dw_over_da_fld, integral_fld);
  }
  /* No Fluid species: fall back to CLP defaults (w0=-1, wa=0) so callers
     (HyRec, etc.) get a sensible w(a) without needing to gate on has_fld. */
  *w_fld          = -1.;
  *dw_over_da_fld = 0.;
  *integral_fld   = 0.;
  return _SUCCESS_;
}

void BackgroundModule::background_idm_drmd(
    double a, double rho_idm_over_rho_idr, double* Rint, double* csp2, double* Gint) const {
  static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"))
      .ComputeIdmDrmd(a, rho_idm_over_rho_idr, Rint, csp2, Gint);
}

double BackgroundModule::f_idr_drmd() const {
  if (!all_species_.count("IDM_DRMD_IDR_DRMD"))
    return 0.;
  return static_cast<const IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"))
      .f_idr_drmd();
}

double BackgroundModule::z_dec_drmd() const {
  if (!all_species_.count("IDM_DRMD_IDR_DRMD"))
    return -1.;
  return static_cast<const IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"))
      .z_dec_drmd();
}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @return the error status
 */

void BackgroundModule::background_init() {
  /** Summary: */

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n", _VERSION_);
    printf("Computing background\n");

    /* below we want to inform the user about ncdm species and/or the total N_eff */
    const double Omega0_idr = all_species_.count("IDM_DR_IDR")
                                  ? static_cast<const IDM_DR_IDR_Species&>(
                                        *all_species_.at("IDM_DR_IDR"))
                                        .idr()
                                        .GetOmega0()
                                  : 0.;
    const double Omega0_ur  = all_species_.count("UR") ? all_species_.at("UR")->GetOmega0() : 0.;
    if (!GetNcdmSpecies(all_species_).empty() || (Omega0_idr != 0.)) {
      /* contribution of ultra-relativistic species _ur to N_eff */
      double Neff = Omega0_ur / 7. * 8. / pow(4. / 11., 4. / 3.) / pba->Omega0_g;

      /* contribution of ncdm species to N_eff*/
      if (!GetNcdmSpecies(all_species_).empty()) {
        for (auto& sp : all_species_) {
          Neff += sp->NeffContribution(0.);
          sp->PrintNeffInfo();
        }
      }

      /* contribution of interacting dark radiation _idr to N_eff */
      if (Omega0_idr != 0.) {
        double N_dark  = Omega0_idr / 7. * 8. / pow(4. / 11., 4. / 3.) / pba->Omega0_g;
        Neff          += N_dark;
        printf(" -> dark radiation Delta Neff %e\n", N_dark);
      }

      printf(
          " -> total N_eff = %g (sumed over ultra-relativistic species, ncdm and dark radiation)\n",
          Neff);
    }
  }

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  background_indices();

  /* fluid equation of state */
  if (all_species_.count("Fluid")) {
    double w_fld, dw_over_da, integral_fld;
    background_w_fld(0., &w_fld, &dw_over_da, &integral_fld);

    class_test(w_fld >= 1. / 3.,
               "Your choice for w(a--->0)=%g is suspicious, since it is bigger than -1/3 there "
               "cannot be radiation domination at early times\n",
               w_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (!GetNcdmSpecies(all_species_).empty())) {
    for (auto& sp : all_species_)
      sp->PrintMassInfo();
  }

  class_test(_Gyr_over_Mpc_ <= 0,
             "_Gyr_over_Mpc = %e instead of strictly positive",
             _Gyr_over_Mpc_);

  /** - this function integrates the background over time, allocates
      and fills the background table */
  background_solve_evolver();

  /** - this function finds and stores a few derived parameters at radiation-matter equality */
  background_find_equality();

  background_output_budget();
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @return the error status
 */

void BackgroundModule::background_indices() {
  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all indices */

  index_bg = 0;

  /* index for scale factor */
  class_define_index(index_bg_a_, _TRUE_, index_bg, 1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(index_bg_H_, _TRUE_, index_bg, 1);
  class_define_index(index_bg_H_prime_, _TRUE_, index_bg, 1);

  /* - end of indices in the short vector of background values */
  bg_size_short_ = index_bg;

  // ── Every species registers its own background indices (flat loop) ──────────
  for (auto& [name, sp] : all_species_)
    sp->RegisterBackgroundIndices(index_bg);

  // ── Module aggregate indices ──────────────────────────────────────────────
  /* - index for total density */
  class_define_index(index_bg_rho_tot_, _TRUE_, index_bg, 1);

  /* - index for total pressure */
  class_define_index(index_bg_p_tot_, _TRUE_, index_bg, 1);

  /* - index for derivative of total pressure */
  class_define_index(index_bg_p_tot_prime_, _TRUE_, index_bg, 1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(index_bg_Omega_r_, _TRUE_, index_bg, 1);

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */

  /* - end of indices in the normal vector of background values */
  bg_size_normal_ = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(index_bg_rho_crit_, _TRUE_, index_bg, 1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(index_bg_Omega_m_, _TRUE_, index_bg, 1);

  /* -> conformal distance */
  class_define_index(index_bg_conf_distance_, _TRUE_, index_bg, 1);

  /* -> angular diameter distance */
  class_define_index(index_bg_ang_distance_, _TRUE_, index_bg, 1);

  /* -> luminosity distance */
  class_define_index(index_bg_lum_distance_, _TRUE_, index_bg, 1);

  /* -> proper time (for age of the Universe) */
  class_define_index(index_bg_time_, _TRUE_, index_bg, 1);

  /* -> conformal sound horizon */
  class_define_index(index_bg_rs_, _TRUE_, index_bg, 1);

  /* -> density growth factor in dust universe */
  class_define_index(index_bg_D_, _TRUE_, index_bg, 1);

  /* -> velocity growth factor in dust universe */
  class_define_index(index_bg_f_, _TRUE_, index_bg, 1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  bg_size_ = index_bg;

  /* - now, indices in vector of variables to integrate.
     First {B} variables, then {C} variables. */

  index_bi = 0;

  /* -> scale factor */
  class_define_index(index_bi_a_, _TRUE_, index_bi, 1);

  /* -> integration indices for all species (each species owns its own offsets) */
  for (auto& [name, sp] : all_species_) {
    sp->RegisterIntegrationIndices(index_bi);
  }

  /* End of {B} variables, now continue with {C} variables */
  bi_B_size_ = index_bi;

  /* -> proper time (for age of the Universe) */
  class_define_index(index_bi_time_, _TRUE_, index_bi, 1);

  /* -> sound horizon */
  class_define_index(index_bi_rs_, _TRUE_, index_bi, 1);

  /* -> Second order equation for growth factor */
  class_define_index(index_bi_D_, _TRUE_, index_bi, 1);
  class_define_index(index_bi_D_prime_, _TRUE_, index_bi, 1);

  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(index_bi_tau_, _TRUE_, index_bi, 1);

  /* -> end of indices in the vector of variables to integrate */
  bi_size_ = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(index_bi_tau_ != index_bi - 1,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  /* Set BackgroundModule pointer on all active species (default is no-op) */
  for (const auto& [name, sp] : all_species_)
    sp->SetBackgroundModule(this);
}

void BackgroundModule::background_solve_evolver() {
  /** Summary: */

  /** - define local variables */

  /* parameters and workspace for the background_derivs_loga / background_derivs evolver functions */
  struct background_parameters_and_workspace bpaw{this};
  /* vector of quantities to be integrated */
  std::vector<double> pvecback_integration(bi_size_);
  /* vector of all background quantities */
  std::vector<double> pvecback(bg_size_);

  bpaw.pvecback = pvecback.data();

  /** - impose initial conditions with background_initial_conditions() */
  background_initial_conditions(pvecback.data(), pvecback_integration.data());

  /** - Determine output vector */
  double loga_ini   = log(pvecback_integration[index_bi_a_]);
  double loga_final = 0.;
  bt_size_          = (loga_final - loga_ini) / ppr->back_integration_stepsize;
  std::vector<double> loga(bt_size_);
  std::vector<int> used_in_output(bt_size_);
  for (int index_loga = 0; index_loga < bt_size_; index_loga++) {
    loga[index_loga]           = loga_ini + index_loga * (loga_final - loga_ini) / (bt_size_ - 1);
    used_in_output[index_loga] = 1;
  }
  // -ffast-math robustness fix: preserve exact boundary values.
  loga.front() = loga_ini;
  loga.back()  = loga_final;

  /** - Remember that we evolve tau at index_bi_a: */
  pvecback_integration[index_bi_a_] = pvecback_integration[index_bi_tau_];

  /** - allocate background tables */
  tau_table_.resize(bt_size_);
  z_table_.resize(bt_size_);
  d2tau_dz2_table_.resize(bt_size_);
  background_table_.resize(bt_size_ * bg_size_);
  d2background_dtau2_table_.resize(bt_size_ * bg_size_);

  auto generic_evolver = &evolver_ndf15;
  if (ppr->evolver == evolver_type::rk) {
    generic_evolver = &evolver_rk;
  }
  else if (ppr->evolver == evolver_type::rkdp45) {
    generic_evolver = &evolver_rkdp45;
  }

  /* Size of vector to integrate is (bi_size_-1) rather than
   * (bi_size_), since a is not integrated.
   */
  generic_evolver(background_derivs_loga,
                  loga_ini,
                  loga_final,
                  pvecback_integration.data(),
                  used_in_output.data(),
                  bi_size_ - 1,
                  &bpaw,
                  ppr->tol_background_integration,
                  ppr->smallest_allowed_variation,
                  background_timescale,
                  ppr->perturb_integration_stepsize,
                  loga.data(),
                  bt_size_,
                  background_add_line_to_bg_table,
                  nullptr);

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  age_ = pvecback_integration[index_bi_time_] / _Gyr_over_Mpc_;
  /* -> conformal age in Mpc. Remember that tau is stored at index_bi_a now */
  conformal_age_ = pvecback_integration[index_bi_a_];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  Omega0_dr_ = 0.;
  for (auto& sp : all_species_)
    Omega0_dr_ += sp->DarkRadiationRhoToday(pvecback_integration.data()) / pba->H0 / pba->H0;
  if (all_species_.count("DCDM_DR")) {
    auto& dcdm_dr = dynamic_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
    Omega0_dcdm_  = pvecback_integration[dcdm_dr.dcdm().bi_rho_index()] / pba->H0 / pba->H0;
  }

  /** Recover some quantities today */
  double D_today = pvecback_integration[index_bi_D_];
  for (int i = 0; i < bt_size_; i++) {
    double* bg_table_row = background_table_.data() + i * bg_size_;
    /** Set cosmological distances */
    double conformal_distance             = conformal_age_ - tau_table_[i];
    bg_table_row[index_bg_conf_distance_] = conformal_distance;
    double comoving_radius                = conformal_distance;
    if (pba->sgnK > 0) {
      comoving_radius = sin(sqrt(pba->K) * conformal_distance) / sqrt(pba->K);
    }
    else if (pba->sgnK < 0) {
      comoving_radius = sinh(sqrt(-pba->K) * conformal_distance) / sqrt(-pba->K);
    }

    bg_table_row[index_bg_ang_distance_] = comoving_radius / (1. + z_table_[i]);
    bg_table_row[index_bg_lum_distance_] = comoving_radius * (1. + z_table_[i]);
    /** Normalise D(z=0)=1 */
    bg_table_row[index_bg_D_] /= D_today;
  }

  /* Let species run table-scope analysis now that the full background table is
     filled (e.g. DRMD decoupling redshift). */
  for (auto& [name, sp] : all_species_)
    sp->ProcessBackgroundTable(background_table_.data(), bt_size_, bg_size_, z_table_.data());

  /** - fill tables of second derivatives (in view of spline interpolation) */
  array_spline_table_lines(z_table_.data(),
                           bt_size_,
                           tau_table_.data(),
                           1,
                           d2tau_dz2_table_.data(),
                           _SPLINE_EST_DERIV_);

  array_spline_table_lines(tau_table_.data(),
                           bt_size_,
                           background_table_.data(),
                           bg_size_,
                           d2background_dtau2_table_.data(),
                           _SPLINE_EST_DERIV_);

  /** - compute remaining "related parameters"
   *     - so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */
  {
    const double* earliest      = background_table_.data();
    const double rho_g_earliest = all_species_.photons().Rho(earliest);
    Neff_ = (background_table_[index_bg_Omega_r_] * background_table_[index_bg_rho_crit_] -
             rho_g_earliest) /
            (7. / 8. * pow(4. / 11., 4. / 3.) * rho_g_earliest);
  }

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n", age_);
    printf(" -> conformal age = %f Mpc\n", conformal_age_);
  }

  if (pba->background_verbose > 2) {
    if (all_species_.count("DCDM_DR")) {
      const auto& dcdm_dr_comp2 = static_cast<const DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n", Omega0_dcdm_);
      printf("     -> Omega0_dr = %f\n", Omega0_dr_);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
             Omega0_dr_ + Omega0_dcdm_,
             dcdm_dr_comp2.dcdm().GetOmega0());
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n",
             dcdm_dr_comp2.dcdm().Omega_ini_dcdm() / pba->Omega0_b);
    }
    if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
      auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
      printf(" -> Dark Radiation Matter Decoupling details: (DRMD)\n");
      printf(
          "     -> values: (initial) Gamma0 = %f 1/Mpc, zstop= %e,f_idr_drmd=%e, and f_idm= %e \n",
          drmd.Gamma0_drmd_ic(),
          drmd.z_stop(),
          drmd.f_idr_drmd(),
          drmd.f_idm_drmd());
      printf("     -> dark radiation Delta N_eff (DRMD) %e\n", drmd.delta_Neff_drmd());

      if (drmd.z_dec_drmd() > 0)
        printf("     -> decoupling occurred at z=%f \n", drmd.z_dec_drmd());
      else
        printf("     -> no decoupling occurred.\n");
    }
    if (all_species_.count("ScalarField")) {
      const auto& scf = static_cast<const ScalarFieldSpecies&>(*all_species_.at("ScalarField"));
      printf("    Scalar field details:\n");
      printf("     -> Omega_scf = %g, wished %g\n",
             scf.Rho(pvecback.data()) / pvecback[index_bg_rho_crit_],
             scf.GetOmega0());
      if (all_species_.count("Lambda"))
        printf("     -> Omega_Lambda = %g, wished %g\n",
               all_species_.at("Lambda")->Rho(pvecback.data()) / pvecback[index_bg_rho_crit_],
               all_species_.at("Lambda")->GetOmega0());
      printf("     -> parameters: [lambda, alpha, A, B] = \n");
      printf("                    [");
      const auto& scf_params = scf.scf_parameters();
      for (size_t i = 0; i < scf_params.size() - 1; i++) {
        printf("%.3f, ", scf_params[i]);
      }
      printf("%.3f]\n", scf_params[scf_params.size() - 1]);
    }
  }

  /**  - total matter, radiation, dark energy today */
  Omega0_m_  = background_table_[(bt_size_ - 1) * bg_size_ + index_bg_Omega_m_];
  Omega0_r_  = background_table_[(bt_size_ - 1) * bg_size_ + index_bg_Omega_r_];
  Omega0_de_ = 1. - (Omega0_m_ + Omega0_r_ + pba->Omega0_k);

  /* -> Backfill each DNCDM matter child's today density fraction (parallels Omega0_dcdm_).
     In combined/initial modes the child only carries `deg`; without this its GetOmega0()
     stays 0, dropping it from fnu (GetOmega0NcdmTot) and the budget-print neutrino line.
     The DNCDM density is not an evolved bi_ variable, so read it from the today table row. */
  const double* bg_today = background_table_.data() + (bt_size_ - 1) * bg_size_;
  for (auto& [key, sp] : all_species_) {
    if (auto* dncdm_dr = dynamic_cast<DNCDM_DR_Species*>(sp.get()))
      dncdm_dr->dncdm().BackfillOmega0FromToday(bg_today, pba->H0, pba->h);
  }
}

/**
 * Assign initial values to background integrated variables.
 *
 * @param pvecback             Input: vector of background quantities used as workspace
 * @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
 * @return the error status
 */

void BackgroundModule::background_initial_conditions(
    double*
        pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
    double*
        pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size bi_size_) */
) {
  /** Summary: */

  /** - fix initial value of \f$ a \f$ */
  double a = ppr->a_ini_over_a_today_default;

  /** Allow any species to pull the integration start earlier than the default
      (e.g. NCDM must be relativistic at a_ini; relevant for some WDM models). */

  for (auto& [name, sp] : all_species_)
    a = sp->BackgroundAIni(a, ppr->tol_ncdm_initial_w);

  pvecback_integration[index_bi_a_] = a;

  /* Set initial values of {B} variables: */
  double Omega_rad = pba->Omega0_g;
  for (auto& [name, sp] : all_species_)
    Omega_rad += sp->GetRadiationOmega0();
  double rho_rad = Omega_rad * pow(pba->H0, 2) / pow(a, 4);
  /* Set initial conditions for all species background ODE variables.  Each
     species owns its own integration offsets and IC math (DCDM, DNCDM, Fluid,
     ScalarField, ...); the module only supplies the shared context. */
  BackgroundICContext ic;
  ic.a_ini                = a;
  ic.rho_rad              = rho_rad;
  ic.pvecback_integration = pvecback_integration;
  for (auto& [name, sp] : all_species_) {
    sp->SetBackgroundInitialConditions(ic);
  }

  /* Infer pvecback from pvecback_integration */
  background_functions(pvecback_integration, pba->normal_info, pvecback);

  /* Just checking that our initial time indeed is deep enough in the radiation
     dominated regime */
  class_test(fabs(pvecback[index_bg_Omega_r_] - 1.) > ppr->tol_initial_Omega_r,
             "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to "
             "start from radiation domination.",
             pvecback[index_bg_Omega_r_]);

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[index_bg_H_] <= 0.,
             "H = %e instead of strictly positive",
             pvecback[index_bg_H_]);

  if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
    static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"))
        .InitializeDrmdBackground(pvecback[index_bg_rho_tot_], pvecback[index_bg_H_], a, pvecback);
  }

  pvecback_integration[index_bi_time_] = 1. / (2. * pvecback[index_bg_H_]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[index_bi_tau_] = 1. / (a * pvecback[index_bg_H_]);

  /** - compute initial sound horizon, assuming \f$ c_s=1/\sqrt{3} \f$ initially */
  pvecback_integration[index_bi_rs_] = pvecback_integration[index_bi_tau_] / sqrt(3.);

  /** - set initial value of D and D' in RD. D will be renormalised later, but D' must be correct. */
  pvecback_integration[index_bi_D_]       = 1;
  pvecback_integration[index_bi_D_prime_] = 2 * a * pvecback[index_bg_H_];
}

/**
 * Find the time of radiation/matter equality and store characteristic
 * quantitites at that time in the background structure..
 *
 * @return the error status
 */

void BackgroundModule::background_find_equality() {
  double Omega_m_over_Omega_r = 0.;

  /* first bracket the right tau value between two consecutive indices in the table */

  int index_tau_plus  = bisect_index(0, bt_size_ - 1, [&](int i) {
    return background_table_[i * bg_size_ + index_bg_Omega_m_] /
               background_table_[i * bg_size_ + index_bg_Omega_r_] >
           1.;
  });
  int index_tau_minus = index_tau_plus - 1;

  /* then get a better estimate within this range */

  double tau_minus = tau_table_[index_tau_minus];
  double tau_plus  = tau_table_[index_tau_plus];
  double tau_mid   = 0.;

  std::vector<double> pvecback(bg_size_);

  /* Refinement bisection kept inline: a_eq_/H_eq_/tau_eq_ must stay consistent
     with pvecback at the final evaluated tau_mid, which bisect_value's post-loop
     midpoint return would not preserve. */
  while ((tau_plus - tau_minus) > ppr->tol_tau_eq) {
    tau_mid = 0.5 * (tau_plus + tau_minus);

    background_at_tau(tau_mid,
                      pba->long_info,
                      pba->inter_closeby,
                      &index_tau_minus,
                      pvecback.data());

    Omega_m_over_Omega_r = pvecback[index_bg_Omega_m_] / pvecback[index_bg_Omega_r_];

    if (Omega_m_over_Omega_r > 1)
      tau_plus = tau_mid;
    else
      tau_minus = tau_mid;
  }

  a_eq_   = pvecback[index_bg_a_];
  H_eq_   = pvecback[index_bg_H_];
  z_eq_   = 1. / a_eq_ - 1.;
  tau_eq_ = tau_mid;

  if (pba->background_verbose > 0) {
    printf(" -> radiation/matter equality at z = %f\n", z_eq_);
    printf("    corresponding to conformal time = %f Mpc\n", tau_eq_);
  }
}

/**
 * Subroutine for formatting background output
 *
 */

void BackgroundModule::background_output_titles(std::string& titles) const {
  // ── Module header (always present) ──────────────────────────────────────
  class_store_columntitle(titles, "z", _TRUE_);
  class_store_columntitle(titles, "proper time [Gyr]", _TRUE_);
  class_store_columntitle(titles, "conf. time [Mpc]", _TRUE_);
  class_store_columntitle(titles, "H [1/Mpc]", _TRUE_);
  class_store_columntitle(titles, "comov. dist.", _TRUE_);
  class_store_columntitle(titles, "ang.diam.dist.", _TRUE_);
  class_store_columntitle(titles, "lum. dist.", _TRUE_);
  class_store_columntitle(titles, "comov.snd.hrz.", _TRUE_);

  // ── Species output — per-species dispatch ───────────────────────────────
  BackgroundColumnWriter writer(titles);
  for (auto& [name, sp] : all_species_)
    sp->WriteBackgroundColumnTitles(writer);

  // ── Module aggregate columns ────────────────────────────────────────────
  class_store_columntitle(titles, "(.)rho_crit", _TRUE_);
  class_store_columntitle(titles, "(.)rho_tot", _TRUE_);
  class_store_columntitle(titles, "(.)p_tot", _TRUE_);
  class_store_columntitle(titles, "(.)p_tot_prime", _TRUE_);
  class_store_columntitle(titles, "gr.fac. D", _TRUE_);
  class_store_columntitle(titles, "gr.fac. f", _TRUE_);
}

void BackgroundModule::background_output_data(int number_of_titles, double* data) const {
  for (int index_tau = 0; index_tau < bt_size_; index_tau++) {
    double* dataptr  = data + index_tau * number_of_titles;
    double* pvecback = const_cast<double*>(background_table_.data()) + index_tau * bg_size_;
    int storeidx     = 0;

    // ── Module header ──────────────────────────────────────────────────────
    class_store_double(dataptr, 1. / pvecback[index_bg_a_] - 1., _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_time_] / _Gyr_over_Mpc_, _TRUE_, storeidx);
    class_store_double(dataptr,
                       conformal_age_ - pvecback[index_bg_conf_distance_],
                       _TRUE_,
                       storeidx);
    class_store_double(dataptr, pvecback[index_bg_H_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_conf_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_ang_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_lum_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rs_], _TRUE_, storeidx);

    // ── Species data — per-species dispatch ───────────────────────────────
    BackgroundColumnWriter writer(dataptr, storeidx);
    for (auto& [name, sp] : all_species_)
      sp->WriteBackgroundData(pvecback, writer);

    // ── Module aggregate columns ──────────────────────────────────────────
    class_store_double(dataptr, pvecback[index_bg_rho_crit_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_tot_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_tot_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_tot_prime_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_D_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_f_], _TRUE_, storeidx);
  }
}

/**
 * Computes background ODE derivatives with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * Called via background_derivs_loga_member during evolver integration.
 * The arguments follow the calling convention of the evolver interface:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector.
 *
 * - errors are not written to a module error_message buffer, but to a generic
 * error_message passed in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int BackgroundModule::background_derivs_member(
    double tau,
    double* y, /* vector with argument y[index_bi] (must be already allocated with size bi_size_) */
    double* dy, /* vector with argument dy[index_bi]
                                                            (must be already allocated with
                                                            size bi_size_) */
    void* parameters_and_workspace) {
  /** Summary: */

  /** - define local variables */

  double *pvecback, a, H;

  background_parameters_and_workspace* pbpaw = static_cast<background_parameters_and_workspace*>(
      parameters_and_workspace);
  pvecback = pbpaw->pvecback;

  /** - calculate functions of \f$ a \f$ with background_functions() */
  background_functions(y, pba->normal_info, pvecback);

  /** - Short hand notation */
  a = y[index_bi_a_];
  H = pvecback[index_bg_H_];

  /** - calculate \f$ a'=a^2 H \f$ */
  dy[index_bi_a_] = y[index_bi_a_] * y[index_bi_a_] * pvecback[index_bg_H_];

  /** - calculate \f$ t' = a \f$ */
  dy[index_bi_time_] = y[index_bi_a_];

  class_test(all_species_.photons().Rho(pvecback) <= 0.,
             "rho_g = %e instead of strictly positive",
             all_species_.photons().Rho(pvecback));

  /** - calculate \f$ rs' = c_s \f$*/
  dy[index_bi_rs_] = 1. /
                     sqrt(3. * (1. + 3. * all_species_.baryons().Rho(pvecback) / 4. /
                                         all_species_.photons().Rho(pvecback))) *
                     sqrt(1. -
                          pba->K * y[index_bi_rs_] * y[index_bi_rs_]);  // TBC: curvature correction

  /** - solve second order growth equation  \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */
  double rho_M = all_species_.baryons().Rho(pvecback);
  if (all_species_.count("CDM"))
    rho_M += all_species_.at("CDM")->Rho(pvecback);
  if (all_species_.count("IDM_DR_IDR")) {
    auto& idm_idr  = static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
    rho_M         += idm_idr.idm_dr().Rho(pvecback);
  }
  if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
    auto& drmd  = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
    rho_M      += drmd.idm_drmd().Rho(pvecback);
  }
  if (all_species_.count("CDM_SCF_Momentum")) {
    auto& t3  = static_cast<Type3Species&>(*all_species_.at("CDM_SCF_Momentum"));
    rho_M    += t3.cdm().Rho(pvecback);
  }

  dy[index_bi_D_]       = y[index_bi_D_prime_];
  dy[index_bi_D_prime_] = -a * H * y[index_bi_D_prime_] + 1.5 * a * a * rho_M * y[index_bi_D_];

  /* Species background ODE contributions (including DCDM_DR and DNCDM_DR composites). */
  for (const auto& [name, sp] : all_species_) {
    sp->BackgroundDerivs(tau, y, dy, pvecback);  // default is no-op
  }

  return _SUCCESS_;
}

/**
 * Function outputting the fractions Omega of the total critical density
 * today, and also the reduced fractions omega=Omega*h*h
 *
 * It also prints the total budgets of non-relativistic, relativistic,
 * and other contents, and of the total
 *
 * @return the error status
 */

void BackgroundModule::background_output_budget() {
  if (pba->background_verbose > 1) {
    double budget_matter    = 0;
    double budget_radiation = 0;
    double budget_other     = 0;
    double budget_neutrino  = 0;

    // Helper: print one species line and accumulate into a budget bucket.
    auto print_one = [&](const char* label, double omega, double& budget) {
      printf("-> %-30s Omega = %-15g , omega = %-15g\n", label, omega, omega * pba->h * pba->h);
      budget += omega;
    };

    printf(" ---------------------------- Budget equation ----------------------- \n");

    printf(" ---> Nonrelativistic Species \n");
    print_one("Bayrons", pba->Omega0_b, budget_matter);
    if (all_species_.count("CDM"))
      print_one("Cold Dark Matter", all_species_.at("CDM")->GetOmega0(), budget_matter);
    if (all_species_.count("IDM_DR_IDR")) {
      const auto& comp = static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
      print_one("Interacting Dark Matter - DR ", comp.idm_dr().GetOmega0(), budget_matter);
    }
    if (all_species_.count("CDM_SCF_Momentum")) {
      const auto& comp = static_cast<const Type3Species&>(*all_species_.at("CDM_SCF_Momentum"));
      print_one("Cold Dark Matter (Type-3 coupled)", comp.cdm().GetOmega0(), budget_matter);
    }
    if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
      const auto& comp = static_cast<const IDM_DRMD_IDR_DRMD_Species&>(
          *all_species_.at("IDM_DRMD_IDR_DRMD"));
      print_one("Interacting DM (DRMD)", comp.idm_drmd().GetOmega0(), budget_matter);
    }
    if (all_species_.count("DCDM_DR")) {
      // Use integration-derived Omega0_dcdm_ for accuracy (set in background_solve_evolver).
      print_one("Decaying Cold Dark Matter", Omega0_dcdm_, budget_matter);
    }

    printf(" ---> Relativistic Species \n");
    print_one("Photons", pba->Omega0_g, budget_radiation);
    if (all_species_.count("UR"))
      print_one("Ultra-relativistic relics", all_species_.at("UR")->GetOmega0(), budget_radiation);
    if (GetNDecayDr() > 0) {
      // Omega0_dr_ aggregates the decay radiation of DCDM_DR and every DNCDM_DR
      // composite (integration-derived, set in background_solve_evolver), so this
      // single line covers all decay-DR channels.
      print_one("Dark Radiation (from decay)", Omega0_dr_, budget_radiation);
    }
    if (all_species_.count("IDM_DR_IDR")) {
      const auto& comp = static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
      print_one("Interacting Dark Radiation", comp.idr().GetOmega0(), budget_radiation);
    }
    if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
      const auto& comp = static_cast<const IDM_DRMD_IDR_DRMD_Species&>(
          *all_species_.at("IDM_DRMD_IDR_DRMD"));
      print_one("Dark Radiation (DRMD)", comp.idr_drmd().GetOmega0(), budget_radiation);
    }

    if (!GetNcdmSpecies(all_species_).empty()) {
      printf(" ---> Massive Neutrino Species \n");
      for (auto* sp : GetNcdmSpecies(all_species_)) {
        sp->PrintOmegaInfo();
        budget_neutrino += sp->GetOmega0();
      }
    }

    if (all_species_.count("Lambda") || all_species_.count("Fluid") ||
        all_species_.count("ScalarField") || pba->sgnK != 0) {
      printf(" ---> Other Content \n");
    }
    if (all_species_.count("Lambda"))
      print_one("Cosmological Constant", all_species_.at("Lambda")->GetOmega0(), budget_other);
    if (all_species_.count("Fluid"))
      print_one("Dark Energy Fluid", all_species_.at("Fluid")->GetOmega0(), budget_other);
    if (all_species_.count("ScalarField"))
      print_one("Scalar Field", all_species_.at("ScalarField")->GetOmega0(), budget_other);
    if (pba->sgnK != 0)
      print_one("Spatial Curvature", pba->Omega0_k, budget_other);

    printf(" ---> Total budgets \n");
    printf(" Radiation                        Omega = %-15g , omega = %-15g \n",
           budget_radiation,
           budget_radiation * pba->h * pba->h);
    printf(" Non-relativistic                 Omega = %-15g , omega = %-15g \n",
           budget_matter,
           budget_matter * pba->h * pba->h);
    if (!GetNcdmSpecies(all_species_).empty()) {
      printf(" Neutrinos                        Omega = %-15g , omega = %-15g \n",
             budget_neutrino,
             budget_neutrino * pba->h * pba->h);
    }
    if (all_species_.count("Lambda") || all_species_.count("Fluid") ||
        all_species_.count("ScalarField") || pba->sgnK != 0) {
      printf(" Other Content                    Omega = %-15g , omega = %-15g \n",
             budget_other,
             budget_other * pba->h * pba->h);
    }
    printf(" TOTAL                            Omega = %-15g , omega = %-15g \n",
           budget_radiation + budget_matter + budget_neutrino + budget_other,
           (budget_radiation + budget_matter + budget_neutrino + budget_other) * pba->h * pba->h);

    printf(" -------------------------------------------------------------------- \n");
  }
}

/**
 * Derivative function passed to the evolver (generic_evolver / evolver_ndf15 /
 * evolver_rk) evaluating the derivative with respect to log(a)
 * of quantities which are integrated (tau, t, etc).
 *
 * The arguments follow the calling convention of the evolver interface:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector.
 *
 * - errors are not written to a module buffer, but to a generic
 * error_message passed in the list of arguments.
 *
 * @param loga                        Input: scale factor
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int BackgroundModule::background_derivs_loga_member(
    double loga,
    double* y, /* vector with argument y[index_bi] (must be already allocated with size bi_size_) */
    double* dy, /* vector with argument dy[index_bi]
                                                                 (must be already allocated with
                                                                 size bi_size_) */
    void* parameters_and_workspace) {
  background_parameters_and_workspace* pbpaw = static_cast<background_parameters_and_workspace*>(
      parameters_and_workspace);
  double* pvecback = pbpaw->pvecback;

  /** Note that we want to reuse as much as possible of the usual tau
      integration, so inside y, index_bi_a is really the index of tau. */
  double a       = exp(loga);
  double tau     = y[index_bi_a_];
  y[index_bi_a_] = a;

  /** Get derivatives w.r.t. conformal time */
  background_derivs_member(tau, y, dy, parameters_and_workspace);

  /** Swap a and tau again */
  y[index_bi_a_] = tau;
  /** Set dtau/dloga: */
  dy[index_bi_a_] = 1.0;

  /** Convert to scale factor derivatives */
  double H = pvecback[index_bg_H_];
  for (int index_bi = 0; index_bi < bi_size_ - 1; index_bi++) {
    dy[index_bi] *= 1. / (a * H);
  }
  return _SUCCESS_;
}

int BackgroundModule::background_add_line_to_bg_table_member(
    double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace) {
  double a       = exp(loga);
  double tau     = y[index_bi_a_];
  y[index_bi_a_] = a;

  z_table_[index_loga]   = std::max(0., 1. / a - 1.);
  tau_table_[index_loga] = tau;

  double* pvecback = background_table_.data() + index_loga * bg_size_;

  // compute quantities depending only on {B} variables.
  background_functions(y, pba->long_info, pvecback);

  pvecback[index_bg_time_] = y[index_bi_time_];
  pvecback[index_bg_rs_]   = y[index_bi_rs_];
  pvecback[index_bg_D_]    = y[index_bi_D_];
  pvecback[index_bg_f_]    = y[index_bi_D_prime_] / (y[index_bi_D_] * a * pvecback[index_bg_H_]);

  //Swap a and tau again
  y[index_bi_a_] = tau;

  return _SUCCESS_;
}

int BackgroundModule::background_print_variables(double loga,
                                                 double* y,
                                                 double* dy,
                                                 void* parameters_and_workspace) {
  background_parameters_and_workspace* pbpaw = static_cast<background_parameters_and_workspace*>(
      parameters_and_workspace);
  double* pvecback     = pbpaw->pvecback;
  BackgroundModule& bm = *(pbpaw->background_module);

  double a          = exp(loga);
  double tau        = y[bm.index_bi_a_];
  y[bm.index_bi_a_] = a;

  /** - calculate functions of \f$ a \f$ with background_functions() */
  bm.background_functions(y, bm.pba->normal_info, pvecback);

  /** Swap a and tau again */
  y[bm.index_bi_a_] = tau;
  //FILE* fid = fopen("tmp.dat", "a");
  //fprintf(fid, "%.3e %.3e %.3e %.3e %.3e %.3e\n", exp(loga), tau, pvecback[bm.index_bg_rho_ncdm1_], pvecback[bm.index_bg_rho_dr_species_], y[bm.index_bi_rho_dr_species_], pvecback[bm.index_bg_lnf_ncdm_decay_dr1_ + 2]);
  //fclose(fid);
  return _SUCCESS_;
}

int BackgroundModule::GetNcdmCount() const {
  return static_cast<int>(GetNcdmSpecies(all_species_).size());
}

double BackgroundModule::GetOmega0NcdmTot() const {
  double total = 0.;
  for (auto& sp : all_species_)
    total += sp->NeutrinoOmega0();
  return total;
}

int BackgroundModule::GetNDecayDr() const {
  int n = (all_species_.count("DCDM_DR") ? 1 : 0);
  for (const auto& sp : all_species_) {
    if (dynamic_cast<DNCDM_DR_Species*>(sp.get()))
      ++n;
  }
  return n;
}
