/** @file background.h Documented includes for background module */

#ifndef __BACKGROUND__
#define __BACKGROUND__

#include <memory>
#include <vector>

#include "arrays.h"
#include "common.h"
#include "dei_rkck.h"
#include "growTable.h"
#include "parser.h"
#include "quadrature.h"

class NonColdDarkMatter;
class DarkRadiation;

/** list of possible parametrisations of the DE equation of state */

enum equation_of_state { CLP, EDE };

/**
 * All background parameters and evolution that other modules need to know.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (except thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */

/** Which of {Lambda, Fluid, ScalarField} (if any) is the budget-closure
 *  species — inferred during input parsing from which corresponding Omega
 *  parameter was left unspecified (Omega_Lambda or Omega_fld) or set to the
 *  negative-sentinel value (Omega_scf < 0), and consumed by
 *  InputModule::ConstructSpecies after all non-closure species are built.
 *  Defined at file scope (not nested in `background`) so that
 *  generate_wrapper.py's struct parser doesn't see the enum's `};` and
 *  prematurely terminate the background struct. */
enum class ClosureSpecies { None, Lambda, Fluid, ScalarField };

struct background {
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' structure)
   *
   * The background cosmological parameters listed here form a parameter
   * basis which is directly usable by the background module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_init()
   * e.g. for dark energy or decaying dark matter). */

  //@{

  double H0; /**< \f$ H_0 \f$: Hubble parameter (in fact, [\f$H_0/c\f$]) in \f$ Mpc^{-1} \f$ */

  double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$: photons */

  double T_cmb = 2.7255; /**< \f$ T_{cmb} \f$: current CMB temperature in Kelvins */

  double Omega0_b; /**< \f$ \Omega_{0 b} \f$: baryons */

  // Per-species Omega0_* fields removed: each species owns its density
  // through its constructor (queryable via species->GetOmega0()).  Only
  // photons (Omega0_g), baryons (Omega0_b), and curvature (Omega0_k) live
  // on pba.

  // Fluid physics params (fluid_equation_of_state, w0_fld, wa_fld, cs2_fld,
  // Omega_EDE, use_ppf, c_gamma_over_c_fld) live on FluidSpecies — see
  // FluidSpecies::CreateAll for the parser and species/fluid.h for the
  // accessors used by cross-module readers.

  // DCDM physics params (Gamma_dcdm, Omega_ini_dcdm) live on DCDMSpecies — see
  // DCDM_DR_Species::CreateAll for the parser and species/dcdm.h for the
  // accessors used by cross-module readers.
  // IDR physics params (T_idr, l_max_idr) live on IDRSpecies — see
  // IDM_DR_IDR_Species::CreateAll for the parser and species/idr.h for the
  // accessors used by cross-module readers.

  double Omega0_k = 0.; /**< \f$ \Omega_{0_k} \f$: curvature contribution */

  /** @name - related parameters */

  //@{

  double h = 0.67556; /**< reduced Hubble parameter */
  double K = 0.;      /**< \f$ K \f$: Curvature parameter \f$ K=-\Omega0_k*a_{today}^2*H_0^2\f$; */
  int sgnK = 0;       /**< K/|K|: -1, 0 or 1 */

  //@}

  /** @name - other background parameters */

  //@{

  double a_today = 1.; /**< scale factor today (arbitrary and irrelevant for most purposes) */

  //@}

  ClosureSpecies closure_species = ClosureSpecies::None;

  /**
   *@name - some flags needed for calling background functions
   */

  //@{

  short short_info;  /**< flag for calling background_at_eta and return little information */
  short normal_info; /**< flag for calling background_at_eta and return medium information */
  short long_info;   /**< flag for calling background_at_eta and return all information */

  short
      inter_normal; /**< flag for calling background_at_eta and find position in interpolation table normally */
  short
      inter_closeby; /**< flag for calling background_at_eta and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  unsigned int number_of_threads;

  short background_verbose =
      0; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}
};

/**
 * @name Some limits on possible background parameters
 */

//@{

#define _h_BIG_ 1.5          /**< maximal \f$ h \f$ */
#define _h_SMALL_ 0.3        /**< minimal \f$ h \f$ */
#define _omegab_BIG_ 0.039   /**< maximal \f$ omega_b \f$ */
#define _omegab_SMALL_ 0.005 /**< minimal \f$ omega_b \f$ */

//@}

/**
 * @name Some limits imposed in other parts of the module:
 */

//@{

#define _SCALE_BACK_ \
  0.1 /**< logarithmic step used when searching
			     for an initial scale factor at which ncdm
			     are still relativistic */

//@}

#endif
/* @endcond */
