/** @file common.h Generic libraries, parameters and functions used in the whole code. */

#include <stdarg.h>

#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "svnversion.h"

#ifdef __cplusplus
#include <type_traits>

#include "exceptions.h"
#define typeof(x) std::remove_reference<decltype((x))>::type
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class InputModule;
class BackgroundModule;
class ThermodynamicsModule;
class PerturbationsModule;
class PrimordialModule;
class NonlinearModule;
class TransferModule;
class SpectraModule;
class LensingModule;
class FileContent;

typedef std::shared_ptr<const InputModule> InputModulePtr;
typedef std::shared_ptr<const BackgroundModule> BackgroundModulePtr;
typedef std::shared_ptr<const ThermodynamicsModule> ThermodynamicsModulePtr;
typedef std::shared_ptr<const PerturbationsModule> PerturbationsModulePtr;
typedef std::shared_ptr<const PrimordialModule> PrimordialModulePtr;
typedef std::shared_ptr<const NonlinearModule> NonlinearModulePtr;
typedef std::shared_ptr<const TransferModule> TransferModulePtr;
typedef std::shared_ptr<const SpectraModule> SpectraModulePtr;
typedef std::shared_ptr<const LensingModule> LensingModulePtr;
#endif

#ifndef __COMMON__
#define __COMMON__

#define _VERSION_ "v2.9.0"
/* @cond INCLUDE_WITH_DOXYGEN */

#define _TRUE_ 1  /**< integer associated to true statement */
#define _FALSE_ 0 /**< integer associated to false statement */

#define _SUCCESS_ 0 /**< integer returned after successful call of a function */
#define _FAILURE_ 1 /**< integer returned after failure in a function */

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg
    [_ERRORMSGSIZE_]; /**< Generic error messages (there is such a field in each structure) */

#define _FILENAMESIZE_ \
  256 /**< size of the string read in each line of the file (extra characters not taken into account) */
typedef char FileName[_FILENAMESIZE_];

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _PIHALF_ 1.57079632679489661923132169164e0 /**< pi divided by 2 */

#define _TWOPI_ 6.283185307179586476925286766559e0 /**< 2 times pi */

#define _SQRT2_ 1.41421356237309504880168872421e0 /** < square root of 2. */

#define _SQRT6_ 2.4494897427831780981972840747059e0 /**< square root of 6. */

#define _SQRT_PI_ 1.77245385090551602729816748334e0 /**< square root of pi. */

#define _E_ \
  2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696 /**< exponential of one */

#define _MAX_IT_ \
  10000 /**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define _QUADRATURE_MAX_ \
  250 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _QUADRATURE_MAX_BG_ \
  800 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _TOLVAR_ \
  100. /**< The minimum allowed variation is the machine precision times this number */

#define _HUGE_ 1.e99

#define _EPSILON_ 1.e-10

#define _OUTPUTPRECISION_ 12 /**< Number of significant digits in some output files */

#define _COLUMNWIDTH_ \
  24 /**< Must be at least _OUTPUTPRECISION_+8 for guaranteed fixed width columns */

#define _MAXTITLESTRINGLENGTH_ 8000 /**< Maximum number of characters in title strings */

#define _DELIMITER_ "\t" /**< character used for delimiting titles in the title strings */

#ifndef __CLASSDIR__
#define __CLASSDIR__ \
  "." /**< The directory of CLASS. This is set to the absolute path to the CLASS directory so this is just a failsafe. */
#endif

#define _LINE_LENGTH_MAX_ \
  1024 /**< size of the string read in each line of the file (extra characters not taken into account) */
#define _ARGUMENT_LENGTH_MAX_ \
  1024 /**< maximum size of each argument (name or value), including the final null character */

typedef char FileArg[_ARGUMENT_LENGTH_MAX_];

/**
 * @name Some conversion factors and fundamental constants needed by background module:
 */

//@{

#define _Mpc_over_m_ 3.085677581282e22 /**< conversion factor from meters to megaparsecs */
/* remark: CAMB uses 3.085678e22: good to know if you want to compare  with high accuracy */

#define _Gyr_over_Mpc_ \
  3.06601394e2               /**< conversion factor from megaparsecs to gigayears
                 (c=1 units, Julian years of 365.25 days) */
#define _c_ 2.99792458e8     /**< c in m/s */
#define _G_ 6.67428e-11      /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19 /**< 1 eV expressed in J */

/* parameters entering in Stefan-Boltzmann constant sigma_B */
#define _k_B_ 1.3806504e-23
#define _h_P_ 6.62606896e-34
/* remark: sigma_B = 2 pi^5 k_B^4 / (15h^3c^2) = 5.670400e-8
                   = Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

//@}

#define MIN(a, b) (((a) < (b)) ? (a) : (b)) /**< the usual "min" function */
#define MAX(a, b) (((a) < (b)) ? (b) : (a)) /**< the usual "max" function */
#define SIGN(a) (((a) > 0) ? 1. : -1.)
#define NRSIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define index_symmetric_matrix(i1, i2, N)            \
  (((i1) <= (i2))                                    \
       ? ((i2) + N * (i1) - ((i1) * ((i1) + 1)) / 2) \
       : ((i1) + N * (i2) -                          \
          ((i2) * ((i2) + 1)) /                      \
              2)) /**< assigns an index from 0 to [N(N+1)/2-1] to the coefficients M_{i1,i2} of an N*N symmetric matrix; useful for converting a symmetric matrix to a vector, without losing or double-counting any information */
/* @endcond */
// needed because of weird openmp bug on macosx lion...

#define class_build_error_string(dest, tmpl, ...)                          \
  {                                                                        \
    ErrorMsg FMsg;                                                         \
    class_protect_sprintf(FMsg, tmpl, __VA_ARGS__);                        \
    class_protect_sprintf(dest, "%s(L:%d) :%s", __func__, __LINE__, FMsg); \
  }

// Error reporting macros

// Call
#define class_call_message(err_out, extra, err_mess) \
  class_build_error_string(err_out, "error in %s;\n=>%s", extra, err_mess);

/* macro for calling function and returning error if it failed */
#define class_call_except(function,                                                     \
                          error_message_from_function,                                  \
                          error_message_output,                                         \
                          list_of_commands)                                             \
  {                                                                                     \
    if (function == _FAILURE_) {                                                        \
      class_call_message(error_message_output, #function, error_message_from_function); \
      list_of_commands;                                                                 \
      return _FAILURE_;                                                                 \
    }                                                                                   \
  }

/* macro for trying to call function */
#define class_call_try(function,                                                        \
                       error_message_from_function,                                     \
                       error_message_output,                                            \
                       list_of_commands)                                                \
  {                                                                                     \
    if (function == _FAILURE_) {                                                        \
      class_call_message(error_message_output, #function, error_message_from_function); \
      list_of_commands;                                                                 \
    }                                                                                   \
  }

/* macro for calling function and returning error if it failed */
#define class_call(function, error_message_from_function, error_message_output) \
  class_call_except(function, error_message_from_function, error_message_output, )

/* same in parallel region */
#define class_call_parallel(function, error_message_from_function, error_message_output)  \
  {                                                                                       \
    if (abort == _FALSE_) {                                                               \
      if (function == _FAILURE_) {                                                        \
        class_call_message(error_message_output, #function, error_message_from_function); \
        abort = _TRUE_;                                                                   \
      }                                                                                   \
    }                                                                                     \
  }

// Testing

#define class_test_message(err_out, extra, args, ...)                                           \
  {                                                                                             \
    ErrorMsg Optional_arguments;                                                                \
    class_protect_sprintf(Optional_arguments, args, ##__VA_ARGS__);                             \
    class_build_error_string(err_out, "condition (%s) is true; %s", extra, Optional_arguments); \
  }

/* macro for testing condition and returning error if condition is true;
   args is a variable list of optional arguments, e.g.: args="x=%d",x
   args cannot be empty, if there is nothing to pass use args="" */
#define class_test_except(condition, error_message_output, list_of_commands, args, ...) \
  {                                                                                     \
    if (condition) {                                                                    \
      class_test_message(error_message_output, #condition, args);                       \
      list_of_commands;                                                                 \
      return _FAILURE_;                                                                 \
    }                                                                                   \
  }

#define class_test(condition, error_message_output, args, ...)    \
  {                                                               \
    if (condition) {                                              \
      class_test_message(error_message_output, #condition, args); \
      return _FAILURE_;                                           \
    }                                                             \
  }

#define class_test_parallel(condition, error_message_output, args, ...) \
  {                                                                     \
    if (abort == _FALSE_) {                                             \
      if (condition) {                                                  \
        class_test_message(error_message_output, #condition, args);     \
        abort = _TRUE_;                                                 \
      }                                                                 \
    }                                                                   \
  }

/* macro for returning error message;
   args is a variable list of optional arguments, e.g.: args="x=%d",x
   args cannot be empty, if there is nothing to pass use args="" */
#define class_stop(error_message_output, args, ...)                                  \
  {                                                                                  \
    ErrorMsg Optional_arguments;                                                     \
    class_protect_sprintf(Optional_arguments, args, ##__VA_ARGS__);                  \
    class_build_error_string(error_message_output, "error; %s", Optional_arguments); \
    return _FAILURE_;                                                                \
  }

// IO
/* macro for opening file and returning error if it failed */
#define class_open(pointer, filename, mode, error_output)                    \
  {                                                                          \
    pointer = fopen(filename, mode);                                         \
    if (pointer == NULL) {                                                   \
      class_build_error_string(error_output,                                 \
                               "could not open %s with name %s and mode %s", \
                               #pointer,                                     \
                               filename,                                     \
                               #mode);                                       \
      return _FAILURE_;                                                      \
    }                                                                        \
  }

/**
 * C++ exception-throwing overrides.
 *
 * In C++ compilation units, the error macros throw std::runtime_error
 * instead of returning _FAILURE_. This enables natural exception
 * propagation through the call chain:
 *
 * - class_test / class_stop: throw at the error origin.
 * - class_call: if the callee returns _FAILURE_ (C function), builds
 *   the error chain and throws. If the callee itself throws (C++
 *   function), the exception propagates naturally.
 * - class_call_except / class_test_except: execute cleanup commands
 *   before re-throwing.
 * - class_open: throws on file open failure.
 *
 * The C definitions above remain active for .c utility files
 * (arrays.c, hyperspherical.c, etc.).
 */
#ifdef __cplusplus

#undef class_call_except
#define class_call_except(function,                                                       \
                          error_message_from_function,                                    \
                          error_message_output,                                           \
                          list_of_commands)                                               \
  {                                                                                       \
    try {                                                                                 \
      if ((function) == _FAILURE_) {                                                      \
        class_call_message(error_message_output, #function, error_message_from_function); \
        throw std::runtime_error(error_message_output);                                   \
      }                                                                                   \
    }                                                                                     \
    catch (...) {                                                                         \
      list_of_commands;                                                                   \
      throw;                                                                              \
    }                                                                                     \
  }

#undef class_call
#define class_call(function, error_message_from_function, error_message_output)         \
  {                                                                                     \
    if ((function) == _FAILURE_) {                                                      \
      class_call_message(error_message_output, #function, error_message_from_function); \
      throw std::runtime_error(error_message_output);                                   \
    }                                                                                   \
  }

#undef class_test_except
#define class_test_except(condition, error_message_output, list_of_commands, args, ...) \
  {                                                                                     \
    if (condition) {                                                                    \
      class_test_message(error_message_output, #condition, args, ##__VA_ARGS__);        \
      list_of_commands;                                                                 \
      throw std::runtime_error(error_message_output);                                   \
    }                                                                                   \
  }

#undef class_test
#define class_test(condition, error_message_output, args, ...)                   \
  {                                                                              \
    if (condition) {                                                             \
      class_test_message(error_message_output, #condition, args, ##__VA_ARGS__); \
      throw std::runtime_error(error_message_output);                            \
    }                                                                            \
  }

#undef class_stop
#define class_stop(error_message_output, args, ...)                                  \
  {                                                                                  \
    ErrorMsg Optional_arguments;                                                     \
    class_protect_sprintf(Optional_arguments, args, ##__VA_ARGS__);                  \
    class_build_error_string(error_message_output, "error; %s", Optional_arguments); \
    throw std::runtime_error(error_message_output);                                  \
  }

#undef class_open
#define class_open(pointer, filename, mode, error_output)                    \
  {                                                                          \
    pointer = fopen(filename, mode);                                         \
    if (pointer == NULL) {                                                   \
      class_build_error_string(error_output,                                 \
                               "could not open %s with name %s and mode %s", \
                               #pointer,                                     \
                               filename,                                     \
                               #mode);                                       \
      throw std::runtime_error(error_output);                                \
    }                                                                        \
  }

#endif /* __cplusplus exception overrides */

/* macro for defining indices (usually one, sometimes a block) */
#define class_define_index(index, condition, running_index, number_of_indices) \
  {                                                                            \
    if (condition) {                                                           \
      index          = running_index;                                          \
      running_index += number_of_indices;                                      \
    }                                                                          \
  }

/* macros for writing formatted output */
#define class_fprintf_double(file, output, condition)                    \
  {                                                                      \
    if (condition == _TRUE_)                                             \
      fprintf(file, "%*.*e ", _COLUMNWIDTH_, _OUTPUTPRECISION_, output); \
  }

#define class_fprintf_double_or_default(file, output, condition, defaultvalue) \
  {                                                                            \
    if (condition == _TRUE_)                                                   \
      fprintf(file, "%*.*e ", _COLUMNWIDTH_, _OUTPUTPRECISION_, output);       \
    else                                                                       \
      fprintf(file, "%*.*e ", _COLUMNWIDTH_, _OUTPUTPRECISION_, defaultvalue); \
  }

#define class_fprintf_int(file, output, condition)           \
  {                                                          \
    if (condition == _TRUE_)                                 \
      fprintf(file,                                          \
              "%*d%*s ",                                     \
              MAX(0, _COLUMNWIDTH_ - _OUTPUTPRECISION_ - 5), \
              output,                                        \
              _OUTPUTPRECISION_ + 5,                         \
              " ");                                          \
  }

#define class_fprintf_columntitle(file, title, condition, colnum)  \
  {                                                                \
    if (condition == _TRUE_)                                       \
      fprintf(file,                                                \
              "%*s%2d:%-*s ",                                      \
              MAX(0,                                               \
                  MIN(_COLUMNWIDTH_ - _OUTPUTPRECISION_ - 6 - 3,   \
                      _COLUMNWIDTH_ - ((int) strlen(title)) - 3)), \
              "",                                                  \
              colnum++,                                            \
              _OUTPUTPRECISION_ + 6,                               \
              title);                                              \
  }

#define class_store_columntitle(titlestring, title, condition) \
  {                                                            \
    if (condition == _TRUE_) {                                 \
      strcat(titlestring, title);                              \
      strcat(titlestring, _DELIMITER_);                        \
    }                                                          \
  }
//,_MAXTITLESTRINGLENGTH_-strlen(titlestring)-1);

#define class_store_double(storage, value, condition, dataindex) \
  {                                                              \
    if (condition == _TRUE_)                                     \
      storage[dataindex++] = value;                              \
  }

#define class_store_double_or_default(storage, value, condition, dataindex, defaultvalue) \
  {                                                                                       \
    if (condition == _TRUE_)                                                              \
      storage[dataindex++] = value;                                                       \
    else                                                                                  \
      storage[dataindex++] = defaultvalue;                                                \
  }

/** parameters related to the precision of the code and to the method of calculation */

/**
 * list of evolver types for integrating perturbations over time
 */
enum evolver_type {
  rk,   /* Runge-Kutta integrator */
  ndf15 /* stiff integrator */
};

/**
 * List of ways in which matter power spectrum P(k) can be defined.
 * The standard definition is the first one (delta_m_squared) but
 * alternative definitions can be useful in some projects.
 *
 */
enum pk_def {
  delta_m_squared, /**< normal definition (delta_m includes all non-relativistic species at late times) */
  delta_tot_squared, /**< delta_tot includes all species contributions to (delta rho), and only non-relativistic contributions to rho */
  delta_bc_squared, /**< delta_bc includes contribution of baryons and cdm only to (delta rho) and to rho */
  delta_tot_from_poisson_squared /**< use delta_tot inferred from gravitational potential through Poisson equation */
};
/**
 * Different ways to present output files
 */

enum file_format { class_format, camb_format };

/** Approximation scheme labels used as precision parameter defaults */

enum tca_method {
  first_order_MB,
  first_order_CAMB,
  first_order_CLASS,
  second_order_CRS,
  second_order_CLASS,
  compromise_CLASS
};
enum rsa_method { rsa_null, rsa_MD, rsa_MD_with_reio, rsa_none };
enum rsa_idr_method { rsa_idr_none, rsa_idr_MD };
enum ufa_method { ufa_mb, ufa_hu, ufa_CLASS, ufa_none };
enum ncdmfa_method { ncdmfa_mb, ncdmfa_hu, ncdmfa_CLASS, ncdmfa_none };

/**
 * All precision parameters.
 *
 * Includes integrations
 * steps, flags telling how the computation is to be performed, etc.
 */
struct precision {
  FileArg class_dir;

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
   */
  double tol_background_integration = 1.e-2;
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
   * non-cold dark matter phase-space distributions.
   */
  double tol_ncdm = 1.e-3;
  /**
   * Tolerance on the relative precision of the integration over
   * non-cold dark matter phase-space distributions in the synchronous gauge.
   */
  double tol_ncdm_synchronous = 1.e-3;
  /**
   * Tolerance on the relative precision of the integration over
   * non-cold dark matter phase-space distributions in the newtonian gauge.
   */
  double tol_ncdm_newtonian = 1.e-5;
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
  FileName sBBN_file;

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
   * Tolerance of the relative value of integral during thermodynamical integration
   */
  double tol_thermo_integration = 1.0e-2;
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
      _TRUE_; /**< from recfast 1.5, specifies how accurate the Hydrogen recombination should be handled */
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
  double recfast_x_He0_trigger2 =
      0.995; /**< Switch for Helium full calculation during reco, for changing Helium flag, raised from 0.985 to same as previous one for smoother Helium */
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

  FileName hyrec_Alpha_inf_file; /**< File containing the alpha parameter of hyrec */
  FileName hyrec_R_inf_file;     /**< File containing the R_inf parameter of hyrec */
  FileName
      hyrec_two_photon_tables_file; /**< File containing the two-photon interaction parameter of hyrec */

  double k_min_tau0 =
      0.1; /**< number defining k_min for the computation of Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */

  double k_max_tau0_over_l_max =
      2.4; /**< number defining k_max for the computation of Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */
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

  int tight_coupling_approximation =
      (int) compromise_CLASS; /**< method for tight coupling approximation */

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
  int radiation_streaming_approximation = rsa_MD_with_reio;

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

  int idr_streaming_approximation =
      rsa_idr_none; /**< method for dark radiation free-streaming approximation */
  double idr_streaming_trigger_tau_over_tau_k =
      50.0; /**< when to switch on dark radiation (idr) free-streaming approximation, first condition */
  double idr_streaming_trigger_tau_c_over_tau =
      10.0; /**< when to switch on dark radiation (idr) free-streaming approximation, second condition */

  int ur_fluid_approximation = ufa_CLASS; /**< method for ultra relativistic fluid approximation */

  /**
   * when to switch off ur (massless neutrinos / ultra-relativistic
   * relics) fluid approximation
   */
  double ur_fluid_trigger_tau_over_tau_k = 30.0;

  int ncdm_fluid_approximation =
      ncdmfa_CLASS; /**< method for non-cold dark matter fluid approximation */

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
   * The type of evolver to use: options are ndf15 or rk
   */
  enum evolver_type evolver = ndf15;

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

  double nonlinear_min_k_max = 20.0; /**< when
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
      _FALSE_; /**< switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles */
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

  /** @name - zone for writing error messages */
  //@{
  ErrorMsg error_message; /**< zone for writing error messages */
  //@}

  /** Parse precision parameters from a configuration file. */
  void parse(const FileContent& fc);
};

#ifdef __cplusplus
extern "C" {
#endif
void class_protect_sprintf(char* dest, const char* tpl, ...);
void class_protect_fprintf(FILE* dest, char* tpl, ...);
void* class_protect_memcpy(void* dest, void* from, size_t sz);
int get_number_of_titles(char* titlestring);
#ifdef __cplusplus
}
#endif

#endif
