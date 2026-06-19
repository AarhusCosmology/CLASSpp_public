/** @file constants.h Physics constants and conversion factors used across CLASS. */
#ifndef CLASS_CONSTANTS_H
#define CLASS_CONSTANTS_H

#include <cmath>  // pow(), used by the Omega <-> T_cmb helpers below

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _PIHALF_ 1.57079632679489661923132169164e0 /**< pi divided by 2 */

#define _TWOPI_ 6.283185307179586476925286766559e0 /**< 2 times pi */

#define _SQRT2_ 1.41421356237309504880168872421e0 /** < square root of 2. */

#define _SQRT6_ 2.4494897427831780981972840747059e0 /**< square root of 6. */

#define _SQRT_PI_ 1.77245385090551602729816748334e0 /**< square root of pi. */

#define _E_ \
  2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696 /**< exponential of one */

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

/** Photon density parameter from CMB temperature.
 *  Omega0_g = (4 sigma_B / c) T^4 / (3 c^2 1e10 h^2 / Mpc_over_m^2 / 8 pi G),
 *  with sigma_B = 2 pi^5 k_B^4 / (15 h_P^3 c^2). Single source of truth for the
 *  T_cmb <-> Omega0_g conversion. Free function (not a PhotonsSpecies method) so
 *  that background.h's default member initializer can call it without the
 *  photons.h <-> background.h include cycle. */
inline double Omega0gFromTcmb(double T_cmb, double h) {
  const double sigma_B = 2. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 2);
  return (4. * sigma_B / _c_ * pow(T_cmb, 4.)) /
         (3. * _c_ * _c_ * 1.e10 * h * h / _Mpc_over_m_ / _Mpc_over_m_ / 8. / _PI_ / _G_);
}

/** Inverse of Omega0gFromTcmb. */
inline double TcmbFromOmega0g(double Omega0_g, double h) {
  const double sigma_B = 2. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 2);
  return pow(Omega0_g *
                 (3. * _c_ * _c_ * 1.e10 * h * h / _Mpc_over_m_ / _Mpc_over_m_ / 8. / _PI_ / _G_) /
                 (4. * sigma_B / _c_),
             0.25);
}

#endif  // CLASS_CONSTANTS_H
