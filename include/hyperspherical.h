/**
 * definitions for module hyperspherical.c
 */

#ifndef __HYPERSPHERICAL__
#define __HYPERSPHERICAL__

#include <algorithm>
#include <vector>

#include "common.h"

#define _HYPER_OVERFLOW_ 1e200
#define _ONE_OVER_HYPER_OVERFLOW_ 1e-200
#define _HYPER_SAFETY_ 1e-5
#define _TRIG_PRECISSION_ 1e-7
#define _HYPER_BLOCK_ 8
#define _HYPER_CHUNK_ 16
#define _TWO_OVER_THREE_ 0.666666666666666666666666666667e0
#define _HIS_BYTE_ALIGNMENT_ 16

typedef struct HypersphericalInterpolationStructure {
  int K;  //Sign of the curvature, (0,-1,1)
  double beta;
  double delta_x;                     //x-spacing. (xvec is uniformly spaced)
  int trig_order;                     //Order of the interpolation formula for SinK and CosK.
  int l_size;                         //Number of l values
  std::vector<int> l;                 //Vector of l values stored
  std::vector<double> chi_at_phimin;  // vector x_min[index-l] below which neglect Bessels
  int x_size;                         //Number of x-values
  std::vector<double> x;              //Pointer to x-values
  std::vector<double> sinK;           //Vector of sin_K(xvec)
  std::vector<double> cotK;           //Vector of cot_K(xvec)
  std::vector<double> phi;            //array of size nl*nx. [y_{l1}(x1) t_{l1}(x2)...]
  std::vector<double> dphi;           //Same as phivec, but containing derivatives.

  HypersphericalInterpolationStructure() = default;
  HypersphericalInterpolationStructure(int K,
                                       double beta,
                                       int nl,
                                       const int* lvec,
                                       double xmin,
                                       double xmax,
                                       double sampling,
                                       int l_WKB,
                                       double phiminabs);
} HyperInterpStruct;

struct WKB_parameters {
  int K;
  int l;
  double beta;
  double phiminabs;
};

int hyperspherical_forwards_recurrence(int K,
                                       int lmax,
                                       double beta,
                                       double x,
                                       double sinK,
                                       double cotK,
                                       double* sqrtK,
                                       double* one_over_sqrtK,
                                       double* PhiL);
int hyperspherical_forwards_recurrence_chunk(int K,
                                             int lmax,
                                             double beta,
                                             double* x,
                                             double* sinK,
                                             double* cotK,
                                             int chunk,
                                             double* sqrtK,
                                             double* one_over_sqrtK,
                                             double* PhiL);
int hyperspherical_backwards_recurrence(int K,
                                        int lmax,
                                        double beta,
                                        double x,
                                        double sinK,
                                        double cotK,
                                        double* sqrtK,
                                        double* one_over_sqrtK,
                                        double* PhiL);

int hyperspherical_backwards_recurrence_chunk(int K,
                                              int lmax,
                                              double beta,
                                              double* x,
                                              double* sinK,
                                              double* cotK,
                                              int chunk,
                                              double* sqrtK,
                                              double* one_over_sqrtK,
                                              double* PhiL);

int hyperspherical_WKB(int K, int l, double beta, double y, double* Phi);
int hyperspherical_bessel_direct_vector(
    int K, double beta, int* lvec, int nl, double* xvec, int nx, double* Phi);
int ClosedModY(int l, int beta, double* y, int* phisign, int* dphisign);
int get_CF1(int K, int l, double beta, double cotK, double* CF, int* isign);
int CF1_from_Gegenbauer(int l, int beta, double sinK, double cotK, double* CF);
double airy_cheb_approx(double z);
double coef1(double z);
double coef2(double z);
double coef3(double z);
double coef4(double z);
double cheb(double x, int n, const double A[]);

double PhiWKB_minus_phiminabs(double x, void* param);

void hyperspherical_get_xmin_from_Airy(
    int K, int l, double beta, double xtol, double phiminabs, double* xmin, int* fevals);

int fzero_ridder(double (*func)(double, void*),
                 double x1,
                 double x2,
                 double xtol,
                 void* param,
                 double* Fx1,
                 double* Fx2,
                 double* xzero,
                 int* fevals);

void hyperspherical_get_xmin_from_approx(
    int K, int l, double nu, double ignore1, double phiminabs, double* xmin, int* ignore2);

/** Hermite interpolation (order 4 or 6) of the hyperspherical Bessel function
    Phi and/or its first two derivatives, evaluated at the nxi points xinterp
    for the multipole pHIS->l[lnum]. Only the outputs selected by the template
    flags are computed and stored (pass nullptr for the others).

    When xinterp is sorted (increasing), computations are reused across
    points; for randomly ordered input the routine is not much slower. For the
    closed case the interpolation structure only covers
    [safety; pi/2-safety]; the calling routine should respect this.

    Each requested output at derivative level m (Phi: m=0, dPhi: m=1,
    d2Phi: m=2) is Hermite-interpolated from levels m..m+Order/2-1 at the two
    surrounding grid points. The tabulated Phi and dPhi are lifted to higher
    derivatives through the radial equation
      Phi'' = -2 cotK Phi' + (l(l+1)/sinK^2 - beta^2 + K) Phi
    differentiated as often as needed. */
template <int Order, bool DoPhi, bool DoDPhi, bool DoD2Phi>
void hyperspherical_Hermite_interpolation(const HyperInterpStruct* pHIS,
                                          int nxi,
                                          int lnum,
                                          const double* xinterp,
                                          double* Phi,
                                          double* dPhi,
                                          double* d2Phi);

// The definition lives in hyperspherical.cpp (same translation unit as
// ClosedModY, which must inline into the closed-case loop); these are the
// combinations the transfer module dispatches to.
extern template void hyperspherical_Hermite_interpolation<4, true, false, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<4, false, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<4, true, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<4, true, false, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<4, true, true, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<6, true, false, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<6, false, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<6, true, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<6, true, false, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
extern template void hyperspherical_Hermite_interpolation<6, true, true, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);

#endif
