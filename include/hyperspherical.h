/**
 * definitions for module hyperspherical.c
 */

#ifndef __HYPERSPHERICAL__
#define __HYPERSPHERICAL__

#include <algorithm>
#include <vector>

#include "common.h"
#include "function_ref.h"

#define _HYPER_OVERFLOW_ 1e200
#define _ONE_OVER_HYPER_OVERFLOW_ 1e-200
#define _HYPER_SAFETY_ 1e-5
#define _TRIG_PRECISSION_ 1e-7
#define _HYPER_BLOCK_ 8
#define _HYPER_CHUNK_ 16
#define _TWO_OVER_THREE_ 0.666666666666666666666666666667e0
#define _HIS_BYTE_ALIGNMENT_ 16

struct HypersphericalInterpolationStructure {
  int K          = 0;  //Sign of the curvature, (0,-1,1)
  double beta    = 0.0;
  double delta_x = 0.0;               //x-spacing. (xvec is uniformly spaced)
  int trig_order = 0;                 //Order of the interpolation formula for SinK and CosK.
  int l_size     = 0;                 //Number of l values
  std::vector<int> l;                 //Vector of l values stored
  std::vector<double> chi_at_phimin;  // vector x_min[index-l] below which neglect Bessels
  int x_size = 0;                     //Number of x-values
  std::vector<double> x;              //Pointer to x-values
  std::vector<double> sinK;           //Vector of sin_K(xvec)
  std::vector<double> cotK;           //Vector of cot_K(xvec)
  std::vector<double> phi;            //array of size nl*nx. [y_{l1}(x1) t_{l1}(x2)...]
  std::vector<double> dphi;           //Same as phivec, but containing derivatives.
  /* Second and third derivatives at the nodes, from the radial equation, in the layout of phi.
     Filled by OnFixedSpacing only, for hyperspherical_Hermite_evaluate; empty otherwise. */
  std::vector<double> d2phi;
  std::vector<double> d3phi;

  HypersphericalInterpolationStructure() = default;
  /** A table on [xmin, xmax] with about `sampling` points per period 2 pi / beta. The spacing
      is fitted so that the last node is exactly xmax, so it changes with xmax. */
  HypersphericalInterpolationStructure(int K,
                                       double beta,
                                       int nl,
                                       const int* lvec,
                                       double xmin,
                                       double xmax,
                                       double sampling,
                                       int l_WKB,
                                       double phiminabs);

  /** The same table, but on the fixed spacing 2 pi / (beta sampling) from xmin, reaching at
      least xmax. The nodes do not depend on xmax, so a table built for a larger xmax holds
      the same nodes with the same values and can stand in for this one (BesselTableCache).
      The node count is padded so that no node falls in a partial recurrence chunk. */
  static HypersphericalInterpolationStructure OnFixedSpacing(int K,
                                                             double beta,
                                                             int nl,
                                                             const int* lvec,
                                                             double xmin,
                                                             double xmax,
                                                             double sampling,
                                                             int l_WKB,
                                                             double phiminabs);

 private:
  struct Grid {
    int nx;
    double deltax;
  };
  HypersphericalInterpolationStructure(int K,
                                       double beta,
                                       int nl,
                                       const int* lvec,
                                       double xmin,
                                       Grid grid,
                                       int l_WKB,
                                       double phiminabs);
};
using HyperInterpStruct = HypersphericalInterpolationStructure;

/** Where arguments x fall in a table's grid, and the cubic Hermite weights there: at x_j,
    f(x_j) = h_left[j] f_i + h_right[j] f_{i+1} + g_left[j] f'_i + g_right[j] f'_{i+1}, i = left[j].
    They depend on the x only, so a caller that interpolates many multipoles at the same arguments
    (the flat-space transfer functions) finds them once and evaluates every multipole from them.
    Separate arrays rather than one struct per point, so that the weights stream through the cache. */
struct HermiteLocations {
  std::vector<int>
      left; /* left node, or -1 when x is outside the table (the functions are 0 there) */
  std::vector<double> h_left;
  std::vector<double> h_right;
  std::vector<double> g_left; /* the derivative weights include the node spacing */
  std::vector<double> g_right;
  void resize(int n) {
    left.resize(n);
    h_left.resize(n);
    h_right.resize(n);
    g_left.resize(n);
    g_right.resize(n);
  }
};

/** Locate nxi arguments in a table (flat: K = 0). */
void hyperspherical_Hermite_locate(const HyperInterpStruct* pHIS,
                                   int nxi,
                                   const double* xinterp,
                                   HermiteLocations& where);

/** Phi, dPhi and/or d2Phi of multipole index lnum at the located arguments first..first+nxi-1:
    the same cubic Hermite interpolants as hyperspherical_Hermite_interpolation<4, ...>, from the
    node tables of a table built by OnFixedSpacing. */
template <bool DoPhi, bool DoDPhi, bool DoD2Phi>
void hyperspherical_Hermite_evaluate(const HyperInterpStruct* pHIS,
                                     int lnum,
                                     const HermiteLocations& where,
                                     int first,
                                     int nxi,
                                     double* Phi,
                                     double* dPhi,
                                     double* d2Phi);

struct WKB_parameters {
  int K;
  int l;
  double beta;
  double phiminabs;
};

void hyperspherical_forwards_recurrence(int K,
                                        int lmax,
                                        double beta,
                                        double x,
                                        double sinK,
                                        double cotK,
                                        const double* sqrtK,
                                        const double* one_over_sqrtK,
                                        double* PhiL);
void hyperspherical_forwards_recurrence_chunk(int K,
                                              int lmax,
                                              double beta,
                                              double* x,
                                              double* sinK,
                                              double* cotK,
                                              int chunk,
                                              const double* sqrtK,
                                              const double* one_over_sqrtK,
                                              double* PhiL);
void hyperspherical_backwards_recurrence(int K,
                                         int lmax,
                                         double beta,
                                         double x,
                                         double sinK,
                                         double cotK,
                                         const double* sqrtK,
                                         const double* one_over_sqrtK,
                                         double* PhiL);

void hyperspherical_backwards_recurrence_chunk(int K,
                                               int lmax,
                                               double beta,
                                               double* x,
                                               double* sinK,
                                               double* cotK,
                                               int chunk,
                                               const double* sqrtK,
                                               const double* one_over_sqrtK,
                                               double* PhiL);

bool hyperspherical_WKB(int K, int l, double beta, double y, double* Phi);
void hyperspherical_bessel_direct_vector(
    int K, double beta, int* lvec, int nl, const double* xvec, int nx, double* Phi);
void ClosedModY(int l, int beta, double* y, int* phisign, int* dphisign);
bool get_CF1(int K, int l, double beta, double cotK, double* CF, int* isign);
bool CF1_from_Gegenbauer(int l, int beta, double sinK, double cotK, double* CF);
double airy_cheb_approx(double z);
double coef1(double z);
double coef2(double z);
double coef3(double z);
double coef4(double z);
double cheb(double x, int n, const double A[]);

double PhiWKB_minus_phiminabs(double x, const WKB_parameters& wkbparam);

void hyperspherical_get_xmin_from_Airy(
    int K, int l, double beta, double xtol, double phiminabs, double* xmin, int* fevals);

bool fzero_ridder(FunctionRef<double(double x)> func,
                  double x1,
                  double x2,
                  double xtol,
                  const double* Fx1,
                  const double* Fx2,
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
extern template void hyperspherical_Hermite_evaluate<true, false, false>(
    const HyperInterpStruct*, int, const HermiteLocations&, int, int, double*, double*, double*);
extern template void hyperspherical_Hermite_evaluate<true, true, true>(
    const HyperInterpStruct*, int, const HermiteLocations&, int, int, double*, double*, double*);
extern template void hyperspherical_Hermite_evaluate<true, true, false>(
    const HyperInterpStruct*, int, const HermiteLocations&, int, int, double*, double*, double*);
extern template void hyperspherical_Hermite_evaluate<true, false, true>(
    const HyperInterpStruct*, int, const HermiteLocations&, int, int, double*, double*, double*);
extern template void hyperspherical_Hermite_evaluate<false, true, false>(
    const HyperInterpStruct*, int, const HermiteLocations&, int, int, double*, double*, double*);
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
