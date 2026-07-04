/** @file hyperspherical.cpp Hyperspherical Bessel-function implementation.
 *
 * Thomas Tram, 11.01.2013
 *
 * This module computes hyperspherical Bessel functions.
 */

#include "hyperspherical.h"

#include <algorithm>
#define REFINE 10

namespace {

double stable_sinc(double arg) {
  double arg2 = arg * arg;
  if (fabs(arg) < 1e-4) {
    double arg4 = arg2 * arg2;
    return 1.0 - arg2 / 6.0 + arg4 / 120.0;
  }
  return sin(arg) / arg;
}

double stable_sinc_minus_cos(double arg) {
  double arg2 = arg * arg;
  if (fabs(arg) < 1e-3) {
    double arg4 = arg2 * arg2;
    return arg2 / 3.0 - arg4 / 30.0 + arg4 * arg2 / 840.0;
  }
  return stable_sinc(arg) - cos(arg);
}

void hyperspherical_phi01_exact(
    int K, double beta, double x, double sinK, double cotK, double* phi0, double* phi1) {
  double beta2 = beta * beta;
  double arg   = beta * x;
  double arg2  = arg * arg;
  double sinc  = stable_sinc(arg);

  if (K == 0) {
    *phi0 = sinc;
    if (fabs(arg) <= 1e-3) {
      double arg4 = arg2 * arg2;
      *phi1       = arg * (1.0 / 3.0 - arg2 / 30.0 + arg4 / 840.0);
    }
    else {
      *phi1 = stable_sinc_minus_cos(arg) / arg;
    }
    return;
  }

  double root1 = sqrt(std::max(0.0, beta2 - (double) K));
  if (fabs(x) < 1e-4) {
    double x2 = x * x;
    double x4 = x2 * x2;
    if (fabs(arg) < 1e-4) {
      *phi0 = 1.0 - x2 * (beta2 - K) / 6.0;
      *phi1 = x * root1 / 3.0 * (1.0 - (3.0 * beta2 - 7.0 * K) * x2 / 30.0);
      return;
    }

    double chi_over_sin = 1.0 + K * x2 / 6.0 + 7.0 * x4 / 360.0;
    double chi_cot_m1   = -K * x2 / 3.0 - x4 / 45.0;
    *phi0               = sinc * chi_over_sin;
    *phi1 = (stable_sinc_minus_cos(arg) + sinc * chi_cot_m1) * chi_over_sin / (root1 * x);
    return;
  }

  double chi_over_sin = x / sinK;
  double chi_cot_m1   = x * cotK - 1.0;
  *phi0               = sinc * chi_over_sin;
  if (root1 > 0.0) {
    *phi1 = (stable_sinc_minus_cos(arg) + sinc * chi_cot_m1) * chi_over_sin / (root1 * x);
  }
  else {
    *phi1 = 0.0;
  }
}

}  // namespace

HypersphericalInterpolationStructure::HypersphericalInterpolationStructure(int K,
                                                                           double beta,
                                                                           int nl,
                                                                           const int* lvec,
                                                                           double xmin,
                                                                           double xmax,
                                                                           double sampling,
                                                                           int l_WKB,
                                                                           double phiminabs) {
  /** Size the std::vector members for this geometry, then compute the values
      of Phi and dPhi to complete the interpolation structure. */
  double deltax, beta2, lambda, x, xfwd;
  int j, k, l, nx, lmax, l_recurrence_max;
  int current_chunk, index_x;

  beta2  = beta * beta;
  lmax   = lvec[nl - 1];
  lambda = 2 * _PI_ / beta;
  nx     = (int) ((xmax - xmin) * sampling / lambda);
  nx     = std::max(nx, 2);
  deltax = (xmax - xmin) / (nx - 1.0);
  //fprintf(stderr,"dx=%e\n",deltax);
  //fprintf(stderr,"%e %e\n",beta,sampling);
  //Set scalar values:
  this->beta    = beta;
  this->delta_x = deltax;
  this->l_size  = nl;
  this->x_size  = nx;
  this->K       = K;
  //Set pointervalues in pHIS:

  this->l.resize(nl);
  this->chi_at_phimin.resize(nl);
  this->x.resize(nx);
  this->sinK.resize(nx);
  this->cotK.resize(nx);
  this->phi.resize(nx * nl);
  this->dphi.resize(nx * nl);

  //Order needed for trig interpolation: (We are using Taylor's remainder theorem)
  if (0.5 * deltax * deltax < _TRIG_PRECISSION_)
    this->trig_order = 1;
  else if ((pow(deltax, 4) / 24.0) < _TRIG_PRECISSION_)
    this->trig_order = 3;
  else
    this->trig_order = 5;

  //Copy lvector:
  for (j = 0; j < nl; j++) {
    this->l[j] = lvec[j];
  }
  //Allocate sqrtK, and PhiL:
  std::vector<double> sqrtK(lmax + 3);
  std::vector<double> one_over_sqrtK(lmax + 3);

  //Find l_WKB_min, the highest l in lvec where l<l_WKB:
  l_recurrence_max         = -10;
  int index_recurrence_max = -10;
  for (k = nl - 1; k >= 0; k--) {
    l = lvec[k];
    if (l < l_WKB) {
      l_recurrence_max     = l;
      index_recurrence_max = k;
      break;
    }
  }

  //Create xvector and set x, cotK, sinK, sqrtK and fwdidx:
  switch (K) {
    case 0:
      xfwd = sqrt(l_recurrence_max * (l_recurrence_max + 1.0)) / beta;
      for (j = 0; j < nx; j++) {
        x             = xmin + j * deltax;
        this->x[j]    = x;
        this->sinK[j] = x;
        this->cotK[j] = 1.0 / x;
      }
      for (l = 0; l <= (lmax + 2); l++) {
        sqrtK[l]          = beta;
        one_over_sqrtK[l] = 1.0 / sqrtK[l];
      }
      break;
    case 1:
      xfwd = asin(sqrt(l_recurrence_max * (l_recurrence_max + 1.0)) / beta);
      for (j = 0; j < nx; j++) {
        x             = xmin + j * deltax;
        this->x[j]    = x;
        this->sinK[j] = sin(x);
        this->cotK[j] = 1.0 / tan(x);
      }
      for (l = 0; l <= (lmax + 2); l++) {
        sqrtK[l]          = sqrt(beta2 - l * l);
        one_over_sqrtK[l] = 1.0 / sqrtK[l];
      }
      break;
    case -1:
      xfwd = asinh(sqrt(l_recurrence_max * (l_recurrence_max + 1.0)) / beta);
      for (j = 0; j < nx; j++) {
        x             = xmin + j * deltax;
        this->x[j]    = x;
        this->sinK[j] = sinh(x);
        this->cotK[j] = 1.0 / tanh(x);
      }
      for (l = 0; l <= (lmax + 2); l++) {
        sqrtK[l]          = sqrt(beta2 + l * l);
        one_over_sqrtK[l] = 1.0 / sqrtK[l];
      }
      break;
    default:
      class_test(1, "K must be -1, 0, or 1 (got %d)", K);
  }

  int xfwdidx = (xfwd - xmin) / deltax;
  // Clamp into [0, nx]. When the classical turning point xfwd lies below xmin
  // (small l: xfwd ~ sqrt(l(l+1))/beta can be < xmin, e.g. 0 for l=0), the raw
  // value is negative and the forwards loop below would start at a negative j,
  // indexing this->x / phi / dphi / cotK out of bounds (an OOB read of x AND an
  // OOB write of phi/dphi 3 elements before the buffer). This is latent in the
  // transfer module (large multipoles keep xfwd > xmin) but is reachable when
  // the constructor is called directly with small l.
  if (xfwdidx < 0)
    xfwdidx = 0;
  if (xfwdidx > nx)
    xfwdidx = nx;
  //Calculate and assign Phi and dPhi values:

  std::vector<double> PhiL((lmax + 2) * _HYPER_CHUNK_);

  int hit_the_ceiling = ((K == 1) && ((int) (beta + 0.2) == (lmax + 1)));
  if (hit_the_ceiling) {
    /** Take care of special case lmax = beta-1.
        The routine below will try to compute
        Phi_{lmax+1} which is not allowed. However,
        the purpose is to calculate the derivative
        Phi'_{lmax}, and the formula is correct if we set Phi_{lmax+1} = 0.
        Since PhiL uses a chunked layout PhiL[l*chunk + index_x], we
        cannot simply set PhiL[lmax+1] = 0 (that was a scalar index).
        Instead we substitute 0.0 at the point of use in the dphi formula.
    */
    lmax--;
  }

  for (j = 0; j < std::min(nx, xfwdidx); j += _HYPER_CHUNK_) {
    current_chunk = std::min(_HYPER_CHUNK_, std::min(nx, xfwdidx) - j);
    //Use backwards method (chunk version for better SIMD utilization):
    hyperspherical_backwards_recurrence_chunk(K,
                                              std::min(l_recurrence_max, lmax) + 1,
                                              beta,
                                              this->x.data() + j,
                                              this->sinK.data() + j,
                                              this->cotK.data() + j,
                                              current_chunk,
                                              sqrtK.data(),
                                              one_over_sqrtK.data(),
                                              PhiL.data());
    //We have now populated PhiL at x, assign Phi and dPhi for all l in lvec:
    for (k = 0; k <= index_recurrence_max; k++) {
      l = lvec[k];
      for (index_x = 0; index_x < current_chunk; index_x++) {
        this->phi[k * nx + j + index_x] = PhiL[l * current_chunk + index_x];
        double PhiL_plus_one =
            ((hit_the_ceiling && l > lmax) ? 0.0 : PhiL[(l + 1) * current_chunk + index_x]);
        this->dphi[k * nx + j + index_x] = l * this->cotK[j + index_x] *
                                               PhiL[l * current_chunk + index_x] -
                                           sqrtK[l + 1] * PhiL_plus_one;
      }
    }
  }

  for (j = xfwdidx; j < nx; j += _HYPER_CHUNK_) {
    //Use forwards method:
    current_chunk = std::min(_HYPER_CHUNK_, nx - j);
    hyperspherical_forwards_recurrence_chunk(K,
                                             std::min(l_recurrence_max, lmax) + 1,
                                             beta,
                                             this->x.data() + j,
                                             this->sinK.data() + j,
                                             this->cotK.data() + j,
                                             current_chunk,
                                             sqrtK.data(),
                                             one_over_sqrtK.data(),
                                             PhiL.data());

    //We have now populated PhiL at x, assign Phi and dPhi for all l in lvec:
    for (k = 0; k <= index_recurrence_max; k++) {
      l = lvec[k];
      for (index_x = 0; index_x < current_chunk; index_x++) {
        this->phi[k * nx + j + index_x] = PhiL[l * current_chunk + index_x];
        double PhiL_plus_one =
            ((hit_the_ceiling && l > lmax) ? 0.0 : PhiL[(l + 1) * current_chunk + index_x]);
        this->dphi[k * nx + j + index_x] = l * this->cotK[j + index_x] *
                                               PhiL[l * current_chunk + index_x] -
                                           sqrtK[l + 1] * PhiL_plus_one;
      }
    }
  }

  for (k = 0; k < nl; k++) {
    hyperspherical_get_xmin_from_approx(K,
                                        lvec[k],
                                        beta,
                                        0.,
                                        phiminabs,
                                        this->chi_at_phimin.data() + k,
                                        nullptr);
  }
}

int hyperspherical_forwards_recurrence(int K,
                                       int lmax,
                                       double beta,
                                       double x,
                                       double sinK,
                                       double cotK,
                                       double* sqrtK,
                                       double* one_over_sqrtK,
                                       double* PhiL) {
  int l;
  hyperspherical_phi01_exact(K, beta, x, sinK, cotK, PhiL, PhiL + 1);
  for (l = 2; l <= lmax; l++) {
    PhiL[l] = ((2 * l - 1) * cotK * PhiL[l - 1] - PhiL[l - 2] * sqrtK[l - 1]) * one_over_sqrtK[l];
  }
  return _SUCCESS_;
}

int hyperspherical_forwards_recurrence_chunk(int K,
                                             int lmax,
                                             double beta,
                                             double* x,
                                             double* sinK,
                                             double* cotK,
                                             int chunk,
                                             double* sqrtK,
                                             double* one_over_sqrtK,
                                             double* PhiL) {
  int l;
  int index_x;
  for (index_x = 0; index_x < chunk; index_x++) {
    hyperspherical_phi01_exact(K,
                               beta,
                               x[index_x],
                               sinK[index_x],
                               cotK[index_x],
                               PhiL + index_x,
                               PhiL + chunk + index_x);
  }
  for (l = 2; l <= lmax; l++) {
    for (index_x = 0; index_x < chunk; index_x++)
      PhiL[l * chunk + index_x] = ((2 * l - 1) * cotK[index_x] * PhiL[(l - 1) * chunk + index_x] -
                                   PhiL[(l - 2) * chunk + index_x] * sqrtK[l - 1]) *
                                  one_over_sqrtK[l];
  }
  return _SUCCESS_;
}

int hyperspherical_backwards_recurrence(int K,
                                        int lmax,
                                        double beta,
                                        double x,
                                        double sinK,
                                        double cotK,
                                        double* sqrtK,
                                        double* one_over_sqrtK,
                                        double* PhiL) {
  double phi0_exact, phi1_exact, phi_top, phipr1 = 0.0, phi, phi_plus_1_times_sqrtK, phi_minus_1,
                                          phi1_down = 0.0, phi0_down, scaling;
  int l, k, isign;
  int funcreturn = _FAILURE_;
  hyperspherical_phi01_exact(K, beta, x, sinK, cotK, &phi0_exact, &phi1_exact);

  //printf("in backwards. x = %g\n",x);
  if (K == 1) {
    if (beta > 1.5 * lmax) {
      funcreturn = get_CF1(K, lmax, beta, cotK, &phipr1, &isign);
    }
    if (funcreturn == _FAILURE_) {
      CF1_from_Gegenbauer(lmax, (int) (beta + 0.2), sinK, cotK, &phipr1);
    }
    phi_top = 1.0;
  }
  else {
    get_CF1(K, lmax, beta, cotK, &phipr1, &isign);
    phi_top  = isign;
    phipr1  *= phi_top;
    //printf("isign = %d, phi_top = %g, phipr1 = %g\n",isign,phi_top,phipr1);
  }

  PhiL[lmax] = phi_top;
  phi        = phi_top;
  //  phi_plus_1 = 1/sqrtK[lmax+1]*(lmax*cotK*phi1-phipr1);
  phi_plus_1_times_sqrtK = lmax * cotK * phi_top - phipr1;

  int l_ini, l_align;
  l_align = lmax - lmax % _HYPER_BLOCK_;

  // Bring l down to _HYPER_BLOCK_ aligned region:
  for (l = lmax; l > l_align; l--) {
    //    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
    phi_minus_1 = ((2 * l + 1) * cotK * phi - phi_plus_1_times_sqrtK) * one_over_sqrtK[l];
    if (l == 1)
      phi1_down = phi;
    phi_plus_1_times_sqrtK = phi * sqrtK[l];
    phi                    = phi_minus_1;
    PhiL[l - 1]            = phi;
  }
  for (l_ini = l_align; l_ini > 0; l_ini -= _HYPER_BLOCK_) {
    for (l = l_ini; l > (l_ini - _HYPER_BLOCK_); l--) {
      //    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
      phi_minus_1 = ((2 * l + 1) * cotK * phi - phi_plus_1_times_sqrtK) * one_over_sqrtK[l];
      if (l == 1)
        phi1_down = phi;
      phi_plus_1_times_sqrtK = phi * sqrtK[l];
      phi                    = phi_minus_1;
      PhiL[l - 1]            = phi;
    }
    if (fabs(phi) > _HYPER_OVERFLOW_) {
      phi                    *= _ONE_OVER_HYPER_OVERFLOW_;
      phi_plus_1_times_sqrtK *= _ONE_OVER_HYPER_OVERFLOW_;
      //Rescale whole Phi vector until this point:
      for (k = l; k <= lmax; k++)
        PhiL[k] *= _ONE_OVER_HYPER_OVERFLOW_;
      phi1_down *= _ONE_OVER_HYPER_OVERFLOW_;
    }
  }

  /**
  for (l=lmax; l>=1; l--){
    //    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )*one_over_sqrtK[l];
    phi_plus_1_times_sqrtK = phi*sqrtK[l];
    phi = phi_minus_1;
    PhiL[l-1] = phi;
    if (fabs(phi)>_HYPER_OVERFLOW_){
      //
      phi *= _ONE_OVER_HYPER_OVERFLOW_;
      phi_plus_1_times_sqrtK *= _ONE_OVER_HYPER_OVERFLOW_;
      //Rescale whole Phi vector until this point:
      for (k=l-1; k<=lmax; k++)
        PhiL[k] *=_ONE_OVER_HYPER_OVERFLOW_;
    }
  }
  */
  phi0_down = phi;
  // Antony Lewis highlighted exact-angle failures at chi = pi/4 and pi/2 in
  // https://github.com/lesgourg/class_public/issues/662. When phi0 is tiny at
  // those points, normalize Miller recurrence with the safer exact phi1 seed.
  if (fabs(phi0_exact) >= fabs(phi1_exact)) {
    if (fabs(phi0_down) > 1e-280) {
      scaling = phi0_exact / phi0_down;
    }
    else {
      scaling = phi1_exact / phi1_down;
    }
  }
  else {
    if (fabs(phi1_down) > 1e-280) {
      scaling = phi1_exact / phi1_down;
    }
    else {
      scaling = phi0_exact / phi0_down;
    }
  }
  for (k = 0; k <= lmax; k++)
    PhiL[k] *= scaling;
  return _SUCCESS_;
}

int hyperspherical_backwards_recurrence_chunk(int K,
                                              int lmax,
                                              double beta,
                                              double* x,
                                              double* sinK,
                                              double* cotK,
                                              int chunk,
                                              double* sqrtK,
                                              double* one_over_sqrtK,
                                              double* PhiL) {
  double phi0_exact, phi1_exact, phi_top, phipr1 = 0.0;
  int l, k, isign;
  int funcreturn = _FAILURE_;
  int index_x;
  double scalevec[_HYPER_CHUNK_] = {0}, phi1_down[_HYPER_CHUNK_] = {0};

  for (index_x = 0; index_x < chunk; index_x++) {
    if (K == 1) {
      if (beta > 1.5 * lmax) {
        funcreturn = get_CF1(K, lmax, beta, cotK[index_x], &phipr1, &isign);
      }
      if (funcreturn == _FAILURE_) {
        CF1_from_Gegenbauer(lmax, (int) (beta + 0.2), sinK[index_x], cotK[index_x], &phipr1);
      }
      phi_top = 1.0;
    }
    else {
      get_CF1(K, lmax, beta, cotK[index_x], &phipr1, &isign);
      phi_top  = isign;
      phipr1  *= phi_top;
    }

    PhiL[lmax * chunk + index_x]       = phi_top;
    PhiL[(lmax - 1) * chunk + index_x] = one_over_sqrtK[lmax] *
                                         ((lmax + 1) * cotK[index_x] * phi_top + phipr1);
  }
  for (l = lmax - 2; l >= 0; l--) {
    //Use recurrence Phi_{l} = --Phi_{l+1} + -- Phi_{l+2}
    for (index_x = 0; index_x < chunk; index_x++) {
      PhiL[l * chunk + index_x] = one_over_sqrtK[l + 1] *
                                  ((2 * l + 3) * cotK[index_x] * PhiL[(l + 1) * chunk + index_x] -
                                   sqrtK[l + 2] * PhiL[(l + 2) * chunk + index_x]);
      if (l == 0)
        phi1_down[index_x] = PhiL[chunk + index_x];
    }

    /* Check if any entry in the chunk has overflowed — the original code only
       checked index_x=0 which caused NaN in closed-model Cl spectra. */
    double maxabs = 0.0;
    for (index_x = 0; index_x < chunk; index_x++) {
      double v = fabs(PhiL[l * chunk + index_x]);
      if (v > maxabs)
        maxabs = v;
    }
    if (maxabs > _HYPER_OVERFLOW_) {
      //Rescale whole Phi vector until this point.
      //Create scale vector:
      for (index_x = 0; index_x < chunk; index_x++)
        scalevec[index_x] = fabs(1.0 / PhiL[l * chunk + index_x]);
      //Now do the scaling: (We do it this way to access elements in order)
      for (k = l; k <= lmax; k++) {
        for (index_x = 0; index_x < chunk; index_x++) {
          PhiL[k * chunk + index_x] *= scalevec[index_x];
        }
      }
      for (index_x = 0; index_x < chunk; index_x++) {
        phi1_down[index_x] *= scalevec[index_x];
      }
    }
  }
  for (index_x = 0; index_x < chunk; index_x++) {
    hyperspherical_phi01_exact(K,
                               beta,
                               x[index_x],
                               sinK[index_x],
                               cotK[index_x],
                               &phi0_exact,
                               &phi1_exact);
    if (fabs(phi0_exact) >= fabs(phi1_exact)) {
      if (fabs(PhiL[index_x]) > 1e-280) {
        scalevec[index_x] = phi0_exact / PhiL[index_x];
      }
      else {
        scalevec[index_x] = phi1_exact / phi1_down[index_x];
      }
    }
    else {
      if (fabs(phi1_down[index_x]) > 1e-280) {
        scalevec[index_x] = phi1_exact / phi1_down[index_x];
      }
      else {
        scalevec[index_x] = phi0_exact / PhiL[index_x];
      }
    }
  }
  for (k = 0; k <= lmax; k++) {
    for (index_x = 0; index_x < chunk; index_x++) {
      PhiL[k * chunk + index_x] *= scalevec[index_x];
    }
  }

  return _SUCCESS_;
}

int get_CF1(int K, int l, double beta, double cotK, double* CF, int* isign) {
  int maxiter   = 1000000;
  double tiny   = 1e-100;
  double reltol = DBL_EPSILON;
  double aj, bj, fj, Cj, Dj, Delj;
  double beta2 = beta * beta;
  double sqrttmp;
  int j;

  if (K == 1)
    maxiter = (int) (beta - l - 10);
  bj     = l * cotK;  //This is b_0
  fj     = bj;
  Cj     = bj;
  Dj     = 0.0;
  *isign = 1;
  for (j = 1; j <= maxiter; j++) {
    sqrttmp = sqrt(beta2 - K * (l + j + 1) * (l + j + 1));
    aj      = -sqrt(beta2 - K * (l + j) * (l + j)) / sqrttmp;
    if (j == 1)
      aj = sqrt(beta2 - K * (l + 1) * (l + 1)) * aj;
    bj = (2 * (l + j) + 1) / sqrttmp * cotK;
    Dj = bj + aj * Dj;
    if (Dj == 0.0)
      Dj = tiny;
    Cj = bj + aj / Cj;
    if (Cj == 0.0)
      Cj = tiny;
    Dj   = 1.0 / Dj;
    Delj = Cj * Dj;
    fj   = fj * Delj;
    if (Dj < 0)
      *isign *= -1;
    if (fabs(Delj - 1.0) < reltol) {
      *CF = fj;
      //printf("iter:%d, %g, %g\n",j,sqrttmp,fj);
      return _SUCCESS_;
    }
  }
  return _FAILURE_;
}

int CF1_from_Gegenbauer(int l, int beta, double sinK, double cotK, double* CF) {
  int n, alpha, k;
  double x, G, dG, Gkm1, Gkm2;
  n = beta - l - 1;
  if (n < 0)
    return _FAILURE_;
  alpha = l + 1;
  x     = sinK * cotK;  //cos(x)
  switch (n) {
    case 0:
      G  = 1;
      dG = 0;
      break;
    case 1:
      G  = 2 * alpha * x;
      dG = 2 * alpha;
      break;
    case 2:
      G  = -alpha + 2 * alpha * (1 + alpha) * x * x;
      dG = 4 * x * alpha * (1 + alpha);
      break;
    case 3:
      G  = -2 * alpha * (1 + alpha) * x + 4.0 / 3.0 * alpha * (1 + alpha) * (2 + alpha) * x * x * x;
      dG = 2 * alpha * (1 + alpha) * (2 * (2 + alpha) * x * x - 1);
      break;
    default:
      G    = 0.0;
      Gkm2 = -alpha + 2 * alpha * (1 + alpha) * x * x;
      Gkm1 = -2 * alpha * (1 + alpha) * x +
             4.0 / 3.0 * alpha * (1 + alpha) * (2 + alpha) * x * x * x;
      for (k = 4; k <= n; k++) {
        G = (2 * (k + alpha - 1) * x * Gkm1 - (k + 2 * alpha - 2) * Gkm2) / k;
        if (fabs(G) > _HYPER_OVERFLOW_) {
          Gkm2 = Gkm1 / _HYPER_OVERFLOW_;
          G    = G / _HYPER_OVERFLOW_;
          Gkm1 = G;
        }
        else {
          Gkm2 = Gkm1;
          Gkm1 = G;
        }
      }
      //Derivative. Gkm2 is in fact Gkm1..
      dG = (-n * x * G + (n + 2 * alpha - 1) * Gkm2) / (1.0 - x * x);
  }
  //%Phi = G;
  //%dPhi = l*coty.*G-siny.*dG;
  *CF = l * cotK - sinK * dG / G;
  return _SUCCESS_;
}

int hyperspherical_WKB(int K, int l, double beta, double y, double* Phi) {
  double e, w, w2, alpha, alpha2, CscK, ytp, t;
  double S      = 0.0, Q, C, argu, Ai;
  int airy_sign = 1, phisign = 1, dphisign = 1, intbeta;
  double ldbl = l;

  if (K == 1) {
    //Limit range to [0; pi/2]:
    intbeta = (int) (beta + 0.4);  //Round to nearest integer (just to be sure)
    ClosedModY(l, intbeta, &y, &phisign, &dphisign);
  }
  e     = 1.0 / sqrt(ldbl * (ldbl + 1.0));
  alpha = beta * e;
  if (K == -1) {
    CscK = 1.0 / sinh(y);
    ytp  = asinh(1.0 / alpha);
  }
  else if (K == 1) {
    CscK = 1.0 / sin(y);
    ytp  = asin(1.0 / alpha);
  }
  else {
    return _FAILURE_;
  }
  w      = alpha / CscK;
  w2     = w * w;
  alpha2 = alpha * alpha;
  if (K == -1) {
    if (y > ytp) {
      S         = alpha * log((sqrt(w2 - 1.0) + sqrt(w2 + alpha2)) / sqrt(1.0 + alpha2)) +
                  atan(1.0 / alpha * sqrt((w2 + alpha2) / (w2 - 1.0))) - 0.5 * _PI_;
      airy_sign = -1;
    }
    else {
      t         = sqrt(1.0 - w2) / sqrt(1.0 + w2 / alpha2);
      S         = atanh(t) - alpha * atan(t / alpha);
      airy_sign = 1;
    }
  }
  else if (K == 1) {
    if (y > ytp) {
      t         = sqrt(1 - w2 / alpha2) / sqrt(w2 - 1.0);
      S         = atan(t) + alpha * atan(1.0 / (t * alpha)) - 0.5 * _PI_;
      airy_sign = -1;
    }
    else {
      S         = atanh(sqrt(1.0 - w2) / sqrt(1.0 - w2 / alpha2)) -
                  alpha * log((sqrt(alpha2 - w2) + sqrt(1.0 - w2)) / sqrt(alpha2 - 1.0));
      airy_sign = 1;
    }
  }
  argu = 3.0 * S / (2.0 * e);
  Q    = CscK * CscK - alpha2;
  C    = 0.5 * sqrt(alpha) / beta;
  Ai   = airy_cheb_approx(airy_sign * pow(argu, 2.0 / 3.0));
  *Phi = phisign * 2.0 * sqrt(_PI_) * C * pow(argu, 1.0 / 6.0) * pow(fabs(Q), -0.25) * Ai * CscK;
  return _SUCCESS_;
}

double airy_cheb_approx(double z) {
  double Ai;
  if (z <= -7) {
    Ai = coef1(z);
    return Ai;
  }
  if (z <= 0) {
    Ai = coef2(z);
    return Ai;
  }
  if (z < 7) {
    Ai = coef3(z);
    return Ai;
  }
  Ai = coef4(z);
  return Ai;
}

double coef1(double z) {
  const double A[5] = {1.1282427601, -0.6803534e-4, 0.16687e-6, -0.128e-8, 0.2e-10};
  const double B[5] = {0.7822108673e-1, -0.6895649e-4, 0.32857e-6, -0.37e-8, 0.7e-10};
  double x, y, t, Ai, zeta, theta, sintheta, costheta, FA, FB;

  x        = -z;
  zeta     = _TWO_OVER_THREE_ * x * sqrt(x);
  theta    = zeta + 0.25 * _PI_;
  sintheta = sin(theta);
  costheta = cos(theta);

  y  = pow(7.0 / x, 3);  //y = (7.0/x)*(7.0/x)*(7.0/x);
  FA = cheb(y, 5, A);
  FB = cheb(y, 5, B) / zeta;

  t  = pow(x, -0.25);
  Ai = t * (sintheta * FA - costheta * FB);
  //Bi = t*(costheta*FA+sintheta*FB);
  return Ai;
}

double coef2(double z) {
  const double A[17] = {0.11535880704,
                        0.6542816649e-1,
                        0.26091774326,
                        0.21959346500,
                        0.12476382168,
                        -0.43476292594,
                        0.28357718605,
                        -0.9751797082e-1,
                        0.2182551823e-1,
                        -0.350454097e-2,
                        0.42778312e-3,
                        -0.4127564e-4,
                        0.323880e-5,
                        -0.21123e-6,
                        0.1165e-7,
                        -0.55e-9,
                        0.2e-10};
  const double B[16] = {0.10888288487,
                        -0.17511655051,
                        0.13887871101,
                        -0.11469998998,
                        0.22377807641,
                        -0.18546243714,
                        0.8063565186e-1,
                        -0.2208847864e-1,
                        0.422444527e-2,
                        -0.60131028e-3,
                        0.6653890e-4,
                        -0.590842e-5,
                        0.43131e-6,
                        -0.2638e-7,
                        0.137e-8,
                        -0.6e-10};
  //Ej = 3^(-j/3)/Gamma(j/3);
  double E1 = 0.355028053887817, E2 = 0.258819403792807;
  double x, FA, FB, Ai;
  x  = -(z / 7.0) * (z / 7.0) * (z / 7.0);
  FA = E1 * cheb(x, 17, A);
  FB = E2 * z * cheb(x, 16, B);
  Ai = FA - FB;
  //Bi = sqrt(3)*(FA+FB);
  return Ai;
}
double coef3(double z) {
  const double A[20] = {1.2974695794,    -0.20230907821,   -0.45786169516,  0.12953331987,
                        0.6983827954e-1, -0.3005148746e-1, -0.493036334e-2, 0.390425474e-2,
                        -0.1546539e-4,   -0.32193999e-3,   0.3812734e-4,    0.1714935e-4,
                        -0.416096e-5,    -0.50623e-6,      0.26346e-6,      -0.281e-8,
                        -0.1122e-7,      0.120e-8,         0.31e-9,         -0.7e-10};
  /**
double B[25]={0.47839902387,-0.6881732880e-1,0.20938146768,
            -0.3988095886e-1,0.4758441683e-1,-0.812296149e-2,
            0.462845913e-2,0.70010098e-3,-0.75611274e-3,
            0.68958657e-3,-0.33621865e-3,0.14501668e-3,-0.4766359e-4,
            0.1251965e-4,-0.193012e-5,-0.19032e-6,0.29390e-6,
            -0.13436e-6,0.4231e-7,-0.967e-8,0.135e-8,0.7e-10,
            -0.12e-9,0.4e-10,-0.1e-10};
  */
  double x, EX, EY, Ai;
  x  = z / 7.0;
  EX = exp(1.75 * z);
  EY = 1.0 / EX;

  Ai = EY * cheb(x, 20, A);
  //Bi = EX*cheb(x,25,&(B[0]));
  return Ai;
}

double coef4(double z) {
  const double A[7] =
      {0.56265126169, -0.76136219e-3, 0.765252e-5, -0.14228e-6, 0.380e-8, -0.13e-9, 0.1e-10};
  /**  double B[7]={1.1316635302,0.166141673e-02,0.1968882e-04,0.47047e-06,
            0.1769e-7,0.94e-9,0.6e-10};
  */
  double x, Y, t, zeta, EX, EY, Ai;

  Y    = z * sqrt(z);
  zeta = 2.0 / 3.0 * Y;
  EX   = exp(zeta);
  EY   = 1.0 / EX;
  x    = 7 * sqrt(7) / Y;
  t    = pow(z, -0.25);

  Ai = t * EY * cheb(x, 7, A);
  //Bi = t*EX*cheb(x, 7, &(B[0]));
  return Ai;
}

double cheb(double x, int n, const double A[]) {
  double b, d, u, y, c, F;
  int j;
  b = 0.0;
  d = A[n - 1];
  u = 2 * x - 1.0;
  y = 2 * u;

  for (j = (n - 1); j > 1; j--) {
    c = b;
    b = d;
    d = y * b - c + A[j - 1];
  }
  F = u * d - b + 0.5 * A[0];
  return F;
}

int ClosedModY(int l, int beta, double* y, int* phisign, int* dphisign) {
  *phisign  = 1;
  *dphisign = 1;

  while (*y > _TWOPI_)
    *y -= _TWOPI_;

  if ((*y) > _PI_) {
    *y = 2.0 * _PI_ - (*y);
    //phisign *= pow(-1,l)
    if (l % 2 == 1)  //l is odd
      *phisign = -*phisign;
    else  //l is even
      *dphisign = -*dphisign;
  }
  if ((*y) > 0.5 * _PI_) {
    *y = _PI_ - (*y);
    //phisign *= pow(-1,beta-l-1)
    if ((beta - l) % 2 == 0)  //beta-l-1 odd
      *phisign = -*phisign;
    else  //beta-l-1 even
      *dphisign = -*dphisign;
  }
  return _SUCCESS_;
}

void hyperspherical_get_xmin_from_Airy(
    int K, int l, double beta, double xtol, double phiminabs, double* xmin, int* fevals) {
  double xold, xtp = 0, xleft, xright, xnew;
  double Fnew, Fold, Fleft, Fright;
  double delx, lambda;
  double AIRY_SAFETY = 1e-6;
  struct WKB_parameters wkbstruct;
  //Start searching from turning point:
  switch (K) {
    case -1:
      xtp = asinh(sqrt(l * (l + 1.)) / beta);
      break;
    case 0:
      xtp = sqrt(l * (l + 1.)) / beta;
      break;
    case 1:
      xtp = asin(sqrt(l * (l + 1.)) / beta);
      break;
  }
  wkbstruct.K         = K;
  wkbstruct.l         = l;
  wkbstruct.beta      = beta;
  wkbstruct.phiminabs = phiminabs;

  xnew = 0.99 * xtp;

  Fnew    = PhiWKB_minus_phiminabs(xnew, &wkbstruct);
  *fevals = (*fevals) + 1;

  lambda = 2 * _PI_ / (beta + 5.0);
  if (Fnew > 0)
    delx = -lambda;
  else
    delx = 0.25 * lambda;

  do {
    //printf("In the loop: xnew = %g, Fnew=%g, Fold=%g\n",xnew,Fnew,Fold);
    xold  = xnew;
    Fold  = Fnew;
    xnew += delx;
    if (xnew < AIRY_SAFETY) {
      xnew    = AIRY_SAFETY;
      Fnew    = PhiWKB_minus_phiminabs(xnew, &wkbstruct);
      *fevals = (*fevals) + 1;
      if (Fnew >= 0.0) {
        *xmin = xnew;
        return;
      }
      else {
        break;
      }
    }
    Fnew    = PhiWKB_minus_phiminabs(xnew, &wkbstruct);
    *fevals = (*fevals) + 1;
  } while (std::copysign(1.0, Fnew) == std::copysign(1.0, Fold));

  if (Fnew <= 0.0) {
    xleft  = xnew;
    Fleft  = Fnew;
    xright = xold;
    Fright = Fold;
  }
  else {
    xleft  = xold;
    Fleft  = Fold;
    xright = xnew;
    Fright = Fnew;
  }

  fzero_ridder(PhiWKB_minus_phiminabs,
               xleft,
               xright,
               xtol,
               &wkbstruct,
               &Fleft,
               &Fright,
               xmin,
               fevals);
}

double PhiWKB_minus_phiminabs(double x, void* param) {
  double phiwkb                   = 0.;
  struct WKB_parameters* wkbparam = (struct WKB_parameters*) param;
  hyperspherical_WKB(wkbparam->K, wkbparam->l, wkbparam->beta, x, &phiwkb);
  return (fabs(phiwkb) - wkbparam->phiminabs);
}

int fzero_ridder(double (*func)(double, void*),
                 double x1,
                 double x2,
                 double xtol,
                 void* param,
                 double* Fx1,
                 double* Fx2,
                 double* xzero,
                 int* fevals) {
  /**Using Ridders' method, return the root of a function func known to
      lie between x1 and x2. The root, returned as zriddr, will be found to
      an approximate accuracy xtol.
   */
  int j, MAXIT = 1000;
  double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;
  if ((Fx1 != nullptr) && (Fx2 != nullptr)) {
    fl = *Fx1;
    fh = *Fx2;
  }
  else {
    fl      = (*func)(x1, param);
    fh      = (*func)(x2, param);
    *fevals = (*fevals) + 2;
  }
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl  = x1;
    xh  = x2;
    ans = -1.11e11;
    for (j = 1; j <= MAXIT; j++) {
      xm      = 0.5 * (xl + xh);
      fm      = (*func)(xm, param);
      *fevals = (*fevals) + 1;
      s       = sqrt(fm * fm - fl * fh);
      if (s == 0.0) {
        *xzero = ans;
        //printf("Success 1\n");
        return _SUCCESS_;
      }
      xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
      if (fabs(xnew - ans) <= xtol) {
        *xzero = ans;
        return _SUCCESS_;
      }
      ans     = xnew;
      fnew    = (*func)(ans, param);
      *fevals = (*fevals) + 1;
      if (fnew == 0.0) {
        *xzero = ans;
        //printf("Success 2, ans=%g\n",ans);
        return _SUCCESS_;
      }
      if (std::copysign(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      }
      else if (std::copysign(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      }
      else if (std::copysign(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      }
      else
        return _FAILURE_;
      if (fabs(xh - xl) <= xtol) {
        *xzero = ans;
        //        printf("Success 3\n");
        return _SUCCESS_;
      }
    }
    fprintf(stderr, "zriddr exceed maximum iterations");
    return _FAILURE_;
  }
  else {
    if (fl == 0.0)
      return x1;
    if (fh == 0.0)
      return x2;
    fprintf(stderr, "root must be bracketed in zriddr.");
    return _FAILURE_;
  }
  return _FAILURE_;
}

void hyperspherical_get_xmin_from_approx(
    int K, int l, double nu, double ignore1, double phiminabs, double* xmin, int* ignore2) {
  double l_plus_half;
  double lhs;
  double alpha;
  double ldbl = l;
  double x;

  l_plus_half = l + 0.5;
  lhs         = 1.0 / l_plus_half * log(2 * phiminabs * l_plus_half);
  //Using Chebyshev cubic root, much cleaner:
  alpha = -2.0 * lhs / 5.0 *
          (1.0 + 2.0 * cosh(1.0 / 3.0 * acosh(1.0 + 375.0 / (16.0 * lhs * lhs))));
  x     = l_plus_half / cosh(alpha) / nu;
  if (K == -1) {
    //%Correct for open case:
    x *= asinh(ldbl / nu) / (ldbl / nu);
    //...and fudge for small nu:
    x *= ((nu + 0.4567) / (nu + 1.24) - 2.209e-3);
  }
  else if (K == 1) {
    //Correct for closed case if possible
    x *= asin(ldbl / nu) / (ldbl / nu);
  }
  *xmin = x;
}

int hyperspherical_bessel_direct_vector(
    int K, double beta, int* lvec, int nl, double* xvec, int nx, double* Phi) {
  /** Evaluate Phi_l^beta(x) by direct forwards/backwards recurrence at each
      requested x, for every requested l. Keeps the numerical considerations
      that are otherwise spread across hyperspherical_HIS_create:
        - closed-case (K=1) symmetry fold into [0, pi/2] via ClosedModY,
        - turning-point stability (backwards below xfwd, forwards above),
        - maximum l in the closed case (l < beta).
      Phi is row-major: Phi[il*nx + ix]. */
  int il, ix, l, lmax, intbeta = 0;
  double beta2 = beta * beta, xfwd, folded_y, sinK, cotK;
  int phisign, dphisign;

  class_test((K != 0) && (K != 1) && (K != -1), "K must be -1, 0, or 1 (got %d)", K);
  class_test(beta <= 0.0, "beta must be positive (got %g)", beta);

  lmax = 0;
  for (il = 0; il < nl; il++) {
    class_test(lvec[il] < 0, "l must be non-negative (got %d)", lvec[il]);
    if (lvec[il] > lmax)
      lmax = lvec[il];
  }

  if (K == 1) {
    intbeta = (int) (beta + 0.2);
    class_test(fabs(beta - intbeta) > 1e-6,
               "closed case (K=1) requires integer beta (got %g)",
               beta);
    class_test(lmax >= intbeta,
               "closed case requires l < beta; got max l=%d, beta=%d",
               lmax,
               intbeta);
  }

  // The forwards recurrence unconditionally computes PhiL[1] (using
  // one_over_sqrtK[1]) even when lmax == 0, so the arrays need at least 2
  // entries. Beyond that, only sqrtK[0..lmax] are read. Sizing to
  // max(lmax,1)+1 (rather than lmax+2/+3) keeps the fill loop below beta in the
  // closed case (lmax = beta-1), avoiding sqrt(beta^2 - l^2) for l >= beta which
  // would be a NaN.
  int nsqrtK = (lmax >= 1) ? lmax : 1;
  std::vector<double> sqrtK(nsqrtK + 1), one_over_sqrtK(nsqrtK + 1);
  for (l = 0; l <= nsqrtK; l++) {
    if (K == 0)
      sqrtK[l] = beta;
    else if (K == 1)
      sqrtK[l] = sqrt(beta2 - (double) l * l);
    else
      sqrtK[l] = sqrt(beta2 + (double) l * l);
    one_over_sqrtK[l] = 1.0 / sqrtK[l];
  }

  if (K == 0)
    xfwd = sqrt(lmax * (lmax + 1.0)) / beta;
  else if (K == 1)
    xfwd = asin(sqrt(lmax * (lmax + 1.0)) / beta);
  else
    xfwd = asinh(sqrt(lmax * (lmax + 1.0)) / beta);

  std::vector<double> PhiL(lmax + 2);

  for (ix = 0; ix < nx; ix++) {
    folded_y = xvec[ix];
    if (K == 1) {
      /* Fold y into [0, pi/2]. The geometric fold is l-independent; the l
         argument to ClosedModY only controls its sign outputs, which are
         discarded here and recomputed per-l below. */
      ClosedModY(0, intbeta, &folded_y, &phisign, &dphisign);
    }

    if (K == 0) {
      sinK = folded_y;
      cotK = 1.0 / folded_y;
    }
    else if (K == 1) {
      sinK = sin(folded_y);
      cotK = 1.0 / tan(folded_y);
    }
    else {
      sinK = sinh(folded_y);
      cotK = 1.0 / tanh(folded_y);
    }

    if (folded_y < xfwd)
      hyperspherical_backwards_recurrence(K,
                                          lmax,
                                          beta,
                                          folded_y,
                                          sinK,
                                          cotK,
                                          sqrtK.data(),
                                          one_over_sqrtK.data(),
                                          PhiL.data());
    else
      hyperspherical_forwards_recurrence(K,
                                         lmax,
                                         beta,
                                         folded_y,
                                         sinK,
                                         cotK,
                                         sqrtK.data(),
                                         one_over_sqrtK.data(),
                                         PhiL.data());

    for (il = 0; il < nl; il++) {
      l           = lvec[il];
      double sign = 1.0;
      if (K == 1) {
        double tmp = xvec[ix];
        phisign    = 1;
        dphisign   = 1;
        ClosedModY(l, intbeta, &tmp, &phisign, &dphisign);
        sign = phisign;
      }
      Phi[il * nx + ix] = sign * PhiL[l];
    }
  }
  return _SUCCESS_;
}

template <int Order, bool DoPhi, bool DoDPhi, bool DoD2Phi>
void hyperspherical_Hermite_interpolation(const HyperInterpStruct* pHIS,
                                          int nxi,
                                          int lnum,
                                          const double* xinterp,
                                          double* Phi,
                                          double* dPhi,
                                          double* d2Phi) {
  static_assert(Order == 4 || Order == 6, "implemented interpolation orders");
  static_assert(DoPhi || DoDPhi || DoD2Phi, "at least one output required");
  // Highest derivative level of the tabulated function entering the fit.
  // Everything below stays in named scalars so the hot loop lives entirely in
  // registers; do not "simplify" to arrays or level-generic helpers (interior
  // pointers defeat scalar replacement and cost ~7% of Transfer).
  constexpr int kDepth = (DoD2Phi ? 2 : DoDPhi ? 1 : 0) + Order / 2 - 1;

  const double* xvec   = pHIS->x.data();
  const double* sinKv  = pHIS->sinK.data();
  const double* cotKv  = pHIS->cotK.data();
  const double beta    = pHIS->beta;
  const double beta2   = beta * beta;
  const double deltax  = pHIS->delta_x;
  const double deltax2 = deltax * deltax;
  const int K          = pHIS->K;
  const int nx         = pHIS->x_size;
  const double* Phi_l  = pHIS->phi.data() + lnum * nx;
  const double* dPhi_l = pHIS->dphi.data() + lnum * nx;
  const int l          = pHIS->l[lnum];
  const double lxlp1   = l * (l + 1.0);
  const double xmin    = xvec[0];
  const double xmax    = xvec[nx - 1];

  // Derivative levels 0..kDepth at the left (..m) and right (..p) border of
  // the current grid interval, and the Hermite coefficients per output
  // (a: Phi, b: dPhi, c: d2Phi). Unused ones are optimized away.
  double ym = 0, dym = 0, d2ym = 0, d3ym = 0, d4ym = 0;
  double yp = 0, dyp = 0, d2yp = 0, d3yp = 0, d4yp = 0;
  double a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
  double b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0;
  double c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0;
  double left_border = xmax, right_border = xmin, next_border = xmin;
  int current_border_idx = 0;
  int phisign = 1, dphisign = 1;

  for (int j = 0; j < nxi; j++) {
    double x = xinterp[j];
    //take advantage of periodicity of functions in closed case
    if (pHIS->K == 1)
      ClosedModY(l, (int) (pHIS->beta + 0.2), &x, &phisign, &dphisign);
    if ((x < xmin) || (x > xmax)) {
      //Outside interpolation region, set to zero.
      if constexpr (DoPhi)
        Phi[j] = 0.0;
      if constexpr (DoDPhi)
        dPhi[j] = 0.0;
      if constexpr (DoD2Phi)
        d2Phi[j] = 0.0;
      continue;
    }
    if ((x > right_border) || (x < left_border)) {
      if ((x > next_border) || (x < left_border)) {
        current_border_idx = ((int) ((x - xmin) / deltax)) + 1;
        //max operation takes care of case x = xmin,
        //min operation takes care of case x = xmax.
        current_border_idx = std::max(1, current_border_idx);
        current_border_idx = std::min(nx - 1, current_border_idx);
        //Calculate left derivatives via the differentiated radial equation:
        ym  = Phi_l[current_border_idx - 1];
        dym = dPhi_l[current_border_idx - 1];
        if constexpr (kDepth >= 2) {
          const double cotKm  = cotKv[current_border_idx - 1];
          const double sinKm  = sinKv[current_border_idx - 1];
          const double sinKm2 = sinKm * sinKm;
          d2ym                = -2 * dym * cotKm + ym * (lxlp1 / sinKm2 - beta2 + K);
          if constexpr (kDepth >= 3)
            d3ym = -2 * cotKm * d2ym - 2 * ym * lxlp1 * cotKm / sinKm2 +
                   dym * (K - beta2 + (2 + lxlp1) / sinKm2);
          if constexpr (kDepth >= 4)
            d4ym = -2 * cotKm * d3ym + d2ym * (K - beta2 + (4 + lxlp1) / sinKm2) +
                   dym * (-4 * (1 + lxlp1) * cotKm / sinKm2) +
                   ym * (2 * lxlp1 / sinKm2 * (2 * cotKm * cotKm + 1 / sinKm2));
        }
      }
      else {
        //x>current_border but not next border: I have moved to next block.
        current_border_idx++;
        //Copy former right derivatives to left derivatives.
        ym  = yp;
        dym = dyp;
        if constexpr (kDepth >= 2)
          d2ym = d2yp;
        if constexpr (kDepth >= 3)
          d3ym = d3yp;
        if constexpr (kDepth >= 4)
          d4ym = d4yp;
      }
      left_border  = xvec[std::max(0, current_border_idx - 1)];
      right_border = xvec[current_border_idx];
      next_border  = xvec[std::min(nx - 1, current_border_idx + 1)];
      //Evaluate right derivatives:
      yp  = Phi_l[current_border_idx];
      dyp = dPhi_l[current_border_idx];
      if constexpr (kDepth >= 2) {
        const double cotKp  = cotKv[current_border_idx];
        const double sinKp  = sinKv[current_border_idx];
        const double sinKp2 = sinKp * sinKp;
        d2yp                = -2 * dyp * cotKp + yp * (lxlp1 / sinKp2 - beta2 + K);
        if constexpr (kDepth >= 3)
          d3yp = -2 * cotKp * d2yp - 2 * yp * lxlp1 * cotKp / sinKp2 +
                 dyp * (K - beta2 + (2 + lxlp1) / sinKp2);
        if constexpr (kDepth >= 4)
          d4yp = -2 * cotKp * d3yp + d2yp * (K - beta2 + (4 + lxlp1) / sinKp2) +
                 dyp * (-4 * (1 + lxlp1) * cotKp / sinKp2) +
                 yp * (2 * lxlp1 / sinKp2 * (2 * cotKp * cotKp + 1 / sinKp2));
      }
      //Calculate coefficients:
      if constexpr (Order == 4) {
        if constexpr (DoPhi) {
          a1 = dym * deltax;
          a2 = -2 * dym * deltax - dyp * deltax - 3 * ym + 3 * yp;
          a3 = dym * deltax + dyp * deltax + 2 * ym - 2 * yp;
        }
        if constexpr (DoDPhi) {
          b1 = d2ym * deltax;
          b2 = -2 * d2ym * deltax - d2yp * deltax - 3 * dym + 3 * dyp;
          b3 = d2ym * deltax + d2yp * deltax + 2 * dym - 2 * dyp;
        }
        if constexpr (DoD2Phi) {
          c1 = d3ym * deltax;
          c2 = -2 * d3ym * deltax - d3yp * deltax - 3 * d2ym + 3 * d2yp;
          c3 = d3ym * deltax + d3yp * deltax + 2 * d2ym - 2 * d2yp;
        }
      }
      else {
        if constexpr (DoPhi) {
          a1 = dym * deltax;
          a2 = 0.5 * d2ym * deltax2;
          a3 = (-1.5 * d2ym + 0.5 * d2yp) * deltax2 - (6 * dym + 4 * dyp) * deltax - 10 * (ym - yp);
          a4 = (1.5 * d2ym - d2yp) * deltax2 + (8 * dym + 7 * dyp) * deltax + 15 * (ym - yp);
          a5 = (-0.5 * d2ym + 0.5 * d2yp) * deltax2 - 3 * (dym + dyp) * deltax - 6 * (ym - yp);
        }
        if constexpr (DoDPhi) {
          b1 = d2ym * deltax;
          b2 = 0.5 * d3ym * deltax2;
          b3 = (-1.5 * d3ym + 0.5 * d3yp) * deltax2 - (6 * d2ym + 4 * d2yp) * deltax -
               10 * (dym - dyp);
          b4 = (1.5 * d3ym - d3yp) * deltax2 + (8 * d2ym + 7 * d2yp) * deltax + 15 * (dym - dyp);
          b5 = (-0.5 * d3ym + 0.5 * d3yp) * deltax2 - 3 * (d2ym + d2yp) * deltax - 6 * (dym - dyp);
        }
        if constexpr (DoD2Phi) {
          c1 = d3ym * deltax;
          c2 = 0.5 * d4ym * deltax2;
          c3 = (-1.5 * d4ym + 0.5 * d4yp) * deltax2 - (6 * d3ym + 4 * d3yp) * deltax -
               10 * (d2ym - d2yp);
          c4 = (1.5 * d4ym - d4yp) * deltax2 + (8 * d3ym + 7 * d3yp) * deltax + 15 * (d2ym - d2yp);
          c5 = (-0.5 * d4ym + 0.5 * d4yp) * deltax2 - 3 * (d3ym + d3yp) * deltax -
               6 * (d2ym - d2yp);
        }
      }
    }
    //Evaluate polynomial:
    const double z  = (x - left_border) / deltax;
    const double z2 = z * z;
    const double z3 = z2 * z;
    if constexpr (Order == 4) {
      if constexpr (DoPhi)
        Phi[j] = (ym + a1 * z + a2 * z2 + a3 * z3) * phisign;
      if constexpr (DoDPhi)
        dPhi[j] = (dym + b1 * z + b2 * z2 + b3 * z3) * dphisign;
      if constexpr (DoD2Phi)
        d2Phi[j] = (d2ym + c1 * z + c2 * z2 + c3 * z3) * phisign;
    }
    else {
      const double z4 = z2 * z2;
      const double z5 = z2 * z3;
      if constexpr (DoPhi)
        Phi[j] = (ym + a1 * z + a2 * z2 + a3 * z3 + a4 * z4 + a5 * z5) * phisign;
      if constexpr (DoDPhi)
        dPhi[j] = (dym + b1 * z + b2 * z2 + b3 * z3 + b4 * z4 + b5 * z5) * dphisign;
      if constexpr (DoD2Phi)
        d2Phi[j] = (d2ym + c1 * z + c2 * z2 + c3 * z3 + c4 * z4 + c5 * z5) * phisign;
    }
  }
}

template void hyperspherical_Hermite_interpolation<4, true, false, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<4, false, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<4, true, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<4, true, false, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<4, true, true, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<6, true, false, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<6, false, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<6, true, true, false>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<6, true, false, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
template void hyperspherical_Hermite_interpolation<6, true, true, true>(
    const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
