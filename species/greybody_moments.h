#pragma once

namespace greybody {

/** Riemann zeta via the Borwein eta (alternating-series) algorithm. Vendored
 *  because libc++ lacks std::riemann_zeta; used uniformly on all compilers so
 *  results are toolchain-independent. Valid for s != 1. */
double riemann_zeta(double s);

/** Grey-body PSD parameters. Construct either directly from (alpha, x, q0) or
 *  by inverting the velocity moments (r, M2, M3). Stores the log-space
 *  combinations used by the PSD and quadrature. */
class GreyBodyParams {
 public:
  static GreyBodyParams FromDirect(double alpha, double x, double q0);
  static GreyBodyParams FromMoments(double r, double M2, double M3);

  double alpha() const {
    return alpha_;
  }
  double x() const {
    return x_;
  }
  double q0() const {
    return q0_;
  }
  double x_times_alpha() const {
    return x_times_alpha_;
  }
  double alpham1_logq0() const {
    return alpham1_logq0_;
  }

  /** Closed-form dimensionless moments of this distribution. */
  void moments(double& M2, double& M3, double& M4) const;

 private:
  GreyBodyParams()      = default;
  double alpha_         = 1.;
  double x_             = 1.;
  double q0_            = 1.;
  double x_times_alpha_ = 1.;
  double alpham1_logq0_ = 0.;
};

}  // namespace greybody
