#include "bbn_weak.h"

#include <cmath>

#include "arrays.h"
#include "bbn_nuclides.h"
#include "errors.h"
#include "quadrature.h"

namespace {

/** Fine-structure constant, for the Coulomb correction. */
constexpr double kAlphaEM = 7.2973525693e-3;

/** Neutron-proton mass difference in units of the electron mass. Built from the
 *  NUCLEAR masses in bbn_nuclides.h -- atomic masses would give the
 *  neutron-to-hydrogen difference, which is wrong by one electron mass here. */
const double kQ = kDeltaMnpMeV / kElectronMassMeV;

/** Fermi-Dirac occupation at scaled energy, written to underflow rather than
 *  overflow: the argument is always non-negative. */
inline double Occupation(double energy_over_T) {
  const double e = std::exp(-energy_over_T);
  return e / (1.0 + e);
}

/** beta * F(eps), where beta = sqrt(eps^2 - 1) and F is the non-relativistic
 *  Fermi function for a lepton of charge `sign` in the field of the proton.
 *
 *  Returned as a PRODUCT because F alone diverges as 1/beta at threshold while
 *  the phase-space factor beta vanishes there; the product is finite and smooth.
 *  Forming it directly avoids both the 0/0 and the overflow that a separate
 *  F would hit for a repelled positron near threshold. */
double BetaTimesFermi(double eps, int sign) {
  const double beta = std::sqrt(std::max(eps * eps - 1.0, 0.0));
  const double y    = 2.0 * M_PI * kAlphaEM * eps; /* = 2 pi eta * beta */

  if (sign > 0) {
    /* Attractive (electron): beta*F = y / (1 - exp(-y/beta)).
       As beta -> 0 the exponential vanishes and this tends to y, finite. */
    if (beta <= 0.0) {
      return y;
    }
    return y / (1.0 - std::exp(-y / beta));
  }

  /* Repulsive (positron): beta*F = y * exp(-y/beta) / (1 - exp(-y/beta)),
     which tends to 0 at threshold -- the Coulomb barrier -- and to beta at
     high energy. Written with exp(-y/beta) so it never overflows. */
  if (beta <= 0.0) {
    return 0.0;
  }
  const double e = std::exp(-y / beta);
  return y * e / (1.0 - e);
}

}  // namespace

BbnWeakRates::BbnWeakRates(double tau_n_seconds, int nodes) : tau_n_(tau_n_seconds) {
  /* Checked BEFORE anything is sized: a negative node count reaches the vector
     constructor as an enormous size_t, and a non-positive lifetime would divide
     into the normalization below and hand back rates of the wrong sign rather
     than failing. Neither is a value a sampler varies -- one is an API argument,
     the other is nonsense rather than an unlikely cosmology -- so both abort. */
  class_test_severe(nodes < 8, "BbnWeakRates needs at least 8 quadrature nodes, got %d", nodes);
  class_test_severe(!(tau_n_seconds > 0.),
                    "BbnWeakRates needs a positive neutron lifetime, got %g s",
                    tau_n_seconds);

  x_.resize(nodes);
  w_.resize(nodes);
  std::vector<double> mu(nodes), weight(nodes);
  quadrature_gauss_legendre(mu.data(), weight.data(), nodes, 1.0e-14);
  for (int i = 0; i < nodes; ++i) {
    x_[i] = 0.5 * (mu[i] + 1.0);
    w_[i] = 0.5 * weight[i];
  }

  /* Normalization. In vacuum every occupation vanishes and only the decay
     channel survives, so
         1/tau_n = K * Int_1^q d(eps) eps * beta F(eps) * (q - eps)^2 .
     Solving for K absorbs G_F, V_ud, g_A and the (1 + 3 g_A^2) factor into one
     measured number, which is why tau_n is the only weak-sector input. */
  const double decay_integral = Integrate(1.0, kQ, [](double eps) {
    const double gap = kQ - eps;
    return eps * BetaTimesFermi(eps, +1) * gap * gap;
  });
  norm_                       = 1.0 / (tau_n_ * decay_integral);
}

void BbnWeakRates::Rates(double T_gamma_mev,
                         double T_nu_mev,
                         double* lambda_np,
                         double* lambda_pn) const {
  /* Inverse temperatures in units of the electron mass. */
  const double z    = kElectronMassMeV / T_gamma_mev;
  const double z_nu = kElectronMassMeV / T_nu_mev;

  /* Occupations. f_e is evaluated at the electron/positron energy in T_gamma;
     f_nu at the neutrino energy in T_nu. Keeping these as two separate lambdas
     makes it impossible to accidentally pass one temperature for the other. */
  const auto f_e  = [z](double eps) { return Occupation(eps * z); };
  const auto f_nu = [z_nu](double w) { return Occupation(w * z_nu); };

  /* Each channel is integrated in the variable that makes its exponential decay
     have unit scale, so one fixed Gauss-Legendre rule serves all of them across
     the whole temperature range. `kTail` is how many decay lengths to keep. */
  constexpr double kTail = 60.0;
  const double nu_span   = kTail / z_nu;
  const double e_span    = kTail / z;

  /* ---- n -> p ------------------------------------------------------------ */

  /* (1) n + nu_e -> p + e-. Integrated over the neutrino energy w = eps - q. */
  const double i1 = Integrate(0.0, nu_span, [&](double w) {
    const double eps = w + kQ;
    return w * w * f_nu(w) * eps * BetaTimesFermi(eps, +1) * (1.0 - f_e(eps));
  });

  /* (2) n + e+ -> p + nubar_e. No Coulomb factor: the positron shares its state
     with a neutron, not with the proton. Integrated over the positron energy. */
  const double i2 = Integrate(1.0, 1.0 + e_span, [&](double eps) {
    const double w = eps + kQ;
    return eps * BetaTimesFermi(eps, -1) * f_e(eps) * w * w * (1.0 - f_nu(w));
  });

  /* (3) n -> p + e- + nubar_e, the decay channel, over the finite Dalitz range. */
  const double i3 = Integrate(1.0, kQ, [&](double eps) {
    const double w = kQ - eps;
    return eps * BetaTimesFermi(eps, +1) * (1.0 - f_e(eps)) * w * w * (1.0 - f_nu(w));
  });

  /* ---- p -> n : the exact reverses ---------------------------------------- */

  /* (1') p + e- -> n + nu_e. */
  const double j1 = Integrate(0.0, e_span, [&](double y) {
    const double eps = kQ + y;
    const double w   = y;
    return eps * BetaTimesFermi(eps, +1) * f_e(eps) * w * w * (1.0 - f_nu(w));
  });

  /* (2') p + nubar_e -> n + e+. */
  const double j2 = Integrate(1.0, 1.0 + nu_span, [&](double eps) {
    const double w = eps + kQ;
    return eps * BetaTimesFermi(eps, -1) * (1.0 - f_e(eps)) * w * w * f_nu(w);
  });

  /* (3') p + e- + nubar_e -> n, the inverse of decay. Tiny, but it is what makes
     the two directions satisfy detailed balance exactly at a common temperature,
     which bbn_weak_test checks. */
  const double j3 = Integrate(1.0, kQ, [&](double eps) {
    const double w = kQ - eps;
    return eps * BetaTimesFermi(eps, +1) * f_e(eps) * w * w * f_nu(w);
  });

  *lambda_np = norm_ * (i1 + i2 + i3);
  *lambda_pn = norm_ * (j1 + j2 + j3);
}

BbnWeakTable::BbnWeakTable(const BbnWeakRates& exact,
                           const std::vector<double>& tau,
                           const std::vector<double>& T_gamma,
                           const std::vector<double>& T_nu) {
  class_test_severe(tau.size() < 4 || tau.size() != T_gamma.size() || tau.size() != T_nu.size(),
                    "BBN weak-rate table needs at least 4 matching trajectory samples, got "
                    "%zu/%zu/%zu",
                    tau.size(),
                    T_gamma.size(),
                    T_nu.size());

  n_         = static_cast<int>(tau.size());
  tau_first_ = tau.front();
  tau_step_  = (tau.back() - tau.front()) / (n_ - 1);
  class_test_severe(!(tau_step_ > 0.), "BBN weak-rate table needs an increasing tau grid");

  ln_rates_.resize(static_cast<std::size_t>(n_) * 2);
  ln_second_.resize(ln_rates_.size());

  for (int i = 0; i < n_; ++i) {
    double np = 0., pn = 0.;
    exact.Rates(T_gamma[i], T_nu[i], &np, &pn);
    /* Clamped rather than guarded with an if: log(0) is -inf, and these
       translation units are built with -ffinite-math-only, where materialising
       an infinity is undefined rather than merely large. */
    ln_rates_[2 * i]     = np > 0. ? std::max(std::log(np), kLogFloor) : kLogFloor;
    ln_rates_[2 * i + 1] = pn > 0. ? std::max(std::log(pn), kLogFloor) : kLogFloor;
  }

  array_spline_table_lines(tau.data(),
                           n_,
                           ln_rates_.data(),
                           2,
                           ln_second_.data(),
                           _SPLINE_EST_DERIV_);
}

void BbnWeakTable::Rates(double tau, double* lambda_np, double* lambda_pn) const {
  /* The grid is uniform by construction, so the bracketing interval is an index
     computation rather than a search -- which is most of why this is faster than
     evaluating the integrals. */
  double position = (tau - tau_first_) / tau_step_;
  if (position < 0.) {
    position = 0.;
  }
  else if (position > n_ - 1) {
    position = n_ - 1;
  }

  int lo = static_cast<int>(position);
  if (lo > n_ - 2) {
    lo = n_ - 2;
  }

  const double b  = position - lo;
  const double a  = 1.0 - b;
  const double h2 = tau_step_ * tau_step_ / 6.0;
  const double ca = (a * a * a - a) * h2;
  const double cb = (b * b * b - b) * h2;

  const std::size_t i0 = static_cast<std::size_t>(lo) * 2;
  const std::size_t i1 = i0 + 2;

  const double ln_np = a * ln_rates_[i0] + b * ln_rates_[i1] + ca * ln_second_[i0] +
                       cb * ln_second_[i1];
  const double ln_pn = a * ln_rates_[i0 + 1] + b * ln_rates_[i1 + 1] + ca * ln_second_[i0 + 1] +
                       cb * ln_second_[i1 + 1];

  *lambda_np = std::exp(ln_np);
  *lambda_pn = std::exp(ln_pn);
}
