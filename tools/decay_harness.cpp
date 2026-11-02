#include "decay_harness.h"

#include <algorithm>
#include <cmath>

#include "common.h"
#include "errors.h"
#include "evolver_ndf15.h"
#include "evolver_options.h"

namespace {
constexpr double kGH   = 2.;  // parent spin states
constexpr double kGL   = 2.;  // fermion daughter
constexpr double kGPhi = 1.;  // scalar
}  // namespace

DecayHarness::Grid DecayHarness::LogGrid(double qmin, double qmax, int n, bool floor_to_zero) {
  // Same convention as the species' daughter grids: log-spaced nodes with
  // trapezoid weights, half cells at the ends.
  class_test_severe(!(qmin > 0.) || !(qmax > qmin) || n < 3,
                    "decay harness grid needs 0 < q_min < q_max and n >= 3 (got %g, %g, %d)",
                    qmin,
                    qmax,
                    n);
  Grid g;
  g.q.resize(n);
  g.dq.resize(n);
  const double h = std::log(qmax / qmin) / (n - 1);
  for (int i = 0; i < n; ++i) {
    g.q[i] = qmin * std::exp(i * h);
  }
  g.dq[0] = 0.5 * (g.q[1] - g.q[0]);
  for (int i = 1; i < n - 1; ++i) {
    g.dq[i] = 0.5 * (g.q[i + 1] - g.q[i - 1]);
  }
  g.dq[n - 1] = 0.5 * (g.q[n - 1] - g.q[n - 2]);
  if (floor_to_zero) {
    const double edge = 0.5 * (g.q[0] + g.q[1]);
    g.dq[0]           = edge * edge * edge / (3. * g.q[0] * g.q[0]);
  }
  return g;
}

DecayHarness::DecayHarness(Config cfg) : cfg_(cfg) {
  class_test_severe(!(cfg_.m_H_mev > 0.) || !(cfg_.tau_s > 0.),
                    "decay harness needs a positive mass and lifetime (got m = %g MeV, "
                    "tau = %g s)",
                    cfg_.m_H_mev,
                    cfg_.tau_s);
  class_test_severe(!(cfg_.T_start_mev > cfg_.T_end_mev) || !(cfg_.T_end_mev > 0.),
                    "decay harness runs from T_start down to T_end > 0 (got %g -> %g MeV)",
                    cfg_.T_start_mev,
                    cfg_.T_end_mev);
  gamma_s_ = 1. / cfg_.tau_s;

  P_ = LogGrid(cfg_.parent_q_min, cfg_.parent_q_max, cfg_.parent_n);
  // Both daughters share one grid. The network's emission saturates once the parent is
  // gone, but MaxDaughterMomentum is a bound at any a, so it is evaluated at the latest
  // one; the log grid makes the extra decades cheap.
  double qd_max = cfg_.daughter_q_max;
  if (!(qd_max > 0.)) {
    qd_max = 1.05 * DecayTransitionKernel::MaxDaughterMomentum(P_.q.back(),
                                                               cfg_.m_H_mev / cfg_.T_end_mev,
                                                               cfg_.r_l,
                                                               cfg_.r_phi);
  }
  D_ = LogGrid(cfg_.daughter_q_min, qd_max, cfg_.daughter_n, true);

  DecayTransitionKernel::Config k = cfg_.kernel;
  k.fermion_mass_ratio            = cfg_.r_l;
  k.boson_mass_ratio              = cfg_.r_phi;
  kernel_                         = std::make_unique<DecayTransitionKernel>(View(P_),
                                                                            View(D_),
                                                                            View(D_),
                                                                            Statistics::Fermion,
                                                                            Statistics::Boson,
                                                                            k);

  nP_ = static_cast<int>(P_.q.size());
  nD_ = static_cast<int>(D_.q.size());
  iT_ = nP_ + 2 * nD_;
  it_ = iT_ + 1;
  iW_ = iT_ + 2;
  ny_ = iT_ + 3;
}

double DecayHarness::BareNumberParent(const double* f) const {
  double N = 0.;
  for (int i = 0; i < nP_; ++i) {
    N += P_.dq[i] * P_.q[i] * P_.q[i] * f[i];
  }
  return N;
}

double DecayHarness::BareNumberDaughter(const double* f) const {
  double N = 0.;
  for (int j = 0; j < nD_; ++j) {
    N += D_.dq[j] * D_.q[j] * D_.q[j] * f[j];
  }
  return N;
}

void DecayHarness::SectorEnergy(
    double a, const double* y, double& rho, double& Ec, double& dEc) const {
  // E_c = a^4 rho = Sum_i g_i/(2 pi^2) Sum q^2 dq eps_i f_i, all dimensionless, and its
  // expansion source d eps/dln a = (a m)^2 / eps per particle.
  const double inv2pi2 = 1. / (2. * _PI_ * _PI_);
  Ec                   = 0.;
  dEc                  = 0.;
  const auto species   = [&](const Grid& g, const double* f, double am, double weight) {
    const double am2 = am * am;
    for (size_t j = 0; j < g.q.size(); ++j) {
      const double w    = weight * inv2pi2 * g.dq[j] * g.q[j] * g.q[j] * f[j];
      const double eps  = std::sqrt(g.q[j] * g.q[j] + am2);
      Ec               += w * eps;
      dEc              += w * am2 / eps;
    }
  };
  const double am = a * cfg_.m_H_mev;
  species(P_, y, am, kGH);
  species(D_, y + nP_, am * cfg_.r_l, kGL);
  species(D_, y + nP_ + nD_, am * cfg_.r_phi, kGPhi);
  rho = Ec / (a * a * a * a);
}

double DecayHarness::Hubble(double a, const double* y, double* rho_sector_out) const {
  double rho_sector = 0.;
  double Ec         = 0.;
  double dEc        = 0.;
  SectorEnergy(a, y, rho_sector, Ec, dEc);
  const BbnEmState em = plasma_.ElectromagneticState(y[iT_]);
  if (rho_sector_out) {
    *rho_sector_out = rho_sector;
  }
  return BbnPlasma::HubbleFromRho(em.rho + BbnPlasma::RhoNeutrino(1. / a, cfg_.n_eff_spectator) +
                                  rho_sector);
}

void DecayHarness::Derivs(double lna, const double* y, double* dy) const {
  const double a = std::exp(lna);
  double rho_sector;
  double Ec;
  double dEc;
  SectorEnergy(a, y, rho_sector, Ec, dEc);
  const double T      = y[iT_];
  const BbnEmState em = plasma_.ElectromagneticState(T);
  const double H      = BbnPlasma::HubbleFromRho(
      em.rho + BbnPlasma::RhoNeutrino(1. / a, cfg_.n_eff_spectator) + rho_sector);

  // The kernel returns a df/dt (it works in conformal time) with K = a^2 m Gamma, so
  // with a in MeV^-1, m in MeV and Gamma in s^-1 the result is in MeV^-1 s^-1.
  kernel_->ComputeBackgroundDerivs(a,
                                   cfg_.m_H_mev,
                                   gamma_s_,
                                   y,
                                   y + nP_,
                                   y + nP_ + nD_,
                                   dy,
                                   dy + nP_,
                                   dy + nP_ + nD_);
  const double to_lna = 1. / (a * H);
  for (int i = 0; i < iT_; ++i) {
    dy[i] *= to_lna;
  }
  // The plasma's own entropy: the decay sector never heats it.
  dy[iT_] = -3. * (em.rho + em.pressure) / em.drho_dT;
  dy[it_] = 1. / H;
  dy[iW_] = dEc;
}

void DecayHarness::Run(const std::vector<double>& T_nu_samples,
                       const std::function<void(const Snapshot&)>& on_sample) {
  const double a_ini  = 1. / cfg_.T_start_mev;
  const double a_end  = 1. / cfg_.T_end_mev;
  const double lna_in = std::log(a_ini);
  const double lna_en = std::log(a_end);

  std::vector<double> y(ny_, 0.);
  const double am = a_ini * cfg_.m_H_mev;
  for (int i = 0; i < nP_; ++i) {
    y[i] = 1. / (std::exp(std::sqrt(P_.q[i] * P_.q[i] + am * am)) + 1.);
  }
  if (cfg_.daughter_thermal) {
    const double aml = am * cfg_.r_l;
    for (int j = 0; j < nD_; ++j) {
      y[nP_ + j] = 1. / (std::exp(std::sqrt(D_.q[j] * D_.q[j] + aml * aml)) + 1.);
    }
  }
  y[iT_] = cfg_.T_start_mev;

  const double I_l_ini = 2. * BareNumberParent(y.data()) + 2. * BareNumberDaughter(y.data() + nP_);
  const double I_phi_ini = 2. * BareNumberParent(y.data()) +
                           BareNumberDaughter(y.data() + nP_ + nD_);
  double rho_ini;
  double Ec_ini;
  double dEc_ini;
  SectorEnergy(a_ini, y.data(), rho_ini, Ec_ini, dEc_ini);

  std::vector<double> xs;
  for (double T : T_nu_samples) {
    class_test_severe(!(T <= cfg_.T_start_mev) || !(T >= cfg_.T_end_mev),
                      "sample T_nu = %g MeV lies outside the run [%g, %g]",
                      T,
                      cfg_.T_end_mev,
                      cfg_.T_start_mev);
    xs.push_back(std::log(1. / T));
  }
  xs.push_back(lna_en);
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  if (xs.front() > lna_in) {
    xs.insert(xs.begin(), lna_in);  // the evolvers want two points, the first at x_ini
  }
  const bool first_is_start =
      (T_nu_samples.end() != std::find(T_nu_samples.begin(), T_nu_samples.end(), cfg_.T_start_mev));

  const std::vector<int> used(ny_, 1);
  const auto derivs = [this](double x, double* yy, double* dyy) { Derivs(x, yy, dyy); };
  const auto output = [&](double x, double* yy, double*, int) {
    if (x == lna_in && !first_is_start) {
      return;
    }
    const double a = std::exp(x);
    Snapshot s;
    s.T_nu    = 1. / a;
    s.T_gamma = yy[iT_];
    s.t       = yy[it_];
    s.H       = Hubble(a, yy);
    s.f_H     = yy;
    s.f_l     = yy + nP_;
    s.f_phi   = yy + nP_ + nD_;
    on_sample(s);
  };
  EvolverOptions opt;
  opt.rtol            = cfg_.rtol;
  opt.x_sampling      = xs.data();
  opt.x_sampling_size = static_cast<int>(xs.size());
  opt.used_in_output  = used.data();
  opt.output          = output;
  diag_.stats         = EvolverStats{};
  opt.stats           = &diag_.stats;
  evolver_ndf15(derivs, lna_in, lna_en, y.data(), ny_, opt);

  const double I_l_end = 2. * BareNumberParent(y.data()) + 2. * BareNumberDaughter(y.data() + nP_);
  const double I_phi_end = 2. * BareNumberParent(y.data()) +
                           BareNumberDaughter(y.data() + nP_ + nD_);
  double rho_end;
  double Ec_end;
  double dEc_end;
  SectorEnergy(a_end, y.data(), rho_end, Ec_end, dEc_end);
  diag_.number_defect_l   = (I_l_end - I_l_ini) / I_l_ini;
  diag_.number_defect_phi = (I_phi_end - I_phi_ini) / I_phi_ini;
  diag_.energy_defect     = (Ec_end - Ec_ini - y[iW_]) / Ec_ini;
}
