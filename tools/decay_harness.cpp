#include "decay_harness.h"

#include <algorithm>
#include <cmath>

#include "common.h"
#include "errors.h"
#include "evolver_ndf15.h"
#include "evolver_options.h"

// ---------------------------------------------------------------------------------------
// DecayNetwork
// ---------------------------------------------------------------------------------------

DecayNetwork::Grid DecayNetwork::LogGrid(double qmin, double qmax, int n, bool floor_to_zero) {
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

double DecayNetwork::BareNumber(const Grid& g, const double* f) {
  double N = 0.;
  for (size_t j = 0; j < g.q.size(); ++j) {
    N += g.dq[j] * g.q[j] * g.q[j] * f[j];
  }
  return N;
}

std::vector<double> DecayNetwork::ThermalFermion(const Grid& g, double a_m) {
  std::vector<double> f(g.q.size());
  for (size_t j = 0; j < g.q.size(); ++j) {
    f[j] = 1. / (std::exp(std::sqrt(g.q[j] * g.q[j] + a_m * a_m)) + 1.);
  }
  return f;
}

DecayNetwork::DecayNetwork(Setup setup) : setup_(std::move(setup)) {
  class_test_severe(!(setup_.T_start_mev > setup_.T_end_mev) || !(setup_.T_end_mev > 0.),
                    "decay harness runs from T_start down to T_end > 0 (got %g -> %g MeV)",
                    setup_.T_start_mev,
                    setup_.T_end_mev);
  const int ns = n_species();
  class_test_severe(ns < 1, "decay network needs at least one species");
  offset_.resize(ns);
  n_occ_ = 0;
  for (int s = 0; s < ns; ++s) {
    const Species& sp = setup_.species[s];
    class_test_severe(sp.grid < 0 || sp.grid >= static_cast<int>(setup_.grids.size()),
                      "decay network species %d names grid %d of %d",
                      s,
                      sp.grid,
                      static_cast<int>(setup_.grids.size()));
    class_test_severe(!(sp.m_mev >= 0.) || !(sp.g > 0.),
                      "decay network species %d needs m >= 0 and g > 0 (got %g, %g)",
                      s,
                      sp.m_mev,
                      sp.g);
    offset_[s]  = n_occ_;
    n_occ_     += static_cast<int>(grid(s).q.size());
  }

  // The grids are not touched again, so the kernels' views into them stay valid.
  int scratch = 0;
  for (const Channel& ch : setup_.channels) {
    const auto valid = [ns](int s) { return s >= 0 && s < ns; };
    class_test_severe(!valid(ch.parent) || !valid(ch.fermion) || !valid(ch.boson) ||
                          ch.parent == ch.fermion || ch.parent == ch.boson ||
                          ch.fermion == ch.boson,
                      "decay network channel needs three distinct species (got %d, %d, %d)",
                      ch.parent,
                      ch.fermion,
                      ch.boson);
    // The kernel's rates are per internal state with the boson leg weighted by
    // g_H/g_phi = 2 (folded into its deposit), so they conserve the energy SectorEnergy
    // books only when g_parent = g_fermion = 2 g_boson. Any other choice would keep the
    // kernel's N/E identities and still break E_c - W.
    {
      const double gH   = setup_.species[ch.parent].g;
      const double gl   = setup_.species[ch.fermion].g;
      const double gphi = setup_.species[ch.boson].g;
      class_test_severe(std::fabs(gl - gH) > 1e-12 * gH || std::fabs(gH - 2. * gphi) > 1e-12 * gH,
                        "decay network channel (%d -> %d + %d) needs g_parent = g_fermion = "
                        "2 g_boson, the kernel's weights (got %g, %g, %g)",
                        ch.parent,
                        ch.fermion,
                        ch.boson,
                        gH,
                        gl,
                        gphi);
    }
    const double mH = setup_.species[ch.parent].m_mev;
    class_test_severe(!(mH > 0.) || !(ch.tau_s > 0.),
                      "decay harness needs a positive mass and lifetime (got m = %g MeV, "
                      "tau = %g s)",
                      mH,
                      ch.tau_s);
    DecayTransitionKernel::Config k = ch.kernel;
    k.fermion_mass_ratio            = setup_.species[ch.fermion].m_mev / mH;
    k.boson_mass_ratio              = setup_.species[ch.boson].m_mev / mH;
    kernels_.push_back(std::make_unique<DecayTransitionKernel>(grid(ch.parent).View(),
                                                               grid(ch.fermion).View(),
                                                               grid(ch.boson).View(),
                                                               Statistics::Fermion,
                                                               Statistics::Boson,
                                                               k));
    gamma_s_.push_back(1. / ch.tau_s);
    scratch_offset_.push_back(scratch);
    scratch += static_cast<int>(grid(ch.parent).q.size() + grid(ch.fermion).q.size() +
                                grid(ch.boson).q.size());
  }
  scratch_.assign(scratch, 0.);

  iT_ = n_occ_;
  it_ = iT_ + 1;
  iW_ = iT_ + 2;
  ny_ = iT_ + 3;
}

void DecayNetwork::SectorEnergy(
    double a, const double* y, double& rho, double& Ec, double& dEc) const {
  // E_c = a^4 rho = Sum_i g_i/(2 pi^2) Sum q^2 dq eps_i f_i, all dimensionless, and its
  // expansion source d eps/dln a = (a m)^2 / eps per particle.
  const double inv2pi2 = 1. / (2. * _PI_ * _PI_);
  Ec                   = 0.;
  dEc                  = 0.;
  for (int s = 0; s < n_species(); ++s) {
    const Grid& g    = grid(s);
    const double* f  = y + offset_[s];
    const double am  = a * setup_.species[s].m_mev;
    const double am2 = am * am;
    const double wg  = setup_.species[s].g * inv2pi2;
    for (size_t j = 0; j < g.q.size(); ++j) {
      const double w    = wg * g.dq[j] * g.q[j] * g.q[j] * f[j];
      const double eps  = std::sqrt(g.q[j] * g.q[j] + am2);
      Ec               += w * eps;
      dEc              += w * am2 / eps;
    }
  }
  rho = Ec / (a * a * a * a);
}

double DecayNetwork::Hubble(double a, const double* y) const {
  double rho_sector = 0.;
  double Ec         = 0.;
  double dEc        = 0.;
  SectorEnergy(a, y, rho_sector, Ec, dEc);
  const BbnEmState em = plasma_.ElectromagneticState(y[iT_]);
  return BbnPlasma::HubbleFromRho(em.rho + BbnPlasma::RhoNeutrino(1. / a, setup_.n_eff_spectator) +
                                  rho_sector);
}

double DecayNetwork::TnuAtTgamma(double T_gamma_mev) const {
  class_test_severe(!(T_gamma_mev > 0.) || !(T_gamma_mev <= setup_.T_start_mev),
                    "T_gamma = %g MeV is outside (0, T_start = %g MeV]",
                    T_gamma_mev,
                    setup_.T_start_mev);
  const BbnEmState e0 = plasma_.ElectromagneticState(setup_.T_start_mev);
  const BbnEmState e1 = plasma_.ElectromagneticState(T_gamma_mev);
  const double s0     = (e0.rho + e0.pressure) / setup_.T_start_mev;
  const double s1     = (e1.rho + e1.pressure) / T_gamma_mev;
  return setup_.T_start_mev * std::cbrt(s1 / s0);
}

void DecayNetwork::Derivs(double lna, const double* y, double* dy) const {
  const double a = std::exp(lna);
  double rho_sector;
  double Ec;
  double dEc;
  SectorEnergy(a, y, rho_sector, Ec, dEc);
  const double T      = y[iT_];
  const BbnEmState em = plasma_.ElectromagneticState(T);
  const double H      = BbnPlasma::HubbleFromRho(
      em.rho + BbnPlasma::RhoNeutrino(1. / a, setup_.n_eff_spectator) + rho_sector);

  // Each kernel returns a df/dt (it works in conformal time) with K = a^2 m Gamma, so
  // with a in MeV^-1, m in MeV and Gamma in s^-1 the result is in MeV^-1 s^-1. Every
  // kernel overwrites its outputs, hence the per-channel scratch and the sum.
  std::fill(dy, dy + n_occ_, 0.);
  for (size_t c = 0; c < kernels_.size(); ++c) {
    const Channel& ch = setup_.channels[c];
    const int nP      = static_cast<int>(grid(ch.parent).q.size());
    const int nL      = static_cast<int>(grid(ch.fermion).q.size());
    const int nB      = static_cast<int>(grid(ch.boson).q.size());
    double* dP        = scratch_.data() + scratch_offset_[c];
    double* dL        = dP + nP;
    double* dB        = dL + nL;
    kernels_[c]->ComputeBackgroundDerivs(a,
                                         setup_.species[ch.parent].m_mev,
                                         gamma_s_[c],
                                         y + offset_[ch.parent],
                                         y + offset_[ch.fermion],
                                         y + offset_[ch.boson],
                                         dP,
                                         dL,
                                         dB);
    const auto add = [dy](int off, const double* src, int n) {
      for (int j = 0; j < n; ++j) {
        dy[off + j] += src[j];
      }
    };
    add(offset_[ch.parent], dP, nP);
    add(offset_[ch.fermion], dL, nL);
    add(offset_[ch.boson], dB, nB);
  }
  const double to_lna = 1. / (a * H);
  for (int i = 0; i < n_occ_; ++i) {
    dy[i] *= to_lna;
  }
  // The plasma's own entropy: the decay sector never heats it.
  dy[iT_] = -3. * (em.rho + em.pressure) / em.drho_dT;
  dy[it_] = 1. / H;
  dy[iW_] = dEc;
}

void DecayNetwork::Run(const std::vector<std::vector<double>>& f_ini,
                       const std::vector<double>& T_nu_samples,
                       const std::function<void(const State&)>& on_sample) {
  const double a_ini  = 1. / setup_.T_start_mev;
  const double a_end  = 1. / setup_.T_end_mev;
  const double lna_in = std::log(a_ini);
  const double lna_en = std::log(a_end);
  const int ns        = n_species();

  class_test_severe(static_cast<int>(f_ini.size()) != ns,
                    "decay network: %d initial occupations for %d species",
                    static_cast<int>(f_ini.size()),
                    ns);
  std::vector<double> y(ny_, 0.);
  for (int s = 0; s < ns; ++s) {
    class_test_severe(f_ini[s].size() != grid(s).q.size(),
                      "decay network: species %d has %d initial occupations for %d bins",
                      s,
                      static_cast<int>(f_ini[s].size()),
                      static_cast<int>(grid(s).q.size()));
    std::copy(f_ini[s].begin(), f_ini[s].end(), y.begin() + offset_[s]);
  }
  y[iT_] = setup_.T_start_mev;

  diag_ = Diagnostics{};
  diag_.N_ini.resize(ns);
  diag_.N_end.resize(ns);
  for (int s = 0; s < ns; ++s) {
    diag_.N_ini[s] = BareNumber(s, y.data());
  }
  double rho;
  double dEc;
  SectorEnergy(a_ini, y.data(), rho, diag_.Ec_ini, dEc);

  std::vector<double> xs;
  for (double T : T_nu_samples) {
    class_test_severe(!(T <= setup_.T_start_mev) || !(T >= setup_.T_end_mev),
                      "sample T_nu = %g MeV lies outside the run [%g, %g]",
                      T,
                      setup_.T_end_mev,
                      setup_.T_start_mev);
    xs.push_back(std::log(1. / T));
  }
  xs.push_back(lna_en);
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  if (xs.front() > lna_in) {
    xs.insert(xs.begin(), lna_in);  // the evolvers want two points, the first at x_ini
  }
  const bool first_is_start = (T_nu_samples.end() != std::find(T_nu_samples.begin(),
                                                               T_nu_samples.end(),
                                                               setup_.T_start_mev));

  const std::vector<int> used(ny_, 1);
  const auto derivs = [this](double x, double* yy, double* dyy) { Derivs(x, yy, dyy); };
  const auto output = [&](double x, double* yy, double*, int) {
    if (x == lna_in && !first_is_start) {
      return;
    }
    const double a = std::exp(x);
    State s;
    s.T_nu    = 1. / a;
    s.T_gamma = yy[iT_];
    s.t       = yy[it_];
    s.H       = Hubble(a, yy);
    s.W       = yy[iW_];
    s.y       = yy;
    on_sample(s);
  };
  EvolverOptions opt;
  opt.rtol            = setup_.rtol;
  opt.x_sampling      = xs.data();
  opt.x_sampling_size = static_cast<int>(xs.size());
  opt.used_in_output  = used.data();
  opt.output          = output;
  opt.stats           = &diag_.stats;
  evolver_ndf15(derivs, lna_in, lna_en, y.data(), ny_, opt);

  for (int s = 0; s < ns; ++s) {
    diag_.N_end[s] = BareNumber(s, y.data());
  }
  SectorEnergy(a_end, y.data(), rho, diag_.Ec_end, dEc);
  diag_.W_end         = y[iW_];
  diag_.energy_defect = (diag_.Ec_end - diag_.Ec_ini - diag_.W_end) / diag_.Ec_ini;
}

// ---------------------------------------------------------------------------------------
// DecayHarness: one channel
// ---------------------------------------------------------------------------------------

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

  DecayNetwork::Setup s;
  s.grids.push_back(DecayNetwork::LogGrid(cfg_.parent_q_min, cfg_.parent_q_max, cfg_.parent_n));
  // Both daughters share one grid. The network's emission saturates once the parent is
  // gone, but MaxDaughterMomentum is a bound at any a, so it is evaluated at the latest
  // one; the log grid makes the extra decades cheap.
  double qd_max = cfg_.daughter_q_max;
  if (!(qd_max > 0.)) {
    qd_max = 1.05 * DecayTransitionKernel::MaxDaughterMomentum(s.grids[0].q.back(),
                                                               cfg_.m_H_mev / cfg_.T_end_mev,
                                                               cfg_.r_l,
                                                               cfg_.r_phi);
  }
  s.grids.push_back(DecayNetwork::LogGrid(cfg_.daughter_q_min, qd_max, cfg_.daughter_n, true));
  s.species = {{cfg_.m_H_mev, 2., 0},
               {cfg_.r_l * cfg_.m_H_mev, 2., 1},
               {cfg_.r_phi * cfg_.m_H_mev, 1., 1}};
  DecayNetwork::Channel ch;
  ch.parent  = 0;
  ch.fermion = 1;
  ch.boson   = 2;
  ch.tau_s   = cfg_.tau_s;
  ch.kernel  = cfg_.kernel;
  s.channels.push_back(ch);
  s.n_eff_spectator = cfg_.n_eff_spectator;
  s.T_start_mev     = cfg_.T_start_mev;
  s.T_end_mev       = cfg_.T_end_mev;
  s.rtol            = cfg_.rtol;
  net_              = std::make_unique<DecayNetwork>(std::move(s));
}

void DecayHarness::Run(const std::vector<double>& T_nu_samples,
                       const std::function<void(const Snapshot&)>& on_sample) {
  const double a_ini = 1. / cfg_.T_start_mev;
  const auto& D      = net_->grid(1);
  std::vector<std::vector<double>> f(3);
  f[0] = DecayNetwork::ThermalFermion(net_->grid(0), a_ini * cfg_.m_H_mev);
  f[1] = cfg_.daughter_thermal ? DecayNetwork::ThermalFermion(D, a_ini * cfg_.m_H_mev * cfg_.r_l)
                               : std::vector<double>(D.q.size(), 0.);
  f[2].assign(D.q.size(), 0.);

  net_->Run(f, T_nu_samples, [&](const DecayNetwork::State& st) {
    Snapshot s;
    s.T_nu    = st.T_nu;
    s.T_gamma = st.T_gamma;
    s.t       = st.t;
    s.H       = st.H;
    s.f_H     = st.y + net_->offset(0);
    s.f_l     = st.y + net_->offset(1);
    s.f_phi   = st.y + net_->offset(2);
    on_sample(s);
  });

  const auto& d           = net_->diagnostics();
  const double I_l_ini    = 2. * d.N_ini[0] + 2. * d.N_ini[1];
  const double I_phi_ini  = 2. * d.N_ini[0] + d.N_ini[2];
  diag_.number_defect_l   = (2. * d.N_end[0] + 2. * d.N_end[1] - I_l_ini) / I_l_ini;
  diag_.number_defect_phi = (2. * d.N_end[0] + d.N_end[2] - I_phi_ini) / I_phi_ini;
  diag_.energy_defect     = d.energy_defect;
  diag_.stats             = d.stats;
}

// ---------------------------------------------------------------------------------------
// DecayCascade: nu_tau -> nu_mu + phi, nu_mu -> nu_e + phi
// ---------------------------------------------------------------------------------------

DecayCascade::DecayCascade(Config cfg) : cfg_(cfg) {
  class_test_severe(!(cfg_.m_tau_mev > cfg_.m_mu_mev + cfg_.m_phi_mev) ||
                        !(cfg_.m_mu_mev > cfg_.m_e_mev + cfg_.m_phi_mev) || !(cfg_.m_e_mev >= 0.) ||
                        !(cfg_.m_phi_mev >= 0.),
                    "decay cascade needs m_tau > m_mu + m_phi and m_mu > m_e + m_phi (got "
                    "%g, %g, %g, %g MeV)",
                    cfg_.m_tau_mev,
                    cfg_.m_mu_mev,
                    cfg_.m_e_mev,
                    cfg_.m_phi_mev);
  const double a_end = 1. / cfg_.T_end_mev;
  const double r1_l  = cfg_.m_mu_mev / cfg_.m_tau_mev;
  const double r1_p  = cfg_.m_phi_mev / cfg_.m_tau_mev;
  const double r2_l  = cfg_.m_e_mev / cfg_.m_mu_mev;
  const double r2_p  = cfg_.m_phi_mev / cfg_.m_mu_mev;

  DecayNetwork::Setup s;
  s.grids.push_back(DecayNetwork::LogGrid(cfg_.tau_q_min, cfg_.tau_q_max, cfg_.tau_n));
  const double k1_max = 1.05 * DecayTransitionKernel::MaxDaughterMomentum(cfg_.tau_q_max,
                                                                          a_end * cfg_.m_tau_mev,
                                                                          r1_l,
                                                                          r1_p);
  // nu_mu holds a thermal population initially, so it reaches at least the parent's top.
  const double mu_max = cfg_.mu_q_max > 0. ? cfg_.mu_q_max : std::fmax(k1_max, cfg_.tau_q_max);
  s.grids.push_back(DecayNetwork::LogGrid(cfg_.mu_q_min, mu_max, cfg_.mu_n, true));
  const double k2_max =
      1.05 * DecayTransitionKernel::MaxDaughterMomentum(mu_max, a_end * cfg_.m_mu_mev, r2_l, r2_p);
  const double e_max = cfg_.e_q_max > 0. ? cfg_.e_q_max : std::fmax(k2_max, cfg_.tau_q_max);
  s.grids.push_back(DecayNetwork::LogGrid(cfg_.e_q_min, e_max, cfg_.e_n, true));
  const double phi_max = cfg_.phi_q_max > 0. ? cfg_.phi_q_max : std::fmax(k1_max, k2_max);
  s.grids.push_back(DecayNetwork::LogGrid(cfg_.phi_q_min, phi_max, cfg_.phi_n, true));

  s.species = {{cfg_.m_tau_mev, 2., kTau},
               {cfg_.m_mu_mev, 2., kMu},
               {cfg_.m_e_mev, 2., kE},
               {cfg_.m_phi_mev, 1., kPhi}};
  DecayNetwork::Channel k1;
  k1.parent  = kTau;
  k1.fermion = kMu;
  k1.boson   = kPhi;
  k1.tau_s   = cfg_.tau_tau_s;
  k1.kernel  = cfg_.kernel;
  DecayNetwork::Channel k2;
  k2.parent         = kMu;
  k2.fermion        = kE;
  k2.boson          = kPhi;
  k2.tau_s          = cfg_.tau_mu_s;
  k2.kernel         = cfg_.kernel;
  s.channels        = {k1, k2};
  s.n_eff_spectator = cfg_.n_eff_spectator;
  s.T_start_mev     = cfg_.T_start_mev;
  s.T_end_mev       = cfg_.T_end_mev;
  s.rtol            = cfg_.rtol;
  net_              = std::make_unique<DecayNetwork>(std::move(s));
}

double DecayCascade::BareEnergy(int species, const double* f, double T_nu_mev) const {
  const auto& g   = net_->grid(species);
  const double am = mass(species) / T_nu_mev;
  double E        = 0.;
  for (size_t j = 0; j < g.q.size(); ++j) {
    E += g.dq[j] * g.q[j] * g.q[j] * std::sqrt(g.q[j] * g.q[j] + am * am) * f[j];
  }
  return E;
}

void DecayCascade::Run(const std::vector<double>& T_nu_samples,
                       const std::function<void(const Snapshot&)>& on_sample) {
  const double a_ini = 1. / cfg_.T_start_mev;
  std::vector<std::vector<double>> f(kNumSpecies);
  f[kTau] = DecayNetwork::ThermalFermion(grid(kTau), a_ini * cfg_.m_tau_mev);
  f[kMu]  = cfg_.mu_thermal ? DecayNetwork::ThermalFermion(grid(kMu), a_ini * cfg_.m_mu_mev)
                            : std::vector<double>(grid(kMu).q.size(), 0.);
  f[kE]   = cfg_.e_thermal ? DecayNetwork::ThermalFermion(grid(kE), a_ini * cfg_.m_e_mev)
                           : std::vector<double>(grid(kE).q.size(), 0.);
  f[kPhi].assign(grid(kPhi).q.size(), 0.);

  net_->Run(f, T_nu_samples, [&](const DecayNetwork::State& st) {
    Snapshot s;
    s.T_nu    = st.T_nu;
    s.T_gamma = st.T_gamma;
    s.t       = st.t;
    s.H       = st.H;
    s.W       = st.W;
    s.f_tau   = st.y + net_->offset(kTau);
    s.f_mu    = st.y + net_->offset(kMu);
    s.f_e     = st.y + net_->offset(kE);
    s.f_phi   = st.y + net_->offset(kPhi);
    on_sample(s);
  });

  const auto& d = net_->diagnostics();
  for (int i = 0; i < kNumSpecies; ++i) {
    diag_.N_ini[i] = d.N_ini[i];
    diag_.N_end[i] = d.N_end[i];
  }
  const auto fermion = [](const double* N) { return 2. * N[kTau] + 2. * N[kMu] + 2. * N[kE]; };
  const auto phi_leg = [](const double* N) { return N[kPhi] + 2. * N[kTau] - 2. * N[kE]; };
  const double scale = fermion(diag_.N_ini);
  diag_.fermion_number_defect = (fermion(diag_.N_end) - scale) / scale;
  diag_.phi_number_defect     = (phi_leg(diag_.N_end) - phi_leg(diag_.N_ini)) / scale;
  diag_.energy_defect         = d.energy_defect;
  diag_.stats                 = d.stats;
}
