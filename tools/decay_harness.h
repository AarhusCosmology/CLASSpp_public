#pragma once

/** @file decay_harness.h Decay channels H <-> l + phi in the MeV-era universe, driven
 *  outside CLASS: one channel (DecayHarness) or Hannestad & Madsen's three-level
 *  cascade nu_tau -> nu_mu + phi, nu_mu -> nu_e + phi (DecayCascade). Both are front
 *  ends of one engine, DecayNetwork, which owns the plasma and the expansion.
 *
 *  WHY OUTSIDE CLASS (issue #389, step 2). Hannestad & Madsen 1997 (astro-ph/9702125)
 *  put MeV-mass parents with second-scale lifetimes in the middle of e+- annihilation.
 *  CLASS's background pins T_gamma = T_cmb/a with a fixed N_eff, and a thermal MeV
 *  neutrino would dominate its budget there, so a CLASS species would be testing the
 *  cosmology rather than the collision operator. This harness owns its own expansion
 *  instead: the photon + e+- plasma from BbnPlasma, spectator neutrinos at T_nu = 1/a,
 *  and the decay sector's own energy density. The thing under test is the kernel.
 *
 *  UNITS. Temperatures and masses in MeV, times in s. The scale factor is normalised
 *  as a = 1/T_nu in MeV^-1, so the kernel's comoving momentum q = a p is p/T_nu and
 *  a*m = m/T_nu. The kernel's right-hand side is a df/dt (conformal time), so the
 *  harness divides by a.
 *
 *  THE PLASMA. The decay sector is decoupled from the plasma, so the plasma's own
 *  entropy is conserved: dT_gamma/dln a = -3 (rho + p) / (drho/dT). T_gamma/T_nu
 *  therefore runs from 1 to (11/4)^{1/3} through annihilation, exactly as in the
 *  instantaneous-decoupling standard model; decay_harness_test pins that limit. The
 *  same fact makes T_gamma(a) a function of the plasma alone, which is what
 *  DecayNetwork::TnuAtTgamma inverts.
 *
 *  WHAT IS INTEGRATED. y = [f_species... | T_gamma | t | W], against ln a.
 *  W = int (d eps/d ln a) weighted by the sector's occupations, i.e. the comoving
 *  energy the expansion hands the sector's massive particles. The collision term
 *  conserves comoving energy exactly, so E_c(a) - E_c(a_ini) - W is a check on the
 *  whole integration, not only on the kernels.
 *
 *  STATISTICAL WEIGHTS. g = 2 for every neutrino, g_phi = 1: the kernel's convention
 *  (the boson's dof ratio is folded into its deposit), and Hannestad's. */

#include <functional>
#include <memory>
#include <vector>

#include "bbn_plasma.h"
#include "decay_transition_kernel.h"
#include "evolver_options.h"

/** Species on momentum grids, coupled by DecayTransitionKernel channels, in the
 *  MeV-era expansion. Each channel is parent <-> fermion + boson with its own
 *  kernel, parent mass and Gamma; a species may be the daughter of one channel and
 *  the parent of the next (the cascade's nu_mu), and two species may share a grid.
 *
 *  ComputeBackgroundDerivs OVERWRITES its outputs, so every channel writes into its
 *  own scratch and the network sums: df_s = sum over the channels that touch s. Each
 *  channel's identities are exact per right-hand side, so every linear combination of
 *  bare numbers that each channel conserves is conserved by the sum. */
class DecayNetwork {
 public:
  struct Grid {
    std::vector<double> q, dq;
    GridView View() const {
      return GridView{q.data(), dq.data(), static_cast<int>(q.size())};
    }
  };

  /** Log grid in q with trapezoid weights, half cells at the ends.
   *
   *  `floor_to_zero`: the first cell owns the whole volume [0, (q0+q1)/2], i.e.
   *  q0^2 dq0 = ((q0+q1)/2)^3 / 3, instead of the trapezoid's half cell. A massive
   *  daughter can be emitted at rest, and everything emitted below q_min is filed into
   *  the first bin (number-conserving; the energy shows in the energy defect), so that
   *  bin IS the reservoir for q < q0 and its f is the mean occupation there. With the
   *  half cell it held ~6x too little volume for a 96-point grid, so a massive daughter
   *  emitted near rest drove it to the Pauli limit and stalled the integrator
   *  (measured at m_F/m_H = 0.1, tau = 0.2 s, from T = 10 MeV). */
  static Grid LogGrid(double qmin, double qmax, int n, bool floor_to_zero = false);

  struct Species {
    double m_mev = 0.; /**< rest mass */
    double g     = 2.; /**< internal dof: 2 for a neutrino, 1 for phi */
    int grid     = 0;  /**< index into Setup::grids */
  };

  /** parent <-> fermion + boson. The kernel's mass ratios are overwritten from the
   *  species' masses, so the energies the kernel conserves are the ones the network
   *  books. */
  struct Channel {
    int parent   = 0;
    int fermion  = 1;
    int boson    = 2;
    double tau_s = 1.; /**< parent rest-frame lifetime */
    DecayTransitionKernel::Config kernel;
  };

  struct Setup {
    std::vector<Grid> grids;
    std::vector<Species> species;
    std::vector<Channel> channels;
    double n_eff_spectator = 0.;  /**< thermal nu + nubar pairs at T_nu, outside the sector */
    double T_start_mev     = 10.; /**< T_gamma = T_nu here */
    double T_end_mev       = 0.01;
    double rtol            = 1e-6;
  };

  /** The state at one sample. Species s's occupations are y + offset(s). */
  struct State {
    double T_nu     = 0.; /**< MeV, = 1/a */
    double T_gamma  = 0.; /**< MeV */
    double t        = 0.; /**< s since T_start */
    double H        = 0.; /**< s^-1 */
    double W        = 0.; /**< comoving energy handed to the sector by expansion so far */
    const double* y = nullptr;
  };

  /** Integration-wide record, from the final state against the initial one. */
  struct Diagnostics {
    std::vector<double> N_ini, N_end; /**< bare numbers per species */
    double Ec_ini = 0.;
    double Ec_end = 0.;
    double W_end  = 0.;
    /** (E_c(end) - E_c(ini) - W) / E_c(ini). Catches energy filed at a grid edge by an
     *  undersized grid (the kernel's clamped_energy_residual) and at a floor cell. */
    double energy_defect = 0.;
    EvolverStats stats; /**< ndf15's own counters, for cost and stiffness studies */
  };

  explicit DecayNetwork(Setup setup);
  DecayNetwork(const DecayNetwork&)            = delete;
  DecayNetwork& operator=(const DecayNetwork&) = delete;

  /** Integrate from T_start to T_end. `f_ini` holds one vector per species, sized to
   *  its grid. `on_sample` is called at each T_nu in `T_nu_samples` (inside
   *  [T_end, T_start]) and at T_end. */
  void Run(const std::vector<std::vector<double>>& f_ini,
           const std::vector<double>& T_nu_samples,
           const std::function<void(const State&)>& on_sample);

  int n_species() const {
    return static_cast<int>(setup_.species.size());
  }
  int offset(int s) const {
    return offset_[s];
  }
  const Grid& grid(int s) const {
    return setup_.grids[setup_.species[s].grid];
  }
  const Setup& setup() const {
    return setup_;
  }
  const Diagnostics& diagnostics() const {
    return diag_;
  }

  /** Bare comoving number Sum q^2 dq f, per dof and without the 1/(2 pi^2): the
   *  moment the kernel's identities are stated in. */
  static double BareNumber(const Grid& g, const double* f);
  double BareNumber(int s, const double* y) const {
    return BareNumber(grid(s), y + offset_[s]);
  }

  /** Comoving energy a^4 rho of the sector (dimensionless) and its expansion source
   *  d(a^4 rho)/dln a; rho in MeV^4. */
  void SectorEnergy(double a, const double* y, double& rho, double& Ec, double& dEc) const;

  /** The T_nu at which the plasma reaches T_gamma (<= T_start): the plasma's entropy
   *  a^3 s(T_gamma) is conserved, so T_nu = T_start (s(T_gamma)/s(T_start))^{1/3}
   *  exactly, independent of the sector. */
  double TnuAtTgamma(double T_gamma_mev) const;

  /** Equilibrium Fermi-Dirac in energy, 1/(exp(sqrt(q^2 + (a m)^2)) + 1) at T = T_nu:
   *  the initial occupation of every neutrino in issue #389. */
  static std::vector<double> ThermalFermion(const Grid& g, double a_m);

 private:
  double Hubble(double a, const double* y) const;
  void Derivs(double lna, const double* y, double* dy) const;

  Setup setup_;
  BbnPlasma plasma_;
  std::vector<std::unique_ptr<DecayTransitionKernel>> kernels_;
  std::vector<double> gamma_s_;          // 1/tau per channel
  std::vector<int> offset_;              // species -> first slot in y
  mutable std::vector<double> scratch_;  // per channel: df_parent | df_fermion | df_boson
  std::vector<int> scratch_offset_;      // channel -> first slot in scratch_
  int n_occ_ = 0;                        // occupation slots
  int iT_    = 0;                        // T_gamma
  int it_    = 0;                        // t
  int iW_    = 0;                        // W
  int ny_    = 0;
  Diagnostics diag_;
};

/** One channel H <-> l + phi (issue #389, step 2): DecayNetwork with three species,
 *  l and phi on one shared daughter grid. */
class DecayHarness {
 public:
  struct Config {
    double m_H_mev = 1.0; /**< parent mass */
    double tau_s   = 1.0; /**< parent rest-frame lifetime */
    /** Daughter masses over the parent's (DecayTransitionKernel::Config). */
    double r_l   = 0.;
    double r_phi = 0.;

    double T_start_mev = 10.;  /**< T_gamma = T_nu here: the start of the run */
    double T_end_mev   = 0.01; /**< T_nu at the end of the run */

    /** Thermal neutrino flavours outside the sector, each nu + nubar at T_nu. With three
     *  flavours and one channel, parent and daughter are two of them, so ONE is left
     *  (Hannestad & Madsen's fig 1: H = nu_tau, F = nu_mu, spectator nu_e). The cascade
     *  has none: all three flavours are in its sector. */
    double n_eff_spectator = 1.;

    /** Initial occupations, at T_start. The parent is kinetically thermal there,
     *  f = 1/(exp(eps/T) + 1) with its mass in eps; the fermion daughter is a thermal
     *  neutrino (mass in eps too) if `daughter_thermal`, else empty; phi starts empty. */
    bool daughter_thermal = true;

    /** Log grids in q = p/T_nu. 0 for the daughters' q_max means "size it": 5% above
     *  the largest momentum the network can emit by T_end (MaxDaughterMomentum), since
     *  undersizing is an energy leak, not a resolution trade-off.
     *
     *  The bottom has no such bound once a daughter is massive: one can be emitted at
     *  rest, so the daughters' first cell is DecayNetwork::LogGrid's floor_to_zero
     *  reservoir. The grid should still start a decade or two below the momenta of
     *  interest -- the dr_q_min convergence question of issue #389. */
    double parent_q_min   = 1e-2;
    double parent_q_max   = 30.;
    int parent_n          = 48;
    double daughter_q_min = 1e-3;
    double daughter_q_max = 0.;
    int daughter_n        = 96;

    /** The kernel's own flags. The default is the scheme the kernel documents as one
     *  piece: inverse decays and quantum statistics on, balanced gather with the exact
     *  split, two Gauss points per emission sub-interval. Mass ratios are overwritten
     *  from r_l / r_phi. */
    DecayTransitionKernel::Config kernel = DefaultKernel();

    double rtol = 1e-6;

    static DecayTransitionKernel::Config DefaultKernel() {
      DecayTransitionKernel::Config k;
      k.inverse_decays     = true;
      k.quantum_statistics = true;
      k.balanced_gather    = true;
      k.lumped_loss        = false;
      k.emission_gauss     = 2;
      // Smooth floor on the balanced gather's log-odds (kernel Config::gather_floor).
      // Without it ndf15 stalls at m_H = 2 MeV, tau = 0.2 s (issue #389, fig 1).
      k.gather_floor = 1e-18;
      return k;
    }
  };

  /** The state at one sample. Occupations are the kernel's, on the harness grids. */
  struct Snapshot {
    double T_nu         = 0.; /**< MeV, = 1/a */
    double T_gamma      = 0.; /**< MeV */
    double t            = 0.; /**< s since T_start */
    double H            = 0.; /**< s^-1 */
    const double* f_H   = nullptr;
    const double* f_l   = nullptr;
    const double* f_phi = nullptr;
  };

  /** Integration-wide checks, all from the final state against the initial one.
   *  Relative to the sector's initial comoving number / energy. */
  struct Diagnostics {
    double number_defect_l   = 0.; /**< (2N_H + 2N_l) drift: the fermion-leg invariant */
    double number_defect_phi = 0.; /**< (2N_H + N_phi) drift: the boson-leg invariant */
    /** E_c(end) - E_c(ini) - W. This also catches energy filed at a grid edge by an
     *  undersized daughter grid (the kernel's clamped_energy_residual), since that is
     *  energy the network did not place where it was emitted. */
    double energy_defect = 0.;
    /** The integrator's own counters (ndf15), for cost and stiffness studies. */
    EvolverStats stats;
  };

  explicit DecayHarness(Config cfg);

  /** Integrate from T_start to T_end. `on_sample` is called at each T_nu in
   *  `T_nu_samples` (descending, inside [T_end, T_start]) and at T_end. */
  void Run(const std::vector<double>& T_nu_samples,
           const std::function<void(const Snapshot&)>& on_sample);

  const std::vector<double>& parent_q() const {
    return net_->grid(0).q;
  }
  const std::vector<double>& daughter_q() const {
    return net_->grid(1).q;
  }
  const Diagnostics& diagnostics() const {
    return diag_;
  }
  const Config& config() const {
    return cfg_;
  }

  /** Bare comoving number Sum q^2 dq f on the parent / daughter grid, per dof and
   *  without the 1/(2 pi^2): the moment the kernel's identities are stated in. */
  double BareNumberParent(const double* f) const {
    return DecayNetwork::BareNumber(net_->grid(0), f);
  }
  double BareNumberDaughter(const double* f) const {
    return DecayNetwork::BareNumber(net_->grid(1), f);
  }

 private:
  Config cfg_;
  std::unique_ptr<DecayNetwork> net_;
  Diagnostics diag_;
};

/** Hannestad & Madsen's three-level cascade (issue #389, step 3):
 *
 *      K1: nu_tau <-> nu_mu + phi,     K2: nu_mu <-> nu_e + phi,
 *
 *  two DecayTransitionKernel instances on ONE nu_mu grid (K1's fermion daughter and
 *  K2's parent) and ONE phi grid (both kernels' boson). Every flavour is in the sector,
 *  so there are no spectators by default.
 *
 *  EXACT INVARIANTS (bare numbers N = Sum q^2 dq f per dof). K1 conserves 2N_tau +
 *  2N_mu and 2N_tau + N_phi; K2 conserves 2N_mu + 2N_e and 2N_mu + N_phi. So the sum
 *  conserves
 *
 *      2N_tau + 2N_mu + 2N_e           (fermion number)
 *      N_phi + 2N_tau - 2N_e           (phi leg)
 *      E_c - W                         (comoving energy, all four species)
 *
 *  and once both parents are gone N_phi = 4 N_tau,0 + 2 N_mu,0 + N_phi,0. */
class DecayCascade {
 public:
  enum Index : int { kTau = 0, kMu = 1, kE = 2, kPhi = 3, kNumSpecies = 4 };

  struct Config {
    double m_tau_mev = 2.;
    double m_mu_mev  = 0.5;
    double m_e_mev   = 0.005;
    double m_phi_mev = 0.005;
    double tau_tau_s = 0.2; /**< nu_tau rest-frame lifetime (K1) */
    double tau_mu_s  = 1.0; /**< nu_mu rest-frame lifetime (K2) */

    double T_start_mev = 5.;   /**< T_gamma = T_nu here (issue #389: 5.00 MeV) */
    double T_end_mev   = 0.01; /**< T_nu at the end of the run */

    /** No flavour is left outside the cascade. */
    double n_eff_spectator = 0.;

    /** Initial occupations at T_start: all three neutrinos equilibrium Fermi-Dirac in
     *  energy (each with its mass in eps), phi empty. */
    bool mu_thermal = true;
    bool e_thermal  = true;

    /** Log grids in q = p/T_nu. 0 for a q_max means "size it" as DecayHarness does:
     *  5% above the largest momentum the channels feeding the species can emit by
     *  T_end. nu_mu, nu_e and phi get the floor_to_zero reservoir cell, since each is a
     *  daughter that can be emitted at (or, for phi, near) rest. phi reaches far down:
     *  the lasing figures put q^2 f_phi ~ 1 at q = 1e-4. */
    double tau_q_min = 1e-2;
    double tau_q_max = 30.;
    int tau_n        = 48;
    double mu_q_min  = 1e-3;
    double mu_q_max  = 0.;
    int mu_n         = 112;
    double e_q_min   = 1e-4;
    double e_q_max   = 0.;
    int e_n          = 128;
    double phi_q_min = 1e-5;
    double phi_q_max = 0.;
    int phi_n        = 160;

    /** Both kernels' flags (mass ratios overwritten from the masses). */
    DecayTransitionKernel::Config kernel = DefaultKernel();

    /** DecayHarness's scheme with a wider gather floor, d = 1e-10 instead of 1e-18.
     *  phi starts empty in every bin, and every bin crosses the floor, where the
     *  gather's slope G/(2d) is largest; below some d ndf15 falls into a regime of
     *  thousands of Jacobians, and where that cliff sits depends on the grid. Measured on
     *  the lasing model of figs 5-7:
     *    - T = 5 -> 4 MeV, coarse grids: d <= 1e-15 needs ~2700 Jacobians (73 s), d >=
     *      1e-12 27 (1 s);
     *    - the whole run on test-decay-cascade's grid (16/36/36/48): d = 1e-12 still
     *      needs 1381 Jacobians (44 s), d = 1e-10 and 1e-9 188 (12 s);
     *    - T = 5 -> 2.3 MeV on grids 25/57/65/81: d = 1e-12 needs 1802 Jacobians.
     *  The physics does not see it: d = 1e-12 -> 1e-10 moves the final cold peak by 3e-4
     *  relative on grids 25/57/65/81 and 49/113/129/161, and d up to 1e-9 moves
     *  f_phi(q -> 0) by 0.5%. DecayHarness keeps 1e-18: at fig 1, 1e-12 raises its
     *  run-level number drift to 1e-10 for no speed. */
    static DecayTransitionKernel::Config DefaultKernel() {
      DecayTransitionKernel::Config k = DecayHarness::Config::DefaultKernel();
      k.gather_floor                  = 1e-10;
      return k;
    }

    double rtol = 1e-6;
  };

  struct Snapshot {
    double T_nu         = 0.;
    double T_gamma      = 0.;
    double t            = 0.;
    double H            = 0.;
    double W            = 0.; /**< DecayNetwork::State::W */
    const double* f_tau = nullptr;
    const double* f_mu  = nullptr;
    const double* f_e   = nullptr;
    const double* f_phi = nullptr;
  };

  /** All three invariants relative to the initial fermion number 2N_tau + 2N_mu +
   *  2N_e (the phi-leg combination starts near zero, so it gets the same scale), and
   *  the energy relative to the initial comoving energy. */
  struct Diagnostics {
    double fermion_number_defect = 0.;
    double phi_number_defect     = 0.;
    double energy_defect         = 0.;
    double N_ini[kNumSpecies]    = {};
    double N_end[kNumSpecies]    = {};
    EvolverStats stats;
  };

  explicit DecayCascade(Config cfg);

  void Run(const std::vector<double>& T_nu_samples,
           const std::function<void(const Snapshot&)>& on_sample);

  const std::vector<double>& q(int species) const {
    return net_->grid(species).q;
  }
  const DecayNetwork::Grid& grid(int species) const {
    return net_->grid(species);
  }
  double BareNumber(int species, const double* f) const {
    return DecayNetwork::BareNumber(net_->grid(species), f);
  }
  /** Bare comoving energy Sum q^2 dq eps f at a = 1/T_nu, per dof. */
  double BareEnergy(int species, const double* f, double T_nu_mev) const;
  double TnuAtTgamma(double T_gamma_mev) const {
    return net_->TnuAtTgamma(T_gamma_mev);
  }
  double mass(int species) const {
    return net_->setup().species[species].m_mev;
  }
  const Diagnostics& diagnostics() const {
    return diag_;
  }
  const Config& config() const {
    return cfg_;
  }
  const DecayNetwork& network() const {
    return *net_;
  }

 private:
  Config cfg_;
  std::unique_ptr<DecayNetwork> net_;
  Diagnostics diag_;
};
