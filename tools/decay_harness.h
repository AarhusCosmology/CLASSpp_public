#pragma once

/** @file decay_harness.h One decay channel H <-> l + phi in the MeV-era universe, driven
 *  outside CLASS.
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
 *  instantaneous-decoupling standard model; decay_harness_test pins that limit.
 *
 *  WHAT IS INTEGRATED. y = [f_H | f_l | f_phi | T_gamma | t | W], against ln a.
 *  W = int (d eps/d ln a) weighted by the sector's occupations, i.e. the comoving
 *  energy the expansion hands the sector's massive particles. The collision term
 *  conserves comoving energy exactly, so E_c(a) - E_c(a_ini) - W is a check on the
 *  whole integration, not only on the kernel (Diagnostics::energy_defect).
 *
 *  STATISTICAL WEIGHTS. g_H = g_l = 2, g_phi = 1, the kernel's convention (the boson's
 *  dof ratio is folded into its deposit), and Hannestad's. */

#include <functional>
#include <memory>
#include <vector>

#include "bbn_plasma.h"
#include "decay_transition_kernel.h"
#include "evolver_options.h"

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

    /** Thermal neutrino species outside the sector, each nu + nubar at T_nu. Two is
     *  nu_e plus the one flavour that is neither parent nor daughter. */
    double n_eff_spectator = 2.;

    /** Initial occupations, at T_start. The parent is kinetically thermal there,
     *  f = 1/(exp(eps/T) + 1) with its mass in eps; the fermion daughter is a thermal
     *  light neutrino if `daughter_thermal`, else empty; phi starts empty. */
    bool daughter_thermal = true;

    /** Log grids in q = p/T_nu. 0 for the daughters' q_max means "size it": 5% above
     *  the largest momentum the network can emit by T_end (MaxDaughterMomentum), since
     *  undersizing is an energy leak, not a resolution trade-off.
     *
     *  The bottom has no such bound once a daughter is massive: one can be emitted at
     *  rest, and everything emitted below daughter_q_min is filed into the first bin
     *  (number-conserving; the energy shows in Diagnostics::energy_defect). So the
     *  daughters' first cell is given the whole phase-space volume [0, (q0+q1)/2]
     *  rather than the trapezoid's half cell: it IS the reservoir for q < q0, and its f
     *  is then the mean occupation there. With the half cell it held ~6x too little
     *  volume for a 96-point grid, so a massive daughter emitted near rest drove it to
     *  the Pauli limit and stalled the integrator (measured at m_F/m_H = 0.1, tau =
     *  0.2 s, from T = 10 MeV). The grid should still start a decade or two below the
     *  momenta of interest -- the dr_q_min convergence question of issue #389. */
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
    return P_.q;
  }
  const std::vector<double>& daughter_q() const {
    return D_.q;
  }
  const Diagnostics& diagnostics() const {
    return diag_;
  }
  const Config& config() const {
    return cfg_;
  }

  /** Bare comoving number Sum q^2 dq f on the parent / daughter grid, per dof and
   *  without the 1/(2 pi^2): the moment the kernel's identities are stated in. */
  double BareNumberParent(const double* f) const;
  double BareNumberDaughter(const double* f) const;

 private:
  struct Grid {
    std::vector<double> q, dq;
  };
  /** `floor_to_zero`: the first cell owns the whole volume [0, (q0+q1)/2], i.e.
   *  q0^2 dq0 = ((q0+q1)/2)^3 / 3, instead of the trapezoid's half cell. */
  static Grid LogGrid(double qmin, double qmax, int n, bool floor_to_zero = false);
  static GridView View(const Grid& g) {
    return GridView{g.q.data(), g.dq.data(), static_cast<int>(g.q.size())};
  }

  /** Physical energy density of the sector, MeV^4, and the comoving energy a^4 rho
   *  (dimensionless) with its expansion source d(a^4 rho)/dln a. */
  void SectorEnergy(double a, const double* y, double& rho, double& Ec, double& dEc) const;
  double Hubble(double a, const double* y, double* rho_sector_out = nullptr) const;
  void Derivs(double lna, const double* y, double* dy) const;

  Config cfg_;
  Grid P_, D_;
  BbnPlasma plasma_;
  std::unique_ptr<DecayTransitionKernel> kernel_;
  int nP_         = 0;
  int nD_         = 0;
  int iT_         = 0;  // T_gamma
  int it_         = 0;  // t
  int iW_         = 0;  // W
  int ny_         = 0;
  double gamma_s_ = 0.;  // 1/tau
  Diagnostics diag_;
};
