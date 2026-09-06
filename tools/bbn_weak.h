#pragma once

/** @file bbn_weak.h Two-temperature n <-> p weak rates, normalized to the neutron lifetime.
 *
 *  Design: docs/superpowers/specs/2026-09-06-bbn-nuclear-network-design.md section 5
 *
 *  These rates, not the nuclear network, are what set Y_p: essentially every
 *  neutron that survives to the deuterium bottleneck ends up in He4, so Y_p is
 *  fixed by the n/p ratio at weak freeze-out. The reaction rates in bbn_rates.h
 *  matter far more for D/H than for Y_p.
 *
 *  Three channels run in each direction, and they are NOT at a common
 *  temperature: the electrons and positrons follow T_gamma while the neutrinos,
 *  decoupled, follow T_nu. Conflating the two is a classic and well-hidden bug,
 *  so both are explicit arguments everywhere in this interface.
 *
 *  What is included: Born phase space, Pauli blocking in every final state, the
 *  Coulomb (Fermi function) correction on the channels where the charged lepton
 *  and the proton coexist, and normalization to the measured neutron lifetime,
 *  which absorbs G_F, V_ud and g_A.
 *
 *  What is NOT included: finite-nucleon-mass (recoil) corrections, thermal
 *  radiative corrections, and any neutrino chemical potential.
 */

#include <vector>

class BbnWeakTable;

class BbnWeakRates {
 public:
  /** @param tau_n_seconds  free neutron lifetime; the single normalization input
   *  @param nodes          Gauss-Legendre nodes per channel integral
   *
   *  These integrals converge algebraically rather than spectrally, because
   *  beta = sqrt(eps^2 - 1) puts a square-root branch point exactly at the
   *  eps = 1 endpoint of four of the six channels. Measured worst-case relative
   *  error against a 2048-node reference, over T9 in [0.1, 300]: 1.5e-5 at 64
   *  nodes, 1.6e-6 at 128, 5.3e-8 at the default 256, 1.4e-8 at 512. The default
   *  propagates to well under 1e-6 in Y_p, far below anything else in the
   *  budget, and the integrals are cheap enough that buying the extra margin
   *  over 128 costs nothing measurable. */
  explicit BbnWeakRates(double tau_n_seconds, int nodes = 256);

  /** Total n -> p and p -> n rates in s^-1, summed over all three channels.
   *
   *  @param T_gamma_mev  photon (= electron/positron) temperature
   *  @param T_nu_mev     neutrino temperature */
  void Rates(double T_gamma_mev, double T_nu_mev, double* lambda_np, double* lambda_pn) const;

  /** The vacuum decay rate this object reproduces, s^-1. Equals 1/tau_n by
   *  construction; exposed so a test can confirm the normalization rather than
   *  trust it. */
  double FreeDecayRate() const {
    return 1.0 / tau_n_;
  }

 private:
  double tau_n_;
  double norm_;           /**< K, s^-1: 1 / (tau_n * vacuum phase-space integral) */
  std::vector<double> x_; /**< Gauss-Legendre abscissae on [0,1] */
  std::vector<double> w_; /**< Gauss-Legendre weights on [0,1] */

  /** Integrate f over [a,b] with the stored rule. */
  template <typename F>
  double Integrate(double a, double b, F&& f) const {
    double sum = 0.0;
    for (std::size_t i = 0; i < x_.size(); ++i) {
      sum += w_[i] * f(a + (b - a) * x_[i]);
    }
    return sum * (b - a);
  }
};

/** The same rates, splined along a known temperature trajectory.
 *
 *  These integrals are ~87% of a right-hand-side evaluation (measured), and the
 *  network solve calls them tens of thousands of times at temperatures it has
 *  already been told. Precomputing them is the single largest saving available.
 *
 *  The trajectory is the point. The rates depend on BOTH temperatures, and T_nu
 *  is a state variable, so they are not a function of T_gamma alone in general.
 *  Along the actual solution they are, because T_nu/T_gamma is fixed by entropy
 *  conservation in the plasma -- which the eta normalization pre-pass already
 *  computes, for free, before the network solve begins. The residual error is the
 *  difference between the pre-pass trajectory and the real one, which is the
 *  ~1e-6 baryon contribution to H; bbn_weak_test measures what that costs.
 *
 *  ln(lambda) is what gets splined, not lambda: lambda_pn falls through seventy
 *  decades over the integration, and no polynomial interpolates that. Values
 *  below kLogFloor are clamped, which puts a kink in the spline in the region
 *  where lambda_pn is ~1e-300 and multiplies an abundance of order unity -- i.e.
 *  where it has been irrelevant for forty decades already.
 */
class BbnWeakTable {
 public:
  /** @param exact    rates to sample
   *  @param tau      UNIFORM ascending grid of the solver's integration coordinate
   *  @param T_gamma  photon temperature at each tau, MeV
   *  @param T_nu     neutrino temperature at each tau, MeV */
  BbnWeakTable(const BbnWeakRates& exact,
               const std::vector<double>& tau,
               const std::vector<double>& T_gamma,
               const std::vector<double>& T_nu);

  /** Interpolated rates, s^-1. Outside the tabulated range the endpoints are
   *  held: the evolver can ask marginally past tau_end on its final step, and
   *  every rate is flat there anyway. */
  void Rates(double tau, double* lambda_np, double* lambda_pn) const;

 private:
  /** ln of the smallest representable rate worth carrying. exp() of this is
   *  ~1e-300, so it is zero for every purpose in the network. */
  static constexpr double kLogFloor = -690.0;

  double tau_first_ = 0.;
  double tau_step_  = 0.;
  int n_            = 0;
  std::vector<double> ln_rates_;  /**< interleaved [ln_np, ln_pn] per node */
  std::vector<double> ln_second_; /**< spline second derivatives, same layout */
};
