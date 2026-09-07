#include "bbn_solver.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "bbn_plasma.h"
#include "bbn_rates.h"
#include "bbn_weak.h"
#include "constants.h"
#include "errors.h"
#include "evolver_ndf15.h"

namespace {

/* State vector layout. The abundances start at kIdxY and run in kNuclides order,
   which is why that enum is documented as part of the interface. */
constexpr int kIdxLnA  = 0;
constexpr int kIdxEta  = 1;
constexpr int kIdxTime = 2;
constexpr int kIdxY    = 3;
constexpr int kNeq     = kIdxY + kBbnNumNuclides;

/** Everything the right-hand side needs, plus the history buffer. */
struct BbnWorkspace {
  const BbnInput* input;
  const BbnRates* rates;
  const BbnWeakRates* weak;
  const BbnWeakTable* weak_table; /**< null -> evaluate the integrals directly */
  const BbnPlasma* plasma;

  double T9_initial;
  double T_nu_initial_mev; /**< neutrino temperature at a = 1 */
  double neff_total;       /**< 3 Standard-Model species plus delta_neff */
  bool freeze_network;     /**< normalization pre-pass: expansion only, no reactions */

  std::vector<double> history;        /**< flattened rows of (T9, t, and the abundances) */
  std::vector<double> trajectory_lna; /**< ln a at each sampled tau, from the pre-pass */
};

/** Photon temperature in MeV at integration coordinate tau = ln(T9_ini / T9). */
inline double TemperatureMeV(double tau, double T9_initial) {
  return MeVFromT9(T9_initial * std::exp(-tau));
}

/** The right-hand side, in the independent variable tau = ln(T9_ini/T9).
 *
 *  Everything physical is computed per unit TIME first and converted at the end
 *  by dividing by dtau/dt. That ordering matters: the temperature equation needs
 *  the nuclear rates, and the nuclear rates need only the temperature, so
 *  computing in time removes what would otherwise look like a circular
 *  dependency between dT/dtau and dY/dtau. */
void BbnDerivs(double tau, double* y, double* dy, void* parameters) {
  BbnWorkspace& ws = *static_cast<BbnWorkspace*>(parameters);

  const double T9  = ws.T9_initial * std::exp(-tau);
  const double T   = MeVFromT9(T9);
  const double a   = std::exp(y[kIdxLnA]);
  const double eta = y[kIdxEta];
  const double* Y  = &y[kIdxY];

  /* Neutrinos are instantaneously decoupled: their temperature is simply 1/a,
     never re-thermalized. Above ~2 MeV they are still coupled, but there the
     photons also go as 1/a because e+- annihilation has not begun, so the two
     descriptions agree and starting at T_nu = T_gamma is consistent. */
  const double T_nu = ws.T_nu_initial_mev / a;

  /* --- densities ------------------------------------------------------------ */
  const double n_gamma_cm3 = BbnPlasma::NumberDensityPhoton(T);
  const double n_b_cm3     = eta * n_gamma_cm3;
  const double n_b_mev3    = n_b_cm3 * BbnPlasma::InvCm3ToMeV3();

  double sum_Y  = 0.0;
  double sum_Ym = 0.0;
  for (int i = 0; i < kBbnNumNuclides; ++i) {
    sum_Y  += Y[i];
    sum_Ym += Y[i] * kNuclides[i].atomic_mass_mev;
  }

  const BbnEmState em    = ws.plasma->ElectromagneticState(T);
  const double rho_b     = n_b_mev3 * (sum_Ym + 1.5 * T * sum_Y);
  const double p_b       = n_b_mev3 * T * sum_Y;
  const double rho_nu    = BbnPlasma::RhoNeutrino(T_nu, ws.neff_total);
  const double rho_extra = ws.input->rho_extra ? ws.input->rho_extra(a) : 0.0;

  const double H = BbnPlasma::HubbleFromRho(em.rho + rho_b + rho_nu + rho_extra);

  /* --- abundance rates, per second ------------------------------------------ */
  double dYdt[kBbnNumNuclides] = {0.0};

  if (!ws.freeze_network) {
    double lambda_np = 0.0, lambda_pn = 0.0;
    if (ws.weak_table != nullptr) {
      ws.weak_table->Rates(tau, &lambda_np, &lambda_pn);
    }
    else {
      ws.weak->Rates(T, T_nu, &lambda_np, &lambda_pn);
    }
    const double weak_flow  = lambda_np * Y[kNucN] - lambda_pn * Y[kNucP];
    dYdt[kNucN]            -= weak_flow;
    dYdt[kNucP]            += weak_flow;

    for (int r = 0; r < kBbnNumReactions; ++r) {
      const BbnReaction& rxn = kReactions[r];

      /* Flows are reactions per baryon per second. The n_b power is one less
         than the number of particles on that side, because dividing the reaction
         density by n_b turns each Y back into a number density except one. */
      /* One fit evaluation per reaction: the reverse rate is that same forward
         rate times a detailed-balance ratio, so calling BbnRates::Reverse would
         re-evaluate the REACLIB exponential sum a second time for nothing. */
      const double lambda_forward = ws.rates->Forward(r, T9);

      double forward = ReactionPrefactor(rxn.reactants, rxn.n_reactants) * lambda_forward;
      for (int i = 0; i < rxn.n_reactants; ++i) {
        forward *= Y[rxn.reactants[i]];
      }
      for (int i = 1; i < rxn.n_reactants; ++i) {
        forward *= n_b_cm3;
      }

      double reverse = ReactionPrefactor(rxn.products, rxn.n_products) * lambda_forward *
                       ws.rates->ReverseOverForward(r, T9);
      for (int i = 0; i < rxn.n_products; ++i) {
        reverse *= Y[rxn.products[i]];
      }
      for (int i = 1; i < rxn.n_products; ++i) {
        reverse *= n_b_cm3;
      }

      const double net = forward - reverse;
      for (int i = 0; i < rxn.n_reactants; ++i) {
        dYdt[rxn.reactants[i]] -= net;
      }
      for (int i = 0; i < rxn.n_products; ++i) {
        dYdt[rxn.products[i]] += net;
      }
    }
  }

  /* --- temperature, from energy conservation of the coupled sector ---------- */
  /* d(rho a^3)/dt = -p d(a^3)/dt over photons + pairs + baryons. Expanding rho as
     a function of (T, n_b, Y_i) and collecting the dT/dt terms gives
         dT/dt = [ -3H(rho+p) + 3H rho_b - n_b sum_i (m_i + 3T/2) dY_i/dt ] / C.
     The sum over dY_i/dt is the nuclear back-heating: because rho_b is built from
     real nuclear masses, binding-energy release needs no separate term. */
  const double heat_capacity = em.drho_dT + 1.5 * n_b_mev3 * sum_Y;

  double nuclear_source = 0.0;
  for (int i = 0; i < kBbnNumNuclides; ++i) {
    nuclear_source += (kNuclides[i].atomic_mass_mev + 1.5 * T) * dYdt[i];
  }

  const double dT_dt = (-3.0 * H * (em.rho + em.pressure + rho_b + p_b) + 3.0 * H * rho_b -
                        n_b_mev3 * nuclear_source) /
                       heat_capacity;

  /* tau = ln(T9_ini/T9) increases as the universe cools, so dtau/dt > 0. */
  const double dtau_dt = -dT_dt / T;
  class_test(!(dtau_dt > 0.0),
             "BBN temperature stopped decreasing at T9=%g (dtau/dt=%g); the "
             "integration variable is no longer monotonic",
             T9,
             dtau_dt);

  /* --- convert to d/dtau ---------------------------------------------------- */
  dy[kIdxLnA] = H / dtau_dt;
  /* eta = n_b/n_gamma with n_b ~ a^-3 and n_gamma ~ T^3. */
  dy[kIdxEta]  = eta * (-3.0 * H - 3.0 * dT_dt / T) / dtau_dt;
  dy[kIdxTime] = 1.0 / dtau_dt;
  for (int i = 0; i < kBbnNumNuclides; ++i) {
    dy[kIdxY + i] = dYdt[i] / dtau_dt;
  }
}

void BbnOutput(double tau, double* y, double* /*dy*/, int index, void* parameters) {
  BbnWorkspace& ws = *static_cast<BbnWorkspace*>(parameters);

  /* The pre-pass records ln a, which is what turns T_nu into a function of the
     integration coordinate and so makes the weak-rate table possible. */
  if (index < static_cast<int>(ws.trajectory_lna.size())) {
    ws.trajectory_lna[index] = y[kIdxLnA];
  }

  if (ws.history.empty()) {
    return;
  }
  const int stride = 2 + kBbnNumNuclides;
  double* row      = &ws.history[static_cast<std::size_t>(index) * stride];
  row[0]           = ws.T9_initial * std::exp(-tau);
  row[1]           = y[kIdxTime];
  for (int i = 0; i < kBbnNumNuclides; ++i) {
    row[2 + i] = y[kIdxY + i];
  }
}

/** Thermal number-density scale (m T / 2 pi hbar^2)^{3/2} in cm^-3, for the Saha
 *  initial conditions. Deliberately the same expression detailed balance uses. */
double ThermalScale(const BbnNuclide& nuc, double T_mev) {
  const double theta = NuclearMassMeV(nuc) * T_mev / (2.0 * M_PI);
  return std::pow(theta, 1.5) / BbnPlasma::InvCm3ToMeV3();
}

/** Nuclear statistical equilibrium abundances at the starting temperature.
 *
 *  n and p come from weak equilibrium using the ACTUAL rates rather than a
 *  Boltzmann factor, so the initial condition is consistent with the same
 *  integrand the solver will use. Everything heavier follows from Saha with the
 *  proton and neutron densities eliminated in favour of their abundances. */
void SahaInitialConditions(
    const BbnWeakRates& weak, double T_mev, double T_nu_mev, double n_b_cm3, double* Y) {
  double lambda_np = 0.0, lambda_pn = 0.0;
  weak.Rates(T_mev, T_nu_mev, &lambda_np, &lambda_pn);

  Y[kNucN] = lambda_pn / (lambda_np + lambda_pn);
  Y[kNucP] = 1.0 - Y[kNucN];

  const double n_n    = Y[kNucN] * n_b_cm3;
  const double n_p    = Y[kNucP] * n_b_cm3;
  const double base_n = n_n / (kNuclides[kNucN].g * ThermalScale(kNuclides[kNucN], T_mev));
  const double base_p = n_p / (kNuclides[kNucP].g * ThermalScale(kNuclides[kNucP], T_mev));

  for (int i = 0; i < kBbnNumNuclides; ++i) {
    if (i == kNucN || i == kNucP) {
      continue;
    }
    const BbnNuclide& nuc = kNuclides[i];
    const int Z           = nuc.Z;
    const int N           = nuc.A - nuc.Z;

    /* Binding energy. Atomic masses are correct here: the Z electrons appear on
       both sides and cancel. */
    const double binding = Z * kNuclides[kNucP].atomic_mass_mev +
                           N * kNuclides[kNucN].atomic_mass_mev - nuc.atomic_mass_mev;

    const double n_i = nuc.g * ThermalScale(nuc, T_mev) * std::pow(base_p, Z) *
                       std::pow(base_n, N) * std::exp(binding / T_mev);
    Y[i]             = n_i / n_b_cm3;
  }
}

}  // namespace

BbnResult SolveBbn(const BbnInput& input) {
  class_test_severe(input.rates_file.empty(), "BBN solve requires a reaction-rate file");
  class_test_severe(!(input.T9_initial > input.T9_final),
                    "BBN needs T9_initial (%g) above T9_final (%g)",
                    input.T9_initial,
                    input.T9_final);
  class_test_severe(!(input.tolerance > 0.),
                    "BBN needs a positive integration tolerance, got %g",
                    input.tolerance);
  /* Either no table at all, or enough nodes to spline: the sampling grid divides
     by (points - 1), so 1 would divide by zero and anything below 4 cannot carry
     a cubic. */
  class_test_severe(input.weak_table_points != 0 && input.weak_table_points < 4,
                    "BBN weak_table_points must be 0 (exact integrals) or at least 4, got %d",
                    input.weak_table_points);
  /* Numeric physical inputs, so these reject the point rather than aborting --
     a sampler can reach them, unlike the configuration checks above. */
  class_test(!(input.omega_b > 0.), "BBN needs a positive omega_b, got %g", input.omega_b);
  class_test(!(input.T_cmb > 0.), "BBN needs a positive T_cmb, got %g K", input.T_cmb);
  class_test(!(input.tau_n > 0.), "BBN needs a positive neutron lifetime, got %g s", input.tau_n);

  const BbnRates rates(input.rates_file);
  const BbnWeakRates weak(input.tau_n);
  const BbnPlasma plasma;

  const double T_initial_mev = MeVFromT9(input.T9_initial);
  const double tau_end       = std::log(input.T9_initial / input.T9_final);

  BbnWorkspace ws;
  ws.input            = &input;
  ws.rates            = &rates;
  ws.weak             = &weak;
  ws.plasma           = &plasma;
  ws.T9_initial       = input.T9_initial;
  ws.T_nu_initial_mev = T_initial_mev;
  /* Three Standard-Model species at the instantaneous-decoupling temperature,
     plus whatever delta_neff asks for as genuinely extra radiation. */
  ws.neff_total     = 3.0 + input.delta_neff;
  ws.freeze_network = false;
  ws.weak_table     = nullptr;

  std::vector<int> used_in_output(kNeq, 1);
  std::vector<double> y(kNeq, 0.0);

  /* omega_b fixes the baryon MASS density; eta is a NUMBER ratio. The conversion
     divides by the mean mass per baryon, sum_i Y_i m_i, which depends on the final
     composition and so depends weakly on the answer -- about 0.2% between a
     helium-free and a helium-rich mixture. Iterating removes that entirely
     instead of hiding it in a fixed omega_b-to-eta coefficient. */
  const double rho_crit_cgs_per_h2 = 1.87847e-29; /* g cm^-3 */
  const double T_cmb_mev           = input.T_cmb * _k_B_ / _eV_ * 1.0e-6;
  const double n_gamma_today       = BbnPlasma::NumberDensityPhoton(T_cmb_mev);
  const double mev_to_gram         = 1.0e6 * _eV_ / (_c_ * _c_) * 1.0e3;

  BbnResult result;
  result.abundances.assign(kBbnNumNuclides, 0.0);

  /* First guess at the mean mass per baryon, for a hydrogen-plus-helium mixture
     with helium mass fraction kGuessYp:
         <m> = 1 / [ (1-Y)/m_H + 4Y/m_He ].
     Starting from pure hydrogen instead is 0.18% out and costs an extra pass; the
     guess below is within ~3e-5 of every composition this solver produces, so one
     correction pass suffices. It is only a starting point -- the loop measures the
     real value and the break criterion, not the guess, sets the accuracy. */
  constexpr double kGuessYp = 0.245;
  double mean_mass_mev      = 1.0 / ((1.0 - kGuessYp) / kNuclides[kNucP].atomic_mass_mev +
                                     4.0 * kGuessYp / kNuclides[kNucHe4].atomic_mass_mev);
  double eta_target         = 0.0;

  /* ndf15 dereferences the sampling vector unconditionally, so every call gets a
     real one; where nothing is being recorded it is just the two endpoints, and
     BbnOutput returns immediately because the buffers are empty.
 
     The pre-pass samples on the weak-rate table's grid instead, because the
     trajectory it records is exactly what the table needs. */
  const bool tabulate_weak  = input.weak_table_points > 0;
  const int prepass_samples = tabulate_weak ? input.weak_table_points : 2;
  std::vector<double> prepass_tau(prepass_samples);
  for (int i = 0; i < prepass_samples; ++i) {
    prepass_tau[i] = tau_end * i / (prepass_samples - 1.0);
  }
  if (tabulate_weak) {
    ws.trajectory_lna.assign(prepass_samples, 0.0);
  }

  /* --- eta normalization pre-pass, run ONCE ---------------------------------
     The boundary condition on eta is at the END of the integration -- the value
     implied by omega_b today -- while the solve needs it at the start. Baryons
     are ~1e-6 of the energy density, so the growth factor eta_final/eta_initial
     is a property of the plasma alone: one cheap frozen-network solve measures
     it. That makes the normalization exact whatever the expansion history does,
     rather than assuming the textbook 11/4 entropy factor.

     Being a property of the plasma is also why this sits OUTSIDE the composition
     loop below. It is insensitive to eta at the 1e-6 level, so re-running it for
     each refinement of the mean baryon mass bought nothing and doubled the cost.

     It does run at a PHYSICAL eta rather than an arbitrary one: the growth factor
     is eta-independent only while baryons are negligible, and eta = 1 would make
     them dominate the energy density outright.

     The frozen network also makes this pass cheap in its own right -- with no
     reactions the weak rates are never evaluated, and those are ~87% of a full
     right-hand side. */
  eta_target = rho_crit_cgs_per_h2 * input.omega_b / (mean_mass_mev * mev_to_gram) / n_gamma_today;

  ws.freeze_network = true;
  ws.history.clear();
  y.assign(kNeq, 0.0);
  y[kIdxEta] = eta_target;
  SahaInitialConditions(weak,
                        T_initial_mev,
                        T_initial_mev,
                        eta_target * BbnPlasma::NumberDensityPhoton(T_initial_mev),
                        &y[kIdxY]);

  EvolverOptions options;
  options.rtol            = input.tolerance;
  options.used_in_output  = used_in_output.data();
  options.output          = BbnOutput;
  options.x_sampling      = prepass_tau.data();
  options.x_sampling_size = prepass_samples;

  evolver_ndf15(BbnDerivs, 0.0, tau_end, y.data(), kNeq, &ws, options);

  const double eta_growth = y[kIdxEta] / eta_target;
  result.aT_ratio         = std::exp(y[kIdxLnA]) * input.T9_final / input.T9_initial;

  /* --- the splined weak rates ----------------------------------------------
     Built from the pre-pass trajectory, which supplies the T_nu that the rates
     need and that the network solve would otherwise have to be trusted to
     reproduce. The residual approximation is the difference between this
     trajectory and the real one -- the ~1e-6 baryon contribution to H -- and
     bbn_solver_test measures what it costs rather than assuming it is nothing. */
  std::unique_ptr<BbnWeakTable> weak_table;
  if (tabulate_weak) {
    std::vector<double> T_gamma(prepass_samples), T_nu(prepass_samples);
    for (int i = 0; i < prepass_samples; ++i) {
      T_gamma[i] = TemperatureMeV(prepass_tau[i], input.T9_initial);
      T_nu[i]    = T_initial_mev * std::exp(-ws.trajectory_lna[i]);
    }
    weak_table    = std::make_unique<BbnWeakTable>(weak, prepass_tau, T_gamma, T_nu);
    ws.weak_table = weak_table.get();
  }
  ws.trajectory_lna.clear();

  /* --- the coupled solve, refining the mean baryon mass --------------------- */
  for (int iteration = 0; iteration < 8; ++iteration) {
    eta_target = rho_crit_cgs_per_h2 * input.omega_b / (mean_mass_mev * mev_to_gram) /
                 n_gamma_today;

    ws.freeze_network = false;
    ws.history.assign(input.history_file.empty()
                          ? 0
                          : static_cast<std::size_t>(input.history_rows) * (2 + kBbnNumNuclides),
                      0.0);

    const int n_samples = ws.history.empty() ? 2 : input.history_rows;
    std::vector<double> tau_samples(n_samples);
    for (int i = 0; i < n_samples; ++i) {
      tau_samples[i] = tau_end * i / (n_samples - 1.0);
    }

    const double eta_initial = eta_target / eta_growth;
    y.assign(kNeq, 0.0);
    y[kIdxEta] = eta_initial;
    SahaInitialConditions(weak,
                          T_initial_mev,
                          T_initial_mev,
                          eta_initial * BbnPlasma::NumberDensityPhoton(T_initial_mev),
                          &y[kIdxY]);

    options.x_sampling      = tau_samples.data();
    options.x_sampling_size = static_cast<int>(tau_samples.size());

    evolver_ndf15(BbnDerivs, 0.0, tau_end, y.data(), kNeq, &ws, options);

    /* The pre-pass promised this; if it is not delivered, the assumption that the
       growth factor is eta-independent has failed and the answer is not the
       requested omega_b. Fail rather than report an abundance for the wrong
       baryon density. */
    class_test(std::fabs(y[kIdxEta] / eta_target - 1.0) > 1.0e-6,
               "BBN eta normalization missed its target: wanted %.8e, got %.8e",
               eta_target,
               y[kIdxEta]);

    double new_mean = 0.0;
    for (int i = 0; i < kBbnNumNuclides; ++i) {
      new_mean += y[kIdxY + i] * kNuclides[i].atomic_mass_mev;
    }
    const double shift = std::fabs(new_mean - mean_mass_mev) / mean_mass_mev;
    mean_mass_mev      = new_mean;
    /* eta is inversely proportional to the mean mass, and D/H -- the most
       eta-sensitive output -- responds as dln(D/H)/dln(eta) ~ -1.6. Converging
       this to 1e-8 therefore pins eta far below every other error in the problem.
       The previous 1e-12 bought nothing and cost two further passes. */
    if (shift < 1.0e-8) {
      break;
    }
  }

  result.eta_final = eta_target;

  /* --- report -------------------------------------------------------------- */
  for (int i = 0; i < kBbnNumNuclides; ++i) {
    result.abundances[i] = y[kIdxY + i];
  }

  const double Y_p   = y[kIdxY + kNucP];
  const double Y_he4 = y[kIdxY + kNucHe4];

  result.Yp_nucleon = 4.0 * Y_he4;
  result.Yp_mass    = Y_he4 * kNuclides[kNucHe4].atomic_mass_mev / mean_mass_mev;
  result.DoverH     = y[kIdxY + kNucD] / Y_p;
  /* Tritium beta-decays to He3, and Be7 electron-captures to Li7, both long after
     BBN freezes out; the observable quantities are the sums. */
  result.He3overH      = (y[kIdxY + kNucHe3] + y[kIdxY + kNucT]) / Y_p;
  result.Li7overH      = (y[kIdxY + kNucLi7] + y[kIdxY + kNucBe7]) / Y_p;
  result.Li6overH      = y[kIdxY + kNucLi6] / Y_p;
  result.eta10         = result.eta_final * 1.0e10;
  result.neutrino_neff = ws.neff_total;
  result.time_final    = y[kIdxTime];

  if (!input.history_file.empty()) {
    FILE* file = fopen(input.history_file.c_str(), "w");
    class_test(file == nullptr,
               "could not open BBN history file '%s' for writing",
               input.history_file.c_str());
    fprintf(file, "# BBN abundance history. Y_i = n_i / n_b.\n");
    fprintf(file,
            "#\n"
            "# Values far below the evolver's absolute error floor (1e-15) are noise and\n"
            "# can come out slightly NEGATIVE -- the solver is not resolving them, and does\n"
            "# not need to. Every abundance is comfortably positive by the time it matters:\n"
            "# the smallest one reported as a result, Li7+Be7, ends near 5e-10.\n");
    fprintf(file, "# %13s %15s", "T9", "t[s]");
    for (int i = 0; i < kBbnNumNuclides; ++i) {
      fprintf(file, " %15s", kNuclides[i].name);
    }
    fprintf(file, "\n");
    const int stride = 2 + kBbnNumNuclides;
    for (int r = 0; r < input.history_rows; ++r) {
      for (int c = 0; c < stride; ++c) {
        fprintf(file, " %15.8e", ws.history[static_cast<std::size_t>(r) * stride + c]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }

  return result;
}
