#pragma once
#include <functional>
#include <utility>
#include <vector>

/**
 * A scalar-field potential bundled with its first two derivatives. Each callable
 * takes (phi, params); the owning ScalarFieldSpecies passes its scf_parameters_
 * as `params`, so the shooting machinery (which tunes params[tuning_index])
 * keeps working unchanged. The potential is evaluated only on the background ODE
 * (never in the per-k perturbation loop), so std::function indirection is free.
 */
struct ScalarFieldPotential {
  std::function<double(double phi, const std::vector<double>& params)> V;
  std::function<double(double phi, const std::vector<double>& params)> dV;
  std::function<double(double phi, const std::vector<double>& params)> ddV;

  // Returns {xguess, dxdy} for params[tuning_index] given the target Omega0_scf
  // and H0. The shooting guess is potential-specific (e.g. the frozen-field 1EXP
  // V0 guess differs from the exponential-quintessence lambda guess), so it lives
  // with V/dV/ddV. ScalarFieldSpecies::ComputeShootingGuess just delegates here.
  std::function<std::pair<double, double>(double omega0_scf,
                                          double H0,
                                          const std::vector<double>& params,
                                          int tuning_index)>
      shooting_guess;
};

/** The historical built-in: V = exp(-lambda*phi) * ((phi-B)^alpha + A),
 *  params = [lambda, alpha, A, B]. */
ScalarFieldPotential DefaultScalarFieldPotential();

/** Axion bundle: V = m^2 f^2 (1 - cos(phi/f))^n, params = [m, f, n, Theta_ini]
 *  (m in 1/Mpc; f, phi in reduced-Planck units; n >= 1; Theta_ini in (0, pi)).
 *  ddV is written with sin^2 = u(2-u) so it stays regular at the minimum u -> 0.
 *  The shooting guess assumes a frozen field today: V(Theta_ini*f) = 3 H0^2 Omega. */
ScalarFieldPotential AxionScalarFieldPotential();
