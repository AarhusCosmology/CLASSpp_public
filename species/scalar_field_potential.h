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
