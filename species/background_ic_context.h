#pragma once

/**
 * Context passed to BaseSpecies::SetBackgroundInitialConditions().
 * a_rel   : a / a_today at the initial time.
 * a_ini   : absolute initial a (= a_rel * a_today); the value the module stores
 *           at index_bi_a_, exposed so species need not know the {B} layout.
 * rho_rad : total radiation density at the initial time (units of H0^2), needed
 *           by ScalarField attractor ICs.
 * pvecback_integration : the integration vector to fill.
 */
struct BackgroundICContext {
  double a_rel                 = 0.;
  double a_ini                 = 0.;
  double rho_rad               = 0.;
  double* pvecback_integration = nullptr;
};
