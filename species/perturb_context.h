#pragma once

/**
 * Pre-computed scalar perturbation quantities, populated by
 * PerturbationsModule::perturb_derivs_member() before calling species.
 * Stored in perturb_workspace so every species can access shared state.
 */
struct PerturbScalarContext {
  double k = 0., k2 = 0.;
  double a = 0., a2 = 0., a_prime_over_a = 0.;
  double metric_continuity = 0.;
  double metric_euler = 0.;
  double metric_shear = 0.;
  double metric_ufa_class = 0.;
  double cotKgen = 0., s2_squared = 1.;
  /** photon delta/theta (RSA-corrected) */
  double delta_g = 0., theta_g = 0.;
  /** baryon delta/theta */
  double delta_b = 0., theta_b = 0.;
  /** IDR delta/theta (RSA-corrected) */
  double delta_idr = 0., theta_idr = 0.;
  /** 4/3 * rho_g / rho_b, photon-baryon momentum ratio */
  double R = 0.;
  /** baryon sound speed squared */
  double cb2 = 0.;
  /** baryon pressure perturbation / rho */
  double delta_p_b_over_rho_b = 0.;
};
