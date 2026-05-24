#pragma once
#include <string>

struct background;

/** Context for ComputeShootingResidual: a freshly-COMPUTED cosmology's "today" state. */
struct ShootingResidualContext {
  const background* pba  = nullptr;
  const double* bg_today = nullptr;  // last row of BackgroundModule::background_table_
};

/** One scalar fittable target a species owns. (All shooting targets are scalar: the
 *  multi-value dncdm input was removed in favour of dot-syntax, one flavor per instance.) */
struct ShootingTarget {
  std::string target_name;    // input vocabulary, e.g. "Omega_dcdmdr"
  std::string unknown_param;  // fc param this target varies, e.g. "Omega_ini_dcdm"
  double target_value = 0.;   // requested value
};
