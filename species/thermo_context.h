#pragma once
#include <vector>

#include "common.h"

class BackgroundModule;
struct precision;

/**
 * Context passed to thermodynamics hooks in BaseSpecies.
 * Bundles the thermodynamics table, coordinate arrays, and related metadata
 * so species can perform multi-row operations (spline, integrate, etc.)
 * without parameter explosion.
 */
struct ThermoTableContext {
  std::vector<double>& table;           /**< thermodynamics_table_ (writable) */
  const std::vector<double>& tau_table; /**< conformal time for each row */
  const std::vector<double>& z_table;   /**< redshift for each row */
  int th_size;                          /**< number of columns (stride) */
  const BackgroundModule* background_module;
  const precision* ppr;
  ErrorMsg& error_message;
};
