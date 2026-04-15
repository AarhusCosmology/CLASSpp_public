#pragma once
#include <vector>

#include "common.h"  // for file_format, _TRUE_, _FALSE_, class_store_* macros

// Forward declarations
struct perturb_workspace;
struct perturb_vector;
struct precision;
class PerturbationsModule;

// ── existing PerturbScalarContext (unchanged) ────────────────────────────────

/**
 * Pre-computed scalar perturbation quantities, populated by
 * PerturbationsModule::perturb_derivs_member() before calling species.
 * Stored in perturb_workspace so every species can access shared state.
 */
struct PerturbScalarContext {
  double k = 0., k2 = 0.;
  double a = 0., a2 = 0., a_prime_over_a = 0.;
  double metric_continuity = 0.;
  double metric_euler      = 0.;
  double metric_shear      = 0.;
  double metric_ufa_class  = 0.;
  double cotKgen = 0., s2_squared = 1.;
  /** photon delta/theta (RSA-corrected) */
  double delta_g = 0., theta_g = 0.;
  /** photon shear (TCA-corrected: 16/45/dkappa*theta_g in Newtonian gauge) */
  double shear_g = 0.;
  /** baryon delta/theta */
  double delta_b = 0., theta_b = 0.;
  /** nature of idr (free-streaming or fluid) */
  int idr_nature = 0;
  /** 4/3 * rho_g / rho_b, photon-baryon momentum ratio */
  double R = 0.;
  /** baryon sound speed squared */
  double cb2 = 0.;
  /** baryon pressure perturbation / rho */
  double delta_p_b_over_rho_b = 0.;
  /** gauge: 0 = newtonian, 1 = synchronous (mirrors possible_gauges enum) */
  int gauge = 0;
};

// ── PerturbColumnWriter ──────────────────────────────────────────────────────
/**
 * Thin helper that handles both title-writing and data-writing in one pass.
 * Construct in title mode or one of two data modes; call Add() once per column.
 *
 * The two Add overloads let species use the same writer for:
 *   - WriteOutputColumns: Add(title, tp_index, active)  — value = tk[tp_index]
 *   - PrintVariables    : Add(title, value,    active)  — value passed directly
 */
class PerturbColumnWriter {
 public:
  // Title mode (used from perturb_output_titles and perturb_prepare_k_output)
  explicit PerturbColumnWriter(char* titles) : titles_(titles) {}

  // Data mode for WriteOutputColumns: value = tk[tp_index]
  PerturbColumnWriter(double* dataptr, const double* tk, int& storeidx)
      : dataptr_(dataptr), tk_(tk), storeidx_(&storeidx) {}

  // Data mode for PrintVariables: value passed directly
  PerturbColumnWriter(double* dataptr, int& storeidx) : dataptr_(dataptr), storeidx_(&storeidx) {}

  bool IsTitleMode() const {
    return titles_ != nullptr;
  }

  // For WriteOutputColumns: species pass a tp_index; writer looks up tk[tp_index]
  void Add(const char* title, int tp_index, bool active);

  // For PrintVariables: species pass the value directly
  void Add(const char* title, double value, bool active);

 private:
  char* titles_     = nullptr;
  double* dataptr_  = nullptr;
  const double* tk_ = nullptr;
  int* storeidx_    = nullptr;
};

// ── PerturbSourceContext ─────────────────────────────────────────────────────
/**
 * Context passed to BaseSpecies::FillSources().
 * Bundles workspace, addressing info, and pre-computed N-body gauge corrections.
 */
struct PerturbSourceContext {
  PerturbationsModule* p_mod = nullptr;
  perturb_workspace* ppw     = nullptr;
  int index_md               = 0;
  int index_ic               = 0;
  int index_k                = 0;
  int index_tau              = 0;
  double k                   = 0.;
  double a_rel               = 0.;
  double a2_rel              = 0.;
  double a_prime_over_a      = 0.;
  // Pre-computed N-body gauge corrections (0 when not has_Nbody_gauge_transfers)
  double theta_over_k2 = 0.;
  double theta_shift   = 0.;
};

// ── PerturbIcContext ─────────────────────────────────────────────────────────
/**
 * Context passed to BaseSpecies::ApplyInitialConditions().
 * The module pre-computes all shared IC ratios and seed values; species
 * just write to y[pv->index_pt_*].
 */
struct PerturbIcContext {
  // Density fractions
  double fracnu = 0., fracg = 0., fracb = 0., fraccdm = 0., fracidm_drmd = 0.;
  double rho_m_over_rho_r = 0., om = 0.;
  // (k*tau)^2, (k*tau)^3
  double ktau_two = 0., ktau_three = 0.;
  // Curvature factor: 1 - 3K/k^2
  double s2_squared = 1.;
  // IC seed values pre-computed by module (photon delta/theta set by module formulas)
  double delta_g_ic = 0., theta_g_ic = 0.;
  // Shared relativistic IC (= delta_g_ic for adiabatic, computed by module for isocurvatures)
  double delta_ur = 0., theta_ur = 0., shear_ur = 0., l3_ur = 0., delta_dr = 0.;
  // Synchronous metric IC
  double eta = 0.;
  // DO NOT read in ApplyInitialConditions() — filled by module AFTER the species loop
  // for the Newtonian gauge transformation. Zero during the species dispatch.
  double alpha = 0., alpha_prime = 0.;
  // Kinematics
  double k = 0., tau = 0., a = 0., a_prime_over_a = 0.;
  int index_ic = 0;
  int gauge    = 0;
  // Pointers
  perturb_workspace* ppw           = nullptr;
  const precision* ppr             = nullptr;
  const PerturbationsModule* p_mod = nullptr;
};
