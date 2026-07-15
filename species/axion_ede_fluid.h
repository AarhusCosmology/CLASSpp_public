#pragma once

#include "fluid.h"

/**
 * Pheno-axion early-dark-energy fluid (AxiCLASS's ede_parametrization = pheno_axion;
 * Poulin et al. 1806.10608, 1811.04083, 1905.12618). Background: sigmoid equation of
 * state w(a) = w_i + (w_f - w_i)/(1 + (a_c/a)^r), r = 3(w_f - w_i)/nu, with a
 * closed-form density integral. Perturbations: true-fluid delta/theta with the
 * k- and a-dependent GDM effective sound speed
 *   cs2(k,a) = [2 a^2 (n-1) wbar^2 + k^2] / [2 a^2 (n+1) wbar^2 + k^2],
 *   wbar(a) = omega_axion * a^(-3(n-1)/(n+1)),
 * where omega_axion is derived from (a_c, Omega_fld_ac, Theta_i, n) at
 * SetBackgroundModule time (needs the full species budget for E(a_c)).
 * PPF does not apply (w > -1 for all a > 0, and the GDM cs2 must enter the true
 * fluid equations), so this class is always built as a non-PPF fluid.
 */
class AxionEDEFluid : public FluidSpecies {
 public:
  AxionEDEFluid(const background& pba,
                double omega0_fld,
                double a_c,
                double n_axion,
                double nu,
                double w_i,
                double w_f,
                double theta_i);

  // ── closed-form pieces (unit-testable without a background pipeline) ──────
  /** w_f = (n-1)/(n+1). */
  static double WFinal(double n);
  /** int_a^1 da' 3(1+w(a'))/a' for the sigmoid w(a). Overflow-safe (log1p form). */
  static double Integral3OnePlusWOverA(double a, double a_c, double nu, double w_i, double w_f);
  /** Omega0_fld from Omega_fld_ac = rho_fld(a_c)/rho_crit0 (and inverse). */
  static double OmegaZeroFromOmegaAc(
      double omega_ac, double a_c, double nu, double w_i, double w_f);
  static double OmegaAcFromOmegaZero(double omega0, double a_c, double nu, double w_i, double w_f);

  void ComputeWFld(double a,
                   double* w_fld,
                   double* dw_over_da_fld,
                   double* integral_fld) const override;
  double Cs2(double k2, double a) const override;
  bool ReachesPhantomDivide() const override {
    // w(a) = w_i + (w_f - w_i)/(1 + (a_c/a)^r) > w_i = -1 strictly for all a > 0;
    // the exact -1 at the a = 0 asymptote lies outside the integration domain.
    return false;
  }
  bool HyrecCplApproximation(double* /*w0*/, double* /*wa*/) const override {
    // No CPL pair represents a frozen->dilution sigmoid: HyRec would extrapolate
    // rho_DE ~ a^(-3(1+w_f)) from today's density to all redshifts and destroy the
    // recombination-era expansion rate. Callers must use RECFAST (reads the true
    // background) or implement true-H(z) feeding for HyRec.
    return false;
  }
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  void SetBackgroundModule(const BackgroundModule* bgm) override;

  double a_c() const {
    return a_c_;
  }
  double n_axion() const {
    return n_axion_;
  }
  double omega_axion() const {
    return omega_axion_;
  }
  /** Unit-test hook: Cs2 needs omega_axion_, normally derived in SetBackgroundModule. */
  void SetOmegaAxionForTest(double omega_axion) {
    omega_axion_ = omega_axion;
  }

 private:
  /** Port of AxiCLASS background.c 982-1023 (Eqs. 27/28/30 of 1806.10608):
   *  derive m_fld, alpha_fld, omega_axion from (a_c, Omega_fld_ac, Theta_i, n)
   *  and the species budget at a_c. Called from SetBackgroundModule. */
  void DeriveAxionScales();

  double a_c_     = 0.;
  double n_axion_ = 3.;
  double nu_      = 1.;
  double w_i_     = -1.;
  double w_f_     = 0.5;
  double theta_i_ = 0.;
  const background& pba_ref_;  // FluidSpecies::pba_ is private; keep our own ref for H0
  double omega_ac_    = 0.;    // rho_fld(a_c)/rho_crit0
  double m_fld_       = 0.;    // axion mass in units of H0 (derived, diagnostic)
  double alpha_fld_   = 0.;    // decay constant in reduced-Planck units (derived, diagnostic)
  double omega_axion_ = 0.;    // characteristic oscillation frequency today [1/Mpc]
};
