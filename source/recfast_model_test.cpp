#include "recfast_model.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>

#include "parser.h"
#include "precision.h"

namespace {

/* Planck-2018-like, as in hyrec_model_test.cpp. */
constexpr double kT_cmb = 2.7255;
constexpr double kYHe   = 0.2454;
constexpr double kNH0   = 0.1904; /* hydrogen nuclei today, m^-3 */

/* The RECFAST constants exactly as thermodynamics_module.cpp derives them. */
recombination MakeRecombination(const precision& ppr) {
  recombination preco    = {};
  preco.YHe              = kYHe;
  preco.fHe              = kYHe / (_not4_ * (1. - kYHe));
  preco.fu               = ppr.recfast_fudge_H + ppr.recfast_delta_fudge_H;
  const double Lalpha    = 1. / _L_H_alpha_;
  const double Lalpha_He = 1. / _L_He_2p_;
  preco.CDB              = _h_P_ * _c_ * (_L_H_ion_ - _L_H_alpha_) / _k_B_;
  preco.CDB_He           = _h_P_ * _c_ * (_L_He1_ion_ - _L_He_2s_) / _k_B_;
  preco.CR               = 2. * _PI_ * (_m_e_ / _h_P_) * (_k_B_ / _h_P_);
  preco.CK               = pow(Lalpha, 3) / (8. * _PI_);
  preco.CK_He            = pow(Lalpha_He, 3) / (8. * _PI_);
  preco.CL               = _c_ * _h_P_ / (_k_B_ * Lalpha);
  preco.CL_He            = _c_ * _h_P_ / (_k_B_ / _L_He_2s_);
  preco.Bfact            = _h_P_ * _c_ * (_L_He_2p_ - _L_He_2s_) / _k_B_;
  return preco;
}

RecombinationState StateAt(double z, double x_H, double x_He, bool helium_ode = true) {
  RecombinationState s = {};
  s.z                  = z;
  s.x_H                = x_H;
  s.x_He               = x_He;
  s.x                  = x_H + kYHe / (_not4_ * (1. - kYHe)) * x_He;
  s.n_H                = kNH0 * (1. + z) * (1. + z) * (1. + z);
  s.T_rad              = kT_cmb * (1. + z);
  s.T_mat              = s.T_rad;
  s.H = 2.2e-18 * std::sqrt(0.31 * std::pow(1. + z, 3) + 9.2e-5 * std::pow(1. + z, 4));
  s.hydrogen_frozen = false;
  s.helium_ode      = helium_ode;
  return s;
}

const EnergyDeposition kNoInjection = {0., 0., 0., 0., 0.};

precision Parse(const std::string& calibration) {
  FileContent fc;
  fc.set("recfast_calibration", calibration);
  precision ppr;
  ppr.parse(fc);
  return ppr;
}

/* CAMB 2's He I rate correction (recfast.f90), written out independently of the model. */
double CambHeRateFactor(double z) {
  const double u = (z - 1968.1134791219417) / 774.0312440582285;
  return 1. + (0.07805480599148856 + 0.899368710079757 * u - 3.063896868095767 * u * u) /
                  (1. + 3.752520187028219 * u * u + 18.126758791787946 * u * u * u * u);
}

/* The calibration sets all seven fit constants and the switch as a group, 'cosmorec'
   restores the compiled defaults, 'planck' gives RECFAST 1.5.2's own values, and an
   unknown name is refused. */
void test_calibration_sets_the_group() {
  const precision defaults;
  const precision cosmorec = Parse("cosmorec");
  assert(cosmorec.recfast_fudge_He == defaults.recfast_fudge_He);
  assert(cosmorec.recfast_AGauss1 == defaults.recfast_AGauss1);
  assert(cosmorec.recfast_AGauss2 == defaults.recfast_AGauss2);
  assert(cosmorec.recfast_zGauss1 == defaults.recfast_zGauss1);
  assert(cosmorec.recfast_zGauss2 == defaults.recfast_zGauss2);
  assert(cosmorec.recfast_wGauss1 == defaults.recfast_wGauss1);
  assert(cosmorec.recfast_wGauss2 == defaults.recfast_wGauss2);
  assert(cosmorec.recfast_He_rate_correction && defaults.recfast_He_rate_correction);
  assert(cosmorec.recfast_fudge_He == 0.8472367977);
  assert(cosmorec.recfast_AGauss1 == -0.1395272483);
  assert(cosmorec.recfast_AGauss2 == 0.0729891952);
  assert(cosmorec.recfast_zGauss1 == 7.2813061282);
  assert(cosmorec.recfast_zGauss2 == 6.7667038679);
  assert(cosmorec.recfast_wGauss1 == 0.1638966410);
  assert(cosmorec.recfast_wGauss2 == 0.2785834127);

  const precision planck = Parse("planck");
  assert(planck.recfast_fudge_He == 0.86);
  assert(planck.recfast_AGauss1 == -0.14);
  assert(planck.recfast_AGauss2 == 0.079);
  assert(planck.recfast_zGauss1 == 7.28);
  assert(planck.recfast_zGauss2 == 6.73);
  assert(planck.recfast_wGauss1 == 0.18);
  assert(planck.recfast_wGauss2 == 0.33);
  assert(!planck.recfast_He_rate_correction);
  /* The hydrogen fudge is 1.125 in both calibrations. */
  assert(planck.recfast_fudge_H == defaults.recfast_fudge_H);
  assert(planck.recfast_delta_fudge_H == defaults.recfast_delta_fudge_H);

  bool threw = false;
  try {
    Parse("CosmoRec2");
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);
}

/* A single constant given alongside the calibration wins over the group. */
void test_individual_keys_override_the_group() {
  FileContent fc;
  fc.set("recfast_calibration", std::string("cosmorec"));
  fc.set("recfast_fudge_He", std::string("0.9"));
  fc.set("recfast_He_rate_correction", std::string("no"));
  precision ppr;
  ppr.parse(fc);
  assert(ppr.recfast_fudge_He == 0.9);
  assert(!ppr.recfast_He_rate_correction);
  assert(ppr.recfast_AGauss1 == -0.1395272483);
}

/* Toggling only the switch scales dx_He/dz by CAMB's F_He inside 1500 < z < 3000 while
   helium is evolved, and nowhere else; hydrogen never sees it. */
void test_he_rate_factor_matches_camb() {
  precision off;
  precision on;
  off.recfast_He_rate_correction = false;
  on.recfast_He_rate_correction  = true;
  const recombination preco      = MakeRecombination(off);
  const RecfastModel model_off(&off, &preco);
  const RecfastModel model_on(&on, &preco);

  for (double z : {1600., 1800., 1968., 2200., 2500., 2900.}) {
    const RecombinationState s    = StateAt(z, 0.9999, 0.6);
    const IonisationDerivatives a = model_off.Derivatives(s, kNoInjection, 0.);
    const IonisationDerivatives b = model_on.Derivatives(s, kNoInjection, 0.);
    assert(a.dx_He_dz != 0.);
    assert(std::fabs(b.dx_He_dz / a.dx_He_dz / CambHeRateFactor(z) - 1.) < 1e-12);
    assert(b.dx_H_dz == a.dx_H_dz);
  }
  /* Outside the redshift window, above x_He = 0.98, or with helium held at Saha. */
  for (const RecombinationState& s : {StateAt(1400., 0.999, 0.05),
                                      StateAt(3100., 1., 0.9),
                                      StateAt(2200., 0.9999, 0.99),
                                      StateAt(2200., 0.9999, 0.6, false)}) {
    assert(model_on.Derivatives(s, kNoInjection, 0.).dx_He_dz ==
           model_off.Derivatives(s, kNoInjection, 0.).dx_He_dz);
  }
}

}  // namespace

int main() {
  test_calibration_sets_the_group();
  test_individual_keys_override_the_group();
  test_he_rate_factor_matches_camb();

  std::printf("recfast model tests passed\n");
  return 0;
}
