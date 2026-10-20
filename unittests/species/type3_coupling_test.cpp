// Unit test for the Type-3 (arXiv:1604.04222) synchronous-gauge momentum-transfer
// coupling, exercised through the two pure-arithmetic free functions that
// AddCouplingDerivs calls:
//
//   Type3CouplingDeltaPhiPrime(beta, phi_prime_bg, theta_cdm)
//   Type3CouplingDeltaThetaCdm(beta, k2, a_prime_over_a, rho_cdm, a2, Zbar,
//                              dV, phi, phi_prime, phi_prime_bg, theta_cdm)
//
// Two guarantees:
//   1. STRUCTURAL: every term carries an explicit `beta` factor, so both
//      functions return exactly 0 at beta = 0 (the coupling is a no-op, and the
//      composite reduces to its free-streaming children).
//   2. NUMERIC: at beta = -0.5 each function equals a value computed BY HAND for
//      a fixed set of inputs (see the derivation in the asserts below).
//
// No perturb_parameters_and_workspace is built here: AddCouplingDerivs is a thin
// wrapper that only reads ctx/pvecback/y slots and forwards them to these
// functions, so testing the functions directly locks the physics.

#include <cassert>
#include <cmath>
#include <cstdio>

#include "type3_species.h"

int main() {
  // ── Guarantee 1: beta = 0 => exactly zero (no coupling) ────────────────────
  {
    const double beta = 0.;
    // DeltaPhiPrime is identically 0 at beta=0 for any inputs.
    assert(Type3CouplingDeltaPhiPrime(beta, /*phi_prime_bg=*/0.4, /*theta_cdm=*/2.0) == 0.);
    // DeltaThetaCdm is identically 0 at beta=0 for any inputs (every term * beta).
    const double dtheta0 = Type3CouplingDeltaThetaCdm(beta,
                                                      /*k2=*/3.0,
                                                      /*a_prime_over_a=*/0.5,
                                                      /*rho_cdm=*/1.0,
                                                      /*a2=*/0.25,
                                                      /*Zbar=*/-0.8,
                                                      /*dV=*/0.2,
                                                      /*phi=*/0.3,
                                                      /*phi_prime=*/0.1,
                                                      /*phi_prime_bg=*/0.4,
                                                      /*theta_cdm=*/2.0);
    assert(dtheta0 == 0.);
  }

  // ── Guarantee 2: beta = -0.5, hand-computed value ──────────────────────────
  // Inputs (chosen so a=0.5 => a2=0.25 and Zbar = -phi_prime_bg/a = -0.8):
  //   beta = -0.5,  phi_prime_bg = 0.4,  theta_cdm = 2.0,  k2 = 3.0,
  //   a_prime_over_a = 0.5,  rho_cdm = 1.0,  a2 = 0.25,  Zbar = -0.8,
  //   dV = 0.2,  phi = 0.3,  phi_prime = 0.1.
  const double beta           = -0.5;
  const double phi_prime_bg   = 0.4;
  const double theta_cdm      = 2.0;
  const double k2             = 3.0;
  const double a_prime_over_a = 0.5;
  const double rho_cdm        = 1.0;
  const double a2             = 0.25;
  const double Zbar           = -0.8;  // -phi_prime_bg / a, with a = 0.5
  const double dV             = 0.2;
  const double phi            = 0.3;
  const double phi_prime      = 0.1;

  // DeltaPhiPrime = 2 * (beta/(1-2 beta)) * phi_prime_bg * theta_cdm
  //   1-2 beta = 2;  beta/(1-2 beta) = -0.25
  //   = 2 * (-0.25) * 0.4 * 2.0 = -0.4   (exact)
  const double dphi_prime = Type3CouplingDeltaPhiPrime(beta, phi_prime_bg, theta_cdm);
  assert(std::fabs(dphi_prime - (-0.4)) < 1e-13);

  // DeltaThetaCdm hand calculation:
  //   one_m_2b = 2,  Zbar^2 = 0.64
  //   D = 3*rho_cdm*a2 - 2*beta*Zbar^2 = 0.75 - (-0.64) = 1.39
  //   term1 = -k2*2beta*(one_m_2b*Zbar^2*phi_prime + dV*phi)/(one_m_2b*D)
  //         = -3*(-1)*(2*0.64*0.1 + 0.2*0.3)/(2*1.39)
  //         = 3*(0.128+0.06)/2.78 = 0.564/2.78 = 141/695
  //   term2 = 4*beta*dV*phi_prime_bg*(k2*phi/phi_prime_bg - 2*beta*theta_cdm)/(one_m_2b*D)
  //         = (-2*0.08)*(2.25 - (-2.0))/2.78 = (-0.16*4.25)/2.78 = -0.68/2.78 = -34/139
  //   term3 = -(6*beta*aH*Zbar^2 + 4*beta*dV*phi_prime_bg)*theta_cdm/D
  //         = -(-0.96 + -0.16)*2.0/1.39 = (1.12*2.0)/1.39 = 2.24/1.39 = 224/139
  //   total = 141/695 - 34/139 + 224/139 = (141 - 170 + 1120)/695 = 1091/695
  //         = 1.5697841726618705...
  const double dtheta          = Type3CouplingDeltaThetaCdm(beta,
                                                            k2,
                                                            a_prime_over_a,
                                                            rho_cdm,
                                                            a2,
                                                            Zbar,
                                                            dV,
                                                            phi,
                                                            phi_prime,
                                                            phi_prime_bg,
                                                            theta_cdm);
  const double dtheta_expected = 1091.0 / 695.0;
  assert(std::fabs(dtheta - dtheta_expected) < 1e-12);

  std::printf(
      "type3 coupling free-function test passed "
      "(beta=0 no-op; beta=-0.5 dphi'=%.15g, dtheta=%.15g)\n",
      dphi_prime,
      dtheta);
  return 0;
}
