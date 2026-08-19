#ifndef __EVOETD__
#define __EVOETD__

#include "common.h"

/**
 * Exponential Rosenbrock (ETDRK4, Cox-Matthews) integrator for systems whose
 * stiffness lives on the Jacobian DIAGONAL.
 *
 * L = diag(J), frozen over the step; N(y) = f(y) - L y, reformed at every stage.
 * The phi functions phi0..phi3 are applied COMPONENTWISE, so each stage is exact for
 * the scalar y' = d y + g at constant g -- the stiff part is integrated rather than
 * resolved.
 *
 * A species offers its diagonal through BaseSpecies::BackgroundDerivsDiagonal or
 * ::PerturbDerivsDiagonal; the default is zero, where phi0 = 1, phi1 = 1, phi2 = 1/2,
 * phi3 = 1/6 and the scheme reduces exactly to its explicit counterpart. So this is
 * not a decaying-neutrino solver -- it is an evolver that exploits a diagonal when
 * one is offered, and costs a little extra per step when none is.
 *
 * Because N is reformed as f - Ly at every stage, the fixed point is exact for ANY
 * L: at f(y*) = 0, N(y*) = -L y* and h phi1(hL) L = e^{hL} - 1, so the update returns
 * y*. A wrong diagonal costs efficiency and conditioning, never which solution is
 * converged to. A wrong SIGN is the dangerous case: e^{hL} then grows and the
 * explicit remainder has to cancel it.
 *
 * TWO THINGS NOT TO "SIMPLIFY".
 *   - Order 4. A two-stage ETD2RK was built and measured an order of magnitude
 *     SLOWER than rkdp45: accuracy-limited everywhere, so it never spent the
 *     stability it was built for.
 *   - The EXPONENTIAL dense output. When a step jumps a relaxation layer of width
 *     1/|diag| -- the whole point of the method -- a cubic Hermite through the two
 *     endpoints cannot see the layer, and the SAMPLED table ends up far less accurate
 *     than the trajectory that produced it. In the perturbations the sampled sources
 *     are the entire deliverable. The extension here interpolates N quadratically
 *     through its own stage samples and integrates that exactly, is exact for
 *     constant N at every theta, meets the accepted step at theta = 1, and costs no
 *     RHS evaluations.
 *
 * The error estimate is an embedded third-order companion reusing the endpoint RHS;
 * the accepted trajectory stays Cox-Matthews order 4.
 *
 * Derivation, the stability/cost measurements behind the order-4 and dense-output
 * choices, and the rejected alternatives (ndf15, phi-lifted DP54, published EPIRK
 * 5(4) pairs) are in docs/superpowers/specs/2026-08-06-etdrk4-perturbations-design.md
 * and docs/superpowers/plans/2026-08-06-etdrk4-perturbations.md.
 *
 * `derivs_diagonal` may be null, in which case the diagonal is taken as zero
 * throughout. minimum_variation, evaluate_timescale and timestep_over_timescale are
 * part of the shared CLASS evolver signature and unused here, as in rkdp45.
 */
void evolver_etd(
    void (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double minimum_variation,
    void (*evaluate_timescale)(double x, void* parameters_and_workspace, double* timescale),
    double timestep_over_timescale,
    double* x_sampling,
    int x_size,
    void (*output)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace),
    void (*print_variables)(double x, double y[], double dy[], void* parameters_and_workspace),
    void (*derivs_diagonal)(double x, double* y, double* diag, void* parameters_and_workspace));

/**************************************************************/

#endif
