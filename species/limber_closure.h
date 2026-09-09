#pragma once
#include <algorithm>
#include <cmath>

/**
 * Closure of a continuously-sourced free-streaming hierarchy at l = l_max.
 *
 * Ma & Bertschinger's truncation asserts the hyperspherical Bessel recurrence
 *
 *     s_{l+1} F_{l+1} + s_l F_{l-1} = (2l+1) [cot_K(chi)/k] F_l
 *
 * evaluated at the CURRENT free-streaming distance chi. That is correct for a
 * single coherent mode released at tau = 0, and wrong for a species injected
 * continuously, which is a superposition of shells free-streaming as
 * j_l(k(chi - chi')) over their own release times. The recurrence holds
 * pointwise in the argument and does not survive that superposition; each
 * multipole is instead dominated by the shell whose age matches its own Bessel
 * peak, so cot_K belongs at the LIMBER distance rather than the current one.
 *
 * The closure's entire content is therefore the argument, and these two
 * functions are that argument. See
 * docs/superpowers/specs/2026-09-07-dr-limber-closure-design.md for the
 * derivation and docs/superpowers/specs/2026-09-09-massive-daughter-limber-closure-design.md
 * for which species it applies to and which it does not.
 */
namespace limber_closure {

/**
 * Order at which the Bessel weight of the shell superposition sits: the exact
 * flat moment ratio
 *   L(l) = Int j_l dy / Int j_l dy/y = 2(l+1)/l [Gamma((l+1)/2)/Gamma(l/2)]^2
 *        = l + 1/2 - 3/(8l) + 3/(16 l^2) + O(l^-3),
 * whose leading term is Limber's l+1/2; the closed form is the extended-Limber
 * correction to it. Depends only on l, so callers cache it on their
 * PerturbLayout rather than recomputing it in the RHS.
 */
inline double Order(int l) {
  return 2. * (l + 1.) / l * std::exp(2. * (std::lgamma((l + 1.) / 2.) - std::lgamma(l / 2.)));
}

/**
 * cot_K evaluated at the Limber distance for order L, in the units of
 * perturb_workspace::cotKgen (that is, cot_K(chi)/k):
 *   sqrt(1 - K(L^2-1)/k^2) / L,
 * the s_l expression at non-integer order L, over L. Derived from the turning
 * points TransferModule::transfer_limber uses, so it covers all three
 * geometries at once. Reduces to 1/L when flat, and clamps to zero above
 * nu = q/sqrt(K) in a closed universe, where the turning point does not exist
 * and s_l clamps too.
 */
inline double Cot(double L, double k, double K) {
  if (K == 0.)
    return 1. / L;
  return std::sqrt(std::max(1. - K * (L * L - 1.) / (k * k), 0.)) / L;
}

}  // namespace limber_closure
