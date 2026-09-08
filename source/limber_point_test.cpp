/**
 * Unit tests for where the Limber delta sits, and for the curvature factor
 * that belongs to it (issue #423).
 *
 * The Limber prescription replaces the radial function by a delta of weight
 * Int_0^inf j_l(x) dx placed at an effective multipole L. Substituting
 * U = Phi * sin_K(chi) turns the hyperspherical radial equation into the flat
 * one with chi -> sin_K(chi),
 *
 *     U'' + [q^2 - L^2 / sin_K^2(chi)] U = 0,
 *
 * so the delta sits at sin_K(chi_L) = L/q, and the WKB weight of a delta at a
 * turning point picks up cos_K(chi_L)^(-1/2) = (1 - K L^2/q^2)^(-1/4) -- the
 * *same* L in both. These tests pin exactly that: one L in all three
 * geometries, and an amplitude that belongs to the point the code returns.
 *
 * L = l + 1/2 is the Langer-corrected multipole, not the classical turning
 * point sqrt(l(l+1)); the two differ at O(1/l) and the choice is fixed by the
 * weight IPhiFlat, which is the flat moment Int j_l dx of that same
 * prescription. See
 * docs/superpowers/specs/2026-09-08-limber-effective-multipole-design.md.
 */

#include <cassert>
#include <cmath>
#include <cstdio>

#include "transfer_module.h"

namespace {

/* The curvature K itself (not sqrt(K)), for a closed model with Omega_k = -0.01
   and h = 0.67: K = (h/2997.9 * sqrt(0.01))^2 Mpc^-2. */
constexpr double kK = 4.9955e-10;

constexpr int kFlat   = 0;
constexpr int kClosed = 1;
constexpr int kOpen   = -1;

double sin_K(double chi, int sgnK, double K) {
  if (sgnK == 0)
    return chi;
  const double sqrt_absK = std::sqrt(sgnK * K);
  return (sgnK == 1) ? std::sin(sqrt_absK * chi) / sqrt_absK
                     : std::sinh(sqrt_absK * chi) / sqrt_absK;
}

double cos_K(double chi, int sgnK, double K) {
  if (sgnK == 0)
    return 1.;
  const double sqrt_absK = std::sqrt(sgnK * K);
  return (sgnK == 1) ? std::cos(sqrt_absK * chi) : std::cosh(sqrt_absK * chi);
}

/* q that puts the delta where sin_K^2(chi_L) * K = s, i.e. K L^2/q^2 = s. */
double q_for_curvature_fraction(double l, double s, int sgnK, double K) {
  const double L = l + 0.5;
  if (sgnK == 0)
    return L / 1.e4;
  return L * std::sqrt(sgnK * K / s);
}

void expect_close(double got, double want, double tol, const char* what) {
  const double rel = std::abs(got - want) / std::abs(want);
  if (!(rel < tol)) {
    std::printf("  FAIL %-46s got %.16e want %.16e  rel %.3e (tol %.1e)\n",
                what,
                got,
                want,
                rel,
                tol);
  }
  assert(rel < tol);
}

/* The delta sits at sin_K(chi_L) q = l + 1/2, in every geometry. This is the
   whole content of the issue: the closed branch used sqrt(l(l+1)) instead. */
void test_evaluation_point_is_the_langer_multipole() {
  for (double l : {40., 100., 400., 2000.}) {
    for (double s : {0.05, 0.2, 0.5}) {
      for (int sgnK : {kFlat, kClosed, kOpen}) {
        const double K                      = (sgnK == 0) ? 0. : sgnK * kK;
        const double q                      = q_for_curvature_fraction(l, s, sgnK, K);
        const TransferModule::LimberPoint p = TransferModule::LimberDeltaPoint(l, q, sgnK, K);
        expect_close(sin_K(p.tau0_minus_tau, sgnK, K) * q,
                     l + 0.5,
                     1.e-12,
                     "sin_K(chi_L) * q == l + 1/2");
      }
    }
  }
}

/* The curvature factor is cos_K evaluated at the delta's own location. This
   holds whichever L one picks, so it is an invariant of the pair and not a
   restatement of the formula: it is what catches an amplitude evaluated at a
   different multipole from the turning point. */
void test_amplitude_belongs_to_its_own_point() {
  for (double l : {40., 100., 400., 2000.}) {
    for (double s : {0.05, 0.2, 0.5}) {
      for (int sgnK : {kFlat, kClosed, kOpen}) {
        const double K                      = (sgnK == 0) ? 0. : sgnK * kK;
        const double q                      = q_for_curvature_fraction(l, s, sgnK, K);
        const TransferModule::LimberPoint p = TransferModule::LimberDeltaPoint(l, q, sgnK, K);
        expect_close(p.amplitude * p.amplitude * cos_K(p.tau0_minus_tau, sgnK, K),
                     1.,
                     1.e-12,
                     "amplitude^2 * cos_K(chi_L) == 1");
      }
    }
  }
  /* Direction, so a factor of one could never pass the invariant by accident. */
  const double q_closed = q_for_curvature_fraction(100., 0.2, kClosed, kK);
  const double q_open   = q_for_curvature_fraction(100., 0.2, kOpen, -kK);
  assert(TransferModule::LimberDeltaPoint(100., q_closed, kClosed, kK).amplitude > 1.);
  assert(TransferModule::LimberDeltaPoint(100., q_open, kOpen, -kK).amplitude < 1.);
  assert(TransferModule::LimberDeltaPoint(100., 100.5 / 1.e4, kFlat, 0.).amplitude == 1.);
}

/* Both curved geometries must reduce to the flat one as K -> 0. Open already
   did; closed did not, and the offset survived the limit -- a discontinuity in
   the observable across Omega_k = 0. */
void test_curved_geometries_are_continuous_with_flat() {
  const double l = 40.;
  const double L = l + 0.5;
  for (double s : {1.e-8, 1.e-10, 1.e-12}) {
    for (int sgnK : {kClosed, kOpen}) {
      const double K                      = sgnK * kK;
      const double q                      = q_for_curvature_fraction(l, s, sgnK, K);
      const double chi_lat                = L / q;  // what the flat branch would return
      const TransferModule::LimberPoint p = TransferModule::LimberDeltaPoint(l, q, sgnK, K);
      /* sin_K(chi) = chi (1 - K chi^2/6 + ...), so the deviation is s/6. */
      expect_close(p.tau0_minus_tau, chi_lat, 10. * s, "chi_L -> flat as K -> 0");
      expect_close(p.amplitude, 1., 10. * s, "amplitude -> 1 as K -> 0");
    }
  }
}

/* In a closed universe nu = q/sqrt(K) is an integer >= l+1, so the smallest
   representable mode has L/nu = (l+1/2)/(l+1) < 1 and the point exists. The
   classical turning point sqrt(l(l+1)) also stayed below nu, so this is a
   property the fix must not break. */
void test_closed_edge_mode_is_representable() {
  for (double l : {2., 40., 2000.}) {
    const double nu                     = l + 1.;
    const double q                      = nu * std::sqrt(kK);
    const TransferModule::LimberPoint p = TransferModule::LimberDeltaPoint(l, q, kClosed, kK);
    assert(std::isfinite(p.tau0_minus_tau));
    assert(std::isfinite(p.amplitude));
    assert(p.tau0_minus_tau > 0.);
    /* Below the equator: asin returns the near turning point. */
    assert(p.tau0_minus_tau < 0.5 * _PI_ / std::sqrt(kK));
    expect_close(sin_K(p.tau0_minus_tau, kClosed, kK) * q, l + 0.5, 1.e-12, "edge mode: L");
  }
}

/* The closed-space delta uses only the near turning point; asin cannot return
   the far one at pi R - chi_L, so the Limber lensing transfer is incomplete once
   the source support crosses the equator. Thresholds below are chi_max read off
   the background table of each model (h = 0.6732, Planck-like), against the
   equator (pi/2) R: -0.2 is inside, -0.3 is not. */
void test_equator_predicate() {
  const double h = 0.6732;
  auto K_of      = [h](double omega_k) {
    const double sqrt_K = h / 2997.9 * std::sqrt(std::abs(omega_k));
    return sqrt_K * sqrt_K;
  };

  struct Case {
    double omega_k;
    double chi_max;
    bool crosses;
  };
  const Case cases[] = {
      {-0.05, 14345., false},
      {-0.1, 14528., false},
      {-0.2, 14928., false},
      {-0.3, 15384., true},
      {-0.5, 15900., true},
  };
  for (const auto& c : cases) {
    const bool got =
        TransferModule::LimberSourceCrossesEquator(c.chi_max, kClosed, K_of(c.omega_k));
    if (got != c.crosses)
      std::printf("  FAIL Omega_k=%.2f chi_max=%.0f: got %d want %d\n",
                  c.omega_k,
                  c.chi_max,
                  static_cast<int>(got),
                  static_cast<int>(c.crosses));
    assert(got == c.crosses);
  }

  /* Flat and open have no equator, at any distance or curvature. */
  assert(!TransferModule::LimberSourceCrossesEquator(1.e9, kFlat, 0.));
  assert(!TransferModule::LimberSourceCrossesEquator(1.e9, kOpen, -K_of(0.5)));
}

/* The full-Limber CMB lensing grid is a pure geometric progression anchored at
   q_min, so where it starts fixes every node. A closed universe has modes down
   to nu = 3 (q = 3 sqrt(K)), and anchoring there made the whole grid slide with
   K and never converge to the flat grid as K -> 0. Below the flat truncation
   scale those nodes carry no integrand anyway -- the Limber transfer vanishes
   for q < l_switch_limber/tau0 -- so the anchor is clamped. */
void test_limber_grid_anchor_converges_to_flat() {
  const double h          = 0.6732;
  const double tau0       = 14160.;      // Mpc, Planck-like
  const double k_min_flat = 0.1 / tau0;  // k_min_tau0/tau0
  auto K_of               = [h](double omega_k) {
    const double sqrt_K = h / 2997.9 * std::sqrt(std::abs(omega_k));
    return sqrt_K * sqrt_K;
  };

  /* Flat is the reference and must not move. */
  assert(TransferModule::LimberGridQMin(k_min_flat, k_min_flat, kFlat, 0.) == k_min_flat);

  /* Curvature large enough to matter: the physical nu = 3 edge still wins, so
     nothing changes for any cosmology anyone runs. */
  for (double omega_k : {-0.01, -0.05, -0.2}) {
    const double K     = K_of(omega_k);
    const double k_min = std::sqrt(8. * K);  // the nu = 3 mode
    expect_close(TransferModule::LimberGridQMin(k_min, k_min_flat, kClosed, K),
                 3. * std::sqrt(K),
                 1.e-12,
                 "closed anchor stays at 3 sqrt(K)");
    assert(TransferModule::LimberGridQMin(k_min, k_min_flat, kClosed, K) > k_min_flat);
  }

  /* Approaching flat: the anchor must land on the flat one, not slide to zero. */
  for (double omega_k : {-1.e-6, -1.e-8, -1.e-10, -1.e-14}) {
    const double K     = K_of(omega_k);
    const double k_min = std::sqrt(8. * K);
    assert(3. * std::sqrt(K) < k_min_flat);  // the case the clamp exists for
    expect_close(TransferModule::LimberGridQMin(k_min, k_min_flat, kClosed, K),
                 k_min_flat,
                 1.e-12,
                 "closed anchor -> flat anchor as K -> 0");
  }

  /* Open already converged on its own: q_min = sqrt(k_min^2 + K) with the open
     k_min = sqrt(-K + k_min_flat^2) collapses to the flat anchor exactly. */
  for (double omega_k : {1.e-6, 1.e-10}) {
    const double K     = -K_of(omega_k);
    const double k_min = std::sqrt(-K + k_min_flat * k_min_flat);
    expect_close(TransferModule::LimberGridQMin(k_min, k_min_flat, kOpen, K),
                 k_min_flat,
                 1.e-12,
                 "open anchor is the flat anchor");
  }
}

}  // namespace

int main() {
  test_evaluation_point_is_the_langer_multipole();
  test_amplitude_belongs_to_its_own_point();
  test_curved_geometries_are_continuous_with_flat();
  test_closed_edge_mode_is_representable();
  test_equator_predicate();
  test_limber_grid_anchor_converges_to_flat();
  std::printf("limber_point_test: all tests passed\n");
  return 0;
}
