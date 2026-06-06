#include "greybody_moments.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>

int main() {
  using greybody::riemann_zeta;
  const double PI = 3.14159265358979323846;
  // Known closed forms.
  assert(std::fabs(riemann_zeta(2.0) - PI * PI / 6.0) < 1e-9);
  assert(std::fabs(riemann_zeta(4.0) - std::pow(PI, 4) / 90.0) < 1e-9);
  // Apery's constant.
  assert(std::fabs(riemann_zeta(3.0) - 1.2020569031595942854) < 1e-9);
  // Non-integer argument vs reference value zeta(2.5) = 1.3414872573...
  assert(std::fabs(riemann_zeta(2.5) - 1.3414872573009741) < 1e-7);
  // Negative argument: exercises the n=21 fallback branch. zeta(-1) = -1/12.
  assert(std::fabs(riemann_zeta(-1.0) - (-1.0 / 12.0)) < 1e-12);
  std::printf("greybody zeta tests passed\n");

  // Round-trip: for each physical (alpha,x,q0), compute its moments, feed
  // (r,M2,M3) back into the solver, and confirm it recovers consistent moments.
  // Two samples at different alpha guard against an alpha-specific shortcut.
  {
    using greybody::GreyBodyParams;
    const double samples[][3] = {{2.5, 0.8, 1.3}, {1.5, 1.2, 0.9}};  // alpha, x, q0
    for (const auto& s : samples) {
      GreyBodyParams forward = GreyBodyParams::FromDirect(s[0], s[1], s[2]);
      double M2, M3, M4;
      forward.moments(M2, M3, M4);
      double r = M2 * M4 / (M3 * M3);

      GreyBodyParams solved = GreyBodyParams::FromMoments(r, M2, M3);
      double m2b, m3b, m4b;
      solved.moments(m2b, m3b, m4b);

      assert(std::fabs(solved.alpha() - s[0]) < 1e-6);
      assert(std::fabs(m2b - M2) / M2 < 1e-8);
      assert(std::fabs(m3b - M3) / M3 < 1e-8);
      assert(std::fabs(m4b - M4) / M4 < 1e-6);
    }
  }
  std::printf("greybody moment round-trip passed\n");

  // The defensive guard: non-positive direct parameters (alpha,x,q0) would
  // otherwise produce inf/nan downstream, so FromDirect must reject them.
  {
    using greybody::GreyBodyParams;
    const double bad[][3] =
        {{0., 1., 1.}, {-1., 1., 1.}, {1., 0., 1.}, {1., -1., 1.}, {1., 1., 0.}, {1., 1., -1.}};
    for (const auto& b : bad) {
      bool threw = false;
      try {
        GreyBodyParams::FromDirect(b[0], b[1], b[2]);
      }
      catch (const std::invalid_argument&) {
        threw = true;
      }
      assert(threw);
    }
  }
  std::printf("greybody direct-param guard passed\n");

  // FromMoments takes log(M2) and divides by M3, so non-positive moments must be
  // rejected (r = 2 is a valid ratio; only the moments are bad here).
  {
    using greybody::GreyBodyParams;
    const double bad[][3] = {{2., 0., 1.}, {2., -1., 1.}, {2., 1., 0.}, {2., 1., -1.}};
    for (const auto& b : bad) {
      bool threw = false;
      try {
        GreyBodyParams::FromMoments(b[0], b[1], b[2]);
      }
      catch (const std::invalid_argument&) {
        threw = true;
      }
      assert(threw);
    }
  }
  std::printf("greybody moment-param guard passed\n");
  return 0;
}
