// The three-zone baryon-clumping solve (H0 Olympics 2107.10291 sec. 2.4.1) must meet
// its three constraints, give a homogeneous plasma at b = 0, and refuse a shape no
// density distribution has.

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "thermodynamics_module.h"

int main() {
  // (f_V^2, Delta_1, Delta_2): the paper's model M1, and a second shape.
  const double shapes[2][3] = {{1. / 3., 0.1, 1.}, {0.5, 0.3, 1.2}};
  for (const auto& shape : shapes) {
    for (double b : {0.1, 0.5, 2.}) {
      double volume = 0., baryons = 0., second_moment = 0.;
      for (const auto& zone : BaryonClumpingZones(b, shape[0], shape[1], shape[2])) {
        assert(zone.f_V >= 0. && zone.Delta >= 0.);
        volume        += zone.f_V;
        baryons       += zone.f_V * zone.Delta;
        second_moment += zone.f_V * zone.Delta * zone.Delta;
      }
      assert(std::fabs(volume - 1.) < 1e-12);
      assert(std::fabs(baryons - 1.) < 1e-12);
      assert(std::fabs(second_moment - 1. - b) < 1e-12);
    }
  }

  // b = 0 in M1: zone 1 empties and the rest sit at the mean density.
  const auto homogeneous = BaryonClumpingZones(0., 1. / 3., 0.1, 1.);
  assert(homogeneous[0].f_V == 0.);
  assert(std::fabs(homogeneous[2].Delta - 1.) < 1e-12);

  // A zone 1 denser than the mean cannot carry a variance with zone 2 at the mean.
  bool refused = false;
  try {
    BaryonClumpingZones(0.5, 1. / 3., 1.5, 1.);
  }
  catch (const std::runtime_error&) {
    refused = true;
  }
  assert(refused);
  return 0;
}
