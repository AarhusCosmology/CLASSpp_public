// The transfer module's Bessel table can be kept between cosmologies only if keeping it cannot
// change a result: a table built for a larger x range must hold the same nodes with the same
// values, bit for bit, as the one a fresh computation would build.
// docs/superpowers/specs/2026-10-26-bessel-table-cache-design.md
#include "bessel_table_cache.h"

#include <cassert>
#include <cstdio>
#include <vector>

#include "hyperspherical.h"

namespace {

/* A multipole list like the transfer module's: every l at first, then growing steps. */
std::vector<int> Multipoles(int l_max) {
  std::vector<int> l;
  for (int ll = 2; ll <= 10; ll++) {
    l.push_back(ll);
  }
  int step = 2;
  while (l.back() + step < l_max) {
    l.push_back(l.back() + step);
    step = std::min(step + 2, 30);
  }
  l.push_back(l_max);
  return l;
}

constexpr double kXMin     = 1e-5;
constexpr double kSampling = 8.;
constexpr double kPhiMin   = 1e-8;

HyperInterpStruct Fixed(const std::vector<int>& l, double x_max) {
  return HyperInterpStruct::OnFixedSpacing(0,
                                           1.,
                                           static_cast<int>(l.size()),
                                           l.data(),
                                           kXMin,
                                           x_max,
                                           kSampling,
                                           l.back() + 1,
                                           kPhiMin);
}

/* Every node of `small` is a node of `large`, with the same Phi and dPhi for every l. */
void AssertSameNodes(const HyperInterpStruct& small, const HyperInterpStruct& large) {
  assert(small.delta_x == large.delta_x);
  assert(small.l == large.l);
  assert(small.chi_at_phimin == large.chi_at_phimin);
  assert(small.x_size <= large.x_size);
  for (int k = 0; k < small.l_size; k++) {
    for (int j = 0; j < small.x_size; j++) {
      assert(small.x[j] == large.x[j]);
      assert(small.phi[k * small.x_size + j] == large.phi[k * large.x_size + j]);
      assert(small.dphi[k * small.x_size + j] == large.dphi[k * large.x_size + j]);
    }
  }
}

}  // namespace

int main() {
  const std::vector<int> l = Multipoles(1200);

  /* Fixed spacing: the nodes do not depend on x_max, and the values at a shared node do not
     depend on where the table ends (the recurrences run over chunks of nodes). */
  const HyperInterpStruct reference = Fixed(l, 2000.);
  assert(reference.x.back() >= 2000.);
  for (double x_max : {2000.3, 2001.7, 2013., 2100., 2600.}) {
    AssertSameNodes(reference, Fixed(l, x_max));
  }

  /* The fitted constructor is unchanged: its last node is x_max itself. */
  const HyperInterpStruct fitted(0,
                                 1.,
                                 static_cast<int>(l.size()),
                                 l.data(),
                                 kXMin,
                                 2000.,
                                 kSampling,
                                 l.back() + 1,
                                 kPhiMin);
  assert(fitted.x.back() == 2000.);

  /* The cache serves a request its table covers, and rebuilds for anything else. */
  BesselTableCache cache;
  const auto first = cache.Get(l, kXMin, 2000., kSampling, l.back() + 1, kPhiMin);
  assert(cache.Get(l, kXMin, 1900., kSampling, l.back() + 1, kPhiMin) == first);
  assert(cache.Get(l, kXMin, 2050., kSampling, l.back() + 1, kPhiMin) == first); /* headroom */
  AssertSameNodes(reference, *first);

  const auto wider = cache.Get(l, kXMin, 2400., kSampling, l.back() + 1, kPhiMin);
  assert(wider != first);
  AssertSameNodes(reference, *wider);

  const std::vector<int> other = Multipoles(1300);
  assert(cache.Get(other, kXMin, 2000., kSampling, other.back() + 1, kPhiMin) != wider);
  assert(cache.Get(l, kXMin, 2000., 9., l.back() + 1, kPhiMin) != wider);

  std::printf("bessel table cache tests passed\n");
  return 0;
}
