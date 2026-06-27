#include <cstdio>
#include <memory>
#include <vector>

#include "background.h"
#include "cdm.h"
#include "composite_species.h"
#include "scalar_field.h"
#include "type3_species.h"

int main() {
  background pba{};
  pba.H0                           = 1e-4;
  const std::vector<double> params = {1.22, 0.0, 0.0, 0.0};
  auto scf = std::make_unique<ScalarFieldSpecies>(pba,
                                                  /*omega0=*/0.7,
                                                  params,
                                                  /*tuning=*/0,
                                                  /*attractor=*/true,
                                                  /*phi_ini=*/1.,
                                                  /*phi_prime_ini=*/1.,
                                                  DefaultScalarFieldPotential(),
                                                  /*beta=*/-2.0);
  Type3Species t3(pba, /*omega0_cdm=*/0.12, std::move(scf));

  auto layout = t3.CreatePerturbLayout();
  auto& comp  = static_cast<CompositeSpecies::PerturbLayout&>(*layout);

  // Generic CreatePerturbLayout builds one sub-layout per child, in children_ order.
  // Direct references (not in assert) so NDEBUG builds still require these to compile.
  const std::size_t n = comp.child_layouts.size();
  if (n != 2) {
    std::fprintf(stderr, "FAIL: expected 2 child_layouts, got %zu\n", n);
    return 1;
  }
  if (dynamic_cast<CDMSpecies::PerturbLayout*>(comp.child_layouts[0].get()) == nullptr) {
    std::fprintf(stderr, "FAIL: child_layouts[0] is not a CDMSpecies::PerturbLayout\n");
    return 1;
  }
  if (dynamic_cast<ScalarFieldSpecies::PerturbLayout*>(comp.child_layouts[1].get()) == nullptr) {
    std::fprintf(stderr, "FAIL: child_layouts[1] is not a ScalarFieldSpecies::PerturbLayout\n");
    return 1;
  }

  std::printf("composite layout test passed\n");
  return 0;
}
