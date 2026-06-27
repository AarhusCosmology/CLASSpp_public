#include <cassert>
#include <cstdio>
#include <memory>

#include "background.h"
#include "cdm.h"
#include "perturbations.h"  // possible_gauges

int main() {
  background pba{};
  pba.H0         = 1e-4;
  const int sync = static_cast<int>(possible_gauges::synchronous);

  // Uncoupled CDM in synchronous gauge: no theta variable.
  {
    CDMSpecies cdm(pba, 0.12, /*coupled=*/false);
    auto layout  = cdm.CreatePerturbLayout();
    int index_pt = 0;
    cdm.RegisterPerturbationIndices(*layout, nullptr, nullptr, index_pt, nullptr, sync);
    const auto& l = static_cast<const CDMSpecies::PerturbLayout&>(*layout);
    assert(l.idx_delta == 0);
    assert(l.idx_theta == -1);  // theta_cdm = 0 by gauge choice
    assert(index_pt == 1);
  }

  // Coupled CDM in synchronous gauge: theta is registered.
  {
    CDMSpecies cdm(pba, 0.12, /*coupled=*/true);
    auto layout  = cdm.CreatePerturbLayout();
    int index_pt = 0;
    cdm.RegisterPerturbationIndices(*layout, nullptr, nullptr, index_pt, nullptr, sync);
    const auto& l = static_cast<const CDMSpecies::PerturbLayout&>(*layout);
    assert(l.idx_delta == 0);
    assert(l.idx_theta == 1);
    assert(index_pt == 2);
  }

  std::printf("cdm coupled-velocity tests passed\n");
  return 0;
}
