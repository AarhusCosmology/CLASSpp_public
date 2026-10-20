#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "parser.h"                         // FileContent
#include "species/species_build_context.h"  // SpeciesBuildContext
#include "species/type3_species.h"          // Type3Species::CreateAll

// Returns true iff CreateAll throws the gauge guard specifically (its message).
static bool RejectsAsNewtonian(const std::string& gauge_value) {
  FileContent fc;
  fc.set("scf_veta", "-0.5");  // coupling active, beta < 1/2
  fc.set("gauge", gauge_value);
  SpeciesBuildContext ctx{};  // only pfc is read before the gauge guard fires
  ctx.pfc = &fc;
  try {
    Type3Species::CreateAll(ctx);
    return false;  // no throw at all
  }
  catch (const std::exception& e) {
    // Distinguish the gauge guard from the later "no scalar field configured" throw
    // that a slipped-through gauge would reach (no Omega_scf set here).
    return std::string(e.what()).find("synchronous gauge only") != std::string::npos;
  }
}

int main() {
  // Every spelling the parser resolves to Newtonian must be rejected by the guard.
  assert(RejectsAsNewtonian("newtonian"));
  assert(RejectsAsNewtonian("Newtonian"));
  assert(RejectsAsNewtonian("new"));   // abbreviated; the find("newton") bug let this slip
  assert(RejectsAsNewtonian("newt"));  // ditto

  // Synchronous gauge must NOT trip the guard (it proceeds past the gauge check;
  // RejectsAsNewtonian is false because the throw it then hits is the unrelated
  // "no scalar field configured" one, not the "synchronous gauge only" guard).
  assert(!RejectsAsNewtonian("synchronous"));
  assert(!RejectsAsNewtonian("sync"));

  std::printf("type3 gauge guard test passed\n");
  return 0;
}
