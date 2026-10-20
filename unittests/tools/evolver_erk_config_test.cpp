// The explicit-RK controller settings live in a process-wide slot, because the
// evolver signature shared with ndf15/rk/etd has no room for them. That is only
// safe if the slot is written by the module that is about to integrate, reading
// its OWN precision struct.
//
// It used to be written at parse time instead, which leaves the value that reaches
// the integrator dependent on how many other cosmologies were parsed in between.
// This pins the property that it does not.
//
// Honest caveat: this test was written to reproduce that hazard and could NOT be
// made to fail against the old arrangement -- reinstating parse-time
// configuration and removing the module-site call still produced the right
// answer here, for a reason not run to ground. So read it as pinning a property
// that must hold, not as the reproduction of a bug that was observed.

#include <cassert>
#include <cstdio>

#include "cosmology.h"
#include "evolver_erk.h"
#include "input_module.h"
#include "perturbations_module.h"

namespace {

// Small enough to run in about a second, big enough that the controller choice
// changes the step count.
FileContent input(int controller) {
  FileContent fc;
  fc.set("output", "tCl");
  fc.set("l_max_scalars", "300");
  fc.set("evolver_perturbations", "2");  // rkdp45: the controller only bites here
  fc.set("rk_controller", std::to_string(controller));
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("perturbations_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

long long steps_for(FileContent fc) {
  Cosmology cosmology{fc};
  /* Step counts come from the module, which sums them over every wavenumber and
     thread. There is no enable/reset dance because there is no shared counter to
     reset: each evolver call fills its own struct. */
  return cosmology.GetPerturbationsModule()->evolver_stats_.steps_accepted;
}

}  // namespace

int main() {
  const long long alone_pi     = steps_for(input(2));
  const long long alone_legacy = steps_for(input(0));
  printf("  PI alone: %lld accepted steps; legacy alone: %lld\n", alone_pi, alone_legacy);
  // If these agreed the test could not detect anything, so check the premise.
  assert(alone_pi != alone_legacy);

  // The ordering that used to break. Constructing a Cosmology parses nothing --
  // it is lazy all the way down -- so the input modules have to be pulled
  // explicitly, which is exactly what any caller inspecting derived parameters
  // before computing does.
  FileContent fc_pi = input(2);
  Cosmology first{fc_pi};
  FileContent fc_legacy = input(0);
  Cosmology second{fc_legacy};
  first.GetInputModule();
  second.GetInputModule();  // the hazard: the LAST parse used to win for everyone

  const long long interleaved = first.GetPerturbationsModule()->evolver_stats_.steps_accepted;

  printf("  PI with a legacy cosmology parsed in between: %lld accepted steps\n", interleaved);
  assert(interleaved == alone_pi);

  /* This hazard is now structurally impossible rather than merely absent: the
     controller travels in EvolverOptions, per call, so there is no process-wide
     setting for a second parse to overwrite. The test is kept because it checks
     the surviving requirement -- that each cosmology's controller reaches its own
     evolver calls -- and because it would catch a reintroduction. */
  printf("evolver_erk_config_test: the controller follows the cosmology, not the last parse\n");
  return 0;
}
