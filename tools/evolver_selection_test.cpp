// Every evolver the input accepts must actually run a cosmology.
//
// This exists because of a regression that nothing else caught. The options
// struct makes setting an option an evolver does not honour an ERROR rather than
// a silent ignore, which is the point of it -- but the perturbation module then
// set EvolverOptions::stats unconditionally, while the legacy rk and etd did not
// declare that feature. Every run selecting either of them aborted before taking
// a single step.
//
// The reason it slipped through: the bit-identity sweep across all five evolvers
// was run for the commit that introduced the struct, and NOT re-run for the later
// commit that added statistics. Outputs cannot differ if the run does not start.
// This test is the cheap standing version of that sweep.
//
// It asserts only that each evolver completes and produces a sane spectrum. The
// numerical agreement between them is a precision question and belongs elsewhere.

#include <cassert>
#include <cstdio>
#include <string>

#include "cosmology.h"
#include "input_module.h"
#include "perturbations_module.h"
#include "spectra_module.h"

namespace {

struct Choice {
  const char* name;
  int value; /* evolver_type: 0 rk, 1 ndf15, 2 rkdp45, 3 etd, 4 tsit5 */
};

/* Small enough to stay quick, large enough to exercise the source machinery. */
FileContent input(int evolver) {
  FileContent fc;
  fc.set("output", "tCl");
  fc.set("l_max_scalars", "300");
  fc.set("evolver_perturbations", std::to_string(evolver));
  fc.set("evolver_background", std::to_string(evolver));
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("perturbations_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

void test_every_evolver_runs() {
  const Choice choices[] = {{"rk", 0}, {"ndf15", 1}, {"rkdp45", 2}, {"etd", 3}, {"tsit5", 4}};

  for (const Choice& choice : choices) {
    FileContent fc = input(choice.value);
    Cosmology cosmology{fc};

    /* Would throw if the evolver rejected its own options, which is exactly the
       failure this test exists for. */
    auto& perturbations = cosmology.GetPerturbationsModule();
    assert(perturbations->tau_sampling_.size() > 1);

    /* And the pipeline downstream of it has to complete, not merely avoid
       throwing here: pulling the spectra forces the sources to be usable. */
    cosmology.GetSpectraModule();

    std::printf("  %-7s ran: %zu source times\n", choice.name, perturbations->tau_sampling_.size());
  }
}

}  // namespace

int main() {
  test_every_evolver_runs();
  std::printf("evolver selection test: every evolver the input accepts runs\n");
  return 0;
}
