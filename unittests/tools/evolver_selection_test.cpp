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

/* ── automatic selection ────────────────────────────────────────────────────
   InputModule::SelectPerturbationEvolver() picks the evolver when the input did
   not. The choice is three-way and the interesting case is that a species can be
   UNSAFE under a plain explicit method and still prefer the exponential one,
   when its stiffness is the Jacobian diagonal etd integrates exactly. These
   assert the resolution rather than the performance. */

FileContent auto_select_input() {
  FileContent fc;
  fc.set("output", "tCl");
  fc.set("l_max_scalars", "300");
  fc.set("H0", "67.32117");
  fc.set("omega_b", "0.02238280");
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("perturbations_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

void add_standard_neutrino(FileContent& fc) {
  fc.set("nu1.type", "ncdm_standard");
  fc.set("nu1.m", "0.06");
}

void add_self_interacting_neutrino(FileContent& fc, const char* log10G_eff) {
  fc.set("nusi.type", "ncdm_self_interacting");
  fc.set("nusi.m", "0.02");
  fc.set("nusi.deg", "2.0");
  fc.set("nusi.log10G_eff", log10G_eff);
}

/* Constructing InputModule is enough: SelectPerturbationEvolver() runs in its
   constructor, after the species are built. */
evolver_type selected(FileContent fc) {
  InputModule input{fc};
  return input.precision_.evolver_perturbations;
}

void test_automatic_selection() {
  /* Ordinary content: every species opts in, so the explicit pair is taken. */
  {
    FileContent fc = auto_select_input();
    fc.set("omega_cdm", "0.1201075");
    fc.set("N_ur", "2.0308");
    add_standard_neutrino(fc);
    assert(selected(std::move(fc)) == evolver_type::rkdp45);
  }

  /* Self-interacting neutrinos: the collision term is a pure diagonal
     relaxation, so etd is preferred even though rkdp45 would not complete. */
  {
    FileContent fc = auto_select_input();
    fc.set("omega_cdm", "0.1201075");
    fc.set("N_ur", "0.0");
    add_standard_neutrino(fc);
    add_self_interacting_neutrino(fc, "-1.35");
    assert(selected(std::move(fc)) == evolver_type::etd);
  }

  /* A species that is neither explicit-safe nor diagonal (dcdm_dr) VETOES etd
     even when something else asks for it: the exponential step is chosen for the
     whole system, so one badly-conditioned sector is enough to fall back. */
  {
    FileContent fc = auto_select_input();
    fc.set("omega_cdm", "0.09");
    fc.set("Omega_ini_dcdm", "0.05");
    fc.set("Gamma_dcdm", "100.0");
    fc.set("N_ur", "0.0");
    add_standard_neutrino(fc);
    add_self_interacting_neutrino(fc, "-1.35");
    assert(selected(std::move(fc)) == evolver_type::ndf15);
  }

  /* An explicit choice in the input always wins over all of the above. */
  {
    FileContent fc = auto_select_input();
    fc.set("omega_cdm", "0.1201075");
    fc.set("N_ur", "0.0");
    add_standard_neutrino(fc);
    add_self_interacting_neutrino(fc, "-1.35");
    fc.set("evolver_perturbations", "1");
    assert(selected(std::move(fc)) == evolver_type::ndf15);
  }
}

}  // namespace

int main() {
  test_every_evolver_runs();
  test_automatic_selection();
  std::printf("evolver selection test: every evolver runs, and automatic selection resolves\n");
  return 0;
}
