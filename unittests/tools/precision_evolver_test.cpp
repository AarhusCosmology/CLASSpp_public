// Evolver input resolution: the blanket `evolver` key XOR the per-module family.
//
// The three callers (background_module.cpp, thermodynamics_module.cpp,
// perturbations_module.cpp) do not want the same integrator for every model --
// the decaying-NCDM inverse-decay sector has a DENSE background Jacobian (where
// ndf15 costs ~84x rkdp45 for the same answer) and a SPARSE, stiffness-limited
// perturbation Jacobian. These tests pin the resolution rules so a later reader
// does not have to infer them from precision::parse.

#include <cassert>
#include <stdexcept>

#include "errors.h"
#include "parser.h"
#include "precision.h"

static precision parsed(const FileContent& fc) {
  precision ppr;
  ppr.parse(fc);
  return ppr;
}

// Every module defaults to ndf15 when nothing is specified.
static void test_default_is_ndf15() {
  FileContent fc;
  const precision p = parsed(fc);
  assert(p.evolver_background == evolver_type::ndf15);
  assert(p.evolver_thermodynamics == evolver_type::ndf15);
  assert(p.evolver_perturbations == evolver_type::ndf15);
}

// The blanket key still sets all three -- existing .ini/.pre files keep working.
static void test_blanket_sets_all() {
  FileContent fc;
  fc.set("evolver", "2");  // rkdp45
  const precision p = parsed(fc);
  assert(p.evolver_background == evolver_type::rkdp45);
  assert(p.evolver_thermodynamics == evolver_type::rkdp45);
  assert(p.evolver_perturbations == evolver_type::rkdp45);
}

// A per-module key sets ONLY its own module; the others keep the default. This
// is the property that matters: adding a fourth caller later must not silently
// require every existing config to name it.
static void test_per_module_is_independent() {
  FileContent fc;
  fc.set("evolver_background", "2");  // rkdp45
  const precision p = parsed(fc);
  assert(p.evolver_background == evolver_type::rkdp45);
  assert(p.evolver_thermodynamics == evolver_type::ndf15);
  assert(p.evolver_perturbations == evolver_type::ndf15);

  FileContent fc2;
  fc2.set("evolver_perturbations", "0");  // rk
  const precision p2 = parsed(fc2);
  assert(p2.evolver_background == evolver_type::ndf15);
  assert(p2.evolver_thermodynamics == evolver_type::ndf15);
  assert(p2.evolver_perturbations == evolver_type::rk);
}

// Naming the whole family is equivalent to the blanket key.
static void test_full_family_matches_blanket() {
  FileContent fc;
  fc.set("evolver_background", "2");
  fc.set("evolver_thermodynamics", "2");
  fc.set("evolver_perturbations", "2");
  const precision p = parsed(fc);
  assert(p.evolver_background == evolver_type::rkdp45);
  assert(p.evolver_thermodynamics == evolver_type::rkdp45);
  assert(p.evolver_perturbations == evolver_type::rkdp45);
}

// The tsit5 value has to survive the parser and reach the right slot. The driver
// tests exercise the integrator directly and would not notice if `evolver = 4`
// were dropped on the floor between the input file and the module dispatch.
static void test_tsit5_is_parsed() {
  FileContent fc;
  fc.set("evolver", "4");
  const precision p = parsed(fc);
  assert(p.evolver_background == evolver_type::tsit5);
  assert(p.evolver_thermodynamics == evolver_type::tsit5);
  assert(p.evolver_perturbations == evolver_type::tsit5);

  FileContent fc2;
  fc2.set("evolver_perturbations", "4");
  const precision p2 = parsed(fc2);
  assert(p2.evolver_perturbations == evolver_type::tsit5);
  assert(p2.evolver_background == evolver_type::ndf15);
}

// Mixing the two spellings is REJECTED rather than resolved by precedence. With
// a precedence rule, a visible `evolver = 2` plus a stale per-module key in an
// inherited .pre file would silently run a different integrator than the file
// appears to ask for. Structural (decidable from which keys are present, not
// from any value a sampler varies) => severe, i.e. invalid_argument.
static void test_mixing_is_rejected() {
  auto rejects = [](const char* per_module_key) {
    FileContent fc;
    fc.set("evolver", "2");
    fc.set(per_module_key, "1");
    try {
      precision ppr;
      ppr.parse(fc);
    }
    catch (const std::invalid_argument&) {
      return true;
    }
    catch (...) {
      return false;
    }
    return false;
  };
  assert(rejects("evolver_background"));
  assert(rejects("evolver_thermodynamics"));
  assert(rejects("evolver_perturbations"));
}

int main() {
  test_default_is_ndf15();
  test_blanket_sets_all();
  test_per_module_is_independent();
  test_full_family_matches_blanket();
  test_tsit5_is_parsed();
  test_mixing_is_rejected();
  std::printf("precision evolver test passed\n");
  return 0;
}
