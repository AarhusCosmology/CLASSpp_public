#include "parser.h"

#include <cassert>

static void test_instances_with_basic() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  fc.set("nu1.m", "0.06");
  fc.set("nu2.type", "ncdm_self_interacting");
  fc.set("nu2.m", "0.07");
  fc.set("nu3.type", "ncdm_standard");

  auto names = fc.instances_with("type", "ncdm_standard");
  assert(names.size() == 2);
  assert(names[0] == "nu1");
  assert(names[1] == "nu3");

  auto none = fc.instances_with("type", "unknown");
  assert(none.empty());

  auto si = fc.instances_with("type", "ncdm_self_interacting");
  assert(si.size() == 1 && si[0] == "nu2");

  // Must NOT mark entries as read.
  assert(!fc.was_read("nu1.type"));
  assert(!fc.was_read("nu3.type"));
}

#include <stdexcept>
#include <string>

#include "species/species_input.h"

static void test_species_input_prefixing() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  fc.set("nu1.m", "0.06");
  fc.set("nu1.T", "0.71611");

  SpeciesInput input(&fc, "nu1");
  assert(input.instance_name() == "nu1");

  double m = 0.;
  assert(input.read_double("m", m));
  assert(m == 0.06);
  assert(fc.was_read("nu1.m"));

  double missing = 0.;
  assert(!input.read_double("deg", missing));
  assert(!fc.was_read("nu1.deg"));

  std::string t;
  assert(input.read_string("type", t));
  assert(t == "ncdm_standard");
}

static void test_species_input_required_throws() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  SpeciesInput input(&fc, "nu1");
  bool threw = false;
  try {
    (void) input.required_double("m");
  }
  catch (const std::invalid_argument& e) {
    threw           = true;
    std::string msg = e.what();
    assert(msg.find("nu1") != std::string::npos);
    assert(msg.find("m") != std::string::npos);
  }
  assert(threw);
}

int main() {
  test_instances_with_basic();
  test_species_input_prefixing();
  test_species_input_required_throws();
  return 0;
}
