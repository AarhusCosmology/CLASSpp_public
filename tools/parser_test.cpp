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

  auto m = input.get<double>("m");
  assert(m && *m == 0.06);
  assert(fc.was_read("nu1.m"));

  assert(!input.get<double>("deg").has_value());
  assert(!fc.was_read("nu1.deg"));

  auto t = input.get<std::string>("type");
  assert(t && *t == "ncdm_standard");
}

static void test_species_input_required_throws() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  SpeciesInput input(&fc, "nu1");
  bool threw = false;
  try {
    (void) input.require<double>("m");
  }
  catch (const std::invalid_argument& e) {
    threw           = true;
    std::string msg = e.what();
    assert(msg.find("nu1") != std::string::npos);
    assert(msg.find("m") != std::string::npos);
  }
  assert(threw);
}

static void test_get_typed() {
  FileContent fc;
  fc.set("an_int", "42");
  fc.set("a_double", "0.067");
  fc.set("a_string", "hello");
  fc.set("a_list", "1.0,2.5,3");

  assert(fc.get<int>("an_int") == 42);
  assert(fc.get<double>("a_double") == 0.067);
  assert(fc.get<std::string>("a_string").value() == "hello");
  auto l = fc.get<std::vector<double>>("a_list");
  assert(l && l->size() == 3 && (*l)[1] == 2.5);

  fc.set("int_list", "10,20,30");
  auto il = fc.get<std::vector<int>>("int_list");
  assert(il && il->size() == 3 && (*il)[2] == 30);

  fc.set("str_list", "a,b,c");
  auto sl = fc.get<std::vector<std::string>>("str_list");
  assert(sl && sl->size() == 3 && (*sl)[0] == "a");

  assert(fc.get_or("missing_int", 42) == 42);  // get_or<int> template path

  // Absent -> nullopt, and absence is not marked read.
  assert(!fc.get<int>("missing").has_value());
  assert(!fc.was_read("missing"));

  // A successful get marks the key read (parity with the old read_*).
  assert(fc.was_read("an_int"));

  // get_or: present returns value, absent returns fallback.
  assert(fc.get_or("a_double", 1.0) == 0.067);
  assert(fc.get_or("missing", 1.0) == 1.0);
  assert(fc.get_or("a_string", "fallback") == "hello");  // const char* overload
  assert(fc.get_or("missing_str", "fallback") == "fallback");
}

static void test_get_parse_error_throws() {
  FileContent fc;
  fc.set("bad_int", "not_a_number");
  bool threw = false;
  try {
    (void) fc.get<int>("bad_int");
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);
}

int main() {
  test_instances_with_basic();
  test_species_input_prefixing();
  test_species_input_required_throws();
  test_get_typed();
  test_get_parse_error_throws();
  return 0;
}
