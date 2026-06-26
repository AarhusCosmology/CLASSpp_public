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

static void test_dot_baryons_translation() {
  FileContent fc;
  fc.set("atom.type", "baryons");
  fc.set("atom.Omega", "0.05");
  TranslateSingleInstanceDotSyntax(&fc);
  auto ob = fc.get<std::string>("Omega_b");
  assert(ob && *ob == "0.05");
  assert(fc.was_read("atom.type"));
  assert(fc.was_read("atom.Omega"));
}

static void test_dot_name_is_ignored() {
  // Any instance name maps to the same legacy key.
  FileContent fc;
  fc.set("whatever.type", "cdm");
  fc.set("whatever.Omega", "0.25");
  TranslateSingleInstanceDotSyntax(&fc);
  auto oc = fc.get<std::string>("Omega_cdm");
  assert(oc && *oc == "0.25");
}

static void test_dot_fluid_full_set() {
  FileContent fc;
  fc.set("de.type", "fluid");
  fc.set("de.Omega", "0.7");
  fc.set("de.w0", "-0.9");
  fc.set("de.wa", "0.1");
  fc.set("de.cs2", "1.0");
  fc.set("de.use_ppf", "no");
  fc.set("de.equation_of_state", "CLP");
  fc.set("de.Omega_EDE", "0.01");
  fc.set("de.c_gamma_over_c", "0.4");
  TranslateSingleInstanceDotSyntax(&fc);
  assert(*fc.get<std::string>("Omega_fld") == "0.7");
  assert(*fc.get<std::string>("w0_fld") == "-0.9");
  assert(*fc.get<std::string>("wa_fld") == "0.1");
  assert(*fc.get<std::string>("cs2_fld") == "1.0");
  assert(*fc.get<std::string>("use_ppf") == "no");
  assert(*fc.get<std::string>("fluid_equation_of_state") == "CLP");
  assert(*fc.get<std::string>("Omega_EDE") == "0.01");
  assert(*fc.get<std::string>("c_gamma_over_c_fld") == "0.4");
}

static void test_dot_ur_forms() {
  FileContent fc;
  fc.set("nu.type", "ur");
  fc.set("nu.N", "3.044");
  TranslateSingleInstanceDotSyntax(&fc);
  assert(*fc.get<std::string>("N_ur") == "3.044");
}

static void test_dot_duplicate_instance_throws() {
  FileContent fc;
  fc.set("a.type", "baryons");
  fc.set("a.Omega", "0.05");
  fc.set("b.type", "baryons");
  fc.set("b.Omega", "0.04");
  bool threw = false;
  try {
    TranslateSingleInstanceDotSyntax(&fc);
  }
  catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

static void test_dot_conflict_with_legacy_throws() {
  FileContent fc;
  fc.set("Omega_b", "0.04");
  fc.set("x.type", "baryons");
  fc.set("x.Omega", "0.05");
  bool threw = false;
  try {
    TranslateSingleInstanceDotSyntax(&fc);
  }
  catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

static void test_dot_identical_legacy_ok() {
  FileContent fc;
  fc.set("Omega_b", "0.05");
  fc.set("x.type", "baryons");
  fc.set("x.Omega", "0.05");
  TranslateSingleInstanceDotSyntax(&fc);  // must not throw
  assert(*fc.get<std::string>("Omega_b") == "0.05");
}

static void test_dot_unknown_field_left_unread() {
  FileContent fc;
  fc.set("x.type", "baryons");
  fc.set("x.Omeega", "0.05");  // typo: not in the table
  TranslateSingleInstanceDotSyntax(&fc);
  assert(!fc.was_read("x.Omeega"));  // stays unread -> warned later
}

int main() {
  test_instances_with_basic();
  test_species_input_prefixing();
  test_species_input_required_throws();
  test_get_typed();
  test_get_parse_error_throws();
  test_dot_baryons_translation();
  test_dot_name_is_ignored();
  test_dot_fluid_full_set();
  test_dot_ur_forms();
  test_dot_duplicate_instance_throws();
  test_dot_conflict_with_legacy_throws();
  test_dot_identical_legacy_ok();
  test_dot_unknown_field_left_unread();
  return 0;
}
