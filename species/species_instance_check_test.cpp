// A dot-syntax species instance that no factory builds is an input error (#430).
//
// `nuM.type = ncdm` (the type is `ncdm_standard`) used to build nothing and say
// nothing: the run completed and silently computed a cosmology without the
// massive neutrino. RejectUnbuiltSpeciesInstances now rejects it, together with
// a type that exists but cannot be declared by dot syntax, an illegal instance
// name, and instance fields with no `.type` line at all.
//
// The unit cases drive the check on a FileContent alone; "read" stands in for a
// factory having built the instance. The InputModule cases run the real
// factories, so they also pin the contract in all_species.h that a factory
// reads N.type for exactly the instances it builds.
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "input_module.h"
#include "parser.h"
#include "species/species_input.h"

namespace {

/* The message of the std::invalid_argument the check throws, or "" if it passed. */
std::string CheckMessage(const FileContent& fc) {
  try {
    RejectUnbuiltSpeciesInstances(fc);
  }
  catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

bool Has(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

FileContent Quiet() {
  FileContent fc;
  fc.set("output", "");
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

/* The message of the std::invalid_argument InputModule throws, or "" if it built. */
std::string BuildMessage(FileContent fc) {
  try {
    InputModule input(fc);
  }
  catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

void test_unit() {
  // ── The #430 reproducer: an unknown type name is rejected, and the message
  //    lists the species types so the fix is on screen ──
  {
    FileContent fc;
    fc.set("nuM.type", "ncdm");
    fc.set("nuM.m", "0.06");
    const std::string msg = CheckMessage(fc);
    assert(Has(msg, "'nuM.type = ncdm'"));
    assert(Has(msg, "'ncdm' is not a species type"));
    assert(Has(msg, "ncdm_standard"));
    // One problem per instance: nuM.m belongs to a declared instance, so it is
    // not also reported as an orphan.
    assert(!Has(msg, "'nuM.m'"));
  }

  // ── A consumed type passes (read = a factory built it) ──
  {
    FileContent fc;
    fc.set("nuM.type", "ncdm_standard");
    fc.set("nuM.m", "0.06");
    (void) fc.get<std::string>("nuM.type");
    assert(CheckMessage(fc).empty());
  }

  // ── Unread fields of a BUILT instance are out of scope for this check ──
  {
    FileContent fc;
    fc.set("nuM.type", "ncdm_standard");
    fc.set("nuM.mass", "0.06");
    (void) fc.get<std::string>("nuM.type");
    assert(CheckMessage(fc).empty());
  }

  // ── A real species type that no factory reads from N.type ──
  {
    FileContent fc;
    fc.set("x.type", "scalar_field");
    const std::string msg = CheckMessage(fc);
    assert(Has(msg, "species type 'scalar_field' cannot be declared"));
    assert(Has(msg, "legacy input keys"));
  }

  // ── An illegal instance name is rejected even if something read the key ──
  {
    FileContent fc;
    fc.set("nu-M.type", "ncdm_standard");
    (void) fc.get<std::string>("nu-M.type");
    assert(Has(CheckMessage(fc), "'nu-M' is not a legal instance name"));
  }

  // ── Fields with no N.type line: a misspelled or forgotten type line.
  //    Keys are grouped per instance, in input order ──
  {
    FileContent fc;
    fc.set("nuM.typ", "ncdm_standard");
    fc.set("nuM.m", "0.06");
    const std::string msg = CheckMessage(fc);
    assert(Has(msg, "'nuM.typ', 'nuM.m' set, but there is no 'nuM.type'"));
  }

  // ── The orphan test must not depend on key order: N.type after its fields ──
  {
    FileContent fc;
    fc.set("nuM.m", "0.06");
    fc.set("nuM.type", "ncdm_standard");
    (void) fc.get<std::string>("nuM.type");
    assert(CheckMessage(fc).empty());
  }

  // ── Every offending instance is named in one message ──
  {
    FileContent fc;
    fc.set("a.type", "ncdm");
    fc.set("b.m", "1");
    const std::string msg = CheckMessage(fc);
    assert(Has(msg, "'a.type = ncdm'"));
    assert(Has(msg, "no 'b.type'"));
  }

  // ── No dot-syntax keys at all: nothing to check ──
  {
    FileContent fc;
    fc.set("N_ur", "3.046");
    fc.set("write parameters", "no");
    assert(CheckMessage(fc).empty());
  }
}

void test_input_module() {
  // ── The #430 input, exactly as filed, now fails instead of running ──
  {
    FileContent fc = Quiet();
    fc.set("nuM.type", "ncdm");
    fc.set("nuM.m", "0.06");
    fc.set("nuM.T", "0.71611");
    fc.set("nuM.deg", "1.0");
    assert(Has(BuildMessage(fc), "'ncdm' is not a species type"));
  }

  // ── ...and with the right type name it builds the species ──
  {
    FileContent fc = Quiet();
    fc.set("nuM.type", "ncdm_standard");
    fc.set("nuM.m", "0.06");
    fc.set("nuM.T", "0.71611");
    fc.set("nuM.deg", "1.0");
    InputModule input(fc);
    assert(input.all_species_.count("nuM") == 1);
  }

  // ── A single-instance type is consumed by the dot-syntax translation ──
  {
    FileContent fc = Quiet();
    fc.set("g.type", "photons");
    fc.set("g.T_cmb", "2.7255");
    assert(BuildMessage(fc).empty());
  }

  // ── The legacy N_ncdm path synthesises ncdm__1.* and must pass the check ──
  {
    FileContent fc = Quiet();
    fc.set("N_ncdm", "1");
    fc.set("m_ncdm", "0.06");
    assert(BuildMessage(fc).empty());
  }

  // ── A species type that no factory reads from N.type ──
  {
    FileContent fc = Quiet();
    fc.set("x.type", "scalar_field");
    assert(Has(BuildMessage(fc), "species type 'scalar_field' cannot be declared"));
  }

  // ── The type line misspelled: previously the species was silently absent ──
  {
    FileContent fc = Quiet();
    fc.set("nuM.typ", "ncdm_standard");
    fc.set("nuM.m", "0.06");
    assert(Has(BuildMessage(fc), "there is no 'nuM.type'"));
  }
}

}  // namespace

int main() {
  test_unit();
  test_input_module();
  std::printf("species instance check: all tests passed\n");
  return 0;
}
