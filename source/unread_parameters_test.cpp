// A parameter that nothing reads is an input error, in the binary as in classy.
//
// `./class` used to ignore it without a word: a misspelt key ran a healthy-looking
// cosmology without the setting it named. Cosmology(FileContent&) now rejects it
// (InputModule::RejectUnreadParameters), unless the input sets
// `allow_unread_parameters = yes`.
//
// The check must see the caller's input only. Shooting and pk_eq build their runs
// from copies, and a copy can leave keys of the original unread: with legacy NCDM
// keys, every 100*theta_s trial leaves `m_ncdm` unread. The last case fails if the
// check ever moves into InputModule.
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "cosmology.h"

namespace {

FileContent Quiet() {
  FileContent fc;
  fc.set("output", "");
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

/* The message of the std::invalid_argument Cosmology throws, or "" if it built. */
std::string BuildMessage(FileContent fc) {
  try {
    Cosmology cosmology{fc};
  }
  catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

bool Has(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

void test_unread_parameter_is_rejected() {
  FileContent fc = Quiet();
  fc.set("omega_cmd", "0.12");  // omega_cdm, misspelt
  const std::string message = BuildMessage(fc);
  assert(Has(message, "did not read input parameter(s): omega_cmd"));
  assert(Has(message, "allow_unread_parameters = yes"));
  std::printf("unread parameter rejected: ok\n");
}

void test_every_unread_parameter_is_named() {
  FileContent fc = Quiet();
  fc.set("omega_cmd", "0.12");
  fc.set("N_eff", "3.044");  // read; a deprecated synonym of N_ur
  fc.set("not_a_parameter", "1");
  const std::string message = BuildMessage(fc);
  assert(Has(message, "omega_cmd"));
  assert(Has(message, "not_a_parameter"));
  assert(!Has(message, "N_eff"));
  std::printf("every unread parameter named: ok\n");
}

void test_allow_unread_parameters() {
  FileContent fc = Quiet();
  fc.set("omega_cmd", "0.12");
  fc.set("allow_unread_parameters", "yes");
  assert(BuildMessage(fc).empty());
  std::printf("allow_unread_parameters = yes runs anyway: ok\n");
}

void test_clean_input_builds() {
  assert(BuildMessage(Quiet()).empty());
  std::printf("clean input builds: ok\n");
}

void test_shooting_copies_are_not_checked() {
  FileContent fc = Quiet();
  fc.set("100*theta_s", "1.041783");
  fc.set("N_ur", "2.0328");
  fc.set("N_ncdm", "1");
  fc.set("m_ncdm", "0.06");
  Cosmology cosmology{fc};
  cosmology.GetInputModule();  // shoots theta_s; its trials leave m_ncdm unread
  std::printf("shooting copies not checked: ok\n");
}

}  // namespace

int main() {
  test_unread_parameter_is_rejected();
  test_every_unread_parameter_is_named();
  test_allow_unread_parameters();
  test_clean_input_builds();
  test_shooting_copies_are_not_checked();
  return 0;
}
