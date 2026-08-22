/* Regression test for #395: a precision-knob pair that puts the start of
   integration after tight coupling has switched off must be rejected as a
   SEVERE (configuration) error, not as an ordinary computation error.

   The distinction is the whole point. Whether perturb_vector_init actually
   trips depends on the cosmology, so an ordinary class_test would reject
   isolated points in parameter space and a sampler would read that as zero
   likelihood -- silently carving a hole out of the posterior. Both operands are
   precision parameters, fixed for a whole run, so the configuration is what is
   wrong and the run must abort. */

#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "cosmology.h"
#include "input_module.h"

namespace {

FileContent input(const char* modes,
                  const char* start_small_k,
                  const char* trigger_key,
                  const char* trigger) {
  FileContent fc;
  fc.set("output", "tCl");
  fc.set("modes", modes);
  if (start_small_k)
    fc.set("start_small_k_at_tau_c_over_tau_h", start_small_k);
  if (trigger_key)
    fc.set(trigger_key, trigger);
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

// Parsing is lazy, so the guard only runs once the input module is pulled.
void parse(FileContent fc) {
  Cosmology cosmology{fc};
  cosmology.GetInputModule();
}

void expect_severe(FileContent fc, const char* named_parameter) {
  bool threw = false;
  try {
    parse(fc);
  }
  catch (const std::invalid_argument& e) {
    threw                 = true;
    const std::string msg = e.what();
    // The message must name the parameters, or there is nothing to act on.
    assert(msg.find("start_small_k_at_tau_c_over_tau_h") != std::string::npos);
    assert(msg.find(named_parameter) != std::string::npos);
  }
  assert(threw && "misordered tight-coupling knobs must raise a SEVERE error");
}

void expect_ok(FileContent fc) {
  parse(fc);  // must not throw
}

}  // namespace

int main() {
  // Scalars: start (0.02) is later than the switch-off trigger (0.015 default).
  expect_severe(input("s", "0.02", nullptr, nullptr), "tight_coupling_trigger_tau_c_over_tau_h");

  // Equality is already too late: at the start time tau_c/tau_H has reached the
  // trigger, so the approximation is off.
  expect_severe(input("s", "0.015", nullptr, nullptr), "tight_coupling_trigger_tau_c_over_tau_h");

  // Tensors carry their own trigger and must be guarded separately.
  expect_severe(input("t", "0.02", nullptr, nullptr),
                "tight_coupling_trigger_tau_c_over_tau_h_ten");

  // A tensor-only knob must not reject a scalar-only run, and vice versa.
  expect_ok(input("s", "0.006", "tight_coupling_trigger_tau_c_over_tau_h_ten", "0.001"));
  expect_ok(input("t", "0.006", "tight_coupling_trigger_tau_c_over_tau_h", "0.001"));

  // Defaults must stay valid.
  expect_ok(input("s", nullptr, nullptr, nullptr));

  std::printf("tca ic guard tests passed\n");
  return 0;
}
