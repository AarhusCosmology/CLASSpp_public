#include "errors.h"

#include <cassert>
#include <stdexcept>
#include <string>

// class_test: runtime_error, message carries func name, condition text, %g-formatted arg.
static void test_class_test_type_and_message() {
  bool threw     = false;
  const double q = 1.5e-12;
  try {
    class_test(q > 0., "q = %g must not be positive", q);
  }
  catch (const std::invalid_argument&) {
    assert(false && "class_test must throw runtime_error, not invalid_argument");
  }
  catch (const std::runtime_error& e) {
    threw                 = true;
    const std::string msg = e.what();
    assert(msg.find("test_class_test_type_and_message") != std::string::npos);
    assert(msg.find("condition (q > 0.) is true") != std::string::npos);
    assert(msg.find("1.5e-12") != std::string::npos);
  }
  assert(threw);
}

// class_test_severe: identical format, but invalid_argument.
static void test_class_test_severe_type_and_message() {
  bool threw = false;
  try {
    class_test_severe(1 == 1, "keys '%s' and '%s' are mutually exclusive", "a", "b");
  }
  catch (const std::invalid_argument& e) {
    threw                 = true;
    const std::string msg = e.what();
    assert(msg.find("condition (1 == 1) is true") != std::string::npos);
    assert(msg.find("keys 'a' and 'b' are mutually exclusive") != std::string::npos);
  }
  assert(threw);
}

// invalid_argument IS-A logic_error, NOT runtime_error: severe must not be
// catchable as runtime_error (that would re-blur the severity channel).
static void test_severe_is_not_runtime_error() {
  bool caught_as_runtime = false;
  try {
    class_stop_severe("structural failure");
  }
  catch (const std::runtime_error&) {
    caught_as_runtime = true;
  }
  catch (const std::exception&) {
  }
  assert(!caught_as_runtime);
}

static void test_class_stop_and_stop_severe() {
  bool threw = false;
  try {
    class_stop("gave up after %d iterations", 20);
  }
  catch (const std::runtime_error& e) {
    threw = true;
    assert(std::string(e.what()).find("error; gave up after 20 iterations") != std::string::npos);
  }
  assert(threw);

  threw = false;
  try {
    class_stop_severe("incomprehensible input '%s'", "spam");
  }
  catch (const std::invalid_argument& e) {
    threw = true;
    assert(std::string(e.what()).find("error; incomprehensible input 'spam'") != std::string::npos);
  }
  assert(threw);
}

static void test_false_condition_does_not_throw() {
  class_test(false, "unreachable");
  class_test_severe(false, "unreachable");
}

int main() {
  test_class_test_type_and_message();
  test_class_test_severe_type_and_message();
  test_severe_is_not_runtime_error();
  test_class_stop_and_stop_severe();
  test_false_condition_does_not_throw();
  return 0;
}
