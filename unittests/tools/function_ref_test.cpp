#include "function_ref.h"

#include <cassert>
#include <cstdio>
#include <optional>
#include <stdexcept>
#include <type_traits>

/* The same test builds against C++26's std::function_ref when CLASS_TEST_STD_FUNCTION_REF is
   defined (g++-16 -std=c++26). That is how the claim in function_ref.h -- what compiles against
   FunctionRef compiles against std::function_ref and does the same -- is checked. The few
   checks that hold for FunctionRef alone are marked. */
#ifdef CLASS_TEST_STD_FUNCTION_REF
#include <functional>
template <class Signature>
using Ref = std::function_ref<Signature>;
#else
template <class Signature>
using Ref = FunctionRef<Signature>;
#endif

namespace {

int add(int a, int b) {
  return a + b;
}

void scale(double factor, double* value) {
  *value *= factor;
}

struct Callable {
  void operator()(int /*value*/) const {}
};

/* Never default constructed, so never null: an absent callback is an empty std::optional. */
static_assert(!std::is_default_constructible_v<Ref<void(int)>>);
/* Two words, copied freely: passing one by value costs what the old (pointer, void*) did. */
static_assert(std::is_trivially_copyable_v<Ref<void(int)>>);
/* Assigning a callable is an error in both, since nothing would keep it alive... */
static_assert(!std::is_assignable_v<Ref<void(int)>&, Callable&>);
/* ...while rebinding to another reference or to a function is allowed. */
static_assert(std::is_assignable_v<Ref<void(int)>&, const Ref<void(int)>&>);
static_assert(std::is_assignable_v<Ref<int(int, int)>&, decltype(&add)>);
/* The signature is checked where the reference is made. */
static_assert(!std::is_constructible_v<Ref<void(int)>, decltype(&scale)>);

#ifndef CLASS_TEST_STD_FUNCTION_REF
/* FunctionRef only: a temporary cannot be bound, so a callback stored in an options struct
   cannot dangle. std::function_ref accepts both. */
static_assert(!std::is_constructible_v<Ref<void(int)>, Callable&&>);
static_assert(!std::is_assignable_v<std::optional<Ref<void(int)>>&, Callable&&>);
static_assert(std::is_assignable_v<std::optional<Ref<void(int)>>&, Callable&>);
static_assert(sizeof(Ref<void(int)>) == 2 * sizeof(void*));
#endif

void test_lambda_with_state() {
  int calls     = 0;
  auto counting = [&calls](int a, int b) {
    ++calls;
    return a * b;
  };
  const Ref<int(int, int)> ref = counting;
  assert(ref(6, 7) == 42);
  assert(ref(2, 3) == 6);
  assert(calls == 2);
}

void test_const_lambda() {
  const int offset   = 10;
  const auto shifted = [offset](int a) { return a + offset; };
  Ref<int(int)> ref  = shifted;
  assert(ref(5) == 15);
}

void test_free_function() {
  const Ref<int(int, int)> by_name    = add;
  const Ref<int(int, int)> by_pointer = &add;
  assert(by_name(2, 3) == 5);
  assert(by_pointer(4, 5) == 9);

  double value                         = 2.;
  const Ref<void(double, double*)> ref = scale;
  ref(3., &value);
  assert(value == 6.);
}

/* A void signature discards whatever the callable returns. */
void test_void_signature_discards_result() {
  int calls     = 0;
  auto counting = [&calls](int a) {
    ++calls;
    return a;
  };
  const Ref<void(int)> ref = counting;
  ref(1);
  assert(calls == 1);
}

/* A copy refers to the same callable; it does not copy it. */
void test_copy_refers_to_the_same_callable() {
  int calls              = 0;
  auto counting          = [&calls]() { ++calls; };
  const Ref<void()> ref  = counting;
  const Ref<void()> copy = ref;
  ref();
  copy();
  assert(calls == 2);
}

void test_rebinding() {
  auto twice             = [](int a, int b) { return 2 * (a + b); };
  Ref<int(int, int)> ref = add;
  assert(ref(1, 2) == 3);
  ref = Ref<int(int, int)>(twice);
  assert(ref(1, 2) == 6);
  ref = &add;
  assert(ref(1, 2) == 3);
}

/* Reference parameters reach the callable as references. */
void test_reference_arguments() {
  auto increment            = [](int& value) { ++value; };
  const Ref<void(int&)> ref = increment;
  int value                 = 1;
  ref(value);
  assert(value == 2);
}

/* class_test throws from inside the derivatives; it has to reach the caller. */
void test_exceptions_propagate() {
  auto failing                = [](double) { throw std::runtime_error("from the callable"); };
  const Ref<void(double)> ref = failing;
  bool caught                 = false;
  try {
    ref(1.);
  }
  catch (const std::runtime_error&) {
    caught = true;
  }
  assert(caught);
}

/* The optional-callback pattern of EvolverOptions. */
void test_optional_callback() {
  std::optional<Ref<void(double, double*)>> callback;
  assert(!callback.has_value());

  auto halve = [](double, double* value) { *value *= 0.5; };
  callback   = halve;
  assert(callback.has_value());
  double value = 4.;
  (*callback)(0., &value);
  assert(value == 2.);

  callback = scale;
  (*callback)(3., &value);
  assert(value == 6.);

  callback.reset();
  assert(!callback.has_value());
}

}  // namespace

int main() {
  test_lambda_with_state();
  test_const_lambda();
  test_free_function();
  test_void_signature_discards_result();
  test_copy_refers_to_the_same_callable();
  test_rebinding();
  test_reference_arguments();
  test_exceptions_propagate();
  test_optional_callback();
  std::printf("function_ref tests passed\n");
  return 0;
}
