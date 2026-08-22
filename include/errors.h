/** @file errors.h Error-reporting (C++ exceptions). */
#ifndef CLASS_ERRORS_H
#define CLASS_ERRORS_H

#include <cstdio>  // fopen / nullptr-compare in class_open

[[noreturn]] void ThrowFormatted(
    const char* func, int line, const char* prefix, const char* fmt, ...);

[[noreturn]] void ThrowFormattedSevere(
    const char* func, int line, const char* prefix, const char* fmt, ...);

/* class_test / class_stop throw std::runtime_error from the error origin
   (computation errors). The severe variants further down throw
   std::invalid_argument; see their comment block for which to use when. */

#define class_test(condition, args, ...)                     \
  {                                                          \
    if (condition) {                                         \
      ThrowFormatted(__func__,                               \
                     __LINE__,                               \
                     "condition (" #condition ") is true; ", \
                     args,                                   \
                     ##__VA_ARGS__);                         \
    }                                                        \
  }

#define class_stop(args, ...)                                           \
  {                                                                     \
    ThrowFormatted(__func__, __LINE__, "error; ", args, ##__VA_ARGS__); \
  }

/* Severe variants: throw std::invalid_argument (classy: CosmoSevereError, the
   sampler ABORTS). Use ONLY for checks that depend exclusively on the
   STRUCTURE of the input: key presence, mutual exclusivity, unparseable
   strings, unknown enum names, list-length/count consistency, API-argument
   validation. A severe check must never depend on possibly-varying (numeric)
   parameters: range checks on cosmological values and all numerical failures
   use class_test/class_stop (std::runtime_error -> CosmoComputationError, the
   sampler rejects the point and the chain survives).

   "Possibly-varying" is the operative test, not "numeric". Precision and
   run-configuration parameters are fixed for a whole run -- a sampler never
   varies them -- so a numeric check depending only on those is a statement
   about the CONFIGURATION and belongs in the severe channel. Putting such a
   check in class_test is actively harmful: whether it fires is then
   cosmology-dependent, so the sampler reads a broken setup as zero likelihood
   and silently carves a hole out of the posterior instead of reporting it
   (#395). The rule is about who can change the operand, not its type. */

#define class_test_severe(condition, args, ...)                    \
  {                                                                \
    if (condition) {                                               \
      ThrowFormattedSevere(__func__,                               \
                           __LINE__,                               \
                           "condition (" #condition ") is true; ", \
                           args,                                   \
                           ##__VA_ARGS__);                         \
    }                                                              \
  }

#define class_stop_severe(args, ...)                                          \
  {                                                                           \
    ThrowFormattedSevere(__func__, __LINE__, "error; ", args, ##__VA_ARGS__); \
  }

#define class_open(pointer, filename, mode)                        \
  {                                                                \
    pointer = fopen(filename, mode);                               \
    if (pointer == nullptr) {                                      \
      ThrowFormatted(__func__,                                     \
                     __LINE__,                                     \
                     "",                                           \
                     "could not open %s with name %s and mode %s", \
                     #pointer,                                     \
                     filename,                                     \
                     #mode);                                       \
    }                                                              \
  }

#endif  // CLASS_ERRORS_H
