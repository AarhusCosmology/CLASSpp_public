/** @file errors.h Error-reporting (C++ exceptions) + status codes. */
#ifndef CLASS_ERRORS_H
#define CLASS_ERRORS_H

#include <cstdio>  // fopen / nullptr-compare in class_open

#define _TRUE_ 1
#define _FALSE_ 0
#define _SUCCESS_ 0
#define _FAILURE_ 1

[[noreturn]] void ThrowFormatted(
    const char* func, int line, const char* prefix, const char* fmt, ...);

/* All error macros throw std::runtime_error from the error origin. */

#define class_call(function)                                         \
  {                                                                  \
    if ((function) == _FAILURE_) {                                   \
      ThrowFormatted(__func__, __LINE__, "error in " #function, ""); \
    }                                                                \
  }

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
