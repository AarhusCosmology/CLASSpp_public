#include "common.h"

#include <cstdarg>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

std::string FormatThrowMessage(
    const char* func, int line, const char* prefix, const char* fmt, va_list args) {
  va_list args_copy;
  va_copy(args_copy, args);
  const int n = std::vsnprintf(nullptr, 0, fmt, args_copy);
  va_end(args_copy);

  std::string body;
  if (n > 0) {
    body.resize(static_cast<std::size_t>(n));
    std::vsnprintf(body.data(), static_cast<std::size_t>(n) + 1, fmt, args);
  }

  // Preserves the historical "<func>(L:<line>) :<message>" format.
  return std::string(func) + "(L:" + std::to_string(line) + ") :" + prefix + body;
}

}  // namespace

[[noreturn]] void ThrowFormatted(
    const char* func, int line, const char* prefix, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  const std::string message = FormatThrowMessage(func, line, prefix, fmt, args);
  va_end(args);
  throw std::runtime_error(message);
}

[[noreturn]] void ThrowFormattedSevere(
    const char* func, int line, const char* prefix, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  const std::string message = FormatThrowMessage(func, line, prefix, fmt, args);
  va_end(args);
  throw std::invalid_argument(message);
}

int get_number_of_titles(const std::string& titlestring) {
  int number_of_titles = 0;
  for (char c : titlestring)
    if (c == '\t')
      number_of_titles++;
  return number_of_titles;
}
