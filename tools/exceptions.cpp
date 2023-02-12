#include "exceptions.h"
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdarg>
#include <typeinfo>

void ThrowRuntimeErrorIf(bool condition, std::string string_for_printf, ...) {
  if (condition) {
    va_list args;
    va_start(args, string_for_printf);
    const size_t error_message_size = 2048;
    char error_message[error_message_size];
    vsnprintf(error_message, error_message_size, string_for_printf.c_str(), args);
    throw std::runtime_error(error_message);
  }
}

void ThrowRuntimeError(std::string string_for_printf, ...) {
  va_list args;
  va_start(args, string_for_printf);
  const size_t error_message_size = 2048;
  char error_message[error_message_size];
  vsnprintf(error_message, error_message_size, string_for_printf.c_str(), args);
  throw std::runtime_error(error_message);
}

std::pair<std::string, std::string> get_my_py_error_message() {
  std::pair<std::string, std::string> res;
  try {
    throw;
  } catch (const std::exception& e) {
    res.first = std::string(typeid(e).name());
    res.second = std::string(e.what());
  }
  return res;
}
