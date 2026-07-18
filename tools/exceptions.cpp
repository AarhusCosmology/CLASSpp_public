#include "exceptions.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>

std::pair<std::string, std::string> get_my_py_error_message() {
  std::pair<std::string, std::string> res;
  try {
    throw;
  }
  catch (const std::exception& e) {
    res.first  = std::string(typeid(e).name());
    res.second = std::string(e.what());
  }
  return res;
}
