#include "common.h"

#include <cstdarg>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>

namespace {

std::string FormatThrowMessage(
    const char* func, int line, const char* prefix, const char* fmt, va_list args) {
  va_list args_copy;
  va_copy(args_copy, args);
  const int n = std::vsnprintf(nullptr, 0, fmt, args_copy);
  // va_copy initialised args_copy; the analyzer misreads va_copy on AArch64, where
  // va_list is a struct.
  va_end(args_copy);  // NOLINT(clang-analyzer-security.VAList)

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
  for (char c : titlestring) {
    if (c == '\t') {
      number_of_titles++;
    }
  }
  return number_of_titles;
}

std::string first_duplicate_title(const std::string& titlestring) {
  std::set<std::string> seen;
  size_t start = 0;
  while (start < titlestring.size()) {
    size_t end = titlestring.find('\t', start);
    if (end == std::string::npos) {
      end = titlestring.size();
    }
    std::string title = titlestring.substr(start, end - start);
    if (!title.empty() && !seen.insert(title).second) {
      return title;
    }
    start = end + 1;
  }
  return {};
}

NamedColumns SplitIntoColumns(const std::string& titlestring,
                              const double* rows,
                              int number_of_rows) {
  // The background and perturbation titles are checked where they are built; this check
  // covers every table, so a repeated title can never hide a column in classy's dicts.
  const std::string duplicate = first_duplicate_title(titlestring);
  class_test_severe(!duplicate.empty(),
                    "two columns of one output table are titled '%s'",
                    duplicate.c_str());
  // Every title ends with the delimiter, so the delimiters count the columns.
  const int number_of_titles = get_number_of_titles(titlestring);
  NamedColumns columns;
  columns.reserve(number_of_titles);
  size_t start = 0;
  for (int index_title = 0; index_title < number_of_titles; index_title++) {
    const size_t end = titlestring.find('\t', start);
    std::vector<double> column(number_of_rows);
    for (int index_row = 0; index_row < number_of_rows; index_row++) {
      column[index_row] = rows[index_row * number_of_titles + index_title];
    }
    columns.emplace_back(titlestring.substr(start, end - start), std::move(column));
    start = end + 1;
  }
  return columns;
}
