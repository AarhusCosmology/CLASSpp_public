#include "parser.h"

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

/*****************************************************************************
 * FileContent – C++ implementation
 *****************************************************************************/

/* static */ FileContent FileContent::from_file(const std::string& filename) {
  // std::ifstream is not guaranteed to set errno on failure, so reset it first
  // and only report a strerror reason when the open actually set one.
  errno = 0;
  std::ifstream in(filename);
  if (!in) {
    const std::string reason = errno ? std::strerror(errno) : "unable to open file";
    throw std::invalid_argument("Cannot open file '" + filename + "': " + reason);
  }
  FileContent fc;
  fc.filename_ = filename;

  std::string line;
  while (std::getline(in, line)) {
    std::string name, value;
    if (!parse_line(line, name, value))
      continue;
    if (fc.params_.count(name)) {
      throw std::invalid_argument("Multiple entries of parameter '" + name + "' in file '" +
                                  filename + "'");
    }
    fc.keys_.push_back(name);
    fc.params_[name] = value;
  }
  if (fc.params_.empty()) {
    throw std::invalid_argument("No readable input in file '" + filename + "'");
  }
  return fc;
}

void FileContent::set(const std::string& name, const std::string& value) {
  if (!params_.count(name)) {
    keys_.push_back(name);
  }
  params_[name] = value;
  read_params_.erase(name);
}

FileContent& FileContent::operator+=(const FileContent& other) {
  for (const auto& key : other.keys_) {
    if (params_.count(key)) {
      throw std::invalid_argument("Multiple entries of parameter '" + key + "' in files '" +
                                  filename_ + "' and '" + other.filename_ + "'");
    }
    keys_.push_back(key);
    params_[key] = other.params_.at(key);
  }
  if (filename_.empty()) {
    filename_ = other.filename_;
  }
  else if (!other.filename_.empty()) {
    filename_ += " or " + other.filename_;
  }
  return *this;
}

std::vector<std::string> FileContent::instances_with(const std::string& field,
                                                     const std::string& value) const {
  std::vector<std::string> out;
  const std::string suffix = "." + field;
  for (const std::string& key : keys_) {
    if (key.size() <= suffix.size())
      continue;
    if (key.compare(key.size() - suffix.size(), suffix.size(), suffix) != 0)
      continue;
    const std::string name = key.substr(0, key.size() - suffix.size());
    if (name.empty())
      continue;
    char c0 = name[0];
    if (!(std::isalpha(static_cast<unsigned char>(c0)) || c0 == '_'))
      continue;
    bool ok = true;
    for (char c : name) {
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_')) {
        ok = false;
        break;
      }
    }
    if (!ok)
      continue;
    auto it = params_.find(key);
    if (it == params_.end())
      continue;
    if (it->second == value)
      out.push_back(name);
  }
  return out;
}

std::vector<std::string> FileContent::unread_parameters() const {
  std::vector<std::string> unread;
  for (const auto& key : keys_) {
    if (!read_params_.count(key))
      unread.push_back(key);
  }
  return unread;
}

void FileContent::for_each(
    const std::function<void(const std::string&, const std::string&, bool)>& fn) const {
  for (const auto& key : keys_) {
    fn(key, params_.at(key), read_params_.count(key) > 0);
  }
}

bool FileContent::read_int(const std::string& name, int& value) const {
  auto it = params_.find(name);
  if (it == params_.end())
    return false;
  if (std::sscanf(it->second.c_str(), "%d", &value) != 1) {
    throw std::invalid_argument("Cannot read integer value of parameter '" + name + "' in file '" +
                                filename_ + "'");
  }
  read_params_.insert(name);
  return true;
}

bool FileContent::read_double(const std::string& name, double& value) const {
  auto it = params_.find(name);
  if (it == params_.end())
    return false;
  if (std::sscanf(it->second.c_str(), "%lg", &value) != 1) {
    throw std::invalid_argument("Cannot read double value of parameter '" + name + "' in file '" +
                                filename_ + "'");
  }
  read_params_.insert(name);
  return true;
}

bool FileContent::read_string(const std::string& name, std::string& value) const {
  auto it = params_.find(name);
  if (it == params_.end())
    return false;
  value = it->second;
  read_params_.insert(name);
  return true;
}

bool FileContent::read_list_of_doubles(const std::string& name, std::vector<double>& values) const {
  auto it = params_.find(name);
  if (it == params_.end())
    return false;

  const auto parts = split_csv(it->second);
  values.resize(parts.size());
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (std::sscanf(parts[i].c_str(), "%lg", &values[i]) != 1) {
      throw std::invalid_argument("Cannot read double entry " + std::to_string(i + 1) +
                                  " of parameter '" + name + "' in file '" + filename_ + "'");
    }
  }
  read_params_.insert(name);
  return true;
}

bool FileContent::read_list_of_integers(const std::string& name, std::vector<int>& values) const {
  auto it = params_.find(name);
  if (it == params_.end())
    return false;

  const auto parts = split_csv(it->second);
  values.resize(parts.size());
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (std::sscanf(parts[i].c_str(), "%d", &values[i]) != 1) {
      throw std::invalid_argument("Cannot read integer entry " + std::to_string(i + 1) +
                                  " of parameter '" + name + "' in file '" + filename_ + "'");
    }
  }
  read_params_.insert(name);
  return true;
}

bool FileContent::read_list_of_strings(const std::string& name,
                                       std::vector<std::string>& values) const {
  auto it = params_.find(name);
  if (it == params_.end())
    return false;
  values = split_csv(it->second);
  read_params_.insert(name);
  return true;
}

/* static */ bool FileContent::parse_line(const std::string& line,
                                          std::string& name,
                                          std::string& value) {
  /* A valid data line must contain '=' */
  const auto eq_pos = line.find('=');
  if (eq_pos == std::string::npos)
    return false;

  /* Ignore the line if '#' appears before (or immediately after) '=' */
  const auto hash_pos = line.find('#');
  if (hash_pos != std::string::npos && hash_pos < eq_pos + 2)
    return false;

  /* Extract the name: trim whitespace and optional surrounding quotes */
  auto trim_quotes = [](std::string s) -> std::string {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\'' || s.front() == '"'))
      s.erase(s.begin());
    while (!s.empty() && (s.back() == ' ' || s.back() == '\'' || s.back() == '"'))
      s.pop_back();
    return s;
  };

  name = trim_quotes(line.substr(0, eq_pos));
  if (name.empty())
    return false;

  /* Extract the value: text after '=', up to '#' or end of line, trimmed */
  const auto value_end  = (hash_pos != std::string::npos) ? hash_pos : line.size();
  std::string raw_value = line.substr(eq_pos + 1, value_end - eq_pos - 1);

  /* Trim trailing whitespace / newline */
  while (!raw_value.empty() && (unsigned char) raw_value.back() <= ' ')
    raw_value.pop_back();
  /* Trim leading whitespace */
  std::size_t start = 0;
  while (start < raw_value.size() && raw_value[start] == ' ')
    ++start;
  raw_value = raw_value.substr(start);

  if (raw_value.empty())
    return false;

  value = std::move(raw_value);
  return true;
}

/* static */ std::vector<std::string> FileContent::split_csv(const std::string& s) {
  std::vector<std::string> parts;
  std::istringstream ss(s);
  std::string token;
  while (std::getline(ss, token, ',')) {
    /* trim whitespace */
    std::size_t start = token.find_first_not_of(' ');
    std::size_t end   = token.find_last_not_of(' ');
    if (start != std::string::npos) {
      parts.push_back(token.substr(start, end - start + 1));
    }
  }
  return parts;
}

/*****************************************************************************
 * Legacy C-style wrapper functions
 * These delegate to the FileContent class so that existing call-sites in
 * input_module.cpp (which use the class_call / class_test macros) keep working
 * without modification.
 *****************************************************************************/

int parser_read_file(const char* filename, FileContent* pfc) {
  try {
    *pfc = FileContent::from_file(filename);
  }
  catch (const std::exception& e) {
    class_stop("%s", e.what());
  }
  return _SUCCESS_;
}

int parser_read_int(FileContent* pfc, const char* name, int* value, int* found) {
  try {
    *found = pfc->read_int(name, *value) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop("%s", e.what());
  }
  return _SUCCESS_;
}

int parser_read_double(FileContent* pfc, const char* name, double* value, int* found) {
  try {
    *found = pfc->read_double(name, *value) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop("%s", e.what());
  }
  return _SUCCESS_;
}

int parser_read_string(FileContent* pfc, const char* name, std::string& value, int* found) {
  try {
    *found = pfc->read_string(name, value) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop("%s", e.what());
  }
  return _SUCCESS_;
}

int parser_cat(const FileContent* pfc1, const FileContent* pfc2, FileContent* pfc3) {
  try {
    *pfc3 = *pfc1 + *pfc2;
  }
  catch (const std::exception& e) {
    class_stop("%s", e.what());
  }
  return _SUCCESS_;
}
