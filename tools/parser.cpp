#include "parser.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>

/*****************************************************************************
 * FileContent – C++ implementation
 *****************************************************************************/

/* static */ FileContent FileContent::from_file(const std::string& filename) {
  FILE* fp = std::fopen(filename.c_str(), "r");
  if (!fp) {
    throw std::invalid_argument("Cannot open file '" + filename + "': " + std::strerror(errno));
  }

  FileContent fc;
  fc.filename_ = filename;

  char line_buf[_LINE_LENGTH_MAX_];
  while (std::fgets(line_buf, _LINE_LENGTH_MAX_, fp) != nullptr) {
    std::string name, value;
    if (!parse_line(line_buf, name, value))
      continue;

    if (fc.params_.count(name)) {
      std::fclose(fp);
      throw std::invalid_argument("Multiple entries of parameter '" + name + "' in file '" +
                                  filename + "'");
    }
    fc.keys_.push_back(name);
    fc.params_[name] = value;
  }

  if (fc.params_.empty()) {
    std::fclose(fp);
    throw std::invalid_argument("No readable input in file '" + filename + "'");
  }

  std::fclose(fp);
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

int parser_read_file(const char* filename, FileContent* pfc, ErrorMsg errmsg) {
  try {
    *pfc = FileContent::from_file(filename);
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

int parser_init(FileContent* pfc, int /*size*/, const char* /*filename*/, ErrorMsg /*errmsg*/) {
  /* With the new FileContent class the size is managed dynamically.
   * Reset to a fresh instance without inserting any placeholder parameter,
   * so the parameter map remains empty until real values are added. */
  *pfc = FileContent();
  return _SUCCESS_;
}

int parser_read_line(char* line, int* is_data, char* name, char* value, ErrorMsg errmsg) {
  std::string sname, svalue;
  if (!FileContent::parse_line(line, sname, svalue)) {
    *is_data = _FALSE_;
    return _SUCCESS_;
  }
  class_test(sname.size() >= _ARGUMENT_LENGTH_MAX_,
             errmsg,
             "name starting by '%s' too long; shorten it or increase _ARGUMENT_LENGTH_MAX_",
             sname.c_str());
  class_test(svalue.size() >= _ARGUMENT_LENGTH_MAX_,
             errmsg,
             "value starting by '%s' too long; shorten it or increase _ARGUMENT_LENGTH_MAX_",
             svalue.c_str());
  std::strncpy(name, sname.c_str(), _ARGUMENT_LENGTH_MAX_ - 1);
  name[_ARGUMENT_LENGTH_MAX_ - 1] = '\0';
  std::strncpy(value, svalue.c_str(), _ARGUMENT_LENGTH_MAX_ - 1);
  value[_ARGUMENT_LENGTH_MAX_ - 1] = '\0';
  *is_data                         = _TRUE_;
  return _SUCCESS_;
}

int parser_read_int(FileContent* pfc, const char* name, int* value, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_int(name, *value) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

int parser_read_double(
    FileContent* pfc, const char* name, double* value, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_double(name, *value) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

int parser_read_string(
    FileContent* pfc, const char* name, FileArg* value, int* found, ErrorMsg errmsg) {
  std::string s;
  try {
    *found = pfc->read_string(name, s) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  if (*found == _TRUE_) {
    class_test(s.size() >= _ARGUMENT_LENGTH_MAX_,
               errmsg,
               "value of '%s' too long; increase _ARGUMENT_LENGTH_MAX_",
               name);
    std::strncpy(*value, s.c_str(), _ARGUMENT_LENGTH_MAX_ - 1);
    (*value)[_ARGUMENT_LENGTH_MAX_ - 1] = '\0';
  }
  return _SUCCESS_;
}

int parser_cat(const FileContent* pfc1,
               const FileContent* pfc2,
               FileContent* pfc3,
               ErrorMsg errmsg) {
  try {
    *pfc3 = *pfc1 + *pfc2;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}
