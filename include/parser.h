#ifndef __PARSER__
#define __PARSER__

#include <functional>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.h"

/**
 * Holds the parsed contents of one or more .ini/.pre files.
 *
 * Internally uses std::map for O(log N) parameter lookup and a std::vector to
 * preserve insertion order.  A std::set tracks which parameters have been
 * consumed so that callers can warn about unused entries.
 */
class FileContent {
 public:
  bool is_shooting = false;

  FileContent()                              = default;
  ~FileContent()                             = default;
  FileContent(const FileContent&)            = default;
  FileContent(FileContent&&)                 = default;
  FileContent& operator=(const FileContent&) = default;
  FileContent& operator=(FileContent&&)      = default;

  /** Load parameters from an .ini/.pre file.  Throws std::invalid_argument on
   *  I/O errors or duplicate keys within the same file. */
  static FileContent from_file(const std::string& filename);

  /** Insert or overwrite a single parameter. Marks the key as unread. */
  void set(const std::string& name, const std::string& value);

  /** Merge another FileContent into this one.  Throws std::invalid_argument if
   *  a key already present in *this also appears in @p other. */
  FileContent& operator+=(const FileContent& other);
  friend FileContent operator+(FileContent lhs, const FileContent& rhs) {
    return lhs += rhs;
  }

  /** Number of stored parameters. */
  int size() const {
    return static_cast<int>(params_.size());
  }

  /** Source filename (or "file1 or file2" after a merge). */
  const std::string& get_filename() const {
    return filename_;
  }

  /** Typed accessor: the value if @p name is present (marking it read), or
   *  std::nullopt if absent. Throws on a present-but-unparseable value.
   *  Supported T: int, double, std::string, std::vector<{int,double,std::string}>. */
  template <class T>
  std::optional<T> get(const std::string& name) const;

  /** Read-with-default: get<T>(name) if present, else @p fallback. */
  template <class T>
  T get_or(const std::string& name, T fallback) const {
    if (auto v = get<T>(name))
      return *v;
    return fallback;
  }
  /** Overload so a string-literal fallback does not deduce T = const char*. */
  std::string get_or(const std::string& name, const char* fallback) const;

  /** Mark every parameter as unread (used before a shooting iteration). */
  void mark_all_unread() const {
    read_params_.clear();
  }

  /** Return true if @p name has been marked as read. */
  bool was_read(const std::string& name) const {
    return read_params_.count(name) > 0;
  }

  /** Return every instance name N such that the entry "N.<field>" has the
   *  given value. The dot is a literal separator; N must match the instance
   *  regex [A-Za-z_][A-Za-z0-9_]*. Results are returned in insertion order.
   *  Does NOT mark anything as read. */
  std::vector<std::string> instances_with(const std::string& field, const std::string& value) const;

  /** Return the names of all parameters that have not yet been read. */
  std::vector<std::string> unread_parameters() const;

 private:
  std::string filename_;
  std::vector<std::string> keys_; /**< insertion-order key list */
  std::map<std::string, std::string> params_;
  mutable std::set<std::string> read_params_;

 public:
  /** Invoke @p fn for every parameter in insertion order.
   *  Arguments: parameter name, value string, whether it was read.
   *  Declared in the second public section so the Cython wrapper generator
   *  (which stops at the first private:) does not attempt to parse it. */
  void for_each(
      const std::function<void(const std::string& name, const std::string& value, bool read)>& fn)
      const;

  /** Parse a single .ini line into name/value.  Returns true when the line
   *  contains a valid key=value pair.  Public for the legacy wrapper. */
  static bool parse_line(const std::string& line, std::string& name, std::string& value);

  /** Split a comma-separated value string into trimmed substrings. */
  static std::vector<std::string> split_csv(const std::string& s);
};

// Explicit specialization declarations — definitions are in parser.cpp.
// These prevent implicit instantiation from the primary template in any TU.
template <>
std::optional<int> FileContent::get<int>(const std::string& name) const;
template <>
std::optional<double> FileContent::get<double>(const std::string& name) const;
template <>
std::optional<std::string> FileContent::get<std::string>(const std::string& name) const;
template <>
std::optional<std::vector<double>> FileContent::get<std::vector<double>>(
    const std::string& name) const;
template <>
std::optional<std::vector<int>> FileContent::get<std::vector<int>>(const std::string& name) const;
template <>
std::optional<std::vector<std::string>> FileContent::get<std::vector<std::string>>(
    const std::string& name) const;

/**************************************************************/
/* Legacy C-style free functions.                             */
/**************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

int parser_read_file(const char* filename, FileContent* pfc);

int parser_cat(const FileContent* pfc1, const FileContent* pfc2, FileContent* pfc3);

#ifdef __cplusplus
}
#endif

#endif
