#ifndef __PARSER__
#define __PARSER__

#include <functional>
#include <map>
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

  /** Read an integer parameter.  Returns true and marks the key as read when
   *  found; returns false when absent.  Throws on parse error. */
  bool read_int(const std::string& name, int& value) const;

  /** Read a double parameter. */
  bool read_double(const std::string& name, double& value) const;

  /** Read a string parameter. */
  bool read_string(const std::string& name, std::string& value) const;

  /** Read a comma-separated list of doubles. Allocates and fills @p values. */
  bool read_list_of_doubles(const std::string& name, std::vector<double>& values) const;

  /** Read a comma-separated list of integers. Allocates and fills @p values. */
  bool read_list_of_integers(const std::string& name, std::vector<int>& values) const;

  /** Read a comma-separated list of strings. Fills @p values. */
  bool read_list_of_strings(const std::string& name, std::vector<std::string>& values) const;

  /** Mark every parameter as unread (used before a shooting iteration). */
  void mark_all_unread() const {
    read_params_.clear();
  }

  /** Return true if @p name has been marked as read. */
  bool was_read(const std::string& name) const {
    return read_params_.count(name) > 0;
  }

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

/**************************************************************/
/* Legacy C-style free functions - thin wrappers kept for     */
/* backward compatibility with existing call-sites.           */
/**************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

int parser_read_file(const char* filename, FileContent* pfc, ErrorMsg errmsg);

int parser_init(FileContent* pfc, int size, const char* filename, ErrorMsg errmsg);

int parser_read_line(char* line, int* is_data, char* name, char* value, ErrorMsg errmsg);

int parser_read_int(FileContent* pfc, const char* name, int* value, int* found, ErrorMsg errmsg);

int parser_read_double(
    FileContent* pfc, const char* name, double* value, int* found, ErrorMsg errmsg);

int parser_read_string(
    FileContent* pfc, const char* name, FileArg* value, int* found, ErrorMsg errmsg);

int parser_cat(const FileContent* pfc1,
               const FileContent* pfc2,
               FileContent* pfc3,
               ErrorMsg errmsg);

#ifdef __cplusplus
}
#endif

#endif
