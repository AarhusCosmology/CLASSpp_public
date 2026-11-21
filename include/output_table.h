#ifndef OUTPUT_TABLE_H
#define OUTPUT_TABLE_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

/** One column of an OutputTable: its values in row order. A view, valid until the table
 *  changes. */
class OutputColumn {
 public:
  /** For Cython, which declares a temporary before assigning a returned column to it. */
  OutputColumn() = default;
  OutputColumn(const double* data, std::size_t size) : data_(data), size_(size) {}

  const double* data() const {
    return data_;
  }
  std::size_t size() const {
    return size_;
  }
  double operator[](std::size_t row) const {
    return data_[row];
  }
  const double* begin() const {
    return data_;
  }
  const double* end() const {
    return data_ + size_;
  }

 private:
  const double* data_ = nullptr;
  std::size_t size_   = 0;
};

/**
 * An output table (background.dat, a transfer function table, ...) as named columns, filled
 * one row at a time with Set(name, value), like a Python dict of columns.
 *
 * Each value travels with its name, so nothing else lists the columns, and nothing can put
 * a value under another column's title. The columns keep the order in which the first row
 * named them: the order of the .dat files and of classy's dicts. A later row may set its
 * names in any order. Every inconsistency is a severe error that names the table: a name set
 * twice in one row, a name the first row did not set, and a row that leaves out a column.
 *
 * Every producer sets the names of a later row in the first row's order, so that is the path
 * made fast. Set compares the name with the column expected next, inline (for a literal
 * name the compiler knows its length), and stores the value in its slot. A row that goes in
 * order from the first column to the last has set every column exactly once, so EndRow only
 * checks that it got to the end. A row that leaves the order takes an out-of-line path: the
 * hash map finds the column, and per-column stamps catch a name set twice and a column left
 * out, for that row only.
 *
 * The values are stored column by column in one buffer with room for capacity_ rows:
 * column c starts at data_[c * capacity_]. A column is then contiguous, as classy wants it,
 * and a value in order costs one store; the room is checked once per row. A producer that
 * knows its row count reserves it, and the buffer doubles otherwise. The first row is held
 * apart until it ends, since the number of columns is not known before.
 *
 * Design: docs/superpowers/specs/2026-11-20-output-table-design.md. This storage came from
 * an independent review of two others (2026-11-21), measured there and in classy.
 */
class OutputTable {
 public:
  /** `what` names the table in error messages: "background", "transfer functions", ... */
  explicit OutputTable(std::string what);
  /** For Cython, which declares a temporary before assigning a returned table to it. */
  OutputTable() = default;

  /** A value of the current row, in the column called `name`. */
  void Set(std::string_view name, double value) {
    if (expected_ < columns_) {
      const std::string& next = names_[expected_];
      if (next.size() == name.size() && std::memcmp(next.data(), name.data(), name.size()) == 0) {
        data_[expected_ * capacity_ + rows_] = value;
        ++expected_;
        return;
      }
    }
    SetOutOfOrder(name, value);
  }
  /** Closes the current row; refuses it if a column got no value. */
  void EndRow() {
    if (expected_ == columns_) {
      ++rows_;
      StartRow();
      return;
    }
    EndRowChecked();
  }
  /** Room for `rows` rows, for a producer that knows how many it writes. */
  void ReserveRows(std::size_t rows);
  /** The rows of `other`, after this table's own. The columns must be the same, in the same
   *  order; a table without columns takes those of `other`. */
  void Append(const OutputTable& other);
  /** A header line for the file; a note given twice is kept once. */
  void AddNote(std::string_view note);

  const std::string& what() const {
    return what_;
  }
  /** The column names, in the order the first row set them. */
  const std::vector<std::string>& names() const {
    return names_;
  }
  OutputColumn column(std::size_t index) const {
    return {data_.data() + index * capacity_, rows_};
  }
  OutputColumn column(const std::string& name) const;
  /** The number of completed rows. */
  std::size_t rows() const {
    return rows_;
  }
  const std::vector<std::string>& notes() const {
    return notes_;
  }

 private:
  // Values of expected_ that send every Set to SetOutOfOrder.
  static constexpr std::size_t kOutOfOrder = SIZE_MAX;      // the first row, or a row out of order
  static constexpr std::size_t kFull       = SIZE_MAX - 1;  // no room for a new row yet

  void SetOutOfOrder(std::string_view name, double value);
  void EndRowChecked();
  void StartRow() {
    expected_ = (rows_ < capacity_ || columns_ == 0) ? 0 : kFull;
  }
  bool RowOpen() const;
  void Grow(std::size_t capacity);

  std::string what_;
  std::vector<std::string> names_;
  std::unordered_map<std::string, std::size_t> index_;
  std::vector<double> data_;       // column c, row r at data_[c * capacity_ + r]
  std::vector<double> first_row_;  // the first row, until it ends
  // Per column, the row (counted from 1) in which a row out of order last set it.
  std::vector<std::size_t> last_row_set_;
  std::size_t columns_  = 0;            // fixed when the first row ends
  std::size_t rows_     = 0;            // completed rows
  std::size_t capacity_ = 0;            // rows data_ has room for
  std::size_t reserved_ = 0;            // rows asked for by ReserveRows
  std::size_t expected_ = kOutOfOrder;  // the column a row in order sets next
  std::vector<std::string> notes_;
};

#endif
