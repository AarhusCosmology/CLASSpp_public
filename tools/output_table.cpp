#include "output_table.h"

#include <algorithm>
#include <utility>

#include "common.h"

namespace {
// The room a table gets when it does not know its row count, before doubling.
constexpr std::size_t kMinimumCapacity = 16;
}  // namespace

OutputTable::OutputTable(std::string what) : what_(std::move(what)) {}

// The first row, a row that left the first row's order, and the first value after the
// buffer filled up.
void OutputTable::SetOutOfOrder(std::string_view name, double value) {
  if (expected_ == kFull) {
    Grow(2 * capacity_);
    expected_ = 0;
    Set(name, value);
    return;
  }
  const std::string key(name);
  const auto found      = index_.find(key);
  const std::size_t row = rows_ + 1;
  if (rows_ == 0) {
    class_test_severe(found != index_.end(),
                      "row 1 of the %s table sets '%s' twice. If a species instance's name "
                      "makes the title, rename the instance so its columns do not repeat "
                      "another's",
                      what_.c_str(),
                      key.c_str());
    index_.emplace(key, names_.size());
    names_.push_back(key);
    first_row_.push_back(value);
    return;
  }
  class_test_severe(found == index_.end(),
                    "'%s' is not a column of the %s table: after its first row, a row may only "
                    "set the columns the first row set",
                    key.c_str(),
                    what_.c_str());
  if (expected_ != kOutOfOrder) {
    // The row leaves the order here: the columns before it got their values in order.
    std::fill_n(last_row_set_.begin(), expected_, row);
    expected_ = kOutOfOrder;
  }
  const std::size_t column = found->second;
  class_test_severe(last_row_set_[column] == row,
                    "row %zu of the %s table sets '%s' twice. If a species instance's name makes "
                    "the title, rename the instance so its columns do not repeat another's",
                    row,
                    what_.c_str(),
                    key.c_str());
  last_row_set_[column]             = row;
  data_[column * capacity_ + rows_] = value;
}

void OutputTable::EndRowChecked() {
  const std::size_t row = rows_ + 1;
  if (rows_ == 0) {
    columns_ = names_.size();
    last_row_set_.assign(columns_, 0);
    capacity_ = std::max(reserved_, kMinimumCapacity);
    data_.assign(columns_ * capacity_, 0.);
    for (std::size_t column = 0; column < columns_; ++column) {
      data_[column * capacity_] = first_row_[column];
    }
    std::vector<double>().swap(first_row_);
  }
  else if (expected_ == kOutOfOrder) {
    for (std::size_t column = 0; column < columns_; ++column) {
      class_test_severe(last_row_set_[column] != row,
                        "row %zu of the %s table gives no value for '%s'",
                        row,
                        what_.c_str(),
                        names_[column].c_str());
    }
  }
  else {
    // A row in order that stopped short, or (kFull) a row that set nothing.
    const std::size_t missing = (expected_ == kFull) ? 0 : expected_;
    class_stop_severe("row %zu of the %s table gives no value for '%s'",
                      row,
                      what_.c_str(),
                      names_[missing].c_str());
  }
  ++rows_;
  StartRow();
}

void OutputTable::Grow(std::size_t capacity) {
  // The completed rows, and the current row's slots when there is room for it.
  const std::size_t keep = std::min(rows_ + 1, capacity_);
  std::vector<double> data(columns_ * capacity, 0.);
  for (std::size_t column = 0; column < columns_; ++column) {
    std::copy_n(data_.begin() + column * capacity_, keep, data.begin() + column * capacity);
  }
  data_.swap(data);
  capacity_ = capacity;
}

void OutputTable::ReserveRows(std::size_t rows) {
  reserved_ = rows;
  if (rows_ > 0 && rows > capacity_) {
    Grow(rows);
    if (expected_ == kFull) {
      expected_ = 0;
    }
  }
}

bool OutputTable::RowOpen() const {
  if (rows_ == 0) {
    return !names_.empty();
  }
  return expected_ != 0 && expected_ != kFull;
}

void OutputTable::Append(const OutputTable& other) {
  class_test_severe(RowOpen() || other.RowOpen(),
                    "a %s table can only be joined to another between rows",
                    what_.c_str());
  if (rows_ == 0) {
    // No columns yet: take those of other, keeping this table's name and notes.
    std::string what               = std::move(what_);
    std::vector<std::string> notes = std::move(notes_);
    *this                          = other;
    what_                          = std::move(what);
    notes_                         = std::move(notes);
  }
  else if (other.rows_ > 0) {
    class_test_severe(other.names_ != names_,
                      "two %s tables with different columns cannot be joined",
                      what_.c_str());
    if (rows_ + other.rows_ > capacity_) {
      Grow(rows_ + other.rows_);
    }
    for (std::size_t column = 0; column < columns_; ++column) {
      std::copy_n(other.data_.begin() + column * other.capacity_,
                  other.rows_,
                  data_.begin() + column * capacity_ + rows_);
    }
    rows_ += other.rows_;
    StartRow();
  }
  for (const std::string& note : other.notes_) {
    AddNote(note);
  }
}

void OutputTable::AddNote(std::string_view note) {
  if (std::find(notes_.begin(), notes_.end(), note) == notes_.end()) {
    notes_.emplace_back(note);
  }
}

OutputColumn OutputTable::column(const std::string& name) const {
  const auto found = index_.find(name);
  class_test_severe(found == index_.end(),
                    "the %s table has no column '%s'",
                    what_.c_str(),
                    name.c_str());
  return column(found->second);
}
