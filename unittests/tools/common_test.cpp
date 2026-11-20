#include "common.h"

#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>

// The rows of a .dat table become one column per title, in the order of the titles.
static void test_split_into_columns_transposes_the_rows() {
  std::string titles;
  class_store_columntitle(titles, "z", true);
  class_store_columntitle(titles, "H [1/Mpc]", true);
  class_store_columntitle(titles, "skipped", false);
  const std::vector<double> rows = {1., 2., 3., 4., 5., 6.};

  const NamedColumns columns = SplitIntoColumns(titles, rows.data(), 3);

  assert(columns.size() == 2);
  assert(columns[0].first == "z");
  assert(columns[1].first == "H [1/Mpc]");
  assert((columns[0].second == std::vector<double>{1., 3., 5.}));
  assert((columns[1].second == std::vector<double>{2., 4., 6.}));
}

// A table without rows still names its columns; one without titles has no columns.
static void test_split_into_columns_of_empty_tables() {
  const NamedColumns no_rows = SplitIntoColumns("a\tb\t", nullptr, 0);
  assert(no_rows.size() == 2);
  assert(no_rows[1].first == "b");
  assert(no_rows[1].second.empty());

  assert(SplitIntoColumns("", nullptr, 0).empty());
}

// A title is a column's key in classy's dicts: a repeated one would hide a column.
static void test_split_into_columns_refuses_repeated_titles() {
  bool refused = false;
  try {
    SplitIntoColumns("z\tH\tz\t", nullptr, 0);
  }
  catch (const std::invalid_argument& error) {
    refused = std::string(error.what()).find("'z'") != std::string::npos;
  }
  assert(refused);
}

int main() {
  test_split_into_columns_transposes_the_rows();
  test_split_into_columns_of_empty_tables();
  test_split_into_columns_refuses_repeated_titles();
  return 0;
}
