// OutputTable: named columns in insertion order, filled one row at a time, with every
// inconsistency between rows a severe error
// (docs/superpowers/specs/2026-11-20-output-table-design.md).
#include "output_table.h"

#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// The message of the severe error `fill` throws, or "" if it throws none.
template <class F>
std::string SevereMessage(F fill) {
  try {
    fill();
  }
  catch (const std::invalid_argument& error) {
    return error.what();
  }
  return "";
}

bool Contains(const std::string& text, const std::string& part) {
  return text.find(part) != std::string::npos;
}

std::vector<double> Values(const OutputColumn& column) {
  return {column.begin(), column.end()};
}

// Rows of z = row, x_e = 2 row, in order.
void FillInOrder(OutputTable& table, int from, int to) {
  for (int row = from; row < to; ++row) {
    table.Set("z", row);
    table.Set("x_e", 2. * row);
    table.EndRow();
  }
}

void ColumnsKeepTheOrderOfTheFirstRow() {
  OutputTable table("background");
  for (int row = 0; row < 3; ++row) {
    table.Set("z", row);
    table.Set("H [1/Mpc]", 10. * row);
    table.Set("a", 0.5 * row);
    table.EndRow();
  }
  assert((table.names() == std::vector<std::string>{"z", "H [1/Mpc]", "a"}));
  assert(table.rows() == 3);
  assert((Values(table.column(1)) == std::vector<double>{0., 10., 20.}));
  assert((Values(table.column("a")) == std::vector<double>{0., 0.5, 1.}));
}

// A value goes to its column by name, so a later row may set its names in any order.
void ALaterRowMaySetItsNamesInAnyOrder() {
  OutputTable table("thermodynamics");
  table.Set("z", 1.);
  table.Set("x_e", 2.);
  table.EndRow();
  table.Set("x_e", 4.);
  table.Set("z", 3.);
  table.EndRow();
  assert((Values(table.column("z")) == std::vector<double>{1., 3.}));
  assert((Values(table.column("x_e")) == std::vector<double>{2., 4.}));
}

void ANameTwiceInOneRowIsAnError() {
  OutputTable first("background");
  assert(Contains(SevereMessage([&] {
                    first.Set("(.)rho_b", 1.);
                    first.Set("(.)rho_b", 2.);
                  }),
                  "'(.)rho_b'"));
  OutputTable later("background");
  later.Set("z", 0.);
  later.Set("a", 1.);
  later.EndRow();
  assert(Contains(SevereMessage([&] {
                    later.Set("z", 1.);
                    later.Set("z", 2.);
                  }),
                  "'z'"));
}

// Set compares a name with the column expected next before it looks the name up. That
// column may already have a value in this row, after a row that went out of order.
void ARepeatOnTheExpectedColumnIsAnError() {
  OutputTable table("background");
  for (const char* name : {"a", "b", "c"}) {
    table.Set(name, 0.);
  }
  table.EndRow();
  assert(Contains(SevereMessage([&] {
                    table.Set("b", 1.);  // looked up; c is expected next
                    table.Set("c", 1.);
                    table.Set("a", 1.);  // looked up; b is expected next, and has its value
                    table.Set("b", 2.);
                  }),
                  "sets 'b' twice"));
}

// The buffer doubles when a table does not know its row count, and reserving room, before
// the first row or after it, changes no value.
void GrowingAndReservingKeepTheValues() {
  OutputTable grown("thermodynamics");
  FillInOrder(grown, 0, 1000);
  OutputTable reserved("thermodynamics");
  reserved.ReserveRows(1000);
  FillInOrder(reserved, 0, 1000);
  OutputTable late("thermodynamics");
  FillInOrder(late, 0, 10);
  late.ReserveRows(1000);
  FillInOrder(late, 10, 1000);
  for (const OutputTable* table : {&grown, &reserved, &late}) {
    assert(table->rows() == 1000);
    for (int row = 0; row < 1000; ++row) {
      assert(table->column(0)[row] == row && table->column(1)[row] == 2. * row);
    }
  }
}

// A row out of order goes through the lookup, and the next row in order is fast again.
void ARowOutOfOrderThenOneInOrder() {
  OutputTable table("background");
  FillInOrder(table, 0, 16);  // fills the first room of 16 rows, so row 17 grows the buffer
  table.Set("x_e", 42.);
  table.Set("z", 16.);
  table.EndRow();
  FillInOrder(table, 17, 40);
  assert(table.rows() == 40);
  assert(table.column("z")[16] == 16. && table.column("x_e")[16] == 42.);
  assert(table.column("z")[39] == 39. && table.column("x_e")[39] == 78.);
}

// A row in order that stops short, a row out of order that misses a column, and a row that
// sets nothing just after the buffer filled up.
void EveryWayOfMissingAColumnIsAnError() {
  OutputTable short_row("background");
  FillInOrder(short_row, 0, 2);
  assert(Contains(SevereMessage([&] {
                    short_row.Set("z", 1.);
                    short_row.EndRow();
                  }),
                  "no value for 'x_e'"));
  OutputTable out_of_order("background");
  FillInOrder(out_of_order, 0, 2);
  assert(Contains(SevereMessage([&] {
                    out_of_order.Set("x_e", 1.);
                    out_of_order.EndRow();
                  }),
                  "no value for 'z'"));
  OutputTable full("background");
  full.ReserveRows(16);
  FillInOrder(full, 0, 16);
  assert(Contains(SevereMessage([&] { full.EndRow(); }), "no value for 'z'"));
}

// The tables of two initial conditions are joined row after row; their columns must agree.
void AppendJoinsTheRows() {
  OutputTable ad("scalar k_output_values");
  FillInOrder(ad, 0, 20);
  OutputTable cdi("scalar k_output_values");
  FillInOrder(cdi, 100, 130);
  OutputTable joined = ad;
  joined.Append(cdi);
  assert(joined.rows() == 50 && joined.names() == ad.names());
  assert(joined.column("z")[19] == 19. && joined.column("z")[20] == 100.);
  assert(joined.column("x_e")[49] == 2. * 129);
  OutputTable empty("scalar k_output_values");
  empty.Append(cdi);
  assert(Values(empty.column("x_e")) == Values(cdi.column("x_e")));
  OutputTable other("scalar k_output_values");
  other.Set("z", 0.);
  other.EndRow();
  assert(Contains(SevereMessage([&] { joined.Append(other); }), "different columns"));
}

void ANewNameAfterTheFirstRowIsAnError() {
  OutputTable table("transfer functions");
  table.Set("k (h/Mpc)", 0.1);
  table.EndRow();
  const std::string message = SevereMessage([&] { table.Set("d_tot", 1.); });
  assert(Contains(message, "'d_tot'") && Contains(message, "transfer functions"));
}

void ARowWithAMissingColumnIsAnError() {
  OutputTable table("perturbations");
  table.Set("tau [Mpc]", 1.);
  table.Set("a", 1e-5);
  table.EndRow();
  table.Set("tau [Mpc]", 2.);
  const std::string message = SevereMessage([&] { table.EndRow(); });
  assert(Contains(message, "'a'") && Contains(message, "perturbations"));
}

void AnEmptyTable() {
  const OutputTable table("primordial");
  assert(table.names().empty() && table.rows() == 0 && table.notes().empty());
}

void ANoteIsKeptOnce() {
  OutputTable table("background");
  table.AddNote("phi_<N> in units of M_pl");
  table.AddNote("rho in units of 3 M_pl^2 / Mpc^2");
  table.AddNote("phi_<N> in units of M_pl");
  assert((table.notes() == std::vector<std::string>{"phi_<N> in units of M_pl",
                                                    "rho in units of 3 M_pl^2 / Mpc^2"}));
}

}  // namespace

int main() {
  ColumnsKeepTheOrderOfTheFirstRow();
  ALaterRowMaySetItsNamesInAnyOrder();
  ANameTwiceInOneRowIsAnError();
  ARepeatOnTheExpectedColumnIsAnError();
  GrowingAndReservingKeepTheValues();
  ARowOutOfOrderThenOneInOrder();
  EveryWayOfMissingAColumnIsAnError();
  AppendJoinsTheRows();
  ANewNameAfterTheFirstRowIsAnError();
  ARowWithAMissingColumnIsAnError();
  AnEmptyTable();
  ANoteIsKeptOnce();
  return 0;
}
