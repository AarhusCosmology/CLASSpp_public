# Design: one insertion-ordered table for every output, filled one named value at a time

**Status (2026-11-21): implemented on branch `output-table`, on master after #504 (squash
1aebd2b4). Approved by Thomas on 2026-11-20.** This PR comes before the nanobind prototype.
Where the implementation differs from the text below:
- **The storage of section 1 is superseded** [M]. Section 1 keeps one vector per column, a
  per-column stamp for each value, and the hash map. Every producer sets a later row in the first
  row's order, so the profile showed no hashing, only the bookkeeping of each value: the stamp,
  the push_back into one of many vectors, and an out-of-line call. An independent review
  (2026-11-21) compared that with a row-major buffer and with design C in a standalone replica.
  C was the fastest, and is what shipped: the values column by column in one buffer with room
  for a number of rows (column c at `data_[c * capacity_]`). A value in order costs an inline
  comparison with the name expected next and one store. EndRow checks once per row that the row
  reached the last column, which, in order, means every column got exactly one value. Stamps and
  the hash map serve only the first row, a row out of order and lookups by name. Columns stay
  contiguous, so classy copies each with one memcpy. `column(i)` returns a view
  (`OutputColumn`), and the tables of several initial conditions are joined by `Append`.
  Producers that know their row count reserve it; the others double. Not done: species that
  build column names per row ("(.)rho_" + name()) could build them once, worth about 0.03 ms
  of `get_background()` (below).
- **k_output_values tables per initial condition** [M]. Section 2 has one table per k and mode.
  Found while building: with `ic = ad,cdi`, k_output_values and 8 threads, #504's `./class` died
  with SIGTRAP in 8 runs out of 8, and ran with 1 thread. perturb_init queues every initial
  condition and k at once, so the ad and cdi integrations of one k ran together and appended to
  one buffer. class_public's loop over initial conditions is outside its parallel region. Each
  initial condition now has its own table, and the file and classy get their rows one after the
  other, in the order of index_ic: what one thread wrote, and what class_public writes. Pinned
  by `unittests/source/k_output_tables_test.cpp`.
- **`AddTensorPerturbationOutput` also takes k** [C]: the workspace holds k only for scalars.
- **d_tot and t_tot** are set under the flags that fill them (their titles were unconditional);
  CAMB's `-T_tot/k2` reads delta_tot only when it is a source [C].
- **A repeated name in a k_output_values table** fails during the integration, in the first
  row (section 5, applied to tables that are built while integrating).
- **The faster path of "Left to the builder" was needed** [M]. With section 1's storage and
  `Set` out of line, `get_thermodynamics()` was 2.35x slower than #504; with the comparison
  inline, 1.9x; with the storage above, 1.05x.

Verification [M], against #504 (e5f07430, the tree of 1aebd2b4):
- .dat files: 128 runs, 991 files, byte-identical, except (i) the 6 `primordial_Pk.dat` of runs
  without perturbations (an empty table has no column titles), (ii) one value of
  `perturbations_k0_t.dat` in the tensors + ncdm run: shear_ncdm, in its last printed digit
  (3.7e-13 relative). With `CLASS_FAST_MATH=OFF` in both builds, all 991 are identical: the loop
  that moved from the module into NCDMSpecies rounds differently under fast-math.
- Decaying NCDM with k_output_values (`benchmarks/output_table/custom_dncdm/`, outside the
  repository): the integrated representation byte-identical; psd with inverse decays and the
  l-moment columns differs in 12 of 6039 x 95 values, by at most 1e-13 of the column's scale
  (the `Mlcancel` ratios), and is identical with fast-math off.
- classy: 151 calls in 11 cosmologies identical (`benchmarks/classy_glue/compare_tables.py`, a
  local harness outside the repository, like the timing script below).
- ctest 73/73; TEST_LEVEL=1 scenarios pass; the other Python tests fail exactly where #504 does
  (the committed goldens of `test_background_columns.py` and `test_transfer_columns.py`).
- Timing, best of 3x15, ms, #504 -> branch, on AC power (`benchmarks/classy_glue/timing_tables.py`,
  one process per build, alternating): `get_background()` 0.37 -> 0.40 (1.07x),
  `get_thermodynamics()` 1.32 -> 1.39 (1.05x), `get_perturbations()` 0.17 -> 0.08 (0.46x),
  `get_transfer()` 0.042 -> 0.038 (0.9x). With section 1's storage they were 1.5x, 1.9x, 0.4x
  and 1.0x.

## Why

Every output table of CLASS++ is built twice: once as a tab-separated string of titles, and
once as a flat buffer of values, row after row. Two code paths, or one path run in two modes,
must then agree on the order and on the conditions of every column:
- `class_store_columntitle` / `class_store_double`, called with the same conditions in the
  same order;
- `BackgroundColumnWriter`, with one virtual for the titles and one for the values;
- `PerturbColumnWriter`, one species function run once in "title mode" and once per row in
  "data mode".

A disagreement shifts every later column silently. Nothing checks it, and reading the code
does not settle it. Two examples from the code today:
- The transfer table writes the title `d_tot` unconditionally, but its value only when
  `has_source_delta_tot_`; `t_tot` likewise with `has_source_theta_tot_`. They are in step only
  because `perturb_indices` sets both flags wherever their transfers are requested
  (perturbations_module.cpp, about line 1066), which is three hundred lines from either.
- The tensor k_output_values table takes its NCDM titles from a species virtual, but the
  module computes their values itself, after a `dynamic_cast<NCDMSpecies>` (the species
  picking STYLE.md forbids).

classy then splits the title string again (`SplitIntoColumns`, #503) into `NamedColumns`. Thomas
(2026-11-20): with the right data structure there would be no problem to begin with.

## Goals

- Each value is stored together with its name, in one call. Nothing else lists the columns.
- The column order is the order in which the names first appear, so the .dat files and
  classy's dicts keep today's order without anyone writing it down.
- A name set twice in one row, a name the table does not have after the first row, and a row
  that misses a column are all severe errors.
- Every .dat file is byte-identical to #504's. classy's tables are identical: same keys, same
  order, same values.
- Gone: `NamedColumns`, `SplitIntoColumns`, the title strings, `get_number_of_titles`,
  `first_duplicate_title`, `class_store_columntitle`, `class_store_double`,
  `class_store_double_or_default`, `BackgroundColumnWriter`, and the title mode of
  `PerturbColumnWriter`.

## Non-goals

- No change to a column name, the column order, a file format or a value.
- No change to the C_l and P(k) files, which have their own writers.

## Anti-goals

- A second list of names to keep in step with the code that makes the columns (for example
  "known names first, the rest after").
- `get_background()`, `get_thermodynamics()` or `get_perturbations()` more than twice as slow
  as on #504.

## Constraints

- STYLE.md: all species are equal; modules loop over species and dispatch, never pick a type.
- The k_output_values tables are filled during the integration, one k per thread.
- `class_fprintf_columntitle` and `class_fprintf_double` fix the .dat format.

## Approaches considered

| | Approach | Trade-off |
|---|---|---|
| A | `std::map<std::string, std::vector<double>>` | Sorts the names: `z` would be the last column of background.dat, while people plot by column number. Also keeps the first of two equal keys silently. Rejected with Thomas, 2026-11-20. |
| B | Columns declared first (`int c = table.AddColumn(name)`), rows filled by handle | Safe order, but every condition is still written twice, once to declare and once to fill. |
| **C** | **Rows filled by `Set(name, value)`; the first row defines the columns, later rows are checked against it** | **One place per column. Costs one hash lookup per value. Chosen.** |

## Design

### 1. `OutputTable`

`include/output_table.h`, `tools/output_table.cpp`. Like a Python dict: the names in a vector,
which keeps their insertion order, backed by an `std::unordered_map` from name to column, so
that every `Set` finds its column in O(1) (Thomas, 2026-11-20, with the per-row `Set`s of the
perturbation output in mind).

```cpp
class OutputTable {
 public:
  explicit OutputTable(std::string what);            // "background", "transfer functions", ...
  void Set(const std::string& name, double value);  // a value of the current row
  void EndRow();                                     // closes the current row
  void AddNote(const std::string& note);             // a header line; a repeated note is kept once
  const std::vector<std::string>& names() const;     // in insertion order
  const std::vector<double>& column(size_t index) const;
  const std::vector<double>& column(const std::string& name) const;  // O(1)
  size_t rows() const;
  const std::vector<std::string>& notes() const;
 private:
  std::vector<std::string> names_;
  std::unordered_map<std::string, size_t> index_;
  std::vector<std::vector<double>> columns_;
  std::vector<size_t> last_row_set_;  // per column: the row it last got a value in
};
```

- **First row:** a new name appends a column.
- **Later rows:** a name the table does not have is severe ("'x' is not a column of the
  background table"). A value goes to its column by name, so a row may set its names in any
  order.
- **Any row:** a name set twice in one row is severe, and so is `EndRow` with a column that got
  no value in that row (the per-column row stamp tells both).
- **Storage:** column by column, so classy and the .dat writer read whole columns.
- **Messages** name the table, so the reader knows which output failed.

### 2. Producers

| Today | Becomes |
|---|---|
| `background_output_titles` + `background_output_data` + `background_output_notes` | `OutputTable BackgroundModule::GetOutputTable() const` |
| `thermodynamics_output_titles` + `thermodynamics_output_data` | `OutputTable ThermodynamicsModule::GetOutputTable() const` |
| `primordial_output_titles` + `primordial_output_data` | `OutputTable PrimordialModule::GetOutputTable() const` |
| `perturb_output_titles` + `perturb_output_data` + `perturb_output_data_at_index_tau` + `perturb_store_columns` | `std::vector<OutputTable> GetTransferTables(file_format, double z)` and `GetTransferTablesAtIndexTau(file_format, int)`, one table per initial condition |
| `perturb_prepare_k_output` (titles at initialisation) + `perturb_print_variables` (rows) | rows appended in `perturb_print_variables` to one `OutputTable` per k_output_values and mode; no titles at initialisation |
| `GetOutputColumns`, `GetTransferColumns*`, `GetPerturbationColumns` (#503) | the tables above |

The CAMB-format `class_store_double_or_default` becomes `Set(name, condition ? value : 0.)`.

### 3. Species

| Today | Becomes |
|---|---|
| `WriteBackgroundColumnTitles(BackgroundColumnWriter&)` + `WriteBackgroundData(pvecback, BackgroundColumnWriter&)` | `virtual void AddBackgroundOutput(OutputTable& row, const double* pvecback) const` |
| `WriteOutputColumns(PerturbColumnWriter&, mod, fmt, section)` | `virtual void AddTransferOutput(OutputTable& row, const double* tk, const PerturbationsModule& mod, file_format fmt, TransferColumnSection section) const` |
| `PrintVariables(PerturbColumnWriter&, layout or null, tau, y, mod, ppw)` | `virtual void AddPerturbationOutput(OutputTable& row, const PerturbLayout& layout, double tau, const double* y, const PerturbationsModule& mod, const perturb_workspace* ppw) const` (no title mode, so the layout is never null) |
| `WriteTensorOutputColumnTitles(std::string&)` + the module's `dynamic_cast<NCDMSpecies>` loop | `virtual void AddTensorPerturbationOutput(OutputTable& row, const PerturbLayout& layout, double a, const double* y, const perturb_workspace* ppw) const`, the NCDM integrals moving into `NCDMSpecies` |

- Composites chain to their children as they do now.
- A species that writes a column only under a condition writes `if (condition) row.Set(...)`.
- Notes: `row.AddNote(...)` replaces `BackgroundColumnWriter::Note`.

### 4. Output and classy

- **.dat files:** `OutputModule::output_print_table(FILE*, const OutputTable&)` replaces
  `output_print_data`, writing the same bytes: `#`, one `class_fprintf_columntitle` per
  name, then one line per row.
- **classy:** the core converts a table into a dict (`_table_dict`) where it now converts
  `NamedColumns` (`_columns_dict`). `generate_wrapper.py` declares `OutputTable` to Cython
  as it declares `NamedColumns` today.

### 5. Repeated names fail when the table is built

Today the background, transfer and k_output_values titles are checked for repeats at
initialisation. With the titles gone, a repeat is found when a table is built: when a .dat
file is written, or classy asks for the table. *My call:* a repeated name matters only to
the table that holds it, so a run that never builds that table is not wrong. The error and
its message ("rename the instance") are unchanged.

## Verification

- **Unit test** `unittests/tools/output_table_test.cpp` (and the repeated-title checks of
  `unittests/species/scalar_field_multi_instance_test.cpp` move onto tables):
  - order of first appearance;
  - a repeated name in the first row;
  - a missing column and an unknown name in a later row, and a later row in another order;
  - an empty table;
  - notes kept once.
- **.dat files, byte for byte against #504:** run `./class` from both trees on a set of .ini
  files and compare every output file with `cmp`. The set covers:
  - background with many species: ncdm, fluid, two scalar fields, idm_dr, dcdm, dncdm,
    greybody;
  - thermodynamics with idm_dr and the damping scale;
  - primordial with scalars and tensors, and tensor-only;
  - mTk and vTk in CLASS and CAMB format, at z = 0 and z > 0;
  - k_output_values for scalars and tensors (with NCDM);
  - ad + cdi.
- **classy:** `benchmarks/classy_glue/compare_step2.py`, run against #504 rather than
  master, must find nothing.
- **Test suites:** ctest; TEST_LEVEL=1 scenarios; `test_class.py`; `test_package.py`; the
  level-2 CI job.
- **Timing:** `get_background()`, `get_thermodynamics()`, `get_perturbations()` and
  `get_transfer()` before and after (anti-goal: at most 2×).

## Left to the builder

- Whether hashing a name per value costs more than the anti-goal allows; if so, a faster path
  for a later row that sets its names in the first row's order (compare with the expected
  name before hashing).
- The order of commits. One per producer is natural, since each keeps the build green.
