# Issue #314 — modern `FileContent` accessor API (`optional<T> get<T>`)

**Date:** 2026-06-18
**Issue:** [#314](https://github.com/AarhusCosmology/CLASSpp/issues/314) — *v4 prep: modern FileContent accessor API to replace `parser_read_*` out-params*
**Status:** design approved (brainstorming), ready for implementation plan

## Problem

`FileContent` is a modern class (std::map storage, read-tracking, `operator+`) but two C-era reading idioms persist on top of it:

1. **Free `parser_read_int/double/string(pfc, name, &val, &flag)`** with `int* found` out-params and the `_TRUE_/_FALSE_` flag convention — **98 call sites, all in `source/input_module.cpp` (93) + `source/input_module.h` (5).** No other module or species uses them.
2. **`input_module.h` macros** `class_read_double/int/string` (133 sites), `class_read_double_one_of_two` (30 sites), and the mutual-exclusion helpers `class_at_least_two_of_three` (3), `class_at_least_two_of_four` (1), `class_none_of_three` (2). All are built on **shared function-top locals** (`param1/param2/param3`, `int1`, `string1`, `flag1/flag2/flag3` — 301 + 256 references) threaded implicitly through the macros. This shared-local coupling is brittle: a macro silently depends on which locals exist at its call site.

A third, half-modern surface already exists: member accessors `FileContent::read_int/read_double/read_string(name, val&) → bool` and `read_list_of_doubles/integers/strings`, used at **45 sites in `species/` and other modules**, defined in `tools/parser.cpp`.

The goal: one modern, composable accessor surface; retire the C idioms and the shared-local coupling entirely.

## Decisions (from brainstorming)

- **Scope:** full modernize. `input_module.cpp` ends with no `class_read_*` macros, no shared `param/int/string/flag` locals, no `_TRUE_/_FALSE_`.
- **API shape:** template `get<T>` returning `std::optional<T>`, plus `get_or<T>`. Single primary surface.
- **Legacy surfaces — full sweep:** `get<T>` becomes the *only* `FileContent` accessor; `read_*` members are deleted and every caller migrates to `get<T>`. The unused free `parser_read_int/double/string` are deleted. `parser_read_file` / `parser_cat` are kept (file loading / concatenation — different concern, still used).
- **`SpeciesInput` reshaped too (full consistency):** the per-instance prefixing wrapper's `read_*` / `required_*` API is replaced by a `get<T>` / `require<T>` surface mirroring `FileContent`, and all 35 `input.read_*` species call sites migrate. (Chosen over keeping its API; accepted that #246 may later rebuild this layer.)
- **Macros:** deleted, not rebodied. The simple reads become one-line `get_or`; the one-of-two / mutual-exclusion patterns become small free helpers in `input_module` (input *policy* does not belong on the generic container).

## The accessor API (`include/parser.h`, `tools/parser.cpp`)

```cpp
class FileContent {
  // ...
  // Returns the value if the parameter is present (marking it read), nullopt if absent.
  // Throws (as today) on a present-but-unparseable value.
  template <class T> std::optional<T> get(const std::string& name) const;

  // Read-with-default convenience: get<T>(name).value_or(fallback).
  template <class T> T get_or(const std::string& name, T fallback) const {
    return get<T>(name).value_or(fallback);
  }
  // Overload so a string-literal fallback does not deduce T = const char*.
  std::string get_or(const std::string& name, const char* fallback) const {
    return get<std::string>(name).value_or(fallback);
  }
};
```

- **Supported `T` (explicit specializations in `tools/parser.cpp`):** `int`, `double`, `std::string`, `std::vector<int>`, `std::vector<double>`, `std::vector<std::string>`. The list specializations subsume the old `read_list_of_*`.
- The **primary template is left undefined** (or `static_assert(sizeof(T)==0, ...)` in a fallback) so an unsupported `T` is a compile/link error, not silent breakage.
- **Read-tracking** (`read_params_`) moves into `get<T>`; `was_read` / `mark_all_unread` / `unread_parameters` / `is_shooting` behave exactly as today. `get_or` marks read on a hit only (same as a present read today); on a miss nothing is marked, matching current `class_read_*` semantics.
- **Parse-error behavior unchanged:** a present value that fails to parse throws the same exception as the current `read_*` / `parser_read_*`.

## `input_module` migration

Two small free helpers (file-local, e.g. anonymous namespace in `input_module.cpp`; declared in `input_module.h` only if needed across its TUs):

```cpp
// Number of the given optionals that hold a value.
template <class... Opt> int n_present(const Opt&... opts) {
  return (0 + ... + (opts.has_value() ? 1 : 0));
}

// Value of whichever of two names is present; class_test error if both present; nullopt if neither.
template <class T>
std::optional<T> read_one_of(const FileContent& fc, const char* n1, const char* n2) {
  auto a = fc.get<T>(n1);
  auto b = fc.get<T>(n2);
  class_test(a && b, "In input file, you can only enter one of %s, %s, choose one", n1, n2);
  return a ? a : b;
}
```

Call-site transformations:

| Old | New |
|-----|-----|
| `class_read_double("x", dest);` | `dest = pfc->get_or("x", dest);` |
| `class_read_int("x", dest);` (int dest) | `dest = pfc->get_or("x", dest);` |
| `class_read_int("extrapolation_method", pnl->extrapolation_method);` (enum dest — the one case) | `pnl->extrapolation_method = (enum source_extrapolation) pfc->get_or<int>("extrapolation_method", pnl->extrapolation_method);` |
| `class_read_string("x", dest);` | `dest = pfc->get_or("x", dest);` |
| `class_read_double_one_of_two(n1, n2, dest);` | `if (auto v = read_one_of<double>(*pfc, n1, n2)) dest = *v;` |
| `parser_read_double(pfc, "N_idr", &param1, &flag1); ... if (flag1 == _TRUE_) ...` | `auto N_idr = pfc->get<double>("N_idr"); ... if (N_idr) ...` |
| `class_test(class_at_least_two_of_three(flag1, flag2, flag3), ...)` | `class_test(n_present(a, b, c) > 1, ...)` |
| `class_none_of_three(flag1, flag2, flag3)` | `n_present(a, b, c) == 0` |

After migration, delete from `input_module.h`: `class_read_double`, `class_read_int`, `class_read_string`, `class_read_double_one_of_two`, `class_at_least_two_of_three`, `class_at_least_two_of_four`, `class_none_of_three`. (`class_any_nonzero_four` / `class_all_nonzero_four` already have 0 uses — delete too.) Delete the shared `param1/param2/param3/int1/string1/flag1/flag2/flag3` locals from every Read function.

## Full sweep — the `read_*` callers

The `read_*` callers fall into three groups (counts are approximate, file-by-file in the plan):

**(a) Direct `FileContent::read_*` callers (~55 sites)** — `ctx.pfc->read_double(...)`, `fc.read_*`, `pfc.read_*` across `input_module.cpp` and 9 species files (`dcdm_dr`, `fluid`, `idm_dr_idr`, `ncdm`, `ncdm_interacting`, `dncdm`, `scalar_field`, `ultra_relativistic`, `greybody_ncdm`). Mechanical migration:
- `bool has_G = ctx.pfc->read_double("G_eff", v);` → `auto g = ctx.pfc->get<double>("G_eff");` then use `g.has_value()` / `*g`. The separate `v` is dropped where it existed only for the out-param.
- `if (pfc.read_string(key, unused))` → `if (pfc.get<std::string>(key))`.
- `ctx.pfc->read_list_of_doubles(name, v)` → `ctx.pfc->get<std::vector<double>>(name)`.

**(b) `input_module.cpp` local helpers** (not the macro smell — keep, rewire bodies):
- The three `read(const FileContent&, name, T& v)` overloads + `read_enum` (175 call sites) become `v = fc.get_or(name, v);` / `v = (E) fc.get_or<int>(name, (int)v);` in their bodies. **Call sites unchanged.**
- The file-local `readDoubleList(pfc, name, values, &found)` adapter → reimplemented on `get<std::vector<double>>` (or removed in favor of `get<>` directly where it is called).

**(c) `SpeciesInput` reshape.** Replace its public surface with templates that prefix and delegate to `FileContent::get<T>`:
```cpp
class SpeciesInput {
  template <class T> std::optional<T> get(const std::string& field) const;        // pfc_->get<T>(qualify(field))
  template <class T> T require(const std::string& field) const;                   // throws std::invalid_argument naming instance+field if absent
};
```
Delete `SpeciesInput::read_double/read_int/read_string/read_list_of_doubles` and `required_double/required_int/required_string`. Migrate all 35 `input.read_*` species call sites to `input.get<T>` and the (test-only) `required_*` to `require<T>`. Update `tools/parser_test.cpp` accordingly and add direct `FileContent::get<T>` coverage.

**Deletions:**
- `FileContent::read_int/read_double/read_string/read_list_of_doubles/read_list_of_integers/read_list_of_strings`.
- Free `parser_read_int/double/string` and their declarations in `parser.h` (the `extern "C"` block keeps only `parser_read_file`; `parser_cat` stays; the non-`extern "C"` `parser_read_string` decl is removed).

**Cython:** confirmed `python/classy.pyx` and `python/cclassy.pxd` reference none of `read_*` / `parser_read_int/double/string` / `get`, so no wrapper changes are needed. The plan still re-greps before deleting (lesson from #320: a Cython reference only surfaces when `classy.cpp` is regenerated).

## Verification

This refactor changes *how* parameter values are read, not *which* values are read or any computed number, so a correct migration is **byte-identical in output**. That makes bit-identity a cheap, strong regression witness here (in contrast to the physics-touching v4 refactors where bit-identity was the wrong bar — see the project norm).

- Build clean (all three manifests if applicable; primary is CMake).
- Full `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` green against the current `classyref`.
- Targeted spot-check on `.ini` inputs exercising the migrated branches: the IDR `N_idr/N_dg/xi_idr` mutual-exclusion, a `class_read_double_one_of_two` parameter, and a list parameter (e.g. `z_pk`), confirming identical parsed values and identical run output.
- Acceptance is correctness; byte-identity is the witness, not the contract.

## Forward-compatibility (#246)

#246 (object-based species input) will rewrite much of the input-parsing surface. This API is the agreed foundation: both `FileContent` and the reshaped `SpeciesInput` expose `get<T>` (+ `get_or` / `require<T>`), never the deleted C shims. `get<T>` returning `optional<T>` composes naturally with per-species input objects (a field is present-or-absent), so #246 builds on this accessor rather than replacing it. Reshaping `SpeciesInput` now (rather than keeping its `read_*`/`required_*` API) was chosen with eyes open: if #246 rebuilds the layer, it starts from the consistent `get<T>` surface.

## Out of scope

- No change to `parser_read_file` / `parser_cat` (file I/O / concatenation).
- No change to the `.ini` file format or parameter names.
- The ~130 hard-coded species *name* branches in modules (`all_species_.count("Fluid")`, …) are a separate concern, not touched here.
