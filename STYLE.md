# STYLE GUIDE

This document defines the coding style and architectural principles for the codebase.
All contributions—human or automated—must follow these rules.

---

# 1. Core Principles

* **RAII everywhere**: All resources are acquired in constructors and released in destructors.
* **No partial states**: Objects are always fully initialized after construction.
* **Immutable shared inputs**: Input data may be normalized during initialization, but must not be modified after that point.
* **Explicit ownership**: Ownership must be expressed using smart pointers.
* **Deterministic behavior**: Avoid hidden state and implicit dependencies.

---

# 2. Architecture

## 2.1 Module Design

Each module is split into two parts:

1. **Input struct**

   * Contains input parameters and directly derived quantities.
   * Plain data only (no complex logic).
   * Mutable only during initialization, including canonicalization of user input into internal form.

2. **Module class**

   * Computes and stores derived quantities.
   * Owns all allocated resources.
   * Enforces invariants.

---

## 2.2 Construction Rules

* Modules must be fully initialized during construction.
* Default constructors are forbidden unless the object is trivially valid.
* The following pattern is **forbidden**:

```cpp
BackgroundModule ba;  // ❌ Forbidden
```

* Construction must always provide required dependencies:

```cpp
BackgroundModule ba(input_module);  // ✅ Required
```

---

## 2.3 RAII

* Resource allocation happens **only in constructors**.
* Resource deallocation happens **only in destructors**.
* Legacy patterns like `init()` / `free()` must not be used in new code.

---

## 2.4 Module Dependencies

* Modules depend on each other via:

```cpp
std::shared_ptr<const T>
```

* Input data is accessed through pointers inherited from `BaseModule`:

  * `ppr`, `pba`, `pth`, `ppt`, etc.

These pointers:

* must not be reassigned
* must not be used to modify input data after initialization

---

# 3. Ownership and Pointers

## 3.1 Rules

* Use `std::shared_ptr<const T>` for shared ownership.
* Use `std::unique_ptr<T>` for exclusive ownership.
* Raw pointers are **non-owning only**.

---

## 3.2 Invariants

* Ownership must always be clear from the type.
* No hidden ownership transfer.
* No manual memory management (`new` / `delete`) outside constructors.

---

# 4. Naming Conventions

All naming is enforced via tooling and must be followed strictly.

| Element             | Style       | Example            |
| ------------------- | ----------- | ------------------ |
| Variables           | snake_case  | `scale_factor`     |
| Function parameters | snake_case  | `omega_b`          |
| Struct members      | snake_case  | `omega_cdm`        |
| Class members       | snake_case_ | `lmax_`            |
| Classes / Structs   | PascalCase  | `BackgroundModule` |
| Methods             | PascalCase  | `ComputeSpectrum`  |
| Constants           | kCamelCase  | `kDaysInAWeek`     |
| Static constants    | sCamelCase  | `sCacheSize`       |

---

## 4.1 Rationale

* Trailing `_` prevents shadowing of class members.
* Naming must be consistent to support automated tooling.

---

# 5. Formatting Rules

Formatting is enforced via `.clang-format`.

Key rules:

* Indentation: **2 spaces**, no tabs
* Braces: K&R style for control statements, with `else` on a new line

```cpp
if (condition) {
  ...
}
else {
  ...
}
```

* Space after commas
* Space around `=` and comparison operators
* Consistent line wrapping

Do not manually format code—run the formatter.

---

# 6. Variables and Scope

## 6.1 Rules

* Variables must be initialized at declaration.
* Variables must be declared in the narrowest possible scope.

---

## 6.2 Forbidden

```cpp
double x;
if (cond) {
  x = compute();
}
```

---

## 6.3 Required

```cpp
if (cond) {
  double x = compute();
}
```

---

## 6.4 Rationale

* Prevents use of uninitialized variables
* Improves readability
* Reduces hidden dependencies

---

# 7. Const Correctness

* Prefer `const` wherever possible.
* Input data must be treated as immutable.
* Use `std::shared_ptr<const T>` for shared data.

---

# 8. Error Handling

Errors are C++ exceptions thrown at the error origin via the macros in
`include/errors.h`. The exception type encodes severity; the Python wrapper
maps it to sampler behaviour (`classy.pyx: raise_my_py_error`):

* `class_test_severe` / `class_stop_severe` → `std::invalid_argument` →
  `CosmoSevereError`: the sampler **aborts**. Use ONLY for checks that depend
  exclusively on the *structure* of the input: key presence, mutual
  exclusivity, unparseable strings, unknown enum names, list-length/count
  consistency, API-argument validation.
* `class_test` / `class_stop` → `std::runtime_error` →
  `CosmoComputationError`: the sampler **rejects the point** and the chain
  survives. Use for everything that depends on a parsed numeric value (range
  checks included) and for every numerical failure, wherever it occurs — the
  input phase and shooting trial builds included.
* The governing rule: a severe check must never depend on possibly-varying
  parameters. Aborting a long chain on one unlucky proposal destroys work;
  rejecting is always the safe direction.
* Programmer-error invariants throw `std::logic_error` (surfaces as
  `CosmoSevereError` via the classy fallback).
* Do not introduce alternative error mechanisms.

Design rationale: `docs/superpowers/specs/2026-07-17-error-severity-conventions-design.md`.

---

# 9. BaseModule

All modules (except `InputModule`) inherit from `BaseModule`.

`BaseModule` provides:

* Access to input structs (`pba`, `pth`, etc.)
* Ownership of the `InputModule`

Rules:

* Do not duplicate this functionality in derived classes.
* Do not pass input pointers manually between functions.

---

# 10. Forbidden Patterns

The following are not allowed:

* Partially initialized objects
* Calling methods on uninitialized modules
* Manual memory management outside RAII
* Mutable shared state without explicit ownership
* Declaring variables far from usage
* Shadowing class members without `_`

---

# 11. Required Patterns

* RAII-based construction
* Narrow variable scope
* Explicit ownership via smart pointers
* Consistent naming conventions
* Use of `BaseModule` for shared functionality

---

# 12. Tooling

All code must pass:

* `.clang-format`
* `.clang-tidy`

These tools are part of the definition of correctness.

---

# 13. For Automated Systems

* Follow all rules strictly.
* Do not introduce alternative styles.
* Do not introduce partially initialized objects.
* Prefer existing patterns over new abstractions.
* When in doubt, prioritize:

  1. RAII
  2. Const correctness
  3. Explicit ownership
  4. Simplicity

---

# 14. Design Documentation

Non-trivial work gets a **design spec**, committed under
`docs/superpowers/specs/` and named `YYYY-MM-DD-<topic>-design.md`. A spec states
the problem, the approach taken, the alternatives rejected and why, and how the
result is verified. Write it to be read a year later by someone who did not do
the work.

Step-by-step **implementation plans are not committed.** Write one if it helps
you execute, and keep it out of the repo. The commit history and the spec are the
durable record; a task list is scaffolding that goes stale the moment execution
diverges from it, and it is never updated to say so.

* Cite the spec from the code whose rationale it carries, by path, at the
  relevant declaration — see the file comments in `include/evolver_etd.h` and
  `species/dncdm_proxy_species.h`. Prose does the same where it needs to: the
  error-handling rules in section 8 close by naming the spec behind them.
* Record measured numbers in the spec, or in a companion
  `YYYY-MM-DD-<topic>-results.md`. Results nobody wrote down get re-measured.

## 14.1 A spec need not match the code

It must be honest about itself. Keeping a design document in sync with an evolving
implementation is a bargain you lose, and specs that pretend otherwise are the ones
nobody trusts. The requirement is narrower: **a spec must never misrepresent what
shipped.** A reader cannot tell a carried-out recommendation from a pending one, or
a dropped section from an implemented one, unless you say which it is.

Two conventions discharge that, both cheap, both already used here:

* **A status header.** One dated block at the top saying what shipped, what did
  not, and what superseded it. Correct an earlier spec *by reference* rather than
  editing it: `2026-08-12-dncdm-reduced-operator-galerkin-results.md` names the
  exact sections of its design note it overrides and the ones that still stand.
  A spec written as recommendations needs this most — an imperative table reads
  as pending work forever until a header says it was executed.
* **Evidence tags.** Mark what you did not verify yourself. The disposition tables
  in `2026-08-17-dncdm-minimal-knob-set.md` tag every row: **[M]** measured here,
  **[C]** stated by the code, **[P]** earlier campaign, not re-measured, **[?]** no
  evidence — do not act on this without checking. The **[?]** on `dr_rate_cap` is
  why that knob is still in the code while fifteen of the knobs tabulated
  alongside it are gone.

What ages well is rationale, rejected alternatives, and dated measurements: those
are facts about a moment and stay true. What ages badly is interface detail — knob
and symbol names, paths — which the code documents better anyway. Cite it lightly.

---

# 15. Guiding Principle

When making design decisions:

> Prefer correctness, clarity, and invariants over flexibility.
