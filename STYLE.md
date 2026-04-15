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

* Errors are stored in `error_message_`.
* This member may be declared `mutable`.
* Do not introduce alternative error handling mechanisms without strong justification.

---

# 9. BaseModule

All modules (except `InputModule`) inherit from `BaseModule`.

`BaseModule` provides:

* Access to input structs (`pba`, `pth`, etc.)
* Ownership of the `InputModule`
* Error handling via `error_message_`

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

# 14. Guiding Principle

When making design decisions:

> Prefer correctness, clarity, and invariants over flexibility.
