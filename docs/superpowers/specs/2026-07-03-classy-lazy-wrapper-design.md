# Design: Lazy dirty-check gate in the classy.pyx wrapper

**Date:** 2026-07-03
**Status:** Approved
**Branch:** `classy-lazy-gate`

## Problem

The C++ core is lazy and constructor-driven: `Cosmology(fc)` parses input cheaply,
`GetInputModule()` runs shooting on first use (`cosmology.cpp`), and every
`GetXxxModule()` builds its module on demand. The Cython wrapper predates this and
still follows the C `set()`/`compute()` lifecycle it inherited for MontePython
compatibility. The mismatch produces three concrete defects in `classy.pyx`:

1. **Silent staleness.** `set()` only flips `parameters_changed`; `_thisptr` still
   points at the old `Cosmology`. Any method called after `set()` but before
   `compute()` silently returns results for the *old* parameters.
2. **Eager shooting at construction.** `__init__` → `reset()` →
   `GetInputModule()` runs shooting inside the constructor, even if no method is
   ever called. The cause is the cache of ten raw struct pointers
   (`self.pr`, `self.ba`, … `self.op`): `DoShooting` *replaces* the `InputModule`
   instance, so the pointers can only be cached post-shooting, which forces
   `reset()` to shoot eagerly.
3. `compute()`'s level dispatch and the `state` property are vestigial shims whose
   semantics are only honest in the one MontePython flow they were shaped around.

## Goal

Make the wrapper as lazy as the C++ core: construction is cheap parse-only
validation, any method call always reflects the current parameter dict, and
`set()`/`compute()` remain as fully compatible legacy shims with unchanged error
timing for MontePython/Cobaya/CosmoHammer.

Non-goals: no C++ changes (the core already supports this), no API removals, no
hand edits to the generated `cclassy.pxd`.

## Design

### 1. The gate

A single accessor becomes the only way to reach the C++ object:

```cython
cdef inline Cosmology* _cosmo(self) except NULL:
    if self.parameters_changed or not self._thisptr:
        self.reset()
    return self._thisptr.get()
```

All ~118 `deref(self._thisptr)` sites become `deref(self._cosmo())`
(mechanical replacement). `_thisptr` is never dereferenced anywhere else.

### 2. Remove the struct-pointer cache

The ten cached members `self.pr, self.ba, self.th, self.pt, self.pm, self.nl,
self.tr, self.sp, self.le, self.op` are deleted and replaced by inline accessors
that fetch through the gate on every use:

```cython
cdef inline const background* ba(self):
    return &deref(deref(self._cosmo()).GetInputModule()).background_
```

Call sites change `self.ba.Omega0_k` → `self.ba().Omega0_k` (Cython auto-derefs
pointer member access; one sed per struct, ~56 sites total). These accessors
trigger shooting via `GetInputModule()` by design: the post-shooting structs are
the only ones safe to read, because `DoShooting` replaces the `InputModule`.

Perf note: each access costs two null-checks instead of a raw cached pointer.
These are Python-facing convenience getters, not hot loops; the cost is
irrelevant. (Considered and rejected: re-caching post-shooting pointers keyed to
the current `Cosmology` instance — complexity without measurable benefit.)

### 3. `reset()` slims down

Keeps: rebuild `_fc` from `_pars`, construct `new Cosmology(self._fc)`, run the
unread-parameters check (all parse-time, ~ms, no shooting), clear
`parameters_changed`. Loses: the `GetInputModule()` call and the ten pointer
assignments.

### 4. Lifecycle methods

- `__init__(params)`: stores `_pars`, calls `reset()` once. Constructing with bad
  parameters still raises immediately (parse-time validation) but no longer
  shoots. "If you can construct it, you can call its methods" holds for real.
- `set()` / `empty()`: unchanged — dict update + dirty flag. Repeated `set()`
  calls stay free.
- `compute(level)`: drops its explicit `if parameters_changed: reset()` (the gate
  handles it); keeps the level dispatch routed through `_cosmo()`. Documented as
  a legacy shim whose remaining purpose is "trigger everything now so errors
  surface here" — preserving MontePython's error timing exactly.
- `struct_cleanup()`: unchanged no-op shim.
- `state`: unchanged (`True`) — now honest, since any access is guaranteed
  current.
- CosmoHammer `__call__`: replace its direct `self.reset()` with setting the
  dirty flag; its `compute()` call gates anyway.

### 5. Error semantics (unchanged where it matters)

| Error class | Constructor-style use | MontePython-style use |
|---|---|---|
| Parse-time (unread param, contradictory inputs) | at `__init__` | at `compute()` / first access after `set()` — same as today |
| Shooting / computation failure | at first method call | at `compute()` — same as today |

## Testing

Add to `python/test_class.py`:

1. **Staleness regression:** `compute()`, read a Cl; `set({'h': …})`; read again
   *without* `compute()`; assert the result changed.
2. **Compute-less use:** constructor then method call with no `compute()` ever.
3. **Error timing:** `set()` with an unread garbage parameter raises at
   `compute()` (or first access), not at `set()`.
4. **MontePython-style flow:** set/set/compute/methods sequence matches current
   behavior.

Existing wrapper tests must pass. Numerical verification against master uses the
usual ~0.1%-tolerance comparison (never bit-identity; Cl^TE zero-crossings
handled).

## Files touched

- `classy.pyx` — all changes live here.
- `python/test_class.py` — new regression tests.
