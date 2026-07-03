# CLASSpp Documentation

CLASSpp is the C++ version of CLASS, the Cosmic Linear Anisotropy Solving
System. This site is meant to be a practical map of the codebase: start from
the conceptual pages below, then use the generated API reference when you need
the exact class, function, or source-level contract.

The code and tests remain the authoritative description of behavior. The
Doxygen manual should explain public contracts, ownership and lifetime rules,
units, numerical assumptions, and the reasons behind non-obvious algorithms.

## Start here

| Need | Go to |
| ---- | ----- |
| Understand what documentation exists | [Documentation map](md_chap2.html) |
| Find a C++ type | [Class list](annotated.html) |
| Inspect inheritance and the species model | [Graphical class hierarchy](inherits.html) |
| Browse source files | [File list](files.html) |
| Regenerate the manual locally | [Updating the manual](md_mod.html) |

## Visual maps

Graphviz diagrams are generated for inheritance and directory dependencies.
These are the best places to see them:

| Graph | Go to |
| ----- | ----- |
| Full graphical inheritance overview | [Graphical class hierarchy](inherits.html) |
| Species inheritance tree | [BaseSpecies](classBaseSpecies.html) |
| Non-cold species inheritance tree | [NCDMBaseSpecies](classNCDMBaseSpecies.html) |
| Species directory dependencies | [species directory](dir_5cf5c0309e5412b6c393145eca6c2225.html) |
| Source directory dependencies | [source directory](dir_b2f33c71d4aa5e7af42a1ca61ff5af1b.html) |

## Architecture map

The high-level computation still follows the CLASS pipeline, but current
CLASSpp code uses C++ ownership and module objects rather than the old C
`*_init()` / `*_free()` lifecycle.

| Area | Entry points |
| ---- | ------------ |
| Input parsing and normalization | [input_module.h](input__module_8h.html), [input_module.cpp](input__module_8cpp.html) |
| Background evolution | [background_module.cpp](background__module_8cpp.html) |
| Thermodynamics | [thermodynamics_module.cpp](thermodynamics__module_8cpp.html) |
| Perturbations and sources | [perturbations_module.cpp](perturbations__module_8cpp.html) |
| Transfers and spectra | [transfer_module.cpp](transfer__module_8cpp.html), [spectra_module.cpp](spectra__module_8cpp.html) |
| Species infrastructure | [SpeciesCollection](classSpeciesCollection.html), [BaseSpecies](classBaseSpecies.html), [SpeciesInput](classSpeciesInput.html) |

## Reading the reference

- Use **Classes** when you know the type you are looking for.
- Use **Files** when you are following the execution path through modules.
- Use search for input names, species names, and physical quantities.
- Treat undocumented generated pages as an index into the code, not as a
  substitute for the implementation.

## Documentation policy

Good Doxygen comments should capture information that the compiler cannot:
units, invariants, valid ranges, ownership rules, algorithmic choices, and
cross-module contracts. Comments that only repeat the function signature or the
next line of code should be removed or replaced with useful context.
