# HYREC-2

The files in this directory are HYREC-2, vendored from upstream with a single
portability patch (see [Local modifications](#local-modifications)).

| | |
|---|---|
| Upstream | https://github.com/nanoomlee/HYREC-2 |
| Commit | `09e8243d0e08edd3603a94dfbc445ae06cafe139` (2023-09-13) |
| Authors | Yacine Ali-Haïmoud, Chris Hirata (2010–17), with contributions from Nanoom Lee (2020) |
| References | arXiv:1011.3758, arXiv:1012.3691, arXiv:2007.14114 |

## What is here, and what is not

Vendored: `energy_injection.{c,h}`, `helium.{c,h}`, `history.{c,h}`, `hydrogen.{c,h}`,
`hyrectools.{c,h}`, and the data tables `Alpha_inf.dat`, `R_inf.dat`, `fit_swift.dat`,
`two_photon_tables.dat`.

Not vendored: `hyrec.c` (the standalone driver's `main`), `input.dat`, `output_xe.dat`,
`two_photon_tables_hires.dat`, `readme.pdf`, and upstream's `Makefile`. CLASS builds these
sources as part of `libclasspp` and drives them through its own wrapper.

## Local modifications

One patch, in `history.c`, for MSVC (the Windows wheel build):

- `#include <unistd.h>` moved inside the `#ifdef CAMB` block, and the unused
  `#include <libgen.h>` deleted. Both headers are POSIX and do not exist in the MSVC C
  runtime, so the Windows wheel build failed with `error C1083: Cannot open include file:
  'unistd.h'`. The only declarations HYREC-2 takes from them are `getcwd()` and `chdir()`,
  used solely by `hyrec_init()` — which lives in the `CAMB` block that CLASS never
  compiles. Nothing from `libgen.h` is referenced at all.

Re-apply this when upgrading; it is otherwise not CLASS-specific and would be worth
sending upstream.

## Keeping the rest unmodified

Nothing else in this directory is edited. Everything CLASS-specific lives in
`source/hyrec_model.{h,cpp}`, behind the `RecombinationModel` interface in
`source/recombination_model.h`; that pair owns the `extern "C"` seam, the `HYREC_DATA`
lifetime, and the translation between CLASS's and HYREC-2's conventions. Upgrading is
therefore a matter of dropping in the new upstream files, re-applying the patch above, and
updating the commit hash.

Two consequences worth knowing when reading the code:

- **HYREC-2 calls `exit(1)` in four places** — `hyrec_xe`, `hyrec_Tm`, `rec_get_cosmoparam`,
  and `interp_Dfnu`. None is reachable from CLASS: we interpolate the history ourselves,
  never read the standalone driver's stdin, and never run `FULL` mode. `test-hyrec` pins
  this.
- **`history.h` sets `MODEL SWIFT`** at compile time. CLASS passes the model as a runtime
  argument instead and uses `SWIFT`, falling back to `PEEBLES` where HYREC-2's own guards
  require it. `FULL` is not available: it needs the radiative-transfer history indexed on
  HYREC-2's uniform `lna` grid, which CLASS's adaptive evolver cannot supply.

## Upstream equivalence

This tree is byte-identical to `class_public` v3.3.4's `external/HyRec2020/` except for
whitespace, one `malloc` added inside a CAMB-only wrapper that neither code calls, and the
`unistd.h` move above. All four data tables match exactly. That makes `class_public` usable
as a parity reference for this integration.

## History

Before 2026-08, CLASS++ shipped the November 2011 release of HyRec, compiled with
`#define MODEL RECFAST` — an effective three-level atom with a fudge factor, not radiative
transfer. See issue #396.
