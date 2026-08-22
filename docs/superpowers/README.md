# Design documentation

`specs/` holds the design specs for non-trivial work in this repo, one per topic,
named `YYYY-MM-DD-<topic>-design.md`. A spec states the problem, the approach
taken, the alternatives rejected and why, and how the result is verified.
Measured numbers live in the spec or in a companion
`YYYY-MM-DD-<topic>-results.md`.

The convention is written up in `STYLE.md` section 14. In short:

* Specs are committed. They are cited by path from the code they explain — see
  `include/evolver_etd.h`, `species/dncdm_proxy_species.h`, `STYLE.md` section 8.
* Step-by-step implementation plans are **not** committed. Write one if it helps
  you execute; keep it out of the repo.
* A spec is not required to match the code, only to be honest about itself — see
  `STYLE.md` section 14.1. Say what shipped in a dated status header, and tag
  claims you did not verify. `2026-08-12-dncdm-reduced-operator-galerkin-results.md`
  and `2026-08-17-dncdm-minimal-knob-set.md` are the worked examples of each.

## Where `plans/` went

This directory used to carry a `plans/` tree: 62 step-by-step task lists, one per
piece of work. They were removed because they had no forward value. Every one of
them was written before execution and committed exactly once — never revisited,
never corrected when execution diverged from the plan. Nothing in the codebase
depended on them: of the ten places that cited this directory, eight cited a
spec, and the two that named a plan were repointed at the corresponding spec when
`plans/` was deleted.

The plans remain in git history. To read one:

```sh
# find the commit that removed them
git log --diff-filter=D --oneline -- docs/superpowers/plans

# list what was there, and print one
git show <removal-commit>^ --stat -- docs/superpowers/plans
git show <removal-commit>^:docs/superpowers/plans/<file>.md
```

A few archived specs still contain inline references to `plans/…` paths. Those
are left as written: they are historical records of what was believed at the
time, and one of them already flags its own references as stale. Follow the
recipe above if you need the target.
